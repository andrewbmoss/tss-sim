"""Layer: L2
WP: WP3
Claims: C1, C4
Status: implemented
Requirements: INS-11, INS-12

The receive chain: gain, noise figure, dynamic range clipping and quantisation.

Simple on purpose.  This is a lumped model of the front end, not a circuit
simulation: one gain, one noise figure, one clipping level and one quantiser.
What it exists to do is make the two things that decide whether the processing
works appear in the samples rather than in a spreadsheet.

**Dynamic range is a real constraint here and not a formality.**  The instrument
transmits while it receives, so the strongest thing at the receiver input is the
transmitter, not the echo.  INS-11's 60 dB instantaneous dynamic range is what
has to hold the leakage and the echo at the same time.  A receiver model that
cannot clip cannot show what happens when it does not, and the whole of C4 is
about what happens when it does not.

**Quantisation sets a floor the processing cannot get below.**  With the
external-noise-limited design of this instrument the right operating point is
the one where the noise itself dithers the quantiser, which needs the noise to
sit a few least significant bits above the quantisation step.  Too much gain and
the leakage clips; too little and the quantiser, not the sky, sets the floor.
Both ends are visible in this model and that is what it is for.

CONVENTIONS
-----------
Samples carry the square root of power in watts, so ``abs(x)**2`` is watts.
That is the same convention :mod:`tss.l2_instrument.waveform` uses.  Levels are
given in dBm because that is how a front end is specified, and the conversion
happens at the boundary of this module rather than being sprinkled through it.

Noise is complex circular Gaussian.  For a complex baseband sample stream at
sample rate ``fs``, the noise power in the sampled band is ``k T fs``, split
equally between the real and imaginary parts.  Randomness comes only from the
:class:`numpy.random.Generator` passed in, per contract rule 4.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import numpy as np

from tss.constants import K_BOLTZ, T0_K
from tss.units import db, dbm_to_w, undb, w_to_dbm

__all__ = [
    "Receiver",
    "noise_temperature_from_nf_db",
    "thermal_noise_power_w",
]


def noise_temperature_from_nf_db(noise_figure_db: float, t0_k: float = T0_K) -> float:
    """Receiver noise temperature, K, from noise figure.  ``(F - 1) T0``."""
    return float((undb(noise_figure_db) - 1.0) * t0_k)


def thermal_noise_power_w(t_sys_k: float, bandwidth_hz: float) -> float:
    """``k T B`` in watts."""
    return float(K_BOLTZ * float(t_sys_k) * float(bandwidth_hz))


class Receiver:
    """A lumped receive chain.

    Parameters
    ----------
    gain_db
        Voltage-referred power gain of the chain ahead of the quantiser.  It
        does not change the signal-to-noise ratio and it is not there to; it is
        there because it decides where the signal sits against the clipping
        level and the quantisation step, and those two do change the answer.
    noise_figure_db
        Noise figure of the chain.  INS-12's threshold is 3 dB.
    instantaneous_dynamic_range_db
        The range the chain covers without changing gain, INS-11.  Used to place
        the clipping level relative to the noise floor when no explicit
        full-scale level is given.
    gain_controlled_dynamic_range_db
        The range covered when gain switching is allowed.  Recorded and reported;
        this model does not switch gain, because when gain switching happens and
        how fast it settles is a design question WP3 has not answered, and a
        model that switched instantaneously would hide the question.
    full_scale_dbm
        Explicit clipping level.  When ``None`` it is placed
        ``instantaneous_dynamic_range_db`` above the noise floor at the time
        :meth:`apply` is called, which is the operating point an automatic level
        control would seek.
    n_bits
        Quantiser resolution.  ``None`` disables quantisation, which is useful
        for isolating whether an effect is quantisation or something else.
    """

    def __init__(
        self,
        *,
        gain_db: float = 0.0,
        noise_figure_db: float = 3.0,
        instantaneous_dynamic_range_db: float = 60.0,
        gain_controlled_dynamic_range_db: float = 90.0,
        full_scale_dbm: float | None = None,
        n_bits: int | None = 14,
        channels: int = 1,
    ) -> None:
        if channels < 1:
            raise ValueError("channels must be at least 1")
        if n_bits is not None and n_bits < 2:
            raise ValueError("n_bits must be at least 2, or None to disable quantisation")
        self.gain_db = float(gain_db)
        self.noise_figure_db = float(noise_figure_db)
        self.instantaneous_dynamic_range_db = float(instantaneous_dynamic_range_db)
        self.gain_controlled_dynamic_range_db = float(gain_controlled_dynamic_range_db)
        self.full_scale_dbm = None if full_scale_dbm is None else float(full_scale_dbm)
        self.n_bits = None if n_bits is None else int(n_bits)
        self.channels = int(channels)

    # -- noise ---------------------------------------------------------------

    @property
    def noise_temperature_k(self) -> float:
        """``(F - 1) T0``, the receiver's own contribution."""
        return noise_temperature_from_nf_db(self.noise_figure_db)

    def system_temperature_k(self, antenna_temperature_k: float) -> float:
        """Antenna temperature plus the receiver's own.

        The antenna temperature is the sky-dominated one from
        :mod:`tss.l2_instrument.antenna`, not a physical temperature.  Passing a
        physical temperature here is the mistake the legacy model made and it is
        worth 20 to 35 dB.
        """
        return float(antenna_temperature_k) + self.noise_temperature_k

    def add_noise(
        self,
        samples: np.ndarray,
        fs_hz: float,
        antenna_temperature_k: float,
        rng: np.random.Generator,
    ) -> np.ndarray:
        """Add complex circular Gaussian noise for the sampled bandwidth.

        Returns a new array; the input is not modified.  The total noise power
        is ``k T_sys fs``, so each of the real and imaginary parts gets half of
        it, which is a standard deviation of ``sqrt(k T_sys fs / 2)``.
        """
        arr = np.asarray(samples)
        power_w = thermal_noise_power_w(self.system_temperature_k(antenna_temperature_k), fs_hz)
        sigma = np.sqrt(power_w / 2.0)
        noise = rng.normal(0.0, sigma, size=arr.shape) + 1j * rng.normal(
            0.0, sigma, size=arr.shape
        )
        return np.asarray(arr, dtype=complex) + noise

    # -- level handling ------------------------------------------------------

    def apply_gain(self, samples: np.ndarray) -> np.ndarray:
        """Scale by the chain gain.  Amplitude scales as the square root."""
        return np.asarray(samples, dtype=complex) * np.sqrt(undb(self.gain_db))

    def clip_level_w(self, noise_floor_w: float | None = None) -> float:
        """Full-scale power in watts, referred to the sample stream.

        Either the explicit ``full_scale_dbm``, or the noise floor raised by the
        instantaneous dynamic range.
        """
        if self.full_scale_dbm is not None:
            return float(dbm_to_w(self.full_scale_dbm))
        if noise_floor_w is None or noise_floor_w <= 0.0:
            raise ValueError(
                "with no full_scale_dbm the clipping level is placed above the "
                "noise floor, so a positive noise_floor_w is required"
            )
        return float(noise_floor_w) * float(undb(self.instantaneous_dynamic_range_db))

    def clip(self, samples: np.ndarray, clip_level_w: float) -> tuple[np.ndarray, int]:
        """Hard-limit the magnitude, preserving phase.

        Returns the clipped array and the number of samples that were clipped.
        Phase is preserved because a saturating amplifier limits amplitude and
        not phase, and because a clipper that also destroyed phase would make
        the de-chirped residual look worse than it is for the wrong reason.
        """
        arr = np.asarray(samples, dtype=complex)
        limit = np.sqrt(float(clip_level_w))
        magnitude = np.abs(arr)
        over = magnitude > limit
        n_clipped = int(np.count_nonzero(over))
        if n_clipped == 0:
            return arr.copy(), 0
        out = arr.copy()
        scale = np.ones_like(magnitude)
        np.divide(limit, magnitude, out=scale, where=over)
        out[over] = arr[over] * scale[over]
        return out, n_clipped

    def quantise(self, samples: np.ndarray, clip_level_w: float) -> np.ndarray:
        """Uniform mid-tread quantisation of the real and imaginary parts.

        The step is set by the full-scale level and the bit depth, with one bit
        given to the sign, so the step is ``2 * sqrt(P_fs) / 2**n_bits``.  A
        sample beyond full scale is limited to the last code, which is what a
        real converter does and is not the same thing as wrapping.
        """
        if self.n_bits is None:
            return np.asarray(samples, dtype=complex).copy()
        arr = np.asarray(samples, dtype=complex)
        full_scale = np.sqrt(float(clip_level_w))
        levels = 2 ** self.n_bits
        step = 2.0 * full_scale / levels
        limit = full_scale - step / 2.0

        def _q(x: np.ndarray) -> np.ndarray:
            return np.clip(np.round(x / step) * step, -limit, limit)

        return _q(arr.real) + 1j * _q(arr.imag)

    # -- the whole chain -----------------------------------------------------

    def apply(
        self,
        samples: np.ndarray,
        fs_hz: float,
        antenna_temperature_k: float,
        rng: np.random.Generator,
        *,
        quantise: bool = True,
    ) -> tuple[np.ndarray, dict[str, Any]]:
        """Noise, gain, clipping and quantisation, in that order.

        The order is the physical one.  Noise enters at the input, gain follows,
        and the converter sees whatever the gain produced.  Adding noise after
        the gain would make the noise figure meaningless and is the easiest way
        to get a receiver model wrong.

        Returns the samples and a flags mapping carrying the noise floor, the
        clipping level, how many samples clipped and the quantisation step, all
        of which go into ``ReceivedSignal.flags`` so that a run records the
        operating point it was actually at.
        """
        t_sys = self.system_temperature_k(antenna_temperature_k)
        noisy = self.add_noise(samples, fs_hz, antenna_temperature_k, rng)
        amplified = self.apply_gain(noisy)

        noise_floor_w = thermal_noise_power_w(t_sys, fs_hz) * float(undb(self.gain_db))
        level_w = self.clip_level_w(noise_floor_w)
        clipped, n_clipped = self.clip(amplified, level_w)
        out = self.quantise(clipped, level_w) if quantise else clipped

        flags: dict[str, Any] = {
            "system_temperature_k": t_sys,
            "noise_floor_dbm": float(w_to_dbm(noise_floor_w)),
            "clip_level_dbm": float(w_to_dbm(level_w)),
            "n_clipped": n_clipped,
            "clipped_fraction": n_clipped / max(int(np.asarray(out).size), 1),
            "n_bits": self.n_bits,
            "quantised": bool(quantise and self.n_bits is not None),
            "peak_dbm": float(w_to_dbm(float(np.max(np.abs(amplified)) ** 2)))
            if np.asarray(amplified).size
            else float("-inf"),
        }
        if self.n_bits is not None:
            flags["quantisation_step_w"] = float((2.0 * np.sqrt(level_w) / 2**self.n_bits) ** 2)
            flags["noise_to_lsb_db"] = float(
                db(noise_floor_w / flags["quantisation_step_w"])
            )
        return out, flags

    def describe(self) -> Mapping[str, Any]:
        """JSON-serialisable record for the run manifest."""
        return {
            "gain_db": self.gain_db,
            "noise_figure_db": self.noise_figure_db,
            "noise_temperature_k": self.noise_temperature_k,
            "instantaneous_dynamic_range_db": self.instantaneous_dynamic_range_db,
            "gain_controlled_dynamic_range_db": self.gain_controlled_dynamic_range_db,
            "full_scale_dbm": self.full_scale_dbm,
            "n_bits": self.n_bits,
            "channels": self.channels,
            "gain_switching": "not modelled; see the class docstring",
        }

    @classmethod
    def from_config(cls, cfg: Mapping[str, Any]) -> Receiver:
        """Build from the ``receiver:`` block of an instrument configuration."""
        return cls(
            gain_db=float(cfg.get("gain_db", 0.0)),
            noise_figure_db=float(cfg.get("noise_figure_db", 3.0)),
            instantaneous_dynamic_range_db=float(
                cfg.get("instantaneous_dynamic_range_db", 60.0)
            ),
            gain_controlled_dynamic_range_db=float(
                cfg.get("gain_controlled_dynamic_range_db", 90.0)
            ),
            full_scale_dbm=cfg.get("full_scale_dbm"),
            n_bits=cfg.get("n_bits", 14),
            channels=int(cfg.get("channels", 1)),
        )
