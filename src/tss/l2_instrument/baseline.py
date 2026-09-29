"""Layer: L2
WP: WP3
Claims: C1, C4
Status: partial
Requirements: INS-11, INS-12, INS-14, INS-21, INS-26

The baseline L2 tier: echoes in, complex baseband samples out.

This is the tier registered as ``l2/baseline`` and named by
``config/default.yaml``.  It composes the other modules in this package rather
than computing anything itself, which is the point: the antenna model, the noise
reference, the sweep law, the receive chain and the oscillator are each testable
on their own, and this module is the wiring.

WHAT IT PRODUCES
----------------
De-chirped complex baseband, one record per sweep.  The instrument de-chirps in
hardware against the same reference that generated the transmission, which is
INS-21, so the samples that leave the receiver are already at baseband and an
echo appears as a tone rather than as a delayed sweep.  Generating the passband
signal and then mixing it down in software would be a more literal simulation
and would cost a factor of several hundred in sample rate for no change in the
answer, so it is not done.

The de-chirped phase of an echo at two-way group delay ``tau`` is exact rather
than approximated:

    phi_b(t) = 2 pi [ f0 tau + d K tau t - d K tau^2 / 2 ]

from subtracting the transmit phase law of :mod:`tss.l2_instrument.waveform` at
``t - tau`` from the same law at ``t``.  The middle term is the beat frequency,
``d K tau``, which is negative for a down-sweep; carrying the sign rather than
the magnitude is what lets an up/down pair be processed as a pair, which C3
needs.  The last term is a constant phase and is the piece a first-order
treatment drops.

WHEN AN ECHO IS PRESENT
-----------------------
An echo at sounding frequency ``f_k`` is only there while the sweep is near
``f_k``.  The tier gives each row a rectangular window centred on the time the
sweep crosses that row's frequency, of width equal to the frequency spacing of
the echo set divided by the sweep rate.  That is the dwell the sounding actually
gets at that frequency, and it is why a fine frequency step costs dwell.

The window is rectangular and not tapered.  A taper would put a more realistic
skirt on each tone, and it would also be a processing choice made at L2 where it
belongs at L3, so the taper is left for the transform to apply.

ABSORPTION IS NOT APPLIED HERE
------------------------------
``EchoSet.amplitude_db`` from the L1 analytic tier already includes the
absorption integral; ``EchoSet.absorption_db`` is carried beside it for
traceability, not as a second loss to apply.  This tier therefore uses
``amplitude_db`` and records ``absorption_db`` in the flags.  Applying both
would double-count the absorption, which at the top of the band near foF2 is
several decibels.

WHAT IS NOT IN IT YET
---------------------
Transmit-to-receive leakage.  :mod:`tss.l2_instrument.leakage` is a stub and it
is on the critical path, so the samples this tier produces have no leakage in
them and the receiver therefore never clips on leakage.  That is the single
largest gap between this tier and the instrument, it is recorded in every
``ReceivedSignal`` this tier returns under the flag ``leakage_modelled``, and no
verdict on INS-11 or INS-24 may be taken from a run made without it.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import numpy as np

from tss.constants import K_BOLTZ
from tss.interfaces import EchoSet, ReceivedSignal, SweepSchedule
from tss.l2_instrument.antenna import ShortDipole, antenna_noise_temperature_k
from tss.l2_instrument.noise import (
    is_extrapolated,
    noise_density_dbm_per_hz,
    sky_temperature,
)
from tss.l2_instrument.oscillator import Oscillator
from tss.l2_instrument.receiver import Receiver
from tss.l2_instrument.waveform import (
    sweep_bandwidth_hz,
    sweep_duration_s,
    validate_schedule,
)
from tss.registry import register
from tss.units import db, undb, w_to_dbm

__all__ = ["BaselineInstrument", "DEFAULT_POLARIZATION_LOSS_DB"]

#: Polarization loss, dB.  -3 dB for the O and X split and -3 dB for an
#: elliptically polarized wave received on a linear antenna.  The same figure
#: the link budget uses as its default.  It is not in the instrument
#: configuration because it is a property of the propagation and the antenna
#: together rather than of the instrument alone.
DEFAULT_POLARIZATION_LOSS_DB = -6.0


@register("l2", "baseline")
class BaselineInstrument:
    """Compose the L2 modules into an :class:`~tss.interfaces.L2Instrument`.

    Parameters
    ----------
    tx_power_w
        Transmit power at the power amplifier output.  INS-14's threshold is 1 W.
    tx_chain_efficiency_db
        Loss between the amplifier and the antenna terminals, negative.  The
        baseline configuration states -20 dB and calls it matching-network
        dominated.  See the note in :meth:`describe` about the risk of
        double-counting this against the antenna radiation efficiency.
    antenna
        A :class:`~tss.l2_instrument.antenna.ShortDipole`, or ``None`` to build
        the default 5 m one.  The same antenna is used for transmit and receive;
        the instrument shares one aperture and a separate receive antenna would
        be a different configuration.
    receiver, oscillator
        The receive chain and the frequency reference.
    stated_efficiency_db
        When given, this radiation efficiency is used instead of the one the
        physical antenna model computes.  The baseline configuration states one,
        and a run against the baseline has to be a run against the threshold the
        requirement set writes down rather than against a physics model that may
        be more generous.
    polarization_loss_db, extra_loss_db
        Applied to every echo.
    oversample
        Complex sample rate as a multiple of the largest beat frequency present.
        2.5 is the default: a complex sampler needs only to exceed the highest
        beat frequency, and the margin is for the window skirts.
    max_samples
        Cap on the record length.  A sweep longer than this is generated only up
        to the cap and the truncation is recorded in the flags, together with
        how many echo rows fell outside the generated window.  A silently
        shortened record would read as an ionogram with a missing top end.
    sky_model
        Passed to :func:`tss.l2_instrument.noise.sky_temperature`.
    quantise
        Whether to run the quantiser.  Off is useful for isolating whether an
        effect is quantisation or something else.
    """

    def __init__(
        self,
        *,
        tx_power_w: float = 1.0,
        tx_chain_efficiency_db: float = 0.0,
        antenna: ShortDipole | None = None,
        receiver: Receiver | None = None,
        oscillator: Oscillator | None = None,
        stated_efficiency_db: float | None = None,
        polarization_loss_db: float = DEFAULT_POLARIZATION_LOSS_DB,
        extra_loss_db: float = 0.0,
        oversample: float = 2.5,
        max_samples: int = 1 << 20,
        sky_model: str = "franklin",
        quantise: bool = True,
        isolation: Mapping[str, Any] | None = None,
        resources: Mapping[str, Any] | None = None,
        configuration_name: str = "",
    ) -> None:
        if tx_power_w <= 0.0:
            raise ValueError("tx_power_w must be positive")
        if oversample <= 1.0:
            raise ValueError("oversample must exceed 1 for a complex sampler")
        if max_samples < 16:
            raise ValueError("max_samples must be at least 16")
        self.tx_power_w = float(tx_power_w)
        self.tx_chain_efficiency_db = float(tx_chain_efficiency_db)
        self.antenna = antenna if antenna is not None else ShortDipole()
        self.receiver = receiver if receiver is not None else Receiver()
        self.oscillator = oscillator if oscillator is not None else Oscillator()
        self.stated_efficiency_db = (
            None if stated_efficiency_db is None else float(stated_efficiency_db)
        )
        self.polarization_loss_db = float(polarization_loss_db)
        self.extra_loss_db = float(extra_loss_db)
        self.oversample = float(oversample)
        self.max_samples = int(max_samples)
        self.sky_model = str(sky_model)
        self.quantise = bool(quantise)
        self.isolation = dict(isolation or {})
        self.resources = dict(resources or {})
        self.configuration_name = str(configuration_name)

    # -- configuration -------------------------------------------------------

    @classmethod
    def from_config(cls, cfg: Mapping[str, Any]) -> BaselineInstrument:
        """Build from the instrument configuration shape.

        The shape is the one in ``config/instrument/baseline_6u.yaml``: blocks
        named ``transmit``, ``antenna``, ``receiver``, ``isolation``,
        ``oscillator`` and ``resources``.  Keys that are not recognised are
        ignored rather than rejected, so that adding a key to the configuration
        for another consumer does not break this one.
        """
        transmit = dict(cfg.get("transmit", {}))
        antenna_cfg = dict(cfg.get("antenna", {}))
        antenna = ShortDipole.from_config(antenna_cfg)
        return cls(
            tx_power_w=float(transmit.get("power_w", 1.0)),
            tx_chain_efficiency_db=float(transmit.get("efficiency_db", 0.0)),
            antenna=antenna,
            receiver=Receiver.from_config(dict(cfg.get("receiver", {}))),
            oscillator=Oscillator.from_config(dict(cfg.get("oscillator", {}))),
            stated_efficiency_db=antenna_cfg.get("radiation_efficiency_db"),
            isolation=dict(cfg.get("isolation", {})),
            resources=dict(cfg.get("resources", {})),
            configuration_name=str(cfg.get("name", "")),
        )

    # -- antenna terms -------------------------------------------------------

    def efficiency_db(self, f_mhz: np.ndarray | float) -> np.ndarray:
        """Radiation efficiency used for this run, dB.

        The stated figure when the configuration gives one, otherwise the
        physical model.  Which was used is in :meth:`describe`.
        """
        if self.stated_efficiency_db is not None:
            return np.full(np.shape(np.asarray(f_mhz, dtype=float)), self.stated_efficiency_db)
        return np.asarray(self.antenna.efficiency_db(f_mhz), dtype=float)

    def gain_dbi(self, f_mhz: np.ndarray | float) -> np.ndarray:
        """``G = D + eta``.  Directivity and efficiency added, never lumped."""
        return self.efficiency_db(f_mhz) + self.antenna.directivity_dbi

    def antenna_temperature_k(self, f_mhz: np.ndarray | float) -> np.ndarray:
        """``T_a = eta T_sky + (1-eta) T_phys`` at the given frequencies."""
        eta = undb(self.efficiency_db(f_mhz))
        t_sky = np.asarray(sky_temperature(f_mhz, self.sky_model), dtype=float)
        return antenna_noise_temperature_k(eta, t_sky, self.antenna.physical_temperature)

    # -- the contract --------------------------------------------------------

    def receive(
        self,
        echoes: EchoSet,
        sweep: SweepSchedule,
        rng: np.random.Generator,
    ) -> ReceivedSignal:
        """De-chirped baseband for one sweep.  Inputs are never modified."""
        validate_schedule(sweep)
        freq = np.asarray(echoes.freq_hz, dtype=float)
        delay = np.asarray(echoes.group_delay_s, dtype=float)
        amp_db = np.asarray(echoes.amplitude_db, dtype=float)
        if not (freq.shape == delay.shape == amp_db.shape):
            raise ValueError(
                "EchoSet.freq_hz, .group_delay_s and .amplitude_db must be parallel; "
                f"got {freq.shape}, {delay.shape}, {amp_db.shape}"
            )

        flags_in = (
            np.asarray(echoes.flags, dtype=object)
            if echoes.flags is not None
            else np.array([""] * freq.size, dtype=object)
        )
        usable = (
            np.isfinite(freq)
            & np.isfinite(delay)
            & np.isfinite(amp_db)
            & (delay > 0.0)
            & np.array([str(f) not in ("evanescent", "penetrating") for f in flags_in.ravel()[: freq.size]])
        )

        duration_s = sweep_duration_s(sweep)
        direction = int(sweep.direction)
        k = float(sweep.k_f_hz_per_s)
        f0 = float(sweep.f_start_hz) if direction > 0 else float(sweep.f_stop_hz)

        max_delay = float(np.max(delay[usable])) if np.any(usable) else 0.0
        beat_max = max(k * max_delay, 1.0)
        fs_hz = self.oversample * beat_max
        n_full = int(np.floor(duration_s * fs_hz))
        n_samples = int(min(max(n_full, 16), self.max_samples))
        truncated = n_samples < n_full
        modelled_s = n_samples / fs_hz

        t = np.arange(n_samples, dtype=float) / fs_hz
        f_tx = f0 + direction * k * t

        # Per-row visit time and dwell window.  The spacing is taken from the
        # echo set itself, because that is the frequency step the sounding
        # actually used, and the dwell that follows from it is the dwell the
        # sounding actually got.
        unique_f = np.unique(freq[usable]) if np.any(usable) else np.array([], dtype=float)
        if unique_f.size > 1:
            df = float(np.median(np.diff(unique_f)))
        else:
            df = sweep_bandwidth_hz(sweep)
        dwell_s = max(abs(df) / k, 2.0 / fs_hz)

        samples = np.zeros(n_samples, dtype=complex)
        rows_used = 0
        rows_outside = 0
        for row in np.nonzero(usable)[0]:
            f_k = float(freq[row])
            t_k = direction * (f_k - f0) / k
            if not (0.0 <= t_k <= modelled_s):
                rows_outside += 1
                continue
            lo = int(np.floor(max(t_k - dwell_s / 2.0, 0.0) * fs_hz))
            hi = int(np.ceil(min(t_k + dwell_s / 2.0, modelled_s) * fs_hz))
            if hi <= lo:
                rows_outside += 1
                continue

            tau = float(delay[row])
            f_mhz = f_k / 1.0e6
            # The same antenna transmits and receives, so its gain appears
            # twice.  amp_db already carries the spreading and the absorption.
            gain_dbi = float(self.gain_dbi(f_mhz))
            power_dbw = (
                db(self.tx_power_w)
                + self.tx_chain_efficiency_db
                + gain_dbi
                + float(amp_db[row])
                + self.polarization_loss_db
                + self.extra_loss_db
                + gain_dbi
            )
            amplitude = np.sqrt(undb(power_dbw))

            tt = t[lo:hi]
            phase = 2.0 * np.pi * (
                f0 * tau + direction * k * tau * tt - direction * k * tau**2 / 2.0
            )
            samples[lo:hi] += amplitude * np.exp(1j * phase)
            rows_used += 1

        # Noise.  The antenna temperature follows the sweep, because the sky is
        # three orders of magnitude colder at the top of the band than at the
        # bottom and a single band-averaged figure would be wrong at both ends.
        t_ant = self.antenna_temperature_k(np.clip(f_tx, 1.0, None) / 1.0e6)
        t_sys = t_ant + self.receiver.noise_temperature_k
        sigma = np.sqrt(K_BOLTZ * t_sys * fs_hz / 2.0)
        samples = samples + rng.normal(0.0, sigma) + 1j * rng.normal(0.0, sigma)

        # Receive chain.  Gain, then clipping against the instantaneous dynamic
        # range referred to the median noise floor, then quantisation.
        amplified = self.receiver.apply_gain(samples)
        noise_floor_w = float(np.median(K_BOLTZ * t_sys * fs_hz)) * float(
            undb(self.receiver.gain_db)
        )
        clip_level_w = self.receiver.clip_level_w(noise_floor_w)
        clipped, n_clipped = self.receiver.clip(amplified, clip_level_w)
        out = self.receiver.quantise(clipped, clip_level_w) if self.quantise else clipped

        channels = max(int(self.receiver.channels), 1)
        rows_out = [out]
        # The vector configuration of C5 is three channels.  Until
        # tss.l2_instrument.vector_array exists the extra channels carry the same
        # signal with an independent noise realisation, which is enough to
        # exercise the interface and is not enough to test anything about C5.
        # The flag ``vector_channels_are_placeholders`` says so on every record.
        for _ in range(channels - 1):
            extra_noise = rng.normal(0.0, sigma) + 1j * rng.normal(0.0, sigma)
            chain = self.receiver.apply_gain(samples + extra_noise)
            chain, _ = self.receiver.clip(chain, clip_level_w)
            rows_out.append(
                self.receiver.quantise(chain, clip_level_w) if self.quantise else chain
            )
        block = np.vstack(rows_out)

        floor_f_mhz = unique_f / 1.0e6 if unique_f.size else np.array([1.0])
        noise_floor = noise_density_dbm_per_hz(
            self.antenna_temperature_k(floor_f_mhz) + self.receiver.noise_temperature_k
        )

        flags: dict[str, Any] = {
            # The samples are de-chirped baseband (see the module docstring).
            # L3 reads this and must not mix them against the sweep again.
            "dechirped": True,
            "fs_hz": fs_hz,
            "n_samples": n_samples,
            "sweep_duration_s": duration_s,
            "modelled_duration_s": modelled_s,
            "truncated": bool(truncated),
            "rows_in": int(freq.size),
            "rows_used": rows_used,
            "rows_unusable": int(freq.size - int(np.count_nonzero(usable))),
            "rows_outside_window": rows_outside,
            "dwell_s": dwell_s,
            "max_beat_frequency_hz": beat_max,
            "clip_level_dbm": float(w_to_dbm(clip_level_w)),
            "n_clipped": n_clipped,
            "quantised": bool(self.quantise and self.receiver.n_bits is not None),
            "n_bits": self.receiver.n_bits,
            "sky_model": self.sky_model,
            "noise_floor_f_mhz": floor_f_mhz.tolist(),
            "noise_reference_extrapolated": bool(np.any(is_extrapolated(floor_f_mhz))),
            "leakage_modelled": False,
            "leakage_note": (
                "tss.l2_instrument.leakage is a stub (WP3, critical path).  No "
                "transmit-to-receive leakage is present in these samples, so the "
                "receiver never clips on leakage and no INS-11 or INS-24 verdict "
                "may be taken from this run."
            ),
            "vector_channels_are_placeholders": channels > 1,
            "absorption_db_carried_not_applied": True,
        }

        return ReceivedSignal(
            samples=block,
            fs_hz=fs_hz,
            sweep=sweep,
            noise_floor_dbm_hz=np.asarray(noise_floor, dtype=float),
            flags=flags,
            source=self.source_label(),
        )

    # -- manifest ------------------------------------------------------------

    def source_label(self) -> str:
        name = self.configuration_name or "baseline"
        return (
            f"l2/baseline({name}, Pt={self.tx_power_w:g} W, "
            f"L={self.antenna.length_m:g} m, NF={self.receiver.noise_figure_db:g} dB)"
        )

    def describe(self) -> Mapping[str, Any]:
        """JSON-serialisable identity and settings, for the run manifest."""
        return {
            "layer": "l2",
            "name": getattr(self, "name", type(self).__name__),
            "class": type(self).__name__,
            "module": type(self).__module__,
            "configuration": self.configuration_name,
            "tx_power_w": self.tx_power_w,
            "tx_chain_efficiency_db": self.tx_chain_efficiency_db,
            "polarization_loss_db": self.polarization_loss_db,
            "extra_loss_db": self.extra_loss_db,
            "oversample": self.oversample,
            "max_samples": self.max_samples,
            "sky_model": self.sky_model,
            "quantise": self.quantise,
            "efficiency_source": (
                "configuration radiation_efficiency_db"
                if self.stated_efficiency_db is not None
                else "ShortDipole physical model"
            ),
            "stated_efficiency_db": self.stated_efficiency_db,
            "antenna": dict(self.antenna.describe()),
            "receiver": dict(self.receiver.describe()),
            "oscillator": dict(self.oscillator.describe()),
            "isolation": dict(self.isolation),
            "resources": dict(self.resources),
            "leakage_modelled": False,
            "open_question_transmit_efficiency": (
                "The baseline configuration states transmit.efficiency_db = -20 dB, "
                "described as matching-network dominated, and separately "
                "antenna.radiation_efficiency_db = -20 dB.  A matching network loss "
                "is part of the radiation efficiency in the ShortDipole model, so "
                "applying both may double-count it.  Both are applied here because "
                "that is what the configuration says; the reading has not been "
                "confirmed against the WP3 antenna trade."
            ),
        }
