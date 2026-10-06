"""Layer: L3
WP: WP4
Claims: C1, C4
Status: implemented
Requirements: INS-21, INS-24, INS-26

De-chirp the received block against the transmitted sweep, and the range
correlation that decides how pure the sweep has to be.

De-chirping
-----------
The instrument transmits ``exp(i phi(t))`` and receives a delayed, Doppler
shifted copy.  Multiplying the received samples by the conjugate of the
reference collapses the sweep and leaves a tone whose frequency is the beat
frequency.  The reference is the sweep that was actually transmitted, carried on
:attr:`tss.interfaces.ReceivedSignal.sweep`, not the one the campaign asked for,
which is the distinction that makes a nonlinearity term meaningful at all.

The mixer product ``echo * conj(reference)`` has instantaneous frequency
``f_rx - f_tx = -K_f tau + f_d`` on an up-sweep, so the echo lands on the
negative frequency axis; on a down-sweep it lands on the positive one.  Carrying
that asymmetry into the transform is a reliable source of sign errors, so
:func:`dechirp` orients the result by conjugating on an up-sweep.  After
orientation the echo is always at positive beat frequency ``K_f tau - direction
* f_d``, which is exactly the apparent range of
:func:`tss.l1_propagation.doppler.apparent_range_bias_m`.  The orientation can
be switched off for anyone who wants the raw mixer product.

Range correlation, and what it does to the waveform purity requirement
----------------------------------------------------------------------
This is the part that changes a requirement, so it is written out.

INS-21 mandates a single reference oscillator shared by the transmitter and the
de-chirp mixer.  The consequence is that oscillator noise is common to both
arms of the mixer, and at the output it survives only to the extent that the
two arms are decorrelated by the round-trip delay ``tau``.  For a sinusoidal
frequency modulation of peak deviation ``df`` at rate ``f_m``, the de-chirped
phase residual has peak amplitude

    beta = (2 * df / f_m) * sin(pi * f_m * tau)

which for ``f_m * tau << 1`` tends to ``2 * pi * df * tau``.  In power, the
sideband at offset ``f_m`` is suppressed by

    S(f_m) = 4 * sin^2(pi * f_m * tau)

relative to what a non-coherent receiver would see.  The same factor applies to
a continuous phase noise spectrum: the de-chirped residual spectral density is
``L(f_m) + 10 log10(4 sin^2(pi f_m tau))``.

Put numbers on it.  A topside echo at 300 km has ``tau = 2 ms``.  Close-in phase
noise at 100 Hz offset is suppressed by ``4 sin^2(pi * 100 * 0.002)`` which is
0.38, about 4 dB.  At 10 Hz offset the suppression is ``4 sin^2(0.0628)``, about
0.0157, or 18 dB of rejection.  At 1 Hz it is 38 dB.  Close-in noise, which is
the expensive part of an oscillator specification, is exactly the part that
range correlation removes.

The consequence for the design is direct.  The de-chirped noise floor demanded
by INS-10 does not translate into a demand for a laboratory-grade reference; it
translates into a much weaker demand that commercial temperature-compensated
parts already meet, with the INS-26 Allan deviation figure covering the slow
drift that range correlation does **not** remove.  That is the honest reading of
C4.  The claim is not that the waveform is pure enough.  The claim is about
transmit-to-receive isolation and receiver dynamic range, because the leakage
path has ``tau`` near zero, gets no range correlation benefit at all, and is
tens of decibels stronger than the echo.  The hard problem moved; it did not go
away.  :mod:`tss.l3_signal.novel.leakage_cancel` is where it sits now.

One caveat on the suppression factor.  It is derived for a two-way delay that
is constant across the coherent interval.  Over the 1.15 s goal interval of
INS-15, the delay changes, and the suppression is then an average over the
delays visited rather than the value at one of them.  For a topside sounding the
delay changes by a few per cent at most, so the correction is small, but it is a
correction and not zero.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

import numpy as np

from tss.constants import C_LIGHT
from tss.interfaces import ReceivedSignal, SweepSchedule
from tss.l3_signal.base import (
    instantaneous_frequency_hz,
    local_rate_hz_per_s,
    sweep_phase_rad,
)

__all__ = [
    "Dechirped",
    "reference_waveform",
    "dechirp",
    "two_way_delay_s",
    "range_correlation_suppression",
    "range_correlation_suppression_db",
    "phase_residual_beta_rad",
    "residual_phase_noise_dbc_hz",
]


@dataclass(frozen=True)
class Dechirped:
    """A de-chirped block, ready for the transform.

    Attributes
    ----------
    samples
        Complex samples, ``(n_channels, n_samples)``.
    fs_hz
        Sample rate, unchanged by de-chirping.
    t_s
        Sample times measured from the start of the sweep this block belongs
        to.  Carried so that the transform can label the block with the
        instantaneous sweep frequency without re-deriving the timebase.
    sweep
        The sweep the block was de-chirped against.
    oriented
        True when the direction-dependent conjugation was applied, so that the
        echo sits at positive beat frequency whichever way the sweep ran.
    centre_frequency_hz
        Instantaneous transmitted frequency at the centre of the block.  This
        is the frequency the resulting ionogram column is filed under.
    k_f_hz_per_s
        Sweep rate magnitude at the centre of the block, the rate that converts
        this block's beat frequencies to range.  For a uniform sweep it is the
        schedule's rate.  For a stepped sweep it is the rate of the step the
        block sits in.  ``None`` means use the schedule's rate, which keeps
        blocks built by hand before this field existed working unchanged.
    rate_uniform
        False when the block spans samples from steps with different rates, in
        which case a single range scale is only approximate for it.
    """

    samples: np.ndarray
    fs_hz: float
    t_s: np.ndarray
    sweep: SweepSchedule
    oriented: bool
    centre_frequency_hz: float
    k_f_hz_per_s: float | None = None
    rate_uniform: bool = True

    def describe(self) -> Mapping[str, Any]:
        """Small record of the block, for diagnostics and manifests."""
        return {
            "n_channels": int(self.samples.shape[0]),
            "n_samples": int(self.samples.shape[-1]),
            "fs_hz": float(self.fs_hz),
            "oriented": bool(self.oriented),
            "centre_frequency_hz": float(self.centre_frequency_hz),
            "k_f_hz_per_s": None if self.k_f_hz_per_s is None else float(self.k_f_hz_per_s),
            "rate_uniform": bool(self.rate_uniform),
        }


def two_way_delay_s(range_m: np.ndarray | float) -> np.ndarray:
    """Round-trip delay for a one-way range.  ``tau = 2 R / c``."""
    return 2.0 * np.asarray(range_m, dtype=float) / C_LIGHT


def reference_waveform(
    sweep: SweepSchedule,
    t_s: np.ndarray,
) -> np.ndarray:
    """The de-chirp reference, ``exp(i phi(t))``, at the given sample times.

    Unit amplitude.  Under INS-21 this is generated from the same reference
    oscillator as the transmitted sweep, which is what makes the range
    correlation of :func:`range_correlation_suppression` legitimate.  A
    simulator that wants to model two independent oscillators has to break that
    assumption deliberately, and it should not do it by quietly passing a
    different ``SweepSchedule`` here.
    """
    return np.exp(1j * sweep_phase_rad(sweep, np.asarray(t_s, dtype=float)))


def dechirp(
    signal: ReceivedSignal,
    *,
    start: int = 0,
    stop: int | None = None,
    t0_s: float = 0.0,
    orient: bool = True,
) -> Dechirped:
    """De-chirp a slice of a received record against its own sweep.

    A record flagged ``flags["dechirped"]`` (see
    :class:`tss.interfaces.ReceivedSignal`) is already at baseband and is
    passed through, oriented, without a second mix.

    Parameters
    ----------
    signal
        The received record.  ``signal.sweep`` is the reference, per INS-21.
    start, stop
        Sample indices bounding the block, as for a slice.  The default takes
        the whole record.
    t0_s
        Time of sample index 0 measured from the start of the sweep.  A block
        taken from the middle of a sweep must say where it starts, or the
        reference is the wrong reference and the beat frequency is wrong by the
        offset.
    orient
        When true, conjugate the product on an up-sweep so that the echo sits
        at positive beat frequency either way.  See the module docstring.

    Returns
    -------
    Dechirped
        A new frozen object.  The input is not modified.
    """
    samples = np.atleast_2d(np.asarray(signal.samples))
    if samples.ndim != 2:
        raise ValueError(
            f"ReceivedSignal.samples must be (n_channels, n_samples); got {samples.shape}"
        )
    stop_idx = samples.shape[1] if stop is None else int(stop)
    block = samples[:, int(start):stop_idx]
    if block.shape[1] == 0:
        raise ValueError("de-chirp block is empty")

    fs = float(signal.fs_hz)
    if fs <= 0.0:
        raise ValueError("ReceivedSignal.fs_hz must be positive")

    t = t0_s + (np.arange(int(start), stop_idx, dtype=float) / fs)
    direction = 1 if int(signal.sweep.direction) >= 0 else -1

    if bool(signal.flags.get("dechirped", False)):
        # Already de-chirped in hardware (INS-21): an echo is a tone at
        # ``direction * K_f * tau``.  Mixing it against the sweep again would
        # turn every tone back into a chirp across the block and smear it over
        # every range bin, which is the fault this branch exists to prevent.
        # Orientation then only has to flip the down-sweep.
        product = np.asarray(block, dtype=complex).copy()
        if orient and direction < 0:
            product = np.conj(product)
    else:
        ref = reference_waveform(signal.sweep, t)
        product = np.asarray(block, dtype=complex) * np.conj(ref)[None, :]
        if orient and direction > 0:
            product = np.conj(product)

    t_mid = float(np.mean(t))
    centre = float(np.atleast_1d(instantaneous_frequency_hz(signal.sweep, t_mid))[0])
    rates = np.atleast_1d(local_rate_hz_per_s(signal.sweep, np.array([t[0], t_mid, t[-1]])))
    k_centre = float(rates[1])
    uniform = bool(np.all(np.abs(rates - k_centre) <= 1e-9 * max(k_centre, 1.0)))

    return Dechirped(
        samples=product,
        fs_hz=fs,
        t_s=t,
        sweep=signal.sweep,
        oriented=bool(orient),
        centre_frequency_hz=centre,
        k_f_hz_per_s=k_centre,
        rate_uniform=uniform,
    )


def range_correlation_suppression(
    f_m_hz: np.ndarray | float,
    tau_s: np.ndarray | float,
) -> np.ndarray:
    """Range-correlation power suppression factor, ``4 sin^2(pi f_m tau)``.

    A linear power ratio, not decibels, and not bounded by one: the factor
    reaches 4, which is 6 dB of *enhancement*, when ``f_m tau`` is a half
    integer.  That is real and it is why the factor is returned as computed
    rather than clipped.  For the offsets that matter in a topside sounding,
    ``f_m tau`` is well below a half and the factor is a suppression.

    Small-argument limit: ``4 (pi f_m tau)^2``, so the suppression improves as
    the square of the offset frequency as the offset falls.  Close-in noise is
    rejected hardest.

    Parameters
    ----------
    f_m_hz
        Offset from the carrier, Hz.
    tau_s
        Two-way delay, seconds.
    """
    f_m = np.asarray(f_m_hz, dtype=float)
    tau = np.asarray(tau_s, dtype=float)
    return 4.0 * np.sin(np.pi * f_m * tau) ** 2


def range_correlation_suppression_db(
    f_m_hz: np.ndarray | float,
    tau_s: np.ndarray | float,
) -> np.ndarray:
    """:func:`range_correlation_suppression` in decibels.

    Negative where the factor suppresses.  Returns ``-inf`` at the exact zeros
    of the sine, which is correct and is also an artefact of assuming a
    perfectly constant delay; see the caveat in the module docstring.
    """
    from tss.units import db

    return db(range_correlation_suppression(f_m_hz, tau_s))


def phase_residual_beta_rad(
    df_hz: np.ndarray | float,
    f_m_hz: np.ndarray | float,
    tau_s: np.ndarray | float,
) -> np.ndarray:
    """Peak de-chirped phase residual from a sinusoidal frequency modulation.

        ``beta = (2 df / f_m) * sin(pi f_m tau)``

    Evaluated as ``2 pi df tau * sinc(f_m tau)`` with the numpy normalisation of
    :func:`numpy.sinc`, which is the same expression rearranged and is finite at
    ``f_m = 0`` instead of dividing by zero there.  The limit it takes is
    ``2 pi df tau``, which is the statement that a very slow frequency error
    turns into a pure range error and range correlation does nothing about it.
    That slow term is what the INS-26 Allan deviation requirement covers.

    Parameters
    ----------
    df_hz
        Peak frequency deviation of the modulation, Hz.
    f_m_hz
        Modulation rate, Hz.
    tau_s
        Two-way delay, seconds.

    Returns
    -------
    numpy.ndarray
        Peak phase residual in radians.
    """
    df = np.asarray(df_hz, dtype=float)
    f_m = np.asarray(f_m_hz, dtype=float)
    tau = np.asarray(tau_s, dtype=float)
    return 2.0 * np.pi * df * tau * np.sinc(f_m * tau)


def residual_phase_noise_dbc_hz(
    phase_noise_dbc_hz: np.ndarray | float,
    f_m_hz: np.ndarray | float,
    tau_s: np.ndarray | float,
) -> np.ndarray:
    """Oscillator phase noise after range correlation, dBc/Hz.

    ``L(f_m) + 10 log10(4 sin^2(pi f_m tau))``.  This is the spectrum that
    actually competes with an echo at delay ``tau`` in a de-chirped FMCW
    receiver with a shared reference, and it is the number that should be
    compared against the INS-10 signal to noise budget.  Comparing the raw
    oscillator specification against that budget overstates the requirement by
    the suppression, which for the close-in region is tens of decibels.

    The leakage path is the case where this gives no relief: its delay is
    essentially zero, the suppression factor is essentially zero as well, and
    the leakage carries the oscillator's full noise into the receiver.  Scoring
    INS-24 against this function with ``tau_s = 0`` returns ``-inf`` dBc/Hz,
    which is not a useful answer and is a sign the wrong question was asked.
    """
    from tss.units import db

    return np.asarray(phase_noise_dbc_hz, dtype=float) + db(
        range_correlation_suppression(f_m_hz, tau_s)
    )
