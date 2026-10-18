"""Layer: L3
WP: WP4
Claims: C1, C4
Status: implemented
Requirements: -

Shared helpers for the L3 signal processing tiers.

Every L3 tier works from the same three things: the transmitted sweep law, the
sample timebase, and the mapping between a de-chirped beat frequency and an
apparent range.  Those live here so that the baseline tier of
:mod:`tss.l3_signal.ionogram` and the novel tiers under
:mod:`tss.l3_signal.novel` cannot disagree about them.

The sweep law
-------------
:class:`~tss.interfaces.SweepSchedule` describes a linear sweep from
``f_start_hz`` to ``f_stop_hz`` at rate ``k_f_hz_per_s``, repeated every
``sweep_period_s``, with ``direction`` +1 for an up-sweep and -1 for a
down-sweep.  The instantaneous frequency at time ``t`` measured from the start
of a sweep is

    f(t) = f_start + direction * K_f * t          (up-sweep)
    f(t) = f_stop  - |direction| * K_f * t        (down-sweep)

which this module writes once, in :func:`instantaneous_frequency_hz`, because
getting the down-sweep convention wrong is the single easiest way to get the
sign of the Doppler bias wrong, and the sign of the Doppler bias is C3.

``k_f_hz_per_s`` is carried as a magnitude and the sign lives in ``direction``.
That is the convention :mod:`tss.l1_propagation.doppler` uses too.

Nonlinearity
------------
``SweepSchedule.nonlinearity_ppm`` is the peak departure of the real sweep from
the ideal linear law, in parts per million of the swept bandwidth.  It is
modelled here as a single quadratic term, which is the shape a directly
synthesised sweep actually shows and is not a general model.  A real
synthesiser's departure has structure at many scales and measuring it is a WP9
bench task.  Until then, a quadratic is an honest placeholder as long as nobody
reads a result from it as a bound.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import numpy as np

from tss import sweep_law
from tss.constants import C_LIGHT, beat_frequency
from tss.interfaces import SweepSchedule

__all__ = [
    "SignalProcessingBase",
    "sweep_duration_s",
    "sweep_bandwidth_hz",
    "instantaneous_frequency_hz",
    "sweep_phase_rad",
    "range_from_beat_hz",
    "beat_from_range_m",
    "local_rate_hz_per_s",
]


def sweep_bandwidth_hz(sweep: SweepSchedule) -> float:
    """Swept bandwidth, always positive."""
    return abs(float(sweep.f_stop_hz) - float(sweep.f_start_hz))


def sweep_duration_s(sweep: SweepSchedule) -> float:
    """Time to cover the band at the stated rate, seconds.

    This is the active sweep time, which is at most ``sweep_period_s``.  A
    schedule whose period is shorter than its sweep duration is inconsistent and
    this raises, because the alternative is a silently truncated sweep and a
    silently wrong range scale.
    """
    k_f = abs(float(sweep.k_f_hz_per_s))
    if k_f <= 0.0:
        raise ValueError("SweepSchedule.k_f_hz_per_s must be non-zero")
    if sweep_law.is_stepped(sweep):
        duration = sweep_law.active_duration_s(sweep)
    else:
        duration = sweep_bandwidth_hz(sweep) / k_f
    period = float(sweep.sweep_period_s)
    if period > 0.0 and duration > period * (1.0 + 1e-9):
        raise ValueError(
            f"sweep needs {duration:.6g} s at {k_f:.6g} Hz/s but the period is "
            f"{period:.6g} s; the schedule is inconsistent"
        )
    return duration


def instantaneous_frequency_hz(
    sweep: SweepSchedule,
    t_s: np.ndarray | float,
) -> np.ndarray:
    """Transmitted frequency at time ``t_s`` from the start of a sweep, Hz.

    Includes the quadratic nonlinearity term when
    ``SweepSchedule.nonlinearity_ppm`` is non-zero.  The term is written so that
    it vanishes at both ends of the sweep and reaches its stated peak departure
    in the middle, which is the parabola a two-point calibration of a real
    synthesiser leaves behind.
    """
    t = np.asarray(t_s, dtype=float)
    if sweep_law.is_stepped(sweep):
        return sweep_law.frequency_hz(sweep, t)
    k_f = abs(float(sweep.k_f_hz_per_s))
    direction = 1 if int(sweep.direction) >= 0 else -1
    f0 = float(sweep.f_start_hz) if direction > 0 else float(sweep.f_stop_hz)
    f_ideal = f0 + direction * k_f * t

    ppm = float(sweep.nonlinearity_ppm)
    if ppm == 0.0:
        return f_ideal

    duration = sweep_duration_s(sweep)
    if duration <= 0.0:
        return f_ideal
    u = np.clip(t / duration, 0.0, 1.0)
    peak_hz = ppm * 1e-6 * sweep_bandwidth_hz(sweep)
    return f_ideal + direction * peak_hz * 4.0 * u * (1.0 - u)


def sweep_phase_rad(
    sweep: SweepSchedule,
    t_s: np.ndarray | float,
) -> np.ndarray:
    """Phase of the transmitted sweep at time ``t_s``, radians.

    The integral of the instantaneous frequency.  For the ideal linear law this
    is the closed form ``2 pi (f0 t + direction K_f t^2 / 2)``.  When a
    nonlinearity is present the integral is taken numerically along whatever
    timebase the caller passed, using the cumulative trapezoidal rule, so a
    caller that passes a sparse timebase gets a phase that is only as accurate
    as that timebase.  That is stated rather than guarded against, because the
    sample timebase is what the receiver actually has.
    """
    t = np.asarray(t_s, dtype=float)
    if sweep_law.is_stepped(sweep):
        # Exact per-step integral, so the processor's reference and the
        # transmitter's waveform carry the same phase.
        return sweep_law.phase_rad(sweep, t)
    k_f = abs(float(sweep.k_f_hz_per_s))
    direction = 1 if int(sweep.direction) >= 0 else -1
    f0 = float(sweep.f_start_hz) if direction > 0 else float(sweep.f_stop_hz)

    if float(sweep.nonlinearity_ppm) == 0.0:
        return 2.0 * np.pi * (f0 * t + direction * k_f * t**2 / 2.0)

    f_inst = instantaneous_frequency_hz(sweep, t)
    phase = np.zeros_like(t)
    if t.size > 1:
        dt = np.diff(t)
        mid = 0.5 * (f_inst[1:] + f_inst[:-1])
        phase[1:] = np.cumsum(mid * dt)
    return 2.0 * np.pi * phase


def beat_from_range_m(
    range_m: np.ndarray | float,
    k_f_hz_per_s: float,
) -> np.ndarray:
    """De-chirped beat frequency for a two-way range.  ``f_b = 2 K_f R / c``.

    A thin wrapper on :func:`tss.constants.beat_frequency`, kept so that L3 code
    reads symmetrically with :func:`range_from_beat_hz` and so that there is one
    obvious place to find the pair.
    """
    return beat_frequency(range_m, abs(float(k_f_hz_per_s)))


def range_from_beat_hz(
    beat_hz: np.ndarray | float,
    k_f_hz_per_s: float,
) -> np.ndarray:
    """Apparent range from a de-chirped beat frequency.  ``R = c f_b / (2 K_f)``.

    The exact inverse of :func:`beat_from_range_m`, and the round trip is
    checked in ``tests/unit/test_l3_signal.py``.  "Apparent" is not a hedge: a
    Doppler shift, a sweep rate error and a sweep nonlinearity all add to the
    beat frequency and this function converts the sum, so the number it returns
    is the range the processor believes in and not the range the echo came from.
    """
    f_b = np.asarray(beat_hz, dtype=float)
    return C_LIGHT * f_b / (2.0 * abs(float(k_f_hz_per_s)))


def local_rate_hz_per_s(sweep: SweepSchedule, t_s: np.ndarray | float) -> np.ndarray:
    """Sweep rate magnitude in use at time ``t_s``.

    ``k_f_hz_per_s`` for a uniform sweep.  For a stepped sweep, the rate of the
    step being transmitted, which is what converts a beat frequency in a block
    from that step into range.
    """
    t = np.asarray(t_s, dtype=float)
    if sweep_law.is_stepped(sweep):
        return sweep_law.rate_hz_per_s(sweep, t)
    return np.full_like(t, abs(float(sweep.k_f_hz_per_s)), dtype=float)


class SignalProcessingBase:
    """Common base for L3 tiers.

    Holds the description machinery and the sweep helpers.  A subclass
    registers itself with :func:`tss.registry.register` and implements
    ``form_ionogram`` to satisfy
    :class:`tss.interfaces.L3SignalProcessing`.
    """

    name: str = "base"
    layer: str = "l3"

    def __init__(self, *, window: str = "hann") -> None:
        self.window = str(window)

    # -- sweep helpers ------------------------------------------------------

    @staticmethod
    def instantaneous_frequency_hz(
        sweep: SweepSchedule, t_s: np.ndarray | float
    ) -> np.ndarray:
        return instantaneous_frequency_hz(sweep, t_s)

    @staticmethod
    def range_from_beat_hz(beat_hz: np.ndarray | float, k_f: float) -> np.ndarray:
        return range_from_beat_hz(beat_hz, k_f)

    @staticmethod
    def beat_from_range_m(range_m: np.ndarray | float, k_f: float) -> np.ndarray:
        return beat_from_range_m(range_m, k_f)

    # -- manifest -----------------------------------------------------------

    def _settings(self) -> dict[str, Any]:
        return {"window": self.window}

    def describe(self) -> Mapping[str, Any]:
        """JSON-serialisable identity and settings, for the run manifest."""
        record: dict[str, Any] = {
            "layer": "l3",
            "name": getattr(self, "name", type(self).__name__),
            "class": type(self).__name__,
            "module": type(self).__module__,
        }
        record.update(self._settings())
        return record
