"""Layer: shared (L2 and L3)
WP: WP3, WP4
Claims: C1, C2
Status: implemented
Requirements: INS-02, INS-05, INS-06, INS-07, INS-27

The stepped sweep law, written once for the transmitter and the processor.

Why a stepped law
-----------------
Requirements Rev G withdraws the uniform 5.6 MHz/s sweep as an ionogram
setting (decision-table row 9, 2026-09-30).  A dwell has to last at least as
long as the round trip to INS-04's 1050 km, 14.0 ms, and at 5.6 MHz/s that
dwell spans 78.4 kHz, which is wider than INS-02's step of ``0.05 f`` below
1.57 MHz.  No single rate satisfies INS-02 at both ends of the band.

The baseline ionogram is therefore the A1 stepped sweep: the band is cut into
steps whose width follows the INS-02 law,

    df_i = min(step_fraction * f_i, step_cap_hz)

with ``f_i`` the lower edge of step ``i``, and each step is swept linearly in
the same dwell ``step_dwell_s``.  The rate within step ``i`` is
``K_i = df_i / step_dwell_s``, so the sweep is continuous in frequency and
piecewise linear in time.  At the A1 settings (0.05, 150 kHz, 26.8 ms) that is
159 steps over 0.33 to 20 MHz, 4.26 s, and a peak rate of 5.6 MHz/s above
3 MHz.  The last step is clipped at ``f_stop_hz`` and keeps the full dwell, so
it sweeps more slowly than the others.

``SweepSchedule.k_f_hz_per_s`` is the peak rate for a stepped schedule.  The
quantities that INS-04 and INS-06 bound are set by the fastest step, so the
peak is the right single number to carry, and :func:`validate_stepped` refuses
a schedule whose ``k_f_hz_per_s`` disagrees with its own step law.

What changes downstream
-----------------------
The range scale of a de-chirped block is ``c f_b / (2 K_i)``, so it now
depends on which step the block came from.  The range cell grows at the bottom
of the band, to ``c / (2 df_i)`` = 9.1 km at 0.33 MHz before windowing, and a
processor that used one rate for the whole sweep would put every echo below
3 MHz at the wrong range.  :func:`rate_hz_per_s` gives the local rate and the
L3 tier reads it per block.

A down-sweep visits the same steps in reverse order, each swept downward, so an
up and down pair is contiguous in frequency as it is for the uniform law.

Nonlinearity
------------
For a stepped sweep ``nonlinearity_ppm`` is applied within each step, as the
same zero-ended parabola the uniform law uses, scaled by that step's width.  A
synthesiser that ramps each step separately shows its departure per ramp, so
this is the stepped counterpart of the uniform placeholder and carries the same
caveat: it is a stand-in with a stated peak, not a measured characteristic.

The law at a step boundary is exact.  In the first ``tau`` of each step the
echo was transmitted during the previous step, so its de-chirped beat ramps
from ``K_(i-1) tau`` to ``K_i tau``.  Below 3 MHz, where ``K_i / K_(i-1) = 1.05``,
the apparent range read with ``K_i`` slides from ``R / 1.05`` to ``R`` across that
interval, which is up to 7.0 ms of the 26.8 ms dwell at INS-04's 1050 km.  Above
3 MHz the rate does not change and there is no transition.  The L2 baseline
builds the de-chirped phase as the exact difference ``phi(t) - phi(t - tau)``
from :func:`phase_rad`, so the transition is in its samples; how L3 treats it
(accept the skirt, gate the first ``tau``, or match-filter it) is a WP4 choice.
"""

from __future__ import annotations

from functools import lru_cache

import numpy as np

from tss.interfaces import SweepSchedule

__all__ = [
    "is_stepped",
    "step_edges_hz",
    "step_count",
    "step_widths_hz",
    "step_rates_hz_per_s",
    "peak_rate_hz_per_s",
    "active_duration_s",
    "validate_stepped",
    "frequency_hz",
    "phase_rad",
    "rate_hz_per_s",
    "time_at_frequency_s",
    "step_index",
    "step_centres_hz",
    "stepped_schedule",
]

_RATE_TOLERANCE = 1.0e-6


def is_stepped(sweep: SweepSchedule) -> bool:
    """True when the schedule uses the stepped law rather than one rate."""
    return float(getattr(sweep, "step_dwell_s", 0.0) or 0.0) > 0.0


@lru_cache(maxsize=64)
def _edges(f_start: float, f_stop: float, fraction: float, cap: float) -> np.ndarray:
    if not (f_stop > f_start > 0.0):
        raise ValueError("a stepped sweep needs 0 < f_start_hz < f_stop_hz")
    if fraction <= 0.0 and cap <= 0.0:
        raise ValueError("a stepped sweep needs step_fraction or step_cap_hz")
    edges = [f_start]
    f = f_start
    limit = 1_000_000
    while f < f_stop:
        widths = []
        if fraction > 0.0:
            widths.append(fraction * f)
        if cap > 0.0:
            widths.append(cap)
        f = min(f + min(widths), f_stop)
        edges.append(f)
        if len(edges) > limit:
            raise ValueError("stepped sweep has more than a million steps")
    out = np.asarray(edges, dtype=float)
    out.setflags(write=False)
    return out


def step_edges_hz(sweep: SweepSchedule) -> np.ndarray:
    """Ascending step edges, ``n_steps + 1`` of them, from f_start to f_stop."""
    return _edges(
        float(sweep.f_start_hz),
        float(sweep.f_stop_hz),
        float(getattr(sweep, "step_fraction", 0.0) or 0.0),
        float(getattr(sweep, "step_cap_hz", 0.0) or 0.0),
    )


def step_count(sweep: SweepSchedule) -> int:
    """Number of steps."""
    return int(step_edges_hz(sweep).size - 1)


def step_widths_hz(sweep: SweepSchedule) -> np.ndarray:
    """Width of each step in ascending frequency order."""
    return np.diff(step_edges_hz(sweep))


def step_centres_hz(sweep: SweepSchedule) -> np.ndarray:
    """Centre frequency of each step in ascending order.

    These are the natural sounding frequencies of a stepped ionogram: one
    column per step, filed at the frequency the step was sweeping through at
    the middle of its dwell.
    """
    e = step_edges_hz(sweep)
    return 0.5 * (e[:-1] + e[1:])


def step_rates_hz_per_s(sweep: SweepSchedule) -> np.ndarray:
    """Rate magnitude of each step in ascending frequency order."""
    return step_widths_hz(sweep) / float(sweep.step_dwell_s)


def peak_rate_hz_per_s(sweep: SweepSchedule) -> float:
    """Fastest step rate, which is what ``k_f_hz_per_s`` carries."""
    return float(np.max(step_rates_hz_per_s(sweep)))


def active_duration_s(sweep: SweepSchedule) -> float:
    """Time to cross the band, ``n_steps * step_dwell_s``."""
    return step_count(sweep) * float(sweep.step_dwell_s)


def validate_stepped(sweep: SweepSchedule) -> None:
    """Raise if a stepped schedule contradicts itself.

    ``k_f_hz_per_s`` must equal the peak step rate.  Carrying a different
    number there would make INS-04 and INS-06 score against a rate the
    instrument never transmits.
    """
    if float(sweep.step_dwell_s) <= 0.0:
        raise ValueError("step_dwell_s must be positive for a stepped sweep")
    if float(getattr(sweep, "step_fraction", 0.0)) < 0.0 or float(
        getattr(sweep, "step_cap_hz", 0.0)
    ) < 0.0:
        raise ValueError("step_fraction and step_cap_hz must not be negative")
    peak = peak_rate_hz_per_s(sweep)
    k = abs(float(sweep.k_f_hz_per_s))
    if abs(k - peak) > _RATE_TOLERANCE * peak:
        raise ValueError(
            f"k_f_hz_per_s is {k:.6g} Hz/s but the step law peaks at {peak:.6g} "
            "Hz/s; for a stepped sweep k_f_hz_per_s is the peak step rate"
        )


def _order(sweep: SweepSchedule) -> tuple[np.ndarray, np.ndarray, int]:
    """Step start frequencies and widths in transmission order, and direction."""
    e = step_edges_hz(sweep)
    widths = np.diff(e)
    direction = 1 if int(sweep.direction) >= 0 else -1
    if direction > 0:
        return e[:-1], widths, 1
    return e[1:][::-1], widths[::-1], -1


@lru_cache(maxsize=64)
def _phase_prefix(
    f_start: float, f_stop: float, fraction: float, cap: float,
    dwell: float, direction: int, eps: float,
) -> np.ndarray:
    """Accumulated phase in cycles at the start of each step."""
    e = _edges(f_start, f_stop, fraction, cap)
    widths = np.diff(e)
    if direction > 0:
        starts = e[:-1]
    else:
        starts, widths = e[1:][::-1], widths[::-1]
    k = widths / dwell
    # Full-step integral: f0 T + d K T^2/2 + d B eps * (2T/3), since
    # int_0^T 4 s(1-s) dt = 4T (1/2 - 1/3) = 2T/3.
    per_step = starts * dwell + direction * k * dwell**2 / 2.0 + direction * widths * eps * (
        2.0 * dwell / 3.0
    )
    out = np.concatenate([[0.0], np.cumsum(per_step)])
    out.setflags(write=False)
    return out


def step_index(sweep: SweepSchedule, t_s: np.ndarray | float) -> np.ndarray:
    """Transmission-order step index for each time, clipped to the sweep."""
    t = np.asarray(t_s, dtype=float)
    n = step_count(sweep)
    idx = np.floor(t / float(sweep.step_dwell_s)).astype(int)
    return np.clip(idx, 0, n - 1)


def _local(sweep: SweepSchedule, t_s):
    t = np.asarray(t_s, dtype=float)
    starts, widths, direction = _order(sweep)
    dwell = float(sweep.step_dwell_s)
    i = step_index(sweep, t)
    u = t - i * dwell
    return t, i, u, starts[i], widths[i], widths[i] / dwell, direction, dwell


def frequency_hz(sweep: SweepSchedule, t_s: np.ndarray | float) -> np.ndarray:
    """Instantaneous transmitted frequency, nonlinearity included."""
    _, _, u, f0, width, k, d, dwell = _local(sweep, t_s)
    f = f0 + d * k * u
    eps = float(sweep.nonlinearity_ppm) * 1.0e-6
    if eps != 0.0:
        s = np.clip(u / dwell, 0.0, 1.0)
        f = f + d * width * eps * 4.0 * s * (1.0 - s)
    return f


def phase_rad(sweep: SweepSchedule, t_s: np.ndarray | float) -> np.ndarray:
    """Accumulated transmitted phase, radians, as an exact integral."""
    _, i, u, f0, width, k, d, dwell = _local(sweep, t_s)
    eps = float(sweep.nonlinearity_ppm) * 1.0e-6
    prefix = _phase_prefix(
        float(sweep.f_start_hz), float(sweep.f_stop_hz),
        float(getattr(sweep, "step_fraction", 0.0) or 0.0),
        float(getattr(sweep, "step_cap_hz", 0.0) or 0.0),
        dwell, d, eps,
    )
    cycles = prefix[i] + f0 * u + d * k * u**2 / 2.0
    if eps != 0.0:
        s = np.clip(u / dwell, 0.0, 1.0)
        cycles = cycles + d * width * eps * 4.0 * dwell * (s**2 / 2.0 - s**3 / 3.0)
    return 2.0 * np.pi * cycles


def rate_hz_per_s(sweep: SweepSchedule, t_s: np.ndarray | float) -> np.ndarray:
    """Local sweep rate magnitude at each time, ``K_i`` of the step in use."""
    _, _, _, _, _, k, _, _ = _local(sweep, t_s)
    return k


def time_at_frequency_s(sweep: SweepSchedule, f_hz: np.ndarray | float) -> np.ndarray:
    """Time from the start of the sweep at which ``f_hz`` is transmitted.

    Ignores the nonlinearity, as the uniform law's visit time does.  Frequencies
    outside the band are mapped by extending the end steps.
    """
    f = np.asarray(f_hz, dtype=float)
    e = step_edges_hz(sweep)
    dwell = float(sweep.step_dwell_s)
    n = e.size - 1
    j = np.clip(np.searchsorted(e, f, side="right") - 1, 0, n - 1)
    width = e[j + 1] - e[j]
    frac = (f - e[j]) / width
    if int(sweep.direction) >= 0:
        return (j + frac) * dwell
    return (n - j - frac) * dwell


def stepped_schedule(
    f_start_hz: float,
    f_stop_hz: float,
    step_dwell_s: float,
    *,
    step_fraction: float = 0.05,
    step_cap_hz: float = 150.0e3,
    sweep_period_s: float = 0.0,
    direction: int = 1,
    nonlinearity_ppm: float = 0.0,
) -> SweepSchedule:
    """Build a stepped schedule with ``k_f_hz_per_s`` set to its peak rate.

    ``sweep_period_s`` of zero means the active duration, with no dead time.
    The defaults are the A1 threshold law of INS-02.
    """
    probe = SweepSchedule(
        f_start_hz=float(f_start_hz), f_stop_hz=float(f_stop_hz), k_f_hz_per_s=1.0,
        sweep_period_s=1.0, direction=int(direction),
        nonlinearity_ppm=float(nonlinearity_ppm), step_dwell_s=float(step_dwell_s),
        step_fraction=float(step_fraction), step_cap_hz=float(step_cap_hz),
    )
    peak = peak_rate_hz_per_s(probe)
    duration = active_duration_s(probe)
    period = float(sweep_period_s) if float(sweep_period_s) > 0.0 else duration
    return SweepSchedule(
        f_start_hz=float(f_start_hz), f_stop_hz=float(f_stop_hz), k_f_hz_per_s=peak,
        sweep_period_s=period, direction=int(direction),
        nonlinearity_ppm=float(nonlinearity_ppm), step_dwell_s=float(step_dwell_s),
        step_fraction=float(step_fraction), step_cap_hz=float(step_cap_hz),
    )
