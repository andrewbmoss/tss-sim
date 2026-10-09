"""Layer: L3
WP: WP5
Claims: C3
Status: implemented (parameter level)
Requirements: SCI-15, SCI-23

The route-two estimator of Technical Volume Section 8.11: drift and tilt from
two nodes that look at the same structure a known interval apart.

One pass measures, for the echo it is tracking, an apparent range ``P'`` and a
line-of-sight rate ``r`` (the phase-path rate, ``c f_D / 2 f``, positive when
closing).  For a reflecting surface that is locally a plane, with vertical
velocity ``w`` and upward unit normal whose horizontal part, as seen at the
spacecraft, is ``s``,

    r / mu_sat = w cos(theta) - v . s                                   (1)

where ``v`` is the platform's velocity relative to the medium and ``theta`` is
the arrival angle at the spacecraft (``sin theta = |s|``).  One pass cannot
divide ``w`` from ``s``.  A second pass over the same point ``dt`` later sees
the surface ``w dt`` nearer, so

    dP' = -w dt cos(theta) / mu                                         (2)

with ``1 / mu`` the mean group index at the two spacecraft.  Equation (2) gives
``w``, and equation (1) then gives the component of ``s`` along each pass's
velocity.  Two passes along the same ground track give the along-track
component.  Two passes along crossing tracks give both components, and so the
full tilt.

What is assumed, and tested elsewhere
-------------------------------------
That the surface is a plane over the footprint and holds still, apart from its
drift, over ``dt``.  The two-node study runs this estimator against a forward
model that assumes neither (:mod:`tss.l1_propagation.surface_ranging` in three
dimensions, :mod:`tss.l1_propagation.plate_echo` in two with the plasma), and
measures what evolution between the passes costs.

The compact facet
-----------------
For a compact facet each pass takes ``P'`` and ``r`` at the peak of the echo
power, where the spacecraft is on the facet's normal.  Route one already
returns ``r / mu_sat`` there, as its drift.  Equations (1) and (2) then hold as
written, and the shift of the power peak between the passes gives the facet's
horizontal drift, which is what SCI-23 names.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from tss.constants import C_LIGHT

__all__ = [
    "PassMeasurement",
    "SecondLookEstimate",
    "group_path_from_two_frequencies",
    "los_rate_from_field",
    "estimate_second_look",
    "horizontal_drift_from_peaks",
]


@dataclass(frozen=True)
class PassMeasurement:
    """What one node measured on one pass.

    Attributes
    ----------
    t_s
        Time of the measurement: the moment the node passed the reference
        point, or the power peak of a compact facet.
    group_path_m
        Apparent range ``P'`` of the echo at ``t_s``.
    los_rate_mps
        Phase-path rate at ``t_s``, positive when closing.
    velocity_mps
        Platform velocity relative to the medium, local ``(x, y)`` components.
    mu_sat
        Refractive index at the spacecraft, from its own measured local plasma
        frequency.  1 in free space.
    """

    t_s: float
    group_path_m: float
    los_rate_mps: float
    velocity_mps: tuple[float, float]
    mu_sat: float = 1.0


@dataclass(frozen=True)
class SecondLookEstimate:
    """Drift and tilt from two passes.

    ``slope_x`` and ``slope_y`` are the components of ``s``, the horizontal part
    of the unit normal at the spacecraft; ``asin(|s|)`` is the arrival angle.
    For a same-track pair only the along-track component is observable and the
    other is NaN.  ``slope_along_1`` and ``slope_along_2`` are the along-track
    components at each pass separately.  The drift from equation (2) is the
    mean over the interval, but each pass's line-of-sight rate is
    instantaneous, so each of these is the tilt at its own pass, in error only
    by the difference between the interval-mean drift and the drift at that
    moment, divided by the platform speed.  When the surface evolves between
    the passes they are the better estimates, and their difference measures
    the evolution.  ``slope_along`` is their mean, or the first alone.  ``doppler_change_mps`` is the change between the passes in
    the line-of-sight rate each implies along the first pass's direction
    (``r2/mu2 - r1/mu1`` for two passes in the same direction).  Over a surface
    that held still it is noise only, so it is the evolution check.  It is NaN
    for a crossing pair, whose two rates differ by design.
    """

    drift_mps: float
    slope_x: float
    slope_y: float
    slope_along: float
    crossing: bool
    doppler_change_mps: float
    flagged: bool
    reason: str
    slope_along_1: float = float("nan")
    slope_along_2: float = float("nan")

    @property
    def tilt_deg(self) -> float:
        s = np.hypot(self.slope_x, self.slope_y) if self.crossing else abs(self.slope_along)
        return float(np.degrees(np.arcsin(min(1.0, s))))


def group_path_from_two_frequencies(field_f: complex, field_f2: complex, df_hz: float) -> float:
    """Apparent range from the phase of one echo at two nearby frequencies.

    With the field's phase ``-4 pi f P / c``, its derivative in frequency is
    ``-4 pi P' / c``, which is the group path.  ``df_hz`` must keep the phase
    step below ``pi``, which means ``P' < c / (4 df)``.
    """
    dphi = float(np.angle(complex(field_f2) * np.conj(complex(field_f))))
    return -C_LIGHT * dphi / (4.0 * np.pi * float(df_hz))


def los_rate_from_field(t_s: np.ndarray, field: np.ndarray, f_hz: float) -> float:
    """Phase-path rate from the power-weighted mean phase step between sweeps."""
    e = np.asarray(field, dtype=complex)
    step = np.sum(e[1:] * np.conj(e[:-1]))
    dt = float(np.mean(np.diff(np.asarray(t_s, dtype=float))))
    return float(C_LIGHT * np.angle(step) / (4.0 * np.pi * float(f_hz) * dt))


def estimate_second_look(
    p1: PassMeasurement,
    p2: PassMeasurement,
    *,
    crossing_min_deg: float = 10.0,
    evolution_tol_mps: float | None = None,
    doppler: str = "mean",
    facet: bool = False,
) -> SecondLookEstimate:
    """Drift and tilt from two passes over the same point.

    Parameters
    ----------
    crossing_min_deg
        The two velocities must differ in direction by at least this much for
        the pair to be solved as crossing tracks.  Below it the pair is solved
        as same-track, along the mean direction.
    evolution_tol_mps
        If given, a same-track pair whose ``|r2/mu2 - r1/mu1|`` exceeds it is
        flagged: the surface changed between the passes, and the estimate is
        of the mean over the interval at best.
    doppler
        ``"mean"`` uses both passes' line-of-sight rates, which refers the
        estimate to the middle of the interval.  ``"first"`` uses only the
        first pass's, which refers the tilt to the first pass.
    facet
        True for a compact facet measured at its power peak, where the
        spacecraft is on the normal and the range is the vertical gap over
        ``cos(theta)`` rather than times it.  The difference is second order in
        the tilt, 0.03% at 2 degrees.
    """
    dt = float(p2.t_s - p1.t_s)
    if dt <= 0.0:
        raise ValueError("the second pass must follow the first")
    v1 = np.asarray(p1.velocity_mps, dtype=float)
    v2 = np.asarray(p2.velocity_mps, dtype=float)
    g = 0.5 * (1.0 / p1.mu_sat + 1.0 / p2.mu_sat)   # mean group index at the spacecraft
    dp = float(p2.group_path_m - p1.group_path_m)
    q1 = p1.los_rate_mps / p1.mu_sat
    q2 = p2.los_rate_mps / p2.mu_sat
    # The angle between the two tracks as lines: an ascending and a descending
    # pass along nearly the same line are a same-track pair, not a crossing.
    cosang = abs(float(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2))))
    crossing = np.degrees(np.arccos(np.clip(cosang, -1.0, 1.0))) >= float(crossing_min_deg)
    u1 = v1 / np.linalg.norm(v1)

    cos_t = 1.0
    sx = sy = s_al = float("nan")
    change = float("nan")
    a1 = a2 = float("nan")
    for _ in range(3):
        w = -dp * cos_t / (dt * g) if facet else -dp / (dt * g * cos_t)
        if crossing:
            a = np.vstack([v1, v2])
            rhs = np.array([w * cos_t - q1, w * cos_t - q2])
            sx, sy = (float(x) for x in np.linalg.solve(a, rhs))
            s_mag = float(np.hypot(sx, sy))
            s_al = float(np.dot([sx, sy], u1))
        else:
            # Each pass's along-track component, both referred to the first
            # pass's direction so that an antiparallel pair is handled.
            a1 = (w * cos_t - q1) / float(np.dot(v1, u1))
            a2 = (w * cos_t - q2) / float(np.dot(v2, u1))
            s_al = 0.5 * (a1 + a2) if doppler == "mean" else a1
            change = (a1 - a2) * float(np.linalg.norm(v1))
            s_mag = abs(s_al)
        cos_t = float(np.sqrt(max(0.0, 1.0 - s_mag * s_mag)))

    flagged = False
    reason = "crossing" if crossing else "same track"
    if (not crossing) and evolution_tol_mps is not None and abs(change) > float(evolution_tol_mps):
        flagged = True
        reason = f"line-of-sight rate changed by {change:.1f} m/s between passes; surface evolved"
    return SecondLookEstimate(
        drift_mps=float(w), slope_x=sx, slope_y=sy, slope_along=float(s_al), crossing=bool(crossing),
        slope_along_1=float(a1), slope_along_2=float(a2),
        doppler_change_mps=float(change), flagged=flagged, reason=reason,
    )


def horizontal_drift_from_peaks(x_peak1_m: float, x_peak2_m: float, dt_s: float,
                                drift_mps: float, slope_along: float) -> float:
    """Along-track drift of a compact facet from the shift of its power peak.

    The spacecraft is on the facet's normal at the power peak, which is
    ``X(q)`` along track from the facet.  A facet that rises by ``w dt`` brings
    that point ``w dt tan(theta)`` nearer, so

        u = (x_pk2 - x_pk1) / dt + w tan(theta).
    """
    tan_t = slope_along / np.sqrt(max(1e-12, 1.0 - slope_along**2))
    return float((x_peak2_m - x_peak1_m) / float(dt_s) + float(drift_mps) * tan_t)
