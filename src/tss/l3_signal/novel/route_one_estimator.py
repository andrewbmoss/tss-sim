"""Layer: L3
WP: WP5
Claims: C3
Status: implemented (parameter level)
Requirements: SCI-15, SCI-23

The route-one estimator of Technical Volume Section 8.6, working on the
complex echo of one detected scatterer, one sample per sweep.

This is the first estimator code of WP5.  It takes the per-sweep complex field
of one echo, as the range-gated output of the baseline chain would give it,
and does three things.

1.  It measures the Doppler history from the phase step between consecutive
    sweeps and the echo power from the magnitude.
2.  It decides whether the echo is separable at all.  An extended layer's
    Doppler history does not sweep and its power has no peak, and the
    estimator must say so rather than return a drift (TV Section 8.6, last
    paragraph of the estimator description).
3.  If the echo is separable it takes the overflight time from the peak of the
    echo power, the zero crossing from the Doppler history with the slope
    imported from the ephemeris, and the drift from

        v_drift = v_s^2 (t_zc - t_pk) / (mu_sat P').

It works at the parameter level: the input is the echo of one range cell, not
a raw record.  The signal-level version, which forms that echo from de-chirped
samples through L2 and L3, is the next step of WP5 and depends on the O/X and
dispersion work of WP4.

The model this estimator assumes is the one TV equation 8.x states.  The
verification study (:mod:`tss.studies.c3_verify`) tests it against a forward
model that does not assume it.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from tss.constants import C_LIGHT

__all__ = [
    "DopplerHistory",
    "RouteOneEstimate",
    "doppler_history_from_field",
    "free_t0_degeneracy",
    "estimate_route_one",
]


@dataclass(frozen=True)
class DopplerHistory:
    """Doppler and power at the midpoints between sweeps."""

    t_s: np.ndarray
    doppler_hz: np.ndarray
    power: np.ndarray


@dataclass(frozen=True)
class RouteOneEstimate:
    """What route one returns for one echo.

    When ``separable`` is false the drift, the overflight time and the zero
    crossing are NaN and only ``los_rate_mps`` is reported.  That is the
    line-of-sight rate the echo carries, ``c f_D / (2 f mu_sat)``, which for
    an extended layer is drift and tilt together and cannot be split by one
    pass.
    """

    separable: bool
    reason: str
    slope_ratio: float
    peak_contrast: float
    t_power_peak_s: float
    t_zero_crossing_s: float
    drift_mps: float
    along_track_m: float
    los_rate_mps: float


def doppler_history_from_field(t_s: np.ndarray, field: np.ndarray, *, main_lobe: float = 0.1) -> DopplerHistory:
    """Doppler from the phase step between sweeps, and power, at midpoints.

    Samples whose power is below ``main_lobe`` of the peak on either side of a
    step are dropped, because across a null of the pattern the phase jumps by
    ``pi`` and the step is not a Doppler shift.
    """
    t = np.asarray(t_s, dtype=float)
    e = np.asarray(field, dtype=complex)
    p = np.abs(e) ** 2
    step = e[1:] * np.conj(e[:-1])
    dt = np.diff(t)
    fd = np.angle(step) / (2.0 * np.pi * dt)
    tm = 0.5 * (t[1:] + t[:-1])
    pm = np.sqrt(p[1:] * p[:-1])
    keep = np.minimum(p[1:], p[:-1]) >= float(main_lobe) * p.max()
    return DopplerHistory(t_s=tm[keep], doppler_hz=fd[keep], power=pm[keep])


def free_t0_degeneracy(
    history: DopplerHistory,
    f_hz: float,
    mu_sat: float,
    slope_hz_per_s: float,
    drifts_mps: np.ndarray,
) -> np.ndarray:
    """Best misfit over ``t0`` at each trial drift, with ``t0`` left free.

    For each trial drift the overflight time that best fits the history is
    found by weighted least squares and the residual sum of squares, in
    Hz^2 per sample, is returned.  If drift and overflight time are not
    separable the result is flat in drift.
    """
    t, fd, w = history.t_s, history.doppler_hz, history.power
    w = w / w.sum()
    out = np.empty(len(drifts_mps))
    for i, v in enumerate(np.asarray(drifts_mps, dtype=float)):
        c = 2.0 * float(f_hz) * float(mu_sat) * v / C_LIGHT
        # f_D = c - S (t - t0)  =>  f_D - c + S t = S t0
        t0 = float(np.sum(w * (fd - c + slope_hz_per_s * t)) / slope_hz_per_s)
        r = fd - (c - slope_hz_per_s * (t - t0))
        out[i] = float(np.sum(w * r * r))
    return out


def _fit_peak(t: np.ndarray, p: np.ndarray) -> tuple[float, float]:
    """Centre and contrast of the power envelope.

    The centre comes from the ``sinc^2`` fit of the route-one study, which is
    the fit TV Section 8.8 sizes.  Contrast is the fitted peak over the larger
    of the two window edges, so a value near one means there is no peak inside
    the window.
    """
    from tss.studies.route_one import fit_overflight_time

    centre, _ = fit_overflight_time(t, p)
    k = max(3, int(round(0.05 * t.size)))
    edge = max(float(p[:k].mean()), float(p[-k:].mean()), 1e-300)
    return centre, float(p.max() / edge)


def estimate_route_one(
    t_s: np.ndarray,
    field: np.ndarray,
    *,
    f_hz: float,
    group_path_m: float,
    mu_sat: float,
    v_s_mps: float = 7600.0,
    slope_tolerance: float = 0.2,
    min_contrast: float = 4.0,
) -> RouteOneEstimate:
    """Route one on the echo of one scatterer.

    Parameters
    ----------
    t_s, field
        Sweep times and the complex echo of the scatterer's range cell.
    f_hz
        Sounding frequency.
    group_path_m
        Apparent range of the echo at its power peak, as measured.
    mu_sat
        Refractive index at the platform.
    slope_tolerance
        The free-fit Doppler slope must lie within this fraction of the
        ephemeris slope ``(2 f / c) v_s^2 / P'`` for the echo to count as
        compact.
    min_contrast
        The power must fall to below ``1 / min_contrast`` of its peak at both
        edges of the window, so that the peak is inside it.
    """
    t = np.asarray(t_s, dtype=float)
    e = np.asarray(field, dtype=complex)
    f = float(f_hz)
    slope = 2.0 * f * float(v_s_mps) ** 2 / (C_LIGHT * float(group_path_m))
    hist = doppler_history_from_field(t, e)
    w = hist.power / hist.power.sum()

    # Free slope, power weighted.
    tm = float(np.sum(w * hist.t_s))
    fm = float(np.sum(w * hist.doppler_hz))
    var = float(np.sum(w * (hist.t_s - tm) ** 2))
    free_slope = -float(np.sum(w * (hist.t_s - tm) * (hist.doppler_hz - fm))) / var if var > 0 else 0.0
    ratio = free_slope / slope

    los = C_LIGHT * fm / (2.0 * f * float(mu_sat))
    power = np.abs(e) ** 2
    t_pk, contrast = _fit_peak(t, power)

    nan = float("nan")
    if abs(ratio - 1.0) > float(slope_tolerance):
        return RouteOneEstimate(
            False, f"Doppler slope is {ratio:.2f} of the ephemeris slope; extended reflector",
            ratio, contrast, nan, nan, nan, nan, los,
        )
    if contrast < float(min_contrast):
        return RouteOneEstimate(
            False, f"power contrast {contrast:.1f} below {min_contrast}; no peak in window",
            ratio, contrast, nan, nan, nan, nan, los,
        )

    # Zero crossing with the slope fixed by the ephemeris: f_D = S (t_zc - t).
    t_zc = float(np.sum(w * (hist.doppler_hz / slope + hist.t_s)))
    drift = float(v_s_mps) ** 2 * (t_zc - t_pk) / (float(mu_sat) * float(group_path_m))
    return RouteOneEstimate(
        True, "compact", ratio, contrast, t_pk, t_zc, drift, float(v_s_mps) * t_pk, los,
    )
