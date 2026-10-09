"""Layer: L3
WP: WP5
Claims: C3
Status: implemented (parameter level)
Requirements: SCI-15

Line-of-sight drift of a field-aligned reflector, with its orientation taken
from the field model.

A field-aligned sheet's normal is perpendicular to B, so the field model fixes
the sheet's orientation except for a rotation ``phi`` about B.  Given the
platform ephemeris and the field line the echo comes from (which the field
model and the measured range identify), the platform's own contribution to the
echo's phase path can be predicted for any assumed ``phi`` and removed before
the Doppler is sampled.  What is left is the drift, plus whatever the assumed
``phi`` gets wrong.

Three ways of choosing ``phi`` are provided, so that a study can compare them:

``assume``
    ``phi = 0``: the sheet's normal is taken to lie in the magnetic meridian.
``fit``
    ``phi`` is chosen to make the compensated Doppler history flattest, on the
    argument that drift is constant over the window and a wrong ``phi`` leaves
    a time-varying residual.
``given``
    ``phi`` is supplied from outside, as a cross-track arrival angle (route
    three) would supply it.

Removing the predicted platform phase first is also what keeps the Doppler
inside the 22.5 Hz Nyquist limit of the 45 Hz drift mode.  Without it, a
field-aligned sheet off the dip equator aliases (the platform term is
``v_s sin I`` along a meridional track, 1.3 km/s at a dip of 10 deg).
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np

from tss.constants import C_LIGHT
from tss.l1_propagation import sheet_echo as se
from tss.l1_propagation.plate_echo import ExponentialTopside
from tss.l3_signal.novel.route_one_estimator import doppler_history_from_field

__all__ = ["PathPrediction", "predict_path", "FieldAlignedEstimate", "estimate_field_aligned"]


@dataclass(frozen=True)
class PathPrediction:
    """Predicted one-way phase path and the ray direction at the platform."""

    path_m: np.ndarray
    ray_at_platform: np.ndarray  # (n, 3), unit, from the reflector toward the platform


def _stratified_path(profile: ExponentialTopside, f_hz: float, platform: np.ndarray,
                     point: np.ndarray, nodes, weights) -> tuple[float, float]:
    """Phase path and Snell invariant of the ray from platform to point.

    Depth is measured below the platform's height at ``t = 0``; the caller
    passes coordinates in the medium frame with that origin at ``z = 0`` of
    the platform start.
    """
    rho = np.array([np.hypot(platform[0] - point[0], platform[1] - point[1])])
    p, _, q = se._rays(rho, float(-platform[2]), np.array([float(-point[2])]), f_hz, profile,
                       nodes, weights, 12, None)
    return float(p[0]), float(q[0])


def predict_path(
    t_s: np.ndarray, platform_positions: np.ndarray, sheet: se.FieldAlignedSheet, *,
    f_hz: float, profile: ExponentialTopside | None = None, height_m: float | None = None,
) -> PathPrediction:
    """Phase path to the specular point of a static sheet, sweep by sweep.

    In free space the specular point is found in closed form
    (:func:`tss.l1_propagation.sheet_echo.specular_point`).  In the plasma it
    is the stationary point of the phase path over the sheet (Fermat), found
    by minimisation over the sheet coordinates.
    """
    t = np.asarray(t_s, dtype=float)
    pos = np.asarray(platform_positions, dtype=float)
    paths = np.empty(t.size)
    rays = np.empty((t.size, 3))
    if profile is None:
        for j in range(t.size):
            sp = se.specular_point(sheet, pos[j])
            paths[j] = sp["range_m"]
            d = pos[j] - sp["point"]
            rays[j] = d / np.linalg.norm(d)
        return PathPrediction(paths, rays)

    from scipy.optimize import minimize

    if height_m is None:
        raise ValueError("height_m is needed in the plasma, to place the depth origin")
    fr = se.sheet_frame(sheet)
    nodes, weights = np.polynomial.legendre.leggauss(48)
    mu_sat = float(np.sqrt(profile.mu_squared(f_hz, 0.0)))
    # Depth origin: the platform start height, which is height_m above the sheet reference.
    z0 = float(height_m)
    guess = np.array([0.0, 0.0])
    for j in range(t.size):
        plat = pos[j].copy()

        def fun(ac, plat=plat):
            line, _ = se._centre_line(sheet, np.array(ac[0]), fr)
            point = line + ac[1] * fr["u"]
            p, _ = _stratified_path(profile, f_hz, plat - [0, 0, z0], point - [0, 0, z0], nodes, weights)
            return p

        simplex = guess + np.array([[0.0, 0.0], [500.0, 0.0], [0.0, 500.0]])
        res = minimize(fun, guess, method="Nelder-Mead",
                       options={"xatol": 1e-3, "fatol": 1e-9, "maxiter": 4000,
                                "initial_simplex": simplex})
        guess = res.x
        paths[j] = float(res.fun)
        line, _ = se._centre_line(sheet, np.array(res.x[0]), fr)
        point = line + res.x[1] * fr["u"]
        _, q = _stratified_path(profile, f_hz, plat - [0, 0, z0], point - [0, 0, z0], nodes, weights)
        hdir = plat[:2] - point[:2]
        hn = np.linalg.norm(hdir)
        hdir = hdir / hn if hn > 0 else np.zeros(2)
        sin_t = q / mu_sat
        rays[j] = np.array([hdir[0] * sin_t, hdir[1] * sin_t, np.sqrt(max(1.0 - sin_t**2, 0.0))])
    return PathPrediction(paths, rays)


@dataclass(frozen=True)
class FieldAlignedEstimate:
    method: str
    rotation_rad: float
    los_mps: float
    residual_rms_hz: float
    max_abs_residual_doppler_hz: float


def _compensated(t, field, f_hz, path):
    k0 = 2.0 * np.pi * f_hz / C_LIGHT
    return doppler_history_from_field(t, field * np.exp(2j * k0 * path))


def estimate_field_aligned(
    t_s: np.ndarray, field: np.ndarray, platform_positions: np.ndarray,
    sheet_model: se.FieldAlignedSheet, *, f_hz: float, mu_sat: float, method: str,
    rotation_rad: float = 0.0, profile: ExponentialTopside | None = None,
    height_m: float | None = None, fit_span_rad: float = np.radians(15.0),
) -> FieldAlignedEstimate:
    """Line-of-sight drift of a field-aligned sheet.

    ``sheet_model`` carries the field-line geometry (dip and curvature); its
    own rotation and velocity are ignored.  ``method`` is ``assume``, ``fit``
    or ``given``; ``rotation_rad`` is used only by ``given``.
    """
    t = np.asarray(t_s, dtype=float)

    def residual(phi: float):
        model = replace(sheet_model, rotation_rad=float(phi), velocity_n_mps=0.0, velocity_u_mps=0.0)
        pred = predict_path(t, platform_positions, model, f_hz=f_hz, profile=profile, height_m=height_m)
        hist = _compensated(t, field, f_hz, pred.path_m)
        w = hist.power / hist.power.sum()
        mean = float(np.sum(w * hist.doppler_hz))
        rms = float(np.sqrt(np.sum(w * (hist.doppler_hz - mean) ** 2)))
        return mean, rms, float(np.max(np.abs(hist.doppler_hz)))

    if method == "assume":
        phi = 0.0
    elif method == "given":
        phi = float(rotation_rad)
    elif method == "fit":
        from scipy.optimize import minimize_scalar

        grid = np.linspace(-fit_span_rad, fit_span_rad, 31)
        vals = [residual(g)[1] for g in grid]
        i = int(np.argmin(vals))
        lo, hi = grid[max(i - 1, 0)], grid[min(i + 1, grid.size - 1)]
        res = minimize_scalar(lambda g: residual(g)[1], bounds=(lo, hi), method="bounded",
                              options={"xatol": 1e-6})
        phi = float(res.x)
    else:
        raise ValueError(f"unknown method {method!r}")

    mean, rms, mx = residual(phi)
    return FieldAlignedEstimate(method, phi, C_LIGHT * mean / (2.0 * f_hz * mu_sat), rms, mx)
