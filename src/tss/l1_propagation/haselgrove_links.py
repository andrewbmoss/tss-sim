"""Layer: L1
WP: WP5, WP7
Claims: C3
Status: implemented (verification bridge)
Requirements: SCI-15, SCI-27

Ray-traced truth for vertical and oblique links through a tilted, drifting
ionosphere, using the independent three-dimensional Haselgrove ray tracer
written for the Chartier reproduction (``nisc-verify``).

That ray tracer shares no formulation with the virtual-mirror model of
:mod:`tss.l1_propagation.bistatic_mirror`: it integrates the ray equations
through a Chapman ionosphere on a sphere, with no mirror, no equivalence
theorem and no secant law.  This module adds three things it needs to serve as
truth for drift and tilt.

1.  **A tilted, drifting medium.**  Every isodensity surface is displaced by
    ``gE x + gN y + w (t - t0)``, in the same gnomonic tangent coordinates
    about a patch centre that the mirror model uses.
2.  **Precise homing.**  Gauss-Newton on launch elevation and bearing until the
    ray crosses the receiver's altitude within a few centimetres of it.  The
    Nelder-Mead homing of the reproduction stops at a kilometre, which is
    enough to count links but tilts the arrival direction by a milliradian,
    seven metres per second of Doppler at orbital speed.
3.  **The line-of-sight rate from the ray itself.**  For a ray from A to B the
    phase path changes at

        dP/dt = -mu_A k_A . v_A + mu_B k_B . v_B - (1/2) integral dX/dt dP'

    where ``k_A`` is the launch direction at A and ``k_B`` the arrival
    direction at B.  By reciprocity in an isotropic medium the arrival direction
    at B is minus the launch direction of the ray homed from B back to A, so
    both are known to the precision of the homing.  The last term is the
    medium's own change along the fixed path, which Fermat's principle allows,
    and written against group path it stays finite through the reflection.

``nisc`` is imported lazily from ``NISC_VERIFY_SRC`` or the reproduction's
default location; :func:`available` says whether it was found.
"""

from __future__ import annotations

import datetime as dt
import os
import sys
from dataclasses import dataclass

import numpy as np

__all__ = ["available", "TiltedDriftingIonosphere", "TracedLink", "trace_link", "home_precise"]

_DEFAULT_PATHS = (
    r"C:\Users\Andrew\PhD\Chartier Verification\nisc-verify\src",
    "/home/claude/work/nisc-verify/src",
)


def _import_nisc():
    for p in [os.environ.get("NISC_VERIFY_SRC", "")] + list(_DEFAULT_PATHS):
        if p and os.path.isdir(p) and p not in sys.path:
            sys.path.insert(0, p)
    import nisc  # noqa: F401
    from nisc.geodesy import relaz2ecef
    from nisc.iono.base import IonosphereModel
    from nisc.iono.chapman import ChapmanIonosphere
    from nisc.raytrace.haselgrove import STATUS_RETURNED, TraceOptions, trace
    from nisc.raytrace.medium import Medium
    return IonosphereModel, ChapmanIonosphere, STATUS_RETURNED, TraceOptions, trace, Medium, relaz2ecef


def available() -> bool:
    try:
        _import_nisc()
        return True
    except Exception:  # pragma: no cover - depends on the machine
        return False


R_E = 6371.0e3
K806 = 80.6163849431291  # f_p^2 = K806 Ne


def _tangent(lat, lon, lat0, lon0):
    la, lo = np.radians(lat), np.radians(lon)
    d = np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], axis=-1)
    la0, lo0 = np.radians(lat0), np.radians(lon0)
    e = np.array([-np.sin(lo0), np.cos(lo0), 0.0])
    n = np.array([-np.sin(la0) * np.cos(lo0), -np.sin(la0) * np.sin(lo0), np.cos(la0)])
    u = np.array([np.cos(la0) * np.cos(lo0), np.cos(la0) * np.sin(lo0), np.sin(la0)])
    s = R_E / (d @ u)
    return (d @ e) * s, (d @ n) * s


def TiltedDriftingIonosphere(base, lat0, lon0, g_east, g_north, w_mps, t0: dt.datetime):
    """A copy of ``base`` with every isodensity surface displaced as a tilted, drifting plane."""
    IonosphereModel = _import_nisc()[0]

    class _Tilted(IonosphereModel):
        name = "tilted_drifting"

        def electron_density(self, lat_deg, lon_deg, alt_m, time):
            lat = np.asarray(lat_deg, dtype=float)
            lon = np.asarray(lon_deg, dtype=float)
            x, y = _tangent(lat, lon, lat0, lon0)
            shift = g_east * x + g_north * y + w_mps * (time - t0).total_seconds()
            return base.electron_density(lat, lon, np.asarray(alt_m, dtype=float) - shift, time)

        def metadata(self):
            return {"name": self.name, "base": base.metadata(), "lat0": lat0, "lon0": lon0,
                    "g_east": g_east, "g_north": g_north, "w_mps": w_mps}

    return _Tilted()


def _options(TraceOptions, alt):
    return TraceOptions(step_init=500.0, step_min=1.0, step_max=5000.0, rtol=1e-10, atol_pos=0.01,
                        return_alt=alt, max_group_path=8000e3, max_steps=60000, record=True)


def _crossing(bundle, i, alt):
    """Landing lat/lon and group path where ray ``i`` crosses ``alt`` going up."""
    pts = np.asarray(bundle.paths[i])
    gp = np.asarray(bundle.path_group[i])
    a1, a2 = pts[-2, 2], pts[-1, 2]
    f = (alt - a1) / (a2 - a1) if a2 != a1 else 1.0
    lat = pts[-2, 0] + f * (pts[-1, 0] - pts[-2, 0])
    lon = pts[-2, 1] + f * (pts[-1, 1] - pts[-2, 1])
    g = gp[-2] + f * (gp[-1] - gp[-2])
    return lat, lon, g


def _local_offset(lat, lon, lat_t, lon_t):
    """East and north offset (m) of (lat, lon) from (lat_t, lon_t), small separations."""
    dn = np.radians(lat - lat_t) * R_E
    de = np.radians(((lon - lon_t + 180.0) % 360.0) - 180.0) * R_E * np.cos(np.radians(lat_t))
    return np.array([de, dn])


def _angles(p, vertical):
    """Launch (elevation, bearing) from the homing variables.

    Oblique rays are homed in elevation and bearing directly.  A near-vertical
    ray is homed in the east and north components of its tilt from the nadir,
    because bearing is undefined at the nadir and the Jacobian would be
    singular there.
    """
    if not vertical:
        return p[0], p[1]
    tilt = np.hypot(p[0], p[1])
    return -(90.0 - tilt), float(np.degrees(np.arctan2(p[0], p[1]))) % 360.0


def home_precise(medium, tx, rx, f_hz, elev0, bearing0, *, tol_m: float = 0.05, iterations: int = 25):
    """Gauss-Newton homing.  Returns (elevation, bearing, the homed ray's bundle)."""
    _, _, RET, TraceOptions, trace, _, _ = _import_nisc()
    opt = _options(TraceOptions, rx[2])
    vertical = abs(tx[0] - rx[0]) < 1e-9 and abs(tx[1] - rx[1]) < 1e-9
    if vertical:
        tilt0 = 90.0 + float(elev0)
        p = np.array([tilt0 * np.sin(np.radians(bearing0)), tilt0 * np.cos(np.radians(bearing0))])
    else:
        p = np.array([float(elev0), float(bearing0)])
    d = 2e-4
    for _ in range(int(iterations)):
        trial = [p, p + [d, 0.0], p + [0.0, d]]
        ang = np.array([_angles(q, vertical) for q in trial])
        b = trace(medium, tx[0], tx[1], tx[2], ang[:, 0], ang[:, 1], f_hz, options=opt)
        if np.any(b.status != RET):
            raise RuntimeError(f"ray did not return at {f_hz/1e6:.2f} MHz (status {b.status})")
        res = []
        for i in range(3):
            la, lo, _ = _crossing(b, i, rx[2])
            res.append(_local_offset(la, lo, rx[0], rx[1]))
        r0 = res[0]
        if np.hypot(*r0) < tol_m:
            el, az = _angles(p, vertical)
            return el, az, b
        jac = np.column_stack([(res[1] - r0) / d, (res[2] - r0) / d])
        step = np.clip(np.linalg.solve(jac, -r0), -2.0, 2.0)
        p = p + step
    raise RuntimeError(f"homing did not converge, miss {np.hypot(*r0):.3f} m")


@dataclass(frozen=True)
class TracedLink:
    group_path_m: float
    los_rate_mps: float
    rate_endpoints_mps: float
    rate_medium_mps: float
    mu_a: float
    mu_b: float
    apex_alt_m: float
    elev_a_deg: float
    elev_b_deg: float


def _launch(relaz2ecef, el, az, lat, lon):
    k = relaz2ecef(el, az, lat, lon)
    return np.asarray(k, float) / np.linalg.norm(k)


def trace_link(iono, time: dt.datetime, tx, rx, va, vb, f_hz, seed_a, seed_b, *, dt_medium_s: float = 1.0) -> TracedLink:
    """Group path and line-of-sight rate of the link tx -> rx at one instant.

    ``tx`` and ``rx`` are (lat, lon, alt); ``va`` and ``vb`` Earth-fixed velocity
    vectors (ECEF, m/s); ``seed_a`` and ``seed_b`` launch (elevation, bearing)
    guesses at each end.
    """
    IonosphereModel, ChapmanIonosphere, RET, TraceOptions, trace, Medium, relaz2ecef = _import_nisc()
    med = Medium(ionosphere=iono, time=time, collisions="none", spherical_earth=True)
    ea, aa, ba = home_precise(med, tx, rx, f_hz, *seed_a)
    vertical = abs(tx[0] - rx[0]) < 1e-9 and abs(tx[1] - rx[1]) < 1e-9
    if vertical:
        eb, ab = ea, aa
    else:
        eb, ab, _ = home_precise(med, rx, tx, f_hz, *seed_b)
    ka = _launch(relaz2ecef, ea, aa, tx[0], tx[1])
    kb = _launch(relaz2ecef, eb, ab, rx[0], rx[1])    # = minus the arrival direction at B
    if vertical:
        # A monostatic echo retraces its own path.
        kb = ka
    ne_a = float(iono.electron_density(np.array([tx[0]]), np.array([tx[1]]), np.array([tx[2]]), time)[0])
    ne_b = float(iono.electron_density(np.array([rx[0]]), np.array([rx[1]]), np.array([rx[2]]), time)[0])
    mu_a = float(np.sqrt(max(1.0 - K806 * ne_a / f_hz**2, 0.0)))
    mu_b = float(np.sqrt(max(1.0 - K806 * ne_b / f_hz**2, 0.0)))
    rate_end = mu_a * float(np.dot(ka, va)) + mu_b * float(np.dot(kb, vb))
    # Medium term along the forward ray.
    pts = np.asarray(ba.paths[0])
    gp = np.asarray(ba.path_group[0])
    la, lo, al = pts[:, 0], pts[:, 1], pts[:, 2]
    tp = time + dt.timedelta(seconds=dt_medium_s)
    tm = time - dt.timedelta(seconds=dt_medium_s)
    dn = (iono.electron_density(la, lo, al, tp) - iono.electron_density(la, lo, al, tm)) / (2 * dt_medium_s)
    dxdt = K806 * dn / f_hz**2
    _, _, g_land = _crossing(ba, 0, rx[2])
    keep = gp <= g_land
    gpk = np.append(gp[keep], g_land)
    dk = np.append(dxdt[keep], dxdt[keep][-1])
    rate_med = 0.5 * float(np.trapezoid(dk, gpk))
    return TracedLink(group_path_m=float(g_land), los_rate_mps=rate_end + rate_med,
                      rate_endpoints_mps=rate_end, rate_medium_mps=rate_med, mu_a=mu_a, mu_b=mu_b,
                      apex_alt_m=float(ba.apex_alt[0]), elev_a_deg=ea, elev_b_deg=eb)
