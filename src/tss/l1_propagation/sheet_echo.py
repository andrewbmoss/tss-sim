"""Layer: L1
WP: WP5
Claims: C3
Status: implemented
Requirements: SCI-15

Physical-optics echo of a field-aligned sheet seen from a moving platform, in
three dimensions, in free space and inside a horizontally stratified plasma.

Why a sheet
-----------
Topside irregularities are field aligned and extend for large distances along
field lines; individual structures are typically about 1.3 km thick across the
field (Dyson 1969).  A reflecting structure of that kind is a sheet that
contains the field direction: thin across B, long along it.  Its normal is
therefore perpendicular to B, which fixes one component of its orientation from
the field model.  The other, a rotation about B, is free.

Geometry
--------
Local frame: ``x`` magnetic north (horizontal), ``y`` magnetic east, ``z`` up.
The field at the sheet's reference point points along
``b = (cos I, 0, -sin I)`` (dip ``I`` positive in the northern hemisphere).
``n0 = (sin I, 0, cos I)`` is the upward unit vector perpendicular to ``b`` in
the magnetic meridian and ``e = (0, 1, 0)``.  A rotation ``phi`` about ``b``
gives the sheet's across-field direction ``u = -sin(phi) n0 + cos(phi) e``.

The sheet follows the field line: its centre line is a circle of radius
``Rc`` in the meridian plane, tangent to ``b`` at the reference point and
curving toward ``-n0`` (a dipole field line is convex upward near the equator).
The sheet is swept by a segment of width ``W`` along ``u`` carried along that
line.  A point is

    X(a, c) = Rc sin(a / Rc) b + Rc (cos(a / Rc) - 1) n0 + c u,

and its unit normal is ``t_a x u`` with ``t_a = cos(a / Rc) b - sin(a / Rc) n0``.

The sheet moves rigidly with velocity ``V``, perpendicular to B.  In the plasma
the profile moves with it (frozen in), which is the bulk-drift case.

The platform moves horizontally at ``v_s`` on heading ``psi`` from magnetic
north.  Its start is placed on the normal through the reference point, at
height ``h`` above it, so the echo is centred on ``t = 0``.

Inside the plasma
-----------------
The profile is field free and horizontally stratified, so each ray lies in the
vertical plane through the platform and the point.  The Snell invariant ``q``
is found by Newton's method from the horizontal distance, exactly as in
:mod:`tss.l1_propagation.plate_echo`, and the phase path follows.

What is left out
----------------
The magnetic field's effect on the wave (the sheet is field aligned but the
propagation is field free), absorption, reflectivity variations across the
sheet, the change of field magnitude along the sheet, and partial reflection
from the sheet's interior.  The sheet is a perfect specular surface.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from tss.constants import C_LIGHT
from tss.l1_propagation.plate_echo import ExponentialTopside

__all__ = [
    "FieldAlignedSheet",
    "SheetEcho",
    "dipole_curvature_radius_m",
    "sheet_frame",
    "platform_start",
    "platform_start_stratified",
    "echo_free_space",
    "echo_stratified",
    "specular_point",
]

EARTH_RADIUS_M = 6371.2e3


def dipole_curvature_radius_m(dip_rad: float, radius_m: float) -> float:
    """Radius of curvature of a dipole field line at the point with dip ``I``.

    The field line is ``r = L sin^2(theta)``.  ``tan I = 2 cot(theta)`` gives
    the colatitude; the curvature is evaluated by differentiating the curve.
    ``radius_m`` is the geocentric distance of the point.  At the dip equator
    the result is ``r / 3``.
    """
    theta = np.arctan2(2.0, np.tan(abs(dip_rad))) if dip_rad != 0 else np.pi / 2
    shell = radius_m / np.sin(theta) ** 2
    st, ct = np.sin(theta), np.cos(theta)
    x1 = 3.0 * shell * st**2 * ct
    z1 = shell * (2.0 * st * ct**2 - st**3)
    x2 = 3.0 * shell * (2.0 * st * ct**2 - st**3)
    z2 = shell * (2.0 * ct**3 - 7.0 * st**2 * ct)
    k = abs(x1 * z2 - z1 * x2) / (x1**2 + z1**2) ** 1.5
    return float(1.0 / k)


@dataclass(frozen=True)
class FieldAlignedSheet:
    """A field-aligned reflecting sheet.  SI units, angles in rad."""

    dip_rad: float = 0.0
    rotation_rad: float = 0.0
    width_m: float = 1.3e3
    length_m: float = 60e3
    curvature_radius_m: float = np.inf
    tukey_fraction: float = 0.3
    velocity_n_mps: float = 0.0
    velocity_u_mps: float = 0.0


def sheet_frame(sheet: FieldAlignedSheet) -> dict[str, np.ndarray]:
    """Unit vectors b, n0, e, u and the reference normal n."""
    i = float(sheet.dip_rad)
    b = np.array([np.cos(i), 0.0, -np.sin(i)])
    n0 = np.array([np.sin(i), 0.0, np.cos(i)])
    e = np.array([0.0, 1.0, 0.0])
    p = float(sheet.rotation_rad)
    u = -np.sin(p) * n0 + np.cos(p) * e
    n = np.cos(p) * n0 + np.sin(p) * e
    return {"b": b, "n0": n0, "e": e, "u": u, "n": n}


def _centre_line(sheet: FieldAlignedSheet, a: np.ndarray, fr: dict) -> tuple[np.ndarray, np.ndarray]:
    """Points and tangents of the centre line at arc length ``a``."""
    a = np.asarray(a, dtype=float)
    rc = float(sheet.curvature_radius_m)
    if np.isinf(rc):
        pts = a[..., None] * fr["b"]
        tan = np.broadcast_to(fr["b"], pts.shape).copy()
    else:
        s, cth = np.sin(a / rc), np.cos(a / rc)
        pts = rc * s[..., None] * fr["b"] + rc * (cth - 1.0)[..., None] * fr["n0"]
        tan = cth[..., None] * fr["b"] - s[..., None] * fr["n0"]
    return pts, tan


def _normal(tangent: np.ndarray, u: np.ndarray) -> np.ndarray:
    nrm = np.cross(tangent, u)
    nrm /= np.linalg.norm(nrm, axis=-1, keepdims=True)
    return np.where(nrm[..., 2:3] < 0, -nrm, nrm)


def sheet_velocity(sheet: FieldAlignedSheet) -> np.ndarray:
    fr = sheet_frame(sheet)
    return sheet.velocity_n_mps * fr["n"] + sheet.velocity_u_mps * fr["u"]


def platform_start(sheet: FieldAlignedSheet, height_m: float) -> np.ndarray:
    """Platform position at ``t = 0``: on the reference normal, ``h`` above."""
    n = sheet_frame(sheet)["n"]
    return (float(height_m) / n[2]) * n


def platform_start_stratified(
    sheet: FieldAlignedSheet, height_m: float, profile: ExponentialTopside, f_hz: float,
) -> np.ndarray:
    """Platform position at ``t = 0`` inside the plasma.

    The ray must meet the sheet's reference point along its normal, so its
    Snell invariant is ``q = mu_f |n_h|`` and its azimuth that of ``n_h``; the
    platform sits ``X(q)`` away horizontally and ``h`` above.  Refraction makes
    this differ from the straight-line start of :func:`platform_start`.
    """
    n = sheet_frame(sheet)["n"]
    nh = np.hypot(n[0], n[1])
    nodes, weights = np.polynomial.legendre.leggauss(64)
    d = 0.5 * height_m * (1.0 + nodes)
    mu2 = profile.mu_squared(f_hz, d)
    mu_f = float(np.sqrt(profile.mu_squared(f_hz, height_m)))
    q = mu_f * nh
    x = float(np.sum(0.5 * height_m * weights * q / np.sqrt(mu2 - q * q)))
    az = np.array([n[0], n[1]]) / nh if nh > 0 else np.zeros(2)
    return np.array([x * az[0], x * az[1], float(height_m)])


def heading_vector(psi_rad: float) -> np.ndarray:
    return np.array([np.cos(psi_rad), np.sin(psi_rad), 0.0])


def sheet_points(sheet: FieldAlignedSheet, f_hz: float) -> tuple[np.ndarray, np.ndarray]:
    """Points (N, 3) and weights.  lambda/2 along B, lambda/4 across."""
    lam = C_LIGHT / float(f_hz)
    fr = sheet_frame(sheet)
    na = max(64, int(np.ceil(sheet.length_m / (lam / 2.0))) | 1)
    nc = max(16, int(np.ceil(sheet.width_m / (lam / 4.0))) | 1)
    a = np.linspace(-sheet.length_m / 2.0, sheet.length_m / 2.0, na)
    c = np.linspace(-sheet.width_m / 2.0, sheet.width_m / 2.0, nc)
    line, _ = _centre_line(sheet, a, fr)
    pts = line[:, None, :] + c[None, :, None] * fr["u"][None, None, :]
    # Tukey taper along B so the ends do not diffract into the window.
    x = (a + sheet.length_m / 2.0) / sheet.length_m
    tf = float(sheet.tukey_fraction)
    wa = np.ones(na)
    lo, hi = x < tf / 2.0, x > 1.0 - tf / 2.0
    wa[lo] = 0.5 * (1.0 - np.cos(2.0 * np.pi * x[lo] / tf))
    wa[hi] = 0.5 * (1.0 - np.cos(2.0 * np.pi * (1.0 - x[hi]) / tf))
    w = np.repeat(wa[:, None], nc, axis=1).ravel()
    return pts.reshape(-1, 3), w / w.sum()


@dataclass(frozen=True)
class SheetEcho:
    t_s: np.ndarray
    field: np.ndarray
    group_path_m: np.ndarray
    mu_sat: float


def echo_free_space(
    t_s: np.ndarray, f_hz: float, sheet: FieldAlignedSheet, *, height_m: float, heading_rad: float,
    v_s_mps: float = 7600.0,
) -> SheetEcho:
    t = np.asarray(t_s, dtype=float)
    k = 2.0 * np.pi * float(f_hz) / C_LIGHT
    pts, w = sheet_points(sheet, f_hz)
    vel = sheet_velocity(sheet)
    p0 = platform_start(sheet, height_m)
    th = heading_vector(heading_rad)
    field = np.empty(t.size, dtype=complex)
    group = np.empty(t.size)
    for j, tj in enumerate(t):
        rel = (p0 + v_s_mps * tj * th) - (pts + vel * tj)
        r = np.sqrt(np.einsum("ij,ij->i", rel, rel))
        field[j] = np.sum(w * np.exp(-2j * k * r))
        sp = specular_point(sheet, p0 + v_s_mps * tj * th - vel * tj)
        group[j] = sp["range_m"]
    return SheetEcho(t_s=t, field=field, group_path_m=group, mu_sat=1.0)


def specular_point(sheet: FieldAlignedSheet, platform: np.ndarray) -> dict:
    """Free-space specular point on the sheet (infinite across B) for a platform.

    Solves ``(P - X(a, 0)) . t_a(a) = 0`` for ``a``, then takes ``c`` as the
    projection on ``u``.  Returns the point, its normal, the range and ``c``.
    """
    from scipy.optimize import brentq

    fr = sheet_frame(sheet)
    p = np.asarray(platform, dtype=float)

    def g(a: float) -> float:
        x, ta = _centre_line(sheet, np.array(a), fr)
        return float(np.dot(p - x, ta))

    lim = sheet.length_m
    a = brentq(g, -lim, lim, xtol=1e-9)
    x, ta = _centre_line(sheet, np.array(a), fr)
    c = float(np.dot(p - x, fr["u"]))
    xs = x + c * fr["u"]
    nrm = _normal(ta, fr["u"])
    return {"a": a, "c": c, "point": xs, "normal": nrm, "range_m": float(np.linalg.norm(p - xs))}


# ---------------------------------------------------------------------------
# Stratified plasma
# ---------------------------------------------------------------------------


def echo_stratified(
    t_s: np.ndarray, f_hz: float, sheet: FieldAlignedSheet, profile: ExponentialTopside, *,
    height_m: float, heading_rad: float, v_s_mps: float = 7600.0, n_nodes: int = 48,
    newton_iterations: int = 6,
) -> SheetEcho:
    """Received field inside a stratified plasma frozen to the moving sheet.

    In the frame of the medium the sheet is still and the platform moves at
    ``v_s t_hat - V``.  Depth is measured below the platform's height at
    ``t = 0``.
    """
    t = np.asarray(t_s, dtype=float)
    f = float(f_hz)
    k0 = 2.0 * np.pi * f / C_LIGHT
    pts, w = sheet_points(sheet, f)
    vel = sheet_velocity(sheet)
    p0 = platform_start_stratified(sheet, height_m, profile, f)
    th = heading_vector(heading_rad)
    nodes, weights = np.polynomial.legendre.leggauss(int(n_nodes))
    mu_sat = float(np.sqrt(profile.mu_squared(f, 0.0)))
    if np.any(profile.mu_squared(f, (p0[2] - pts[:, 2]).max() + 5e3) <= 0.05):
        raise ValueError("the sheet reaches the reflection level of this frequency")
    field = np.empty(t.size, dtype=complex)
    group = np.empty(t.size)
    q_prev = None
    for j, tj in enumerate(t):
        plat = p0 + (v_s_mps * th - vel) * tj  # medium frame
        d_top = float(p0[2] - plat[2])  # depth of the platform below its t=0 height
        d_pts = p0[2] - pts[:, 2]
        rho = np.hypot(plat[0] - pts[:, 0], plat[1] - pts[:, 1])
        iters = newton_iterations if q_prev is None else max(3, newton_iterations // 2)
        phase, pg, q_prev = _rays(rho, d_top, d_pts, f, profile, nodes, weights, iters, q_prev)
        field[j] = np.sum(w * np.exp(-2j * k0 * phase))
        group[j] = float(np.min(pg))
    return SheetEcho(t_s=t, field=field, group_path_m=group, mu_sat=mu_sat)


def _rays(rho, d_top, d_bot, f, profile, nodes, weights, iterations, q0):
    span = d_bot - d_top
    centre = 0.5 * (d_bot + d_top)
    d = centre[:, None] + 0.5 * span[:, None] * nodes[None, :]
    half = 0.5 * span[:, None] * weights[None, :]
    mu2 = profile.mu_squared(f, d)
    if q0 is None:
        q = rho / np.sum(half / np.sqrt(mu2), axis=1)
    else:
        q = q0
    cap = 0.95 * np.sqrt(mu2.min(axis=1))
    for _ in range(int(iterations)):
        root = np.sqrt(mu2 - (q * q)[:, None])
        x = np.sum(half * q[:, None] / root, axis=1)
        dx = np.sum(half * mu2 / root**3, axis=1)
        q = np.clip(q - (x - rho) / dx, 0.0, cap)
    root = np.sqrt(mu2 - (q * q)[:, None])
    p = np.sum(half * mu2 / root, axis=1)
    pg = np.sum(half / root, axis=1)
    return p, pg, q
