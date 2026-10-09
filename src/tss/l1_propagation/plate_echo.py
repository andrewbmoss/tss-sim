"""Layer: L1
WP: WP5
Claims: C3
Status: implemented
Requirements: SCI-15, SCI-23

Physical-optics echo of a flat reflecting plate seen from a moving platform,
in free space and inside a horizontally stratified plasma.

This is the forward model of the C3 verification study
(:mod:`tss.studies.c3_verify`, Technical Volume Section 8.15).  It is written
so that it shares no equation with the route-one estimator it is used to test.
The estimator models a linearised Doppler history and a ``sinc^2`` power
envelope.  This module models neither.  It sums the two-way phase of many
points along a plate and returns the complex field the platform would receive
on each sweep, and the estimator has to find the Doppler history and the power
envelope in that field for itself.

The plate
---------
A plate is a straight segment in the along-track and vertical plane.  Its
centre is at along-track position ``x_m`` and depth ``depth_m`` below the
platform at ``t = 0``.  Its normal is tilted by ``tilt_rad`` toward ``+x``, so
the unit normal is ``(sin a, cos a)`` in ``(x, up)``.  A compact facet is a
short plate.  An extended layer is a long one with tapered ends, so that its
ends do not diffract into the window.

Drift is a rigid vertical motion of the plate and, inside the plasma, of the
whole profile with it, positive upward toward the platform.  That is the
"strata moving as a whole" case of Technical Volume Section 8.2.  The plate
may also move along track, at ``drift_x_mps``.  The profile is horizontally
stratified, so moving it along track changes nothing, and only the plate
moves.

Free space
----------
The field on a sweep at time ``t`` is

    E(t) = sum_i w_i exp(-j 2 k R_i(t)),

with ``R_i`` the exact distance from the platform at ``(v_s t, 0)`` to point
``i`` and ``k = 2 pi f / c``.  Amplitude factors that vary slowly across the
plate (spreading, obliquity) are left out, so the result is the scalar
Kirchhoff sum for a plate of uniform reflectivity.

Inside the plasma
-----------------
The profile is field free and horizontally stratified with an exponential
topside, ``f_N^2(d) = f_N,sat^2 exp(d / H)`` at depth ``d`` below the platform's
height at ``t = 0``.  For each point the ray from the platform is found from
Snell's law, with ``q = mu sin(theta)`` held constant along it, by solving

    X(q) = int q / sqrt(mu^2 - q^2) dd = x_platform - x_i

for ``q`` with Newton's method.  The phase path and group path of that ray are

    P(q)  = int mu^2 / sqrt(mu^2 - q^2) dd,
    P'(q) = int 1 / sqrt(mu^2 - q^2) dd,

both taken from the platform's depth in the moving medium to the point's
depth.  The integrals use fixed Gauss-Legendre nodes, which is exact to
rounding for an integrand this smooth.  The field is the same sum with
``k_0 P_i`` in place of ``k R_i``.  Nothing here assumes the ray is short, that
the angle is small, or that the history is linear in time.

What is left out
----------------
The magnetic field, absorption, any variation of reflectivity across the
plate, and the spreading loss.  None of them moves the time of the power peak
or the zero crossing of the Doppler history for a plate a few kilometres long,
which are the two quantities the verification study scores.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from tss.constants import C_LIGHT

__all__ = [
    "Plate",
    "ExponentialTopside",
    "PlateEcho",
    "plate_points",
    "echo_free_space",
    "echo_stratified",
]


@dataclass(frozen=True)
class Plate:
    """A flat reflecting plate.  Lengths in m, angles in rad, speeds in m/s.

    Attributes
    ----------
    x_m
        Along-track position of the centre at ``t = 0``.
    depth_m
        Depth of the centre below the platform at ``t = 0``.
    tilt_rad
        Tilt of the normal from vertical, positive toward ``+x``.
    length_m
        Length along the plate.
    drift_mps
        Rigid vertical velocity, positive upward toward the platform.
    drift_x_mps
        Along-track velocity of the plate, positive toward ``+x``.  Added for
        the two-node study, where a facet's horizontal drift between passes is
        what SCI-23 scores.  Zero by default, which leaves every earlier result
        unchanged.
    taper
        ``"uniform"`` for a facet, whose pattern is then a ``sinc``, or
        ``"tukey"`` for an extended layer, whose ends must not diffract.
    tukey_fraction
        Fraction of the length given to the cosine ends when ``taper`` is
        ``"tukey"``.
    """

    x_m: float = 0.0
    depth_m: float = 500e3
    tilt_rad: float = 0.0
    length_m: float = 2.0e3
    drift_mps: float = 0.0
    taper: str = "uniform"
    tukey_fraction: float = 0.25
    drift_x_mps: float = 0.0


@dataclass(frozen=True)
class ExponentialTopside:
    """Field-free exponential topside, ``f_N^2 = f_N,sat^2 exp(d / H)``.

    ``d`` is depth below the platform's height at ``t = 0``, in the frame of
    the medium.
    """

    f_n_sat_hz: float = 1.5e6
    scale_height_m: float = 200e3

    def mu_squared(self, f_hz: float, depth_m: np.ndarray) -> np.ndarray:
        f_n2 = self.f_n_sat_hz**2 * np.exp(np.asarray(depth_m, dtype=float) / self.scale_height_m)
        return 1.0 - f_n2 / float(f_hz) ** 2

    def depth_for_mu(self, f_hz: float, mu: float) -> float:
        """Depth at which the field-free index falls to ``mu``."""
        f_n2 = (1.0 - float(mu) ** 2) * float(f_hz) ** 2
        return float(self.scale_height_m * np.log(f_n2 / self.f_n_sat_hz**2))


@dataclass(frozen=True)
class PlateEcho:
    """What the platform receives from a plate, one value per sweep.

    ``field`` is complex and unnormalised in phase convention
    ``exp(-j 2 k_0 P)``, so a closing path gives a positive Doppler shift.
    ``group_path_m`` is the group path of the ray to the plate's centre, which
    is the apparent range the sounder would measure.  ``mu_sat`` is the index
    at the platform (1 in free space).
    """

    t_s: np.ndarray
    field: np.ndarray
    group_path_m: np.ndarray
    mu_sat: float


def plate_points(plate: Plate, f_hz: float) -> tuple[np.ndarray, np.ndarray]:
    """Offsets along the plate and their weights.

    The spacing is a quarter of the vacuum wavelength, which keeps the two-way
    phase step between neighbouring points below ``pi/2`` at any angle and so
    keeps the sum a faithful quadrature of the Kirchhoff integral.
    """
    lam = C_LIGHT / float(f_hz)
    n = max(64, int(np.ceil(plate.length_m / (lam / 4.0))) | 1)
    s = np.linspace(-plate.length_m / 2.0, plate.length_m / 2.0, n)
    if plate.taper == "uniform":
        w = np.ones(n)
    elif plate.taper == "tukey":
        a = float(plate.tukey_fraction)
        u = (s + plate.length_m / 2.0) / plate.length_m
        w = np.ones(n)
        lo, hi = u < a / 2.0, u > 1.0 - a / 2.0
        w[lo] = 0.5 * (1.0 - np.cos(2.0 * np.pi * u[lo] / a))
        w[hi] = 0.5 * (1.0 - np.cos(2.0 * np.pi * (1.0 - u[hi]) / a))
    else:
        raise ValueError(f"unknown taper {plate.taper!r}")
    return s, w / w.sum()


def _point_coordinates(plate: Plate, s: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Along-track position and depth of each point at ``t = 0``.

    The tangent is ``(cos a, -sin a)`` in ``(x, up)``, perpendicular to the
    normal ``(sin a, cos a)``.  Depth is positive downward.
    """
    a = float(plate.tilt_rad)
    return plate.x_m + s * np.cos(a), plate.depth_m + s * np.sin(a)


def echo_free_space(
    t_s: np.ndarray,
    f_hz: float,
    plate: Plate,
    *,
    v_s_mps: float = 7600.0,
) -> PlateEcho:
    """Received field from a plate in free space, one complex value per sweep."""
    t = np.asarray(t_s, dtype=float)
    k = 2.0 * np.pi * float(f_hz) / C_LIGHT
    s, w = plate_points(plate, f_hz)
    xi, di = _point_coordinates(plate, s)
    field = np.empty(t.size, dtype=complex)
    group = np.empty(t.size, dtype=float)
    for j, tj in enumerate(t):
        dx = float(v_s_mps) * tj - xi - float(plate.drift_x_mps) * tj
        dz = di - float(plate.drift_mps) * tj
        r = np.hypot(dx, dz)
        field[j] = np.sum(w * np.exp(-2j * k * r))
        group[j] = float(np.hypot(
            float(v_s_mps) * tj - plate.x_m - plate.drift_x_mps * tj, plate.depth_m - plate.drift_mps * tj
        ))
    return PlateEcho(t_s=t, field=field, group_path_m=group, mu_sat=1.0)


def _ray_integrals(
    q: np.ndarray, mu2: np.ndarray, half: np.ndarray
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """X, dX/dq, P and P' for rays ``q`` through nodes ``mu2`` of weight ``half``.

    ``mu2`` has shape (points, nodes) and already carries the interval mapping
    in ``half``, which is the Gauss-Legendre weight times half the interval.
    """
    qq = q[:, None]
    root = np.sqrt(mu2 - qq * qq)
    x = np.sum(half * qq / root, axis=1)
    dx = np.sum(half * mu2 / root**3, axis=1)
    p = np.sum(half * mu2 / root, axis=1)
    pg = np.sum(half / root, axis=1)
    return x, dx, p, pg


def echo_stratified(
    t_s: np.ndarray,
    f_hz: float,
    plate: Plate,
    profile: ExponentialTopside,
    *,
    v_s_mps: float = 7600.0,
    n_nodes: int = 64,
    newton_iterations: int = 8,
) -> PlateEcho:
    """Received field from a plate inside a stratified plasma.

    Raises
    ------
    ValueError
        If the plate lies at or below the reflection level of the sounding
        frequency, where the ordinary ray does not reach it.
    """
    t = np.asarray(t_s, dtype=float)
    f = float(f_hz)
    k0 = 2.0 * np.pi * f / C_LIGHT
    s, w = plate_points(plate, f)
    xi, di = _point_coordinates(plate, s)
    if np.any(profile.mu_squared(f, di.max() + abs(plate.drift_mps) * float(np.ptp(t))) <= 0.05):
        raise ValueError("the plate reaches the reflection level of this frequency")
    nodes, weights = np.polynomial.legendre.leggauss(int(n_nodes))
    mu_sat = float(np.sqrt(profile.mu_squared(f, 0.0)))

    # The centre ray, for the apparent range.
    xc, dc = np.array([plate.x_m]), np.array([plate.depth_m])

    field = np.empty(t.size, dtype=complex)
    group = np.empty(t.size, dtype=float)
    for j, tj in enumerate(t):
        d_sat = float(plate.drift_mps) * tj  # the medium has risen past the platform
        u_t = float(plate.drift_x_mps) * tj
        phase, group_c = _solve(
            float(v_s_mps) * tj - xi - u_t, d_sat, di, f, profile, nodes, weights, newton_iterations
        )
        field[j] = np.sum(w * np.exp(-2j * k0 * phase))
        _, gc = _solve(
            np.array([float(v_s_mps) * tj - xc[0] - u_t]), d_sat, dc, f, profile, nodes, weights,
            newton_iterations,
        )
        group[j] = float(gc[0])
    return PlateEcho(t_s=t, field=field, group_path_m=group, mu_sat=mu_sat)


def _solve(
    offset_m: np.ndarray,
    d_top: float,
    d_bottom: np.ndarray,
    f_hz: float,
    profile: ExponentialTopside,
    nodes: np.ndarray,
    weights: np.ndarray,
    iterations: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Phase path and group path of the ray from depth ``d_top`` to each point."""
    span = d_bottom - d_top
    centre = 0.5 * (d_bottom + d_top)
    d = centre[:, None] + 0.5 * span[:, None] * nodes[None, :]
    half = 0.5 * span[:, None] * weights[None, :]
    mu2 = profile.mu_squared(f_hz, d)
    # Start from the paraxial ray, X = q P'.
    pg0 = np.sum(half / np.sqrt(mu2), axis=1)
    q = offset_m / pg0
    cap = 0.95 * np.sqrt(mu2.min(axis=1))
    for _ in range(int(iterations)):
        x, dx, _, _ = _ray_integrals(q, mu2, half)
        q = np.clip(q - (x - offset_m) / dx, -cap, cap)
    _, _, p, pg = _ray_integrals(q, mu2, half)
    return p, pg
