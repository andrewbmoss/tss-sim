"""Layer: L1
WP: WP5
Claims: C3
Status: partial
Requirements: SCI-15, SCI-23

Rays from a platform to many scattering points in a spherically stratified,
field-free plasma: the three-dimensional curved-Earth geometry of the facet
tier (WP5 Phase 3, Andrew 2026-11-25, decisions (a) and (b)).

WHY ONE PLANE PER RAY
---------------------
In a medium whose index depends on the radius alone, Bouguer's law
``r mu sin(theta) = const`` holds along every ray and the ray stays in the plane
through the Earth's centre, its start and its end.  So a ray from the platform
``P`` to a point ``X`` is a two-dimensional problem in the plane ``(O, P, X)``,
whatever the three-dimensional layout of the points.  In that plane the
Earth-flattening map about a reference radius ``R``

    x = R psi,    d = R ln(R / r),    mu'^2(d) = mu^2(r) (r / R)^2,

turns the sphere into a flat stratified medium in which Snell's invariant
``q = mu' sin(theta')`` is Bouguer's constant over ``R``.  The map is conformal,
so the ray, its angular span ``psi`` and its phase path ``int mu ds`` are those
of the sphere exactly (Kerr 1951).  The group path of a field-free plasma
(``mu_g = 1 / mu``) is ``int ds / mu``, which in the flat frame is

    P'(q) = int (r / R)^2 / sqrt(mu'^2 - q^2) dd,

because ``ds = (r / R) ds'`` and ``mu = mu' R / r``.  The phase path is
``int mu'^2 / sqrt(mu'^2 - q^2) dd`` and the horizontal reach
``int q / sqrt(mu'^2 - q^2) dd = R psi`` fixes ``q`` by Newton's method, as in
:mod:`tss.l1_propagation.plate_echo`.

FRAME AND PRECISION
-------------------
Points are given in a local Cartesian frame (east, north, up) whose origin is
on the sphere of radius ``r0`` and whose Earth centre is ``(0, 0, -r0)``.
Radii and angles are formed from local offsets, never by subtracting two
geocentric vectors, so the frame keeps full precision for any ``r0``; a very
large ``r0`` is the flat limit, which is how the module is held to the flat
plate and sheet models (``tests/unit/test_l1_sphere.py``).

THE CHECK THAT DOES NOT USE THE MAP
-----------------------------------
:func:`bouguer_paths` integrates the same ray directly in polar coordinates,

    psi = int c dr / (r sqrt(r^2 mu^2 - c^2)),
    P   = int r mu^2 dr / sqrt(r^2 mu^2 - c^2),
    P'  = int r dr / sqrt(r^2 mu^2 - c^2),

with Bouguer's constant ``c`` found by bisection.  It shares no step with the
flattening and is slow; it exists to test :func:`ray_paths`.

MOTION
------
Drift is the strata moving as a whole: the profile rises radially by ``w t``
(``mu^2(r) -> mu^2(r - w t)``) and the scatterers move rigidly with their own
velocity.  This is :class:`SphericalMedium`'s ``rise_m``.

What is left out: the magnetic field's effect on the wave, absorption,
horizontal gradients.  A ray that turns (a point beyond the ray's reach)
is not handled; points must lie below the platform within the ray's reach,
which holds for the near-nadir geometry of the drift mode.
"""

from __future__ import annotations

import math
from collections.abc import Callable
from dataclasses import dataclass
from typing import Any

import numpy as np

__all__ = [
    "LocalFrame",
    "SphericalMedium",
    "radius_and_angle",
    "ray_paths",
    "bouguer_paths",
]


@dataclass(frozen=True)
class LocalFrame:
    """East-north-up frame on the sphere of radius ``r0`` at (``lat_deg``, ``lon_deg``).

    Spherical (geocentric) latitude, as :mod:`tss.platform.orbits` reports it.
    """

    lat_deg: float
    lon_deg: float
    r0_m: float

    def axes(self) -> np.ndarray:
        """Rows: east, north, up, as Earth-fixed unit vectors."""
        la, lo = math.radians(self.lat_deg), math.radians(self.lon_deg)
        sl, cl, so, co = math.sin(la), math.cos(la), math.sin(lo), math.cos(lo)
        return np.array([[-so, co, 0.0], [-sl * co, -sl * so, cl], [cl * co, cl * so, sl]])

    def origin_ecef(self) -> np.ndarray:
        return self.r0_m * self.axes()[2]

    def from_llh(self, lat_deg: float, lon_deg: float, r_m: float) -> np.ndarray:
        """Local position of the point at spherical latitude, longitude and radius."""
        la, lo = math.radians(lat_deg), math.radians(lon_deg)
        u = np.array([math.cos(la) * math.cos(lo), math.cos(la) * math.sin(lo), math.sin(la)])
        ax = self.axes()
        # r u - r0 up, formed as (r - r0) u + r0 (u - up) to keep precision.
        up = ax[2]
        diff = (float(r_m) - self.r0_m) * u + self.r0_m * (u - up)
        return ax @ diff

    def vector_from_ned(self, v_ned: np.ndarray, lat_deg: float, lon_deg: float) -> np.ndarray:
        """A north-east-down vector at another point, in this frame's axes."""
        la, lo = math.radians(lat_deg), math.radians(lon_deg)
        sl, cl, so, co = math.sin(la), math.cos(la), math.sin(lo), math.cos(lo)
        north = np.array([-sl * co, -sl * so, cl])
        east = np.array([-so, co, 0.0])
        down = np.array([-cl * co, -cl * so, -sl])
        v = np.asarray(v_ned, dtype=float)
        return self.axes() @ (v[0] * north + v[1] * east + v[2] * down)


def radius_and_angle(p: np.ndarray, x: np.ndarray, r0_m: float) -> tuple[float, np.ndarray, np.ndarray]:
    """Heights above ``r0`` of ``p`` and of the points ``x`` (N, 3), and the angles at the centre.

    Local coordinates, centre at ``(0, 0, -r0)``.  The radius is ``r0 + delta``
    with ``delta = (2 r0 z + |x|^2) / (r0 + r)``, returned as ``delta`` so that
    nothing is lost to the size of ``r0``; ``psi = atan2(|P x X|, P . X)`` with
    ``P x X = p x x + r0 z_hat x (x - p)``.
    """
    p = np.asarray(p, dtype=float)
    x = np.atleast_2d(np.asarray(x, dtype=float))
    r0 = float(r0_m)

    def _height(v: np.ndarray) -> np.ndarray:
        s2 = np.sum(v * v, axis=-1)
        rough = np.sqrt(s2 + 2.0 * r0 * v[..., 2] + r0 * r0)
        return (2.0 * r0 * v[..., 2] + s2) / (r0 + rough)

    rp = float(_height(p))
    rx = _height(x)
    d = x - p
    cross = np.cross(p, x) + r0 * np.stack([-d[:, 1], d[:, 0], np.zeros(len(x))], axis=1)
    dot = x @ p + r0 * (p[2] + x[:, 2]) + r0 * r0
    psi = np.arctan2(np.linalg.norm(cross, axis=1), dot)
    return rp, rx, psi


@dataclass(frozen=True)
class SphericalMedium:
    """A field-free index against radius, from a profile against depth.

    ``base.mu_squared(f, depth)`` is the index at ``depth = r_base - r``; it is
    the exponential topside or the tabulated L0 profile of the facet tier.
    ``base`` None is free space.  ``rise_m`` lifts the whole profile.  Radii
    are given as heights above the frame radius ``r0_m`` (see
    :func:`radius_and_angle`), so ``r_base_m`` is too.
    """

    base: Any
    r_base_m: float
    rise_m: float = 0.0

    def mu_squared(self, f_hz: float, h_m: np.ndarray) -> np.ndarray:
        """Index squared at height ``h_m`` above the frame radius."""
        h = np.asarray(h_m, dtype=float)
        if self.base is None:
            return np.ones_like(h)
        return self.base.mu_squared(f_hz, self.r_base_m - h + self.rise_m)

    def shifted(self, rise_m: float) -> SphericalMedium:
        return SphericalMedium(self.base, self.r_base_m, self.rise_m + float(rise_m))


def ray_paths(
    p: np.ndarray,
    x: np.ndarray,
    f_hz: float,
    medium: SphericalMedium,
    r0_m: float,
    *,
    n_nodes: int = 48,
    newton_iterations: int = 8,
    q0: np.ndarray | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Phase path, group path and Snell invariant of the rays from ``p`` to each ``x``.

    Flattened about ``R = r_p``, the platform's own radius, so its flat depth
    is zero.  Free space (``medium.base`` None) returns the chord for both
    paths.  ``q0`` warm-starts Newton's method (the previous sweep's ``q``).
    """
    hp, hx, psi = radius_and_angle(p, x, r0_m)
    if medium.base is None:
        chord = np.linalg.norm(np.atleast_2d(x) - np.asarray(p, dtype=float), axis=1)
        return chord, chord.copy(), np.zeros_like(chord)
    return _flat_rays(hp, hx, psi, f_hz, medium, n_nodes, newton_iterations, q0, r0_m=r0_m)


def _flat_rays(hp, hx, psi, f_hz, medium, n_nodes, iterations, q0, *, r0_m):
    """The flattened rays; heights ``hp``, ``hx`` above ``r0_m`` (module docstring)."""
    big_r = float(r0_m) + float(hp)
    # Flat depths; d = -R log1p((r - R) / R), with r - R = hx - hp exactly.
    d_bot = -big_r * np.log1p((hx - float(hp)) / big_r)
    if np.any(d_bot <= 0.0):
        raise ValueError("a point is not below the platform")
    offset = big_r * psi
    nodes, weights = np.polynomial.legendre.leggauss(int(n_nodes))
    d = 0.5 * d_bot[:, None] * (1.0 + nodes[None, :])
    half = 0.5 * d_bot[:, None] * weights[None, :]
    ratio = np.exp(-d / big_r)  # r / R
    # Height above r0 of each node: hp + R (ratio - 1), with expm1 for precision.
    h = float(hp) + big_r * np.expm1(-d / big_r)
    mu2 = medium.mu_squared(f_hz, h) * ratio * ratio
    if np.any(mu2 <= 0.0):
        raise ValueError("a ray meets the reflection level before its point")
    if q0 is None or np.shape(q0) != offset.shape:
        q = offset / np.sum(half / np.sqrt(mu2), axis=1)
    else:
        q = np.asarray(q0, dtype=float).copy()
    cap = 0.95 * np.sqrt(mu2.min(axis=1))
    for _ in range(int(iterations)):
        root = np.sqrt(mu2 - (q * q)[:, None])
        xr = np.sum(half * q[:, None] / root, axis=1)
        dx = np.sum(half * mu2 / root**3, axis=1)
        q = np.clip(q - (xr - offset) / dx, 0.0, cap)
    root = np.sqrt(mu2 - (q * q)[:, None])
    phase = np.sum(half * mu2 / root, axis=1)
    group = np.sum(half * ratio * ratio / root, axis=1)
    return phase, group, q


def bouguer_paths(
    rp: float,
    rx: float,
    psi: float,
    f_hz: float,
    mu_squared_of_r: Callable[[float, np.ndarray], np.ndarray],
    *,
    n_nodes: int = 400,
) -> tuple[float, float]:
    """Phase and group path of one ray by Bouguer's law in polar coordinates.

    Independent of the flattening (module docstring).  ``rp`` and ``rx`` are
    true radii; ``mu_squared_of_r(f, r)`` the index at radius ``r``.  The
    integrals in ``r`` run from ``rx`` up to ``rp``; ``c`` is found by
    bisection on ``psi``.
    """
    nodes, weights = np.polynomial.legendre.leggauss(int(n_nodes))
    r = 0.5 * (rp + rx) + 0.5 * (rp - rx) * nodes
    w = 0.5 * (rp - rx) * weights
    m2 = mu_squared_of_r(f_hz, r)
    if np.any(m2 <= 0.0):
        raise ValueError("the ray meets the reflection level")
    a2 = r * r * m2

    def reach(c: float) -> float:
        return float(np.sum(w * c / (r * np.sqrt(a2 - c * c))))

    lo, hi = 0.0, math.sqrt(float(a2.min())) * (1.0 - 1e-12)
    if psi <= 0.0:
        c = 0.0
    else:
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if reach(mid) < psi:
                lo = mid
            else:
                hi = mid
        c = 0.5 * (lo + hi)
    root = np.sqrt(a2 - c * c)
    return float(np.sum(w * r * m2 / root)), float(np.sum(w * r / root))
