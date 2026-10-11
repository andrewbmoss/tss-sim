"""Layer: L1
WP: WP5, WP7
Claims: C3
Status: superseded for oblique links (see stratified_link)
Requirements: SCI-15, SCI-23, SCI-27

The virtual-mirror model of an oblique (bistatic) link between two sounders, and
of a vertical (monostatic) sounding, over a locally tilted, drifting layer on a
spherical Earth.

Why a mirror
------------
The Breit-Tuve theorem makes the group path of a ray reflected from a
stratified layer equal to the length of the straight two-legged path to a
virtual mirror, with the legs at the angle whose sine is the ray's horizontal
wavenumber over the free-space one.  So the measured group path of each link
fixes a virtual mirror height, and the horizontal part of the platform
velocities enters the Doppler exactly as it would off that mirror.  What the
mirror does not carry is the refractive index at the satellites, which scales
the vertical part: a rigid vertical drift ``w`` of the medium changes the phase
path at each end by ``w sqrt(mu^2 - q^2)``, where ``q`` is the sine of the
mirror leg's zenith angle, rather than by the free-space ``w sqrt(1 - q^2)``.
:func:`link_state` applies that factor.

Superseded for oblique links
----------------------------
Checked against the independent Haselgrove ray tracer of the Chartier
reproduction, the mirror gets the horizontal wavenumber of a curved-Earth
oblique ray about 2% wrong (a 1000 km link at 10 MHz: -4528 against -4636
m/s of line-of-sight rate), and an estimator built on it read a level, still
ionosphere as a 53 m/s drift and a 0.47 degree tilt.  The network study uses
:mod:`tss.l1_propagation.stratified_link` instead, which traces the ray through
the profile.  This module is kept for the geometry helpers it shares and as
the record of why the mirror was dropped.

The surface
-----------
A patch is centred on a point ``c`` of the sphere.  Points near it are written
in gnomonic tangent coordinates ``(x, y)``, east and north, so that the
direction of a surface point is ``normalize(R c + x e + y n)``.  The mirror
height is

    h(x, y, t) = h0 + gE x + gN y + w t,

so ``(gE, gN)`` is the gradient of the reflecting surface's height, in m/m,
and ``w`` its vertical drift.  The tilt of the surface normal is
``atan(|g|)`` toward ``-g``.  The same definition is used by the ray-traced
truth model of the two-node network study, so the two can be compared without
a change of variable.

Conventions
-----------
Positions and velocities are Earth-fixed Cartesian (ECEF) on a sphere of radius
6371 km.  The ionosphere co-rotates, so Earth-fixed velocities are velocities
relative to the medium.  The line-of-sight rate ``r`` is minus the rate of
change of the phase path, positive when the path is closing, as for the
two-node estimator of :mod:`tss.l3_signal.novel.second_look_estimator`.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np

__all__ = [
    "R_EARTH_M",
    "Patch",
    "LinkGeometry",
    "llh_to_ecef",
    "ecef_to_llh",
    "enu_basis",
    "specular_point",
    "link_state",
]

R_EARTH_M = 6371.0e3


def llh_to_ecef(lat_deg, lon_deg, alt_m) -> np.ndarray:
    """Spherical-Earth geodetic to ECEF, shape (..., 3)."""
    lat = np.radians(np.asarray(lat_deg, dtype=float))
    lon = np.radians(np.asarray(lon_deg, dtype=float))
    r = R_EARTH_M + np.asarray(alt_m, dtype=float)
    return np.stack([r * np.cos(lat) * np.cos(lon), r * np.cos(lat) * np.sin(lon), r * np.sin(lat)], axis=-1)


def ecef_to_llh(x: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    x = np.asarray(x, dtype=float)
    r = np.linalg.norm(x, axis=-1)
    lat = np.degrees(np.arcsin(x[..., 2] / r))
    lon = np.degrees(np.arctan2(x[..., 1], x[..., 0]))
    return lat, lon, r - R_EARTH_M


def enu_basis(lat_deg: float, lon_deg: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Unit east, north and up at a point."""
    la, lo = np.radians(lat_deg), np.radians(lon_deg)
    e = np.array([-np.sin(lo), np.cos(lo), 0.0])
    n = np.array([-np.sin(la) * np.cos(lo), -np.sin(la) * np.sin(lo), np.cos(la)])
    u = np.array([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])
    return e, n, u


@dataclass(frozen=True)
class Patch:
    """A locally tilted, drifting mirror.  Heights in m, gradients in m/m."""

    lat_deg: float
    lon_deg: float
    h0_m: float = 300e3
    g_east: float = 0.0
    g_north: float = 0.0
    w_mps: float = 0.0

    def basis(self):
        return enu_basis(self.lat_deg, self.lon_deg)

    def height(self, x: float, y: float, t: float) -> float:
        return self.h0_m + self.g_east * x + self.g_north * y + self.w_mps * t

    def point(self, x: float, y: float, t: float) -> np.ndarray:
        e, n, u = self.basis()
        d = R_EARTH_M * u + x * e + y * n
        d = d / np.linalg.norm(d)
        return (R_EARTH_M + self.height(x, y, t)) * d

    def tangent_coords(self, p: np.ndarray) -> tuple[float, float]:
        """Gnomonic (x, y) of the direction of ``p``."""
        e, n, u = self.basis()
        d = np.asarray(p, dtype=float)
        s = R_EARTH_M / float(np.dot(d, u))
        return float(np.dot(d, e) * s), float(np.dot(d, n) * s)


def specular_point(patch: Patch, a: np.ndarray, b: np.ndarray, t: float = 0.0, *,
                   iterations: int = 40) -> tuple[np.ndarray, float, float, float]:
    """The mirror point that makes ``|a - p| + |p - b|`` stationary.

    Returns the point, the path length and its tangent coordinates.  Newton's
    method on the two tangent coordinates with a finite-difference Hessian; the
    length is a smooth function of them, so it converges in a few steps from the
    great-circle midpoint.
    """
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)

    def length(x, y):
        p = patch.point(x, y, t)
        return float(np.linalg.norm(a - p) + np.linalg.norm(p - b))

    mid = 0.5 * (a + b)
    x, y = patch.tangent_coords(mid)
    h = 50.0
    for _ in range(int(iterations)):
        f0 = length(x, y)
        fxp, fxm = length(x + h, y), length(x - h, y)
        fyp, fym = length(x, y + h), length(x, y - h)
        fpp, fpm = length(x + h, y + h), length(x + h, y - h)
        fmp, fmm = length(x - h, y + h), length(x - h, y - h)
        gx, gy = (fxp - fxm) / (2 * h), (fyp - fym) / (2 * h)
        hxx = (fxp - 2 * f0 + fxm) / h**2
        hyy = (fyp - 2 * f0 + fym) / h**2
        hxy = (fpp - fpm - fmp + fmm) / (4 * h * h)
        hess = np.array([[hxx, hxy], [hxy, hyy]])
        step = np.linalg.solve(hess, np.array([gx, gy]))
        x -= float(step[0])
        y -= float(step[1])
        if np.hypot(step[0], step[1]) < 1e-4:
            break
    p = patch.point(x, y, t)
    return p, length(x, y), x, y


@dataclass(frozen=True)
class LinkGeometry:
    """Two ends of a link: positions and Earth-fixed velocities, ECEF, SI.

    ``mu_a`` and ``mu_b`` are the refractive indices at the two satellites at
    the link frequency, from each satellite's own measured local plasma
    frequency.  ``a`` and ``b`` coincide for a vertical sounding.
    """

    a: np.ndarray
    b: np.ndarray
    va: np.ndarray
    vb: np.ndarray
    mu_a: float = 1.0
    mu_b: float = 1.0


@dataclass(frozen=True)
class LinkState:
    group_path_m: float
    los_rate_mps: float
    point: np.ndarray
    sin_zenith_a: float
    sin_zenith_b: float
    drift_coefficient: float


def _zenith_sines(p: np.ndarray, a: np.ndarray, b: np.ndarray) -> tuple[float, float]:
    out = []
    for s in (a, b):
        leg = p - s
        leg = leg / np.linalg.norm(leg)
        up = s / np.linalg.norm(s)
        out.append(float(np.sqrt(max(0.0, 1.0 - np.dot(leg, -up) ** 2))))
    return out[0], out[1]


def link_state(patch: Patch, link: LinkGeometry, t: float = 0.0, *, dt_s: float = 0.05) -> LinkState:
    """Group path and line-of-sight rate of a link off the patch's mirror.

    The rate is the central difference of the mirror path over ``2 dt_s`` with
    the satellites moving and the mirror drifting, plus the refractive-index
    correction to the drift term described in the module docstring.
    """
    a0, b0 = np.asarray(link.a, float), np.asarray(link.b, float)
    va, vb = np.asarray(link.va, float), np.asarray(link.vb, float)
    p, length, _, _ = specular_point(patch, a0, b0, t)
    _, lp, _, _ = specular_point(patch, a0 + va * dt_s, b0 + vb * dt_s, t + dt_s)
    _, lm, _, _ = specular_point(patch, a0 - va * dt_s, b0 - vb * dt_s, t - dt_s)
    r_free = -(lp - lm) / (2.0 * dt_s)
    # The free-space drift coefficient, from the same difference with the
    # satellites held still.
    _, wp, _, _ = specular_point(patch, a0, b0, t + dt_s)
    _, wm, _, _ = specular_point(patch, a0, b0, t - dt_s)
    if patch.w_mps != 0.0:
        c_free = -(wp - wm) / (2.0 * dt_s) / patch.w_mps
    else:
        p2 = replace(patch, w_mps=1.0)
        _, wp, _, _ = specular_point(p2, a0, b0, t + dt_s)
        _, wm, _, _ = specular_point(p2, a0, b0, t - dt_s)
        c_free = -(wp - wm) / (2.0 * dt_s)
    qa, qb = _zenith_sines(p, a0, b0)
    fa = np.sqrt(max(link.mu_a**2 - qa * qa, 0.0)) / np.sqrt(max(1.0 - qa * qa, 1e-12))
    fb = np.sqrt(max(link.mu_b**2 - qb * qb, 0.0)) / np.sqrt(max(1.0 - qb * qb, 1e-12))
    f_mu = 0.5 * (fa + fb)
    c = c_free * f_mu
    r = r_free + (f_mu - 1.0) * c_free * patch.w_mps
    return LinkState(group_path_m=length, los_rate_mps=r, point=p, sin_zenith_a=qa, sin_zenith_b=qb,
                     drift_coefficient=c)
