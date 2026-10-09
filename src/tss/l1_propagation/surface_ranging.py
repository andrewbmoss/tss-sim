"""Layer: L1
WP: WP5
Claims: C3
Status: implemented
Requirements: SCI-15, SCI-23

Exact free-space ranging to a smooth drifting surface in three dimensions.

This is the three-dimensional forward model of the two-node study
(:mod:`tss.studies.route_two`).  It is written so that it shares no equation
with the second-look estimator it tests
(:mod:`tss.l3_signal.novel.second_look_estimator`).  The estimator assumes the
surface is a plane over the footprint and that the range and line-of-sight rate
follow from the plane's normal and vertical motion.  This module assumes
neither.  It finds the specular point by minimising the distance from the
platform to the surface numerically, and takes the line-of-sight rate as the
finite difference of that distance in time, so a curved or evolving surface is
ranged exactly.

The surface
-----------
In a local frame with ``x`` east, ``y`` north and ``z`` up, at the platform's
altitude, the surface height is

    z_s(x, y, t) = z0 + a x + b y + w t + A sin(omega t - Kx x - Ky y + phi)

which is a plane with gradient ``(a, b)`` drifting vertically at ``w``, carrying
a travelling ionospheric disturbance of vertical amplitude ``A``.  ``z0`` is
negative: the surface lies below the platform.  The disturbance makes the
surface's tilt and vertical velocity at a given point change between two
passes, which is the evolution the second look has to survive.

What is left out
----------------
The plasma.  The two-dimensional physical-optics model
(:mod:`tss.l1_propagation.plate_echo`) carries the refractive index at the
spacecraft and at the reflector, and the two-node study verifies the plasma
form of the estimator with it.  This module carries what that one cannot, which
is the third dimension: the cross-track tilt, the track offset that Earth's
rotation introduces between two nodes, and tracks that cross.  Earth's
curvature is also left out.  Over the few hundred kilometres between two
specular points it changes the vertical by less than the platform's own
altitude knowledge would, and both nodes see it alike.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = ["DriftingSurface", "SurfaceSample", "specular_point", "sample_surface"]


@dataclass(frozen=True)
class DriftingSurface:
    """A tilted plane drifting vertically, with an optional travelling wave.

    Lengths in m, times in s, angles in rad.  ``a`` and ``b`` are the slopes
    ``dz/dx`` and ``dz/dy`` of the mean plane, not angles.  The upward unit
    normal of the mean plane is ``(-a, -b, 1) / sqrt(1 + a^2 + b^2)``.
    """

    z0_m: float = -500e3
    a: float = 0.0
    b: float = 0.0
    w_mps: float = 0.0
    amp_m: float = 0.0
    omega_rad_s: float = 0.0
    kx_rad_m: float = 0.0
    ky_rad_m: float = 0.0
    phase_rad: float = 0.0

    def height(self, x, y, t):
        arg = self.omega_rad_s * t - self.kx_rad_m * x - self.ky_rad_m * y + self.phase_rad
        return self.z0_m + self.a * x + self.b * y + self.w_mps * t + self.amp_m * np.sin(arg)

    def gradient(self, x, y, t) -> tuple[float, float]:
        arg = self.omega_rad_s * t - self.kx_rad_m * x - self.ky_rad_m * y + self.phase_rad
        c = self.amp_m * np.cos(arg)
        return self.a - self.kx_rad_m * c, self.b - self.ky_rad_m * c

    def hessian(self, x, y, t) -> tuple[float, float, float]:
        arg = self.omega_rad_s * t - self.kx_rad_m * x - self.ky_rad_m * y + self.phase_rad
        s = -self.amp_m * np.sin(arg)
        return self.kx_rad_m**2 * s, self.kx_rad_m * self.ky_rad_m * s, self.ky_rad_m**2 * s

    def vertical_velocity(self, x, y, t) -> float:
        """``dz_s/dt`` at a fixed point: the drift the second look measures."""
        arg = self.omega_rad_s * t - self.kx_rad_m * x - self.ky_rad_m * y + self.phase_rad
        return float(self.w_mps + self.amp_m * self.omega_rad_s * np.cos(arg))

    def unit_normal(self, x, y, t) -> np.ndarray:
        gx, gy = self.gradient(x, y, t)
        n = np.array([-gx, -gy, 1.0])
        return n / np.linalg.norm(n)


def specular_point(surface: DriftingSurface, p: np.ndarray, t: float, *,
                   iterations: int = 30, tol_m: float = 1e-7) -> tuple[np.ndarray, float]:
    """The surface point nearest the platform at ``p``, and the distance to it.

    Newton's method on the squared distance, started from the foot of the
    perpendicular to the mean plane.  The surface is smooth on the scale of the
    footprint, so the minimum is unique and Newton converges in a few steps.
    """
    p = np.asarray(p, dtype=float)
    # Foot of the perpendicular to the mean plane, ignoring the wave.
    n0 = np.array([-surface.a, -surface.b, 1.0])
    n0 = n0 / np.linalg.norm(n0)
    plane_at_p = surface.z0_m + surface.a * p[0] + surface.b * p[1] + surface.w_mps * t
    d0 = (p[2] - plane_at_p) * n0[2]
    q = p - d0 * n0
    x, y = float(q[0]), float(q[1])
    for _ in range(int(iterations)):
        z = float(surface.height(x, y, t))
        gx, gy = surface.gradient(x, y, t)
        hxx, hxy, hyy = surface.hessian(x, y, t)
        dz = z - p[2]
        g = np.array([2.0 * (x - p[0]) + 2.0 * dz * gx, 2.0 * (y - p[1]) + 2.0 * dz * gy])
        h = np.array([
            [2.0 + 2.0 * gx * gx + 2.0 * dz * hxx, 2.0 * gx * gy + 2.0 * dz * hxy],
            [2.0 * gx * gy + 2.0 * dz * hxy, 2.0 + 2.0 * gy * gy + 2.0 * dz * hyy],
        ])
        step = np.linalg.solve(h, g)
        x -= float(step[0])
        y -= float(step[1])
        if np.hypot(step[0], step[1]) < tol_m:
            break
    foot = np.array([x, y, float(surface.height(x, y, t))])
    return foot, float(np.linalg.norm(p - foot))


@dataclass(frozen=True)
class SurfaceSample:
    """What a platform measures from the surface at one instant, in free space.

    ``range_m`` is the distance to the specular point and ``los_rate_mps`` its
    rate of decrease, positive when closing, so it has the sign of the Doppler.
    ``foot`` is the specular point, where the truth is evaluated.
    """

    t_s: float
    range_m: float
    los_rate_mps: float
    foot: np.ndarray


def sample_surface(surface: DriftingSurface, p0: np.ndarray, v: np.ndarray, t0: float,
                   *, dt_s: float = 0.05) -> SurfaceSample:
    """Range and line-of-sight rate for a platform at ``p0`` at time ``t0``.

    The platform moves at the constant velocity ``v`` (local frame, relative to
    the medium).  The rate is a central difference of the exact range over
    ``2 dt_s``; with the range exact to well below a millimetre the difference
    is accurate to about a millimetre per second.
    """
    p0 = np.asarray(p0, dtype=float)
    v = np.asarray(v, dtype=float)
    if v.size == 2:
        v = np.array([v[0], v[1], 0.0])
    foot, r0 = specular_point(surface, p0, t0)
    _, rp = specular_point(surface, p0 + v * dt_s, t0 + dt_s)
    _, rm = specular_point(surface, p0 - v * dt_s, t0 - dt_s)
    return SurfaceSample(t_s=float(t0), range_m=r0, los_rate_mps=-(rp - rm) / (2.0 * dt_s), foot=foot)
