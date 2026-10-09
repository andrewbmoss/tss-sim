"""Layer: platform
WP: WP5, WP7
Claims: C3
Status: implemented
Requirements: MIS-10, SCI-15, SCI-23

The geometry of two nodes that look at the same point a known interval apart.

Route two (Technical Volume Section 8.11) needs a second node to pass over the
point the first one sounded.  How nearly it does so is set by how the two are
phased, and there are three ways to do it.

``same_plane``
    The second node trails the first in the same orbit plane by ``dt``.  It
    occupies the first node's inertial position ``dt`` later, but the Earth and
    the co-rotating ionosphere have turned east under it by ``omega_E dt``, so
    in the Earth-fixed frame its ground track lies west of the first by
    ``omega_E r cos(lat) dt``.  The part of that offset across the track is
    what matters.  This is what one launch and no propulsion gives.
``same_ground_track``
    The second node's plane is turned east by ``(omega_E - dOmega/dt) dt``, so
    its Earth-fixed track repeats the first node's exactly.  Both planes
    precess at the same J2 rate, so the arrangement holds once it is set up.
    It costs a plane change, or a period of differential nodal precession.
``crossing``
    The second node crosses the point on a descending pass while the first
    crosses it ascending, or in a plane of another inclination.  The two tracks
    then cross at an angle, which is what makes the cross-track tilt
    observable.

All three are computed from :class:`tss.platform.orbits.CircularOrbit`, so the
offsets and velocities come from the orbit and not from a flat-Earth sketch.
Velocities are Earth-fixed, which is the velocity relative to an ionosphere
that co-rotates with the Earth, and are resolved into local east and north at
the point.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize_scalar

from tss.platform.orbits import EARTH_ROTATION_RAD_PER_S, MU, CircularOrbit

__all__ = ["NodePass", "FormationGeometry", "two_node_geometry", "MODES"]

MODES = ("same_plane", "same_ground_track", "crossing")


@dataclass(frozen=True)
class NodePass:
    """One node's pass over the reference point.

    ``offset_en_m`` is the node's position at closest approach relative to the
    point, east and north, at the platform's altitude.  ``velocity_en_mps`` is
    its Earth-fixed velocity in east and north.  ``heading_deg`` is the
    direction of that velocity, clockwise from north.
    """

    t_s: float
    offset_en_m: np.ndarray
    velocity_en_mps: np.ndarray
    heading_deg: float


@dataclass(frozen=True)
class FormationGeometry:
    """Two passes over one point, and what the arrangement costs.

    ``plane_angle_deg`` is the angle between the two orbit planes;
    ``plane_change_dv_mps`` is the impulsive velocity change that would rotate
    one into the other.  ``j2_drift_days`` is how long a temporary 20 km
    altitude offset takes to open the required node difference by differential
    J2 precession instead, with ``j2_drift_dv_mps`` for raising and lowering.
    Both are zero for ``same_plane`` and are not computed for ``crossing``,
    which needs a second plane chosen at launch.
    """

    mode: str
    inclination_deg: float
    inclination2_deg: float
    latitude_deg: float
    dt_s: float
    node1: NodePass
    node2: NodePass
    crossing_angle_deg: float
    cross_track_offset_m: float
    along_track_offset_m: float
    raan_offset_deg: float
    plane_angle_deg: float
    plane_change_dv_mps: float
    j2_drift_days: float
    j2_drift_dv_mps: float


def _enu(up: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    east = np.cross([0.0, 0.0, 1.0], up)
    east = east / np.linalg.norm(east)
    north = np.cross(up, east)
    return east, north


def _pass(orbit: CircularOrbit, t_guess: float, up: np.ndarray, east: np.ndarray,
          north: np.ndarray) -> NodePass:
    """Closest approach of ``orbit`` to the direction ``up`` near ``t_guess``."""
    def ang(t: float) -> float:
        r, _ = orbit.ecef_state(t)
        return -float(np.dot(r / np.linalg.norm(r), up))

    res = minimize_scalar(ang, bounds=(t_guess - 900.0, t_guess + 900.0), method="bounded",
                          options={"xatol": 1e-4})
    t = float(res.x)
    r, v = orbit.ecef_state(t)
    a = np.linalg.norm(r)
    rel = r / a - up
    off = np.array([np.dot(rel, east), np.dot(rel, north)]) * a
    vel = np.array([np.dot(v, east), np.dot(v, north)])
    heading = float(np.degrees(np.arctan2(vel[0], vel[1])))
    return NodePass(t_s=t, offset_en_m=off, velocity_en_mps=vel, heading_deg=heading)


def _normal(orbit: CircularOrbit) -> np.ndarray:
    i = math.radians(orbit.inclination_deg)
    o = math.radians(orbit.raan0_deg)
    return np.array([math.sin(o) * math.sin(i), -math.cos(o) * math.sin(i), math.cos(i)])


def two_node_geometry(
    inclination_deg: float,
    latitude_deg: float,
    dt_s: float,
    mode: str,
    *,
    altitude_m: float = 500e3,
    inclination2_deg: float | None = None,
    j2_offset_m: float = 20e3,
) -> FormationGeometry:
    """Two passes over the point the first node crosses, ascending, at ``latitude_deg``."""
    if mode not in MODES:
        raise ValueError(f"mode must be one of {MODES}")
    inc = float(inclination_deg)
    lat = math.radians(float(latitude_deg))
    if abs(latitude_deg) > min(inc, 180.0 - inc):
        raise ValueError("the orbit does not reach that latitude")
    o1 = CircularOrbit(altitude_m, inc, raan0_deg=0.0, u0_deg=0.0)
    n = o1.mean_motion
    u1 = math.asin(math.sin(lat) / math.sin(math.radians(inc)))
    t1_guess = u1 / n
    r1, _ = o1.ecef_state(t1_guess)
    up = r1 / np.linalg.norm(r1)
    east, north = _enu(up)
    p1 = _pass(o1, t1_guess, up, east, north)
    # Re-centre on the exact closest approach so node 1 passes over the point.
    r1, _ = o1.ecef_state(p1.t_s)
    up = r1 / np.linalg.norm(r1)
    east, north = _enu(up)
    p1 = _pass(o1, p1.t_s, up, east, north)

    inc2 = inc if inclination2_deg is None else float(inclination2_deg)
    raan_off = 0.0
    if mode == "same_plane":
        o2 = CircularOrbit(altitude_m, inc, raan0_deg=0.0, u0_deg=-math.degrees(n * dt_s))
    elif mode == "same_ground_track":
        raan_off = math.degrees((EARTH_ROTATION_RAD_PER_S - o1.raan_rate) * dt_s)
        o2 = CircularOrbit(altitude_m, inc, raan0_deg=raan_off, u0_deg=-math.degrees(n * dt_s))
    else:
        t2 = p1.t_s + dt_s
        tmp = CircularOrbit(altitude_m, inc2)
        n2 = tmp.mean_motion
        if abs(latitude_deg) > min(inc2, 180.0 - inc2):
            raise ValueError("the second orbit does not reach that latitude")
        u2 = math.pi - math.asin(math.sin(lat) / math.sin(math.radians(inc2)))
        lon_p = math.atan2(up[1], up[0])
        ra_p = lon_p + EARTH_ROTATION_RAD_PER_S * t2
        ci2 = math.cos(math.radians(inc2))
        node_at_t2 = ra_p - math.atan2(ci2 * math.sin(u2), math.cos(u2))
        raan0 = node_at_t2 - tmp.raan_rate * t2
        o2 = CircularOrbit(altitude_m, inc2, raan0_deg=math.degrees(raan0),
                           u0_deg=math.degrees(u2 - n2 * t2))
        raan_off = math.degrees(raan0)
    p2 = _pass(o2, p1.t_s + dt_s, up, east, north)

    v1, v2 = p1.velocity_en_mps, p2.velocity_en_mps
    cosang = float(np.dot(v1, v2) / (np.linalg.norm(v1) * np.linalg.norm(v2)))
    cross_ang = float(np.degrees(np.arccos(np.clip(cosang, -1.0, 1.0))))
    u_hat = v1 / np.linalg.norm(v1)
    c_hat = np.array([u_hat[1], -u_hat[0]])
    cross_off = float(np.dot(p2.offset_en_m - p1.offset_en_m, c_hat))
    along_off = float(np.dot(p2.offset_en_m - p1.offset_en_m, u_hat))

    plane_angle = float(np.degrees(np.arccos(np.clip(np.dot(_normal(o1), _normal(o2)), -1.0, 1.0))))
    v_orb = math.sqrt(MU / o1.a)
    if mode == "same_ground_track":
        dv = 2.0 * v_orb * math.sin(math.radians(plane_angle) / 2.0)
        rate = 3.5 * abs(o1.raan_rate) * j2_offset_m / o1.a           # rad/s
        days = math.radians(raan_off) / rate / 86400.0 if rate > 0 else float("inf")
        dv_j2 = v_orb * j2_offset_m / o1.a                           # raise, then lower
    elif mode == "same_plane":
        dv, days, dv_j2 = 0.0, 0.0, 0.0
    else:
        dv, days, dv_j2 = float("nan"), float("nan"), float("nan")

    return FormationGeometry(
        mode=mode, inclination_deg=inc, inclination2_deg=inc2, latitude_deg=float(latitude_deg),
        dt_s=float(p2.t_s - p1.t_s), node1=p1, node2=p2, crossing_angle_deg=cross_ang,
        cross_track_offset_m=cross_off, along_track_offset_m=along_off, raan_offset_deg=raan_off,
        plane_angle_deg=plane_angle, plane_change_dv_mps=dv, j2_drift_days=days, j2_drift_dv_mps=dv_j2,
    )
