"""Layer: platform
WP: WP0, WP1
Claims: C3
Status: implemented
Requirements: MIS-03, SCI-02

Circular-orbit sampling and the platform state the rest of the simulator works
from.

WHY THE ORBIT AND NOT A GRID
----------------------------
The band-edge study asks what the local plasma frequency at the spacecraft is,
as a distribution.  That distribution depends on where the spacecraft actually
goes, so the sample space has to be the orbit track rather than a uniform
latitude and longitude grid.  A uniform grid over-weights the poles by roughly
``1/cos(lat)``, and the polar topside is thin, so a grid-sampled answer is
optimistic in exactly the tail the requirement is set from.

WHAT IS MODELLED AND WHAT IS NOT
--------------------------------
Only the secular J2 nodal precession.  That is what sweeps the orbit through
local time, and local time is the dominant control on topside density, so it is
the term that matters.  Drag, higher zonals and lunisolar terms change where the
satellite is on a given day but not the statistics of where it goes.

The Earth is treated as a sphere for the sub-satellite point, with the sidereal
day taken as 86164.0905 s.  Geodetic against geocentric latitude differs by up
to about 0.19 degrees, which is far below anything the topside statistics
resolve, and the difference is stated here rather than hidden.

WHY THE ECEF VELOCITY IS POPULATED
----------------------------------
C3 is the claim that platform motion can be separated from range using a single
transmit and receive pair.  The whole claim rests on the motion being large and
precisely known, so the velocity has to be a first-class field on
:class:`~tss.interfaces.PlatformState` rather than something reconstructed
downstream from two positions and a time difference.
:func:`platform_state_at` therefore returns it, computed analytically from the
orbit rather than by differencing, so it carries no finite-difference error and
does not depend on the step the caller happened to choose.

The velocity is computed in Earth-centred Earth-fixed components: the inertial
velocity rotated into the rotating frame, less the frame rotation term
``omega x r``, so at the equator it is about 465 m/s smaller than the inertial
speed for a prograde orbit and about 465 m/s larger for a retrograde one.  That
difference is roughly 6 per cent of the orbital speed and is far larger than the
Doppler precision C3 needs, so getting the frame wrong here would be a real
error and not a rounding one.

It is delivered on :class:`~tss.interfaces.PlatformState` in local
north-east-down components (:func:`ned_from_ecef`), because that is the frame
every consumer works in: L1's line of sight for a vertical sounding is local
down, and L5's along-track spacing is the horizontal speed.  The local frame is
taken at the spherical sub-satellite point, the same one the latitude and
longitude are, so a circular orbit has exactly zero vertical velocity.  Until
2026-10-17 the Earth-fixed components were delivered as they were and read as
north-east-down, which at the equator put the northward speed, 7.5 km/s, into
the vertical.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

import numpy as np

from tss.interfaces import PlatformState

__all__ = [
    "MU",
    "RE",
    "J2",
    "SOLAR_DAY",
    "SIDEREAL_DAY",
    "SUN_DRIFT_DEG_PER_DAY",
    "EARTH_ROTATION_RAD_PER_S",
    "CircularOrbit",
    "orbit_classes",
    "sun_synchronous_inclination",
    "platform_state_at",
]

MU = 3.986004418e14        # m^3/s^2
RE = 6378137.0             # m
J2 = 1.08262668e-3
SOLAR_DAY = 86400.0
SIDEREAL_DAY = 86164.0905
SUN_DRIFT_DEG_PER_DAY = 0.98565    # apparent motion of the sun in right ascension
EARTH_ROTATION_RAD_PER_S = 2.0 * math.pi / SIDEREAL_DAY


@dataclass
class CircularOrbit:
    """A circular orbit described by altitude and inclination."""

    altitude: float            # m
    inclination_deg: float
    raan0_deg: float = 0.0
    u0_deg: float = 0.0        # argument of latitude at epoch
    name: str = ""

    @property
    def a(self):
        return RE + self.altitude

    @property
    def mean_motion(self):
        """rad/s"""
        return math.sqrt(MU / self.a**3)

    @property
    def period(self):
        return 2.0 * math.pi / self.mean_motion

    @property
    def raan_rate(self):
        """Secular J2 nodal precession, rad/s.  Negative for prograde orbits."""
        return (
            -1.5
            * self.mean_motion
            * J2
            * (RE / self.a) ** 2
            * math.cos(math.radians(self.inclination_deg))
        )

    @property
    def lt_drift_hours_per_day(self):
        """
        Drift rate of the local time of the ascending node.  Zero means
        sun-synchronous, which MIS-03 disqualifies because the mission would
        never sample the local-time dependence the science requires.
        """
        d_deg_per_day = math.degrees(self.raan_rate) * SOLAR_DAY
        return (d_deg_per_day - SUN_DRIFT_DEG_PER_DAY) / 15.0

    @property
    def days_for_full_lt_coverage(self):
        r = self.lt_drift_hours_per_day
        return float("inf") if r == 0 else abs(24.0 / r)

    @property
    def speed(self):
        """Inertial speed of a circular orbit, m/s.  ``sqrt(mu/a)``.

        Added in the migration.  It is the number the Doppler bias of C3 scales
        with, so having it on the orbit rather than recomputed at each call site
        keeps a factor of ``a`` from being dropped somewhere.
        """
        return math.sqrt(MU / self.a)

    def sample(self, n_samples, duration_days=365.25, seed=0, rng=None):
        """
        Draw ``n_samples`` epochs uniformly at random over ``duration_days``.

        Random sampling rather than a fixed cadence is deliberate: an orbital
        period near 100 minutes will alias against any round-numbered cadence,
        and against the solar day.  Uniform random epochs are unbiased by
        construction and need no argument about commensurability.

        Returns a dict of arrays: t_days, lat_deg, lon_deg, local_time_hours,
        and the argument of latitude.

        ``rng`` was added in the migration and takes precedence over ``seed``.
        Contract rule 4 in :mod:`tss.interfaces` says a tier draws only from an
        injected generator; ``seed`` is kept so that the study scripts written
        against the original signature still run, and it builds a generator the
        same way the original did, so the same seed gives the same sample.
        """
        generator = np.random.default_rng(seed) if rng is None else rng
        t = generator.uniform(0.0, duration_days * SOLAR_DAY, int(n_samples))   # s

        inc = math.radians(self.inclination_deg)
        u = math.radians(self.u0_deg) + self.mean_motion * t
        raan = math.radians(self.raan0_deg) + self.raan_rate * t

        # Sub-satellite geodetic latitude of a circular orbit.
        lat = np.arcsin(np.sin(inc) * np.sin(u))

        # Right ascension of the sub-satellite point.
        ra = raan + np.arctan2(np.cos(inc) * np.sin(u), np.cos(u))

        # Earth rotation.  The epoch offset is arbitrary for a statistical
        # study; only the relative rotation matters.
        gmst = 2.0 * math.pi * (t / SIDEREAL_DAY)          # sidereal day
        lon = np.degrees(ra - gmst)
        lon = (lon + 180.0) % 360.0 - 180.0

        # Mean solar time, which is what IRI and NeQuick take.
        ut_hours = (t / 3600.0) % 24.0
        lt = (ut_hours + lon / 15.0) % 24.0

        return dict(
            t_days=t / SOLAR_DAY,
            t_seconds=t,
            lat_deg=np.degrees(lat),
            lon_deg=lon,
            ut_hours=ut_hours,
            local_time_hours=lt,
            arg_lat_deg=np.degrees(u) % 360.0,
            altitude_m=np.full(int(n_samples), self.altitude),
        )

    # -- state at one instant -----------------------------------------------

    def eci_state(self, t_s: float) -> tuple[np.ndarray, np.ndarray]:
        """Inertial position and velocity at ``t_s`` seconds after epoch, SI.

        The position of a circular orbit in the inertial frame is

            r = a [ cos O cos u - sin O sin u cos i,
                    sin O cos u + cos O sin u cos i,
                    sin u sin i ]

        and the velocity is its exact time derivative with ``u`` advancing at
        the mean motion and ``O`` at the J2 nodal rate.  The nodal rate is
        included in the derivative rather than treated as frozen, because
        leaving it out puts a systematic error of order ``a * dO/dt``, which is
        about half a metre per second, straight into the along-track velocity.
        That is small against 7.5 km/s and is not small against what C3 needs to
        resolve.
        """
        t = float(t_s)
        inc = math.radians(self.inclination_deg)
        u = math.radians(self.u0_deg) + self.mean_motion * t
        node = math.radians(self.raan0_deg) + self.raan_rate * t
        n = self.mean_motion
        node_dot = self.raan_rate
        a = self.a

        cu, su = math.cos(u), math.sin(u)
        co, so = math.cos(node), math.sin(node)
        ci, si = math.cos(inc), math.sin(inc)

        r = a * np.array(
            [co * cu - so * su * ci, so * cu + co * su * ci, su * si], dtype=float
        )
        # d/dt with both u and the node advancing.
        v = a * np.array(
            [
                -node_dot * (so * cu + co * su * ci) + n * (-co * su - so * cu * ci),
                node_dot * (co * cu - so * su * ci) + n * (-so * su + co * cu * ci),
                n * cu * si,
            ],
            dtype=float,
        )
        return r, v

    def ecef_state(self, t_s: float, gmst0_rad: float = 0.0) -> tuple[np.ndarray, np.ndarray]:
        """Earth-fixed position and velocity at ``t_s``, SI.

        The inertial state rotated by the Greenwich hour angle, with the frame
        rotation removed from the velocity:

            v_ecef = R(-theta) v_eci - omega x r_ecef

        ``gmst0_rad`` is the hour angle at epoch.  It defaults to zero, matching
        :meth:`sample`, which notes that the epoch offset is arbitrary for a
        statistical study.  A campaign that cares about real longitudes has to
        set it.
        """
        r_eci, v_eci = self.eci_state(t_s)
        theta = float(gmst0_rad) + EARTH_ROTATION_RAD_PER_S * float(t_s)
        c, s = math.cos(theta), math.sin(theta)
        rot = np.array([[c, s, 0.0], [-s, c, 0.0], [0.0, 0.0, 1.0]], dtype=float)
        r_ecef = rot @ r_eci
        omega = np.array([0.0, 0.0, EARTH_ROTATION_RAD_PER_S], dtype=float)
        v_ecef = rot @ v_eci - np.cross(omega, r_ecef)
        return r_ecef, v_ecef

    def describe(self) -> Mapping[str, Any]:
        """JSON-serialisable record for the run manifest."""
        return {
            "name": self.name,
            "altitude_m": float(self.altitude),
            "inclination_deg": float(self.inclination_deg),
            "raan0_deg": float(self.raan0_deg),
            "u0_deg": float(self.u0_deg),
            "period_s": float(self.period),
            "speed_mps": float(self.speed),
            "raan_rate_rad_per_s": float(self.raan_rate),
            "lt_drift_hours_per_day": float(self.lt_drift_hours_per_day),
            "days_for_full_lt_coverage": float(self.days_for_full_lt_coverage),
            "perturbations": "secular J2 nodal precession only",
        }


#: The orbit classes defined in TSS_Science_Requirements_RevB section 6.3.
def orbit_classes(altitude):
    """Return the O1 and O2 class orbits at a given altitude."""
    return {
        "O1": CircularOrbit(altitude, 80.0, name="O1 global morphology"),
        "O2": CircularOrbit(altitude, 30.0, name="O2 bubble tracking"),
    }


def sun_synchronous_inclination(altitude):
    """Inclination that makes the node track the sun -- the case MIS-03 excludes."""
    a = RE + altitude
    n = math.sqrt(MU / a**3)
    ci = (math.radians(SUN_DRIFT_DEG_PER_DAY) / SOLAR_DAY) / (-1.5 * n * J2 * (RE / a) ** 2)
    return math.degrees(math.acos(ci))


def _geocentric_lat_lon_alt(r_ecef: np.ndarray) -> tuple[float, float, float]:
    """Spherical sub-satellite latitude, longitude and altitude above ``RE``.

    Spherical rather than geodetic, to match :meth:`CircularOrbit.sample`.  The
    two differ by up to about 0.19 degrees in latitude and about 21 km in
    altitude at mid-latitudes, and the altitude difference is the one that would
    matter to a profile lookup, so a campaign that needs geodetic altitudes has
    to convert here rather than downstream.  That limitation is stated in the
    module docstring and is not hidden in this function.
    """
    x, y, z = (float(v) for v in r_ecef)
    radius = math.sqrt(x * x + y * y + z * z)
    lat = math.degrees(math.asin(np.clip(z / radius, -1.0, 1.0)))
    lon = (math.degrees(math.atan2(y, x)) + 180.0) % 360.0 - 180.0
    return lat, lon, radius - RE


def ned_from_ecef(v_ecef: np.ndarray, lat_deg: float, lon_deg: float) -> np.ndarray:
    """Earth-fixed Cartesian components to local north-east-down, at a point.

    ``lat_deg`` and ``lon_deg`` fix the local frame.  Pass the spherical
    latitude :func:`platform_state_at` reports, so that "down" is toward the
    centre of the Earth the orbit is computed around.  The conversion is a
    rotation, so it preserves the speed.
    """
    v = np.asarray(v_ecef, dtype=float)
    if v.shape != (3,):
        raise ValueError(f"an Earth-fixed velocity is a 3-vector; got shape {v.shape}")
    lat, lon = math.radians(float(lat_deg)), math.radians(float(lon_deg))
    sl, cl = math.sin(lat), math.cos(lat)
    so, co = math.sin(lon), math.cos(lon)
    north = np.array([-sl * co, -sl * so, cl])
    east = np.array([-so, co, 0.0])
    down = np.array([-cl * co, -cl * so, -sl])
    return np.array([north @ v, east @ v, down @ v], dtype=float)


def platform_state_at(
    orbit: CircularOrbit,
    t_s: float,
    *,
    epoch: str = "",
    gmst0_rad: float = 0.0,
    attitude_quat: np.ndarray | None = None,
) -> PlatformState:
    """The :class:`~tss.interfaces.PlatformState` at ``t_s`` seconds after epoch.

    Parameters
    ----------
    orbit
        The orbit.  Not modified.
    t_s
        Seconds after the orbit's epoch.
    epoch
        ISO 8601 UTC timestamp for the returned state, carried through
        unchanged.  This function does not do calendar arithmetic; the caller
        owns the mapping from ``t_s`` to a timestamp, because that mapping needs
        an epoch the orbit does not carry.
    gmst0_rad
        Greenwich hour angle at the orbit's epoch.  Zero by default, matching
        :meth:`CircularOrbit.sample`.
    attitude_quat
        Attitude, passed straight through.  Copied, so the returned state does
        not hold a reference to the caller's array.

    The velocity is populated because C3 depends on it; see the module
    docstring.  It is analytic, not a finite difference, and it is delivered in
    local north-east-down components.
    """
    r_ecef, v_ecef = orbit.ecef_state(float(t_s), gmst0_rad=gmst0_rad)
    lat, lon, alt = _geocentric_lat_lon_alt(r_ecef)
    return PlatformState(
        alt_m=float(alt),
        lat_deg=float(lat),
        lon_deg=float(lon),
        velocity_ned_mps=ned_from_ecef(v_ecef, lat, lon),
        attitude_quat=None
        if attitude_quat is None
        else np.array(attitude_quat, dtype=float, copy=True),
        epoch=str(epoch),
    )
