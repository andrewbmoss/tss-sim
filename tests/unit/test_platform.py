"""Unit tests for the platform package.

Two things are checked against numbers that exist outside this codebase, which
is the point of a unit test on an orbit: the circular orbit period against the
closed-form Keplerian value, and the sun-synchronous inclination against the
published figures for the altitudes mission designers quote them at.

Everything else here holds a property the rest of the simulator relies on.  The
Earth-fixed velocity is the one that matters most, because C3 is the claim that
platform motion can be separated from range, and the claim is only as good as
the velocity the simulator says the platform has.
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[2]
SRC = REPO / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from tss.interfaces import PlatformState  # noqa: E402
from tss.platform import (  # noqa: E402
    MU,
    RE,
    CircularOrbit,
    orbit_classes,
    platform_state_at,
    sun_synchronous_inclination,
)


# ---------------------------------------------------------------------------
# Period against the analytic value
# ---------------------------------------------------------------------------
def test_circular_orbit_period_matches_the_analytic_value():
    for altitude in (400.0e3, 500.0e3, 600.0e3, 800.0e3, 1200.0e3):
        orbit = CircularOrbit(altitude, 80.0)
        analytic = 2.0 * math.pi * math.sqrt((RE + altitude) ** 3 / MU)
        assert orbit.period == pytest.approx(analytic, rel=1e-12)


def test_period_and_speed_are_consistent_with_each_other():
    """``v T = 2 pi a`` for a circular orbit.  Two properties, one geometry."""
    orbit = CircularOrbit(600.0e3, 80.0)
    assert orbit.speed * orbit.period == pytest.approx(2.0 * math.pi * orbit.a, rel=1e-12)


def test_period_at_500_km_is_about_ninety_five_minutes():
    """A sanity anchor a reader can check without a calculator."""
    orbit = CircularOrbit(500.0e3, 80.0)
    assert 94.0 < orbit.period / 60.0 < 95.5


# ---------------------------------------------------------------------------
# Sun-synchronous inclination against known cases
# ---------------------------------------------------------------------------
def test_sun_synchronous_inclination_at_known_altitudes():
    """The textbook figures: about 97.4 degrees at 500 km, 98.6 at 800 km.

    These are the numbers quoted in every mission design reference, and they
    come out of the same secular J2 nodal rate this module implements, so
    agreeing with them checks the rate and not only the arithmetic around it.
    """
    assert sun_synchronous_inclination(500.0e3) == pytest.approx(97.4, abs=0.15)
    assert sun_synchronous_inclination(800.0e3) == pytest.approx(98.6, abs=0.15)
    assert sun_synchronous_inclination(600.0e3) == pytest.approx(97.79, abs=0.15)


def test_a_sun_synchronous_orbit_has_no_local_time_drift():
    """Which is exactly why MIS-03 excludes it.

    A mission that never drifts in local time never samples the local-time
    dependence the science requires, so the case is disqualified rather than
    merely inconvenient.
    """
    altitude = 600.0e3
    orbit = CircularOrbit(altitude, sun_synchronous_inclination(altitude))
    assert abs(orbit.lt_drift_hours_per_day) < 1.0e-3
    assert orbit.days_for_full_lt_coverage > 1.0e3


def test_the_orbit_classes_do_drift_in_local_time():
    for orbit in orbit_classes(600.0e3).values():
        assert abs(orbit.lt_drift_hours_per_day) > 0.05
        assert orbit.days_for_full_lt_coverage < 500.0


def test_nodal_precession_is_retrograde_for_a_prograde_orbit():
    assert CircularOrbit(600.0e3, 30.0).raan_rate < 0.0
    assert CircularOrbit(600.0e3, 120.0).raan_rate > 0.0


# ---------------------------------------------------------------------------
# The platform state, and the velocity C3 depends on
# ---------------------------------------------------------------------------
def test_platform_state_is_the_interface_object():
    orbit = CircularOrbit(600.0e3, 80.0)
    state = platform_state_at(orbit, 0.0, epoch="2026-01-01T00:00:00")
    assert isinstance(state, PlatformState)
    assert state.epoch == "2026-01-01T00:00:00"
    assert state.alt_m == pytest.approx(600.0e3, abs=1.0)
    assert state.velocity_ned_mps is not None
    assert state.velocity_ned_mps.shape == (3,)


def test_ecef_velocity_matches_a_finite_difference_of_ecef_position():
    """The velocity is analytic, so it has to agree with differencing.

    C3 rests on the motion being precisely known, so an analytic velocity that
    did not match the position it is the derivative of would be a defect that
    nothing else in the package would catch.
    """
    orbit = CircularOrbit(600.0e3, 80.0)
    dt = 1.0e-3
    for t in (0.0, 1234.5, 3000.0):
        r0, v0 = orbit.ecef_state(t)
        r1, _ = orbit.ecef_state(t + dt)
        fd = (r1 - r0) / dt
        assert np.allclose(fd, v0, rtol=0.0, atol=1.0e-2)


def test_ecef_speed_is_below_the_inertial_speed_for_a_prograde_orbit():
    """The frame rotation term is about 465 m/s at the equator and is real.

    Reporting the inertial speed in a field named ECEF would be a several per
    cent error in the along-track velocity, which is far larger than the Doppler
    precision C3 needs.
    """
    orbit = CircularOrbit(600.0e3, 30.0)
    state = platform_state_at(orbit, 0.0)
    ecef_speed = float(np.linalg.norm(state.velocity_ned_mps))
    assert ecef_speed < orbit.speed
    assert orbit.speed - ecef_speed == pytest.approx(465.0, abs=60.0)


def test_ned_from_ecef_has_the_right_axes_at_the_origin():
    """At latitude 0, longitude 0 the Earth-fixed x is up, y east and z north."""
    from tss.platform.orbits import ned_from_ecef

    assert np.allclose(ned_from_ecef(np.array([1.0, 0.0, 0.0]), 0.0, 0.0), [0.0, 0.0, -1.0])
    assert np.allclose(ned_from_ecef(np.array([0.0, 1.0, 0.0]), 0.0, 0.0), [0.0, 1.0, 0.0])
    assert np.allclose(ned_from_ecef(np.array([0.0, 0.0, 1.0]), 0.0, 0.0), [1.0, 0.0, 0.0])
    # At the north pole, Earth-fixed z is up.
    assert np.allclose(ned_from_ecef(np.array([0.0, 0.0, 1.0]), 90.0, 0.0), [0.0, 0.0, -1.0])


@pytest.mark.parametrize("t_s", [0.0, 700.0, 2100.0, 4000.0])
@pytest.mark.parametrize("inclination", [30.0, 80.0, 97.4])
def test_a_circular_orbit_neither_climbs_nor_falls(t_s, inclination):
    """The vertical component of the delivered velocity is zero for a circular
    orbit, and the horizontal speed is the Earth-fixed speed.  Until 2026-10-17
    the Earth-fixed z (northward at the equator) was read as the vertical."""
    orbit = CircularOrbit(600.0e3, inclination)
    state = platform_state_at(orbit, t_s)
    _, v_ecef = orbit.ecef_state(t_s)
    v = state.velocity_ned_mps
    assert abs(v[2]) < 1.0e-6 * orbit.speed
    assert np.hypot(v[0], v[1]) == pytest.approx(np.linalg.norm(v_ecef), rel=1e-12)


def test_altitude_stays_constant_around_a_circular_orbit():
    orbit = CircularOrbit(600.0e3, 80.0)
    altitudes = [
        platform_state_at(orbit, t).alt_m
        for t in np.linspace(0.0, orbit.period, 25)
    ]
    assert max(altitudes) - min(altitudes) < 1.0


def test_latitude_never_exceeds_the_inclination():
    orbit = CircularOrbit(600.0e3, 30.0)
    lats = [platform_state_at(orbit, t).lat_deg for t in np.linspace(0.0, orbit.period, 200)]
    assert max(abs(lat) for lat in lats) <= 30.0 + 1.0e-6


def test_platform_state_does_not_alias_the_callers_attitude_array():
    orbit = CircularOrbit(600.0e3, 80.0)
    quat = np.array([1.0, 0.0, 0.0, 0.0])
    before = quat.copy()
    state = platform_state_at(orbit, 0.0, attitude_quat=quat)
    assert state.attitude_quat is not quat
    state.attitude_quat[0] = 0.5
    assert np.array_equal(quat, before)


# ---------------------------------------------------------------------------
# Sampling
# ---------------------------------------------------------------------------
def test_sample_reproduces_from_a_seed_and_from_an_injected_generator():
    orbit = CircularOrbit(600.0e3, 80.0)
    a = orbit.sample(500, duration_days=30.0, seed=7)
    b = orbit.sample(500, duration_days=30.0, seed=7)
    assert np.array_equal(a["lat_deg"], b["lat_deg"])

    c = orbit.sample(500, duration_days=30.0, rng=np.random.default_rng(7))
    assert np.array_equal(a["t_seconds"], c["t_seconds"])


def test_sample_stays_inside_the_inclination_and_the_longitude_range():
    orbit = CircularOrbit(600.0e3, 80.0)
    out = orbit.sample(2000, duration_days=90.0, seed=1)
    assert np.all(np.abs(out["lat_deg"]) <= 80.0 + 1.0e-9)
    assert np.all(out["lon_deg"] >= -180.0) and np.all(out["lon_deg"] < 180.0)
    assert np.all((out["local_time_hours"] >= 0.0) & (out["local_time_hours"] < 24.0))
    assert np.all(out["altitude_m"] == 600.0e3)
