"""WP5 Phase 2: the facet tier on an orbit, in the L0 profile (Andrew, 2026-11-25).

New, unverified geometry: the plate model is flat, so the orbit's curvature
enters through the Earth-flattening map (``x = R psi``, ``d = R ln(R/r)``,
``mu'^2 = mu^2 (r/R)^2``), which is conformal and keeps the phase path of every
ray exactly.  The group path is ``d(f P)/df``.  These tests hold the mapping to
the spherical geometry it stands for.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from tss.constants import C_LIGHT
from tss.interfaces import IonosphereState, PlatformState
from tss.l1_propagation import plate_echo as pe
from tss.l1_propagation.doppler import doppler_shift_hz
from tss.platform.orbits import RE, CircularOrbit, platform_state_at
from tss.registry import build

F = 5.0e6
ORBIT = CircularOrbit(600e3, 60.0, u0_deg=20.0)
FACET_ALT = 300e3


def _plat(t):
    p = platform_state_at(ORBIT, t)
    return PlatformState(alt_m=p.alt_m, lat_deg=p.lat_deg, lon_deg=p.lon_deg,
                         velocity_ned_mps=p.velocity_ned_mps, t_s=t)


def _state(f_n_sat=1.5e6, h=200e3):
    alt = np.linspace(50e3, 700e3, 1301)
    fn2 = f_n_sat**2 * np.exp((600e3 - alt) / h)
    return IonosphereState(alt_m=alt, ne=fn2 / 80.616)


def _tier(medium="plasma", profile="exponential", drift=0.0, length=2.0e3):
    p0 = _plat(0.0)
    return build("l1", "facet", geometry="orbit", medium=medium, profile=profile, ref_alt_m=600e3,
                 facets=[{"id": "o1", "lat_deg": p0.lat_deg, "lon_deg": p0.lon_deg,
                          "alt_m": FACET_ALT, "drift_mps": drift, "length_m": length}])


def test_free_space_range_is_the_spherical_chord():
    tier = _tier(medium="free")
    for t in (0.0, 0.4, -0.8):
        p = _plat(t)
        e = tier.echoes(_state(), p, np.array([F]), np.random.default_rng(0))
        # Exact chord from the platform to the facet's centre.
        f0 = _plat(0.0)
        lat1, lon1, lat2, lon2 = map(math.radians, (p.lat_deg, p.lon_deg, f0.lat_deg, f0.lon_deg))
        psi = math.acos(math.sin(lat1) * math.sin(lat2)
                        + math.cos(lat1) * math.cos(lat2) * math.cos(lon1 - lon2))
        rp, rf = RE + p.alt_m, RE + FACET_ALT
        chord = math.sqrt(rp**2 + rf**2 - 2 * rp * rf * math.cos(psi))
        assert e.group_delay_s[0] == pytest.approx(2 * chord / C_LIGHT, rel=2e-6)


def test_doppler_slope_carries_the_curvature():
    """Near overflight d^2 r/dt^2 = v^2 (r_f / r_p) / D, about 7 % below the flat v^2 / D.

    A 50 m plate, far inside the Fresnel zone, so it scatters as a point; a
    2 km plate is part mirror and its slope is smaller (route one's envelope).
    """
    tier = _tier(medium="free", length=50.0)
    h = 0.2
    fd = [tier.echoes(_state(), _plat(t), np.array([F]), np.random.default_rng(0)).doppler_hz[0]
          for t in (-h, h)]
    slope = (fd[1] - fd[0]) / (2 * h)  # Hz/s, negative (approach then recede)
    p = _plat(0.0)
    v = float(np.hypot(*p.velocity_ned_mps[:2]))
    rp, rf = RE + p.alt_m, RE + FACET_ALT
    d = rp - rf
    spherical = -2 * F / C_LIGHT * v**2 * (rf / rp) / d
    flat = -2 * F / C_LIGHT * v**2 / d
    assert slope == pytest.approx(spherical, rel=5e-3)
    assert slope / flat == pytest.approx(rf / rp, rel=5e-3)
    assert 1 - rf / rp > 0.04


def test_overflight_doppler_is_the_drift_doppler():
    tier = _tier(drift=25.0)
    e = tier.echoes(_state(), _plat(0.0), np.array([F]), np.random.default_rng(0))
    mu = math.sqrt(1 - (1.5e6 / F) ** 2)
    assert e.doppler_hz[0] == pytest.approx(float(doppler_shift_hz(25.0, F, mu_sat=mu)), rel=2e-3)


def test_l0_profile_matches_the_exponential_on_a_fine_grid():
    """The tabulated L0 profile gives the exponential's phase history: same
    phase delay to well under a cycle, and the same second difference across
    sweeps (the quantity route one reads)."""
    a, b = _tier(profile="exponential"), _tier(profile="l0")
    st = _state()
    ts = np.arange(-5, 6) / 45.0
    pa = np.array([a.echoes(st, _plat(t), np.array([F]), np.random.default_rng(0)).phase_delay_s[0] for t in ts])
    pb = np.array([b.echoes(st, _plat(t), np.array([F]), np.random.default_rng(0)).phase_delay_s[0] for t in ts])
    assert np.max(np.abs(pa - pb)) * F < 0.01
    d2a, d2b = np.diff(pa, 2) * F * 2 * np.pi, np.diff(pb, 2) * F * 2 * np.pi
    assert np.max(np.abs(d2a - d2b)) < 1e-3


def test_cross_track_offset_is_recorded():
    tier = _tier(medium="free")
    tier.echoes(_state(), _plat(1.0), np.array([F]), np.random.default_rng(0))
    assert 0.0 <= tier.max_cross_track_m < 50.0


def test_orbit_mode_needs_the_facet_position():
    with pytest.raises(ValueError):
        build("l1", "facet", geometry="orbit", ref_alt_m=600e3, facets=[{"id": "x", "depth_m": 1e5}])


def test_straight_mode_unchanged():
    t = build("l1", "facet", facets=[{"id": "f1", "depth_m": 300e3}])
    e = t.echoes(_state(), PlatformState(alt_m=600e3, lat_deg=0, lon_deg=0, t_s=0.0), np.array([F]),
                 np.random.default_rng(0))
    ref = pe.echo_track(np.array([0.0]), F, pe.Plate(depth_m=300e3), pe.ExponentialTopside())
    assert e.group_delay_s[0] == pytest.approx(2 * ref.group_path_m[0] / C_LIGHT, rel=1e-12)
