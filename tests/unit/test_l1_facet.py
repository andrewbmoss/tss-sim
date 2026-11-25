"""WP5 Phase 2: the l1/facet tier (D4: plate_echo wrapped per sweep)."""

from __future__ import annotations

import numpy as np
import pytest

from tss.constants import C_LIGHT
from tss.interfaces import IonosphereState, PlatformState
from tss.l1_propagation import plate_echo as pe
from tss.l1_propagation.doppler import doppler_shift_hz
from tss.registry import build

PROF = pe.ExponentialTopside()
F = 5.0e6
DEPTH = PROF.depth_for_mu(F, 0.5)


def _state():
    alt = np.linspace(100e3, 600e3, 51)
    return IonosphereState(alt_m=alt, ne=np.full(alt.size, 1e11))


def _plat(t):
    return PlatformState(alt_m=600e3, lat_deg=0.0, lon_deg=0.0, t_s=t)


def _tier(**kw):
    facets = kw.pop("facets", [{"id": "f1", "depth_m": DEPTH, "drift_mps": 25.0}])
    return build("l1", "facet", facets=facets, **kw)


def test_registered_and_describes_itself():
    t = _tier()
    assert t.layer == "l1" and t.name == "facet"
    d = t.describe()
    assert d["geometry"] == "straight" and d["facets"][0]["id"] == "f1"
    assert t.group_index_used == "field_free"


def test_rows_carry_the_plate_echo():
    t = _tier()
    e = t.echoes(_state(), _plat(0.3), np.array([F]), np.random.default_rng(0))
    assert list(e.flags) == ["facet:f1"]
    ref = pe.echo_track(np.array([0.3]), F, pe.Plate(depth_m=DEPTH, drift_mps=25.0), PROF)
    assert e.group_delay_s[0] == pytest.approx(2 * ref.group_path_m[0] / C_LIGHT, rel=1e-12)
    # The carrier phase 2 pi f tau_p is the field's phase, in L2's sign.
    want = -np.angle(ref.field[0])
    got = 2 * np.pi * F * e.phase_delay_s[0]
    assert abs(np.angle(np.exp(1j * (got - want)))) < 1e-9
    assert e.phase_delay_s[0] < e.group_delay_s[0]
    lam = C_LIGHT / F
    spread = 20 * np.log10(lam / (4 * np.pi * 2 * ref.group_path_m[0]))
    assert e.amplitude_db[0] == pytest.approx(spread + 20 * np.log10(abs(ref.field[0])), abs=1e-9)


def test_doppler_at_overflight_is_the_drift_doppler():
    t = _tier()
    e = t.echoes(_state(), _plat(0.0), np.array([F]), np.random.default_rng(0))
    mu = float(np.sqrt(PROF.mu_squared(F, 0.0)))
    want = float(doppler_shift_hz(25.0, F, mu_sat=mu))
    assert want > 0.0  # closing (drift upward toward the platform)
    assert e.doppler_hz[0] == pytest.approx(want, rel=1e-3)


def test_doppler_sign_follows_the_platform_motion():
    t = _tier(facets=[{"id": "f1", "depth_m": DEPTH}])
    before = t.echoes(_state(), _plat(-1.0), np.array([F]), np.random.default_rng(0))
    after = t.echoes(_state(), _plat(1.0), np.array([F]), np.random.default_rng(0))
    assert before.doppler_hz[0] > 0.0 > after.doppler_hz[0]  # approaching, then receding


def test_a_facet_beyond_reflection_is_masked_not_raised():
    t = _tier()
    e = t.echoes(_state(), _plat(0.0), np.array([1.0e6, 2.0e6, F]), np.random.default_rng(0))
    assert np.isnan(e.group_delay_s[0]) and np.isnan(e.group_delay_s[1])
    assert str(e.flags[0]).startswith("facet:f1") and "beyond_reflection" in str(e.flags[0])
    assert np.isfinite(e.group_delay_s[2])


def test_free_space_medium():
    t = _tier(medium="free", facets=[{"id": "a", "depth_m": 500e3}])
    e = t.echoes(_state(), _plat(0.0), np.array([F]), np.random.default_rng(0))
    assert e.group_delay_s[0] == pytest.approx(2 * 500e3 / C_LIGHT, rel=1e-12)
    # The summed field of a plate smaller than the Fresnel zone is not the
    # centre ray's phase, but it is within a fraction of a cycle of it.
    assert abs(e.phase_delay_s[0] - e.group_delay_s[0]) * F < 0.125


def test_no_draws_from_the_generator():
    rng = np.random.default_rng(5)
    before = rng.bit_generator.state
    _tier().echoes(_state(), _plat(0.0), np.array([F]), rng)
    assert rng.bit_generator.state == before
