"""WP5 Phase 3 part 3, decision (e): the diffuse halo around a facet.

Andrew, 2026-11-27 (option A): weak point scatterers, equal in all directions,
placed at random along the track over ``ratio`` facet lengths, moving with the
medium, their power about constant through the pass; K is the facet's power
over the halo's mean power at overflight.  ``K = 0`` leaves the facet out and
keeps the halo at the power the facet would have had.  The decorrelation time
is the uniform patch's of ``l2_instrument.fading`` (TV equation c3taud).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from tss import sweep_law
from tss.constants import C_LIGHT
from tss.interfaces import PlatformState
from tss.l1_propagation import plate_echo as pe
from tss.l2_instrument.baseline import BaselineInstrument
from tss.l2_instrument.fading import FadingModel, halo_decorrelation_s
from tss.l3_signal.drift import DriftProcessing, OnboardPrior
from tss.pipeline import Pipeline
from tss.platform.orbits import CircularOrbit, platform_state_at
from tss.registry import build

F = 5.0e6
PROF = pe.ExponentialTopside()
DEPTH = PROF.depth_for_mu(F, 0.5)
FACET = {"id": "f1", "depth_m": DEPTH, "drift_mps": 25.0}
PRIOR = OnboardPrior.exponential(f_n_sat_hz=1.5e6, scale_height_m=200e3, ref_alt_m=600e3, source="truth")


def _echoes(facet, t=0.0):
    tier = build("l1", "facet", facets=[facet])
    plat = PlatformState(alt_m=600e3, lat_deg=0.0, lon_deg=0.0, t_s=t)
    alt = np.linspace(100e3, 600e3, 51)
    state = build("l0", "chapman").state_at(0.0, 0.0, "2026-01-01T12:00:00Z", alt)
    return tier, tier.echoes(state, plat, np.array([F]), np.random.default_rng(0))


def _rows(e, prefix):
    return np.array([str(f).startswith(prefix) for f in e.flags], dtype=bool)


def _field(e, prefix):
    """The summed complex field of the rows, at the carrier (amplitude and phase delay)."""
    m = _rows(e, prefix)
    return complex(np.sum(10 ** (e.amplitude_db[m] / 20) * np.exp(-2j * np.pi * F * e.phase_delay_s[m])))


def _halo(k, ratio=2.0, n=200, seed=0):
    return {**FACET, "halo": {"k": k, "ratio": ratio, "n_points": n, "seed": seed}}


def test_without_a_halo_or_with_k_infinite_the_rows_are_the_facets_alone():
    _, plain = _echoes(FACET, 0.3)
    _, inf = _echoes(_halo(math.inf), 0.3)
    for name in ("amplitude_db", "group_delay_s", "phase_delay_s", "doppler_hz", "freq_hz"):
        np.testing.assert_array_equal(getattr(plain, name), getattr(inf, name))
    assert list(plain.flags) == list(inf.flags) == ["facet:f1"]


def test_halo_rows_are_flagged_apart_and_lie_around_the_facet():
    tier, e = _echoes(_halo(3.0, ratio=4.0, n=50))
    assert int(_rows(e, "facet:f1").sum()) == 1 and int(_rows(e, "halo:f1").sum()) == 50
    d = tier.describe()["facets"][0]["halo"]
    assert d["length_m"] == 8.0e3 and d["k"] == 3.0 and d["n_points"] == 50 and d["seed"] == 0
    # Points up to 4 km to either side.  In the stratified medium a ray to x has
    # the group path P'(0) + x^2 / (2 P'^2) times the integral of dz / mu^3.
    z = (np.arange(4000) + 0.5) * DEPTH / 4000
    mu = np.sqrt(PROF.mu_squared(F, z))
    dz = DEPTH / 4000
    at_edge = 4.0e3**2 * float(np.sum(dz / mu**3)) / (2 * float(np.sum(dz / mu)) ** 2)
    extra = (e.group_delay_s[_rows(e, "halo:")] - e.group_delay_s[_rows(e, "facet:")][0]) * C_LIGHT / 2
    assert np.all(extra > -1.0) and 0.5 * at_edge < extra.max() < 1.02 * at_edge + 1.0


@pytest.mark.parametrize("k", [10.0, 1.0])
def test_each_scatterer_carries_an_equal_share_of_one_over_k(k):
    """Amplitude a with N a^2 = |E_facet(overflight)|^2 / K, spreading apart."""
    _, e = _echoes(_halo(k, n=100))
    spread = 20 * np.log10(C_LIGHT / F / (4 * np.pi * C_LIGHT * e.group_delay_s))
    level = e.amplitude_db - spread
    facet = level[_rows(e, "facet:")][0]
    halo = level[_rows(e, "halo:")]
    assert halo.size == 100
    np.testing.assert_allclose(halo, facet - 10 * math.log10(k * 100), atol=1e-9)


@pytest.mark.parametrize(("k", "t"), [(3.0, 0.0), (3.0, 1.5), (0.0, 0.0)])
def test_the_halos_mean_power_is_the_facets_over_k_all_through_the_pass(k, t):
    """Over 60 draws; at K = 0 the facet is absent and the halo has its overflight power."""
    _, ref = _echoes(FACET, 0.0)
    p_ref = abs(_field(ref, "facet:")) ** 2
    power = []
    for seed in range(60):
        _, e = _echoes(_halo(k, n=100, seed=seed), t)
        assert int(_rows(e, "facet:").sum()) == (0 if k == 0.0 else 1)
        power.append(abs(_field(e, "halo:")) ** 2)
    assert np.mean(power) / p_ref == pytest.approx(1.0 if k == 0.0 else 1.0 / k, rel=0.3)


def test_the_halo_decorrelates_as_a_uniform_patch_does():
    """At the patch's 1/e time (TV equation c3taud) the field correlation over
    200 draws is 1/e."""
    _, e0 = _echoes(FACET, 0.0)
    p_group = float(e0.group_delay_s[0]) * C_LIGHT / 2
    ratio = 4.0
    tau = halo_decorrelation_s(F, p_group, ratio * 2.0e3)
    cross, power = 0j, 0.0
    for seed in range(200):
        facet = {**FACET, "drift_mps": 0.0, "halo": {"k": 1.0, "ratio": ratio, "n_points": 100, "seed": seed}}
        a = _field(_echoes(facet, 0.0)[1], "halo:")
        b = _field(_echoes(facet, tau)[1], "halo:")
        cross += b * a.conjugate()
        power += 0.5 * (abs(a) ** 2 + abs(b) ** 2)
    assert abs(cross) / power == pytest.approx(math.exp(-1), abs=0.12)


def test_the_halo_moves_with_the_medium():
    """Points gathered at the facet's centre carry its drift's Doppler."""
    _, e = _echoes({**FACET, "halo": {"k": 3.0, "ratio": 1e-4, "n_points": 3, "seed": 1}})
    mu_sat = math.sqrt(float(PROF.mu_squared(F, 0.0)))
    want = 2 * F * mu_sat * 25.0 / C_LIGHT
    assert _rows(e, "halo:").sum() == 3
    np.testing.assert_allclose(e.doppler_hz[_rows(e, "halo:")], want, atol=0.01)


def test_the_draw_follows_its_seed():
    a, b, c = (_echoes(_halo(3.0, n=20, seed=s))[1] for s in (5, 5, 6))
    np.testing.assert_array_equal(a.phase_delay_s, b.phase_delay_s)
    assert not np.allclose(a.phase_delay_s, c.phase_delay_s)


@pytest.mark.parametrize("halo", [{"ratio": 2.0}, {"k": -1.0}, {"k": 3.0, "ratio": 0.0},
                                  {"k": 3.0, "n_points": 0}, {"k": 3.0, "radius": 1.0}])
def test_a_bad_halo_is_refused(halo):
    with pytest.raises(ValueError):
        build("l1", "facet", facets=[{**FACET, "halo": halo}])


def test_a_halo_is_refused_where_it_is_not_modelled():
    with pytest.raises(ValueError, match="halo"):
        build("l1", "facet", geometry="orbit", plate_model="2d", ref_alt_m=600e3,
              facets=[{"id": "o", "lat_deg": 0.0, "lon_deg": 0.0, "alt_m": 300e3, "halo": {"k": 3.0}}])
    with pytest.raises(ValueError, match="halo"):
        build("l1", "facet", facets=[{"id": "s", "kind": "sheet", "height_m": 200e3, "halo": {"k": 3.0}}])


def test_a_3d_orbit_plate_carries_its_halo():
    orbit = CircularOrbit(600e3, 60.0, u0_deg=20.0)
    p0 = platform_state_at(orbit, 0.0)
    facet = {"id": "o", "lat_deg": p0.lat_deg, "lon_deg": p0.lon_deg, "alt_m": 300e3, "drift_mps": 50.0,
             "halo": {"k": 3.0, "ratio": 2.0, "n_points": 40, "seed": 2}}
    tier = build("l1", "facet", geometry="orbit", ref_alt_m=600e3, facets=[facet])
    alt = np.linspace(50e3, 700e3, 201)
    state = build("l0", "chapman").state_at(p0.lat_deg, p0.lon_deg, "2026-01-01T12:00:00Z", alt)
    e = tier.echoes(state, p0, np.array([F]), np.random.default_rng(0))
    assert int(_rows(e, "facet:o").sum()) == 1 and int(_rows(e, "halo:o").sum()) == 40
    spread = 20 * np.log10(C_LIGHT / F / (4 * np.pi * C_LIGHT * e.group_delay_s))
    level = e.amplitude_db - spread
    # The platform is overhead at t = 0, so the reference is this facet row.
    np.testing.assert_allclose(level[_rows(e, "halo:")], level[_rows(e, "facet:")][0] - 10 * math.log10(120),
                               atol=0.5)
    assert np.all(np.isfinite(e.doppler_hz)) and np.all(np.isfinite(e.phase_delay_s))


def test_the_facets_fading_gain_does_not_touch_the_halo():
    _, e = _echoes(_halo(3.0, n=10))
    l2 = BaselineInstrument(fading=FadingModel(k_factor=math.inf, specular="rayleigh", tau_i_s=0.5,
                                               tau_log_sigma=0.0))
    gains = l2.fading_gains([e], np.array([0.0]), np.random.default_rng(3))
    assert list(gains) == ["f1"]
    out = l2._apply_gains(e, gains, 0)
    halo = _rows(e, "halo:")
    assert halo.sum() == 10
    np.testing.assert_array_equal(out.amplitude_db[halo], e.amplitude_db[halo])
    np.testing.assert_array_equal(out.phase_delay_s[halo], e.phase_delay_s[halo])
    assert out.amplitude_db[0] != e.amplitude_db[0]


def _sequence(facet, n=91):
    tiers = {"l0": build("l0", "chapman"), "l1": build("l1", "facet", facets=[facet]),
             "l2": BaselineInstrument(quantise=False, thermal_noise=False)}
    t = sweep_law.drift_mode_times(n, t0_s=-(n // 2) / 45 - 1 / 90)
    plats = [PlatformState(alt_m=600e3, lat_deg=0.0, lon_deg=0.0) for _ in t]
    return Pipeline(tiers, seed=11).run_sequence("d", sweep=sweep_law.drift_mode(F), t_s=t, platforms=plats,
                                                  alt_m=np.linspace(100e3, 600e3, 101),
                                                  epoch="2026-01-01T12:00:00Z").sequence


def test_through_the_chain_a_weak_halo_perturbs_the_drift_and_does_not_replace_it():
    plain = _sequence(FACET)
    weak = _sequence(_halo(1000.0))
    a, b = plain.records[45].samples, weak.records[45].samples
    # A halo 30 dB down changes the samples by about 3 per cent.
    assert 3e-3 < np.max(np.abs(a - b)) / np.max(np.abs(a)) < 0.3
    p = DriftProcessing(prior=PRIOR).process_sequence(weak).products
    assert p["echo_class"][0] == "compact_facet"
    assert p["drift_mps"][0] == pytest.approx(25.0, abs=3.0)


def test_route_ones_slope_test_does_not_reject_a_random_patch():
    """Measured 2026-11-27 (docs/drift-tier.md section 14), from route one's own
    Doppler law with no L1 or L2: every point of a randomly phased patch sweeps
    at the ephemeris slope, so the free-slope test reads about 1 and only the
    contrast test can decline it.  A tapered specular plate reads a slope near 0
    and is declined; a rough one is not."""
    from tss.l3_signal.novel.route_one_estimator import estimate_route_one

    f, p_group, mu, v_s = 4.988e6, 529.8e3, 0.9537, 7600.0
    t = (np.arange(209) - 104) / 45.0
    slope = 2 * f * v_s**2 / (C_LIGHT * p_group)
    kappa = mu * p_group / v_s**2
    ratios, formed = [], 0
    for seed in range(60):
        rng = np.random.default_rng(seed)
        ti = rng.uniform(-2000.0, 2000.0, 200) / v_s
        theta = rng.uniform(-math.pi, math.pi, 200)
        u = t[:, None] - ti[None, :]
        y = np.exp(1j * (theta + 2 * math.pi * slope * (kappa * 25.0 * u - 0.5 * u * u))).sum(axis=1)
        e = estimate_route_one(t, y / np.abs(y).max(), f_hz=f, group_path_m=p_group, mu_sat=mu, v_s_mps=v_s)
        ratios.append(e.slope_ratio)
        formed += e.separable
    assert np.mean(np.abs(np.array(ratios) - 1.0) <= 0.2) > 0.9
    assert formed >= 6

