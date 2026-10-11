"""Route two, phase 2: stratified links, the network estimator, the link table and the study."""

from __future__ import annotations

import math
from dataclasses import replace

import numpy as np
import pytest

from tss.l1_propagation.bistatic_mirror import LinkGeometry, enu_basis, llh_to_ecef
from tss.l1_propagation.stratified_link import (
    QUIET_DAYTIME,
    PatchState,
    link_forward,
    link_sensitivity,
    max_central_angle,
    solve_ray,
)

ALT = 890e3
LAT0, LON0 = 20.0, -140.0


def _vel(lat, lon, head, v=7400.0):
    e, n, _ = enu_basis(lat, lon)
    h = math.radians(head)
    return v * (math.sin(h) * e + math.cos(h) * n)


def _link(la1, lo1, h1, la2, lo2, h2):
    return (llh_to_ecef(la1, lo1, ALT), llh_to_ecef(la2, lo2, ALT), _vel(la1, lo1, h1), _vel(la2, lo2, h2))


NS = _link(15.5, -140.0, 120.0, 24.5, -140.0, 60.0)           # 1000 km north-south, moving apart east
DIAG = _link(16.8, -143.4, 150.0, 23.2, -136.6, 30.0)


# -- the level ray --------------------------------------------------------------------------------

def test_level_rate_matches_phase_path_difference():
    """-K dTheta/dt (the analytic rate) equals minus the time difference of the traced phase path."""
    a, b, va, vb = NS
    s = link_sensitivity(QUIET_DAYTIME, LAT0, LON0, a, b, va, vb, 10e6)
    dt = 0.05
    pp = solve_ray(QUIET_DAYTIME, a + va * dt, b + vb * dt, 10e6).phase_path_m
    pm = solve_ray(QUIET_DAYTIME, a - va * dt, b - vb * dt, 10e6).phase_path_m
    assert s.rate_mps == pytest.approx(-(pp - pm) / (2 * dt), abs=0.05)
    # The ray-traced value of the same link (Haselgrove, nisc-verify) was -4636.08 m/s.
    assert s.rate_mps == pytest.approx(-4636.08, abs=0.1)


def test_vertical_group_path_and_drift_coefficient():
    a = llh_to_ecef(LAT0, LON0, ALT)
    ray = solve_ray(QUIET_DAYTIME, a, a, 6e6)
    assert ray.group_path_m == pytest.approx(1100424.1, abs=1.0)       # Haselgrove: 1100424.11 m
    s = link_sensitivity(QUIET_DAYTIME, LAT0, LON0, a, a, _vel(LAT0, LON0, 30), _vel(LAT0, LON0, 30), 6e6)
    mu = math.sqrt(1.0 - 80.6163849431291 * QUIET_DAYTIME.density(ALT) / 36e12)
    assert s.d_rate[0] == pytest.approx(2.0 * mu, rel=2e-3)          # a vertical echo sees 2 mu w


def test_no_ray_beyond_the_maximum_range():
    r = 6371e3 + ALT
    tmax = max_central_angle(QUIET_DAYTIME, r, 10e6)
    a = llh_to_ecef(20.0, -150.0, ALT)
    for frac, closes in ((0.95, True), (1.05, False)):
        lon = -150.0 + math.degrees(tmax * frac) / math.cos(math.radians(20.0))
        b = llh_to_ecef(20.0, lon, ALT)
        th = math.acos(np.dot(a, b) / r**2)
        if abs(th / tmax - frac) > 0.03:   # keep the test about the threshold, not the parallel
            continue
        assert (solve_ray(QUIET_DAYTIME, a, b, 10e6) is not None) is closes


def test_line_of_sight_ray_is_not_a_reflection():
    """At 16 MHz and 963 km only the nearly straight topside ray exists; it must not count."""
    a, b = llh_to_ecef(30.0, -150.0, ALT), llh_to_ecef(30.0, -140.0, ALT)
    assert solve_ray(QUIET_DAYTIME, a, b, 16e6) is None
    assert solve_ray(QUIET_DAYTIME, a, b, 10e6) is not None


# -- first order against the shifted sphere ---------------------------------------------------------

@pytest.mark.parametrize("link", [NS, DIAG])
def test_first_order_sensitivities_match_the_exact_forward_model(link):
    a, b, va, vb = link
    s = link_sensitivity(QUIET_DAYTIME, LAT0, LON0, a, b, va, vb, 10e6)
    lev = solve_ray(QUIET_DAYTIME, a, b, 10e6)
    # Central differences cancel the second-order terms (a tilt shortens every echo
    # by an amount quadratic in the tilt).  Tolerances are set against the largest
    # coefficient of the row, the scale at which the differences are resolved.
    tilt_scale = max(abs(s.d_rate[1]), abs(s.d_rate[2]))
    group_scale = max(abs(s.d_group[1]), abs(s.d_group[2]), 1.0)
    for k, (field, step) in enumerate((("w_mps", 1.0), ("g_east", 1e-3), ("g_north", 1e-3), ("shift_m", 10.0))):
        fw = [link_forward(QUIET_DAYTIME, PatchState(LAT0, LON0, **{field: sg * step}), a, b, va, vb, 10e6,
                           r_ref=lev.turning_radius, k_guess=lev.k) for sg in (1.0, -1.0)]
        got = (fw[0][1] - fw[1][1]) / (2 * step)
        if k in (1, 2):
            assert abs(got - s.d_rate[k]) < 2e-3 * tilt_scale
        else:
            assert got == pytest.approx(s.d_rate[k], rel=0.02, abs=2e-3)
        if k > 0:
            gg = (fw[0][0] - fw[1][0]) / (2 * step)
            tol = 2e-2 * group_scale if k in (1, 2) else 0.02 * abs(s.d_group[k]) + 0.02
            assert abs(gg - s.d_group[k]) < tol


# -- estimator ------------------------------------------------------------------------------------

def test_estimator_recovers_drift_tilt_and_profile_from_exact_data():
    from tss.l3_signal.novel.network_tilt_estimator import (
        PROFILE_PARAMS,
        LinkMeasurement,
        estimate_patch,
    )
    truth_prof = replace(QUIET_DAYTIME, scale=1.02, top_scale=0.97)
    st = PatchState(LAT0, LON0, w_mps=30.0, g_east=math.tan(math.radians(0.5)), g_north=-math.tan(math.radians(0.3)),
                    shift_m=200.0)
    links = [(NS, 10e6), (DIAG, 10e6), (_link(20.0, -144.8, 30.0, 20.0, -135.2, 150.0), 10e6),
             (_link(20.0, -147.2, 60.0, 20.0, -132.8, 120.0), 12e6)]
    a0 = llh_to_ecef(LAT0, LON0, ALT)
    v0 = _vel(LAT0, LON0, 30.0)
    links += [((a0, a0, v0, v0), f) for f in (4e6, 6e6, 8e6, 9e6)]
    meas = []
    for (a, b, va, vb), f in links:
        lev = solve_ray(truth_prof, a, b, f)
        fw = link_forward(truth_prof, st, a, b, va, vb, f, r_ref=lev.turning_radius, k_guess=lev.k or None)
        meas.append(LinkMeasurement(LinkGeometry(a, b, va, vb), f, fw[0], fw[1]))
    est = estimate_patch(QUIET_DAYTIME, LAT0, LON0, meas, profile_params=PROFILE_PARAMS)
    assert est.w_mps == pytest.approx(30.0, abs=0.05)
    assert math.degrees(est.g_east) == pytest.approx(0.5, abs=1e-3)
    assert math.degrees(est.g_north) == pytest.approx(-0.3, abs=1e-3)
    assert est.profile_corrections[0] == pytest.approx(0.02, abs=2e-4)
    assert est.profile_corrections[1] == pytest.approx(-0.03, abs=2e-4)


# -- the link table -------------------------------------------------------------------------------

@pytest.fixture(scope="module")
def table():
    from tss.l1_propagation.link_table import LinkTable
    return LinkTable.build(QUIET_DAYTIME, 6371e3 + ALT, [6e6, 10e6, 12e6])


@pytest.mark.parametrize("link,f,patch", [(NS, 10e6, (LAT0, LON0)), (DIAG, 10e6, (21.0, -139.0)),
                                          (_link(20.0, -147.2, 60.0, 20.0, -132.8, 120.0), 12e6, (LAT0, LON0))])
def test_table_matches_ray_by_ray(table, link, f, patch):
    from tss.l1_propagation.link_table import rows_for_patch
    a, b, va, vb = link
    ev = table.evaluate(f, a[None], b[None], va[None], vb[None])
    assert ev["ok"][0]
    grp, rate = rows_for_patch(ev, [0], *patch)
    s = link_sensitivity(QUIET_DAYTIME, *patch, a, b, va, vb, f)
    assert ev["group"][0] == pytest.approx(s.group_path_m, abs=1.0)
    assert ev["rate"][0] == pytest.approx(s.rate_mps, abs=0.05)
    for k in range(4):
        scale = max(abs(s.d_rate[k]), 1e-2 * abs(s.d_rate[1]) + 1e-2)
        assert abs(rate[0, k] - s.d_rate[k]) < 0.01 * scale + 1e-3
        gscale = max(abs(s.d_group[k]), 1.0)
        assert abs(grp[0, k] - s.d_group[k]) < 0.01 * gscale + 0.05


# -- the study ------------------------------------------------------------------------------------

def test_patch_grid_is_equal_area_and_round_trips():
    from tss.studies.route_two_network import EARTH_AREA_KM2, PatchGrid
    g = PatchGrid(2.0)
    areas = g.area_km2(np.arange(g.n))
    assert areas.sum() == pytest.approx(EARTH_AREA_KM2, rel=1e-9)
    assert areas[g.n // 2] == pytest.approx(222.4**2, rel=0.03)
    rng = np.random.default_rng(1)
    lat = np.degrees(np.arcsin(rng.uniform(-1, 1, 500)))
    lon = rng.uniform(-180, 180, 500)
    pid = g.index(lat, lon)
    clat, clon = g.centre(pid)
    assert np.all(np.abs(clat - lat) <= 1.0 + 1e-9)
    assert np.all(g.index(clat, clon) == pid)


def test_walker_states_are_circular():
    from tss.studies.route_two_network import MU, OMEGA_E, R_EARTH_M, walker
    pos, vel = walker(8, 16, 86.4, 1000e3, 180.0, 4, 1234.0)
    r = np.linalg.norm(pos, axis=1)
    assert np.allclose(r, R_EARTH_M + 1000e3)
    vin = vel + np.cross([0.0, 0.0, OMEGA_E], pos)
    assert np.allclose(np.linalg.norm(vin, axis=1), math.sqrt(MU / (R_EARTH_M + 1000e3)), rtol=1e-9)
    assert np.allclose(np.sum(pos * vin, axis=1), 0.0, atol=1e-3)


def test_one_pair_cannot_separate_drift_from_tilt(table):
    """Every link sees w + g . (velocity of its reflection point): one pair alone is rank deficient,
    while two pairs crossing the same patch with different motion resolve it."""
    from tss.studies.route_two_network import PatchGrid, _sigmas, snapshot
    g = PatchGrid(2.0)
    a, b, va, vb = NS
    out = snapshot(table, np.array([a, b]), np.array([va, vb]), g)
    obl = {p: e for p, e in out["fim"].items() if e["n_obl"]}
    assert obl and all(not np.isfinite(_sigmas(e["all"])[0]) for e in obl.values())
    # Add a second pair crossing at the same reflection point, moving the other way.
    c, d, vc, vd = _link(20.0, -144.8, 0.0, 20.0, -135.2, 180.0)
    out = snapshot(table, np.array([a, b, c, d]), np.array([va, vb, vc, vd]), g)
    best = min((_sigmas(e["all"] + np.diag([0, 0, 0, 0, 0.03**-2, 0.05**-2]))[0]
                for e in out["fim"].values() if e["n_obl"] >= 2), default=np.inf)
    assert np.isfinite(best)


def test_haselgrove_truth_for_a_vertical_link_if_available():
    from tss.l1_propagation import haselgrove_links as hl
    if not hl.available():
        pytest.skip("nisc-verify is not on this machine")
    import datetime as dt
    _, chap, *_ = hl._import_nisc()
    t0 = dt.datetime(2019, 9, 21, 12)
    iono = hl.TiltedDriftingIonosphere(chap.quiet_daytime(), LAT0, LON0, 0.0, 0.0, 0.0, t0)
    v = _vel(LAT0, LON0, 30.0)
    tr = hl.trace_link(iono, t0, (LAT0, LON0, ALT), (LAT0, LON0, ALT), v, v, 6e6, (-89.9, 0.0), None)
    a = llh_to_ecef(LAT0, LON0, ALT)
    assert solve_ray(QUIET_DAYTIME, a, a, 6e6).group_path_m == pytest.approx(tr.group_path_m, abs=1.0)
