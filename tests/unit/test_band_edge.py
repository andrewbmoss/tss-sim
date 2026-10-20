"""The band-edge study (min_frequency): the MinFreqStudy harness, ported.

The first group is the harness's own 18 checks (MinFreqStudy/test_minfreq_study.py),
run against the analytic stand-in, so they verify the sampling, the metrics and
the percentile logic and say nothing about IRI's topside.  Then the checks the
port adds: the ensemble direction, the orbit sampling against a fixture drawn
with the harness's own orbits.py, the plan loader, and the study end to end on
the stand-in.  The ``pharlap``-marked test at the end reruns two cases of the
Rev G run through MATLAB and PHaRLAP and holds them to the reference summary.
"""

from __future__ import annotations

import csv
import datetime as dt
import json
from pathlib import Path

import numpy as np
import pytest

from tss.config import repo_root
from tss.l0_ionosphere.backends import EXTRA_COLUMNS, AnalyticBackend
from tss.platform import orbits
from tss.studies import min_frequency as ms

ALT = np.arange(150e3, 1000e3 + 1, 10e3)
REF = repo_root() / "data" / "reference" / "band_edge"


def _case(alt=800e3, inc=80.0, n=4000, f107=120.0):
    orb = orbits.CircularOrbit(alt, inc)
    b = AnalyticBackend(f107=f107)
    return ms.run_case(b, orb, dt.datetime(2024, 1, 1), n, ALT, seed=3)


# --- orbit sampling ---------------------------------------------------------
def test_latitude_bounded_by_inclination():
    for inc in (30.0, 80.0):
        s = orbits.CircularOrbit(800e3, inc).sample(20000, seed=1)
        assert np.abs(s["lat_deg"]).max() <= inc + 1e-6


def test_latitude_density_peaks_at_turning_points():
    s = orbits.CircularOrbit(800e3, 80.0).sample(200000, seed=1)
    a = np.abs(s["lat_deg"])
    assert np.mean((a > 70) & (a <= 80)) > 1.5 * np.mean(a <= 10)


def test_local_time_is_swept_uniformly_over_a_year():
    s = orbits.CircularOrbit(800e3, 80.0).sample(200000, seed=1)
    h, _ = np.histogram(s["local_time_hours"], bins=12, range=(0, 24))
    frac = h / h.sum()
    assert frac.max() < 0.11 and frac.min() > 0.06, frac


def test_sun_synchronous_orbit_does_not_sweep_local_time():
    sso = orbits.CircularOrbit(800e3, orbits.sun_synchronous_inclination(800e3))
    assert abs(sso.lt_drift_hours_per_day) < 1e-3
    assert sso.days_for_full_lt_coverage > 3650


def test_orbit_classes_sweep_local_time_within_the_mission():
    for k, orb in orbits.orbit_classes(800e3).items():
        assert orb.days_for_full_lt_coverage < 365.25, (k, orb.days_for_full_lt_coverage)


def test_nodal_precession_sign_and_magnitude():
    assert orbits.CircularOrbit(800e3, 30.0).raan_rate < 0
    assert orbits.CircularOrbit(800e3, 120.0).raan_rate > 0
    assert 95.0 < orbits.sun_synchronous_inclination(800e3) < 102.0


# --- metric definitions -----------------------------------------------------
def test_coverage_is_the_lower_tail_not_the_upper():
    fn = np.linspace(0.5e6, 3.0e6, 10001)
    f95 = ms.f_min_for_coverage(fn, 0.95)
    assert abs(ms.coverage_fraction(fn, f95) - 0.95) < 0.01
    assert f95 < np.median(fn)


def test_tighter_coverage_demands_a_lower_f_min():
    _, _, p = _case()
    fn = p["fn_sat"]
    assert ms.f_min_for_coverage(fn, 0.90) > ms.f_min_for_coverage(fn, 0.95) > \
        ms.f_min_for_coverage(fn, 0.99)


def test_completeness_is_unity_when_f_min_is_below_the_local_plasma_frequency():
    _, ne, p = _case()
    c = ms.completeness(1.0, p, ne, ALT, 800e3)
    assert np.allclose(c["altitude_completeness"], 1.0)
    assert np.allclose(c["density_completeness"], 1.0)


def test_completeness_falls_monotonically_as_f_min_rises():
    _, ne, p = _case()
    prev_a = prev_d = 1.0
    for fm in (0.3e6, 1.0e6, 2.0e6, 4.0e6, 8.0e6):
        c = ms.completeness(fm, p, ne, ALT, 800e3)
        a = float(np.mean(c["altitude_completeness"]))
        d = float(np.mean(c["density_completeness"]))
        assert a <= prev_a + 1e-9 and d <= prev_d + 1e-9, (fm, a, d)
        prev_a, prev_d = a, d


def test_completeness_is_zero_above_the_critical_frequency():
    _, ne, p = _case()
    c = ms.completeness(float(np.max(p["fof2"])) * 1.05, p, ne, ALT, 800e3)
    assert np.allclose(c["altitude_completeness"], 0.0)


def test_start_altitude_lies_between_the_peak_and_the_spacecraft():
    _, ne, p = _case()
    c = ms.completeness(1.5e6, p, ne, ALT, 800e3)
    ok = np.isfinite(c["h_start"])
    assert (c["h_start"][ok] <= 800e3 + 1.0).all()
    assert (c["h_start"][ok] >= p["hmf2"][ok] - 1e4).all()


# --- physical direction of the answer --------------------------------------
def test_higher_orbit_lowers_the_required_minimum_frequency():
    vals = [ms.f_min_for_coverage(_case(alt=a)[2]["fn_sat"], 0.95)
            for a in (500e3, 600e3, 800e3, 900e3)]
    assert all(a > b for a, b in zip(vals, vals[1:], strict=False)), vals


def test_solar_maximum_raises_the_required_minimum_frequency():
    lo = _case(f107=70.0)[2]["fn_sat"]
    hi = _case(f107=180.0)[2]["fn_sat"]
    assert ms.f_min_for_coverage(hi, 0.95) > ms.f_min_for_coverage(lo, 0.95)


def test_low_inclination_orbit_sees_a_denser_topside():
    assert np.median(_case(inc=30.0)[2]["fn_sat"]) > np.median(_case(inc=80.0)[2]["fn_sat"])


class _ExtraBackend(AnalyticBackend):
    """Analytic profiles plus a per-call ``extra``, like the MATLAB backend."""

    def profiles(self, lat_deg, lon_deg, when, alt_grid_m, **opts):
        ne = super().profiles(lat_deg, lon_deg, when, alt_grid_m, **opts)
        n = ne.shape[0]
        self.extra = np.zeros((n, len(EXTRA_COLUMNS)))
        self.extra[:, 0] = ne.max(axis=1)
        self.extra[:, 1] = 300.0 + np.arange(n) % 7
        return ne


def test_extra_is_collected_across_chunks_not_overwritten():
    n, chunk = 5000, 2000
    orb = orbits.CircularOrbit(800e3, 80.0)
    _, ne, p = ms.run_case(_ExtraBackend(f107=120.0), orb, dt.datetime(2024, 1, 1), n, ALT,
                           seed=3, chunk=chunk)
    assert ne.shape[0] == n and p["hmf2"].shape[0] == n and p["extra"].shape[0] == n
    want = np.concatenate([300.0 + np.arange(min(chunk, n - a)) % 7 for a in range(0, n, chunk)])
    assert np.allclose(p["hmf2"] / 1e3, want)


def test_mismatched_extra_is_rejected_not_broadcast():
    ne = np.full((10, len(ALT)), 1e11)
    with pytest.raises(ValueError, match="rows"):
        ms.profile_parameters(ne, ALT, 800e3, extra=np.zeros((4, 10)))


def test_backend_records_provenance():
    assert "NOT a science result" in AnalyticBackend().provenance


# --- what the port adds -------------------------------------------------------
def test_the_ensemble_takes_the_minimum_not_the_maximum():
    """The direction error that gave 1.9 MHz and 53.4 per cent coverage."""
    assert ms.combine_ensemble(np.array([0.334, 0.477, 0.354])) == 0.334
    rows = [dict(orbit_class="O1", altitude_km=900.0, solar="min", topside=t,
                 f_min_for_95pct_coverage_MHz=v, f_min_for_99pct_coverage_MHz=v / 2)
            for t, v in (("NeQuick", 0.334), ("IRI-2001", 0.477), ("IRI-2001 corrected", 0.354))]
    (e,) = ms.ensemble_table(rows)
    assert e["f_min_95_ensemble_MHz"] == 0.334 and e["governing_model_95"] == "NeQuick"
    assert e["f_min_95_unsafe_max_MHz"] == 0.477


def test_orbit_sampling_draws_what_the_harness_drew():
    """tss.platform.orbits against MinFreqStudy/orbits.py, same seed, same samples."""
    fix = json.loads((REF / "orbit_sample_fixture.json").read_text())
    idx = [0, 1, 2, 9999, 19999]
    for cls, orb in orbits.orbit_classes(900e3).items():
        s = orb.sample(20000, seed=17)
        for key in ("t_seconds", "lat_deg", "lon_deg", "local_time_hours"):
            assert np.allclose(s[key][idx], fix[cls][key], rtol=0, atol=1e-9), (cls, key)


def test_the_reference_plan_is_the_rev_g_run():
    plan = ms.load_plan({})
    assert plan["name"] == "band_edge_reference"
    assert plan["altitudes_km"] == [500, 600, 700, 800, 900]
    assert plan["samples"] == 20000 and plan["seed"] == 17
    assert [s["r12"] for s in plan["solar"]] == [10, 120]
    assert len(plan["topsides"]) == 3


def test_a_plan_whose_grid_stops_below_the_spacecraft_is_refused():
    with pytest.raises(ValueError, match="grid stops"):
        ms.load_plan({"min_frequency": {"altitudes_km": [1400], "alt_grid_km": [150, 1000, 10]}})


def test_the_reference_summary_is_complete():
    with open(REF / "minfreq_summary.csv", newline="") as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == 60
    assert tuple(rows[0].keys()) == ms.SUMMARY_COLUMNS
    gov = [r for r in rows if r["orbit_class"] == "O1" and r["altitude_km"] == "900.0"
           and r["solar"] == "min" and r["topside"] == "NeQuick"]
    assert abs(float(gov[0]["f_min_for_95pct_coverage_MHz"]) - 0.334) < 0.001


def test_the_study_runs_end_to_end_on_the_stand_in(tmp_path: Path):
    from tss.provenance import Manifest

    config = {"min_frequency": {"backend": "analytic", "samples": 400, "altitudes_km": [600, 900],
                                "topsides": ["NeQuick"]}}
    result = ms.run(tmp_path, config=config, manifest=Manifest(study="min_frequency", config=config))
    with open(tmp_path / "minfreq_summary.csv", newline="") as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == 2 * 2 * 2       # altitudes x classes x solar cases
    assert set(ms.SUMMARY_COLUMNS) <= set(rows[0])
    z = np.load(tmp_path / rows[0]["samples_file"])
    assert z["fn_sat_hz"].shape == (400,)
    assert result.metrics["combination"] == "minimum of per-model 5th percentiles"
    assert (tmp_path / "ensemble.csv").is_file()


# --- through MATLAB and PHaRLAP -------------------------------------------------
@pytest.mark.pharlap
@pytest.mark.parametrize("cls,alt_km,solar,r12,topside", [
    ("O1", 900, "min", 10.0, "NeQuick"),
    ("O2", 500, "max", 120.0, "IRI-2001"),
])
def test_two_cases_reproduce_the_rev_g_run(tmp_path, cls, alt_km, solar, r12, topside):
    with open(REF / "minfreq_summary.csv", newline="") as fh:
        ref = [r for r in csv.DictReader(fh) if r["orbit_class"] == cls
               and float(r["altitude_km"]) == alt_km and r["solar"] == solar
               and r["topside"] == topside][0]
    plan = ms.load_plan({})
    backend = ms.build_backend(plan, topside, {"name": solar, "mode": "r12", "r12": r12},
                               tmp_path / "work", probe=True)
    orb = orbits.orbit_classes(alt_km * 1e3)[cls]
    grid = np.arange(150e3, 1000e3 + 1, 10e3)
    s, ne, p = ms.run_case(backend, orb, dt.datetime(2015, 1, 1), 20000, grid, seed=17)
    summ, _ = ms.summarize(p, ne, grid, alt_km * 1e3)
    for key in ("f_min_for_90pct_coverage_MHz", "f_min_for_95pct_coverage_MHz",
                "f_min_for_99pct_coverage_MHz", "fn_sat_median_MHz", "fof2_median_MHz",
                "hmf2_median_km"):
        assert summ[key] == pytest.approx(float(ref[key]), rel=1e-4), key


def _chapman_case(n_scale, h_sat, alt):
    """Profiles whose density at h_sat scales by ``n_scale`` per sample."""
    from tss.studies.min_frequency import profile_parameters

    base = 1e12 * np.exp(-(alt - 300e3) / 150e3)
    base[alt < 300e3] = 1e12 * np.exp(-((alt[alt < 300e3] - 300e3) / 60e3) ** 2)
    ne = np.asarray(n_scale)[:, None] * base[None, :]
    return ne, profile_parameters(ne, alt, h_sat)


def test_goal_metrics_completeness_and_transition():
    """Phase 4 additions: completeness at a fixed start and the O+/H+ transition."""
    from tss.l0_ionosphere.backends import EXTRA_COLUMNS
    from tss.studies.min_frequency import INS01A_THRESHOLD_HZ, goal_metrics

    alt = np.arange(150e3, 2001e3, 10e3)
    h_sat = 1400e3
    ne, params = _chapman_case([0.25, 0.8, 1.25, 4.0], h_sat, alt)
    f_start = float(np.median(params["fn_sat"]))
    extra = np.full((4, len(EXTRA_COLUMNS)), np.nan)
    k = list(EXTRA_COLUMNS).index("h_oplus_hplus_km")
    extra[:, k] = [800.0, 1200.0, 1600.0, np.nan]
    params["extra"] = extra
    g = goal_metrics(params, ne, alt, h_sat, f_start_hz=f_start)
    assert g["coverage_at_ins01a"] == pytest.approx(0.5)
    assert 0.0 < g["mean_alt_completeness_at_ins01a"] < 1.0
    assert g["frac_oh_transition_valid"] == pytest.approx(0.75)
    assert g["frac_oh_transition_below_hsat"] == pytest.approx(0.5)   # 800 and 1200 km
    assert g["oh_transition_p50_km"] == pytest.approx(1200.0)
    # A start below every fN_sat retrieves the whole topside.
    g0 = goal_metrics(params, ne, alt, h_sat, f_start_hz=1.0)
    assert g0["coverage_at_ins01a"] == 1.0
    assert g0["p05_alt_completeness_at_ins01a"] == pytest.approx(1.0)
    assert INS01A_THRESHOLD_HZ == 0.33e6


def test_goal_metrics_without_transition_column():
    from tss.studies.min_frequency import goal_metrics

    alt = np.arange(150e3, 1001e3, 10e3)
    ne, params = _chapman_case([1.0, 2.0], 900e3, alt)
    params["extra"] = None
    g = goal_metrics(params, ne, alt, 900e3)
    assert np.isnan(g["oh_transition_p50_km"])
    assert np.isnan(g["frac_oh_transition_below_hsat"])
