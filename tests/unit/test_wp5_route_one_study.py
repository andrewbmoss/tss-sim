"""WP5 Phase 3: the study ``wp5_route_one``, route one at signal level over the Phase 0 scenes.

The layout is the one Andrew approved on 2026-11-27 (all scenes noise free; nine
scenes against SNR and fading; the oscillator's goal; the halo cases; orbit
scenes against the prior's scale height).  These tests run single small jobs;
the study itself takes hours.
"""

from __future__ import annotations

import io
import json
import math

import numpy as np
import pytest

from tss.registry import build
from tss.studies import STUDIES, run_study
from tss.studies import wp5_route_one as W
from tss.studies.wp5_campaign import SCENES

FAST_SCENE = "A_free_10MHz_500km_x0_v25"
BY_ID = {s.scene_id: s for s in SCENES}


def test_the_study_is_registered():
    assert STUDIES["wp5_route_one"].status == "implemented"


def test_the_default_layout_is_the_approved_one():
    s = W.settings({})
    jobs = W.jobs(s)
    count = {b: sum(j.block == b for j in jobs) for b in W.BLOCKS}
    assert count == {"scenes": len(SCENES), "monte_carlo": 9 * 4 * 100, "oscillator": 9 * 100,
                     "halo": 6 * 100, "orbit": 3 * 100}
    assert len({j.job_id for j in jobs}) == len(jobs)
    assert s.levels == (0.0, 5.0, 10.0, 15.0, 20.0, "link") and s.tau_c_s == (None, 0.1, 0.5, 1.0)
    assert s.oscillator_main == 5.8e-8 and s.oscillator_second == 3.2e-9
    assert s.halo_cases == {"none": None, "K10": 10.0, "K3": 3.0, "K0": 0.0} and s.halo_ratios == (2.0, 4.0)
    assert {(j.halo, j.halo_ratio) for j in jobs if j.block == "halo"} == {
        ("none", None), ("K10", 2.0), ("K10", 4.0), ("K3", 2.0), ("K3", 4.0), ("K0", 2.0)}
    assert s.halo_levels == (10.0, 20.0) and s.orbit_scale_heights_km == (100.0, 200.0, 300.0)
    assert s.estimators == ("pinned", "joint", "joint_fading")


def test_the_monte_carlo_scenes_are_the_seven_environments_and_two_budget_points():
    s = W.settings({})
    got = [BY_ID[i] for i in s.mc_scenes]
    assert sum(g.group == "A_untilted" and g.x_m == 0.0 and g.drift_mps == 25.0 for g in got) == 7
    assert sum(g.group == "E_budget" for g in got) == 2 and len(got) == 9
    assert len({g.environment for g in got}) == 9


@pytest.mark.parametrize("scene_id", ["A_plasma_5MHz_mu0.5_x5_v25", "B_free_5MHz_500km_t1_v25",
                                      "C_plasma_5MHz_mu0.5_t0.5_v25", "D_free_5MHz_500km_L32",
                                      "S_free_dip10_psi15_phi5",
                                      pytest.param("S_plasma_dip10_psi0_phi5", marks=pytest.mark.slow),
                                      "E_plasma_5MHz_P500km_v25"])
def test_every_kind_of_scene_becomes_a_facet_the_tier_builds(scene_id):
    setup = W.scene_setup(BY_ID[scene_id])
    tier = build("l1", "facet", facets=[setup.facet], **setup.l1_kwargs)
    assert tier.facets and setup.t_s.size > 20
    assert np.all(np.diff(setup.t_s) == pytest.approx(1 / 45))
    if BY_ID[scene_id].group == "C_extended":
        assert setup.facet["taper"] == "tukey" and setup.facet["length_m"] == 80e3
    if BY_ID[scene_id].group == "S_sheet":
        assert setup.facet["kind"] == "sheet" and setup.facet["velocity_u_mps"] == 30.0


def test_a_noise_free_scene_reads_its_drift():
    s = W.settings({})
    rows = W.run_job(W.Job("scenes", FAST_SCENE), s, {})
    echo = [r for r in rows if r["facet_echo"]]
    assert [r["estimator"] for r in echo] == list(s.estimators)
    for r in echo:
        assert set(r) == set(W.ROW_FIELDS)
        assert r["detected"] and r["formed"] and abs(r["err_mps"]) < 1.0
        assert r["level_db"] == "none" and r["oscillator"] == "none" and r["halo"] == "none"
        assert abs(r["param_err_mps"]) < 1.0 and math.isfinite(r["peak_contrast"])


@pytest.fixture(scope="module")
def fast_settings():
    return W.settings({"wp5_route_one": {"seeds": 1, "mc_scenes": [FAST_SCENE], "levels": [20.0, "link"],
                                         "tau_c_s": [None], "halo_points": 12, "halo_levels": [20.0]}})


@pytest.fixture(scope="module")
def cache(fast_settings):
    return {key: W.prepare(key, fast_settings) for key in W.needed_caches(
        [j for j in W.jobs(fast_settings) if j.block in ("monte_carlo", "oscillator")], fast_settings)}


def test_a_monte_carlo_job_sets_the_snr_and_records_the_measures(fast_settings, cache):
    job = next(j for j in W.jobs(fast_settings) if j.block == "monte_carlo")
    rows = [r for r in W.run_job(job, fast_settings, cache) if r["facet_echo"]]
    assert {(r["level_db"], r["estimator"]) for r in rows} == {
        (lv, e) for lv in (20.0, "link") for e in fast_settings.estimators}
    for r in rows:
        assert r["oscillator"] == "threshold" and math.isnan(r["tau_c_s"])
        assert r["detected"]
        assert math.isfinite(r["snr_db"]) and math.isfinite(r["edge_power_over_noise"])
        assert math.isfinite(r["err_mps"]) == r["formed"]
    at20 = [r["snr_db"] for r in rows if r["level_db"] == 20.0]
    assert np.mean(at20) == pytest.approx(20.0, abs=1.5)


def test_the_same_job_gives_the_same_rows(fast_settings, cache):
    job = next(j for j in W.jobs(fast_settings) if j.block == "oscillator")
    a, b = W.run_job(job, fast_settings, cache), W.run_job(job, fast_settings, cache)
    assert json.dumps(a, default=str) == json.dumps(b, default=str)
    assert a[0]["oscillator"] == "goal"


def test_a_halo_job_carries_its_case_and_length(fast_settings, monkeypatch):
    jobs = [j for j in W.jobs(fast_settings) if j.block == "halo"]
    assert [(j.halo, j.halo_ratio) for j in jobs] == [
        ("none", None), ("K10", 2.0), ("K10", 4.0), ("K3", 2.0), ("K3", 4.0), ("K0", 2.0)]
    got = W.needed_caches(jobs, fast_settings)
    cache = {key: W.prepare(key, fast_settings) for key in got}
    built, straight = [], W._straight
    monkeypatch.setattr(W, "_straight", lambda setup: built.append(setup.facet["halo"]) or straight(setup))
    rows = [r for r in W.run_job(jobs[4], fast_settings, cache) if r["facet_echo"]]
    assert [(h["k"], h["ratio"]) for h in built] == [(3.0, 4.0)]
    assert {(r["halo"], r["halo_ratio"]) for r in rows} == {("K3", 4.0)}
    assert {r["level_db"] for r in rows} == {20.0}


def test_the_halo_lengths_are_set_from_the_configuration():
    s = W.settings({"wp5_route_one": {"halo_ratios": [2, 1], "seeds": 1}})
    assert s.halo_ratios == (2.0, 1.0)
    assert [j.halo_ratio for j in W.jobs(s) if j.halo == "K3"] == [2.0, 1.0]
    assert [j.halo_ratio for j in W.jobs(s) if j.halo == "K0"] == [2.0]
    with pytest.raises(ValueError, match="halo_ratios"):
        W.settings({"wp5_route_one": {"halo_ratios": []}})


@pytest.mark.slow
def test_an_orbit_job_sweeps_the_scale_height():
    s = W.settings({"wp5_route_one": {"seeds": 1, "orbit_levels": [20.0]}})
    job = next(j for j in W.jobs(s) if j.block == "orbit" and "5MHz" in j.scene_id)
    cache = {key: W.prepare(key, s) for key in W.needed_caches([job], s)}
    rows = [r for r in W.run_job(job, s, cache) if r["facet_echo"] and r["estimator"] == "pinned"]
    assert sorted(r["h_assumed_km"] for r in rows) == [100.0, 200.0, 300.0]
    assert all(r["formed"] and abs(r["err_mps"]) < 15.0 for r in rows)


def test_the_summary_groups_the_facet_echoes_and_counts_the_formed():
    base = {n: math.nan for n in W.ROW_FIELDS}
    base.update({"block": "monte_carlo", "scene_id": FAST_SCENE, "level_db": 10.0, "tau_c_s": math.nan,
                 "oscillator": "threshold", "halo": "none", "h_assumed_km": math.nan, "estimator": "pinned",
                 "facet_echo": True, "detected": True, "w_s": 0.45, "p_group_m": 500e3, "mu_sat": 1.0,
                 "f_mhz": 10.0, "n_sweeps": 101, "snr_db": 10.0, "echo_class": "compact_facet", "group": "A_untilted",
                 "echo_class_truth": "compact_facet"})
    rows = [{**base, "formed": True, "err_mps": e} for e in (1.0, -2.0, 50.0)]
    rows.append({**base, "formed": False})
    rows.append({**base, "facet_echo": False, "formed": True, "err_mps": 900.0})
    (got,) = W.summarise(rows)
    assert got["n"] == 4 and got["detected"] == 4 and got["formed"] == 3
    assert got["rms_mps"] == pytest.approx(math.sqrt((1 + 4 + 2500) / 3))
    assert got["bound_mps"] > 0 and got["outliers"] == pytest.approx(1 / 3)
    assert got["other_echoes_formed"] == 1
    halo = [{**base, "block": "halo", "halo": "K3", "halo_ratio": r, "formed": True, "err_mps": 1.0}
            for r in (2.0, 4.0)]
    assert [g["halo_ratio"] for g in W.summarise(halo)] == [2.0, 4.0]


def test_a_run_writes_its_outputs_resumes_and_refuses_another_commits_parts(tmp_path):
    overrides = [f'wp5_route_one.scene_ids=["{FAST_SCENE}"]', 'wp5_route_one.blocks=["scenes"]',
                 "wp5_route_one.processes=1"]
    out = tmp_path / "run"
    assert run_study("wp5_route_one", out, overrides, stream=io.StringIO()) == 0
    for name in ("wp5_route_one_echoes.csv.gz", "wp5_route_one_summary.csv", "wp5_route_one.json",
                 "summary.md", "manifest.json"):
        assert (out / name).exists(), name
    part = next((out / "parts").glob("*.json"))
    stamp = part.stat().st_mtime_ns
    assert run_study("wp5_route_one", out, overrides, stream=io.StringIO()) == 0
    assert part.stat().st_mtime_ns == stamp  # reused
    meta = json.loads((out / "parts_meta.json").read_text(encoding="utf-8"))
    meta["commit"] = "0" * 40
    (out / "parts_meta.json").write_text(json.dumps(meta), encoding="utf-8")
    assert run_study("wp5_route_one", out, overrides, stream=io.StringIO()) == 6
