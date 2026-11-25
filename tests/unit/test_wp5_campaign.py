"""WP5 Phase 0: the evaluation harness (tss.studies.wp5_campaign).

The full grids (c3_verify's 107 trials, the 1000-trial budget, the plasma
sheets) run in the study itself and in the reference run
``data/reference/wp5_baseline``.  These tests use reduced grids so the default
suite stays fast.
"""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from tss.studies import wp5_campaign as w


def _by_id():
    return {s.scene_id: s for s in w.SCENES}


def _same(a: dict, b: dict) -> bool:
    """Row equality with NaN equal to NaN (a declined facet has NaN drift)."""
    if a.keys() != b.keys():
        return False
    for k in a:
        x, y = a[k], b[k]
        if isinstance(x, float) and isinstance(y, float) and np.isnan(x) and np.isnan(y):
            continue
        if x != y:
            return False
    return True


def test_scene_list_covers_the_verified_grids():
    groups = {}
    for s in w.SCENES:
        groups.setdefault(s.group, []).append(s)
    # c3_verify's grid: A 7 environments x 2 positions x 4 drifts, B 4 x 4 tilts
    # x 2 drifts, C 2 x 3 tilts x 2 drifts, D 7 lengths.
    assert len(groups["A_untilted"]) == 56
    assert len(groups["B_tilted"]) == 32
    assert len(groups["C_extended"]) == 12
    assert len(groups["D_length"]) == 7
    sheets = groups["S_sheet"]
    assert {s.heading_deg for s in sheets} == {0.0, 15.0, 45.0}
    assert any(s.medium == "plasma" for s in sheets)
    assert {s.f_mhz for s in w.SCENES if s.group == "A_untilted"} == {3.0, 5.0, 10.0}
    drifts = {s.drift_mps for s in w.SCENES if s.group == "A_untilted"}
    assert min(drifts) == -50.0 and max(drifts) == 100.0
    tilts = {s.tilt_deg for s in groups["B_tilted"]}
    assert min(tilts) == 0.1 and max(tilts) == 1.0
    assert all(s.plate_m == 80e3 for s in groups["C_extended"])
    ids = [s.scene_id for s in w.SCENES]
    assert len(ids) == len(set(ids))
    assert {s.echo_class_truth for s in w.SCENES} == {"compact_facet", "extended", "field_aligned"}


def test_scene_seed_depends_only_on_the_scene_and_the_campaign_seed():
    s = w.SCENES[3]
    assert w.scene_seed(s, 7) == w.scene_seed(s, 7)
    assert w.scene_seed(s, 7) != w.scene_seed(s, 8)
    assert w.scene_seed(w.SCENES[3], 7) != w.scene_seed(w.SCENES[4], 7)


def test_a_scene_gives_the_same_row_alone_or_in_a_set():
    by = _by_id()
    pick = [by[i] for i in ("A_free_5MHz_500km_x0_v25", "C_free_5MHz_500km_t0.5_v25")]
    alone = w.run_scene(pick[1], campaign_seed=11)
    together = [w.run_scene(s, campaign_seed=11) for s in pick]
    assert _same(alone, together[1])
    assert np.isnan(alone["drift_mps"])  # the extended plate is declined


def test_reduced_grid_reproduces_each_kind_of_c3_verify_result():
    by = _by_id()
    a = w.run_scene(by["A_plasma_5MHz_mu0.5_x5_v-50"], campaign_seed=1)
    assert a["verdict"] == "pass" and abs(a["drift_err_mps"]) <= 1.0
    assert a["echo_class"] == "compact_facet"
    b = w.run_scene(by["B_free_5MHz_500km_t1_v25"], campaign_seed=1)
    assert b["echo_class"] == "compact_facet" and "TILT_UNFLAGGED" in b["flags"]
    assert 0.99 <= b["drift_err_mps"] / b["predicted_bias_mps"] <= 1.01
    # Free space: the plasma extended plates take about 20 s each (exact rays to
    # every point of an 80 km plate) and run only in the study.
    c = w.run_scene(by["C_free_5MHz_500km_t1_v25"], campaign_seed=1)
    assert c["echo_class"] == "extended" and "TILT_DEGENERATE" in c["flags"]
    assert np.isnan(c["drift_mps"])
    d4 = w.run_scene(by["D_free_5MHz_500km_L4"], campaign_seed=1)
    d1 = w.run_scene(by["D_free_5MHz_500km_L1"], campaign_seed=1)
    assert d1["echo_class"] == "compact_facet" and d4["echo_class"] == "extended"


def test_rows_use_the_contract_names_and_score_against_the_plasma_rate():
    row = w.run_scene(_by_id()["A_free_3MHz_500km_x0_v100"], campaign_seed=1)
    for k in ("scene_id", "echo_class", "echo_class_truth", "truth_mps", "drift_mps", "drift_err_mps",
              "los_rate_mps", "flags", "not_evaluated", "seed"):
        assert k in row
    # Truth is the plasma's own rate along the wave normal, here the vertical drift.
    assert row["truth_mps"] == 100.0
    assert "tau_c_s" in row["not_evaluated"] and "drift_mf_mps" in row["not_evaluated"]


def test_free_sheet_reproduces_the_closed_form_error():
    row = w.run_scene(_by_id()["S_free_dip10_psi15_phi5"], campaign_seed=1)
    assert abs(row["drift_err_mps"] - w.EXPECTED["sheet_assume_err_psi15_phi5_mps"]) < 1.0
    assert abs(row["drift_given_err_mps"]) <= w.EXPECTED["sheet_given_max_abs_err_free_mps"]


def test_expected_values_are_the_tv_numbers():
    e = w.EXPECTED
    assert (e["A_pass"], e["A_n"], e["B_flagged"], e["B_n"], e["C_declined"], e["C_n"]) == (56, 56, 0, 32, 12, 12)
    assert e["budget_total_mps"] == {
        (3.0, 500.0): (12.2, 33.7), (3.0, 1000.0): (11.3, 30.9), (5.0, 500.0): (10.8, 23.1),
        (5.0, 1000.0): (10.4, 20.1), (10.0, 500.0): (10.6, 14.8), (10.0, 1000.0): (10.5, 14.5),
    }


def test_budget_check_holds_to_five_percent():
    exp = w.EXPECTED["budget_total_mps"]
    good = [{"f_mhz": f, "p_km": p, "fading": j, "total_mps": v[j] * 1.04} for (f, p), v in exp.items() for j in (0, 1)]
    assert all(r["within_5pct"] for r in w.check_budget(good))
    bad = [dict(good[0], total_mps=exp[(3.0, 500.0)][0] * 1.06)] + good[1:]
    assert not w.check_budget(bad)[0]["within_5pct"]


def test_run_records_every_seed_in_the_manifest(tmp_path):
    cfg = {"wp5_campaign": {"groups": ["D_length"], "budget_cases": [[5.0, 500.0]], "budget_trials": 20,
                            "plasma_sheets": False, "sheets": False}}
    manifest = SimpleNamespace(seeds={})
    res = w.run(tmp_path, config=cfg, manifest=manifest)
    for key in ("campaign", "scenes", "budget", "budget_trials"):
        assert key in manifest.seeds
    assert len(manifest.seeds["scenes"]) == 7
    assert manifest.seeds["budget"] == {"5MHz_500km_fading0": 20262026, "5MHz_500km_fading1": 20262027}
    assert (tmp_path / "wp5_campaign_scenes.csv").is_file()
    assert (tmp_path / "wp5_campaign_budget.csv").is_file()
    # The plate-length switch reproduces; 20 trials is too few to hold the budget to 5 %.
    assert not any(m.startswith("D:") for m in res.metrics["mismatches"])


@pytest.mark.parametrize("bad", [{"groups": ["Z_nothing"]}])
def test_unknown_group_is_refused(bad, tmp_path):
    with pytest.raises(ValueError, match="group"):
        w.run(tmp_path, config={"wp5_campaign": bad}, manifest=SimpleNamespace(seeds={}))


# --------------------------------------------------------------------------
# Budget-point facets (WP5 P1, decision of 2026-11-24 evening): only the two
# budget points a facet can reach in the reference profile are added.
# --------------------------------------------------------------------------


def test_budget_point_scenes_are_the_two_reachable_points():
    e = [s for s in w.SCENES if s.group == "E_budget"]
    assert [(s.f_mhz, s.p_target_km) for s in e] == [(5.0, 500.0), (10.0, 500.0)]
    for s in e:
        assert s.medium == "plasma" and s.echo_class_truth == "compact_facet"
        assert s.plate_m == 2.0e3 and s.drift_mps == 25.0 and s.x_m == 0.0
    assert e[0].mu_facet == pytest.approx(0.554, abs=0.002)
    assert e[1].mu_facet == pytest.approx(0.870, abs=0.002)
    assert "E_budget" in w.GROUPS


def test_the_other_budget_points_are_out_of_reach_of_the_profile():
    assert w.BUDGET_POINTS_UNREACHABLE == ((3.0, 500.0), (3.0, 1000.0), (5.0, 1000.0), (10.0, 1000.0))
    assert w.max_facet_group_path_m(3.0e6) == pytest.approx(527e3, rel=0.01)
    assert w.max_facet_group_path_m(5.0e6) == pytest.approx(750e3, rel=0.01)
    assert w.max_facet_group_path_m(10.0e6) == pytest.approx(1034e3, rel=0.01)
    for f, p in w.BUDGET_POINTS_UNREACHABLE:
        assert w.budget_point_mu(f * 1e6, p * 1e3) is None
    for f, p in w.BUDGET_POINTS_REACHABLE:
        assert w.budget_point_mu(f * 1e6, p * 1e3) >= w.MIN_MU_FACET


def test_budget_point_scene_runs_at_its_apparent_range_and_is_scored():
    e = [s for s in w.SCENES if s.group == "E_budget"]
    rows = [w.run_scene(s, campaign_seed=7) for s in e]
    for s, r in zip(e, rows, strict=True):
        assert r["group_path_km"] == pytest.approx(s.p_target_km, abs=1.0)
        assert r["verdict"] == "pass" and abs(r["drift_err_mps"]) <= 1.0
    m = w.check_reproduction(rows, [])
    assert (m["E_pass"], m["E_n"]) == (2, 2) and m["reproduced"]
    assert m["E_max_abs_group_path_dev_km"] <= 1.0
