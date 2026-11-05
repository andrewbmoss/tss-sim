"""WP4 Phase 5: the ISIS-1 retrieval study (tss.studies.isis1_retrieval)."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np

from tss.studies import isis1_retrieval as ir


def _prof(alt_km, ne_m3, lat=40.0, lt=12.0, rz=100.0):
    a, n = np.asarray(alt_km, float), np.asarray(ne_m3, float)
    return SimpleNamespace(lat_deg=lat, extra={"lt_hour": lt, "rz12": rz},
                           ascending=lambda: (a * 1e3, n))


def test_span_measures_error_and_reach():
    prof = _prof([400, 500, 600, 700], [4e11, 2e11, 1e11, 5e10])
    ret = SimpleNamespace(alt_m=np.array([450e3, 550e3, 650e3, 700e3]),
                          ne=np.array([2.828e11, 1.414e11, 0.707e11, 0.5e11]) * 1.1)
    err, reach = ir.span_measures(ret, prof)
    assert abs(err - 0.1) < 0.01
    assert abs(reach - 250.0 / 300.0) < 1e-9


def test_span_measures_without_retrieval():
    prof = _prof([400, 500, 600], [3e11, 2e11, 1e11])
    err, reach = ir.span_measures(None, prof)
    assert np.isnan(err) and reach == 0.0


def test_stratum_labels():
    assert ir.stratum(_prof([400], [1e11], lat=-10, lt=2, rz=30)) == "low/night/low_sun"
    assert ir.stratum(_prof([400], [1e11], lat=70, lt=12, rz=120)) == "high/day/high_sun"


def test_block_counts_failures_and_verdicts():
    rows = [{"error": "", "sci05": 0.05, "sci05_verdict": "goal", "sci07_m": 5e3, "sci07_verdict": "goal",
             "ins10_verdict": "goal", "run_s": 1.0},
            {"error": "ValueError: x", "run_s": 0.5}]
    b = ir.block(rows)
    assert b["n"] == 2 and b["n_failed"] == 1
    assert b["sci05_verdicts"] == {"goal": 1}
    assert b["sci05_median"] == 0.05


def test_thresholds_match_the_requirements_file():
    from pathlib import Path

    import yaml

    root = Path(__file__).resolve().parents[2]
    req = yaml.safe_load((root / "config" / "requirements" / "requirements.yaml").read_text(encoding="utf-8"))
    flat = req.get("requirements", req)
    assert (flat["SCI-05"]["threshold"], flat["SCI-05"]["goal"]) == ir.SCI05
    assert (flat["SCI-07"]["threshold"], flat["SCI-07"]["goal"]) == ir.SCI07_M


def test_verdict():
    assert ir.verdict(0.05, ir.SCI05) == "goal"
    assert ir.verdict(0.12, ir.SCI05) == "threshold"
    assert ir.verdict(0.2, ir.SCI05) == "fail"
    assert ir.verdict(float("nan"), ir.SCI05) == "not_evaluated"


def test_span_scores_ignore_the_cap():
    alt = np.arange(200e3, 1001e3, 5e3)
    truth_ne = 1e12 * np.exp(-(alt - 300e3) / 150e3)
    truth_ne[alt < 300e3] = 1e12 * np.exp(-(300e3 - alt[alt < 300e3]) / 50e3)
    truth = SimpleNamespace(alt_m=alt, ne=truth_ne)
    ret_ne = truth_ne.copy()
    ret_ne[alt < 400e3] *= 1.5            # wrong only below the measured span
    ret = SimpleNamespace(alt_m=alt, ne=ret_ne, flags={})
    s05, s07 = ir.span_scores(ret, truth, 400e3, 1000e3)
    assert s05 < 1e-6 and s07 < 1.0


def test_m2_rule():
    good = {"error": "", "sci05_span": 0.08, "sci07_span_m": 8e3, "sci05": 0.3, "sci07_m": 30e3}
    bad = {"error": "", "sci05_span": 0.30, "sci07_span_m": 30e3, "sci05": 0.3, "sci07_m": 30e3}
    out = ir.m2_rule([good] * 8 + [bad] * 2)
    assert out["span/SCI-05"]["fraction_at_threshold_or_better"] == 0.8
    assert out["pass"] is True
    assert out["full/SCI-05"]["pass"] is False
    assert ir.m2_rule([good] * 7 + [bad] * 3)["pass"] is False
