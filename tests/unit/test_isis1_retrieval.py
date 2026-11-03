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
