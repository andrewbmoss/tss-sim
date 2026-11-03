"""WP4 Phase 5: the FM-CW acceptance study (tss.studies.fmcw_acceptance)."""

from __future__ import annotations

import numpy as np

from tss.studies import fmcw_acceptance as fa


def test_parabolic_layer_phase_height_integrates_to_group_height():
    phase, group = fa.parabolic_layer(7.2e6, 0.0, 250e3, 100e3)
    f = np.linspace(5.5e6, 6.9e6, 4001)
    fp = f * phase(f, 0.0)
    deriv = np.gradient(fp, f)
    assert np.max(np.abs(deriv[5:-5] - group(f, 0.0)[5:-5])) < 50.0


def test_fc_schedule_gives_constant_group_velocity():
    sched = fa._fc_schedule(6.4e6, 250e3, 100e3, 376e3, -50.0)
    t = np.array([0.0, 100.0, 200.0])
    x = 6.4e6 / sched(t)
    h = 250e3 + 50e3 * x * np.log((1 + x) / (1 - x))
    assert np.allclose(h, 376e3 - 50.0 * t, atol=5.0)


def test_items_score_as_recorded():
    res = fa.compute({"fmcw_acceptance": {"mc_trials": 200}})
    out = {r["item"]: r["outcome"] for r in res["rows"]}
    for item in ("A1", "A2", "A3", "A4a", "A4b", "A4c", "A5b", "A9", "A10", "A11", "A12"):
        assert out[item] == "reproduced", item
    assert out["A5a"] == "difference explained"
    assert out["A8"] == "reported"
    assert all(out[f"A{i}"] == "not simulated" for i in range(13, 18))
    assert "not reproduced, open" not in out.values()
    text, metrics = fa._summary(res)
    assert "| A9 |" in text and metrics["items"]["A2"] == "reproduced"
