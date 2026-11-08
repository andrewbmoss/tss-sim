"""WP4 Phase 5: the FM-CW acceptance study (tss.studies.fmcw_acceptance)."""

from __future__ import annotations

import numpy as np
import pytest

from tss.studies import fmcw_acceptance as fa


def _flat(f, t):
    return np.full_like(f, 3.0e5)


def test_a_noisy_cell_without_a_generator_is_refused():
    """Code review F8: noise from an unseeded generator would make a study irreproducible."""
    with pytest.raises(ValueError, match="rng"):
        fa.cell(5.0e6, 0.0, fa.K_P, fa.T_P, fa.FS_P, _flat, snr_db=20.0)


def test_a_noisy_cell_is_reproducible_from_its_generator():
    a = fa.cell(5.0e6, 0.0, fa.K_P, fa.T_P, fa.FS_P, _flat, snr_db=20.0, rng=np.random.default_rng(7))
    b = fa.cell(5.0e6, 0.0, fa.K_P, fa.T_P, fa.FS_P, _flat, snr_db=20.0, rng=np.random.default_rng(7))
    assert np.array_equal(a, b)


def test_a_noiseless_cell_needs_no_generator():
    x = fa.cell(5.0e6, 0.0, fa.K_P, fa.T_P, fa.FS_P, _flat)
    assert x.size == int(round(fa.T_P * fa.FS_P))


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
