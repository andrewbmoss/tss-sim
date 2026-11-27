"""WP5 Phase 3 part 3, decision (n): the power cut of route one's Doppler history.

``main_lobe`` is the fraction of the peak power below which a sweep pair is left
out of the Doppler history, which feeds route one's free-slope test and its zero
crossing.  The experiment (``scripts/wp5_power_cut.py``) compares higher cuts
with the one in use; nothing is adopted, so the default must not move.
"""

from __future__ import annotations

import importlib.util
import math
from pathlib import Path

import numpy as np
import pytest

from tss.l3_signal.novel.route_one_estimator import doppler_history_from_field, estimate_route_one
from tss.l3_signal.novel.route_one_joint import joint_model
from tss.studies import c3_verify as cv

KW = {"f_hz": 4.988e6, "group_path_m": 529.8e3, "mu_sat": 0.9537, "v_s_mps": 7600.0}
T = (np.arange(209) - 104) / 45.0


def _noisy(seed, snr_db=10.0):
    rng = np.random.default_rng(seed)
    s = joint_model(T, t0_s=0.0, drift_mps=25.0, width_s=0.926, **KW)
    n = (rng.standard_normal(T.size) + 1j * rng.standard_normal(T.size)) / math.sqrt(2)
    return 10 ** (snr_db / 20) * s + n


def test_the_default_cut_is_a_tenth_of_the_peak():
    for seed in range(5):
        y = _noisy(seed)
        a = estimate_route_one(T, y, **KW)
        b = estimate_route_one(T, y, main_lobe=0.1, **KW)
        assert a.slope_ratio == b.slope_ratio and a.los_rate_mps == b.los_rate_mps
        assert a.separable == b.separable
        assert a.drift_mps == b.drift_mps or (math.isnan(a.drift_mps) and math.isnan(b.drift_mps))


def test_a_higher_cut_reaches_route_ones_history_and_keeps_fewer_samples():
    y = _noisy(3)
    sizes = [doppler_history_from_field(T, y, main_lobe=c).t_s.size for c in (0.1, 0.25, 0.5)]
    assert sizes[0] > sizes[1] > sizes[2] > 0
    ratios = {estimate_route_one(T, y, main_lobe=c, **KW).slope_ratio for c in (0.1, 0.25, 0.5)}
    assert len(ratios) == 3


def test_a_higher_cut_lets_a_4_km_plate_pass_as_compact(monkeypatch):
    """The guard that argues against adopting a higher cut: c3_verify case D's
    4 km plate, declined noise-free at 0.1, is accepted at 0.25."""
    seen = []

    def capture(t, field, **kw):
        seen.append((np.array(t), np.array(field), kw))
        return estimate_route_one(t, field, **kw)

    monkeypatch.setattr(cv, "estimate_route_one", capture)
    env = {e.name: e for e in cv.ENVIRONMENTS}["free_5MHz_500km"]
    cv.run_trial("D_length", env, plate_m=4.0e3, drift_mps=25.0, half_window_s=2.5 * 2.0,
                 rng=np.random.default_rng(0))
    t, field, kw = seen[-1]
    assert not estimate_route_one(t, field, **kw).separable
    assert estimate_route_one(t, field, main_lobe=0.25, **kw).separable


@pytest.fixture(scope="module")
def script():
    path = Path(__file__).resolve().parents[2] / "scripts" / "wp5_power_cut.py"
    spec = importlib.util.spec_from_file_location("wp5_power_cut", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_the_script_runs_one_echo_at_every_cut(script):
    row = script.trial((15.0, None, 7), fit_fading=True)
    assert row["snr_db"] == 15.0 and row["tau_c_s"] is None and row["seed"] == 7
    for c in script.CUTS:
        assert row[f"formed_{c:g}"] in (True, False)
        assert math.isfinite(row[f"pinned_err_{c:g}"]) == row[f"formed_{c:g}"]
    assert math.isfinite(row["fading_err"]) == any(row[f"formed_{c:g}"] for c in script.CUTS)
    assert script.CUTS[0] == 0.1


def test_the_scripts_summary_counts_availability_and_a_common_set(script):
    rows = [script.trial((15.0, None, s), fit_fading=False) for s in range(4)]
    out = script.summarise(rows)
    case = out[(15.0, None)]
    assert case["n"] == 4 and 0 <= case["common"] <= 4
    for c in script.CUTS:
        assert case["cuts"][c]["formed"] == sum(r[f"formed_{c:g}"] for r in rows)
