"""WP3 Phase 5: the RF canceller and the leakage budget."""

from __future__ import annotations

import copy
import math

import numpy as np
import pytest

from tss import sweep_law
from tss.config import load_yaml, resolve_path
from tss.l2_instrument import leakage as lk
from tss.l2_instrument.baseline import BaselineInstrument
from tss.studies.data_volume import A1


def _cfg(canceller=None, **imp):
    cfg = copy.deepcopy(load_yaml(resolve_path("instrument/baseline_6u.yaml")))
    cfg["transmit_impairments"].update(imp)
    cfg["canceller"] = {"enabled": False} if canceller is None else {"enabled": True, **canceller}
    return cfg


def _no_echoes():
    from tss.studies.leakage_residual import _no_echoes as ne

    return ne()


# ----------------------------------------------------------------------------
# The canceller
# ----------------------------------------------------------------------------


def test_the_depth_of_a_weight_error_is_its_vector_distance_from_one():
    c = lk.RFCanceller(amplitude_error_db=0.1, phase_error_deg=1.0)
    ea, ep = 10 ** (0.1 / 20) - 1, math.radians(1.0)
    assert c.nominal_depth_db() == pytest.approx(-10 * math.log10(ea**2 + ep**2))
    assert c.nominal_depth_db() == pytest.approx(33.6, abs=0.1)
    w = lk.RFCanceller(amplitude_error_db=0.5, phase_error_deg=3.0, error_draw="worst")
    assert w.nominal_depth_db() == pytest.approx(-20 * math.log10(abs(1 - (1 + 10 ** (0.5 / 20) - 1) * np.exp(1j * math.radians(3.0)))))
    assert lk.RFCanceller(ideal_depth_db=25.0).nominal_depth_db() == 25.0
    with pytest.raises(ValueError):
        lk.RFCanceller(retune="hourly")


def test_a_biased_canceller_leaves_the_bias_everywhere_on_a_flat_coupling():
    model = lk.LeakageModel(isolation_db=40.0, coupling_delay_s=3e-9)
    sweep = A1.schedule()
    t = np.linspace(0.0, 4.2, 5000)
    c = lk.RFCanceller(amplitude_error_db=0.2, phase_error_deg=2.0, error_draw="worst")
    g = c.residual_factor(sweep, t, model)
    assert np.allclose(np.abs(g) ** 2, 10 ** (-c.nominal_depth_db() / 10), rtol=1e-9)


def test_a_delay_mismatch_limits_the_depth_away_from_the_calibration():
    model = lk.LeakageModel(isolation_db=40.0, coupling_delay_s=3e-9)
    sweep = A1.schedule()
    c = lk.RFCanceller(retune="sweep", calibration_hz=1e6, delay_mismatch_s=1e-9, error_draw="worst")
    f = np.asarray(sweep_law.frequency_hz(sweep, np.array([0.5, 4.2])))
    t = np.array([0.5, 4.2]) + 3e-9
    g = c.residual_factor(sweep, t, model)
    want = np.abs(1 - np.exp(-2j * np.pi * (f - 1e6) * 1e-9))
    assert np.allclose(np.abs(g), want, rtol=1e-3)


def test_tracking_moves_the_weight_continuously_and_stepping_does_not():
    model = lk.mom_anchor_model("bus_crossed_ends_offset", "simultaneous", coupling_delay_s=3e-9)
    sweep = A1.schedule()
    fs = 78.4e3
    t = np.arange(int(20 * sweep.step_dwell_s * fs)) / fs
    kw = dict(amplitude_error_db=0.1, phase_error_deg=1.0, error_correlation_steps=0.0)
    gt = lk.RFCanceller(retune="track", **kw).residual_factor(sweep, t, model, np.random.default_rng(1))
    gs = lk.RFCanceller(retune="step", **kw).residual_factor(sweep, t, model, np.random.default_rng(1))
    assert np.max(np.abs(np.diff(gt))) < 1e-4
    assert np.max(np.abs(np.diff(gs))) > 3e-3
    # Same calibrations: the two agree at the step centres.
    centres = ((np.arange(1, 19) + 0.5) * sweep.step_dwell_s * fs).astype(int)
    assert np.allclose(gt[centres], gs[centres], atol=2e-4)


def test_correlated_errors_keep_their_rms_and_are_smooth():
    c = lk.RFCanceller(amplitude_error_db=1.0, error_correlation_steps=10.0)
    z = np.concatenate([c._correlated(np.random.default_rng(s).normal(size=159)) for s in range(200)])
    assert np.std(z) == pytest.approx(1.0, abs=0.05)
    one = c._correlated(np.random.default_rng(0).normal(size=159))
    assert np.corrcoef(one[:-1], one[1:])[0, 1] > 0.98


def test_the_configuration_block_reads():
    inst = BaselineInstrument.from_config(load_yaml(resolve_path("instrument/baseline_6u.yaml")))
    assert inst.canceller is not None and inst.canceller.retune == "track"   # INS-31 at threshold
    assert lk.RFCanceller.from_config({"enabled": False}) is None
    c = lk.RFCanceller.from_config({"enabled": True, "retune": "step", "added_noise_dbm_hz": -158})
    assert c.retune == "step" and c.added_noise_dbm_hz == -158.0


# ----------------------------------------------------------------------------
# In the L2 record
# ----------------------------------------------------------------------------


def test_the_canceller_takes_the_am_down_with_the_carrier():
    sweep = A1.schedule()
    off = BaselineInstrument.from_config(_cfg()).receive(_no_echoes(), sweep, np.random.default_rng(0))
    on = BaselineInstrument.from_config(_cfg({"ideal_depth_db": 30.0})).receive(
        _no_echoes(), sweep, np.random.default_rng(0))
    d = np.asarray(off.flags["leakage_per_step"]["residual_pedestal_dbw"]) - np.asarray(
        on.flags["leakage_per_step"]["residual_pedestal_dbw"])
    assert np.median(d) == pytest.approx(30.0, abs=0.3)
    m = on.flags["l2_metrics"]
    assert m["rf_canceller_nominal_depth_db"] == 30.0
    assert m["leakage_residual_margin_db"] == pytest.approx(off.flags["l2_metrics"]["leakage_residual_margin_db"] + 30.0, abs=1.0)


def test_an_active_canceller_adds_its_own_noise():
    from tss.l2_instrument.reference_echo import reference_levels

    sweep = A1.schedule()
    inst = BaselineInstrument.from_config(_cfg({"ideal_depth_db": 30.0}, am_noise_dbc_hz=None))
    noisy = copy.deepcopy(_cfg({"amplitude_error_db": 0.0, "phase_error_deg": 0.0, "error_draw": "worst",
                                "delay_mismatch_s": 0.0, "added_noise_dbm_hz": -150.0}, am_noise_dbc_hz=None))
    r = BaselineInstrument.from_config(noisy).receive(_no_echoes(), sweep, np.random.default_rng(0))
    lev = reference_levels(inst, sweep)
    want = -180.0 + 10 * np.log10(lev["cell_bandwidth_hz"])
    got = np.asarray(r.flags["leakage_per_step"]["residual_pedestal_dbw"])
    # The pedestal is the largest of up to four part means, a few tenths high.
    assert np.median(got - want) == pytest.approx(0.3, abs=0.5)


def test_the_default_carries_ins30_and_ins31_at_threshold_and_meets_ins24():
    from tss.l5_assessment.requirements_check import check_requirements

    inst = BaselineInstrument.from_config(load_yaml(resolve_path("instrument/baseline_6u.yaml")))
    assert inst.impairments.am_noise_dbc_hz == -140.0
    assert inst.canceller is not None and inst.canceller.ideal_depth_db == 26.0
    m = inst.receive(_no_echoes(), A1.schedule(), np.random.default_rng(3)).flags["l2_metrics"]
    # The leakage noise is 6.3 dB under the system noise at the band bottom.
    assert m["leakage_noise_over_system_noise_db"] == pytest.approx(-6.3, abs=1.0)
    assert m["leakage_residual_margin_db"] > 10.0
    assert m["rf_canceller_min_depth_db"] == pytest.approx(26.0, abs=1e-6)
    rec = check_requirements(m)
    assert rec["INS-24"]["verdict"] == "threshold"
    assert rec["INS-29"]["verdict"] == "threshold"
    assert rec["INS-30"]["verdict"] == "threshold" and rec["INS-30"]["basis"].startswith("declared")
    assert rec["INS-31"]["verdict"] == "threshold"
    assert rec["INS-23"]["verdict"] == "threshold"      # 36.6 dB against 35 (2026-10-29)
    assert m["isolation_db_min"] < 40.0


def test_without_a_canceller_ins31_says_why_it_is_not_evaluated():
    from tss.l5_assessment.requirements_check import check_requirements

    m = BaselineInstrument.from_config(_cfg(None)).receive(
        _no_echoes(), A1.schedule(), np.random.default_rng(3)).flags["l2_metrics"]
    rec = check_requirements(m)
    assert rec["INS-31"]["verdict"] == "not_evaluated" and "canceller" in rec["INS-31"]["reason"]


# ----------------------------------------------------------------------------
# The budget arithmetic
# ----------------------------------------------------------------------------


def test_the_am_needed_is_where_the_worst_step_meets_ins24():
    from tss.studies import leakage_budget as lb

    ch = {
        "isolation_db": np.array([40.0, 40.0]),
        "X": np.array([1e-13, 2e-13]),
        "Y": np.array([1e-18, 1e-18]),
        "weakest_echo_dbw": np.array([-160.0, -150.0]),
        "noise_cell_dbw": np.array([-170.0, -160.0]),
        "cell_bandwidth_hz": np.array([56.0, 56.0]),
    }
    need = lb.am_needed_dbc_hz(ch, iso_name="mom_bus", power_w=1.0)
    pr = lb.predict(ch, iso_name="mom_bus", am=need, power_w=1.0)
    assert float(np.min(pr["margin_db"])) == pytest.approx(10.0, abs=1e-6)
    # Terms that do not scale with AM can already use up the allowance.
    ch2 = dict(ch, Y=np.array([1e-9, 1e-18]))
    assert lb.am_needed_dbc_hz(ch2, iso_name="mom_bus", power_w=1.0) is None
    # Halving the power buys 3 dB of AM.
    assert lb.am_needed_dbc_hz(ch, iso_name="mom_bus", power_w=0.5) == pytest.approx(need + 3.0103, abs=1e-3)
