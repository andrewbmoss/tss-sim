"""WP3 Phase 4: the residual leakage, its score, and the L2 evidence L5 scores."""

from __future__ import annotations

import copy
import math

import numpy as np
import pytest

from tss import sweep_law
from tss.config import load_yaml, resolve_path
from tss.constants import C_LIGHT
from tss.interfaces import EchoSet
from tss.l2_instrument import leakage as lk
from tss.l2_instrument import receive_chain as rc
from tss.l2_instrument.baseline import BaselineInstrument
from tss.l2_instrument.reference_echo import reference_levels
from tss.studies.data_volume import A1


def _cfg(**leak):
    """The Phase 4 default this file was written against: the masks' AM and no
    canceller (Phase 5 moved the default to INS-30 and INS-31)."""
    cfg = copy.deepcopy(load_yaml(resolve_path("instrument/baseline_6u.yaml")))
    cfg["leakage"].update(leak)
    cfg["transmit_impairments"]["am_noise_dbc_hz"] = -133.0
    cfg["canceller"] = {"enabled": False}
    return cfg


def _echo(f_hz=5.0e6, range_m=6.0e5, amp_db=-95.0):
    return EchoSet(
        freq_hz=np.array([f_hz]), mode=np.array(["O"], dtype=object),
        group_delay_s=np.array([2.0 * range_m / C_LIGHT]), amplitude_db=np.array([amp_db]),
        absorption_db=np.array([0.0]), flags=np.array([""], dtype=object), source="test",
    )


# ----------------------------------------------------------------------------
# The primed start
# ----------------------------------------------------------------------------


def test_the_settled_high_pass_passes_nothing_of_a_carrier_it_has_absorbed():
    sweep = A1.schedule()
    fs = 78.4e3
    hp = rc.LeakageHighpass()
    u0 = 3e-3 * np.exp(0.7j)
    x = np.full(4000, u0)
    st = hp.steady_state(u0, sweep, fs)
    warm = hp.apply(x, sweep, fs, state=st)
    cold = hp.apply(x, sweep, fs)
    assert np.max(np.abs(warm)) < 1e-9 * abs(u0)
    assert np.max(np.abs(cold)) > 0.5 * abs(u0)


def test_primed_leakage_is_there_from_the_first_sample():
    sweep = A1.schedule()
    fs = 78.4e3
    model = lk.LeakageModel(isolation_db=40.0, coupling_delay_s=3e-9)
    kw = dict(n_samples=2000, t0_s=0.0, transmit_power_w=1.0, tx_chain_efficiency_db=0.0)
    cold = lk.leakage_signal(sweep, fs, model, np.random.default_rng(0), **kw)
    warm = lk.leakage_signal(sweep, fs, model, np.random.default_rng(0), primed=True, **kw)
    assert cold[0] == 0 and abs(warm[0]) == pytest.approx(1e-2, rel=1e-6)
    assert np.allclose(warm[1:], cold[1:])


# ----------------------------------------------------------------------------
# Reading the residual
# ----------------------------------------------------------------------------


def _levels():
    inst = BaselineInstrument.from_config(_cfg())
    sweep = A1.schedule()
    return inst, sweep, reference_levels(inst, sweep), inst.sample_rate_for(sweep)


@pytest.mark.parametrize("seed", [0, 1])
def test_a_white_residual_reads_as_its_density_in_a_cell_and_shows_no_lines(seed):
    inst, sweep, lev, fs = _levels()
    n = int(sweep_law.step_count(sweep) * sweep.step_dwell_s * fs)
    n0 = 1e-18                                          # W/Hz
    rng = np.random.default_rng(seed)
    x = math.sqrt(n0 * fs / 2.0) * (rng.normal(size=n) + 1j * rng.normal(size=n))
    r = lk.dechirped_residual(x, sweep, fs, levels=lev)
    want = 10 * np.log10(n0 * lev["cell_bandwidth_hz"])
    assert np.nanmean(r["pedestal_mean_dbw"] - want) == pytest.approx(0.0, abs=0.1)
    # The largest of up to four part means reads a little high on average.
    assert np.all(np.abs(r["pedestal_dbw"] - want) < 2.0)
    assert 0.0 < np.nanmean(r["pedestal_dbw"] - want) < 0.6
    assert np.all(np.isnan(r["line_dbw"]))
    assert np.array_equal(r["level_dbw"], r["pedestal_dbw"])


def test_a_transmitter_line_is_found_at_its_level_above_the_pedestal():
    inst, sweep, lev, fs = _levels()
    n = int(sweep_law.step_count(sweep) * sweep.step_dwell_s * fs)
    rng = np.random.default_rng(3)
    n0 = 1e-18
    t = np.arange(n) / fs
    p_line = 1e-14                                      # 23 dB above the pedestal in a cell
    x = math.sqrt(n0 * fs / 2.0) * (rng.normal(size=n) + 1j * rng.normal(size=n))
    x = x + math.sqrt(p_line) * np.exp(2j * np.pi * 21.0e3 * t)
    r = lk.dechirped_residual(x, sweep, fs, levels=lev)
    inband = (lev["beat_min_hz"] < 20.5e3) & (lev["beat_max_hz"] > 21.5e3)
    assert inband.sum() > 20
    found = r["line_dbw"][inband]
    assert np.all(np.isfinite(found))
    # Hann scalloping is at most 1.42 dB; the noise in the bin adds a few tenths.
    want = 10 * np.log10(p_line)
    assert np.all((found <= want + 0.7) & (found >= want - 2.5))
    assert np.allclose(r["line_offset_hz"][inband], 21.0e3, atol=60.0)
    assert np.all(np.isnan(r["line_dbw"][~inband]))


def test_the_score_is_the_worst_step_against_the_weakest_echo():
    lev = {
        "f_hz": np.array([1e6, 2e6, 3e6]),
        "weakest_echo_dbw": np.array([-160.0, -155.0, -150.0]),
        "echo_rmax_dbw": np.array([-150.0, -140.0, -130.0]),
    }
    res = {"level_dbw": np.array([-175.0, -160.0, np.nan])}
    canc = {"level_dbw": np.array([-205.0, -190.0, np.nan])}
    model = lk.LeakageModel(isolation_db=40.0)
    s = lk.score_ins23_ins24(lev["f_hz"], model, res, levels=lev, canceller_residual=canc)
    assert s["isolation_db_min"] == 40.0 and s["ins23_met"]
    assert s["ins24_margin_db"] == pytest.approx(5.0)
    assert s["ins24_worst_freq_hz"] == 2e6
    assert s["ins24_margin_db_rmax"] == pytest.approx(20.0)
    assert s["ins24_steps_failing"] == 1 and s["ins24_steps_scored"] == 2
    assert s["ins24_margin_db_canceller"] == pytest.approx(35.0)
    assert s["ins24_steps_failing_canceller"] == 0


# ----------------------------------------------------------------------------
# The L2 evidence
# ----------------------------------------------------------------------------


def test_the_phase4_instrument_carries_its_allocations_and_the_leakage():
    inst = BaselineInstrument.from_config(_cfg())
    assert inst.leakage_enabled and inst.leakage_primed and inst.canceller_bound_db == 30.0
    imp = inst.impairments
    assert imp.oscillator is not None
    assert imp.phase_noise_floor_dbc_hz == -123.0 and imp.am_noise_dbc_hz == -133.0
    assert sorted(s.offset_hz for s in imp.spurs) == [16100.0, 39200.0]
    assert all(s.level_dbc == -82.0 and s.kind == "pm" for s in imp.spurs)


def test_the_am_pedestal_is_where_the_allocation_puts_it():
    """With the front end linear, the residual pedestal is the leakage power
    less the isolation plus the AM density in a cell: AM is not taken off by
    the de-chirp, so the -133 dBc/Hz allocation sits in every cell."""
    inst = BaselineInstrument.from_config(_cfg(), front_end_enabled=False)
    sweep = A1.schedule()
    r = inst.receive(_echo(), sweep, np.random.default_rng(5))
    ps = r.flags["leakage_per_step"]
    lev = reference_levels(inst, sweep)
    want = (-np.asarray(ps["isolation_db"]) - 133.0 + 10 * np.log10(lev["cell_bandwidth_hz"]))
    got = np.asarray(ps["residual_pedestal_dbw"])
    upper = lev["f_hz"] > 2e6                       # bands of several hundred bins
    assert np.mean(got[upper] - want[upper]) == pytest.approx(0.0, abs=0.5)
    # Behind the ideal canceller the residual is the bound lower, step by step.
    canc = np.asarray(ps["residual_level_dbw_canceller"])
    assert np.median(np.asarray(ps["residual_level_dbw"]) - canc) == pytest.approx(30.0, abs=0.5)


def test_the_l2_metrics_score_ins11_ins23_ins24_and_declare_ins14_ins26():
    inst = BaselineInstrument.from_config(_cfg())
    r = inst.receive(_echo(), A1.schedule(), np.random.default_rng(7))
    m = r.flags["l2_metrics"]
    assert m["transmit_power_w"] == 1.0 and m["transmit_power_w_basis"].startswith("declared")
    assert m["allan_deviation_composite"] == pytest.approx(2e-9, rel=1e-6)
    assert m["allan_deviation_composite_basis"].startswith("declared")
    assert m["isolation_db_min"] == pytest.approx(36.6, abs=0.1)
    # The AM allocation puts the pedestal above the weakest echo at the bottom
    # of the band: INS-24 fails without a canceller and is met behind one.
    assert m["leakage_residual_margin_db"] < 10.0
    assert m["leakage_residual_margin_db_canceller"] > 10.0
    assert m["leakage_residual_steps_scored"] == sweep_law.step_count(A1.schedule())
    assert m["dynamic_range_db_at_digitiser"] == 60.0
    assert m["ins11_rejection_margin_db"] > 0 and m["ins11_residual_peak_over_clip_db"] < 0
    assert r.flags["n_clipped"] == 0 and r.flags["leakage_primed"]
    from tss.l5_assessment.requirements_check import L2_EVIDENCE_METRICS

    # Every one but the canceller's, which this (Phase 4) instrument does not fit.
    assert set(L2_EVIDENCE_METRICS) - {"rf_canceller_min_depth_db"} <= set(m)


def test_without_the_high_pass_the_record_does_not_support_the_configured_range():
    inst = BaselineInstrument.from_config(_cfg(), highpass_enabled=False)
    m = inst.receive(_echo(), A1.schedule(), np.random.default_rng(7)).flags["l2_metrics"]
    assert m["ins11_residual_peak_over_clip_db"] > 0
    assert m["dynamic_range_db_at_digitiser"] < 60.0 - m["ins11_residual_peak_over_clip_db"] + 1e-9


def test_with_the_leakage_off_only_the_declared_values_are_reported():
    inst = BaselineInstrument.from_config(_cfg(enabled=False))
    m = inst.receive(_echo(), A1.schedule(), np.random.default_rng(7)).flags["l2_metrics"]
    assert set(m) >= {"transmit_power_w", "allan_deviation_composite"}
    assert "isolation_db_min" not in m and "leakage_residual_margin_db" not in m
    assert "dynamic_range_db_at_digitiser" not in m


class _Osc:
    def __init__(self, f):
        self.f = f

    def sigma_y(self, tau):
        return self.f(tau)


@pytest.mark.parametrize(
    ("f", "want"),
    [
        (lambda tau: 2e-9, 2e-9),                                  # goal
        (lambda tau: 1e-8, 1e-8),                                  # threshold
        (lambda tau: 1e-7, 1e-7),                                  # fail
        (lambda tau: 2e-9 * (tau / 0.05) ** 0.5, None),            # goal at 0.05 s only
    ],
)
def test_the_ins26_composite_scores_each_level_at_its_own_interval(f, want):
    from tss.l5_assessment.requirements_check import check_requirements

    inst = BaselineInstrument.from_config(_cfg(enabled=False))
    inst.oscillator = _Osc(f)
    m = inst._declared_metrics()
    if want is not None:
        assert m["allan_deviation_composite"] == pytest.approx(want)
    verdict = check_requirements(m)["INS-26"]["verdict"]
    expected = {2e-9: "goal", 1e-8: "threshold", 1e-7: "fail", None: "threshold"}[want]
    assert verdict == expected


def test_l5_scores_the_l2_evidence_and_says_what_kind_it_is():
    from tss.l5_assessment.requirements_check import check_requirements

    inst = BaselineInstrument.from_config(_cfg())
    m = inst.receive(_echo(), A1.schedule(), np.random.default_rng(7)).flags["l2_metrics"]
    rec = check_requirements(m)
    assert rec["INS-11"]["verdict"] == "threshold" and "verified" in rec["INS-11"]["basis"]
    assert rec["INS-14"]["verdict"] == "threshold" and rec["INS-14"]["basis"].startswith("declared")
    assert rec["INS-23"]["verdict"] == "fail"
    assert rec["INS-24"]["verdict"] == "fail"
    assert rec["INS-26"]["verdict"] == "goal" and rec["INS-26"]["basis"].startswith("declared")
    off = BaselineInstrument.from_config(_cfg(enabled=False))
    m_off = off.receive(_echo(), A1.schedule(), np.random.default_rng(7)).flags["l2_metrics"]
    rec_off = check_requirements(m_off)
    for rid in ("INS-11", "INS-23", "INS-24"):
        assert rec_off[rid]["verdict"] == "not_evaluated"
        assert "leakage is off" in rec_off[rid]["reason"]


def test_the_pipeline_hands_the_l2_evidence_only_to_a_tier_that_takes_it():
    from tss.l5_assessment.standard import StandardAssessment
    from tss.pipeline import _accepts

    assert _accepts(StandardAssessment().score, "l2_metrics")

    def old(retrieved, truth, ionogram, configuration, *, platform=None, truth_echo_snr=None):
        return None

    assert not _accepts(old, "l2_metrics")


# ----------------------------------------------------------------------------
# The echo windows
# ----------------------------------------------------------------------------


def test_the_echo_window_is_the_exact_delay_to_a_constant():
    inst = BaselineInstrument.from_config(_cfg())
    sweep = A1.schedule()
    fs = inst.sample_rate_for(sweep)
    n = 60000
    rz = inst.impairments.realize(sweep, fs, n, np.random.default_rng(2), max_delay_s=8e-3)
    for tau in (0.7e-3, 4.123456e-3):
        lo, hi = 20000, 22100
        exact = rz.tx_phase(tau)[lo:hi]
        win = rz.tx_phase_window(tau, lo, hi)
        # The two integrate the frequency noise from different starts, so they
        # differ by a constant.  Otherwise they agree below 0.8 of the Nyquist
        # frequency, which is where the windowed sinc is exact; above it the
        # white floors differ by some per cent, which no echo can see.
        lp_e, lp_w = _lowpass(exact), _lowpass(win)
        assert np.std(lp_w - lp_e) < 2e-3 * np.std(lp_e)
        assert np.std(win - exact) < 0.1 * np.std(exact)
        a_exact = rz.tx_amplitude(tau)[lo:hi] - 1.0
        a_win = rz.tx_amplitude_window(tau, lo, hi) - 1.0
        assert np.std(_lowpass(a_win) - _lowpass(a_exact)) < 2e-3 * np.std(_lowpass(a_exact))
        assert np.std(a_win - a_exact) < 0.1 * np.std(a_exact)


def _lowpass(x, keep=0.8):
    x = np.asarray(x, dtype=float) - np.mean(x)
    spec = np.fft.rfft(x * np.hanning(x.size))
    spec[int(keep * spec.size):] = 0.0
    return np.fft.irfft(spec, n=x.size)[200:-200]
