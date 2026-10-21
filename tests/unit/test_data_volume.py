"""Unit tests for the INS-13 data-volume study (WP0 Phase 7)."""

from __future__ import annotations

import json
from dataclasses import replace

import pytest

from tss.studies import data_volume as dv

DAY = 86_400.0


def test_a1_reproduces_rev_g_ins13():
    """159 steps x 105 bins x 8 bit = 16.7 kB; 215 MB/day at 6.7 s, 339 at 4.26 s."""
    env = dv.THRESHOLD_ENVELOPE
    assert dv.A1.steps == 159
    assert dv.A1.sweep_s == pytest.approx(4.2612)
    assert env.r_max_km == pytest.approx(1050.0)
    assert dv.range_bins(env) == 105
    m = dv.l1a_matrix_bytes(dv.A1, env)
    assert m == 16_695
    assert m * DAY / 6.7 / 1e6 == pytest.approx(215.3, abs=0.1)
    assert m * DAY / dv.A1.sweep_s / 1e6 == pytest.approx(338.5, abs=0.1)


def test_trace_extraction_reproduces_rev_g():
    """About 2.5 kB of traces per ionogram, 32 to 51 MB/day (Rev G INS-13)."""
    p = dv.DEFAULT_PARAMS
    traces = dv.l1b_bytes(dv.A1) - p.header_bytes
    assert traces == 2544
    assert traces * DAY / 6.7 / 1e6 == pytest.approx(32.8, abs=0.1)
    assert traces * DAY / dv.A1.sweep_s / 1e6 == pytest.approx(51.6, abs=0.1)


def test_raw_samples_and_ins06_rate():
    env = dv.THRESHOLD_ENVELOPE
    # INS-06: f_s >= 4 K_f R_max / c = 78.4 kHz at A1's 5.6 MHz/s and 1050 km.
    fs = dv.raw_rate_bytes_per_s(dv.A1, env) / 2
    assert fs == pytest.approx(78.4e3, rel=2e-3)
    assert dv.raw_rate_bytes_per_s(dv.A1, env) * DAY / 1e9 == pytest.approx(13.5, abs=0.1)


def test_goal_altitude_ionogram():
    """Phase 4 knock-on: 180 bins over 1800 km, 28.6 kB, 580 MB/day back to back."""
    env = dv.GOAL_ENVELOPE
    assert env.r_max_km == pytest.approx(1800.0)
    m = dv.l1a_matrix_bytes(dv.A1, env)
    assert m == 159 * 180
    assert m * DAY / dv.A1.sweep_s / 1e6 == pytest.approx(580.3, abs=0.1)


def test_back_to_back_volume_does_not_depend_on_step_count():
    """Each step adds a column and a dwell, so a lower start costs volume only at a fixed cadence."""
    env = dv.GOAL_ENVELOPE
    lo = replace(dv.A1, name="lo", f_start_hz=0.164e6)
    per_s = [dv.l1a_matrix_bytes(c, env) / c.sweep_s for c in (dv.A1, lo)]
    assert lo.steps == 173
    assert per_s[1] == pytest.approx(per_s[0], rel=1e-9)


def test_vector_configuration_ins19():
    """Three channels triple the raw rate; the delivered trace product about doubles."""
    vf = dv.vector_factor()
    assert vf["raw"] == pytest.approx(3.0)
    assert vf["l1b"] == pytest.approx(2.0, abs=0.05)
    assert 2.0 < vf["threshold_set_O1"] < 2.5
    assert vf["goal_set_O1"] < 3.0


def test_a2_is_the_vector_configuration_without_the_cap():
    assert dv.A2.n_rx == 3
    assert dv.A2.steps == 85
    assert dv.A2.sweep_s == pytest.approx(2.278)


def test_product_sets_and_duty_cycles():
    env = dv.THRESHOLD_ENVELOPE
    for duty in dv.DUTY_CYCLES.values():
        tot = duty.ionogram_fraction + duty.doppler_fraction + duty.drift_fraction
        assert tot == pytest.approx(1.0)
    r = dv.daily_volume(dv.A1, env, dv.DUTY_CYCLES["O1"])
    v, s = r["bytes_per_day"], r["sets_bytes_per_day"]
    assert r["ionograms_per_day"] == pytest.approx(DAY / dv.A1.sweep_s)
    assert v["doppler_threshold"] == 0.0 and v["drift"] == 0.0
    assert s["goal"] - s["threshold"] == pytest.approx(
        v["l1a_all"] - v["l1a_subset"] + v["l2_goal"] - v["l2_threshold"])
    assert s["raw_reference"] > 30 * s["goal"]
    o2 = dv.daily_volume(dv.A1, env, dv.DUTY_CYCLES["O2"])
    # 18-03 LT is 9 h of 24; Doppler takes 20 % of it, drift 10 %.
    assert o2["bytes_per_day"]["doppler_threshold"] > 0
    assert o2["ionograms_per_day"] == pytest.approx(r["ionograms_per_day"] * (1 - 0.375 * 0.3))


def test_bad_duty_cycles_are_refused():
    env = dv.THRESHOLD_ENVELOPE
    with pytest.raises(ValueError):
        dv.daily_volume(dv.A1, env, dv.DutyCycle("fast", cadence_s=3.0))
    with pytest.raises(ValueError):
        dv.daily_volume(dv.A1, env, dv.DutyCycle("over", ionogram_fraction=0.9, doppler_fraction=0.2))


def test_doppler_and_drift_dwells():
    env = dv.THRESHOLD_ENVELOPE
    size_t, t_t = dv.doppler_dwell(env)
    size_g, t_g = dv.doppler_dwell(env, goal=True)
    assert t_t == pytest.approx(20 * 18e-3) and t_g == pytest.approx(60 * 18e-3)
    assert size_t == 105 * 20 + 8 * 8 + 64
    size, t = dv.drift_dwell()
    assert t == 2.0 and size == 90 * 4 * 5 + 64


def test_study_writes_its_outputs(tmp_path):
    res = dv.run(tmp_path, config={}, manifest=None)
    names = {p.name for p in res.outputs}
    assert {"data_volume.csv", "data_volume.json", "data_volume.md"} <= names
    j = json.loads((tmp_path / "data_volume.json").read_text())
    assert len(j["rows"]) == 12
    assert res.metrics["l1a_matrix_bytes_A1_900km"] == 16_695
