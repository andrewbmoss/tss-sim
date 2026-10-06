"""Route one of C3: the model, the overflight-time Monte Carlo and the budget.

The relations are the plasma forms of Technical Volume Section 8.2.  These
tests hold the study to them and to the reference case the chapter quotes, a
plasma frequency at the spacecraft of 1.5 MHz, so that a number in the chapter
can be traced to a function here.
"""

from __future__ import annotations

import numpy as np
import pytest

from tss.constants import C_LIGHT
from tss.studies import route_one as r1

F5, P500 = 5.0e6, 5.0e5


def test_reference_case_puts_3_mhz_at_the_requirement_floor():
    assert float(r1.mu_sat(3.0e6)) == pytest.approx(np.sqrt(0.75))
    assert float(r1.mu_sat(5.0e6)) == pytest.approx(0.9539, abs=1e-4)
    assert np.isnan(float(r1.mu_sat(1.0e6)))


def test_slope_uses_the_group_path_and_no_index():
    assert float(r1.doppler_slope_hz_per_s(F5, P500)) == pytest.approx(3.853, abs=2e-3)
    assert float(r1.doppler_slope_hz_per_s(2.0e7, P500)) == pytest.approx(15.41, abs=1e-2)


def test_drift_and_overflight_time_give_the_same_history():
    mu = float(r1.mu_sat(F5))
    shift = float(r1.zero_crossing_shift_s(25.0, P500, mu))
    assert shift == pytest.approx(0.2065, abs=5e-4)
    t = np.linspace(-0.5, 0.5, 51)
    drifting = r1.doppler_history_hz(t, F5, P500, v_drift_mps=25.0, t0_s=0.0, mu=mu)
    displaced = r1.doppler_history_hz(t, F5, P500, v_drift_mps=0.0, t0_s=shift, mu=mu)
    assert drifting == pytest.approx(displaced, abs=1e-12)
    # The zero crossing is where the drifting history vanishes.
    assert float(r1.doppler_history_hz(np.array([shift]), F5, P500, v_drift_mps=25.0, mu=mu)[0]) == (
        pytest.approx(0.0, abs=1e-12)
    )
    assert float(r1.drift_from_zero_crossing_mps(shift, 0.0, P500, mu)) == pytest.approx(25.0)


def test_free_space_estimator_overstates_the_drift_by_the_index():
    mu = float(r1.mu_sat(3.0e6))
    shift = float(r1.zero_crossing_shift_s(40.0, P500, mu))
    assert float(r1.drift_from_zero_crossing_mps(shift, 0.0, P500, 1.0)) == pytest.approx(mu * 40.0)


def test_envelope_width_is_the_half_power_width():
    width = float(r1.envelope_width_s(F5, P500, 2.0e3))
    assert width == pytest.approx(0.443 * (C_LIGHT / F5) * P500 / (2.0e3 * 7600.0), rel=1e-3)
    assert width == pytest.approx(0.874, abs=2e-3)
    assert r1.envelope_power(np.array([-width / 2, 0.0, width / 2]), 0.0, width) == pytest.approx(
        [0.5, 1.0, 0.5], abs=1e-9
    )


def test_overflight_time_tolerance():
    assert float(r1.overflight_time_tolerance_s(50.0, P500, float(r1.mu_sat(F5)))) == pytest.approx(
        0.413, abs=1e-3
    )
    assert float(r1.overflight_time_tolerance_s(50.0, P500, 1.0)) == pytest.approx(0.433, abs=1e-3)


def test_fit_recovers_a_clean_envelope():
    t = np.arange(-98, 99) / 45.0
    power = 10.0 * r1.envelope_power(t, 0.137, 0.87) + 1.0
    centre, width = r1.fit_overflight_time(t, power)
    assert centre == pytest.approx(0.137, abs=1e-4)
    assert width == pytest.approx(0.87, abs=1e-3)


def test_monte_carlo_is_repeatable_and_fading_costs():
    width = float(r1.envelope_width_s(F5, P500, 2.0e3))
    clear = r1.overflight_time_sigma_s(width, n_trials=120, seed=7)
    assert clear == r1.overflight_time_sigma_s(width, n_trials=120, seed=7)
    assert 0.015 < clear < 0.06
    faded = r1.overflight_time_sigma_s(width, n_trials=120, seed=7, fading_tau_s=0.5)
    assert faded > 2.0 * clear


def test_budget_terms_at_the_floor():
    mu = float(r1.mu_sat(3.0e6))
    # INS-26 is one tenth of 100 m/s at this index, so the oscillator is 10 m/s.
    assert r1.oscillator_term_mps(mu) == pytest.approx(10.0, abs=1e-9)
    assert r1.oscillator_term_mps(1.0) == pytest.approx(8.66, abs=5e-3)
    assert r1.doppler_precision_mps(3.0e6, 1.0) == pytest.approx(2.60, abs=1e-2)
    assert r1.doppler_precision_mps(3.0e6, mu) == pytest.approx(2.60 / mu, abs=1e-2)
    assert r1.ephemeris_term_mps(P500, 1.0) == pytest.approx(0.18, abs=1e-2)


def test_budget_row_is_the_root_sum_square_of_its_terms():
    row = r1.route_one_budget(3.0e6, 5.0e5, n_trials=60, seed=3)
    parts = (row.overflight_mps, row.doppler_mps, row.oscillator_mps, row.ephemeris_mps,
             row.range_scale_mps)
    assert row.total_mps == pytest.approx(float(np.sqrt(np.sum(np.square(parts)))))
    assert row.overflight_mps == pytest.approx(7600.0**2 * row.sigma_t0_s / (row.mu_sat * 5.0e5))
    assert row.range_scale_mps == pytest.approx(2.5)


def test_range_bias_and_tilt_tables():
    bias = r1.range_bias_table_m(np.array([25.0, 133.0]), np.array([1.0e6, 2.0e6, 5.0e6, 2.0e7]))
    assert np.isnan(bias[0, 0])
    assert bias[0, 1] == pytest.approx(25.0 * 0.536 * float(r1.mu_sat(2.0e6)), rel=1e-3)
    assert bias[0, 2] == pytest.approx(21.3, abs=0.05)
    assert bias[1, 3] == pytest.approx(473.7, abs=0.5)
    v_r, f_d = r1.tilt_doppler_table(np.array([1.0]), np.array([5.0e6, 2.0e7]))
    assert v_r[0] == pytest.approx(132.6, abs=0.05)
    assert f_d[0] == pytest.approx([4.22, 17.65], abs=0.01)


def test_study_is_registered_and_writes_its_tables(tmp_path):
    from tss.studies import available_studies

    assert available_studies()["route_one"]["status"] == "implemented"

    class _Manifest:
        outputs: list = []

    result = r1.run(tmp_path, config={"route_one": {"n_trials": 12, "seed": 1}}, manifest=_Manifest())
    assert {p.name for p in result.outputs} == {"route_one_overflight_mc.csv", "route_one_budget.csv"}
    budget = np.loadtxt(tmp_path / "route_one_budget.csv", delimiter=",", skiprows=1)
    assert budget.shape == (12, 12)
    assert result.metrics["mu_sat_3mhz"] == pytest.approx(np.sqrt(0.75))
