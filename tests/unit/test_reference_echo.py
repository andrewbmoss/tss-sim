"""The INS-10/INS-24 reference echo (WP3 Phase 0).

The reference is only useful if it is the level the chain actually measures,
so the last tests put a synthetic echo through the L2 baseline and the L3
measurement and hold the two together.
"""

from __future__ import annotations

import copy
import math

import numpy as np
import pytest

from tss import sweep_law
from tss.config import load_yaml, resolve_path
from tss.constants import C_LIGHT, K_BOLTZ
from tss.l2_instrument import reference_echo as rx
from tss.l2_instrument.antenna import ShortDipole
from tss.l2_instrument.baseline import BaselineInstrument
from tss.studies.data_volume import A1, A2


@pytest.fixture(scope="module")
def inst():
    # The chain's calibration is checked on an ideal transmitter with the
    # leakage off: the default since WP3 Phase 4 carries the AM allocation,
    # whose pedestal is INS-24's business (test_leakage_scoring.py).
    cfg = copy.deepcopy(load_yaml(resolve_path("instrument/baseline_6u.yaml")))
    cfg["leakage"]["enabled"] = False
    cfg["transmit_impairments"] = {}
    return BaselineInstrument.from_config(cfg)


@pytest.fixture(scope="module")
def a1(inst):
    return rx.reference_levels(inst, A1.schedule())


def _at(lev, f_mhz):
    i = int(np.argmin(np.abs(lev["f_hz"] / 1e6 - f_mhz)))
    return {k: float(np.asarray(v)[i]) for k, v in lev.items()}


def test_the_old_122_db_is_the_inherited_budget_at_300_km():
    legacy = rx.legacy_echo_db(5.0e6)
    assert legacy["true_range_m"] == 300e3
    assert legacy["echo_below_tx_db"] == pytest.approx(
        rx.LEGACY_LEAKAGE_OVER_ECHO_PLUS_ISOLATION_DB, abs=1.0
    )


def test_the_free_space_echo_is_the_instruments(inst, a1):
    """1 W, the 5 m Balanis dipole at both ends, 2 x 1050 km, -6 dB polarization."""
    r = _at(a1, 5.0)
    f = r["f_hz"]
    g = float(ShortDipole(length_m=5.0).efficiency_db(f / 1e6)) + 1.76
    spreading = 20 * math.log10(C_LIGHT / f / (4 * math.pi * 2 * 1.05e6))
    assert r["gain_dbi"] == pytest.approx(g, abs=1e-6)
    assert r["echo_rmax_dbw"] == pytest.approx(0.0 + 2 * g + spreading - 6.0, abs=1e-6)
    # The Balanis dipole at 5 MHz, -10.9 dB efficiency: Rev G §5.8.
    assert g - 1.76 == pytest.approx(-10.9, abs=0.1)


def test_the_noise_in_one_cell_is_k_tsys_times_the_hann_cell_bandwidth(inst, a1):
    r = _at(a1, 5.0)
    t_sys = float(inst.antenna_temperature_k(r["f_hz"] / 1e6)) + inst.receiver.noise_temperature_k
    b = 1.5 / A1.dwell_s
    assert r["cell_bandwidth_hz"] == pytest.approx(b, rel=1e-3)
    assert r["noise_cell_dbw"] == pytest.approx(10 * math.log10(K_BOLTZ * t_sys * b), abs=0.01)
    assert r["weakest_echo_dbw"] == pytest.approx(r["noise_cell_dbw"] + 10.0)


def test_the_leakage_ratios_move_one_for_one_with_isolation(inst):
    sweep = A1.schedule()
    l40 = rx.reference_levels(inst, sweep)
    l30 = rx.reference_levels(inst, sweep, isolation_db=30.0)
    for key in ("leakage_over_weakest_db", "leakage_over_rmax_echo_db", "leakage_over_noise_db"):
        assert np.allclose(l30[key] - l40[key], 10.0)
    assert np.allclose(l40["leakage_dbw"], 0.0 - 40.0)


def test_the_echo_band_is_100_to_1050_km_at_each_steps_rate(a1):
    r = _at(a1, 5.0)
    assert r["step_rate_hz_per_s"] == pytest.approx(5.597e6, rel=1e-3)
    assert r["beat_max_hz"] == pytest.approx(2 * r["step_rate_hz_per_s"] * 1.05e6 / C_LIGHT)
    assert r["beat_min_hz"] == pytest.approx(2 * r["step_rate_hz_per_s"] * 1.0e5 / C_LIGHT)
    low = _at(a1, 0.33)
    assert low["step_rate_hz_per_s"] < 0.15 * r["step_rate_hz_per_s"]


def test_a2_has_the_same_cell_levels_and_a_wider_echo_band(inst):
    l1 = rx.reference_levels(inst, A1.schedule(), f_hz=np.array([5.0e6]))
    l2 = rx.reference_levels(inst, A2.schedule(), f_hz=np.array([5.0e6]))
    assert l2["weakest_echo_dbw"][0] == pytest.approx(l1["weakest_echo_dbw"][0])
    assert l2["beat_max_hz"][0] > 1.5 * l1["beat_max_hz"][0]


def test_the_spur_rule_reproduces_srm_section_16():
    """-(82 + 10) - 4 sin^2(pi 39227 Hz 33.4 ns) = -50.3 dBc (SRM D188)."""
    levels = {"leakage_over_rmax_echo_db": np.array([82.0]),
              "leakage_over_weakest_db": np.array([82.0]),
              "beat_max_hz": np.array([2 * 5.6e6 * 1.05e6 / C_LIGHT])}
    spur = rx.required_spur_dbc(levels, coupling_path_m=10.0, reference="rmax")
    assert spur[0] == pytest.approx(-50.30, abs=0.01)


def test_a_uniform_sweep_is_refused(inst):
    from tss.interfaces import SweepSchedule

    with pytest.raises(ValueError, match="stepped"):
        rx.reference_levels(inst, SweepSchedule(f_start_hz=1e6, f_stop_hz=2e6, k_f_hz_per_s=1e6,
                                                sweep_period_s=1.0))


@pytest.mark.parametrize("f_mhz", [1.0, 5.0])
def test_the_chain_measures_the_predicted_snr(inst, f_mhz):
    """L2 baseline plus L3 snr_at agree with the prediction to within 1.5 dB."""
    sweep = A1.schedule()
    centres = sweep_law.step_centres_hz(sweep)
    fc = float(centres[int(np.argmin(np.abs(centres - f_mhz * 1e6)))])
    out = rx.processing_loss_db(inst, sweep, fc, 1.05e6, snr_db=20.0, seeds=(11, 12, 13, 14))
    assert abs(out["loss_db"]) < 1.5
