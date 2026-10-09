"""Arrival-angle options: the bounds, the plasma fields and the requirement mapping."""

from __future__ import annotations

import numpy as np
import pytest

from tss.studies import aoa_options as ao


@pytest.mark.parametrize("p", [(2.9, 0.3, 0.4, 0.7), (2.5, 1.9, 0.1, -1.0), (2.8, 2.0, 0.0, 0.0)])
def test_six_component_bound_is_one_over_root_n_snr(p):
    # Independent of direction and polarization, linear polarization included.
    got = ao.direction_crb_deg(*p, 100.0, 4, ao.CHANNELS["V6"])
    assert got == pytest.approx(np.degrees(1.0 / np.sqrt(4 * 100.0)), rel=1e-6)


def test_three_element_sensors_lose_direction_for_linear_polarization():
    for opt in ("V3E", "V3H"):
        assert not np.isfinite(ao.direction_crb_deg(2.8, 2.0, 0.0, 0.0, 100.0, 1, ao.CHANNELS[opt]))


def test_required_snr_inverts_the_bound():
    snr_db = ao.required_snr_per_sweep_db(0.5, 45)
    got = ao.direction_crb_deg(2.9, 0.3, 0.4, 0.7, 10 ** (snr_db / 10), 45, ao.CHANNELS["V6"])
    assert got == pytest.approx(0.5, rel=1e-6)


def test_plasma_wave_h_is_transverse_and_e_is_not():
    b = np.array([np.cos(np.radians(30)), 0.0, -np.sin(np.radians(30))])
    th, ph = 2.9, 0.4
    k = np.array([np.sin(th) * np.cos(ph), np.sin(th) * np.sin(ph), np.cos(th)])
    f = ao.plasma_fields(th, ph, 0.09, 0.2, b, "O")
    assert abs(np.dot(f[3:], k)) < 1e-12
    assert abs(np.dot(f[:3], k)) > 1e-3
    vac = ao.plasma_fields(th, ph, 1e-9, 1e-9, b, "O")
    assert abs(np.dot(vac[:3], k)) < 1e-9


def test_loops_measure_the_wave_normal_without_plasma_bias():
    b = np.array([np.cos(np.radians(30)), 0.0, -np.sin(np.radians(30))])
    r = ao.plasma_bias_deg("V3H", 0.09, 0.2, b, "O", rng=np.random.default_rng(1), n_draws=4)
    assert r["max_deg"] < 1e-4
    r6 = ao.plasma_bias_deg("V6", 0.09, 0.2, b, "O", rng=np.random.default_rng(1), n_draws=4)
    assert 0.2 < r6["rms_deg"] < 1.0


def test_requirement_mapping():
    assert ao.angle_requirement_deg(50.0, 90.0) == pytest.approx(np.degrees(38.0 / 7600.0))
    assert ao.angle_requirement_deg(50.0, 15.0) > ao.angle_requirement_deg(50.0, 90.0) * 3.8


def test_interferometer_baseline_matches_the_tv_table():
    # TV Table c3baseline: 2 deg phase, 0.38 deg at 5 MHz -> 51 m.
    assert ao.interferometer_baseline_m(5e6, 0.38, 2.0) == pytest.approx(50.5, abs=1.0)


def test_bigger_loops_are_sky_limited_lower():
    assert ao.loop_sky_limited_floor_mhz(1.0) < ao.loop_sky_limited_floor_mhz(0.5) < ao.loop_sky_limited_floor_mhz(0.3)
