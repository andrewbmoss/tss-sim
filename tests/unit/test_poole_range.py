"""Poole's (1985) three-cell FM-CW estimator (tss.l3_signal.poole_range)."""

from __future__ import annotations

import math

import numpy as np
import pytest

from tss.l3_signal import poole_range as pr

C = 299_792_458.0
K, T, FS = 100.0e3, 0.5, 2048.0


def _cell(f0, t0, rstar):
    n = int(T * FS)
    t = np.arange(n) / FS
    f = f0 + K * t
    r = rstar(f, t0 + t)
    return np.exp(2j * np.pi * (2 * f * r / C - 2 * K * (r / C) ** 2))


def _three(f0, df, rstar):
    return [_cell(f0, 0.0, rstar), _cell(f0, T, rstar), _cell(f0 + df, 2 * T, rstar)]


def test_closed_forms():
    assert math.isclose(pr.ambiguity_interval_m(5e3), C / 1e4)
    assert pr.range_frequency_hz(K, 5e3) == 20.0
    assert math.isclose(pr.max_unambiguous_velocity_mps(K, 5e3, 5e6), C * 10 / 1e7, rel_tol=1e-12)
    assert math.isclose(math.degrees(pr.phase_error_rms_rad(20.0)), 4.056, abs_tol=1e-3)


@pytest.mark.parametrize("r0, v", [(300e3, 10.0), (312.3e3, 0.0), (450e3, 250.0)])
def test_three_cell_recovers_range_and_velocity(r0, v):
    f0 = 5.9e6
    res = pr.three_cell(_three(f0, 5e3, lambda f, t: r0 + v * t), FS, k_hz_per_s=K, f0_hz=f0, delta_f_hz=5e3)
    # Poole drops the range change across the cells: phi_R carries 5 df V* T / c,
    # so R' reads 2.5 V* T long (module docstring)
    assert abs(res["group_range_m"] - r0 - 2.5 * v * T) < 50.0
    assert abs(res["phase_velocity_mps"] - v) < 0.05


def test_dispersive_range_is_group_range():
    # a linear phase range R*(f) = a + b f has group range a + 2 b f (d(f R*)/df)
    a, b, f0 = 280e3, 2e-3, 5.9e6
    res = pr.three_cell(_three(f0, 5e3, lambda f, t: a + b * f), FS, k_hz_per_s=K, f0_hz=f0, delta_f_hz=5e3)
    assert abs(res["group_range_m"] - (a + 2 * b * (f0 + 2.5e3 + K * T / 2))) < 1.0e3


def test_velocity_ambiguity_wraps_n():
    f0, r0 = 5.0e6 - K * T / 2, 300e3 - (300e3 % (C / (2 * K * T)))
    ok = pr.three_cell(_three(f0, 5e3, lambda f, t: r0 + 250.0 * t), FS, k_hz_per_s=K, f0_hz=f0, delta_f_hz=5e3)
    bad = pr.three_cell(_three(f0, 5e3, lambda f, t: r0 + 400.0 * t), FS, k_hz_per_s=K, f0_hz=f0, delta_f_hz=5e3)
    assert abs(ok["group_range_m"] - r0) < 1e3
    assert abs(bad["group_range_m"] - r0 - C / 1e4) < 2e3
