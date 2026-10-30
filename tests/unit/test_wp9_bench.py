"""Tests for the WP9 bench predictions (WP3 Phase 6)."""

from __future__ import annotations

import math

import numpy as np
import pytest

from tss.studies import wp9_bench as wb


def _ch(n: int = 4) -> dict[str, np.ndarray]:
    return {
        "X": np.full(n, 1e-13),
        "Y": np.full(n, 1e-18),
        "cell_bandwidth_hz": np.full(n, 56.0),
        "noise_cell_dbw": np.full(n, -150.0),
    }


def test_coax_delay_set_matches_free_space_paths():
    # 0.66 m of RG-58 at vf 0.66 is 1 m of free space.
    from tss.constants import C_LIGHT

    for coax, free in zip(wb.DELAY_SET_M, (1.0, 3.0, 10.0, 30.0), strict=True):
        assert wb.coax_delay_s(coax) == pytest.approx(free / C_LIGHT, rel=1e-12)


def test_paths_are_isolation_plus_canceller():
    assert wb.BENCH_PATHS_DB["E1"] == 35.0 + 26.0
    assert wb.BENCH_PATHS_DB["E2"] == 35.0 + 34.0
    # E0 is the most leakage a level-7 mixer takes at WP9-2's margin.
    assert wb.P_TX_DBM - wb.BENCH_PATHS_DB["E0"] == wb.MIXER_P1DB_DBM[0] - wb.WP9_2_MARGIN_DB


def test_leakage_noise_scales_with_path_and_am():
    ch = _ch()
    ch["Y"][:] = 0.0
    a = wb.leakage_noise_dbm_hz(ch, 39.0, -140.0)
    b = wb.leakage_noise_dbm_hz(ch, 61.0, -140.0)
    c = wb.leakage_noise_dbm_hz(ch, 61.0, -150.0)
    assert np.allclose(a - b, 22.0)
    assert np.allclose(b - c, 10.0)


def test_am_allowed_inverts_the_noise():
    ch = _ch()
    ch["noise_cell_dbw"] = np.array([-160.0, -150.0, -145.0, -140.0])
    allowed = wb.am_allowed_dbc_hz(ch, 61.0)
    for i, am in enumerate(allowed):
        n = wb.leakage_noise_dbm_hz(ch, 61.0, am)[i]
        assert n == pytest.approx(wb.kt_sys_dbm_hz(ch)[i], abs=1e-6)


def test_floor_rise():
    assert wb.floor_rise_db(-165.0, -165.0) == pytest.approx(10 * math.log10(2))
    assert wb.floor_rise_db(-175.0, -165.0) == pytest.approx(0.414, abs=1e-3)


def test_am_measurement_floor():
    # -165 dBm/Hz floor, -9 dBm carrier: N0/(2 Pc) = -159 dBc/Hz.
    assert wb.am_measurement_floor_dbc_hz() == pytest.approx(-165.0 + 9.0 - 3.0103, abs=1e-3)
    assert wb.am_measurement_floor_dbc_hz(averages=1000) == pytest.approx(-159.0103 - 15.0, abs=1e-3)


def test_ins31_jump_limit_is_minus_79_dbc():
    assert wb.INS31_JUMP_DBC == pytest.approx(-78.8, abs=0.05)
    # 0.0066 degrees of phase is the same fraction of the leakage.
    assert 20 * math.log10(math.radians(0.0066)) == pytest.approx(wb.INS31_JUMP_DBC, abs=0.05)


def test_stray_limits_and_dds_echo():
    st = wb.stray_limits(56.0)
    assert st["weakest_echo_dbm"] == pytest.approx(-165.0 + 10 * math.log10(56.0) + 10.0)
    assert st["tx_stray_carrier_db"] == 79.0
    # Behind E1 the leakage is -31 dBm: echo path at most 0 + 31 + 65 - 10 = 86 dB.
    assert st["dds_echo_path_max_db_E1"] == pytest.approx(86.0)
    assert st["dds_echo_min_dbm_E1"] == pytest.approx(-86.0)
    # A DDS-channel echo cannot be a weak echo while leakage is present.
    assert st["dds_echo_min_snr_db_E1"] > 40.0
