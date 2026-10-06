"""The short-dipole terminal reactance, held against Balanis and a method of moments.

The model carried twice the terminal capacitance until 2026-10-05.  These tests
hold the corrected form so the factor of two cannot come back unnoticed.
"""

from __future__ import annotations

import numpy as np
import pytest

from tss.constants import EPS0
from tss.l2_instrument.antenna import short_dipole_efficiency


def test_low_frequency_limit_is_the_half_length_capacitance() -> None:
    length, radius, f = 10.0, 5.0e-3, 0.33e6
    _, terms = short_dipole_efficiency(f, length, radius)
    c_a = np.pi * EPS0 * (length / 2.0) / (np.log(length / (2.0 * radius)) - 1.0)
    assert float(terms["x_a"]) == pytest.approx(1.0 / (2.0 * np.pi * f * c_a), rel=1e-3)


@pytest.mark.parametrize(
    ("length_m", "radius_m", "f_hz", "x_mom_ohm"),
    [
        # MATLAB Antenna Toolbox, strip-equivalent dipole, run 2026-10-05.
        (10.0, 5.0e-3, 0.33e6, 19999.0),
        (5.0, 4.0e-3, 5.0e6, 2360.5),
        (5.0, 4.0e-3, 10.0e6, 1092.3),
        (5.0, 4.0e-3, 15.0e6, 625.1),
    ],
)
def test_reactance_agrees_with_method_of_moments(
    length_m: float, radius_m: float, f_hz: float, x_mom_ohm: float
) -> None:
    _, terms = short_dipole_efficiency(f_hz, length_m, radius_m)
    assert float(terms["x_a"]) == pytest.approx(x_mom_ohm, rel=0.06)


def test_band_edge_efficiencies() -> None:
    expected_db = {5.0: -46.06, 10.0: -37.57, 20.0: -29.02, 40.0: -20.44}
    for length, want in expected_db.items():
        eta, _ = short_dipole_efficiency(0.33e6, length)
        assert 10.0 * np.log10(float(eta)) == pytest.approx(want, abs=0.02)


def test_reactance_is_zero_at_and_above_resonance() -> None:
    _, terms = short_dipole_efficiency(np.array([20.0e6]), 40.0)
    assert float(terms["x_a"][0]) == 0.0
