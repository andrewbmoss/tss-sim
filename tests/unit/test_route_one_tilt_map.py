"""Route-one tilt map: closed forms and one forward-model spot check."""

from __future__ import annotations

import numpy as np
import pytest

from tss.studies import route_one_tilt_map as tm


def test_bias_is_vs_sin_tilt_scaled_by_index_ratio():
    assert tm.tilt_bias_mps(1.0, 0.5, 1.0) == pytest.approx(0.5 * 7600.0 * np.sin(np.radians(1.0)))
    assert tm.tilt_bias_mps(1.0, 0.9, 0.9) == pytest.approx(7600.0 * np.sin(np.radians(1.0)))


def test_allowance_inverts_the_bias():
    a = tm.tilt_allowance_deg(50.0, 10.0, 0.5, 0.95, linear=True)
    assert tm.tilt_bias_mps(a, 0.5, 0.95) + 10.0 == pytest.approx(50.0)
    assert tm.tilt_allowance_deg(50.0, 60.0, 0.5, 0.95, linear=True) == 0.0
    assert np.isnan(tm.tilt_allowance_deg(50.0, float("nan"), 0.5, 0.95, linear=True))


def test_group_path_reduces_to_depth_without_plasma():
    vac = tm.ExponentialTopside(f_n_sat_hz=1.0, scale_height_m=1e9)
    assert tm.vertical_group_path_m(vac, 5e6, 400e3) == pytest.approx(400e3, rel=1e-9)


def test_nyquist_limit_falls_with_frequency():
    prof = tm.PROFILES["ref_1p5MHz"]
    lim = []
    for f in (5e6, 20e6):
        mu_s = float(np.sqrt(prof.mu_squared(f, 0.0)))
        pg = tm.vertical_group_path_m(prof, f, prof.depth_for_mu(f, 0.5))
        lim.append(tm._nyquist_tilt_deg(f, mu_s, pg))
    assert lim[0] > lim[1] > 0.0
