"""WP5 Phase 2: the general track entry of plate_echo reduces to the verified ones.

D4: ``echo_free_space`` and ``echo_stratified`` are left unedited.
``echo_track`` takes the platform's along-track position and its extra depth
per time; on a straight level track it must give the same numbers exactly.
"""

from __future__ import annotations

import numpy as np
import pytest

from tss.l1_propagation import plate_echo as pe

T = np.linspace(-1.5, 1.5, 41)


@pytest.mark.parametrize("tilt", [0.0, np.radians(0.5)])
def test_free_space_reduces_exactly(tilt):
    plate = pe.Plate(x_m=1.0e3, depth_m=500e3, tilt_rad=tilt, drift_mps=25.0, drift_x_mps=3.0)
    a = pe.echo_free_space(T, 5.0e6, plate)
    b = pe.echo_track(T, 5.0e6, plate, None)
    assert np.array_equal(a.field, b.field)
    assert np.array_equal(a.group_path_m, b.group_path_m)
    assert np.array_equal(b.phase_path_m, b.group_path_m)
    assert b.mu_sat == 1.0


def test_stratified_reduces_exactly():
    prof = pe.ExponentialTopside()
    plate = pe.Plate(depth_m=prof.depth_for_mu(5.0e6, 0.5), drift_mps=25.0)
    a = pe.echo_stratified(T, 5.0e6, plate, prof)
    b = pe.echo_track(T, 5.0e6, plate, prof)
    assert np.array_equal(a.field, b.field)
    assert np.array_equal(a.group_path_m, b.group_path_m)
    assert a.mu_sat == b.mu_sat
    # In a plasma the phase path is shorter than the group path.
    assert np.all(b.phase_path_m < b.group_path_m)


def test_explicit_track_equals_the_straight_one():
    prof = pe.ExponentialTopside()
    plate = pe.Plate(depth_m=prof.depth_for_mu(5.0e6, 0.5), drift_mps=25.0)
    a = pe.echo_track(T, 5.0e6, plate, prof)
    b = pe.echo_track(T, 5.0e6, plate, prof, along_m=7600.0 * T, extra_depth_m=np.zeros_like(T))
    assert np.array_equal(a.field, b.field)


def test_descending_platform_lengthens_the_path():
    prof = pe.ExponentialTopside()
    plate = pe.Plate(depth_m=prof.depth_for_mu(5.0e6, 0.5))
    a = pe.echo_track(np.array([0.0]), 5.0e6, plate, prof)
    b = pe.echo_track(np.array([0.0]), 5.0e6, plate, prof, extra_depth_m=np.array([-100.0]))
    # A platform 100 m higher (extra depth -100 m) sees a longer path.
    assert b.group_path_m[0] > a.group_path_m[0]
    assert b.phase_path_m[0] - a.phase_path_m[0] == pytest.approx(100.0 * a.mu_sat, rel=1e-3)
