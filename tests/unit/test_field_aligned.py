"""Field-aligned sheet: forward model checks and the separation result."""

from __future__ import annotations

import numpy as np
import pytest

from tss.l1_propagation import sheet_echo as se
from tss.l1_propagation.plate_echo import ExponentialTopside
from tss.l3_signal.novel import field_aligned_estimator as fe
from tss.studies import field_aligned_verify as fv


def test_dipole_curvature_is_a_third_of_the_radius_at_the_equator():
    r = se.EARTH_RADIUS_M + 350e3
    assert se.dipole_curvature_radius_m(0.0, r) == pytest.approx(r / 3.0, rel=1e-9)
    assert se.dipole_curvature_radius_m(np.radians(30), r) > r / 3.0


def test_sheet_normal_is_perpendicular_to_the_field():
    fr = se.sheet_frame(se.FieldAlignedSheet(dip_rad=np.radians(20), rotation_rad=np.radians(7)))
    assert abs(np.dot(fr["n"], fr["b"])) < 1e-12
    assert abs(np.dot(fr["u"], fr["b"])) < 1e-12


def test_stratified_sheet_echo_reduces_to_free_space_in_vacuum():
    t = np.arange(-0.2, 0.2, 1 / 45.0)
    sheet = se.FieldAlignedSheet(dip_rad=np.radians(10), rotation_rad=np.radians(5), length_m=15e3)
    vac = ExponentialTopside(f_n_sat_hz=1.0, scale_height_m=1e9)
    a = se.echo_free_space(t, 5e6, sheet, height_m=424e3, heading_rad=np.radians(15))
    b = se.echo_stratified(t, 5e6, sheet, vac, height_m=424e3, heading_rad=np.radians(15))
    assert np.max(np.abs(np.angle(a.field * np.conj(b.field)))) < 1e-6


def test_plasma_path_prediction_matches_free_space_in_vacuum():
    t = np.arange(-0.3, 0.3, 1 / 15.0)
    sheet = se.FieldAlignedSheet(dip_rad=np.radians(10), rotation_rad=np.radians(5), length_m=30e3)
    p0 = se.platform_start(sheet, 424e3)
    pos = p0 + 7600.0 * t[:, None] * se.heading_vector(np.radians(15))
    vac = ExponentialTopside(f_n_sat_hz=1.0, scale_height_m=1e9)
    a = fe.predict_path(t, pos, sheet, f_hz=5e6)
    b = fe.predict_path(t, pos, sheet, f_hz=5e6, profile=vac, height_m=424e3)
    assert np.max(np.abs(a.path_m - b.path_m)) < 1e-4


def test_meridional_track_needs_no_rotation():
    r = fv.run_sheet_trial(10.0, 0.0, 5.0, half_window_s=1.0)
    assert abs(r.err_given_mps) < 0.5
    # Assumed rotation costs only the second-order term v_s sin I (1 - cos phi).
    assert r.err_assume_mps == pytest.approx(7600 * np.sin(np.radians(10)) * (1 - np.cos(np.radians(5))), abs=0.5)


def test_oblique_track_reads_rotation_as_drift_unless_given():
    r = fv.run_sheet_trial(10.0, 15.0, 5.0, half_window_s=1.0)
    assert abs(r.err_given_mps) < 0.5
    assert r.err_assume_mps == pytest.approx(r.predicted_assume_err_mps, abs=1.0)
    assert abs(r.err_assume_mps) > 100.0
