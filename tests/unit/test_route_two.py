"""Route two: the second-look estimator, the 3-D forward model and the formation geometry."""

from __future__ import annotations

import math

import numpy as np
import pytest

from tss.l1_propagation.plate_echo import Plate, echo_free_space
from tss.l1_propagation.surface_ranging import DriftingSurface, sample_surface, specular_point
from tss.l3_signal.novel.second_look_estimator import (
    PassMeasurement,
    estimate_second_look,
    group_path_from_two_frequencies,
)
from tss.platform.formation import two_node_geometry
from tss.platform.orbits import EARTH_ROTATION_RAD_PER_S, RE
from tss.studies import route_two as rt
from tss.studies.c3_verify import ENVIRONMENTS

ENV = {e.name: e for e in ENVIRONMENTS}


def _passes(surface, v1, v2, dt, off2=(0.0, 0.0)):
    a = sample_surface(surface, np.zeros(3), v1, 0.0)
    b = sample_surface(surface, np.array([off2[0], off2[1], 0.0]), v2, dt)
    return (PassMeasurement(0.0, a.range_m, a.los_rate_mps, tuple(v1)),
            PassMeasurement(dt, b.range_m, b.los_rate_mps, tuple(v2)))


def test_specular_point_of_a_plane_is_the_perpendicular_foot():
    s = DriftingSurface(z0_m=-500e3, a=0.02, b=-0.01)
    foot, r = specular_point(s, np.array([1e3, 2e3, 0.0]), 0.0)
    n = s.unit_normal(0.0, 0.0, 0.0)
    expect = abs(0.0 - s.height(1e3, 2e3, 0.0)) * n[2]
    assert r == pytest.approx(expect, abs=1e-6)


def test_same_track_recovers_drift_and_along_track_tilt():
    s = DriftingSurface(z0_m=-500e3, a=math.tan(math.radians(0.7)), b=0.0, w_mps=25.0)
    v = np.array([7200.0, 1500.0])
    m1, m2 = _passes(s, v, v, 264.0)
    est = estimate_second_look(m1, m2)
    n = s.unit_normal(0, 0, 0)
    assert not est.crossing
    assert est.drift_mps == pytest.approx(25.0, abs=1e-3)
    assert est.slope_along == pytest.approx(float(np.dot(n[:2], v / np.linalg.norm(v))), abs=1e-7)


def test_crossing_tracks_recover_the_full_tilt():
    s = DriftingSurface(z0_m=-500e3, a=0.012, b=-0.021, w_mps=-40.0)
    m1, m2 = _passes(s, np.array([7000.0, 2500.0]), np.array([7000.0, -2500.0]), 264.0)
    est = estimate_second_look(m1, m2)
    n = s.unit_normal(0, 0, 0)
    assert est.crossing
    assert est.drift_mps == pytest.approx(-40.0, abs=1e-6)
    assert (est.slope_x, est.slope_y) == pytest.approx((n[0], n[1]), abs=1e-9)


def test_antiparallel_passes_are_a_same_track_pair():
    s = DriftingSurface(z0_m=-500e3, a=0.01, w_mps=10.0)
    m1, m2 = _passes(s, np.array([7500.0, 0.0]), np.array([-7500.0, 0.0]), 120.0)
    est = estimate_second_look(m1, m2)
    assert not est.crossing
    assert est.slope_along == pytest.approx(s.unit_normal(0, 0, 0)[0], abs=1e-7)


def test_two_frequency_range_is_the_distance_in_free_space():
    plate = Plate(depth_m=500e3, length_m=1.0)
    e1 = echo_free_space(np.array([0.0]), 5e6, plate)
    e2 = echo_free_space(np.array([0.0]), 5e6 + rt.DF_HZ, plate)
    assert group_path_from_two_frequencies(e1.field[0], e2.field[0], rt.DF_HZ) == pytest.approx(500e3, abs=0.5)


@pytest.mark.parametrize("name", ["free_5MHz_500km", "plasma_3MHz_mu0.5"])
def test_physical_optics_extended_plate_passes(name):
    r = rt.po_two_pass(ENV[name], extended=True, tilt_deg=1.0, drift_mps=-50.0, dt_s=264.0)
    assert r.verdict == "pass"
    assert abs(r.drift_err_mps) < 0.05


def test_physical_optics_compact_facet_gives_horizontal_drift():
    r = rt.po_two_pass(ENV["plasma_5MHz_mu0.5"], extended=False, tilt_deg=1.0, drift_mps=25.0,
                       hdrift_mps=100.0)
    assert r.verdict == "pass"
    assert abs(r.hdrift_err_mps) < 0.1


def test_same_ground_track_repeats_and_same_plane_drifts_west():
    dt = 264.0
    g = two_node_geometry(30.0, 10.0, dt, "same_ground_track")
    assert abs(g.cross_track_offset_m) < 1.0
    sp = two_node_geometry(30.0, 10.0, dt, "same_plane")
    # Westward shift omega_E r cos(lat) dt, of which the cross-track part is cos(heading).
    shift = EARTH_ROTATION_RAD_PER_S * (RE + 500e3) * math.cos(math.radians(10.0)) * dt
    expect = shift * math.cos(math.radians(sp.node1.heading_deg))
    assert abs(sp.cross_track_offset_m) == pytest.approx(expect, rel=0.03)


def test_same_plane_cross_tilt_bias_is_offset_over_interval_times_slope():
    g = two_node_geometry(45.0, 0.0, 264.0, "same_plane")
    u = g.node1.velocity_en_mps / np.linalg.norm(g.node1.velocity_en_mps)
    c = np.array([u[1], -u[0]])
    slope = math.tan(math.radians(2.0))
    s = DriftingSurface(z0_m=-500e3, a=float(slope * c[0]), b=float(slope * c[1]), w_mps=25.0)
    m1, m2, _ = rt.surface_passes(s, g)
    est = estimate_second_look(m1, m2)
    bias = slope * g.cross_track_offset_m / (m2.t_s - m1.t_s)
    assert est.drift_mps - 25.0 == pytest.approx(bias, rel=0.02)


def test_evolution_is_flagged_and_per_pass_tilt_beats_the_mean():
    k = 2.0 * math.pi / 300e3
    s = DriftingSurface(z0_m=-500e3, amp_m=math.tan(math.radians(2.0)) / k, omega_rad_s=2 * math.pi / 900.0,
                        kx_rad_m=k, phase_rad=0.3)
    g = two_node_geometry(30.0, 10.0, 360.0, "same_ground_track")
    m1, m2, feet = rt.surface_passes(s, g)
    est = estimate_second_look(m1, m2, evolution_tol_mps=42.0)
    u = g.node1.velocity_en_mps / np.linalg.norm(g.node1.velocity_en_mps)
    dt = m2.t_s - m1.t_s
    s1 = float(np.dot(s.unit_normal(feet[0, 0], feet[0, 1], 0.0)[:2], u))
    s2 = float(np.dot(s.unit_normal(feet[1, 0], feet[1, 1], dt)[:2], u))
    assert est.flagged
    assert abs(est.slope_along_1 - s1) < 0.3 * abs(est.slope_along - s1)
    assert abs(est.slope_along_2 - s2) < 0.3 * abs(est.slope_along - s2)


def test_steady_surface_is_not_flagged():
    s = DriftingSurface(z0_m=-500e3, a=0.01, w_mps=30.0)
    g = two_node_geometry(30.0, 10.0, 264.0, "same_ground_track")
    m1, m2, _ = rt.surface_passes(s, g)
    assert not estimate_second_look(m1, m2, evolution_tol_mps=42.0).flagged


def test_compact_facet_is_outside_a_same_plane_followers_beam():
    half = rt.facet_beam_half_width_m(5e6, 500e3, 2e3)
    assert 3e3 < half < 4e3
    g = two_node_geometry(30.0, 0.0, 264.0, "same_plane")
    assert abs(g.cross_track_offset_m) > 10 * half
