"""WP5 Phase 3: rays in a spherically stratified plasma (``tss.l1_propagation.sphere``).

The flattening per ray plane is held to a calculation that does not use it:
Bouguer's law integrated directly in radius (advisor, 2026-11-25: P2's orbit
result was only self-consistent).
"""

from __future__ import annotations

import numpy as np
import pytest

from tss.l1_propagation.plate_echo import ExponentialTopside
from tss.l1_propagation.sphere import (
    LocalFrame,
    SphericalMedium,
    bouguer_paths,
    radius_and_angle,
    ray_paths,
)
from tss.platform.orbits import RE


class _Unity:
    def mu_squared(self, f_hz, depth_m):
        return np.ones_like(np.asarray(depth_m, dtype=float))


FRAME = LocalFrame(10.0, 20.0, RE + 300e3)
PLAT = FRAME.from_llh(10.0, 20.0, RE + 600e3)
POINTS = np.array([[e, n, z] for e in (0.0, 5e3, 20e3) for n in (0.0, 3e3, -40e3) for z in (0.0, 1e3)])
# Heights are above the frame radius; the profile's depth is counted from the platform.
MEDIUM = SphericalMedium(ExponentialTopside(), 300e3)


def _mu2_of_r(f, r):
    return MEDIUM.mu_squared(f, np.asarray(r) - FRAME.r0_m)


def test_platform_above_the_frame_origin_is_on_the_up_axis():
    assert PLAT == pytest.approx([0.0, 0.0, 300e3], abs=1e-6)


def test_flattening_in_free_space_is_the_chord():
    """With mu = 1 the flattened index is (r/R)^2, not 1: a real test of the map."""
    med = SphericalMedium(_Unity(), 300e3)
    hp, hx, psi = radius_and_angle(PLAT, POINTS, FRAME.r0_m)
    from tss.l1_propagation.sphere import _flat_rays

    ph, gr, _ = _flat_rays(hp, hx, psi, 5e6, med, 48, 8, None, r0_m=FRAME.r0_m)
    chord = np.linalg.norm(POINTS - PLAT, axis=1)
    assert np.max(np.abs(ph - chord)) < 1e-6
    assert np.max(np.abs(gr - chord)) < 1e-6


@pytest.mark.parametrize("f", [3.6e6, 5e6, 10e6])
def test_flattening_matches_bouguer_in_polar_coordinates(f):
    ph, gr, _ = ray_paths(PLAT, POINTS, f, MEDIUM, FRAME.r0_m)
    hp, hx, psi = radius_and_angle(PLAT, POINTS, FRAME.r0_m)
    r0 = FRAME.r0_m
    for i in range(len(POINTS)):
        bp, bg = bouguer_paths(r0 + hp, r0 + hx[i], psi[i], f, _mu2_of_r)
        assert ph[i] == pytest.approx(bp, abs=1e-6)
        assert gr[i] == pytest.approx(bg, abs=1e-6)


def test_group_path_is_d_fp_df():
    f, df = 5e6, 500.0
    lo = ray_paths(PLAT, POINTS, f - df, MEDIUM, FRAME.r0_m)[0]
    hi = ray_paths(PLAT, POINTS, f + df, MEDIUM, FRAME.r0_m)[0]
    gr = ray_paths(PLAT, POINTS, f, MEDIUM, FRAME.r0_m)[1]
    assert np.max(np.abs(((f + df) * hi - (f - df) * lo) / (2 * df) - gr)) < 5e-3


def test_strata_rising_with_the_point_shorten_the_path_by_mu_sat():
    """Vertically, lifting the profile by w changes int mu dr by mu_x - mu_p; the
    point rising with it takes off mu_x w more, so the frozen-in drift changes
    the phase path by -mu_sat w, the Doppler of the plate model."""
    w, f = 1.0, 5e6
    x = np.array([[0.0, 0.0, 0.0]])
    a = ray_paths(PLAT, x, f, MEDIUM, FRAME.r0_m)[0][0]
    b = ray_paths(PLAT, x + [0.0, 0.0, w], f, MEDIUM.shifted(w), FRAME.r0_m)[0][0]
    mu_sat = float(np.sqrt(MEDIUM.mu_squared(f, np.array([300e3]))[0]))
    assert b - a == pytest.approx(-mu_sat * w, rel=1e-4)


@pytest.mark.parametrize("r0", [1e9, 1e11, 1e13])
def test_precision_holds_for_a_very_large_sphere(r0):
    p = np.array([0.0, 0.0, 300e3])
    x = np.array([[1e3, 2e3, 0.0]])
    hp, hx, psi = radius_and_angle(p, x, r0)
    assert hp == pytest.approx(300e3, abs=1e-6)
    assert hx[0] == pytest.approx(5e6 / (2 * r0), abs=1e-9)  # the sphere's drop
    assert psi[0] * r0 == pytest.approx(np.hypot(1e3, 2e3), rel=1e-9)


def test_frame_round_trip_against_geocentric_vectors():
    fr = LocalFrame(-35.0, 140.0, RE + 250e3)
    q = fr.from_llh(-34.9, 140.2, RE + 610e3)
    ecef = fr.origin_ecef() + fr.axes().T @ q
    la, lo = np.radians([-34.9, 140.2])
    want = (RE + 610e3) * np.array([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)])
    assert ecef == pytest.approx(want, abs=1e-6)


def test_points_above_the_platform_are_refused():
    with pytest.raises(ValueError):
        ray_paths(PLAT, np.array([[0.0, 0.0, 400e3]]), 5e6, MEDIUM, FRAME.r0_m)
