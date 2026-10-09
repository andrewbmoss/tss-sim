"""C3 verification (TV Section 8.15): forward model, estimator and the three cases.

The forward model is checked against closed forms it does not use internally,
the estimator against the forward model, and the study's headline results are
held so that a change to either cannot silently move them.
"""

from __future__ import annotations

import numpy as np
import pytest

from tss.l1_propagation import plate_echo as pe
from tss.l3_signal.novel.route_one_estimator import estimate_route_one
from tss.studies import c3_verify as cv

ENV = {e.name: e for e in cv.ENVIRONMENTS}


def test_stratified_reduces_to_free_space_without_plasma():
    t = np.arange(-1.0, 1.0, 1 / 45.0)
    plate = pe.Plate(x_m=3e3, tilt_rad=np.radians(0.3), drift_mps=25.0)
    vacuum = pe.ExponentialTopside(f_n_sat_hz=1.0, scale_height_m=1e9)
    a = pe.echo_free_space(t, 5e6, plate)
    b = pe.echo_stratified(t, 5e6, plate, vacuum)
    assert np.max(np.abs(np.angle(a.field * np.conj(b.field)))) < 1e-6
    assert np.max(np.abs(a.group_path_m - b.group_path_m)) < 1e-4


def test_ray_solution_matches_bisection():
    from scipy.optimize import brentq

    prof, f = pe.ExponentialTopside(), 5e6
    depth = prof.depth_for_mu(f, 0.5)
    nodes, weights = np.polynomial.legendre.leggauss(64)
    offsets = np.array([-30e3, 12e3])
    phase, _ = pe._solve(offsets, 0.0, np.full(2, depth), f, prof, nodes, weights, 8)
    d = 0.5 * depth * (1 + nodes)
    mu2 = prof.mu_squared(f, d)
    for off, p in zip(offsets, phase, strict=True):
        q = brentq(lambda q, o=off: np.sum(0.5 * depth * weights * q / np.sqrt(mu2 - q * q)) - o,
                   -0.4, 0.4, xtol=1e-15)
        assert p == pytest.approx(np.sum(0.5 * depth * weights * mu2 / np.sqrt(mu2 - q * q)), abs=1e-6)


def test_plate_below_reflection_is_refused():
    prof = pe.ExponentialTopside()
    deep = pe.Plate(depth_m=prof.depth_for_mu(3e6, 0.0) + 1e3)
    with pytest.raises(ValueError):
        pe.echo_stratified(np.array([0.0]), 3e6, deep, prof)


@pytest.mark.parametrize("env", ["free_5MHz_500km", "plasma_5MHz_mu0.5"])
@pytest.mark.parametrize("x_m", [0.0, 5e3])
def test_untilted_facet_recovers_drift_and_position(env, x_m):
    r = cv.run_trial("A", ENV[env], x_m=x_m, drift_mps=25.0, degeneracy=True)
    assert r.separable
    assert abs(r.drift_err_mps) < cv.DRIFT_TOL_MPS
    assert abs(r.angle_err_deg) < cv.ANGLE_TOL_DEG
    # With t0 free every drift fits equally well.
    assert r.degeneracy_spread_hz2 < 1e-6 * r.degeneracy_reference_hz2


@pytest.mark.parametrize("env", ["free_5MHz_500km", "plasma_3MHz_mu0.5"])
def test_tilted_facet_tilt_is_read_as_drift(env):
    r = cv.run_trial("B", ENV[env], tilt_deg=0.5, drift_mps=25.0)
    assert r.separable  # nothing flags it
    assert r.drift_err_mps == pytest.approx(r.predicted_bias_mps, rel=5e-3)
    assert r.verdict.startswith("fail")


def test_free_space_tilt_bias_is_vs_sin_tilt():
    assert cv.predicted_tilt_bias_mps(ENV["free_5MHz_500km"], np.radians(1.0)) == pytest.approx(
        -7600.0 * np.sin(np.radians(1.0))
    )


@pytest.mark.parametrize("env", ["free_5MHz_500km", "plasma_5MHz_mu0.5"])
def test_extended_plate_is_declined(env):
    r = cv.run_trial("C_extended", ENV[env], plate_m=80e3, tilt_deg=0.5, drift_mps=25.0,
                     extended=True)
    assert not r.separable
    assert abs(r.slope_ratio) < 0.01
    # What it does report is drift and tilt together.
    assert r.los_rate_mps == pytest.approx(r.los_rate_predicted_mps, abs=0.1)


def test_estimator_declines_a_constant_doppler_history():
    t = np.arange(-2.0, 2.0, 1 / 45.0)
    field = np.exp(2j * np.pi * 0.8 * t)
    est = estimate_route_one(t, field, f_hz=5e6, group_path_m=5e5, mu_sat=1.0)
    assert not est.separable
    assert est.los_rate_mps == pytest.approx(0.8 * 299792458.0 / (2 * 5e6), rel=1e-6)
