"""The Doppler relations carry the refractive index at the platform.

A topside sounder is inside the plasma, so the wave it launches has the local
wavelength and not the free-space one.  The shift from a closing velocity is
``2 f mu_sat v / c``.  These tests hold the forward model to that: the factor
is what differentiating the phase path gives, it reaches every relation that
turns a velocity into a frequency or a range bias, the free-space relations are
what a caller gets by saying nothing, and the analytic tier applies it.
"""

from __future__ import annotations

import numpy as np
import pytest

from tss.constants import C_LIGHT, density_from_plasma_frequency
from tss.interfaces import IonosphereState, PlatformState, SweepSchedule
from tss.l1_propagation.analytic import AnalyticPropagation
from tss.l1_propagation.doppler import (
    apparent_range_bias_m,
    coherent_interval_for_velocity_resolution_s,
    doppler_geometry,
    doppler_shift_hz,
    refractive_index_at_platform,
    residual_bias_second_order_m,
    velocity_from_doppler_mps,
)

ALT_SAT_M = 8.0e5
F_N_SAT_HZ = 2.0e6
SCALE_HEIGHT_M = 1.5e5


def _f_n(h_m):
    """Plasma frequency of an exponential topside, 2 MHz at the spacecraft."""
    return F_N_SAT_HZ * np.exp((ALT_SAT_M - np.asarray(h_m, dtype=float)) / (2.0 * SCALE_HEIGHT_M))


def _phase_path_m(f_hz, h_sat_m, n=200_000):
    """Integral of the refractive index from reflection to ``h_sat_m``.

    Budden (1961) eq. 10.14, evaluated with ``h = h_r + (h_sat - h_r) s^2`` so
    that the square-root zero at the reflection level is integrated cleanly.
    """
    h_refl = ALT_SAT_M - 2.0 * SCALE_HEIGHT_M * np.log(f_hz / F_N_SAT_HZ)
    s = (np.arange(n) + 0.5) / n
    h = h_refl + (h_sat_m - h_refl) * s**2
    mu = np.sqrt(np.clip(1.0 - (_f_n(h) / f_hz) ** 2, 0.0, None))
    return float(np.sum(mu * 2.0 * (h_sat_m - h_refl) * s) / n)


def test_index_at_platform_values_and_band_edge():
    mu = refractive_index_at_platform([4.0e6, 3.0e6, 2.2e6], F_N_SAT_HZ)
    assert mu == pytest.approx([np.sqrt(0.75), np.sqrt(5.0) / 3.0, np.sqrt(1.0 - 1.0 / 1.21)])
    # No plasma at the platform is free space.
    assert float(refractive_index_at_platform(5.0e6, 0.0)) == 1.0
    # At and below the plasma frequency the wave is evanescent: no index.
    assert np.all(np.isnan(refractive_index_at_platform([2.0e6, 1.0e6], F_N_SAT_HZ)))


@pytest.mark.parametrize("f_hz", [2.2e6, 3.0e6, 6.0e6])
def test_factor_is_the_derivative_of_the_phase_path(f_hz):
    """Moving the platform end of the phase path by dz changes it by mu_sat dz."""
    dz = 50.0
    slope = (_phase_path_m(f_hz, ALT_SAT_M + dz) - _phase_path_m(f_hz, ALT_SAT_M - dz)) / (2 * dz)
    assert slope == pytest.approx(float(refractive_index_at_platform(f_hz, F_N_SAT_HZ)), rel=1e-4)


def test_default_is_the_free_space_relation():
    assert float(doppler_shift_hz(100.0, 5.0e6)) == pytest.approx(2 * 100.0 * 5.0e6 / C_LIGHT)
    assert float(apparent_range_bias_m(100.0, 5.0e6, 5.6e6)) == pytest.approx(-100.0 * 5.0e6 / 5.6e6)


def test_every_relation_carries_the_factor():
    f, v, k, mu = 3.0e6, 250.0, 5.6e6, float(refractive_index_at_platform(3.0e6, F_N_SAT_HZ))
    assert mu == pytest.approx(0.745356, rel=1e-6)
    assert float(doppler_shift_hz(v, f, mu_sat=mu)) == pytest.approx(mu * float(doppler_shift_hz(v, f)))
    for direction in (1, -1):
        assert float(apparent_range_bias_m(v, f, k, direction=direction, mu_sat=mu)) == pytest.approx(
            mu * float(apparent_range_bias_m(v, f, k, direction=direction))
        )
    assert float(residual_bias_second_order_m(0.5, f, k, 4.26, mu_sat=mu)) == pytest.approx(
        mu * float(residual_bias_second_order_m(0.5, f, k, 4.26))
    )
    assert float(coherent_interval_for_velocity_resolution_s(50.0, f, mu_sat=mu)) == pytest.approx(
        float(coherent_interval_for_velocity_resolution_s(50.0, f)) / mu
    )


def test_velocity_round_trip_and_the_error_of_omitting_the_factor():
    f, v, mu = 3.0e6, 40.0, float(refractive_index_at_platform(3.0e6, F_N_SAT_HZ))
    f_d = doppler_shift_hz(v, f, mu_sat=mu)
    assert float(velocity_from_doppler_mps(f_d, f, mu_sat=mu)) == pytest.approx(v)
    # An estimator that assumes free space reads the velocity low by mu_sat.
    assert float(velocity_from_doppler_mps(f_d, f)) == pytest.approx(mu * v)


@pytest.mark.parametrize("bad", [0.0, -0.2, 1.01])
def test_index_outside_the_physical_range_is_refused(bad):
    with pytest.raises(ValueError):
        doppler_shift_hz(10.0, 5.0e6, mu_sat=bad)


def _state_and_platform():
    alt = np.linspace(2.0e5, 1.0e6, 1601)
    state = IonosphereState(alt_m=alt, ne=density_from_plasma_frequency(_f_n(alt)), nu=np.full(alt.shape, 100.0))
    platform = PlatformState(
        alt_m=ALT_SAT_M,
        lat_deg=0.0,
        lon_deg=0.0,
        velocity_ecef_mps=np.array([7600.0, 0.0, 30.0]),
        epoch="2026-01-01T00:00:00Z",
    )
    return state, platform


def test_analytic_tier_applies_the_factor():
    state, platform = _state_and_platform()
    freqs = np.array([1.5e6, 2.5e6, 4.0e6, 8.0e6])
    echoes = AnalyticPropagation(modes=("O",)).echoes(
        state, platform, freqs, np.random.default_rng(0)
    )
    f_n_sat = float(AnalyticPropagation().plasma_frequency_at(state, ALT_SAT_M))
    assert f_n_sat == pytest.approx(F_N_SAT_HZ, rel=1e-3)
    # Below the plasma frequency at the spacecraft there is no echo and no shift.
    assert np.isnan(echoes.doppler_hz[0])
    mu = np.sqrt(1.0 - (f_n_sat / freqs[1:]) ** 2)
    assert echoes.doppler_hz[1:] == pytest.approx(2.0 * mu * 30.0 * freqs[1:] / C_LIGHT, rel=1e-9)
    # The factor tends to 1 up the band but is 0.6 at 1.25 times f_N.
    assert mu[0] == pytest.approx(0.6, rel=2e-3)


def test_geometry_records_the_index_it_used():
    _, platform = _state_and_platform()
    sweep = SweepSchedule(
        f_start_hz=3.3e5, f_stop_hz=2.0e7, k_f_hz_per_s=5.6e6, sweep_period_s=4.0,
        direction=1,
    )
    plain = doppler_geometry(platform, [3.0e6], sweep)[0]
    inside = doppler_geometry(platform, [3.0e6, 1.0e6], sweep, f_n_sat_hz=F_N_SAT_HZ)
    assert plain.mu_sat == 1.0
    assert inside[0].mu_sat == pytest.approx(0.745356, rel=1e-6)
    assert inside[0].doppler_hz == pytest.approx(0.745356 * plain.doppler_hz, rel=1e-6)
    assert inside[0].range_bias_m == pytest.approx(0.745356 * plain.range_bias_m, rel=1e-6)
    assert np.isnan(inside[1].doppler_hz) and np.isnan(inside[1].range_bias_m)
