"""WP5 Phase 1: the echo fading model (tss.l2_instrument.fading).

Decisions of 2026-11-24 (gate G1), recorded in docs/decorrelation.md:
the intrinsic decorrelation time comes from ground-based literature, held at
low confidence and high variance; the diffuse part's decorrelation is tied to
the facet through the halo extent; the intrinsic times carried are 0.1, 0.5 and
1.0 s.  The model is standalone until P2 wires it into the echo chain.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from tss.l2_instrument import fading as fd

DT = 1.0 / 45.0  # drift-mode sweep period, s


def _t(n: int = 2048) -> np.ndarray:
    return np.arange(n) * DT


def _acf_estimate(gains: np.ndarray, max_lag: int) -> np.ndarray:
    """Normalised field ACF averaged over echoes (rows) at lags 0..max_lag."""
    g = gains - 0.0
    p0 = np.mean(np.abs(g) ** 2)
    out = np.empty(max_lag + 1)
    for k in range(max_lag + 1):
        out[k] = np.real(np.mean(g[:, k:] * np.conj(g[:, : g.shape[1] - k]))) / p0
    return out


def _many(model: fd.FadingModel, n_echo: int, seed: int, n: int = 2048) -> np.ndarray:
    streams = fd.echo_streams(np.random.default_rng(seed), n_echo)
    return np.array([model.realise(_t(n), s) for s in streams])


# --------------------------------------------------------------------------
# Decision record
# --------------------------------------------------------------------------


def test_g1_decisions_are_the_module_defaults():
    assert fd.CARRY_TAU_I_S == (0.1, 0.5, 1.0)
    assert fd.CARRY_HALO_RATIOS == (1.0, 2.0, 4.0)
    assert fd.CARRY_K == (math.inf, 10.0, 3.0, 1.0, 0.0)
    # Low confidence, high variance: the intrinsic time is drawn per echo.
    assert fd.LOW_CONFIDENCE_LOG_SIGMA > 0.0
    assert fd.FadingModel(tau_i_s=0.5).tau_log_sigma == fd.LOW_CONFIDENCE_LOG_SIGMA


# --------------------------------------------------------------------------
# Correlation in time
# --------------------------------------------------------------------------


@pytest.mark.parametrize("kind", ["exp", "gauss"])
def test_time_acf_matches_the_closed_form_at_lags_1_to_50(kind):
    m = fd.FadingModel(tau_i_s=0.5, kind=kind, tau_log_sigma=0.0)
    est = _acf_estimate(_many(m, 200, 11), 50)
    lags = np.arange(51) * DT
    want = m.acf(lags)
    assert want[0] == pytest.approx(1.0)
    assert np.max(np.abs(est[1:] - want[1:])) < 0.05


def test_closed_forms():
    m = fd.FadingModel(tau_i_s=0.5, kind="exp", tau_log_sigma=0.0)
    assert m.acf(0.5) == pytest.approx(math.exp(-1.0))
    g = fd.FadingModel(tau_i_s=0.5, kind="gauss", tau_log_sigma=0.0)
    assert g.acf(0.5) == pytest.approx(math.exp(-1.0))
    # Two time scales: K / (K+1) of the power is specular and slow.
    two = fd.FadingModel(tau_i_s=1.0, k_factor=1.0, tau_d_s=0.05, kind="exp", tau_log_sigma=0.0)
    assert two.acf(0.2) == pytest.approx(0.5 * math.exp(-0.2) + 0.5 * math.exp(-4.0))


def test_two_time_scale_acf_is_realised():
    m = fd.FadingModel(tau_i_s=1.0, k_factor=1.0, tau_d_s=0.05, kind="exp", tau_log_sigma=0.0)
    est = _acf_estimate(_many(m, 200, 12), 50)
    assert np.max(np.abs(est[1:] - m.acf(np.arange(1, 51) * DT))) < 0.05


def test_diffuse_part_shorter_than_a_sweep_is_white():
    for kind in ("exp", "gauss"):
        m = fd.FadingModel(k_factor=0.0, tau_d_s=0.005, kind=kind, tau_log_sigma=0.0)
        est = _acf_estimate(_many(m, 100, 13, n=1024), 3)
        assert abs(est[1]) < 0.05


# --------------------------------------------------------------------------
# Power and amplitude statistics
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    "kw",
    [
        dict(tau_i_s=0.5, kind="exp"),
        dict(tau_i_s=0.5, kind="gauss"),
        dict(k_factor=0.0, tau_d_s=0.1),
        dict(tau_i_s=1.0, k_factor=3.0, tau_d_s=0.05),
        dict(k_factor=10.0, tau_d_s=0.02),
    ],
)
def test_unit_mean_power(kw):
    g = _many(fd.FadingModel(tau_log_sigma=0.0, **kw), 200, 14)
    assert np.mean(np.abs(g) ** 2) == pytest.approx(1.0, abs=0.05)


def test_no_fading_is_unity():
    g = fd.FadingModel(tau_i_s=None, k_factor=math.inf).realise(_t(64), np.random.default_rng(1))
    assert np.allclose(g, 1.0)


def test_k_zero_is_rayleigh():
    g = _many(fd.FadingModel(k_factor=0.0, tau_d_s=None, tau_log_sigma=0.0), 50, 15)
    p = np.abs(g) ** 2
    assert np.mean(p) == pytest.approx(1.0, abs=0.03)
    assert np.mean(p**2) == pytest.approx(2.0, abs=0.1)


@pytest.mark.parametrize("k", [1.0, 3.0, 10.0])
def test_steady_specular_with_diffuse_is_rician(k):
    g = _many(fd.FadingModel(tau_i_s=None, k_factor=k, tau_d_s=None, tau_log_sigma=0.0), 50, 16)
    p = np.abs(g) ** 2
    want = (k * k + 4.0 * k + 2.0) / (k + 1.0) ** 2
    assert np.mean(p**2) == pytest.approx(want, rel=0.05)


# --------------------------------------------------------------------------
# Random streams
# --------------------------------------------------------------------------


def test_same_stream_reproduces_and_echo_streams_are_independent():
    m = fd.FadingModel(tau_i_s=0.5, k_factor=3.0, tau_d_s=0.05)
    a = m.realise(_t(512), np.random.default_rng(99))
    b = m.realise(_t(512), np.random.default_rng(99))
    assert np.array_equal(a, b)
    s1, s2 = fd.echo_streams(np.random.default_rng(5), 2)
    x = fd.FadingModel(tau_i_s=None, k_factor=0.0, tau_d_s=None).realise(_t(4096), s1)
    y = fd.FadingModel(tau_i_s=None, k_factor=0.0, tau_d_s=None).realise(_t(4096), s2)
    assert abs(np.mean(x * np.conj(y))) < 0.05


def test_switching_the_diffuse_part_leaves_the_specular_draw_alone():
    a = fd.FadingModel(tau_i_s=0.5, k_factor=math.inf).realise_parts(_t(256), np.random.default_rng(3))
    b = fd.FadingModel(tau_i_s=0.5, k_factor=1.0, tau_d_s=0.05).realise_parts(
        _t(256), np.random.default_rng(3)
    )
    assert np.array_equal(a.specular, b.specular)
    assert a.tau_i_s == b.tau_i_s


def test_intrinsic_time_is_drawn_per_echo_with_the_stated_spread():
    m = fd.FadingModel(tau_i_s=0.5, tau_log_sigma=1.0)
    taus = np.array([m.draw_tau_i_s(s) for s in fd.echo_streams(np.random.default_rng(8), 4000)])
    assert np.mean(np.log(taus)) == pytest.approx(math.log(0.5), abs=0.05)
    assert np.std(np.log(taus)) == pytest.approx(1.0, abs=0.05)
    assert fd.FadingModel(tau_i_s=0.5, tau_log_sigma=0.0).draw_tau_i_s(np.random.default_rng(1)) == 0.5


# --------------------------------------------------------------------------
# Correlation across frequency
# --------------------------------------------------------------------------


def test_frequency_correlation_follows_the_coherence_bandwidth():
    f = 5.0e6 + np.arange(4) * 50.0e3
    m = fd.FadingModel(k_factor=0.0, tau_d_s=0.1, coherence_bw_hz=100.0e3, tau_log_sigma=0.0)
    gs = [m.realise(_t(512), s, f_hz=f) for s in fd.echo_streams(np.random.default_rng(21), 200)]
    g = np.array(gs)  # echo, frequency, time
    assert g.shape == (200, 4, 512)
    for j, df in ((1, 50.0e3), (2, 100.0e3), (3, 150.0e3)):
        c = np.real(np.mean(g[:, 0, :] * np.conj(g[:, j, :])))
        assert c == pytest.approx(math.exp(-df / 100.0e3), abs=0.05)
        assert m.acf(0.0, df) == pytest.approx(math.exp(-df / 100.0e3))


def test_no_coherence_bandwidth_means_one_draw_for_all_frequencies():
    f = np.array([3.0e6, 3.1e6, 3.2e6])
    g = fd.FadingModel(k_factor=0.0, tau_d_s=0.1).realise(_t(64), np.random.default_rng(2), f_hz=f)
    assert np.array_equal(g[0], g[1]) and np.array_equal(g[0], g[2])


def test_frequency_order_does_not_matter():
    f = np.array([5.1e6, 5.0e6, 5.2e6])
    m = fd.FadingModel(k_factor=0.0, tau_d_s=0.1, coherence_bw_hz=1.0e5)
    a = m.realise(_t(64), np.random.default_rng(4), f_hz=f)
    b = m.realise(_t(64), np.random.default_rng(4), f_hz=np.sort(f))
    assert np.array_equal(a[1], b[0]) and np.array_equal(a[0], b[1]) and np.array_equal(a[2], b[2])


# --------------------------------------------------------------------------
# The diffuse part is tied to the facet (decision 2 of G1)
# --------------------------------------------------------------------------


def test_halo_decorrelation_reference_values():
    # 5 MHz, 500 km: 0.5 s needs a patch of about 2.8 km, 0.1 s about 14 km.
    assert fd.halo_extent_m(0.5, 5.0e6, 5.0e5) == pytest.approx(2.76e3, rel=0.01)
    assert fd.halo_extent_m(0.1, 5.0e6, 5.0e5) == pytest.approx(13.8e3, rel=0.01)
    L = 4.0e3
    assert fd.halo_extent_m(fd.halo_decorrelation_s(5.0e6, 5.0e5, L), 5.0e6, 5.0e5) == pytest.approx(L)


@pytest.mark.parametrize("f_hz,p_m", [(3.0e6, 5.0e5), (5.0e6, 1.0e6), (10.0e6, 5.0e5)])
def test_diffuse_time_against_the_facet_envelope_is_0p79_d_over_l(f_hz, p_m):
    from tss.studies.route_one import envelope_width_s

    for d, L in ((2.0e3, 2.0e3), (2.0e3, 8.0e3), (1.0e3, 0.5e3)):
        ratio = fd.halo_decorrelation_s(f_hz, p_m, L) / float(envelope_width_s(f_hz, p_m, d))
        assert ratio == pytest.approx(fd.diffuse_to_envelope_ratio(d, L), rel=1e-9)
        assert ratio == pytest.approx(0.789 * d / L, rel=0.002)


def test_from_halo_ties_the_diffuse_time_and_bandwidth_to_the_geometry():
    m = fd.FadingModel.from_halo(
        f_hz=5.0e6, p_group_m=5.0e5, facet_m=2.0e3, halo_ratio=2.0, tau_i_s=0.5, k_factor=3.0
    )
    assert m.tau_d_s == pytest.approx(fd.halo_decorrelation_s(5.0e6, 5.0e5, 4.0e3))
    assert m.coherence_bw_hz == pytest.approx(fd.halo_coherence_bandwidth_hz(5.0e5, 4.0e3))
    assert m.k_factor == 3.0 and m.tau_i_s == 0.5
    # A wider halo decorrelates faster and over less bandwidth.
    assert fd.halo_coherence_bandwidth_hz(5.0e5, 14e3) < fd.halo_coherence_bandwidth_hz(5.0e5, 2.8e3)


def test_extended_class_threshold():
    # The 0.1 s case moved to the extended reflector: its halo is wider than
    # the 4 km at which the D_length scenes become extended.
    assert fd.halo_extent_m(0.1, 5.0e6, 5.0e5) > 4.0e3
    assert fd.halo_extent_m(0.5, 5.0e6, 5.0e5) < 4.0e3


# --------------------------------------------------------------------------
# Validation
# --------------------------------------------------------------------------


@pytest.mark.parametrize(
    "kw",
    [
        dict(tau_i_s=-1.0),
        dict(k_factor=-0.5),
        dict(tau_d_s=0.0),
        dict(kind="lorentz"),
        dict(coherence_bw_hz=-1.0),
        dict(tau_log_sigma=-0.1),
    ],
)
def test_bad_parameters_are_refused(kw):
    with pytest.raises(ValueError):
        fd.FadingModel(**kw)


def test_gaussian_needs_a_uniform_grid():
    t = np.array([0.0, 0.02, 0.05, 0.07])
    with pytest.raises(ValueError):
        fd.FadingModel(tau_i_s=0.5, kind="gauss").realise(t, np.random.default_rng(0))
    # The exponential kind takes any increasing grid.
    g = fd.FadingModel(tau_i_s=0.5, kind="exp").realise(t, np.random.default_rng(0))
    assert g.shape == (4,)
