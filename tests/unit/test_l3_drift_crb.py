"""WP5 Phase 3 part 2: the Cramer-Rao harness for route one (l3_signal/novel/drift_crb.py).

Anchors: the Doppler-sample model of the degeneracy note (0.42 m/s with t0
known, singular with it free); the tone bound with t0 known against Stein's
FDOA accuracy as Moss (2014) states it; the information matrices against
finite differences of the full model (the phase terms that the harness
leaves to phi0 kept in).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from tss.constants import C_LIGHT
from tss.l3_signal.novel import drift_crb as D
from tss.l3_signal.novel.route_one_estimator import HALF_POWER_X, envelope_power

F, P, MU, W, VS = 4.988e6, 529.8e3, 0.9537, 0.926, 7600.0
T = (np.arange(209) - 104) / 45.0
KW = {"width_s": W, "f_hz": F, "mu_sat": MU, "p_prime_m": P, "v_s_mps": VS}


def test_doppler_samples_reproduce_the_degeneracy_note():
    t = np.linspace(-2.5, 2.5, 50)
    known = D.route_one_crb(D.fisher_doppler_samples(t, sigma_fd_hz=0.1, f_hz=5e6, height_m=500e3, t0_known=True))
    assert known.sigma["v"] == pytest.approx(C_LIGHT / (2 * 5e6) * 0.1 / math.sqrt(50), rel=1e-9)
    assert known.sigma["v"] == pytest.approx(0.42, abs=0.005)
    free = D.route_one_crb(D.fisher_doppler_samples(t, sigma_fd_hz=0.1, f_hz=5e6, height_m=500e3))
    assert free.condition > 1e12 and math.isinf(free.sigma["v"])


@pytest.mark.parametrize("snr_db", [10.0, 20.0])
def test_with_t0_known_the_drift_bound_is_steins_tone_bound(snr_db):
    r = D.route_one_crb(D.fisher_route_one(T, snr_db=snr_db, t0_known=True, width_known=True, **KW))
    amp = 10 ** (snr_db / 20) * np.sqrt(envelope_power(T, 0.0, W))
    want = C_LIGHT / (2 * F * MU) * D.stein_fdoa_sigma_hz(T, amp)
    assert r.sigma["v"] == pytest.approx(want, rel=1e-6)


def test_steins_bound_for_a_flat_envelope_is_the_budgets_closed_form():
    """sqrt(6) / (2 pi T sqrt(rho N)) for N samples of SNR rho over T (large N)."""
    t = np.arange(2001) / 45.0
    rho = 10.0
    got = D.stein_fdoa_sigma_hz(t, np.full(t.size, math.sqrt(rho)))
    span = t.size / 45.0
    assert got == pytest.approx(math.sqrt(6) / (2 * math.pi * span * math.sqrt(rho * t.size)), rel=2e-3)


def test_phase_alone_cannot_separate_t0_and_drift():
    """An envelope far wider than the window carries no t0 information."""
    r = D.route_one_crb(D.fisher_route_one(T, snr_db=20.0, width_s=1e6, f_hz=F, mu_sat=MU, p_prime_m=P))
    assert r.condition > 1e12 and math.isinf(r.sigma["v"])
    r2 = D.route_one_crb(D.fisher_route_one(T, snr_db=20.0, **KW))
    assert r2.condition < 1e4 and r2.correlation_t0_v < -0.9


def _full_mean(theta):
    a_, phi0, t0, v, w = theta
    s = 2 * F * VS**2 / (C_LIGHT * P)
    tzc = t0 + v * MU * P / VS**2
    x = 2 * HALF_POWER_X * (T - t0) / w
    return a_ * np.sinc(x / np.pi) * np.exp(1j * (phi0 + 2 * math.pi * s * (tzc * T - 0.5 * T * T)))


def test_deterministic_matrix_against_finite_differences_of_the_full_model():
    theta = np.array([math.sqrt(10.0), 0.3, 0.05, 25.0, W])
    steps = np.array([1e-6, 1e-6, 1e-7, 1e-5, 1e-7])
    jac = np.stack([(_full_mean(theta + np.eye(5)[i] * steps[i]) - _full_mean(theta - np.eye(5)[i] * steps[i]))
                    / (2 * steps[i]) for i in range(5)], axis=1)
    fis = 2 * np.real(jac.conj().T @ jac)
    want = math.sqrt(np.linalg.inv(fis)[3, 3])
    got = D.route_one_crb(D.fisher_route_one(T, snr_db=10.0, t0_s=0.05, drift_mps=25.0, **KW))
    assert got.sigma["v"] == pytest.approx(want, rel=1e-4)


def test_fading_matrix_against_finite_differences_of_the_full_model():
    tau = 0.5
    rr = np.exp(-np.abs(T[:, None] - T[None, :]) / tau)

    def cov(theta):
        p, t0, v, w = theta
        d = _full_mean((math.sqrt(p), 0.0, t0, v, w))
        return d[:, None] * rr * d.conj()[None, :] + np.eye(T.size)

    theta = np.array([10.0, 0.05, 25.0, W])
    steps = np.array([1e-5, 1e-7, 1e-5, 1e-7])
    c0inv = np.linalg.inv(cov(theta))
    dc = [c0inv @ ((cov(theta + np.eye(4)[i] * steps[i]) - cov(theta - np.eye(4)[i] * steps[i])) / (2 * steps[i]))
          for i in range(4)]
    fis = np.array([[np.real(np.trace(dc[i] @ dc[j])) for j in range(4)] for i in range(4)])
    want = math.sqrt(np.linalg.inv(fis)[2, 2])
    got = D.route_one_crb(D.fisher_route_one(T, snr_db=10.0, t0_s=0.05, drift_mps=25.0, tau_c_s=tau, **KW))
    assert got.sigma["v"] == pytest.approx(want, rel=1e-3)


def test_fading_raises_the_bound_and_a_shorter_correlation_raises_it_more():
    det = D.route_one_crb(D.fisher_route_one(T, snr_db=10.0, **KW)).sigma["v"]
    slow = D.route_one_crb(D.fisher_route_one(T, snr_db=10.0, tau_c_s=2.0, **KW)).sigma["v"]
    fast = D.route_one_crb(D.fisher_route_one(T, snr_db=10.0, tau_c_s=0.1, **KW)).sigma["v"]
    assert det < slow < fast
