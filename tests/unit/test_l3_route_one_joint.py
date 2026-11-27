"""WP5 Phase 3 part 3: the joint route-one estimators (l3_signal/novel/route_one_joint.py).

Parameter level, on the straight scene of ``docs/drift-tier.md`` section 9 and
on the ``c3_verify`` fields.  The deterministic fit is held to the
deterministic Cramer-Rao bound, the fading-aware fit to the Rayleigh-fading
bound (``novel/drift_crb.py``).  The tier's use of them is tested in
``test_l3_drift.py``.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from tss.l3_signal.novel import drift_crb as D
from tss.l3_signal.novel import route_one_joint as J
from tss.l3_signal.novel.route_one_estimator import estimate_route_one
from tss.studies import c3_verify as cv

F, P, MU, W, VS = 4.988e6, 529.8e3, 0.9537, 0.926, 7600.0
T = (np.arange(209) - 104) / 45.0
KW = {"f_hz": F, "group_path_m": P, "mu_sat": MU, "v_s_mps": VS}
DRIFT = 25.0


def _noise(rng, n):
    """Complex white noise of unit variance."""
    return (rng.standard_normal(n) + 1j * rng.standard_normal(n)) / math.sqrt(2)


def _bound(snr_db, tau_c_s=None):
    return D.route_one_crb(D.fisher_route_one(T, snr_db=snr_db, width_s=W, f_hz=F, mu_sat=MU, p_prime_m=P,
                                              v_s_mps=VS, drift_mps=DRIFT, tau_c_s=tau_c_s)).sigma["v"]


def _fit(kind, t, y, e, noise_var=1.0, **kw):
    start = {"t0_start_s": e.t_power_peak_s, "drift_start_mps": e.drift_mps}
    if kind == "joint":
        return J.estimate_route_one_joint(t, y, **start, **kw)
    return J.estimate_route_one_joint_fading(t, y, noise_var=noise_var, **start, **kw)


def _monte_carlo(kind, snr_db, tau_c_s, seeds):
    """Fits on the echoes the pinned estimator accepts (decision (s))."""
    amp = 10 ** (snr_db / 20)
    lc = None
    if tau_c_s is not None:
        lag = np.abs(T[:, None] - T[None, :])
        lc = np.linalg.cholesky(np.exp(-lag / tau_c_s) + 1e-12 * np.eye(T.size))
    rows = []
    for seed in seeds:
        rng = np.random.default_rng(seed)
        g = 1.0 if lc is None else lc @ _noise(rng, T.size)
        y = amp * J.joint_model(T, t0_s=0.0, drift_mps=DRIFT, width_s=W, **KW) * g + _noise(rng, T.size)
        e = estimate_route_one(T, y, **KW)
        if not e.separable:
            continue
        r = _fit(kind, T, y, e, **KW)
        assert r.ok, r.reason
        rows.append((r.drift_mps - DRIFT, r.drift_sigma_mps, e.drift_mps - DRIFT, r.tau_c_fit_s))
    return np.array(rows)


def _rms(x):
    return float(np.sqrt(np.mean(np.square(x))))


def _robust_sd(x):
    return float(1.4826 * np.median(np.abs(x - np.median(x))))


def test_the_model_is_the_bounds_model():
    """``joint_model`` has drift_crb's envelope and, up to a constant phase, its phase."""
    s = J.joint_model(T, t0_s=0.05, drift_mps=DRIFT, width_s=W, **KW)
    m = D._model(T, amp=1.0, t0=0.05, v=DRIFT, width=W,
                 s_rate=2 * F * VS**2 / (299792458.0 * P), kappa=MU * P / VS**2)
    np.testing.assert_allclose(s, m["a"] * np.exp(1j * m["phi"]), atol=1e-12)


@pytest.mark.parametrize("kind", ["joint", "joint_fading"])
def test_noise_free_model_field_is_recovered(kind):
    y = (0.8 - 0.3j) * J.joint_model(T, t0_s=0.07, drift_mps=-40.0, width_s=1.1, **KW)
    e = estimate_route_one(T, y, **KW)
    r = _fit(kind, T, y, e, noise_var=1e-12, **KW)
    assert r.ok, r.reason
    assert r.drift_mps == pytest.approx(-40.0, abs=1e-3 if kind == "joint" else 0.5)
    assert r.t0_s == pytest.approx(0.07, abs=1e-5 if kind == "joint" else 5e-3)
    assert r.width_s == pytest.approx(1.1, abs=1e-4 if kind == "joint" else 0.05)
    assert r.t_zc_s == pytest.approx(r.t0_s + r.drift_mps * MU * P / VS**2, abs=1e-12)
    if kind == "joint":
        assert r.amplitude == pytest.approx(0.8 - 0.3j, abs=1e-5)


def _c3_fields(monkeypatch, group, every=True):
    """``(trial, t, field, kwargs, pinned)`` for c3_verify's trials of one case.

    With ``every`` false, one trial per environment (A) or per environment and
    tilt at its first drift (B, 0.5 degrees only): the fading fit takes seconds
    on the long windows, and the full grids run in the slow suite.
    """
    seen = []

    def capture(t, field, **kw):
        est = estimate_route_one(t, field, **kw)
        seen.append((np.array(t), np.array(field), kw, est))
        return est

    monkeypatch.setattr(cv, "estimate_route_one", capture)
    out = []
    env = {e.name: e for e in cv.ENVIRONMENTS}
    rng = np.random.default_rng(0)
    if group == "A":
        cases = [(e.name, 0.0, x, v) for e in cv.ENVIRONMENTS for x in (0.0, 5.0e3)
                 for v in (-50.0, 0.0, 25.0, 100.0)]
    else:
        cases = [(name, tilt, 0.0, v)
                 for name in ("free_5MHz_500km", "free_3MHz_500km", "plasma_5MHz_mu0.5", "plasma_3MHz_mu0.5")
                 for tilt in (0.1, 0.25, 0.5, 1.0) for v in (0.0, 25.0)]
    if not every:
        keep = (lambda c: c[2] == 5.0e3 and c[3] == 25.0) if group == "A" else (
            lambda c: c[1] == 0.5 and c[3] == 25.0)
        cases = [c for c in cases if keep(c)]
    for name, tilt, x, v in cases:
        trial = cv.run_trial(group, env[name], tilt_deg=tilt, x_m=x, drift_mps=v, rng=rng)
        out.append((trial, *seen[-1]))
    return out


def _grid(n_fast, n_full):
    """The deterministic fit on the whole grid; the fading fit on a subset, and whole when slow."""
    return [pytest.param("joint", True, n_full, id="joint"),
            pytest.param("joint_fading", False, n_fast, id="joint_fading"),
            pytest.param("joint_fading", True, n_full, id="joint_fading-all", marks=pytest.mark.slow)]


@pytest.mark.parametrize(("kind", "every", "count"), _grid(7, 56))
def test_c3_verify_untilted_facets_within_1_mps(monkeypatch, kind, every, count):
    """c3_verify case A, noise-free: the forward model does not assume the
    ``sinc^2`` envelope's exactness; 56 of 56 within 1 m/s."""
    trials = _c3_fields(monkeypatch, "A", every)
    assert len(trials) == count
    err = []
    for trial, t, field, kw, est in trials:
        assert est.separable
        peak = float(np.max(np.abs(field)))
        r = _fit(kind, t, field / peak, est, noise_var=1e-12, **kw)
        assert r.ok, (trial.environment, r.reason)
        err.append(r.drift_mps - trial.drift_true_mps)
    assert np.max(np.abs(err)) < 1.0


@pytest.mark.parametrize(("kind", "every", "count"), _grid(4, 32))
def test_c3_verify_tilt_bias_is_reproduced_within_1_percent(monkeypatch, kind, every, count):
    """c3_verify case B: a tilted facet reads ``-v_s mu_f sin(a) / mu_sat`` too."""
    trials = _c3_fields(monkeypatch, "B", every)
    assert len(trials) == count
    for trial, t, field, kw, est in trials:
        assert est.separable
        r = _fit(kind, t, field / np.max(np.abs(field)), est, noise_var=1e-12, **kw)
        assert r.ok, (trial.environment, r.reason)
        got = r.drift_mps - trial.drift_true_mps
        assert got == pytest.approx(trial.predicted_bias_mps, rel=0.01), (trial.environment, trial.tilt_deg)


@pytest.mark.parametrize("snr_db", [20.0, 10.0])
def test_deterministic_fit_meets_the_bound_and_states_its_spread(snr_db):
    """60 seeds: RMS error within 1.25 of the deterministic bound (exit: 1.2,
    with the sampling error of 60 seeds), well under the pinned estimator's,
    and the stated sigma matching the spread."""
    r = _monte_carlo("joint", snr_db, None, range(100, 160))
    err, sig, pinned = r[:, 0], r[:, 1], r[:, 2]
    assert err.size >= 20
    assert 0.75 < _rms(err) / _bound(snr_db) < 1.25
    assert _rms(err) < 0.85 * _rms(pinned)
    assert 0.75 < np.median(sig) / np.std(err) < 1.3


@pytest.mark.parametrize("tau_c_s", [0.5, 2.0])
def test_fading_fit_is_near_the_fading_bound_at_20_db(tau_c_s):
    """30 seeds, Rayleigh fading: the robust spread within 1.5 of the fading
    bound, the RMS under half the pinned estimator's, tau_c recovered to a
    factor of two in the median, and the stated sigma near the spread."""
    r = _monte_carlo("joint_fading", 20.0, tau_c_s, range(200, 230))
    err, sig, pinned, tau = r[:, 0], r[:, 1], r[:, 2], r[:, 3]
    assert err.size >= 15
    bound = _bound(20.0, tau_c_s)
    assert _robust_sd(err) < 1.5 * bound
    assert _rms(err) < 0.5 * _rms(pinned)
    assert 0.5 < np.median(tau) / tau_c_s < 2.0
    assert 0.6 < np.nanmedian(sig) / _robust_sd(err) < 1.6


def test_the_fading_fit_floors_the_noise_at_the_models_own_error():
    y = 2.0 * J.joint_model(T, t0_s=0.0, drift_mps=DRIFT, width_s=W, **KW)
    low = J.estimate_route_one_joint_fading(T, y, noise_var=1e-12, **KW)
    assert low.noise_var == pytest.approx(J.MODEL_ERROR_FLOOR * np.max(np.abs(y) ** 2), rel=1e-12)
    high = J.estimate_route_one_joint_fading(T, y, noise_var=0.5, **KW)
    assert high.noise_var == 0.5


def test_fading_likelihood_is_the_dense_gaussian_one():
    """The O(N) form (the exponential correlation's inverse is tridiagonal)
    against ``log det C + y^H C^-1 y`` formed directly, uneven times included."""
    rng = np.random.default_rng(4)
    t = np.sort(T + rng.uniform(-0.004, 0.004, T.size))
    for tau, power, nv in ((0.05, 10.0, 1.0), (0.5, 100.0, 0.7), (30.0, 3.0, 0.01)):
        y = 3.0 * _noise(rng, t.size)
        s = J.joint_model(t, t0_s=0.1, drift_mps=60.0, width_s=0.8, **KW)
        c = power * s[:, None] * np.exp(-np.abs(t[:, None] - t[None, :]) / tau) * np.conj(s)[None, :]
        c = c + nv * np.eye(t.size)
        dense = np.linalg.slogdet(c)[1] + np.real(np.vdot(y, np.linalg.solve(c, y)))
        fast = J.fading_nll(t, y, power=power, t0_s=0.1, drift_mps=60.0, width_s=0.8, tau_c_s=tau,
                            noise_var=nv, **KW)
        assert fast == pytest.approx(dense, rel=1e-9)


def test_noise_variance_from_the_detectors_bin_and_the_gated_sum():
    """Noise only: a Hann-windowed bin's mean power over ``sum h^2`` is the
    sample variance, and the sum of N gated samples has N times it."""
    rng = np.random.default_rng(8)
    sig2, m, n_gate = 3.7, 800, 600
    x = math.sqrt(sig2) * np.stack([_noise(rng, m) for _ in range(400)])
    per_bin = float(np.median(np.abs(np.fft.fft(x * np.hanning(m), axis=1)) ** 2)) / math.log(2.0)
    got = J.sample_noise_variance(per_bin, m)
    assert got == pytest.approx(sig2, rel=0.02)
    z = x[:, :n_gate] @ np.exp(-2j * math.pi * 0.013 * np.arange(n_gate))
    assert np.mean(np.abs(z) ** 2) == pytest.approx(n_gate * got, rel=0.15)


def test_two_echoes_with_different_drifts_are_ambiguous():
    """Two facets of equal strength apart in time: the likelihood has two
    well-separated maxima of nearly the same height, and the fit says so."""
    t = T
    a = J.joint_model(t, t0_s=-1.3, drift_mps=60.0, width_s=0.5, **KW)
    b = J.joint_model(t, t0_s=1.3, drift_mps=-60.0, width_s=0.5, **KW)
    y = 10.0 * (a + b) + _noise(np.random.default_rng(1), t.size)
    r = J.estimate_route_one_joint_fading(t, y, noise_var=1.0, **KW)
    assert r.ok, r.reason
    assert 0.0 <= r.nll_margin < 5.0
    assert abs(r.drift_alt_mps - r.drift_mps) > 3.0 * r.drift_sigma_mps


def test_one_strong_echo_is_not_ambiguous():
    y = 10.0 * J.joint_model(T, t0_s=0.1, drift_mps=DRIFT, width_s=W, **KW)
    y = y + _noise(np.random.default_rng(2), T.size)
    r = J.estimate_route_one_joint_fading(T, y, noise_var=1.0, **KW)
    assert r.ok, r.reason
    assert r.nll_margin > 10.0
    # No second maximum is reported as a drift when there is none.
    assert math.isfinite(r.drift_alt_mps) == math.isfinite(r.nll_margin)


def test_the_separation_that_makes_a_second_maximum_is_stated_in_sigmas():
    t = T
    a = J.joint_model(t, t0_s=-1.3, drift_mps=60.0, width_s=0.5, **KW)
    b = J.joint_model(t, t0_s=1.3, drift_mps=-60.0, width_s=0.5, **KW)
    y = 10.0 * (a + b) + _noise(np.random.default_rng(1), t.size)
    far = J.estimate_route_one_joint_fading(t, y, noise_var=1.0, separation_sigmas=1e6, **KW)
    assert math.isinf(far.nll_margin) and math.isnan(far.drift_alt_mps)


def test_the_fits_say_when_they_cannot_run():
    y = J.joint_model(T, t0_s=0.0, drift_mps=DRIFT, width_s=W, **KW)
    r = J.estimate_route_one_joint(T, y, t0_start_s=math.nan, drift_start_mps=0.0, **KW)
    assert not r.ok and math.isnan(r.drift_mps) and r.reason
    r = J.estimate_route_one_joint_fading(T, y, noise_var=math.nan, **KW)
    assert not r.ok and math.isnan(r.drift_mps) and "noise" in r.reason
    far = J.joint_model(T, t0_s=0.0, drift_mps=900.0, width_s=W, **KW)
    r = J.estimate_route_one_joint_fading(T, 10 * far, noise_var=1.0, v_max_mps=100.0, **KW)
    assert not r.ok and "bound" in r.reason
