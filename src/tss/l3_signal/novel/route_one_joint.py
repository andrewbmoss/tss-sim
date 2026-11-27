"""Layer: L3
WP: WP5
Claims: C3
Status: partial (parameter level and the drift tier; not yet the default)
Requirements: SCI-15

Joint estimators for route one (WP5 Phase 3 part 3, decisions (o), (r) and
(s), 2026-11-27, and the fading-aware fit added the same day).

WHY
---
The pinned estimator (:func:`estimate_route_one`) splits the problem: the
overflight time ``t_pk`` from a fit to the power alone, the zero crossing
``t_zc`` from a power-weighted mean of phase steps.  Its RMS drift error was
1.6 and 2.8 times the Cramer-Rao bound at 20 and 10 dB per sweep
(``docs/drift-tier.md`` section 9).  The power fit ignores the phase, and the
phase steps ignore what the envelope says about ``t0``, which is what makes
drift and overflight time separable.

THE MODEL
---------
That of :mod:`tss.l3_signal.novel.drift_crb`:

    y_n = c s_n g_n + w_n,
    s_n = sinc(x_n) exp(j 2 pi S [kappa v u_n - u_n^2 / 2]),
    u_n = t_n - t0,   x_n = 2 X u_n / W,
    S = 2 f v_s^2 / (c P'),   kappa = mu_sat P' / v_s^2,

so ``t_zc - t0 = kappa v`` and ``v`` is the drift.  ``S`` is the ephemeris
slope (decision (s)): the pinned free-slope test stays the gate, and both
estimators here form a drift only on echoes the pinned one accepted.

DETERMINISTIC (:func:`estimate_route_one_joint`)
------------------------------------------------
``g_n = 1`` and ``w`` white.  Maximum likelihood is least squares on the
complex samples.  The complex amplitude ``c`` is linear and solved in closed
form at each trial ``(t0, v, W)`` (variable projection), so the search has
three parameters, started from the pinned estimate.  The covariance is
``s^2 (J^T J)^-1`` with ``J`` the Jacobian of the projected residual and
``s^2`` the residual variance per real component.  It is the covariance under
this model: under fading the model is wrong and the sigma is too small.

FADING-AWARE (:func:`estimate_route_one_joint_fading`)
------------------------------------------------------
``g`` is zero-mean complex Gaussian with ``E g_n g_m* = exp(-|t_n - t_m| /
tau_c)``, the Rayleigh gain of the fading bound, so

    y ~ CN(0, P D R D^H + sigma^2 I),   D = diag(s_n),

with ``P``, ``t0``, ``v``, ``W`` and ``tau_c`` fitted by maximum likelihood.

* ``sigma^2`` is **given**, not fitted: the drift tier derives it from the
  detector's noise reference (:func:`sample_noise_variance`).  With it free, a
  fit can explain a fading spike as a narrow facet plus more noise.  It is not
  taken below :data:`MODEL_ERROR_FLOOR` of the largest sample power, the
  envelope model's own error: on an echo that does not fade, a much smaller
  value makes the fit explain that error as fading.
* The exponential correlation is that of a first-order Markov process, so
  ``R^-1`` is tridiagonal and each likelihood costs ``O(N)``
  (:func:`fading_nll`), not ``O(N^3)``.
* The likelihood has several minima.  It is searched on a grid of ``(t0, v, W,
  tau_c)`` over the whole window and the drift prior, and a local fit is run
  from the best ``starts_per_width`` distinct cells at each width and from the
  pinned estimate; the lowest wins.  The global maximum itself can still be far
  from the truth, which no search removes (``docs/drift-tier.md`` sections 11
  and 12).
* The covariance is the inverse of the observed information (the Hessian of
  the negative log likelihood) over the parameters not held at a bound.
* **A second maximum** (the ambiguity check, decision 2B of 2026-11-27).  Among
  the local fits' answers, those whose drift is more than ``separation_sigmas``
  (3) of the best fit's sigma from the best drift are other maxima of the
  likelihood.  ``nll_margin`` is how far the lowest of them lies above the best
  in negative log likelihood (a margin of 2 means ``e^-2`` as likely), and
  ``drift_alt_mps`` its drift; infinite and NaN when there is none, and NaN
  both when the best fit has no sigma.  It sees only the maxima the starts
  reach.  Whether a small margin declines the echo is the caller's choice
  (``l3/drift``, ``ambiguity_margin``).

LIMITS
------
The envelope is the ``sinc^2`` power envelope that the pinned fit and the
bound assume; a facet whose envelope differs is a model error shared by all
three.  The fading is Rayleigh (no specular part) with an exponential
correlation; a Rician gain or another correlation shape is a mismatch that has
not been measured.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from tss.constants import C_LIGHT
from tss.l3_signal.novel.route_one_estimator import HALF_POWER_X, fit_overflight_time

__all__ = [
    "JointEstimate",
    "joint_model",
    "sample_noise_variance",
    "fading_nll",
    "estimate_route_one_joint",
    "estimate_route_one_joint_fading",
]

#: Least noise variance the fading fit assumes, relative to the largest sample
#: power: the ``sinc^2`` envelope model's own error.  On the noise-free
#: ``c3_verify`` fields the deterministic fit leaves a residual of 1.6e-4 to
#: 2.7e-3 of the peak.  Told that much less is noise, the fading fit explains
#: the misfit as fading and moves ``t0`` and the width: at 1e-6, 45 of 56
#: untilted facets miss by more than 1 m/s; at 1e-4 one tilted facet of 32
#: misses its bias by 30 %; from 1e-3 up all 88 read as they should
#: (``docs/drift-tier.md`` section 11).
MODEL_ERROR_FLOOR = 3.0e-3

#: The search grid: largest ``t0`` step, number of widths, decorrelation times,
#: and the drift step as a fraction of the Doppler resolution.
GRID_T0_STEP_S = 0.1
GRID_WIDTHS = 6
GRID_TAUS_S = (0.1, 0.5, 2.0)
GRID_V_FRACTION = 0.5


@dataclass(frozen=True)
class JointEstimate:
    """What a joint fit returns for one echo.

    ``ok`` is false when no fit converged; the numbers are then NaN and
    ``reason`` says why.  ``noise_var`` is ``E |w|^2`` per sample, in the
    field's units: estimated from the residual (deterministic) or as given
    (fading).  ``power``, ``tau_c_fit_s`` (the fitted time of the exponential
    correlation), ``nll``, ``nll_margin`` and ``drift_alt_mps`` (module
    docstring, a second maximum) belong to the fading fit.
    ``at_bound`` names the parameters the fading fit left on a bound; their
    uncertainty is not in the sigmas.
    """

    ok: bool
    reason: str
    drift_mps: float = math.nan
    drift_sigma_mps: float = math.nan
    t0_s: float = math.nan
    t0_sigma_s: float = math.nan
    t_zc_s: float = math.nan
    width_s: float = math.nan
    amplitude: complex = complex(math.nan, math.nan)
    noise_var: float = math.nan
    corr_t0_v: float = math.nan
    power: float = math.nan
    tau_c_fit_s: float = math.nan
    nll: float = math.nan
    at_bound: tuple[str, ...] = ()
    nll_margin: float = math.nan
    drift_alt_mps: float = math.nan


def _rates(f_hz: float, group_path_m: float, mu_sat: float, v_s_mps: float) -> tuple[float, float]:
    """``S`` (Hz/s) and ``kappa`` (s per m/s) of the module docstring."""
    s = 2.0 * float(f_hz) * float(v_s_mps) ** 2 / (C_LIGHT * float(group_path_m))
    return s, float(mu_sat) * float(group_path_m) / float(v_s_mps) ** 2


def joint_model(t_s: np.ndarray, *, t0_s: float, drift_mps: float, width_s: float, f_hz: float,
                group_path_m: float, mu_sat: float, v_s_mps: float = 7600.0) -> np.ndarray:
    """``s_n``, the unit-amplitude model of the module docstring."""
    s_rate, kappa = _rates(f_hz, group_path_m, mu_sat, v_s_mps)
    u = np.asarray(t_s, dtype=float) - float(t0_s)
    x = 2.0 * HALF_POWER_X * u / float(width_s)
    phase = 2.0 * math.pi * s_rate * (kappa * float(drift_mps) * u - 0.5 * u * u)
    return np.sinc(x / math.pi) * np.exp(1j * phase)


def sample_noise_variance(noise_power_per_bin: float, segment_samples: int) -> float:
    """``E |w|^2`` per raw sample from the windowed detector's noise reference.

    The detector's bin is ``sum_k h_k w_k e^{-j ...}`` over ``segment_samples``
    with a Hann window ``h``, so for white noise its mean power is the sample
    variance times ``sum h^2``.  A sum of ``N`` unweighted gated samples (the
    drift tier's one value per sweep) then has variance ``N`` times this.
    """
    h = np.hanning(int(segment_samples))
    return float(noise_power_per_bin) / float(np.sum(h * h))


# -- deterministic -------------------------------------------------------------------

def estimate_route_one_joint(
    t_s: np.ndarray,
    field: np.ndarray,
    *,
    f_hz: float,
    group_path_m: float,
    mu_sat: float,
    v_s_mps: float = 7600.0,
    t0_start_s: float,
    drift_start_mps: float,
    width_start_s: float | None = None,
) -> JointEstimate:
    """The deterministic fit of ``(c, t0, v, W)`` (module docstring).

    ``t0_start_s`` and ``drift_start_mps`` come from the pinned estimator;
    ``width_start_s`` from the ``sinc^2`` power fit if not given.
    """
    from scipy.optimize import least_squares

    t = np.asarray(t_s, dtype=float)
    y = np.asarray(field, dtype=complex)
    dt = float(np.median(np.diff(t)))
    span = float(t[-1] - t[0])
    if width_start_s is None:
        _, width_start_s = fit_overflight_time(t, np.abs(y) ** 2)
    x0 = np.array([float(t0_start_s), float(drift_start_mps), float(width_start_s)])
    if not np.all(np.isfinite(x0)):
        return JointEstimate(False, "no finite start from the pinned estimator")
    lo = np.array([t[0], -np.inf, 2.0 * dt])
    hi = np.array([t[-1], np.inf, span])
    x0 = np.clip(x0, lo + 1e-9, hi - 1e-9)
    kw = {"f_hz": f_hz, "group_path_m": group_path_m, "mu_sat": mu_sat, "v_s_mps": v_s_mps}

    def project(q: np.ndarray) -> tuple[complex, np.ndarray]:
        s = joint_model(t, t0_s=q[0], drift_mps=q[1], width_s=q[2], **kw)
        ss = float(np.real(np.vdot(s, s)))
        c = complex(np.vdot(s, y) / ss) if ss > 0.0 else 0j
        return c, y - c * s

    def resid(q: np.ndarray) -> np.ndarray:
        r = project(q)[1]
        return np.concatenate([r.real, r.imag])

    # Steps and tolerances on the parameters' natural scales, not absolute ones
    # (the lesson of part 2's envelope fit, docs/drift-tier.md section 6).
    sol = least_squares(resid, x0=x0, bounds=(lo, hi), x_scale=np.array([0.05, 5.0, 0.1]), jac="3-point",
                        xtol=1e-10, ftol=1e-12, gtol=1e-12, method="trf")
    if not sol.success:
        return JointEstimate(False, f"fit did not converge: {sol.message}")
    t0, v, w = (float(a) for a in sol.x)
    c, r = project(sol.x)
    s2 = float(np.sum(np.abs(r) ** 2) / (2 * t.size - 5))  # per real component; c is two parameters
    jtj = sol.jac.T @ sol.jac
    try:
        cov = s2 * np.linalg.inv(jtj)
    except np.linalg.LinAlgError:
        return JointEstimate(False, "singular information matrix")
    if not np.all(np.isfinite(cov)) or np.any(np.diag(cov) < 0.0):
        return JointEstimate(False, "covariance not positive")
    sig_t0, sig_v = math.sqrt(cov[0, 0]), math.sqrt(cov[1, 1])
    _, kappa = _rates(f_hz, group_path_m, mu_sat, v_s_mps)
    corr = float(cov[0, 1] / (sig_t0 * sig_v)) if sig_t0 > 0.0 and sig_v > 0.0 else math.nan
    return JointEstimate(
        True, "converged", drift_mps=v, drift_sigma_mps=sig_v, t0_s=t0, t0_sigma_s=sig_t0,
        t_zc_s=t0 + kappa * v, width_s=w, amplitude=c, noise_var=2.0 * s2, corr_t0_v=corr,
    )


# -- fading-aware ----------------------------------------------------------------------

def _markov_precision(dt: np.ndarray, tau_c_s: float) -> tuple[np.ndarray, np.ndarray, float]:
    """Diagonal, off-diagonal and log determinant for ``R_nm = exp(-|t_n - t_m| / tau)``.

    ``R`` is the covariance of a first-order Markov process, so ``R^-1`` is
    tridiagonal.  With ``r_k = exp(-dt_k / tau)`` and ``q_k = 1 - r_k^2``:
    ``det R = prod q_k``.
    """
    r = np.exp(-dt / float(tau_c_s))
    q = -np.expm1(-2.0 * dt / float(tau_c_s))
    diag = np.ones(dt.size + 1)
    diag[1:] = 1.0 / q
    diag[:-1] += r * r / q
    return diag, -r / q, float(np.sum(np.log(q)))


def _nll(dt: np.ndarray, a: np.ndarray, z: np.ndarray, power: float, tau_c_s: float,
         noise_var: float, energy: float) -> np.ndarray:
    """``log det C + y^H C^-1 y`` for each column of ``z`` (module docstring).

    ``a`` is the real envelope ``sinc(x_n)``, ``z = e^{-j phi} y`` (``N`` or
    ``N x K``, one column per trial drift) and ``energy = sum |y|^2``.  With
    ``M = sigma^2 I + P A R A``, ``T = R^-1 / P + A^2 / sigma^2``:
    ``log det M = N log(sigma^2 P) + log det R + log det T`` and
    ``y^H C^-1 y = |y|^2 / sigma^2 - b^H T^-1 b / sigma^4``, ``b = A z``.
    """
    from scipy.linalg import cho_solve_banded, cholesky_banded

    diag, off, logdet_r = _markov_precision(dt, tau_c_s)
    ab = np.empty((2, a.size))
    ab[0, 0] = 0.0
    ab[0, 1:] = off / power
    ab[1] = diag / power + a * a / noise_var
    u = cholesky_banded(ab, check_finite=False)
    b = (a * z.T).T.reshape(a.size, -1)
    # T is real: solve the real and imaginary parts as real right-hand sides.
    parts = np.concatenate([b.real, b.imag], axis=1)
    solved = cho_solve_banded((u, False), parts, check_finite=False)
    bhtb = np.sum(parts * solved, axis=0).reshape(2, -1).sum(axis=0).reshape(z.shape[1:])
    quad = energy / noise_var - bhtb / noise_var**2
    logdet = a.size * math.log(noise_var * power) + logdet_r + 2.0 * float(np.sum(np.log(u[1])))
    return logdet + quad


def fading_nll(t_s: np.ndarray, field: np.ndarray, *, power: float, t0_s: float, drift_mps: float,
               width_s: float, tau_c_s: float, noise_var: float, f_hz: float, group_path_m: float,
               mu_sat: float, v_s_mps: float = 7600.0) -> float:
    """The fading model's negative log likelihood, up to the constant ``N log pi``."""
    t = np.asarray(t_s, dtype=float)
    y = np.asarray(field, dtype=complex)
    s_rate, kappa = _rates(f_hz, group_path_m, mu_sat, v_s_mps)
    u = t - float(t0_s)
    a = np.sinc(2.0 * HALF_POWER_X * u / float(width_s) / math.pi)
    z = y * np.exp(-2j * math.pi * s_rate * (kappa * float(drift_mps) * u - 0.5 * u * u))
    return float(_nll(np.diff(t), a, z, float(power), float(tau_c_s), float(noise_var),
                      float(np.sum(np.abs(y) ** 2))))


def _grid_starts(t: np.ndarray, y: np.ndarray, *, s_rate: float, kappa: float, noise_var: float,
                 v_max: float, starts_per_width: int) -> list[tuple[float, float, float, float, float]]:
    """The best distinct grid cells at each width, as ``(t0, v, W, P, tau_c)``.

    The likelihood's minima differ in width as much as in time, so each width
    contributes its own starts: the best cells ranked over the whole grid alone
    are nearly always wide ones, and the local fit from them stays wide.

    ``t0`` over the whole window, ``v`` over ``+/- v_max``, the widths and
    decorrelation times of the ``GRID_*`` constants; ``P`` by weighted least
    squares on the power.  The ``v`` step is ``GRID_V_FRACTION`` of the Doppler
    resolution of 3 s or of the window, whichever is shorter.
    """
    dt = np.diff(t)
    step_t = float(np.median(dt))
    span = float(t[-1] - t[0])
    energy = float(np.sum(np.abs(y) ** 2))
    t0s = np.linspace(t[0], t[-1], int(math.ceil(span / GRID_T0_STEP_S)) + 1)
    hz_per_mps = s_rate * kappa  # 2 f mu_sat / c
    dv = float(np.clip(GRID_V_FRACTION / (hz_per_mps * min(span, 3.0)), 0.5, 10.0))
    vs = np.linspace(-v_max, v_max, 2 * int(math.ceil(v_max / dv)) + 1)
    widths = np.geomspace(max(4.0 * step_t, 0.1 * min(span, 3.0)), min(span, 3.0), GRID_WIDTHS)
    p2 = np.abs(y) ** 2 - noise_var
    cells = []
    for t0 in t0s:
        u = t - t0
        z = (y * np.exp(1j * math.pi * s_rate * u * u))[:, None] * np.exp(
            -2j * math.pi * hz_per_mps * np.outer(u, vs))
        for w in widths:
            a = np.sinc(2.0 * HALF_POWER_X * u / w / math.pi)
            a2 = a * a
            power = max(float(np.dot(a2, p2) / np.dot(a2, a2)), noise_var)
            for tau in GRID_TAUS_S:
                f = _nll(dt, a, z, power, tau, noise_var, energy)
                j = int(np.argmin(f))
                cells.append((float(f[j]), float(t0), float(vs[j]), float(w), power, float(tau)))
    cells.sort()
    out: list[tuple[float, float, float, float, float]] = []
    for w_i in widths:
        picked: list[tuple[float, float, float, float, float]] = []
        for _, t0, v, w, power, tau in cells:
            if len(picked) >= starts_per_width:
                break
            # Distinct: not the same hollow in (t0, zero crossing).
            if w == w_i and all(abs(t0 - q[0]) > 0.25 or abs((t0 + kappa * v) - (q[0] + kappa * q[1])) > 0.25
                                for q in picked):
                picked.append((t0, v, w, power, tau))
        out.extend(picked)
    return out


def _hessian(fun, q: np.ndarray, steps: np.ndarray) -> np.ndarray:
    """Central second differences."""
    n = q.size
    h = np.empty((n, n))
    f0 = fun(q)
    for i in range(n):
        ei = np.zeros(n)
        ei[i] = steps[i]
        h[i, i] = (fun(q + ei) - 2.0 * f0 + fun(q - ei)) / steps[i] ** 2
        for j in range(i + 1, n):
            ej = np.zeros(n)
            ej[j] = steps[j]
            h[i, j] = h[j, i] = (fun(q + ei + ej) - fun(q + ei - ej) - fun(q - ei + ej)
                                 + fun(q - ei - ej)) / (4.0 * steps[i] * steps[j])
    return h


def estimate_route_one_joint_fading(
    t_s: np.ndarray,
    field: np.ndarray,
    *,
    f_hz: float,
    group_path_m: float,
    mu_sat: float,
    noise_var: float,
    v_s_mps: float = 7600.0,
    t0_start_s: float | None = None,
    drift_start_mps: float | None = None,
    width_start_s: float | None = None,
    v_max_mps: float = 300.0,
    starts_per_width: int = 6,
    separation_sigmas: float = 3.0,
) -> JointEstimate:
    """Maximum likelihood of ``(P, t0, v, W, tau_c)`` under fading (module docstring).

    ``noise_var`` is ``E |w|^2`` of the field's samples, measured, not fitted;
    below ``MODEL_ERROR_FLOOR`` of the largest sample power the floor is used,
    and the result's ``noise_var`` says what was.
    The pinned estimate, when given, is one more starting point.  ``v_max_mps``
    bounds the drift searched (the drift prior is +/- 250 m/s,
    ``docs/drift-mode.md``).  ``separation_sigmas`` says how far in drift
    another maximum must be to count as one (module docstring).
    """
    from scipy.optimize import minimize

    t = np.asarray(t_s, dtype=float)
    y = np.asarray(field, dtype=complex)
    if not (math.isfinite(noise_var) and noise_var > 0.0):
        return JointEstimate(False, "no noise variance")
    if not np.all(np.isfinite(y)) or t.size < 8:
        return JointEstimate(False, "too few finite samples")
    dt = np.diff(t)
    step_t = float(np.median(dt))
    span = float(t[-1] - t[0])
    s_rate, kappa = _rates(f_hz, group_path_m, mu_sat, v_s_mps)
    energy = float(np.sum(np.abs(y) ** 2))
    noise_var = max(float(noise_var), MODEL_ERROR_FLOOR * float(np.max(np.abs(y) ** 2)))
    peak = max(float(np.max(np.abs(y) ** 2)), noise_var)

    names = ("t0", "v", "log_W", "log_P", "log_tau_c")
    scale = np.array([0.05, 5.0, 0.1, 0.3, 0.3])
    lo = np.array([t[0], -v_max_mps, math.log(2.0 * step_t), math.log(1e-3 * noise_var), math.log(0.5 * step_t)])
    hi = np.array([t[-1], v_max_mps, math.log(span), math.log(100.0 * peak), math.log(100.0)])

    def nll(q: np.ndarray) -> float:
        u = t - q[0]
        a = np.sinc(2.0 * HALF_POWER_X * u / math.exp(q[2]) / math.pi)
        z = y * np.exp(-2j * math.pi * s_rate * (kappa * q[1] * u - 0.5 * u * u))
        return float(_nll(dt, a, z, math.exp(q[3]), math.exp(q[4]), noise_var, energy))

    starts = _grid_starts(t, y, s_rate=s_rate, kappa=kappa, noise_var=noise_var, v_max=v_max_mps,
                          starts_per_width=starts_per_width)
    if t0_start_s is not None and drift_start_mps is not None and math.isfinite(t0_start_s + drift_start_mps):
        if width_start_s is None:
            _, width_start_s = fit_overflight_time(t, np.abs(y) ** 2)
        a2 = np.sinc(2.0 * HALF_POWER_X * (t - t0_start_s) / width_start_s / math.pi) ** 2
        p0 = max(float(np.dot(a2, np.abs(y) ** 2 - noise_var) / np.dot(a2, a2)), noise_var)
        starts.append((float(t0_start_s), float(drift_start_mps), float(width_start_s), p0, 0.5))

    best = None
    found: list[tuple[float, float]] = []  # (nll, drift) of every local fit
    for t0, v, w, power, tau in starts:
        q0 = np.clip(np.array([t0, v, math.log(w), math.log(power), math.log(tau)]), lo, hi)
        sol = minimize(lambda r: nll(r * scale), q0 / scale, method="L-BFGS-B",
                       bounds=list(zip(lo / scale, hi / scale, strict=True)),
                       options={"maxiter": 500, "ftol": 1e-12, "gtol": 1e-6})
        if math.isfinite(sol.fun):
            found.append((float(sol.fun), float(np.clip(sol.x[1] * scale[1], lo[1], hi[1]))))
            if best is None or sol.fun < best.fun:
                best = sol
    if best is None:
        return JointEstimate(False, "no local fit returned a finite likelihood")
    q = np.clip(best.x * scale, lo, hi)
    t0, v = float(q[0]), float(q[1])
    w, power, tau = (float(math.exp(a)) for a in q[2:])

    margin = 1e-6 * (hi - lo)
    free = (q > lo + margin) & (q < hi - margin)
    at_bound = tuple(n for n, f in zip(names, free, strict=True) if not f)
    if "t0" in at_bound or "v" in at_bound:
        return JointEstimate(False, f"the fit ran to a bound in {', '.join(at_bound)}", noise_var=noise_var,
                             nll=float(best.fun), at_bound=at_bound)
    idx = np.nonzero(free)[0]
    steps = np.array([1e-4, 2e-2, 1e-3, 3e-3, 3e-3])[idx]

    def sub(r: np.ndarray) -> float:
        full = q.copy()
        full[idx] = r
        return nll(full)

    h = _hessian(sub, q[idx].copy(), steps)
    sig_t0 = sig_v = corr = math.nan
    reason = "converged"
    eig = np.linalg.eigvalsh(h) if np.all(np.isfinite(h)) else np.array([-1.0])
    if np.all(eig > 0.0):
        cov = np.linalg.inv(h)
        sig_t0, sig_v = math.sqrt(cov[0, 0]), math.sqrt(cov[1, 1])
        corr = float(cov[0, 1] / (sig_t0 * sig_v))
    else:
        reason = "converged; the observed information is not positive definite, no sigma"
    margin_nll = drift_alt = math.nan
    if math.isfinite(sig_v):
        others = [(f, u) for f, u in found if abs(u - v) > float(separation_sigmas) * sig_v]
        margin_nll = math.inf
        if others:
            f_alt, drift_alt = min(others)
            margin_nll = f_alt - float(best.fun)
    return JointEstimate(
        True, reason, drift_mps=v, drift_sigma_mps=sig_v, t0_s=t0, t0_sigma_s=sig_t0, t_zc_s=t0 + kappa * v,
        width_s=w, noise_var=noise_var, corr_t0_v=corr, power=power, tau_c_fit_s=tau, nll=float(best.fun),
        at_bound=at_bound, nll_margin=margin_nll, drift_alt_mps=drift_alt,
    )
