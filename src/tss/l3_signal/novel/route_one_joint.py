"""Layer: L3
WP: WP5
Claims: C3
Status: partial
Requirements: SCI-15

The joint route-one estimator (WP5 Phase 3 part 3, decisions (o), (r) and (s),
Andrew 2026-11-27).

WHY
---
The pinned estimator (:func:`tss.l3_signal.novel.route_one_estimator.estimate_route_one`)
splits the problem in two.  It takes the overflight time ``t_pk`` from a fit to
the power alone and the zero crossing ``t_zc`` from a power-weighted mean of
phase steps.  In part 2 its RMS drift error was 1.6 and 2.8 times the
Cramer-Rao bound at 20 and 10 dB per sweep (``docs/drift-tier.md`` section 9).
Both halves discard information: the power fit ignores the phase, and the
phase steps ignore the envelope's dependence on ``t0``, which is what makes
drift and overflight time separable.

THE MODEL
---------
The deterministic model of :mod:`tss.l3_signal.novel.drift_crb` (decision (r)):

    y_n = c s_n(t0, v, W) + w_n,
    s_n = sinc(x_n) exp(j 2 pi S [kappa v u_n - u_n^2 / 2]),
    u_n = t_n - t0,   x_n = 2 X u_n / W,
    S = 2 f v_s^2 / (c P'),   kappa = mu_sat P' / v_s^2,

with ``c = A exp(j phi0)`` complex and ``w_n`` white complex Gaussian noise.
The slope ``S`` is the ephemeris slope (decision (s)): the pinned free-slope
test stays the gate, so both estimators form a drift on the same echoes.
``t_zc - t0 = kappa v``, so ``v`` is the drift directly.

Under white Gaussian noise the maximum-likelihood estimate is the least-squares
fit of the complex samples.  The complex amplitude ``c`` enters linearly and is
solved for in closed form at each trial ``(t0, v, W)`` (variable projection), so
the nonlinear search has three parameters.  It starts from the pinned
estimate.  The covariance of ``(t0, v, W)`` is ``sigma^2 (J^T J)^-1`` from the
Jacobian of the projected residual, with ``sigma^2`` estimated from the
residual.  That gives ``drift_sigma_mps``.

LIMITS
------
The model is mismatched under fading, since a fading gain is not a constant
``c``.  It is then not efficient, and it is compared with the fading bound,
not the deterministic one.  The envelope is the ``sinc^2`` power envelope that
the pinned fit and the bound assume, so a facet whose envelope differs from it
is a model error shared by all three.
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

from tss.constants import C_LIGHT
from tss.l3_signal.novel.route_one_estimator import HALF_POWER_X, fit_overflight_time

__all__ = ["JointEstimate", "joint_model", "estimate_route_one_joint"]


@dataclass(frozen=True)
class JointEstimate:
    """What the joint fit returns for one echo.

    ``ok`` is false when the fit did not converge or its covariance is not
    usable; the other fields are then NaN and ``reason`` says why.
    """

    ok: bool
    reason: str
    drift_mps: float
    drift_sigma_mps: float
    t0_s: float
    t0_sigma_s: float
    t_zc_s: float
    width_s: float
    amplitude: complex
    noise_var: float
    corr_t0_v: float
    n_eval: int


def joint_model(t_s: np.ndarray, *, t0_s: float, drift_mps: float, width_s: float, slope_hz_per_s: float,
                kappa_s2_per_m: float) -> np.ndarray:
    """``s_n``, the unit-amplitude model of the module docstring."""
    u = np.asarray(t_s, dtype=float) - float(t0_s)
    x = 2.0 * HALF_POWER_X * u / float(width_s)
    phase = 2.0 * math.pi * float(slope_hz_per_s) * (float(kappa_s2_per_m) * float(drift_mps) * u - 0.5 * u * u)
    return np.sinc(x / math.pi) * np.exp(1j * phase)


def _project(y: np.ndarray, s: np.ndarray) -> tuple[complex, np.ndarray]:
    """The least-squares complex amplitude and the residual."""
    ss = float(np.real(np.vdot(s, s)))
    c = complex(np.vdot(s, y) / ss) if ss > 0.0 else 0.0j
    return c, y - c * s


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
    """The joint maximum-likelihood fit of ``(c, t0, v, W)`` (module docstring).

    ``t0_start_s`` and ``drift_start_mps`` come from the pinned estimator;
    ``width_start_s`` from the ``sinc^2`` power fit if not given.
    """
    from scipy.optimize import least_squares

    t = np.asarray(t_s, dtype=float)
    y = np.asarray(field, dtype=complex)
    s_rate = 2.0 * float(f_hz) * float(v_s_mps) ** 2 / (C_LIGHT * float(group_path_m))
    kappa = float(mu_sat) * float(group_path_m) / float(v_s_mps) ** 2
    dt = float(np.median(np.diff(t)))
    span = float(t[-1] - t[0])
    if width_start_s is None:
        _, width_start_s = fit_overflight_time(t, np.abs(y) ** 2)
    nan = math.nan

    def fail(reason: str, n: int = 0) -> JointEstimate:
        return JointEstimate(False, reason, nan, nan, nan, nan, nan, nan, complex(nan, nan), nan, nan, n)

    def resid(q: np.ndarray) -> np.ndarray:
        s = joint_model(t, t0_s=q[0], drift_mps=q[1], width_s=q[2], slope_hz_per_s=s_rate, kappa_s2_per_m=kappa)
        _, r = _project(y, s)
        return np.concatenate([r.real, r.imag])

    x0 = np.array([float(t0_start_s), float(drift_start_mps),
                   float(np.clip(width_start_s, 2.0 * dt, span))])
    if not np.all(np.isfinite(x0)):
        return fail("no finite start from the pinned estimator")
    lo = np.array([t[0], -np.inf, 2.0 * dt])
    hi = np.array([t[-1], np.inf, span])
    x0 = np.clip(x0, lo + 1e-12, hi - 1e-12)
    # Tolerances relative to the parameters' natural scales, not absolute
    # (the lesson of part 2's envelope fit, docs/drift-tier.md section 6).
    sol = least_squares(resid, x0=x0, bounds=(lo, hi), x_scale=np.array([dt, 1.0, dt]),
                        xtol=1e-10, ftol=1e-12, gtol=1e-10, method="trf")
    if not sol.success:
        return fail(f"fit did not converge: {sol.message}", int(sol.nfev))
    t0, v, w = (float(a) for a in sol.x)
    s = joint_model(t, t0_s=t0, drift_mps=v, width_s=w, slope_hz_per_s=s_rate, kappa_s2_per_m=kappa)
    c, r = _project(y, s)
    n_real = 2 * t.size
    n_par = 5  # c (two real), t0, v, W
    noise_var = float(np.sum(np.abs(r) ** 2) / (n_real - n_par))  # per real component
    jac = np.asarray(sol.jac, dtype=float)
    jtj = jac.T @ jac
    try:
        cov = noise_var * np.linalg.inv(jtj)
    except np.linalg.LinAlgError:
        return fail("singular information matrix", int(sol.nfev))
    if not np.all(np.isfinite(cov)) or np.any(np.diag(cov) <= 0.0):
        return fail("covariance not positive", int(sol.nfev))
    sig_t0, sig_v = math.sqrt(cov[0, 0]), math.sqrt(cov[1, 1])
    return JointEstimate(
        ok=True, reason="converged", drift_mps=v, drift_sigma_mps=sig_v, t0_s=t0, t0_sigma_s=sig_t0,
        t_zc_s=t0 + kappa * v, width_s=w, amplitude=c, noise_var=noise_var,
        corr_t0_v=float(cov[0, 1] / (sig_t0 * sig_v)), n_eval=int(sol.nfev),
    )
