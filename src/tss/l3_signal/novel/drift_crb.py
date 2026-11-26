"""Layer: L3
WP: WP5
Claims: C3
Status: partial
Requirements: SCI-15

Cramer-Rao bounds for route one (WP5 Phase 3, TV Section 8.21 item 2).

THE MODEL
---------
One complex value per sweep for one compact facet, as the drift tier forms it
(:mod:`tss.l3_signal.drift`, step 3), at times ``t_n``:

    y_n = a(t_n) exp(j phi(t_n)) g_n + w_n,      w_n ~ CN(0, 1),

with the facet's amplitude envelope ``a(t) = A sinc(x)``, ``x = 2 X (t - t0) / W``
(``X`` = :data:`HALF_POWER_X`, so ``a^2`` is the ``sinc^2`` power envelope of
half-power width ``W`` that route one fits), and the phase of route one's
Doppler law ``f_D = S (t_zc - t)``:

    phi(t) = phi0 + 2 pi S [kappa v (t - t0) - (t - t0)^2 / 2],
    S = 2 f v_s^2 / (c P'),   kappa = mu_sat P' / v_s^2,

so ``t_zc - t0 = kappa v`` and ``v`` is the drift.  The terms of the phase that
depend on ``t0`` and ``v`` but not on ``t`` are absorbed into ``phi0``.  The
SNR per sweep is ``A^2`` over the unit noise, the route-one budget's
definition (the per-sweep value at the envelope's peak).

``g_n`` is the fading gain: 1 for the deterministic bound, or a zero-mean
complex Gaussian process of unit power and correlation ``exp(-|tau| / tau_c)``
(``fading_kind="exp"``; ``"gauss"`` for ``exp(-tau^2 / tau_c^2)``) for the
Gaussian-covariance bound, the Rayleigh gain of the route-one budget
(:mod:`tss.l2_instrument.fading`, ``specular="rayleigh"``).

Parameters: ``A`` (or the power ``A^2`` with fading), ``phi0`` (deterministic
only: a common phase does not change a zero-mean process's covariance),
``t0``, ``v``, and ``W`` unless it is given as known.

WHAT IT SHOWS
-------------
``d phi / d v = 2 pi S kappa (t - t0)`` and ``d phi / d t0 = 2 pi S (t - t0)``
up to a constant: the phase alone cannot tell ``v`` from ``t0`` (the
degeneracy of ``claude/c3-t0-degeneracy.md``, TV Section 8.7), and the
envelope's dependence on ``t0`` is what makes the information matrix regular.
With ``t0`` known, ``v`` is a frequency offset ``(2 f mu_sat / c) v`` on a
tone of envelope ``a(t)``, and its bound is the tone's: Stein's FDOA accuracy
``sigma_f = 1 / (T_e sqrt(B T gamma))`` with ``T_e`` the RMS duration times
``2 pi``, as Moss (2014, eqs. 2.25 to 2.28) states it, in the limit of one
noisy channel (``gamma_2 -> inf``, so ``gamma = 2 gamma_1``); here
``B T gamma = 2 E``, ``E = sum a^2``.  :func:`stein_fdoa_sigma_hz` is that
closed form.

:func:`fisher_doppler_samples` is the Doppler-sample model of the degeneracy
note (``f_D`` measured with known variance), which gave 0.42 m/s with ``t0``
known and a singular matrix with ``t0`` free.

NOT DONE
--------
The oscillator, the range scale and the ephemeris terms of the budget; a
fading envelope other than Rayleigh (the Rician ``K`` of G1); the tier's
estimator is not efficient and is compared with the bound, not built from it.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from tss.constants import C_LIGHT
from tss.l3_signal.novel.route_one_estimator import HALF_POWER_X

__all__ = [
    "CrbResult",
    "fisher_route_one",
    "route_one_crb",
    "fisher_doppler_samples",
    "stein_fdoa_sigma_hz",
]


@dataclass(frozen=True)
class CrbResult:
    """An information matrix, its parameters, and the bounds read from it."""

    names: tuple[str, ...]
    fisher: np.ndarray
    sigma: dict[str, float] = field(default_factory=dict)
    condition: float = math.nan
    correlation_t0_v: float = math.nan


def _sinc_and_slope(x: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """``sin x / x`` and its derivative."""
    s = np.sinc(x / np.pi)
    with np.errstate(invalid="ignore", divide="ignore"):
        d = np.where(np.abs(x) < 1e-6, -x / 3.0, (x * np.cos(x) - np.sin(x)) / (x * x))
    return s, d


def _model(t: np.ndarray, *, amp: float, t0: float, v: float, width: float, s_rate: float,
           kappa: float) -> dict[str, np.ndarray]:
    """Envelope, phase and their derivatives (module docstring)."""
    u = t - t0
    x = 2.0 * HALF_POWER_X * u / width
    sn, ds = _sinc_and_slope(x)
    a = amp * sn
    return {
        "a": a, "phi": 2.0 * math.pi * s_rate * (kappa * v * u - 0.5 * u * u),
        "a_A": sn, "a_t0": -amp * ds * 2.0 * HALF_POWER_X / width, "a_W": -amp * ds * x / width,
        "phi_t0": 2.0 * math.pi * s_rate * u,  # up to a constant, which phi0 absorbs
        "phi_v": 2.0 * math.pi * s_rate * kappa * u,
    }


def fisher_route_one(
    t_s: np.ndarray,
    *,
    snr_db: float,
    width_s: float,
    f_hz: float,
    mu_sat: float,
    p_prime_m: float,
    v_s_mps: float = 7600.0,
    t0_s: float = 0.0,
    drift_mps: float = 0.0,
    tau_c_s: float | None = None,
    fading_kind: str = "exp",
    width_known: bool = False,
    t0_known: bool = False,
) -> CrbResult:
    """The information matrix of route one's per-sweep samples (module docstring).

    Deterministic (``tau_c_s`` None): ``I_ij = 2 sum (a_i a_j + a^2 phi_i phi_j)``
    over ``(A, phi0, t0, v[, W])``.  With fading: the zero-mean complex
    Gaussian form ``I_ij = tr(C^-1 C_i C^-1 C_j)``, ``C = D R D^H + I``,
    ``D = diag(a e^{j phi})``, over ``(P, t0, v[, W])``.  ``t0_known`` drops
    ``t0``; ``width_known`` drops ``W``.
    """
    t = np.asarray(t_s, dtype=float)
    s_rate = 2.0 * float(f_hz) * float(v_s_mps) ** 2 / (C_LIGHT * float(p_prime_m))
    kappa = float(mu_sat) * float(p_prime_m) / float(v_s_mps) ** 2
    amp = 10.0 ** (float(snr_db) / 20.0)
    m = _model(t, amp=amp, t0=float(t0_s), v=float(drift_mps), width=float(width_s), s_rate=s_rate,
               kappa=kappa)
    a, phi = m["a"], m["phi"]
    zeros = np.zeros_like(t)
    if tau_c_s is None:
        derivs = {"A": (m["a_A"], zeros), "phi0": (zeros, np.ones_like(t)),
                  "t0": (m["a_t0"], m["phi_t0"]), "v": (zeros, m["phi_v"]), "W": (m["a_W"], zeros)}
        names = [n for n in ("A", "phi0", "t0", "v", "W")
                 if not (n == "t0" and t0_known) and not (n == "W" and width_known)]
        fis = np.empty((len(names), len(names)))
        for i, ni in enumerate(names):
            ai, pi_ = derivs[ni]
            for j, nj in enumerate(names):
                aj, pj = derivs[nj]
                fis[i, j] = 2.0 * float(np.sum(ai * aj + a * a * pi_ * pj))
        return CrbResult(names=tuple(names), fisher=fis)

    if fading_kind not in ("exp", "gauss"):
        raise ValueError("fading_kind must be 'exp' or 'gauss'")
    lag = np.abs(t[:, None] - t[None, :]) / float(tau_c_s)
    r = np.exp(-lag) if fading_kind == "exp" else np.exp(-lag * lag)
    e = np.exp(1j * phi)
    dvec = a * e
    # Power P = A^2: d a / d P = a / (2 P).
    p = amp * amp
    dd = {
        "P": (a / (2.0 * p)) * e,
        "t0": (m["a_t0"] + 1j * a * m["phi_t0"]) * e,
        "v": (1j * a * m["phi_v"]) * e,
        "W": m["a_W"] * e,
    }
    names = [n for n in ("P", "t0", "v", "W") if not (n == "t0" and t0_known) and not (n == "W" and width_known)]
    cov = (dvec[:, None] * r * np.conj(dvec)[None, :]) + np.eye(t.size)
    cinv = np.linalg.inv(cov)
    dc = {}
    for n in names:
        x = dd[n][:, None] * r * np.conj(dvec)[None, :]
        dc[n] = cinv @ (x + x.conj().T)
    fis = np.empty((len(names), len(names)))
    for i, ni in enumerate(names):
        for j, nj in enumerate(names):
            fis[i, j] = float(np.real(np.sum(dc[ni] * dc[nj].T)))
    return CrbResult(names=tuple(names), fisher=fis)


def route_one_crb(result: CrbResult) -> CrbResult:
    """The bounds (square roots of the inverse's diagonal) and the conditioning.

    ``condition`` is the information matrix's condition number after scaling
    each parameter to unit diagonal; above about 1e12 the inverse is not
    meaningful and the bounds are reported as infinite.
    """
    f = np.asarray(result.fisher, dtype=float)
    dscale = np.sqrt(np.clip(np.diag(f), 1e-300, None))
    fn = f / np.outer(dscale, dscale)
    cond = float(np.linalg.cond(fn))
    sigma: dict[str, float] = {}
    corr = math.nan
    if not cond < 1e12:
        sigma = {n: math.inf for n in result.names}
    else:
        inv = np.linalg.inv(fn) / np.outer(dscale, dscale)
        for i, n in enumerate(result.names):
            sigma[n] = float(math.sqrt(max(inv[i, i], 0.0)))
        if "t0" in result.names and "v" in result.names:
            i, j = result.names.index("t0"), result.names.index("v")
            corr = float(inv[i, j] / math.sqrt(inv[i, i] * inv[j, j]))
    return CrbResult(names=result.names, fisher=f, sigma=sigma, condition=cond, correlation_t0_v=corr)


def fisher_doppler_samples(
    t_s: np.ndarray,
    *,
    sigma_fd_hz: float,
    f_hz: float,
    height_m: float,
    v_s_mps: float = 7600.0,
    t0_known: bool = False,
) -> CrbResult:
    """The Doppler-sample model of the degeneracy note.

    ``f_D(t) = (2 f / c) v - k (t - t0)``, ``k = 2 f v_s^2 / (c h)``, each
    sample with variance ``sigma_fd^2``.  Parameters ``(v, t0)`` or ``(v,)``.
    """
    t = np.asarray(t_s, dtype=float)
    k = 2.0 * float(f_hz) * float(v_s_mps) ** 2 / (C_LIGHT * float(height_m))
    g_v = np.full_like(t, 2.0 * float(f_hz) / C_LIGHT)
    g_t0 = np.full_like(t, k)
    cols = [g_v] if t0_known else [g_v, g_t0]
    j = np.stack(cols, axis=1)
    return CrbResult(names=("v",) if t0_known else ("v", "t0"),
                     fisher=(j.T @ j) / float(sigma_fd_hz) ** 2)


def stein_fdoa_sigma_hz(t_s: np.ndarray, amplitude: np.ndarray) -> float:
    """Stein's FDOA accuracy for one noisy channel (Moss 2014, eqs. 2.25 to 2.28).

    ``sigma_f = 1 / (T_e sqrt(B T gamma))`` with ``T_e = 2 pi t_rms`` (the RMS
    duration of ``|u|^2`` about its centroid) and, with the second channel
    noise free, ``gamma = 2 gamma_1``, so ``B T gamma = 2 E`` for samples of
    unit noise and energy ``E = sum |a|^2``.
    """
    t = np.asarray(t_s, dtype=float)
    p = np.abs(np.asarray(amplitude)) ** 2
    e = float(p.sum())
    tc = float(np.sum(p * t) / e)
    t_rms = math.sqrt(float(np.sum(p * (t - tc) ** 2) / e))
    return 1.0 / (2.0 * math.pi * t_rms * math.sqrt(2.0 * e))
