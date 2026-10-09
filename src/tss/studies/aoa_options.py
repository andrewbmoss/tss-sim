"""Arrival-angle options for route three, graded against the drift requirement.

Layer:        study
WP:           WP11
Claims:       C3, C5
Status:       implemented
Requirements: SCI-15, SCI-27, INS-18

What route three must deliver
-----------------------------
The C3 studies fixed what an arrival angle buys.  Route one's tilt error is
``v_s sin(theta_sat)`` with ``theta_sat`` the along-track arrival angle at the
platform when the echo power peaks (:mod:`tss.studies.c3_verify`), so a measured
angle removes it with ``sigma_v = v_s sigma_theta`` whatever the plasma.  For a
field-aligned sheet only the rotation about B is unknown, and the drift error
from an error in it scales as ``v_s sin(psi)`` with ``psi`` the track heading
from magnetic north (:mod:`tss.studies.field_aligned_verify`).  The angle
requirement is therefore

    sigma_theta <= (L - sigma_r) / (v_s sin psi_eff),

with ``L`` the SCI-15 accuracy (50 m/s goal, 100 m/s threshold), ``sigma_r``
the rest of the route-one budget (about 12 m/s without fading, from
:mod:`tss.studies.route_one`) added linearly because a tilt bias does not
average down, and ``psi_eff = 90 deg`` for a compact facet.

The options
-----------
V6   co-located six-component vector sensor, three dipoles and three loops
     (AERO-VISTA geometry);
V3E  three orthogonal dipoles (RPI geometry);
V3H  three orthogonal loops;
WW   two crossed loops and a sense dipole (Watson-Watt);
IF   two-element interferometer on booms;
BI   a second spacecraft receiving the same echo (time difference of arrival);
P1   Poole's single-antenna phase gradient.

What is computed
----------------
1.  The Cramer-Rao bound on direction for each co-located sensor, from the
    deterministic-signal Fisher information with the polarization and the
    complex amplitude as nuisance parameters, averaged over direction and
    polarization.  Vacuum steering vectors at the platform.
2.  The bias from instrument errors, by Monte Carlo: channel gain and phase
    mismatch, and errors in the effective antenna directions (the effect
    Cassini calibrated in flight).
3.  Two plasma effects at the platform, from Appleton-Hartree: the angle
    between ray (Poynting) and wave normal, which biases any estimator built on
    ``E x H``, and the polarization axial ratio, which decides whether three
    dipoles alone can find a direction at all.
4.  The frequency below which a loop of given size is no longer sky-noise
    limited, against the Cane and Dulk background.
5.  The interferometer baseline each requirement needs.

Mass, volume and mission constraints are not computed; they are tabulated in
the summary from the heritage cited there and are estimates.
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from tss.constants import C_LIGHT
from tss.l1_propagation.magnetoionic import refractive_index_squared
from tss.l2_instrument.noise import cane_dulk_temperature_k
from tss.studies import StudyResult, register_study

__all__ = [
    "steering",
    "direction_crb_deg",
    "mean_crb_deg",
    "systematic_bias_deg",
    "ray_normal_offset_deg",
    "plasma_fields",
    "plasma_bias_deg",
    "plasma_model_bias_deg",
    "axial_ratio",
    "loop_sky_limited_floor_mhz",
    "interferometer_baseline_m",
    "plasma_crb_deg",
    "required_snr_per_sweep_db",
    "angle_requirement_deg",
    "CHANNELS",
]

V_S = 7600.0
SCI15 = {"goal": 50.0, "threshold": 100.0}
SIGMA_REST_MPS = 12.0

#: Which field components each co-located option measures: E then H, x y z.
CHANNELS = {
    "V6": np.array([1, 1, 1, 1, 1, 1], dtype=bool),
    "V3E": np.array([1, 1, 1, 0, 0, 0], dtype=bool),
    "V3H": np.array([0, 0, 0, 1, 1, 1], dtype=bool),
    "WW": np.array([0, 0, 1, 1, 1, 0], dtype=bool),  # loops x, y; sense dipole z
}


def steering(theta: float, phi: float, gamma: float, eta: float) -> np.ndarray:
    """Six-component response to a unit plane wave (Nehorai and Paldi 1994).

    ``theta, phi`` give the propagation direction ``k``; ``gamma, eta`` the
    polarization (auxiliary angle and phase difference).  E is in the plane
    perpendicular to k and H = k x E, in units of the free-space impedance, so
    every channel has the same noise scale.
    """
    st, ct, sp, cp = np.sin(theta), np.cos(theta), np.sin(phi), np.cos(phi)
    k = np.array([st * cp, st * sp, ct])
    e_th = np.array([ct * cp, ct * sp, -st])
    e_ph = np.array([-sp, cp, 0.0])
    pol = np.array([np.cos(gamma), np.sin(gamma) * np.exp(1j * eta)])
    e = e_ph * pol[0] + (-e_th) * pol[1]
    h = np.cross(k, e)
    return np.concatenate([e, h])


def _derivs(p: np.ndarray, mask: np.ndarray, step: float = 1e-6) -> np.ndarray:
    cols = []
    for i in range(4):
        d = np.zeros(4)
        d[i] = step
        cols.append((steering(*(p + d)) - steering(*(p - d)))[mask] / (2 * step))
    return np.array(cols).T


def direction_crb_deg(theta, phi, gamma, eta, snr_linear: float, n_snapshots: int,
                      mask: np.ndarray) -> float:
    """Root mean-square angular error bound, deg, for one co-located sensor.

    Deterministic CRB with an unknown complex amplitude per snapshot:
    ``FIM = (2 N SNR) Re(D^H P_perp D)``, ``P_perp`` the projector orthogonal
    to the steering vector.  ``SNR`` is the signal power for a unit E field
    over the noise power per channel.  The angular error combines the two
    direction parameters as ``sqrt(C_tt + sin^2(theta) C_pp)``.  Returns inf
    where the direction is not identifiable.
    """
    p = np.array([theta, phi, gamma, eta], dtype=float)
    a = steering(*p)[mask]
    d = _derivs(p, mask)
    proj = np.eye(a.size) - np.outer(a, a.conj()) / np.vdot(a, a).real
    fim = 2.0 * n_snapshots * snr_linear * np.real(d.conj().T @ proj @ d)
    # Direction block after removing the polarization nuisance (Schur
    # complement, pseudo-inverse so that a parameter that drops out, such as
    # eta for a linearly polarized wave, does not make the bound infinite).
    f_dd, f_dn, f_nn = fim[:2, :2], fim[:2, 2:], fim[2:, 2:]
    eff = f_dd - f_dn @ np.linalg.pinv(f_nn, rcond=1e-10) @ f_dn.T
    if np.linalg.cond(eff) > 1e10:
        return float("inf")
    crb = np.linalg.inv(eff)
    msae = crb[0, 0] + np.sin(theta) ** 2 * crb[1, 1]
    if not np.isfinite(msae) or msae <= 0:
        return float("inf")
    return float(np.degrees(np.sqrt(msae)))


def mean_crb_deg(option: str, snr_db: float, n_snapshots: int, *, rng, n_draws: int = 400,
                 near_nadir_deg: float | None = 10.0) -> dict:
    """Median and 90th percentile of the bound over random direction and polarization.

    Directions are drawn within ``near_nadir_deg`` of nadir (echoes from below),
    or over the sphere when ``None``.  Polarization is uniform on the Poincare
    sphere.
    """
    mask = CHANNELS[option]
    snr = 10.0 ** (snr_db / 10.0)
    vals = []
    for _ in range(n_draws):
        if near_nadir_deg is None:
            th = np.arccos(rng.uniform(-1, 1))
        else:
            th = np.pi - np.radians(near_nadir_deg) * np.sqrt(rng.uniform())
        ph = rng.uniform(0, 2 * np.pi)
        # Uniform on the Poincare sphere: 2*gamma has cos uniform.
        g = 0.5 * np.arccos(rng.uniform(-1, 1))
        e = rng.uniform(-np.pi, np.pi)
        vals.append(direction_crb_deg(th, ph, g, e, snr, n_snapshots, mask))
    v = np.array(vals)
    fin = v[np.isfinite(v)]
    return {
        "median_deg": float(np.median(fin)) if fin.size else float("inf"),
        "p90_deg": float(np.percentile(fin, 90)) if fin.size else float("inf"),
        "fraction_unidentifiable": float(1.0 - fin.size / v.size),
    }


def plasma_crb_deg(option: str, f_hz: float, dip_deg: float, *, f_n_sat_hz: float = 1.5e6,
                   f_h_hz: float = 1.0e6, snr_db: float, n_snapshots: int, mode: str = "O", rng,
                   n_draws: int = 200, near_nadir_deg: float = 10.0) -> dict:
    """Bound when the echo carries the magnetoionic polarization of its mode.

    The plasma wave's E is projected on the transverse basis to give the
    polarization angles, and the vacuum-manifold bound is evaluated there.  This
    shows where three dipoles or three loops lose the direction because the
    wave is close to linear.
    """
    mask = CHANNELS[option]
    snr = 10.0 ** (snr_db / 10.0)
    x, y = (f_n_sat_hz / f_hz) ** 2, f_h_hz / f_hz
    d = np.radians(dip_deg)
    b = np.array([np.cos(d), 0.0, -np.sin(d)])
    vals = []
    for _ in range(n_draws):
        th = np.pi - np.radians(near_nadir_deg) * np.sqrt(rng.uniform())
        ph = rng.uniform(0, 2 * np.pi)
        e = plasma_fields(th, ph, x, y, b, mode)[:3]
        st, ct, sp, cp = np.sin(th), np.cos(th), np.sin(ph), np.cos(ph)
        e_th = np.array([ct * cp, ct * sp, -st])
        e_ph = np.array([-sp, cp, 0.0])
        c1, c2 = np.vdot(e_ph, e), np.vdot(-e_th, e)
        nrm = np.hypot(abs(c1), abs(c2))
        gam = np.arctan2(abs(c2), abs(c1))
        eta = np.angle(c2) - np.angle(c1) if nrm > 0 else 0.0
        vals.append(direction_crb_deg(th, ph, gam, eta, snr, n_snapshots, mask))
    v = np.array(vals)
    fin = v[np.isfinite(v)]
    return {"median_deg": float(np.median(fin)) if fin.size else float("inf"),
            "p90_deg": float(np.percentile(fin, 90)) if fin.size else float("inf"),
            "fraction_unidentifiable": float(1 - fin.size / v.size)}


def required_snr_per_sweep_db(sigma_deg: float, n_snapshots: int) -> float:
    """Per-sweep SNR at which the six-component bound equals ``sigma_deg``.

    For the six-component sensor the bound is exactly ``1 / sqrt(N SNR)`` rad,
    independent of direction and polarization (checked numerically in the
    unit tests).
    """
    return float(10 * np.log10(1.0 / (n_snapshots * np.radians(sigma_deg) ** 2)))


def _estimate(y: np.ndarray, mask: np.ndarray, start: np.ndarray) -> np.ndarray:
    """Least-squares direction and polarization, amplitude projected out."""
    from scipy.optimize import minimize

    def cost(p):
        a = steering(*p)[mask]
        return -abs(np.vdot(a, y)) ** 2 / np.vdot(a, a).real

    res = minimize(cost, start, method="Nelder-Mead",
                   options={"xatol": 1e-9, "fatol": 1e-14, "maxiter": 6000})
    return res.x


def _angle_between(t1, p1, t2, p2) -> float:
    """Angle between two directions, folded onto one hemisphere.

    Three dipoles or three loops give the wave normal only up to sign, from
    the plane of the polarization ellipse.  The sounder knows its echoes come
    from below, so the sign is resolved and the fold is legitimate.  The
    six-component sensor never needs it.
    """
    k1 = np.array([np.sin(t1) * np.cos(p1), np.sin(t1) * np.sin(p1), np.cos(t1)])
    k2 = np.array([np.sin(t2) * np.cos(p2), np.sin(t2) * np.sin(p2), np.cos(t2)])
    return float(np.degrees(np.arccos(np.clip(abs(np.dot(k1, k2)), 0, 1))))


def systematic_bias_deg(option: str, *, gain_db: float = 0.0, phase_deg: float = 0.0,
                        direction_deg: float = 0.0, rng, n_draws: int = 60,
                        near_nadir_deg: float = 10.0) -> dict:
    """RMS direction bias, deg, from instrument errors, noise free.

    Each channel gets an independent gain error (dB rms), phase error (deg
    rms), and its antenna axis an independent pointing error (deg rms, the
    Cassini effect).  The estimator uses the ideal model.
    """
    mask = CHANNELS[option]
    idx = np.flatnonzero(mask)
    errs = []
    for _ in range(n_draws):
        th = np.pi - np.radians(near_nadir_deg) * np.sqrt(rng.uniform())
        ph = rng.uniform(0, 2 * np.pi)
        g = 0.5 * np.arccos(rng.uniform(-1, 1))
        e = rng.uniform(-np.pi, np.pi)
        full = steering(th, ph, g, e)
        y = np.empty(idx.size, dtype=complex)
        for j, ch in enumerate(idx):
            axis = np.zeros(3)
            axis[ch % 3] = 1.0
            if direction_deg > 0:
                pert = rng.normal(0, np.radians(direction_deg), 3)
                axis = axis + pert - np.dot(pert, axis) * axis
                axis /= np.linalg.norm(axis)
            field = full[:3] if ch < 3 else full[3:]
            amp = 10 ** (rng.normal(0, gain_db) / 20) * np.exp(1j * np.radians(rng.normal(0, phase_deg)))
            y[j] = amp * np.dot(axis, field)
        est = _estimate(y, mask, np.array([th, ph, g, e]))
        errs.append(_angle_between(th, ph, est[0], est[1]))
    errs = np.array(errs)
    return {"rms_deg": float(np.sqrt(np.mean(errs**2))), "p90_deg": float(np.percentile(errs, 90))}


def ray_normal_offset_deg(x: float, y: float, theta_b_rad: float, mode: str) -> float:
    """Angle between ray (Poynting) and wave normal, deg (Budden 1961 ch. 13).

    ``tan(alpha) = -(1/mu) d mu / d theta``.
    """
    h = 1e-5
    mu = np.sqrt(refractive_index_squared(x, y, theta_b_rad, mode))
    mp = np.sqrt(refractive_index_squared(x, y, theta_b_rad + h, mode))
    mm = np.sqrt(refractive_index_squared(x, y, theta_b_rad - h, mode))
    dmu = (mp - mm) / (2 * h)
    return float(np.degrees(np.arctan(abs(dmu / mu))))


def axial_ratio(x: float, y: float, theta_b_rad: float, mode: str) -> float:
    """Minor/major axis ratio of the wave polarization ellipse (0 linear, 1 circular).

    From the Appleton-Hartree polarization ``R = E_y / E_x`` (Budden 1961
    eq. 12.4 form): ``R = -i [Y_T^2 / (2 (1 - X)) -/+ sqrt(Y_T^4/(4 (1-X)^2) + Y_L^2)] / Y_L``,
    upper sign the ordinary wave.  For field-free or perpendicular propagation
    the wave is linear.
    """
    yt, yl = y * np.sin(theta_b_rad), y * np.cos(theta_b_rad)
    if abs(yl) < 1e-12:
        return 0.0
    a = yt**2 / (2 * (1 - x))
    s = np.sqrt(a**2 + yl**2)
    r = abs((a - s) / yl) if mode.upper() == "O" else abs((a + s) / yl)
    return float(min(r, 1.0 / r))


def plasma_fields(theta, phi, x: float, y: float, b_dir: np.ndarray, mode: str) -> np.ndarray:
    """E and H of a cold magnetoionic characteristic wave, six components.

    The wave normal is ``k`` (from ``theta, phi``) and the field direction
    ``b_dir``.  The dielectric tensor in the field frame is Stix's
    ``[[S, -iD, 0], [iD, S, 0], [0, 0, P]]`` with ``S = 1 - X/(1 - Y^2)``,
    ``D = -X Y/(1 - Y^2)`` and ``P = 1 - X``; the refractive index comes from
    Appleton-Hartree and E is the null vector of
    ``n^2 (k k^T - I) + epsilon``.  ``H = n k x E`` in free-space-impedance
    units, so H is exactly perpendicular to k while E is not.
    """
    k = np.array([np.sin(theta) * np.cos(phi), np.sin(theta) * np.sin(phi), np.cos(theta)])
    b = np.asarray(b_dir, dtype=float) / np.linalg.norm(b_dir)
    cos_tb = float(np.clip(np.dot(k, b), -1, 1))
    tb = np.arccos(cos_tb)
    n2 = float(refractive_index_squared(x, y, tb, mode))
    n = np.sqrt(n2)
    s_ = 1 - x / (1 - y * y)
    d_ = -x * y / (1 - y * y)
    p_ = 1 - x
    # Field frame: z along b.
    zb = b
    xb = np.cross(zb, [0.0, 0.0, 1.0]) if abs(zb[2]) < 0.99 else np.cross(zb, [1.0, 0.0, 0.0])
    xb /= np.linalg.norm(xb)
    yb = np.cross(zb, xb)
    rot = np.array([xb, yb, zb])  # rows: field-frame axes in lab coordinates
    eps_b = np.array([[s_, -1j * d_, 0], [1j * d_, s_, 0], [0, 0, p_]])
    eps = rot.T @ eps_b @ rot
    m = n2 * (np.outer(k, k) - np.eye(3)) + eps
    _, sv, vh = np.linalg.svd(m)
    e = vh[-1].conj()
    e /= np.linalg.norm(e)
    h = n * np.cross(k, e)
    return np.concatenate([e, h])


def plasma_bias_deg(option: str, x: float, y: float, b_dir: np.ndarray, mode: str, *, rng,
                    n_draws: int = 40, near_nadir_deg: float = 10.0) -> dict:
    """Direction bias, deg, when a vacuum-model estimator sees a magnetoionic wave."""
    mask = CHANNELS[option]
    errs = []
    for _ in range(n_draws):
        th = np.pi - np.radians(near_nadir_deg) * np.sqrt(rng.uniform())
        ph = rng.uniform(0, 2 * np.pi)
        yv = plasma_fields(th, ph, x, y, b_dir, mode)[mask]
        best = None
        for g0 in (0.3, 0.8, 1.2):
            for e0 in (-1.5, 0.0, 1.5):
                est = _estimate(yv, mask, np.array([th, ph, g0, e0]))
                a = steering(*est)[mask]
                c = -abs(np.vdot(a, yv)) ** 2 / np.vdot(a, a).real
                if best is None or c < best[0]:
                    best = (c, est)
        errs.append(_angle_between(th, ph, best[1][0], best[1][1]))
    errs = np.array(errs)
    return {"rms_deg": float(np.sqrt(np.mean(errs**2))), "max_deg": float(errs.max())}


def plasma_model_bias_deg(option: str, x: float, y: float, b_dir: np.ndarray, mode: str, *,
                          dx_frac: float = 0.0, dy_frac: float = 0.0, db_deg: float = 0.0,
                          rng, n_draws: int = 30, near_nadir_deg: float = 10.0) -> dict:
    """Bias, deg, when the estimator uses the magnetoionic model with wrong inputs.

    The data are the exact characteristic wave.  The estimator fits only the
    direction (the mode's polarization comes from the model) with ``X``, ``Y``
    and the field direction perturbed by the given fractions and angle, in a
    random sense on each draw.  This is the residual after the 0.5 deg
    vacuum-model bias has been modelled out.
    """
    from scipy.optimize import minimize

    mask = CHANNELS[option]
    errs = []
    for _ in range(n_draws):
        th = np.pi - np.radians(near_nadir_deg) * np.sqrt(rng.uniform())
        ph = rng.uniform(0, 2 * np.pi)
        yv = plasma_fields(th, ph, x, y, b_dir, mode)[mask]
        xm = x * (1 + dx_frac * rng.choice([-1, 1]))
        ym = y * (1 + dy_frac * rng.choice([-1, 1]))
        ax = np.cross(b_dir, rng.normal(size=3))
        ax /= np.linalg.norm(ax)
        ang = np.radians(db_deg)
        bm = b_dir * np.cos(ang) + np.cross(ax, b_dir) * np.sin(ang)

        def cost(p, xm=xm, ym=ym, bm=bm, yv=yv):
            a = plasma_fields(p[0], p[1], xm, ym, bm, mode)[mask]
            return -abs(np.vdot(a, yv)) ** 2 / np.vdot(a, a).real

        res = minimize(cost, [th, ph], method="Nelder-Mead", options={"xatol": 1e-10, "fatol": 1e-15})
        errs.append(_angle_between(th, ph, res.x[0], res.x[1]))
    errs = np.array(errs)
    return {"rms_deg": float(np.sqrt(np.mean(errs**2))), "max_deg": float(errs.max())}


def loop_sky_limited_floor_mhz(side_m: float, *, tube_diameter_m: float = 0.01,
                               threshold_k: float = 580.0) -> float:
    """Lowest frequency at which a square single-turn copper loop is sky-noise limited.

    Efficiency ``R_rad / (R_rad + R_loss)`` with ``R_rad = 31171 (A / lambda^2)^2``
    and skin-effect loss; sky limited when ``eta T_sky > threshold`` (the
    criterion of the project's band-edge analysis).  Cane and Dulk background.
    """
    area = side_m**2
    perim = 4 * side_m
    sigma = 5.8e7
    f = np.logspace(5, np.log10(3e7), 800)
    lam = C_LIGHT / f
    rrad = 31171.0 * (area / lam**2) ** 2
    rs = np.sqrt(np.pi * f * 4e-7 * np.pi / sigma)
    rloss = perim / (np.pi * tube_diameter_m) * rs
    eta = rrad / (rrad + rloss)
    ok = eta * cane_dulk_temperature_k(f / 1e6) > threshold_k
    return float(f[np.argmax(ok)] / 1e6) if ok.any() else float("inf")


def interferometer_baseline_m(f_hz: float, sigma_theta_deg: float, sigma_phase_deg: float) -> float:
    """``d = lambda sigma_phi / (2 pi sigma_theta)`` near broadside."""
    return float((C_LIGHT / f_hz) * np.radians(sigma_phase_deg) / (2 * np.pi * np.radians(sigma_theta_deg)))


def angle_requirement_deg(limit_mps: float, psi_eff_deg: float = 90.0,
                          sigma_rest_mps: float = SIGMA_REST_MPS) -> float:
    """Arrival-angle accuracy that keeps the drift inside ``limit``, linear sum."""
    room = limit_mps - sigma_rest_mps
    s = np.sin(np.radians(psi_eff_deg))
    return float(np.degrees(room / (V_S * s))) if s > 0 else float("inf")


@dataclass(frozen=True)
class _Rows:
    crb: list
    sys: list
    plasma: list
    loops: list
    ifm: list
    req: list
    pbias: list
    pcrb: list
    pmodel: list


def _compute(seed: int, n_draws: int) -> _Rows:
    rng = np.random.default_rng(seed)
    crb = []
    for opt in ("V6", "V3E", "V3H", "WW"):
        for snr_db in (10.0, 20.0, 30.0, 40.0):
            for n in (1, 45):
                r = mean_crb_deg(opt, snr_db, n, rng=rng, n_draws=n_draws)
                crb.append((opt, snr_db, n, r["median_deg"], r["p90_deg"], r["fraction_unidentifiable"]))
    sysr = []
    for opt in ("V6", "V3E", "V3H"):
        for kind, kw in (("phase_1deg", {"phase_deg": 1.0}), ("phase_2deg", {"phase_deg": 2.0}),
                         ("gain_0p2dB", {"gain_db": 0.2}), ("gain_0p5dB", {"gain_db": 0.5}),
                         ("axis_0p1deg", {"direction_deg": 0.1}), ("axis_0p2deg", {"direction_deg": 0.2}),
                         ("axis_0p5deg", {"direction_deg": 0.5}), ("axis_1deg", {"direction_deg": 1.0}),
                         ("axis_8deg", {"direction_deg": 8.0}),
                         ("cal_tight", {"direction_deg": 0.1, "gain_db": 0.05, "phase_deg": 1.0}),
                         ("cal_good", {"direction_deg": 0.2, "gain_db": 0.1, "phase_deg": 1.0}),
                         ("cal_cassini", {"direction_deg": 0.5, "gain_db": 0.2, "phase_deg": 2.0})):
            r = systematic_bias_deg(opt, rng=rng, n_draws=max(30, n_draws // 8), **kw)
            sysr.append((opt, kind, r["rms_deg"], r["p90_deg"]))
    plasma = []
    for f_mhz in (3.0, 5.0, 10.0, 20.0):
        x = (1.5 / f_mhz) ** 2
        for fh in (0.8, 1.4):
            y = fh / f_mhz
            for thb in (10.0, 45.0, 80.0):
                for mode in ("O", "X"):
                    plasma.append((f_mhz, fh, thb, mode,
                                   ray_normal_offset_deg(x, y, np.radians(thb), mode),
                                   axial_ratio(x, y, np.radians(thb), mode)))
    loops = [(s, loop_sky_limited_floor_mhz(s)) for s in (0.3, 0.5, 1.0, 2.0)]
    ifm = []
    for f in (1e6, 3e6, 5e6, 10e6, 20e6):
        for req in (0.29, 0.66, 1.1, 3.3):
            for ph in (1.0, 2.0):
                ifm.append((f / 1e6, req, ph, interferometer_baseline_m(f, req, ph)))
    req = []
    for psi in (90.0, 45.0, 15.0, 5.0, 2.0):
        g, t = angle_requirement_deg(SCI15["goal"], psi), angle_requirement_deg(SCI15["threshold"], psi)
        req.append((psi, g, t, required_snr_per_sweep_db(g, 45), required_snr_per_sweep_db(t, 45),
                    required_snr_per_sweep_db(g, 1), required_snr_per_sweep_db(t, 1)))
    pbias = []
    for f_mhz in (3.0, 5.0, 10.0):
        x, y = (1.5 / f_mhz) ** 2, 1.0 / f_mhz
        for dip in (10.0, 30.0, 60.0):
            d = np.radians(dip)
            b = np.array([np.cos(d), 0.0, -np.sin(d)])
            for opt in ("V6", "V3E", "V3H"):
                for mode in ("O", "X"):
                    r = plasma_bias_deg(opt, x, y, b, mode, rng=rng, n_draws=max(10, n_draws // 40))
                    pbias.append((f_mhz, dip, opt, mode, r["rms_deg"], r["max_deg"]))
    pcrb = []
    for dip in (0.0, 10.0, 30.0, 60.0):
        for opt in ("V6", "V3E", "V3H", "WW"):
            r = plasma_crb_deg(opt, 5e6, dip, snr_db=30.0, n_snapshots=45, rng=rng, n_draws=max(50, n_draws // 4))
            pcrb.append((dip, opt, r["median_deg"], r["p90_deg"], r["fraction_unidentifiable"]))
    pmodel = []
    x, y = (1.5 / 5) ** 2, 1.0 / 5
    b = np.array([np.cos(np.radians(30)), 0.0, -np.sin(np.radians(30))])
    for label, kw in (("X_10pct", {"dx_frac": 0.1}), ("X_1pct", {"dx_frac": 0.01}),
                      ("Y_5pct", {"dy_frac": 0.05}), ("Y_1pct", {"dy_frac": 0.01}),
                      ("B_0p5deg", {"db_deg": 0.5}), ("B_0p1deg", {"db_deg": 0.1})):
        for mode in ("O", "X"):
            r = plasma_model_bias_deg("V6", x, y, b, mode, rng=rng, n_draws=max(10, n_draws // 40), **kw)
            pmodel.append((label, mode, r["rms_deg"], r["max_deg"]))
    return _Rows(crb, sysr, plasma, loops, ifm, req, pbias, pcrb, pmodel)


def _write(path: Path, header: str, rows: list) -> None:
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(header + "\n")
        for r in rows:
            fh.write(",".join(f"{v:.6g}" if isinstance(v, float) else str(v) for v in r) + "\n")


CAL_LEVELS = {"cal_tight": "0.1 deg axes, 0.05 dB, 1 deg", "cal_good": "0.2 deg axes, 0.1 dB, 1 deg",
              "cal_cassini": "0.5 deg axes, 0.2 dB, 2 deg"}


def _figure(rows: _Rows, out_dir: Path) -> list[Path]:
    try:
        import matplotlib  # noqa: F401
    except ImportError:  # pragma: no cover
        return []
    from tss.figures import style

    style.use()
    fig, axes = style.figure(1, 3, figsize=(11.5, 3.6))

    # (a) Random error against per-sweep SNR at 45 sweeps, with requirement lines.
    ax = axes[0]
    for k, opt in enumerate(("V6", "V3E", "V3H", "WW")):
        r = sorted((x for x in rows.crb if x[0] == opt and x[2] == 45), key=lambda x: x[1])
        ax.semilogy([x[1] for x in r], [x[3] for x in r], "o-", ms=3, color=style.SERIES[k], label=f"{opt} median")
        if opt != "V6":
            ax.semilogy([x[1] for x in r], [x[4] for x in r], ":", lw=0.8, color=style.SERIES[k])
    for psi, g, _t, *_ in rows.req:
        if psi in (90.0, 15.0):
            ax.axhline(g, color=style.INK_MUTED, ls="--", lw=0.6)
            ax.text(40, g * 1.05, f"goal, psi {psi:g}", fontsize=6, ha="right", color=style.INK_MUTED)
    ax.set_xlabel("SNR per sweep (dB), 45 sweeps")
    ax.set_ylabel("arrival-angle error bound (deg)")
    ax.set_title("(a) Random error (dotted: 90th pct)", fontsize=9)
    ax.legend(fontsize=6, frameon=False)

    # (b) Systematic error per calibration level for V6, against requirements.
    ax = axes[1]
    labels = ["phase_2deg", "gain_0p2dB", "axis_0p2deg", "axis_0p5deg", "cal_tight", "cal_good", "cal_cassini"]
    vals = [next(x[2] for x in rows.sys if x[0] == "V6" and x[1] == lab) for lab in labels]
    ax.barh(range(len(labels)), vals, color=style.SERIES[0])
    ax.set_yticks(range(len(labels)))
    ax.set_yticklabels(labels, fontsize=7)
    for psi, g, t, *_ in rows.req:
        if psi in (90.0, 45.0):
            ax.axvline(g, color=style.INK_MUTED, ls="--", lw=0.6)
            ax.axvline(t, color=style.INK_MUTED, ls=":", lw=0.6)
    ax.set_xlabel("V6 direction bias, rms (deg)")
    ax.set_title("(b) Calibration (dashed goal, dotted threshold)", fontsize=9)

    # (c) Interferometer baseline against frequency for the requirement levels.
    ax = axes[2]
    for k, req in enumerate((0.29, 0.66, 1.1, 3.3)):
        r = sorted((x for x in rows.ifm if x[1] == req and x[2] == 2.0), key=lambda x: x[0])
        ax.loglog([x[0] for x in r], [x[3] for x in r], "o-", ms=3, color=style.SERIES[k], label=f"{req:g} deg")
    ax.axhspan(0.1, 3.0, color=style.SERIES[2], alpha=0.12, lw=0)
    ax.set_xlabel("frequency (MHz)")
    ax.set_ylabel("baseline for 2 deg channel phase (m)")
    ax.set_title("(c) Interferometer; shaded: fits a 6U", fontsize=9)
    ax.legend(fontsize=6, frameon=False)

    fig.tight_layout()
    return style.savefig(fig, out_dir / "aoa_options", formats=("png", "pdf"))


@register_study(
    "aoa_options",
    description="Arrival-angle options for route three, against the drift requirement",
    work_packages=["WP11"],
    claims=["C3", "C5"],
    status="implemented",
    requirements=["SCI-15", "SCI-27", "INS-18"],
)
def run(out_dir: Path, *, config: dict[str, Any], manifest: Any) -> StudyResult:
    opts = dict(config.get("aoa_options", {}))
    seed = int(opts.get("seed", 20261009))
    n_draws = int(opts.get("n_draws", 400))
    rows = _compute(seed, n_draws)
    outs = []
    for name, header, data in (
        ("aoa_crb.csv", "option,snr_db_per_sweep,n_sweeps,median_deg,p90_deg,frac_unidentifiable", rows.crb),
        ("aoa_systematic.csv", "option,error,rms_deg,p90_deg", rows.sys),
        ("aoa_plasma.csv", "f_mhz,f_h_mhz,theta_b_deg,mode,ray_normal_offset_deg,axial_ratio", rows.plasma),
        ("aoa_loops.csv", "loop_side_m,sky_limited_above_mhz", rows.loops),
        ("aoa_interferometer.csv", "f_mhz,required_deg,phase_error_deg,baseline_m", rows.ifm),
        ("aoa_requirement.csv",
         "psi_eff_deg,goal_deg,threshold_deg,snr_goal_n45_db,snr_thr_n45_db,snr_goal_n1_db,snr_thr_n1_db",
         rows.req),
        ("aoa_plasma_bias.csv", "f_mhz,dip_deg,option,mode,rms_deg,max_deg", rows.pbias),
        ("aoa_plasma_crb.csv", "dip_deg,option,median_deg,p90_deg,frac_unidentifiable", rows.pcrb),
        ("aoa_plasma_model_residual.csv", "error,mode,rms_deg,max_deg", rows.pmodel),
    ):
        p = out_dir / name
        _write(p, header, data)
        outs.append(p)
    outs += _figure(rows, out_dir)
    v6 = {x[1]: x[2] for x in rows.sys if x[0] == "V6"}
    metrics = {
        "seed": seed,
        "n_draws": n_draws,
        "requirement_deg": {f"{x[0]:g}": (round(x[1], 3), round(x[2], 3)) for x in rows.req},
        "snr_per_sweep_db_n45": {f"{x[0]:g}": (round(x[3], 1), round(x[4], 1)) for x in rows.req},
        "v6_calibration_rms_deg": {k: round(v6[k], 3) for k in CAL_LEVELS},
        "loop_sky_limited_above_mhz": {f"{s_:g} m": round(f_, 2) for s_, f_ in rows.loops},
    }
    summary = (
        "Arrival-angle options for route three.  Requirement (deg, goal/threshold) by effective "
        f"heading: {metrics['requirement_deg']}.  Per-sweep SNR the six-component sensor needs at "
        f"45 sweeps (dB): {metrics['snr_per_sweep_db_n45']}.  V6 calibration bias (deg rms): "
        f"{metrics['v6_calibration_rms_deg']}.  Loops sky-noise limited above (MHz): "
        f"{metrics['loop_sky_limited_above_mhz']}."
    )
    return StudyResult(outputs=outs, summary=summary, metrics=metrics)
