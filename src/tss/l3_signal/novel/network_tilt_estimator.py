"""Layer: L3
WP: WP5, WP7
Claims: C3
Status: implemented (parameter level)
Requirements: SCI-15, SCI-23, SCI-27

Drift and tilt of one patch of the ionosphere from every link that reflects in
it, vertical and oblique, at one instant.

A networked constellation sounds a patch along many paths at once: the vertical
sounding of a satellite overhead, and the oblique links between pairs whose
reflection points fall in the patch.  Each link measures a group path and a
line-of-sight rate.  Each satellite also measures the profile, from its own
ionogram, so the level, still ionosphere's group path and rate of every link
are known (:mod:`tss.l1_propagation.stratified_link`).  What is left over is
the patch's departure from level and still,

    dh = s + gE x + gN y + w t,

which every link sees through the same four numbers: a height shift ``s``
that absorbs an error in the inverted profile, the two components of the
surface gradient and the vertical drift.  To first order a link's rate moves
by

    dr = c_w w + c_E gE + c_N gN + c_s s,

where ``c_w`` is about ``2 mu`` for a vertical echo and falls with obliquity,
and ``(c_E, c_N)`` is about ``-D`` times the velocity of the link's reflection
point over the ground (``D`` the link's drift coefficient).  A vertical
sounding therefore sees only the gradient along its own track, mixed with
drift; oblique links reflecting at the same place but with different
obliquity and with reflection points moving in different directions separate
all three.  The group paths add little on tilt but fix ``s``.  Nothing waits
for a second pass, so the evolution of the structure costs nothing.

:func:`estimate_patch` solves for the four by Gauss-Newton, with the exact
shifted-sphere forward model of :func:`~tss.l1_propagation.stratified_link.link_forward`
and the first-order Jacobian, and returns the covariance.  :func:`patch_crb`
gives the covariance from the Jacobian alone, for sizing.

An earlier version of this module used the virtual-mirror model; checked
against the Haselgrove ray tracer it read a level, still ionosphere as a
50 m/s drift and a 0.5 degree tilt, because a mirror on a sphere gets the
ray's horizontal wavenumber a couple of per cent wrong.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np

from tss.l1_propagation.bistatic_mirror import LinkGeometry
from tss.l1_propagation.stratified_link import (
    ChapmanProfile,
    LinkSensitivity,
    PatchState,
    link_forward,
    link_sensitivity,
    solve_ray,
)

__all__ = ["LinkMeasurement", "PatchEstimate", "estimate_patch", "patch_crb", "weighted_rows", "PARAMS",
           "PROFILE_PARAMS"]

PARAMS = ("w", "gE", "gN", "s")
#: Profile parameters that can be solved for with the patch: fractional changes
#: of every density and of every topside scale height.
PROFILE_PARAMS = ("scale", "top_scale")


def _perturbed(profile: ChapmanProfile, names, values) -> ChapmanProfile:
    kw = {}
    for nm, v in zip(names, values, strict=True):
        kw[nm] = getattr(profile, nm) * (1.0 + v)
    return replace(profile, **kw) if kw else profile


@dataclass(frozen=True)
class LinkMeasurement:
    """One link's measurement at one frequency, with its one-sigma errors."""

    link: LinkGeometry
    f_hz: float
    group_path_m: float
    los_rate_mps: float
    sigma_group_m: float = 300.0
    sigma_rate_mps: float = 10.0


@dataclass(frozen=True)
class PatchEstimate:
    w_mps: float
    g_east: float
    g_north: float
    shift_m: float
    covariance: np.ndarray        # of the solved-for parameters, in PARAMS order
    iterations: int
    residual_rms: float           # normalised residual, about 1 when the errors are as stated
    n_links: int
    profile_corrections: tuple = ()

    @property
    def tilt_deg(self) -> float:
        return float(np.degrees(np.arctan(np.hypot(self.g_east, self.g_north))))

    @property
    def sigma(self) -> tuple[float, float, float]:
        """One-sigma (w m/s, gE deg, gN deg)."""
        s = np.sqrt(np.diag(self.covariance))
        return float(s[0]), float(np.degrees(s[1])), float(np.degrees(s[2]))


def weighted_rows(sens: LinkSensitivity, sigma_group_m: float, sigma_rate_mps: float, *,
                  solve_shift: bool = True) -> np.ndarray:
    """The link's two rows of the whitened Jacobian in ``(w, gE, gN[, s])``."""
    n = 4 if solve_shift else 3
    return np.vstack([sens.d_group[:n] / sigma_group_m, sens.d_rate[:n] / sigma_rate_mps])


def patch_crb(profile: ChapmanProfile, lat: float, lon: float, links: list[tuple[LinkGeometry, float]], *,
              sigma_group_m: float = 300.0, sigma_rate_mps: float = 10.0,
              solve_shift: bool = True) -> np.ndarray | None:
    """Covariance bound of ``(w, gE, gN[, s])`` from the links that close, or ``None``."""
    rows = []
    for lk, f in links:
        s = link_sensitivity(profile, lat, lon, lk.a, lk.b, lk.va, lk.vb, f)
        if s is not None:
            rows.append(weighted_rows(s, sigma_group_m, sigma_rate_mps, solve_shift=solve_shift))
    if not rows:
        return None
    j = np.vstack(rows)
    return np.linalg.pinv(j.T @ j)


def estimate_patch(profile: ChapmanProfile, lat: float, lon: float, meas: list[LinkMeasurement], *,
                   solve_shift: bool = True, profile_params: tuple[str, ...] = (),
                   iterations: int = 12) -> PatchEstimate:
    """Gauss-Newton solve for ``(w, gE, gN[, s])`` from all the links on a patch.

    ``profile_params`` adds fractional corrections to the named profile
    parameters (see :data:`PROFILE_PARAMS`) to the solve, so that the
    profile is fitted to the same group paths and rates rather than taken as
    given.  The covariance then covers them too, after the patch parameters.
    """
    n = 4 if solve_shift else 3
    npar = n + len(profile_params)
    rows, refs, keep = [], [], []
    for m in meas:
        lk = m.link
        lev = solve_ray(profile, lk.a, lk.b, m.f_hz)
        if lev is None:
            continue
        kg = lev.k if lev.k > 0 else None
        s = link_sensitivity(profile, lat, lon, lk.a, lk.b, lk.va, lk.vb, m.f_hz, k_guess=kg)
        if s is None:
            continue
        row = np.zeros((2, npar))
        row[:, :n] = weighted_rows(s, m.sigma_group_m, m.sigma_rate_mps, solve_shift=solve_shift)
        for c, nm in enumerate(profile_params):
            out = []
            for sg in (+1.0, -1.0):
                pr = _perturbed(profile, (nm,), (sg * 1e-3,))
                fw = link_forward(pr, PatchState(lat, lon), lk.a, lk.b, lk.va, lk.vb, m.f_hz,
                                  r_ref=lev.turning_radius, k_guess=kg)
                out.append(None if fw is None else np.array(fw[:2]))
            if out[0] is None or out[1] is None:
                break
            d = (out[0] - out[1]) / 2e-3
            row[:, n + c] = d / np.array([m.sigma_group_m, m.sigma_rate_mps])
        else:
            rows.append(row)
            refs.append((lev.turning_radius, kg))
            keep.append(m)
    if not keep:
        raise ValueError("no link closes through the profile")
    jw = np.vstack(rows)
    y = np.ravel([[m.group_path_m / m.sigma_group_m, m.los_rate_mps / m.sigma_rate_mps] for m in keep])

    def predict(theta):
        full = np.zeros(4)
        full[:n] = theta[:n]
        st = PatchState(lat, lon, *full)
        pr = _perturbed(profile, profile_params, theta[n:])
        out = []
        for m, (r_ref, kg) in zip(keep, refs, strict=True):
            lk = m.link
            fw = link_forward(pr, st, lk.a, lk.b, lk.va, lk.vb, m.f_hz, r_ref=r_ref, k_guess=kg)
            if fw is None:
                raise ValueError("a link stopped closing during the solve")
            out.append([fw[0] / m.sigma_group_m, fw[1] / m.sigma_rate_mps])
        return np.ravel(out)

    theta = np.zeros(npar)
    tol = np.concatenate([[1e-3, 1e-7, 1e-7, 1e-1][:n], np.full(len(profile_params), 1e-7)])
    it = 0
    for it in range(1, int(iterations) + 1):  # noqa: B007 - the count is reported
        res = y - predict(theta)
        step, *_ = np.linalg.lstsq(jw, res, rcond=None)
        theta = theta + step
        if np.all(np.abs(step) < tol):
            break
    res = y - predict(theta)
    cov = np.linalg.pinv(jw.T @ jw)
    full = np.zeros(4)
    full[:n] = theta[:n]
    return PatchEstimate(w_mps=float(full[0]), g_east=float(full[1]), g_north=float(full[2]),
                         shift_m=float(full[3]), covariance=cov, iterations=it,
                         residual_rms=float(np.sqrt(np.mean(res**2))), n_links=len(keep),
                         profile_corrections=tuple(float(v) for v in theta[n:]))
