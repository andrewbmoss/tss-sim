"""Layer: L1
WP: WP5, WP7
Claims: C3
Status: implemented
Requirements: SCI-15, SCI-23, SCI-27

Tabulated link sensitivities for a constellation whose satellites share one
circular altitude, so that a network study can evaluate thousands of links
per snapshot.

Between two satellites at the same radius in a spherically stratified
ionosphere, the reflected ray depends only on the frequency and the central
angle ``Theta`` between them.  So does everything
:func:`tss.l1_propagation.stratified_link.link_sensitivity` computes, apart from
where the ray sits and how fast it moves.  For each frequency this module
tabulates, along the reflected (rising) branch of ``Theta(K)``:

* ``K``, the Bouguer invariant, so the level rate is ``-K dTheta/dt``;
* the group path ``P'``;
* the drift coefficient ``D = (1/2) int X' dP'`` (negative), the first-order
  phase-path change per metre that the profile rises;
* ``D_g = d(f D)/df``, the same for the group path;
* the derivatives of ``K`` and ``P'`` with respect to a fractional change in
  every density and in every topside scale height, at fixed ``Theta``.

A symmetric ray samples a linear height change ``g . x`` at its two legs
equally either side of the reflection point, so to first order the phase path
changes by ``D (s + g . x_c + w t)``, with ``x_c`` the gnomonic position of the
reflection point.  The sensitivities of one link in ``(w, gE, gN, s, scale,
top_scale)`` are then

    rate:   -D,  -(D xdot_c + Ddot x_c)_E,  -(D xdot_c + Ddot x_c)_N,  -Ddot,
            -(dK/dscale) Thetadot,  -(dK/dtop) Thetadot
    group:   0,   D_g x_c,E,   D_g x_c,N,   D_g,   dP'/dscale,   dP'/dtop

with ``Ddot = (dD/dTheta) Thetadot``.  A vertical sounding is the case
``Theta = 0``.  The table is checked against the ray-by-ray computation in the
unit tests.

Links within ``margin`` of either end of the reflected branch are not used:
there the two neighbouring rays merge (the skip distance and the maximum
range), the ray focuses, and neither the first-order model nor a single echo
can be relied on.
"""

from __future__ import annotations

from dataclasses import dataclass, replace

import numpy as np
from scipy.interpolate import PchipInterpolator

from tss.l1_propagation.bistatic_mirror import R_EARTH_M
from tss.l1_propagation.stratified_link import ChapmanProfile, _leg, _turning_radius, _x, _xr

__all__ = ["FrequencyTable", "LinkTable", "build_table", "PARAMS6", "rows_for_links", "rows_for_patch"]

PARAMS6 = ("w", "gE", "gN", "s", "scale", "top_scale")


def _ray(profile, k, f, r):
    """(Theta, P', D) of the symmetric ray with invariant ``k`` between two satellites at ``r``."""
    rt = _turning_radius(profile, k, f, r)
    if rt is None:
        return np.nan, np.nan, np.nan
    th, _, pg, rn, wn = _leg(profile, k, f, rt, r)
    d = 0.5 * float(np.sum(_xr(profile, rn, f) * wn))
    return 2.0 * th, 2.0 * pg, 2.0 * d


def _branch(profile, f, r, n_scan, n_fine):
    """K, Theta, P', D sampled along the reflected branch, or None."""
    mu2 = 1.0 - float(_x(profile, r, f))
    if mu2 <= 0.0:
        return None
    k_hi = np.sqrt(mu2) * r * (1.0 - 1e-9)
    ks = k_hi * np.sin(np.linspace(0.0, 0.5 * np.pi, n_scan))
    th = np.array([_ray(profile, k, f, r)[0] for k in ks])
    fin = np.isfinite(th)
    if not np.any(fin):
        return None
    i0 = int(np.argmax(fin))
    # First rising run of finite values at least four scan points long; a
    # shorter one is quadrature noise on a falling branch, which is what a
    # ray skimming the peak looks like near penetration.
    i = i0
    start = end = i0
    while i + 1 < n_scan:
        while i + 1 < n_scan and np.isfinite(th[i + 1]) and th[i + 1] <= th[i]:
            i += 1                 # Pedersen (falling) part, if any
        start = i
        while i + 1 < n_scan and np.isfinite(th[i + 1]) and th[i + 1] > th[i]:
            i += 1
        end = i
        if end - start >= 3:
            break
        if i + 1 < n_scan and not np.isfinite(th[i + 1]):
            i += 1
            while i + 1 < n_scan and not np.isfinite(th[i]):
                i += 1
        if end - start < 3 and i + 1 >= n_scan:
            return None
    if end - start < 3:
        return None
    k_lo = ks[max(start - 1, i0)] if start > i0 else ks[start]
    k_up = ks[min(end + 1, n_scan - 1)]
    # Fine sampling, dense at both ends, keeping the strictly rising part.
    u = 0.5 * (1.0 - np.cos(np.linspace(0.0, np.pi, n_fine)))
    kf = k_lo + (k_up - k_lo) * u
    vals = np.array([_ray(profile, k, f, r) for k in kf])
    ok = np.isfinite(vals[:, 0])
    kf, vals = kf[ok], vals[ok]
    j = int(np.argmax(vals[:, 0]))
    kf, vals = kf[:j + 1], vals[:j + 1]
    # Drop any leading falling part (Pedersen) so Theta increases throughout.
    if len(kf) > 2:
        dth = np.diff(vals[:, 0])
        neg = np.flatnonzero(dth <= 0.0)
        if neg.size:
            s0 = int(neg[-1]) + 1
            kf, vals = kf[s0:], vals[s0:]
    if len(kf) < 4:
        return None
    return kf, vals[:, 0], vals[:, 1], vals[:, 2]


def _interp(theta, y):
    return PchipInterpolator(theta, y, extrapolate=False)


@dataclass
class FrequencyTable:
    f_hz: float
    theta_lo: float
    theta_hi: float
    k: PchipInterpolator
    group: PchipInterpolator
    d: PchipInterpolator
    dd_dtheta: PchipInterpolator
    d_g: PchipInterpolator
    dk_dp: tuple
    dgroup_dp: tuple

    def usable(self, theta) -> np.ndarray:
        return (np.asarray(theta) >= self.theta_lo) & (np.asarray(theta) <= self.theta_hi)


def build_table(profile: ChapmanProfile, r_sat: float, f_hz: float, *, n_scan: int = 240, n_fine: int = 160,
                margin: float = 0.01, df_hz: float = 20e3, dp: float = 0.01,
                theta_cap: float = 0.6) -> FrequencyTable | None:
    """Tabulate one frequency.  ``margin`` is the fraction of the branch's span cut from each end.

    Branches that start beyond ``theta_cap`` (0.6 rad, 3800 km) are dropped:
    they are rays skimming the peak far beyond any usable signal, and near
    penetration the quadrature cannot tell their rising and falling parts apart.
    """
    base = _branch(profile, f_hz, r_sat, n_scan, n_fine)
    if base is None:
        return None
    kf, th, pg, d = base
    if th[0] > theta_cap:
        return None
    lo, hi = float(th[0]), float(th[-1])
    span = hi - lo
    # A branch that starts at the vertical needs no lower margin.
    th_lo = lo if lo < 1e-9 else lo + margin * span
    th_hi = hi - margin * span

    def table_at(prof, f):
        b = _branch(prof, f, r_sat, n_scan, n_fine)
        if b is None:
            return None
        return b

    pert = {}
    for key, prof, f in (("f", profile, f_hz + df_hz),
                         ("scale", replace(profile, scale=profile.scale * (1.0 + dp)), f_hz),
                         ("top_scale", replace(profile, top_scale=profile.top_scale * (1.0 + dp)), f_hz)):
        b = table_at(prof, f)
        if b is None:
            return None
        pert[key] = b
        # Restrict the usable range to where every table is defined.
        th_lo = max(th_lo, float(b[1][0]) if b[1][0] > 1e-9 else th_lo)
        th_hi = min(th_hi, float(b[1][-1]) - margin * span)
    if th_hi <= th_lo:
        return None
    k_i, g_i, d_i = _interp(th, kf), _interp(th, pg), _interp(th, d)
    kf2, th2, _, d2 = pert["f"]
    d_f = _interp(th2, d2)
    grid = np.linspace(th_lo, th_hi, 200) if th_lo > 0 else np.concatenate([[0.0], np.linspace(th_lo + 1e-6, th_hi, 199)])
    dg = (( (f_hz + df_hz) * d_f(grid) - f_hz * d_i(grid)) / df_hz)
    dk_dp, dgr_dp = [], []
    for key in ("scale", "top_scale"):
        kp, thp, gp, _ = pert[key]
        dk_dp.append(_interp(grid, (_interp(thp, kp)(grid) - k_i(grid)) / dp))
        dgr_dp.append(_interp(grid, (_interp(thp, gp)(grid) - g_i(grid)) / dp))
    return FrequencyTable(f_hz=f_hz, theta_lo=th_lo, theta_hi=th_hi, k=k_i, group=g_i, d=d_i,
                          dd_dtheta=d_i.derivative(), d_g=_interp(grid, dg), dk_dp=tuple(dk_dp),
                          dgroup_dp=tuple(dgr_dp))


@dataclass
class LinkTable:
    """All frequencies for one profile and one satellite radius."""

    profile: ChapmanProfile
    r_sat: float
    tables: dict

    @classmethod
    def build(cls, profile: ChapmanProfile, r_sat: float, freqs_hz, **kw) -> LinkTable:
        out = {}
        for f in freqs_hz:
            t = build_table(profile, r_sat, float(f), **kw)
            if t is not None:
                out[float(f)] = t
        return cls(profile=profile, r_sat=r_sat, tables=out)

    @property
    def theta_max(self) -> float:
        return max((t.theta_hi for t in self.tables.values()), default=0.0)

    def evaluate(self, f_hz: float, a, b, va, vb, *, dt_s: float = 0.5):
        """Vectorised over links.  ``a, b, va, vb`` are (n, 3) ECEF arrays.

        Returns a dict with ``ok`` (closes and usable), ``theta``, ``group``,
        ``rate`` (level), ``mid`` (unit direction of the reflection point),
        ``mid_dot`` (its rate) and the per-link quantities needed to form the
        sensitivity rows against any patch centre (:func:`rows_for_patch`).
        """
        t = self.tables.get(float(f_hz))
        a, b, va, vb = (np.atleast_2d(np.asarray(z, float)) for z in (a, b, va, vb))

        def geom(sign):
            pa, pb = a + sign * dt_s * va, b + sign * dt_s * vb
            ah = pa / np.linalg.norm(pa, axis=1, keepdims=True)
            bh = pb / np.linalg.norm(pb, axis=1, keepdims=True)
            th = np.arccos(np.clip(np.sum(ah * bh, axis=1), -1.0, 1.0))
            m = ah + bh
            m = m / np.linalg.norm(m, axis=1, keepdims=True)
            return th, m

        th, mid = geom(0.0)
        thp, mp = geom(+1.0)
        thm, mm = geom(-1.0)
        th_dot = (thp - thm) / (2.0 * dt_s)
        n = len(th)
        if t is None:
            return {"ok": np.zeros(n, bool)}
        ok = t.usable(th)
        thc = np.clip(th, t.theta_lo, t.theta_hi)
        k = t.k(thc)
        d = t.d(thc)
        out = {
            "ok": ok & np.isfinite(k) & np.isfinite(d),
            "theta": th, "theta_dot": th_dot,
            "group": t.group(thc), "rate": -k * th_dot,
            "d": d, "d_dot": t.dd_dtheta(thc) * th_dot, "d_g": t.d_g(thc),
            "dk_dp": np.stack([g(thc) for g in t.dk_dp], axis=1),
            "dgroup_dp": np.stack([g(thc) for g in t.dgroup_dp], axis=1),
            "mid": mid, "mid_p": mp, "mid_m": mm, "dt_s": dt_s,
        }
        return out


def rows_for_links(ev: dict, idx, lat0, lon0):
    """Unwhitened (group, rate) sensitivity rows in PARAMS6 for links ``idx``.

    ``lat0`` and ``lon0`` are the patch centres, one per link or one for all.
    Returns two (m, 6) arrays.
    """
    idx = np.asarray(idx)
    m = len(idx)
    la = np.radians(np.broadcast_to(np.asarray(lat0, float), (m,)))
    lo = np.radians(np.broadcast_to(np.asarray(lon0, float), (m,)))
    e = np.stack([-np.sin(lo), np.cos(lo), np.zeros(m)], axis=1)
    n = np.stack([-np.sin(la) * np.cos(lo), -np.sin(la) * np.sin(lo), np.cos(la)], axis=1)
    u = np.stack([np.cos(la) * np.cos(lo), np.cos(la) * np.sin(lo), np.sin(la)], axis=1)

    def gn(dirs):
        s_ = R_EARTH_M / np.sum(dirs * u, axis=1)
        return np.stack([np.sum(dirs * e, axis=1) * s_, np.sum(dirs * n, axis=1) * s_], axis=1)

    xc = gn(ev["mid"][idx])
    xdot = (gn(ev["mid_p"][idx]) - gn(ev["mid_m"][idx])) / (2.0 * ev["dt_s"])
    d = ev["d"][idx][:, None]
    ddot = ev["d_dot"][idx][:, None]
    dg = ev["d_g"][idx][:, None]
    th_dot = ev["theta_dot"][idx][:, None]
    rate = np.zeros((m, 6))
    grp = np.zeros((m, 6))
    rate[:, 0] = -d[:, 0]
    rate[:, 1:3] = -(d * xdot + ddot * xc)
    rate[:, 3] = -ddot[:, 0]
    rate[:, 4:6] = -ev["dk_dp"][idx] * th_dot
    grp[:, 1:3] = dg * xc
    grp[:, 3] = dg[:, 0]
    grp[:, 4:6] = ev["dgroup_dp"][idx]
    return grp, rate


def rows_for_patch(ev: dict, idx, lat0: float, lon0: float):
    """:func:`rows_for_links` against one patch centre."""
    return rows_for_links(ev, idx, lat0, lon0)
