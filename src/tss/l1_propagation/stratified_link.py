"""Layer: L1
WP: WP5, WP7
Claims: C3
Status: implemented
Requirements: SCI-15, SCI-23, SCI-27

Vertical and oblique links through a spherically stratified ionosphere, with
the first-order effect of a tilt, a vertical drift and a uniform height shift.

Why not the virtual mirror alone
--------------------------------
The two satellites of an oblique link usually move at kilometres per second
relative to each other, so the line-of-sight rate of the link is dominated by
the rate at which their separation changes, multiplied by the ray's
horizontal wavenumber.  A drift of tens of metres per second rides on that.
The virtual mirror of :mod:`tss.l1_propagation.bistatic_mirror` gives the
horizontal wavenumber exactly only for a flat Earth; on a sphere, with the
mirror hundreds of kilometres below where the ray actually turns, it is off by
a couple of per cent, which is a hundred metres per second of rate.  Checked
against the Haselgrove ray tracer, a mirror-based estimator read a level,
still ionosphere as a 50 m/s drift and a 0.5 degree tilt.  The rate of the
level ionosphere therefore has to come from a ray through the profile itself,
which each sounder measures with its own ionogram.

The level ray
-------------
For a spherically stratified, field-free medium a ray conserves Bouguer's
invariant ``K = mu r sin(zeta)``.  Between two satellites at radii ``rA`` and
``rB``, ``K`` is fixed by the central angle between them,

    Theta(K) = sum over legs of  int_{r_t}^{r_s} K dr / (r sqrt(mu^2 r^2 - K^2)),

and the phase and group paths are

    P  = sum  int mu^2 r dr / sqrt(mu^2 r^2 - K^2),
    P' = sum  int      r dr / sqrt(mu^2 r^2 - K^2),

with the turning point ``r_t`` where ``mu r = K``.  The substitution
``r = r_t + tau^2`` removes the square-root singularity there, so plain
Gauss-Legendre quadrature is accurate.  ``Theta(K)`` rises from zero for the
vertical ray to a maximum and falls again for rays launched nearly
horizontally, which turn just below the satellite; the link is the
lower-``K`` (deeper) ray on the rising branch, the one a homing ray tracer
finds from a mirror seed.  Where ``Theta`` falls with ``K`` the solution is
either the Pedersen ray or, past the maximum, the nearly straight
line-of-sight ray between the satellites, whose perigee the same integrals
describe; neither is counted.  No reflected ray exists beyond the maximum of
the rising branch, which makes the closure test.

Because the phase path's gradient at each end is ``mu k``, the level rate
needs no differencing of the ray:

    -dP/dt = -( K dTheta/dt + sqrt(mu_A^2 - K^2/rA^2) drA/dt
                            + sqrt(mu_B^2 - K^2/rB^2) drB/dt ).

Tilt, drift and height shift
----------------------------
They enter by Fermat's principle at first order.  If every isodensity surface
is displaced radially by ``dh``, then ``X(r) -> X(r - dh)`` and the phase path
along the unperturbed ray changes by

    dP = -(1/2) int (dX/dr) dh dP'          (dX/dr < 0 on the topside),

which stays finite through the turning point because it is written against
group path.  With ``dh = s + gE x + gN y + w t`` in the gnomonic tangent
coordinates of :class:`~tss.l1_propagation.bistatic_mirror.Patch`, ``dP`` is
linear in ``(w, gE, gN, s)``.  Its rate as the satellites move is the tilt's
share of the Doppler; its frequency derivative ``d(f dP)/df`` is the group
path's share.  :func:`link_sensitivity` returns both.

The sign of ``dP`` deserves a word: raising every surface by ``dh`` moves the
reflection up, toward a satellite above it, so the path shortens.  In the
quadrature ``dX/dr`` is the derivative of ``X`` with height, negative above
the peak, and ``dP = (1/2) int X'(r) dh dP'`` comes out negative as it should;
the code uses that form.

Beyond first order: the shifted sphere
--------------------------------------
First order is good to a fraction of a metre per second in rate for tilts up
to a degree, but a tilt also shortens every echo by a second-order amount,
``h alpha^2 / 2`` for a vertical echo, which at a degree of tilt is a few
hundred metres of group path: as large as the measurement error.  The forward
model the estimator iterates on, :func:`link_forward`, therefore traces the
link exactly through the profile with its centre moved by

    d = R (gE e + gN n) + (s + w t + R^2 |g|^2 / (2 r_ref)) u,

where ``(e, n, u)`` are east, north and up at the patch centre.  Near the
patch the surfaces of a sphere whose centre has moved sideways by ``R g``
rise by ``g . x`` to first order in ``x / R``, exactly the gnomonic tilt, and
sink everywhere by ``|d|^2 / 2r``; the last term of ``d`` puts that back at
the link's own reflection radius ``r_ref``.  What remains differs from a
tilted gnomonic plane by terms cubic in ``x / R``, a few tens of metres at the
far ends of a 1500 km link and nothing at the reflection point.  The ray
through the shifted sphere is a level ray between the shifted satellite
positions, so it costs one call to :func:`solve_ray`.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

from tss.l1_propagation.bistatic_mirror import R_EARTH_M, Patch

__all__ = ["ChapmanProfile", "StratifiedRay", "solve_ray", "link_sensitivity", "LinkSensitivity",
           "QUIET_DAYTIME", "max_central_angle", "PatchState", "link_forward"]

K806 = 80.6163849431291
_NODES, _WEIGHTS = np.polynomial.legendre.leggauss(64)
_SUB_X = 0.5 * (_NODES + 1.0)          # sub-quadrature nodes on [0, 1]
_SUB_W = 0.5 * _WEIGHTS


@dataclass(frozen=True)
class ChapmanProfile:
    """Sum of Chapman layers, ``N = Nm exp(c (1 - z - exp(-z)))``, ``z = (h - hm)/H``.

    ``layers`` holds ``(nm, hm, H_bottom, H_top, c)``.  ``scale`` multiplies
    every density and ``top_scale`` every topside scale height, for
    profile-error studies.
    """

    layers: tuple = ()
    scale: float = 1.0
    top_scale: float = 1.0

    def _terms(self, h):
        for nm, hm, hb, ht, c in self.layers:
            hh = np.where(h >= hm, ht * self.top_scale, hb)
            z = np.clip((h - hm) / hh, -30.0, 60.0)
            n = nm * np.exp(c * (1.0 - z - np.exp(-z)))
            yield n, c, z, hh

    def density(self, alt_m):
        h = np.asarray(alt_m, dtype=float)
        out = np.zeros(h.shape)
        for n, *_ in self._terms(h):
            out = out + n
        return out * self.scale

    def d_density(self, alt_m):
        """Analytic derivative with height, m^-4."""
        h = np.asarray(alt_m, dtype=float)
        out = np.zeros(h.shape)
        for n, c, z, hh in self._terms(h):
            out = out + n * c * (np.exp(-z) - 1.0) / hh
        return out * self.scale


#: The quiet daytime profile of the Chartier reproduction (foF2 about 10.2 MHz, hmF2 300 km).
QUIET_DAYTIME = ChapmanProfile(layers=(
    (1.3e12, 300e3, 45e3, 90e3, 0.5),
    (2.2e11, 180e3, 30e3, 30e3, 0.5),
    (1.4e11, 110e3, 10e3, 10e3, 0.5),
))


def _x(profile, r, f):
    return K806 * profile.density(np.asarray(r) - R_EARTH_M) / f**2


def _xr(profile, r, f):
    return K806 * profile.d_density(np.asarray(r) - R_EARTH_M) / f**2


_SCAN = 1.0 - np.linspace(0.0, 1.0, 1200) ** 1.5      # fractions of the scan depth, top first


def _turning_radius(profile, k, f, r_top):
    """Largest r below r_top where mu r = K, or None if the ray reaches 60 km."""
    def g(r):
        return (1.0 - _x(profile, r, f)) * r * r - k * k
    r_floor = R_EARTH_M + 60e3
    rs = r_floor + (r_top - r_floor) * _SCAN
    gv = g(rs)
    idx = np.flatnonzero(gv <= 0.0)
    if idx.size == 0:
        return None
    i = int(idx[0])
    if i == 0:
        return float(r_top)
    return float(brentq(g, rs[i], rs[i - 1], xtol=1e-7, rtol=1e-15))


def _root(profile, k, f, r_t, tau):
    """sqrt(mu^2 r^2 - K^2) at r = r_t + tau^2, regular through the turning point.

    Near the turning point the difference cancels to the precision of ``r_t``,
    so below a millimetre it is replaced by its linear expansion
    ``tau sqrt(G'(r_t))``.
    """
    r = r_t + tau * tau
    g = (1.0 - _x(profile, r, f)) * r * r - k * k
    gp = 2.0 * r_t * (1.0 - _x(profile, r_t, f)) - r_t * r_t * _xr(profile, r_t, f)
    small = tau * tau < 1e-3
    return np.where(small, tau * np.sqrt(max(float(gp), 1e-300)), np.sqrt(np.maximum(g, 1e-300)))


def _leg(profile, k, f, r_t, r_s, *, cumulative: bool = False):
    """Central angle, phase path, group path, node radii, node dP' weights [, cumulative angle]."""
    tau_max = np.sqrt(max(r_s - r_t, 0.0))
    tau = tau_max * _SUB_X
    wt = tau_max * _SUB_W
    r = r_t + tau * tau
    mu2 = 1.0 - _x(profile, r, f)
    root = _root(profile, k, f, r_t, tau)
    jac = 2.0 * tau
    dtheta = k / (r * root) * jac
    dp = mu2 * r / root * jac
    dpg = r / root * jac
    out = (float(np.sum(wt * dtheta)), float(np.sum(wt * dp)), float(np.sum(wt * dpg)), r, wt * dpg)
    if not cumulative:
        return out
    if k == 0.0:
        return out + (np.zeros(r.size),)
    # Angle from the turning point to each node: the same quadrature on [0, tau_i].
    tt = tau[:, None] * _SUB_X[None, :]
    rr = r_t + tt * tt
    rt_ = _root(profile, k, f, r_t, tt)
    cum = np.sum(_SUB_W[None, :] * (k / (rr * rt_) * 2.0 * tt), axis=1) * tau
    return out + (cum,)


@dataclass(frozen=True)
class StratifiedRay:
    k: float
    turning_radius: float
    theta: float
    phase_path_m: float
    group_path_m: float
    r_nodes: np.ndarray          # radii of the quadrature nodes, leg A then leg B
    w_nodes: np.ndarray          # dP' weight of each node
    phi_nodes: np.ndarray        # central angle of each node from A along the great circle
    theta_max: float | None = None


def _theta(profile, k, f, ra, rb):
    rt = _turning_radius(profile, k, f, min(ra, rb))
    if rt is None:
        return float("nan")
    return _leg(profile, k, f, rt, ra)[0] + _leg(profile, k, f, rt, rb)[0]


def _k_grid(k_hi, n=48):
    return k_hi * np.sin(np.linspace(0.0, 0.5 * np.pi, n))


def max_central_angle(profile: ChapmanProfile, r_sat: float, f_hz: float, *, n: int = 48) -> float:
    """Largest central angle a reflected link between two satellites at ``r_sat`` can span."""
    mu2 = 1.0 - float(_x(profile, r_sat, f_hz))
    if mu2 <= 0.0:
        return 0.0
    ks = _k_grid(np.sqrt(mu2) * r_sat * (1.0 - 1e-9), n)
    th = np.array([_theta(profile, k, f_hz, r_sat, r_sat) for k in ks])
    if not np.any(np.isfinite(th)):
        return 0.0
    i = int(np.nanargmax(th))
    lo, hi = ks[max(i - 1, 0)], ks[min(i + 1, n - 1)]
    # Golden-section refinement of the maximum.
    gr = 0.5 * (np.sqrt(5.0) - 1.0)
    c, d = hi - gr * (hi - lo), lo + gr * (hi - lo)
    fc, fd = _theta(profile, c, f_hz, r_sat, r_sat), _theta(profile, d, f_hz, r_sat, r_sat)
    for _ in range(30):
        if fc > fd:
            hi, d, fd = d, c, fc
            c = hi - gr * (hi - lo)
            fc = _theta(profile, c, f_hz, r_sat, r_sat)
        else:
            lo, c, fc = c, d, fd
            d = lo + gr * (hi - lo)
            fd = _theta(profile, d, f_hz, r_sat, r_sat)
    return float(max(np.nanmax(th), fc, fd))


def solve_ray(profile: ChapmanProfile, a: np.ndarray, b: np.ndarray, f_hz: float, *,
              k_guess: float | None = None) -> StratifiedRay | None:
    """The deeper reflected ray from ``a`` to ``b``, or ``None`` if there is none.

    ``k_guess``, from a ray at a nearby geometry or frequency, skips the scan.
    """
    a, b = np.asarray(a, float), np.asarray(b, float)
    ra, rb = float(np.linalg.norm(a)), float(np.linalg.norm(b))
    theta = float(np.arccos(np.clip(np.dot(a, b) / (ra * rb), -1.0, 1.0)))
    r_top = min(ra, rb)
    mu_top2 = 1.0 - float(_x(profile, r_top, f_hz))
    if mu_top2 <= 0.0:
        return None
    k_hi = np.sqrt(mu_top2) * r_top * (1.0 - 1e-9)
    theta_max = None
    if theta < 1e-7:
        k = 0.0
    else:
        def resid(kk):
            return _theta(profile, kk, f_hz, ra, rb) - theta
        bracket = None
        if k_guess is not None and 0.0 < k_guess < k_hi:
            for d in (1e-4, 1e-3, 1e-2):
                lo, hi = k_guess * (1.0 - d), min(k_guess * (1.0 + d), k_hi)
                flo, fhi = resid(lo), resid(hi)
                if np.isfinite(flo) and np.isfinite(fhi) and flo <= 0.0 <= fhi:
                    bracket = (lo, hi)
                    break
        if bracket is None:
            ks = _k_grid(k_hi)
            res = np.array([resid(kk) for kk in ks])
            theta_max = float(np.nanmax(res) + theta) if np.any(np.isfinite(res)) else 0.0
            # Only the rising branch is a reflection: where Theta falls with K
            # the ray is either the Pedersen ray, or, beyond the maximum, the
            # nearly straight line-of-sight ray between the satellites, whose
            # perigee the same integrals also describe.
            ok = (np.isfinite(res[:-1]) & np.isfinite(res[1:]) & (res[:-1] <= 0.0) & (res[1:] >= 0.0))
            hits = np.flatnonzero(ok)
            if hits.size == 0:
                return None
            i = int(hits[0])
            bracket = (ks[i], ks[i + 1])
        k = float(brentq(resid, *bracket, xtol=1e-10 * k_hi, rtol=1e-14))
    rt = _turning_radius(profile, k, f_hz, r_top)
    if rt is None:
        return None
    ta, pa, ga, r_a, w_a, c_a = _leg(profile, k, f_hz, rt, ra, cumulative=True)
    tb, pb, gb, r_b, w_b, c_b = _leg(profile, k, f_hz, rt, rb, cumulative=True)
    return StratifiedRay(k=k, turning_radius=rt, theta=ta + tb, phase_path_m=pa + pb, group_path_m=ga + gb,
                         r_nodes=np.concatenate([r_a, r_b]), w_nodes=np.concatenate([w_a, w_b]),
                         phi_nodes=np.concatenate([ta - c_a, ta + c_b]), theta_max=theta_max)


def _node_dirs(ray: StratifiedRay, a, b):
    ah = a / np.linalg.norm(a)
    bh = b / np.linalg.norm(b)
    perp = bh - np.dot(ah, bh) * ah
    n = np.linalg.norm(perp)
    if n < 1e-9:
        return np.tile(ah, (ray.r_nodes.size, 1))
    perp = perp / n
    return np.outer(np.cos(ray.phi_nodes), ah) + np.outer(np.sin(ray.phi_nodes), perp)


def _moments(profile, ray: StratifiedRay, a, b, f, patch: Patch):
    """First-order phase-path change per unit (w t, gE, gN, s): (D, D x, D y, D) with D = (1/2) int X' dP'."""
    xr = _xr(profile, ray.r_nodes, f)
    wgt = 0.5 * xr * ray.w_nodes
    dirs = _node_dirs(ray, a, b)
    e, n, u = patch.basis()
    s = R_EARTH_M / (dirs @ u)
    x = (dirs @ e) * s
    y = (dirs @ n) * s
    d = float(np.sum(wgt))
    return np.array([d, float(np.sum(wgt * x)), float(np.sum(wgt * y)), d])


@dataclass(frozen=True)
class LinkSensitivity:
    """Level-ionosphere group path and rate of a link, and their sensitivities.

    ``rate_mps`` is the line-of-sight rate (minus the phase-path rate) of the
    level, still ionosphere.  ``d_rate`` and ``d_group`` are derivatives with
    respect to ``(w, gE, gN, s)``: vertical drift (m/s), surface gradient
    east and north (m/m) and a uniform height shift (m).
    """

    group_path_m: float
    rate_mps: float
    d_rate: np.ndarray
    d_group: np.ndarray
    turning_alt_m: float
    k: float
    theta: float
    sin_zenith_a: float
    reflection_dir: np.ndarray


def _level_rate(profile, ray, a, b, va, vb, f, dt_s):
    ra, rb = np.linalg.norm(a), np.linalg.norm(b)

    def ang(sign):
        pa, pb = a + sign * dt_s * va, b + sign * dt_s * vb
        return float(np.arccos(np.clip(np.dot(pa, pb) / (np.linalg.norm(pa) * np.linalg.norm(pb)), -1.0, 1.0)))

    th_dot = (ang(+1) - ang(-1)) / (2.0 * dt_s)
    mu2a = 1.0 - float(_x(profile, ra, f))
    mu2b = 1.0 - float(_x(profile, rb, f))
    ca = np.sqrt(max(mu2a - (ray.k / ra) ** 2, 0.0))
    cb = np.sqrt(max(mu2b - (ray.k / rb) ** 2, 0.0))
    ra_dot = float(np.dot(a, va)) / ra
    rb_dot = float(np.dot(b, vb)) / rb
    return -(ray.k * th_dot + ca * ra_dot + cb * rb_dot)


def link_sensitivity(profile: ChapmanProfile, lat0: float, lon0: float, a, b, va, vb, f_hz: float, *,
                     dt_s: float = 0.5, df_hz: float = 5e3, k_guess: float | None = None) -> LinkSensitivity | None:
    """Group path, rate and their sensitivities for one link, or ``None`` if it does not close."""
    a, b = np.asarray(a, float), np.asarray(b, float)
    va, vb = np.asarray(va, float), np.asarray(vb, float)
    ray = solve_ray(profile, a, b, f_hz, k_guess=k_guess)
    if ray is None:
        return None
    patch = Patch(lat0, lon0)
    rate = _level_rate(profile, ray, a, b, va, vb, f_hz, dt_s)
    m0 = _moments(profile, ray, a, b, f_hz, patch)
    # Rate of the first-order path change as the satellites move.
    mt = []
    for sgn in (+1, -1):
        ap, bp = a + sgn * dt_s * va, b + sgn * dt_s * vb
        rr = solve_ray(profile, ap, bp, f_hz, k_guess=ray.k if ray.k > 0 else None)
        if rr is None:
            return None
        mt.append(_moments(profile, rr, ap, bp, f_hz, patch))
    m_dot = (mt[0] - mt[1]) / (2.0 * dt_s)
    # dP = m . (w t, gE, gN, s): the rate of w t is w itself.
    d_rate = -np.array([m0[0], m_dot[1], m_dot[2], m_dot[3]])
    # Group path: d(f dP)/df, at t = 0 so drift does not enter.
    mf = []
    for fs in (f_hz + df_hz, f_hz - df_hz):
        rr = solve_ray(profile, a, b, fs, k_guess=ray.k if ray.k > 0 else None)
        if rr is None:
            return None
        mf.append(fs * _moments(profile, rr, a, b, fs, patch))
    d_group = (mf[0] - mf[1]) / (2.0 * df_hz)
    d_group[0] = 0.0
    mu_a = np.sqrt(max(1.0 - float(_x(profile, np.linalg.norm(a), f_hz)), 0.0))
    mid = _node_dirs(ray, a, b)[np.argmin(ray.r_nodes)]
    return LinkSensitivity(group_path_m=ray.group_path_m, rate_mps=rate, d_rate=d_rate, d_group=d_group,
                           turning_alt_m=ray.turning_radius - R_EARTH_M, k=ray.k, theta=ray.theta,
                           sin_zenith_a=ray.k / (mu_a * np.linalg.norm(a)) if mu_a > 0 else float("nan"),
                           reflection_dir=mid)


@dataclass(frozen=True)
class PatchState:
    """Drift, tilt and height shift of a patch: ``dh = s + gE x + gN y + w t``."""

    lat_deg: float
    lon_deg: float
    w_mps: float = 0.0
    g_east: float = 0.0
    g_north: float = 0.0
    shift_m: float = 0.0

    def vector(self) -> np.ndarray:
        return np.array([self.w_mps, self.g_east, self.g_north, self.shift_m])

    def centre_offset(self, r_ref: float, t: float = 0.0) -> tuple[np.ndarray, np.ndarray]:
        """Offset of the profile's centre (ECEF, m) and its rate (m/s)."""
        e, n, u = Patch(self.lat_deg, self.lon_deg).basis()
        dh = R_EARTH_M * (self.g_east * e + self.g_north * n)
        up = self.shift_m + self.w_mps * t + float(np.dot(dh, dh)) / (2.0 * r_ref)
        return dh + up * u, self.w_mps * u


def link_forward(profile: ChapmanProfile, state: PatchState, a, b, va, vb, f_hz: float, *,
                 r_ref: float | None = None, k_guess: float | None = None, dt_s: float = 0.5):
    """Group path and line-of-sight rate through the tilted, drifting profile, or ``None``.

    Exact for the shifted sphere described in the module docstring.  ``r_ref``
    defaults to the turning radius of the level ray.  Returns
    ``(group_path_m, rate_mps, ray)``.
    """
    a, b = np.asarray(a, float), np.asarray(b, float)
    va, vb = np.asarray(va, float), np.asarray(vb, float)
    if r_ref is None:
        level = solve_ray(profile, a, b, f_hz, k_guess=k_guess)
        if level is None:
            return None
        r_ref, k_guess = level.turning_radius, (level.k if level.k > 0 else None)
    d, d_dot = state.centre_offset(r_ref)
    ap, bp = a - d, b - d
    ray = solve_ray(profile, ap, bp, f_hz, k_guess=k_guess)
    if ray is None:
        return None
    rate = _level_rate(profile, ray, ap, bp, va - d_dot, vb - d_dot, f_hz, dt_s)
    return ray.group_path_m, rate, ray
