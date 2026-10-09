"""Route two: drift and tilt from two nodes, and which geometries work.

Layer:        study
WP:           WP5, WP7
Claims:       C3
Status:       implemented
Requirements: SCI-15, SCI-23, MIS-10

Technical Volume Section 8.11 proposes that a second node, passing over the
same point a known interval ``dt`` after the first, separates drift from tilt:
the change in range between the passes gives the vertical drift, and the
line-of-sight rate either node measured then gives the tilt.  Section 8.21
(step 3 of the verification sequence) asks for the case to be simulated.  This
study does that, and asks which arrangements of the two nodes work best.

What is simulated
-----------------
1.  **Forward-model verification (two dimensions, free space and plasma).**
    Two passes of the physical-optics plate model
    (:mod:`tss.l1_propagation.plate_echo`) over an extended tilted plate and
    over a compact tilted facet that drifts vertically and, for the facet,
    along track.  The plasma rises with the plate, so the second node sits in
    denser plasma than the first.  The range is measured from the echo's phase
    at two nearby frequencies and the line-of-sight rate from its phase step
    between sweeps; the compact facet's power peak comes from route one.
    Pass criteria fixed before the run: drift within 0.5 m/s, arrival angle
    within 0.01 deg, horizontal drift within 1 m/s.
2.  **Geometry (three dimensions, free space).**  The three ways of phasing a
    second node (:mod:`tss.platform.formation`): trailing in the same plane,
    on the same ground track, and on a crossing track.  Earth's rotation puts
    a same-plane follower tens of kilometres off the first track; what that
    costs is measured on a plane with a cross-track tilt, and for a compact
    facet as whether the second node is inside the facet's scattered beam at
    all.
3.  **Noise.**  Range and line-of-sight rate noise on each pass, against the
    interval between passes.
4.  **Evolution.**  A travelling ionospheric disturbance carried on the
    surface, so that the tilt and drift at the point change between the
    passes.  The truth is the instantaneous value at the middle of the
    interval (and, for the first-pass estimator, at the first pass).
5.  **Crossing tracks.**  The full tilt, both components, from an ascending
    and a descending pass, against the crossing angle.
6.  **Compact facets.**  Drift, tilt and horizontal drift from the route-one
    measurements of two passes.

The three-dimensional forward model
(:mod:`tss.l1_propagation.surface_ranging`) finds the specular point
numerically on a surface that need not be a plane, so the evolution error is
measured, not assumed.  The estimator is
:mod:`tss.l3_signal.novel.second_look_estimator`.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Any

import numpy as np

from tss.constants import C_LIGHT
from tss.l1_propagation.plate_echo import (
    ExponentialTopside,
    Plate,
    PlateEcho,
    echo_free_space,
    echo_stratified,
)
from tss.l1_propagation.surface_ranging import DriftingSurface, sample_surface
from tss.l3_signal.novel.route_one_estimator import estimate_route_one
from tss.l3_signal.novel.second_look_estimator import (
    PassMeasurement,
    estimate_second_look,
    group_path_from_two_frequencies,
    horizontal_drift_from_peaks,
    los_rate_from_field,
)
from tss.platform.formation import two_node_geometry
from tss.studies import StudyResult, register_study
from tss.studies.c3_verify import ENVIRONMENTS, PROFILE, Environment
from tss.studies.route_one import PLATFORM_SPEED_MPS, envelope_width_s

__all__ = [
    "PoResult",
    "po_two_pass",
    "surface_passes",
    "facet_beam_half_width_m",
    "DRIFT_TOL_MPS",
    "ANGLE_TOL_DEG",
    "HDRIFT_TOL_MPS",
]

DRIFT_TOL_MPS = 0.5
ANGLE_TOL_DEG = 0.01
HDRIFT_TOL_MPS = 1.0
REPETITION_HZ = 45.0
DF_HZ = 50.0                 # frequency step for the two-frequency range, P' < 1500 km
MIS10_DT_S = 264.0           # 4.4 minutes, the earliest phasing MIS-10 allows
SIGMA_RATE_MPS = 10.0        # per-pass line-of-sight rate error, route-one budget without fading
SIGMA_RANGE_M = (30.0, 100.0, 300.0)
LAGS_MIN = (1.0, 2.0, 4.4, 6.0, 8.7, 10.0, 15.0, 20.0)
REF_INC_DEG = 30.0
REF_LAT_DEG = 10.0


# ---------------------------------------------------------------------------
# 1. Physical-optics two-pass verification
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class PoResult:
    kind: str
    environment: str
    f_mhz: float
    tilt_deg: float
    drift_true_mps: float
    hdrift_true_mps: float
    dt_s: float
    mu_sat_1: float
    mu_sat_2: float
    mu_facet: float
    drift_est_mps: float
    drift_err_mps: float
    angle_true_deg: float
    angle_est_deg: float
    angle_err_deg: float
    hdrift_est_mps: float
    hdrift_err_mps: float
    doppler_change_mps: float
    verdict: str


def _profile_at(drift_mps: float, t_s: float) -> ExponentialTopside:
    """The profile in the frame of a node at absolute time ``t_s``.

    The medium has risen ``w t`` past the platform's height, so the plasma at
    the platform is the plasma that was ``w t`` deeper at ``t = 0``.
    """
    scale = np.exp(drift_mps * t_s / (2.0 * PROFILE.scale_height_m))
    return ExponentialTopside(PROFILE.f_n_sat_hz * float(scale), PROFILE.scale_height_m)


def _plate_at(plate: Plate, t_s: float) -> Plate:
    return replace(plate, x_m=plate.x_m + plate.drift_x_mps * t_s,
                   depth_m=plate.depth_m - plate.drift_mps * t_s)


def _echo(env: Environment, plate: Plate, prof: ExponentialTopside, t: np.ndarray, f: float) -> PlateEcho:
    if env.plasma:
        return echo_stratified(t, f, plate, prof, v_s_mps=PLATFORM_SPEED_MPS)
    return echo_free_space(t, f, plate, v_s_mps=PLATFORM_SPEED_MPS)


def _normal_offset(env: Environment, depth: float, prof: ExponentialTopside, tilt: float) -> float:
    if not env.plasma:
        return float(depth * np.tan(tilt))
    nodes, weights = np.polynomial.legendre.leggauss(64)
    d = 0.5 * depth * (1.0 + nodes)
    mu2 = prof.mu_squared(env.f_hz, d)
    q = float(np.sqrt(prof.mu_squared(env.f_hz, depth))) * np.sin(tilt)
    return float(np.sum(0.5 * depth * weights * q / np.sqrt(mu2 - q * q)))


def _po_pass_extended(env: Environment, plate: Plate, t_abs: float) -> tuple[PassMeasurement, float]:
    p = _plate_at(plate, t_abs)
    prof = _profile_at(plate.drift_mps, t_abs)
    tau = (np.arange(23) - 11) / REPETITION_HZ
    e = _echo(env, p, prof, tau, env.f_hz)
    rate = los_rate_from_field(tau, e.field, env.f_hz)
    e1 = _echo(env, p, prof, np.array([0.0]), env.f_hz)
    e2 = _echo(env, p, prof, np.array([0.0]), env.f_hz + DF_HZ)
    pg = group_path_from_two_frequencies(e1.field[0], e2.field[0], DF_HZ)
    return PassMeasurement(t_abs, pg, rate, (PLATFORM_SPEED_MPS, 0.0), e.mu_sat), e.mu_sat


def _po_pass_compact(env: Environment, plate: Plate, t_abs: float) -> tuple[PassMeasurement, float, float]:
    p = _plate_at(plate, t_abs)
    prof = _profile_at(plate.drift_mps, t_abs)
    mu_s = float(np.sqrt(prof.mu_squared(env.f_hz, 0.0))) if env.plasma else 1.0
    centre = (p.x_m + _normal_offset(env, p.depth_m, prof, p.tilt_rad)) / (PLATFORM_SPEED_MPS - p.drift_x_mps)
    probe = _echo(env, replace(p, length_m=1.0, tilt_rad=0.0, drift_mps=0.0, drift_x_mps=0.0), prof,
                  np.array([0.0]), env.f_hz)
    w_env = float(envelope_width_s(env.f_hz * mu_s, float(probe.group_path_m[0]), p.length_m))
    dtau = 1.0 / REPETITION_HZ
    n = int(np.ceil(2.5 * w_env / dtau))
    tau = centre + 0.37 * dtau + np.arange(-n, n + 1) * dtau
    e = _echo(env, p, prof, tau, env.f_hz)
    i_pk = int(np.argmax(np.abs(e.field)))
    e1 = _echo(env, p, prof, tau[i_pk:i_pk + 1], env.f_hz)
    e2 = _echo(env, p, prof, tau[i_pk:i_pk + 1], env.f_hz + DF_HZ)
    pg_pk = group_path_from_two_frequencies(e1.field[0], e2.field[0], DF_HZ)
    est = estimate_route_one(tau, e.field, f_hz=env.f_hz, group_path_m=pg_pk, mu_sat=e.mu_sat,
                             v_s_mps=PLATFORM_SPEED_MPS)
    if not est.separable:
        raise RuntimeError(f"route one declined a compact facet: {est.reason}")
    t_pk = est.t_power_peak_s
    f1 = _echo(env, p, prof, np.array([t_pk]), env.f_hz)
    f2 = _echo(env, p, prof, np.array([t_pk]), env.f_hz + DF_HZ)
    pg = group_path_from_two_frequencies(f1.field[0], f2.field[0], DF_HZ)
    rate = est.drift_mps * e.mu_sat     # route one's drift is r / mu_sat at the power peak
    x_pk = PLATFORM_SPEED_MPS * t_pk    # Earth-fixed: the node passes x = 0 at tau = 0
    return PassMeasurement(t_abs + t_pk, pg, rate, (PLATFORM_SPEED_MPS, 0.0), e.mu_sat), e.mu_sat, x_pk


def po_two_pass(env: Environment, *, extended: bool, tilt_deg: float, drift_mps: float,
                hdrift_mps: float = 0.0, dt_s: float = MIS10_DT_S) -> PoResult:
    """Two passes of the physical-optics model and the second-look estimate."""
    a = np.radians(tilt_deg)
    depth = env.plate_depth()
    plate = Plate(x_m=0.0, depth_m=depth, tilt_rad=a, length_m=80e3 if extended else 2.0e3,
                  drift_mps=drift_mps, drift_x_mps=0.0 if extended else hdrift_mps,
                  taper="tukey" if extended else "uniform")
    mu_f = float(np.sqrt(PROFILE.mu_squared(env.f_hz, depth))) if env.plasma else 1.0
    if extended:
        m1, mu1 = _po_pass_extended(env, plate, 0.0)
        m2, mu2 = _po_pass_extended(env, plate, dt_s)
        est = estimate_second_look(m1, m2)
        u_est = float("nan")
    else:
        m1, mu1, x1 = _po_pass_compact(env, plate, 0.0)
        m2, mu2, x2 = _po_pass_compact(env, plate, dt_s)
        est = estimate_second_look(m1, m2, facet=True)
        u_est = horizontal_drift_from_peaks(x1, x2, m2.t_s - m1.t_s, est.drift_mps, est.slope_along)
        # The facet drifts with the medium, so the platform's speed relative to
        # it is v - u.  Re-solve once with the measured u.
        v_rel = (PLATFORM_SPEED_MPS - u_est, 0.0)
        est = estimate_second_look(replace(m1, velocity_mps=v_rel), replace(m2, velocity_mps=v_rel), facet=True)
        u_est = horizontal_drift_from_peaks(x1, x2, m2.t_s - m1.t_s, est.drift_mps, est.slope_along)
    # Truth: the arrival angle at the spacecraft has sin = mu_f sin(a) / mu_sat,
    # which differs between the passes because mu_sat does; the estimate is
    # their mean.
    s_true = 0.5 * (mu_f * np.sin(a) / mu1 + mu_f * np.sin(a) / mu2)
    ang_true = float(np.degrees(np.arcsin(s_true)))
    ang_est = float(np.degrees(np.arcsin(est.slope_along)))
    d_err = est.drift_mps - drift_mps
    a_err = ang_est - ang_true
    u_err = u_est - hdrift_mps if not extended else float("nan")
    ok = abs(d_err) <= DRIFT_TOL_MPS and abs(a_err) <= ANGLE_TOL_DEG
    if not extended:
        ok = ok and abs(u_err) <= HDRIFT_TOL_MPS
    return PoResult(
        kind="extended" if extended else "compact", environment=env.name, f_mhz=env.f_hz / 1e6,
        tilt_deg=tilt_deg, drift_true_mps=drift_mps, hdrift_true_mps=0.0 if extended else hdrift_mps,
        dt_s=float(m2.t_s - m1.t_s), mu_sat_1=mu1, mu_sat_2=mu2, mu_facet=mu_f,
        drift_est_mps=est.drift_mps, drift_err_mps=d_err, angle_true_deg=ang_true,
        angle_est_deg=ang_est, angle_err_deg=a_err, hdrift_est_mps=u_est, hdrift_err_mps=u_err,
        doppler_change_mps=est.doppler_change_mps, verdict="pass" if ok else "fail",
    )


def _po_checks() -> list[PoResult]:
    env = {e.name: e for e in ENVIRONMENTS}
    names = ("free_5MHz_500km", "plasma_5MHz_mu0.5", "plasma_3MHz_mu0.5")
    out: list[PoResult] = []
    for name in names:
        for tilt in (0.0, 0.5, 1.0):
            for w in (25.0, -50.0):
                for dt in (MIS10_DT_S, 600.0):
                    out.append(po_two_pass(env[name], extended=True, tilt_deg=tilt, drift_mps=w, dt_s=dt))
    for name in names:
        for tilt in (0.25, 0.5, 1.0):
            for u in (0.0, 100.0):
                out.append(po_two_pass(env[name], extended=False, tilt_deg=tilt, drift_mps=25.0,
                                       hdrift_mps=u))
    return out


# ---------------------------------------------------------------------------
# 2-6. Three-dimensional geometry, noise and evolution
# ---------------------------------------------------------------------------


def surface_passes(surface: DriftingSurface, geom) -> tuple[PassMeasurement, PassMeasurement, np.ndarray]:
    """Noise-free measurements of both nodes of a formation over a surface.

    The point is the origin of the local frame; node 1 passes over it at
    ``t = 0`` and node 2 at ``geom.dt_s`` from its own closest approach.
    """
    out = []
    feet = []
    for k, node in enumerate((geom.node1, geom.node2)):
        t = 0.0 if k == 0 else geom.node2.t_s - geom.node1.t_s
        off = node.offset_en_m - geom.node1.offset_en_m
        p0 = np.array([off[0], off[1], 0.0])
        smp = sample_surface(surface, p0, node.velocity_en_mps, t)
        out.append(PassMeasurement(t, smp.range_m, smp.los_rate_mps, tuple(node.velocity_en_mps)))
        feet.append(smp.foot)
    return out[0], out[1], np.array(feet)


def facet_beam_half_width_m(f_hz: float, group_path_m: float, facet_m: float) -> float:
    """Half the half-power width of a compact facet's scattered beam, at the spacecraft."""
    return 0.5 * 0.443 * (C_LIGHT / f_hz) * group_path_m / facet_m


def _noise_stats(m1: PassMeasurement, m2: PassMeasurement, w_true: float, s_true: np.ndarray,
                 sigma_r_m: float, sigma_rate: float, rng, n: int = 4000, crossing: bool = False,
                 facet: bool = False) -> dict[str, float]:
    """Rms errors of drift and tilt with Gaussian range and rate noise."""
    e_w, e_al, e_x, e_y, e_tilt, e_p1 = [], [], [], [], [], []
    u1 = np.asarray(m1.velocity_mps) / np.linalg.norm(m1.velocity_mps)
    for _ in range(n):
        a = replace(m1, group_path_m=m1.group_path_m + rng.normal(0, sigma_r_m),
                    los_rate_mps=m1.los_rate_mps + rng.normal(0, sigma_rate))
        b = replace(m2, group_path_m=m2.group_path_m + rng.normal(0, sigma_r_m),
                    los_rate_mps=m2.los_rate_mps + rng.normal(0, sigma_rate))
        est = estimate_second_look(a, b, facet=facet)
        e_w.append(est.drift_mps - w_true)
        e_al.append(np.degrees(est.slope_along - float(np.dot(s_true, u1))))
        if not crossing:
            e_p1.append(np.degrees(est.slope_along_1 - float(np.dot(s_true, u1))))
        if crossing:
            e_x.append(np.degrees(est.slope_x - s_true[0]))
            e_y.append(np.degrees(est.slope_y - s_true[1]))
            e_tilt.append(np.degrees(np.hypot(est.slope_x, est.slope_y) - np.hypot(*s_true)))
    rms = lambda v: float(np.sqrt(np.mean(np.square(v)))) if v else float("nan")  # noqa: E731
    return {"drift_rms_mps": rms(e_w), "along_rms_deg": rms(e_al), "east_rms_deg": rms(e_x),
            "north_rms_deg": rms(e_y), "tilt_rms_deg": rms(e_tilt), "per_pass_rms_deg": rms(e_p1)}


def _geometry_table() -> list[dict[str, Any]]:
    rows = []
    cases = [(30.0, lat) for lat in (0.0, 10.0, 20.0)] + [(45.0, lat) for lat in (0.0, 20.0, 40.0)] \
        + [(60.0, lat) for lat in (0.0, 30.0, 50.0)] + [(97.4, lat) for lat in (0.0, 30.0, 60.0)]
    cross_tilt = np.tan(np.radians(1.0))
    for inc, lat in cases:
        for mode in ("same_plane", "same_ground_track", "crossing"):
            g = two_node_geometry(inc, lat, MIS10_DT_S, mode)
            # A plane tilted 1 degree across node 1's track, drifting at 25 m/s.
            u = g.node1.velocity_en_mps / np.linalg.norm(g.node1.velocity_en_mps)
            c = np.array([u[1], -u[0]])
            grad = cross_tilt * c
            surf = DriftingSurface(z0_m=-500e3, a=float(grad[0]), b=float(grad[1]), w_mps=25.0)
            m1, m2, _ = surface_passes(surf, g)
            est = estimate_second_look(m1, m2)
            n = surf.unit_normal(0.0, 0.0, 0.0)
            rows.append({
                "inclination_deg": inc, "latitude_deg": lat, "mode": mode,
                "dt_s": g.dt_s, "heading1_deg": g.node1.heading_deg, "heading2_deg": g.node2.heading_deg,
                "speed_mps": float(np.linalg.norm(g.node1.velocity_en_mps)),
                "line_crossing_deg": float(min(g.crossing_angle_deg, 180.0 - g.crossing_angle_deg)),
                "cross_track_offset_km": g.cross_track_offset_m / 1e3,
                "cross_track_rate_km_per_min": g.cross_track_offset_m / 1e3 / (MIS10_DT_S / 60.0),
                "raan_offset_deg": g.raan_offset_deg if mode == "same_ground_track" else 0.0,
                "plane_change_dv_mps": g.plane_change_dv_mps,
                "j2_drift_days_20km": g.j2_drift_days, "j2_drift_dv_mps": g.j2_drift_dv_mps,
                "drift_err_1deg_cross_tilt_mps": est.drift_mps - 25.0,
                "along_err_1deg_cross_tilt_deg": float(np.degrees(est.slope_along - np.dot(n[:2], u))),
                "facet_visible_5MHz_2km": abs(g.cross_track_offset_m)
                <= facet_beam_half_width_m(5e6, 500e3, 2e3),
            })
    return rows


def _phasing_cost_table() -> list[dict[str, Any]]:
    """What a same-ground-track follower costs, against the interval."""
    rows = []
    for inc in (30.0, 45.0, 60.0, 97.4):
        for lag in (1.0, 2.0, 4.4, 8.7):
            g = two_node_geometry(inc, 0.0, lag * 60.0, "same_ground_track")
            sp = two_node_geometry(inc, 0.0, lag * 60.0, "same_plane")
            rows.append({"inclination_deg": inc, "lag_min": lag, "raan_offset_deg": g.raan_offset_deg,
                         "plane_angle_deg": g.plane_angle_deg, "plane_change_dv_mps": g.plane_change_dv_mps,
                         "j2_drift_days_20km": g.j2_drift_days, "j2_drift_dv_mps": g.j2_drift_dv_mps,
                         "same_plane_cross_offset_equator_km": sp.cross_track_offset_m / 1e3})
    return rows


def _lag_table(rng) -> list[dict[str, Any]]:
    rows = []
    surf = DriftingSurface(z0_m=-500e3, a=np.tan(np.radians(0.5)), b=0.0, w_mps=25.0)
    n = surf.unit_normal(0.0, 0.0, 0.0)
    for lag in LAGS_MIN:
        g = two_node_geometry(REF_INC_DEG, REF_LAT_DEG, lag * 60.0, "same_ground_track")
        m1, m2, _ = surface_passes(surf, g)
        for sr in SIGMA_RANGE_M:
            st = _noise_stats(m1, m2, 25.0, n[:2], sr, SIGMA_RATE_MPS, rng, n=3000)
            rows.append({"lag_min": lag, "sigma_range_m": sr, "sigma_rate_mps": SIGMA_RATE_MPS,
                         "drift_rms_mps": st["drift_rms_mps"], "along_tilt_rms_deg": st["along_rms_deg"],
                         "along_tilt_per_pass_rms_deg": st["per_pass_rms_deg"],
                         "drift_closed_form_mps": np.sqrt(2.0) * sr / (lag * 60.0)})
    return rows


def _tid_table(rng) -> list[dict[str, Any]]:
    """Evolution error: a travelling disturbance carried on the surface."""
    rows = []
    lam_h = 300e3
    k = 2.0 * np.pi / lam_h
    tilt_amp = np.radians(2.0)
    amp = np.tan(tilt_amp) / k
    tol = 3.0 * np.sqrt(2.0) * SIGMA_RATE_MPS
    for period_min in (15.0, 30.0, 60.0):
        om = 2.0 * np.pi / (period_min * 60.0)
        for lag in LAGS_MIN[:-1]:
            g = two_node_geometry(REF_INC_DEG, REF_LAT_DEG, lag * 60.0, "same_ground_track")
            u1 = g.node1.velocity_en_mps / np.linalg.norm(g.node1.velocity_en_mps)
            ew_mean, es_mean, ew_first, es_pass, flags, ew_interval = [], [], [], [], [], []
            for az in np.radians(np.arange(0.0, 360.0, 45.0)):
                for ph in np.linspace(0.0, 2.0 * np.pi, 16, endpoint=False):
                    surf = DriftingSurface(z0_m=-500e3, w_mps=0.0, amp_m=amp, omega_rad_s=om,
                                           kx_rad_m=k * np.sin(az), ky_rad_m=k * np.cos(az), phase_rad=ph)
                    m1, m2, feet = surface_passes(surf, g)
                    dt = m2.t_s - m1.t_s
                    tm = 0.5 * dt
                    xf, yf = 0.5 * (feet[0, 0] + feet[1, 0]), 0.5 * (feet[0, 1] + feet[1, 1])
                    w_mid = surf.vertical_velocity(xf, yf, tm)
                    s_mid = float(np.dot(surf.unit_normal(xf, yf, tm)[:2], u1))
                    w_t1 = surf.vertical_velocity(feet[0, 0], feet[0, 1], 0.0)
                    s_t1 = float(np.dot(surf.unit_normal(feet[0, 0], feet[0, 1], 0.0)[:2], u1))
                    s_t2 = float(np.dot(surf.unit_normal(feet[1, 0], feet[1, 1], dt)[:2], u1))
                    w_int = float((surf.height(xf, yf, dt) - surf.height(xf, yf, 0.0)) / dt)
                    e_m = estimate_second_look(m1, m2, evolution_tol_mps=tol)
                    ew_mean.append(e_m.drift_mps - w_mid)
                    es_mean.append(np.degrees(e_m.slope_along - s_mid))
                    ew_first.append(e_m.drift_mps - w_t1)
                    es_pass.append(np.degrees(e_m.slope_along_1 - s_t1))
                    es_pass.append(np.degrees(e_m.slope_along_2 - s_t2))
                    ew_interval.append(e_m.drift_mps - w_int)
                    flags.append(e_m.flagged)
            rms = lambda v: float(np.sqrt(np.mean(np.square(v))))  # noqa: E731
            rows.append({
                "period_min": period_min, "lag_min": lag, "tilt_amplitude_deg": 2.0,
                "wavelength_km": lam_h / 1e3, "drift_amplitude_mps": amp * om,
                "drift_rms_vs_interval_mean_mps": rms(ew_interval),
                "drift_rms_vs_mid_mps": rms(ew_mean), "drift_rms_vs_first_pass_mps": rms(ew_first),
                "mean_tilt_rms_vs_mid_deg": rms(es_mean), "per_pass_tilt_rms_deg": rms(es_pass),
                "flagged_fraction": float(np.mean(flags)),
            })
    return rows


def _crossing_table(rng) -> list[dict[str, Any]]:
    rows = []
    tilt = np.radians(1.0)
    az = np.radians(30.0)
    grad = -np.tan(tilt) * np.array([np.sin(az), np.cos(az)])   # normal leans toward azimuth 30 deg
    surf = DriftingSurface(z0_m=-500e3, a=float(grad[0]), b=float(grad[1]), w_mps=25.0)
    n = surf.unit_normal(0.0, 0.0, 0.0)
    cases = [(30.0, None, lat) for lat in (0.0, 10.0, 20.0, 25.0)] \
        + [(45.0, None, lat) for lat in (0.0, 20.0, 40.0)] \
        + [(60.0, None, lat) for lat in (0.0, 30.0, 50.0)] \
        + [(97.4, None, lat) for lat in (0.0, 30.0, 60.0)] \
        + [(30.0, 97.4, lat) for lat in (0.0, 10.0, 20.0)]
    for inc, inc2, lat in cases:
        g = two_node_geometry(inc, lat, MIS10_DT_S, "crossing", inclination2_deg=inc2)
        m1, m2, _ = surface_passes(surf, g)
        est = estimate_second_look(m1, m2)
        st = _noise_stats(m1, m2, 25.0, n[:2], 100.0, SIGMA_RATE_MPS, rng, n=3000, crossing=True)
        u1 = g.node1.velocity_en_mps / np.linalg.norm(g.node1.velocity_en_mps)
        c1 = np.array([u1[1], -u1[0]])
        # Rms of the component across node 1's track, from the east/north errors.
        rows.append({
            "inclination_deg": inc, "inclination2_deg": inc if inc2 is None else inc2, "latitude_deg": lat,
            "line_crossing_deg": float(min(g.crossing_angle_deg, 180.0 - g.crossing_angle_deg)),
            "noise_free_drift_err_mps": est.drift_mps - 25.0,
            "noise_free_tilt_err_deg": float(np.degrees(np.hypot(est.slope_x, est.slope_y) - np.hypot(*n[:2]))),
            "drift_rms_mps": st["drift_rms_mps"], "east_rms_deg": st["east_rms_deg"],
            "north_rms_deg": st["north_rms_deg"], "tilt_rms_deg": st["tilt_rms_deg"],
            "worst_axis_rms_deg": float(max(st["east_rms_deg"], st["north_rms_deg"])),
            "cross_axis_dir_deg": float(np.degrees(np.arctan2(c1[0], c1[1]))),
        })
    return rows


def _compact_table(rng) -> list[dict[str, Any]]:
    """Compact facet budget from the route-one per-pass errors (TV Table 8.14)."""
    rows = []
    cases = {
        # (sigma of route-one rate r/mu per pass, sigma of power-peak time), 5 MHz, 500 km, 2 km facet
        "no_fading": (10.8, 0.033),
        "fading_0.5s": (23.1, 0.172),
    }
    for label, (s_rate, s_tpk) in cases.items():
        sx = PLATFORM_SPEED_MPS * s_tpk
        for lag in LAGS_MIN:
            dt = lag * 60.0
            for sr in SIGMA_RANGE_M:
                w = np.sqrt(2.0) * sr / dt
                tilt = np.degrees(np.sqrt(s_rate**2 / 2.0 + w**2) / PLATFORM_SPEED_MPS)
                u = np.sqrt(2.0) * sx / dt
                rows.append({"fading": label, "lag_min": lag, "sigma_range_m": sr,
                             "sigma_rate_mps": s_rate, "sigma_peak_time_s": s_tpk,
                             "drift_rms_mps": w, "along_tilt_rms_deg": tilt, "hdrift_rms_mps": u,
                             "match_window_along_km_at_250mps": 250.0 * dt / 1e3,
                             "match_window_range_km_at_100mps": 100.0 * dt / 1e3})
    return rows


def _write_csv(path: Path, rows: list[dict[str, Any]], header: str = "") -> Path:
    fields = list(rows[0].keys())
    with open(path, "w", encoding="utf-8") as fh:
        if header:
            fh.write(f"# {header}\n")
        fh.write(",".join(fields) + "\n")
        for r in rows:
            vals = []
            for k in fields:
                v = r[k]
                if isinstance(v, (bool, np.bool_)):
                    vals.append(str(bool(v)))
                elif isinstance(v, str):
                    vals.append(f'"{v}"')
                elif isinstance(v, (float, np.floating)):
                    vals.append(f"{float(v):.6g}")
                else:
                    vals.append(str(v))
            fh.write(",".join(vals) + "\n")
    return path


def _figure(lag, tid, geom, cross, out_dir: Path) -> list[Path]:
    try:
        import matplotlib  # noqa: F401
    except ImportError:  # pragma: no cover
        return []
    from tss.figures import style

    style.use()
    fig, axes = style.figure(2, 2, figsize=(9.6, 7.0))
    ax = axes[0, 0]
    for k, sr in enumerate(SIGMA_RANGE_M):
        r = [x for x in lag if x["sigma_range_m"] == sr]
        ax.loglog([x["lag_min"] for x in r], [x["drift_rms_mps"] for x in r], "o-", ms=3,
                  color=style.SERIES[k], label=f"noise, range {sr:g} m")
    for k, per in enumerate((15.0, 30.0, 60.0)):
        r = [x for x in tid if x["period_min"] == per]
        ax.loglog([x["lag_min"] for x in r], [x["drift_rms_vs_mid_mps"] for x in r], "s--", ms=3,
                  color=style.SERIES[3 + k], label=f"TID {per:g} min, against mid-interval drift")
    ax.axvline(MIS10_DT_S / 60.0, color=style.INK_MUTED, ls=":", lw=0.8)
    ax.set_xlabel("interval between passes (min)")
    ax.set_ylabel("drift error, rms (m/s)")
    ax.set_title("(a) Drift: noise against evolution", fontsize=9)
    ax.legend(fontsize=6, frameon=False)

    ax = axes[0, 1]
    r = [x for x in lag if x["sigma_range_m"] == 100.0]
    ax.loglog([x["lag_min"] for x in r], [x["along_tilt_per_pass_rms_deg"] for x in r], "o-", ms=3,
              color=style.SERIES[1], label="noise, per pass (range 100 m, rate 10 m/s)")
    for k, per in enumerate((15.0, 30.0, 60.0)):
        r = [x for x in tid if x["period_min"] == per]
        ax.loglog([x["lag_min"] for x in r], [x["mean_tilt_rms_vs_mid_deg"] for x in r], "s--", ms=3,
                  color=style.SERIES[3 + k], label=f"TID {per:g} min, mean of the passes")
        ax.loglog([x["lag_min"] for x in r], [x["per_pass_tilt_rms_deg"] for x in r], "^:", ms=3,
                  color=style.SERIES[3 + k], label=f"TID {per:g} min, each pass")
    ax.axhline(0.29, color=style.INK_MUTED, ls="--", lw=0.8)
    ax.text(19.0, 0.32, "route-three goal, 0.29 deg", fontsize=6, color=style.INK_MUTED, ha="right")
    ax.axvline(MIS10_DT_S / 60.0, color=style.INK_MUTED, ls=":", lw=0.8)
    ax.set_xlabel("interval between passes (min)")
    ax.set_ylabel("along-track tilt error, rms (deg)")
    ax.set_title("(b) Tilt: noise against evolution", fontsize=9)
    ax.set_ylim(5e-4, 8.0)
    ax.legend(fontsize=5.5, frameon=False, ncol=2, loc="upper left")

    ax = axes[1, 0]
    sp = [x for x in geom if x["mode"] == "same_plane"]
    for k, inc in enumerate((30.0, 45.0, 60.0, 97.4)):
        r = [x for x in sp if x["inclination_deg"] == inc]
        ax.plot([x["latitude_deg"] for x in r], [abs(x["cross_track_offset_km"]) for x in r], "o-", ms=3,
                color=style.SERIES[k], label=f"i = {inc:g} deg")
    ax.axhline(facet_beam_half_width_m(5e6, 500e3, 2e3) / 1e3, color=style.INK_MUTED, ls="--", lw=0.8)
    ax.text(1.0, facet_beam_half_width_m(5e6, 500e3, 2e3) / 1e3 + 3.0,
            "2 km facet's beam, 5 MHz", fontsize=6, color=style.INK_MUTED)
    ax.set_xlabel("latitude of the point (deg)")
    ax.set_ylabel("cross-track offset after 4.4 min (km)")
    ax.set_title("(c) Same-plane follower: Earth-rotation offset", fontsize=9)
    ax.legend(fontsize=6, frameon=False)

    ax = axes[1, 1]
    same_inc = [x for x in cross if x["inclination2_deg"] == x["inclination_deg"]]
    ax.plot([x["line_crossing_deg"] for x in same_inc], [x["worst_axis_rms_deg"] for x in same_inc], "o",
            ms=4, color=style.SERIES[0], label="ascending x descending")
    mixed = [x for x in cross if x["inclination2_deg"] != x["inclination_deg"]]
    ax.plot([x["line_crossing_deg"] for x in mixed], [x["worst_axis_rms_deg"] for x in mixed], "s",
            ms=4, color=style.SERIES[2], label="30 deg x 97.4 deg")
    gam = np.linspace(5.0, 90.0, 100)
    base = np.degrees(SIGMA_RATE_MPS / PLATFORM_SPEED_MPS)
    ax.plot(gam, base / (np.sqrt(2.0) * np.sin(np.radians(gam) / 2.0)), "-", color=style.INK_MUTED, lw=0.8,
            label="worst principal axis, sigma_r / (sqrt2 v sin(gamma/2))")
    ax.set_xlabel("angle between the tracks (deg)")
    ax.set_ylabel("worst tilt component, rms (deg)")
    ax.set_yscale("log")
    ax.set_title("(d) Crossing tracks: full tilt", fontsize=9)
    ax.legend(fontsize=6, frameon=False)

    fig.tight_layout()
    return style.savefig(fig, out_dir / "route_two", formats=("png", "pdf"))


@register_study(
    "route_two",
    description="Route two: drift and tilt from two nodes, by geometry (TV 8.11, 8.21 step 3)",
    work_packages=["WP5", "WP7"],
    claims=["C3"],
    status="implemented",
    requirements=["SCI-15", "SCI-23", "MIS-10"],
)
def run(out_dir: Path, *, config: dict[str, Any], manifest: Any) -> StudyResult:
    opts = dict(config.get("route_two", {}))
    seed = int(opts.get("seed", 20261009))
    rng = np.random.default_rng(seed)

    po = _po_checks()
    geom = _geometry_table()
    phasing = _phasing_cost_table()
    lag = _lag_table(rng)
    tid = _tid_table(rng)
    cross = _crossing_table(rng)
    compact = _compact_table(rng)

    outputs = [
        _write_csv(out_dir / "route_two_po_checks.csv", [asdict(r) for r in po],
                   f"pass criteria: drift within {DRIFT_TOL_MPS} m/s, angle within {ANGLE_TOL_DEG} deg, "
                   f"horizontal drift within {HDRIFT_TOL_MPS} m/s"),
        _write_csv(out_dir / "route_two_geometry.csv", geom, "dt 4.4 min; plane tilted 1 deg across node 1's track"),
        _write_csv(out_dir / "route_two_phasing_cost.csv", phasing,
                   "same-ground-track follower: node offset, plane change, or J2 drift with a 20 km offset"),
        _write_csv(out_dir / "route_two_lag.csv", lag, f"same ground track, i {REF_INC_DEG}, lat {REF_LAT_DEG}"),
        _write_csv(out_dir / "route_two_tid.csv", tid, "TID 300 km, 2 deg tilt amplitude; no noise"),
        _write_csv(out_dir / "route_two_crossing.csv", cross, "plane 1 deg toward az 30, 25 m/s; noise 100 m, 10 m/s"),
        _write_csv(out_dir / "route_two_compact.csv", compact, "5 MHz, 500 km, 2 km facet; TV Table 8.14 inputs"),
    ]
    outputs += _figure(lag, tid, geom, cross, out_dir)

    ext = [r for r in po if r.kind == "extended"]
    cmp_ = [r for r in po if r.kind == "compact"]
    lag100 = {x["lag_min"]: x for x in lag if x["sigma_range_m"] == 100.0}
    sp = [x for x in geom if x["mode"] == "same_plane"]
    metrics = {
        "seed": seed,
        "po_extended_pass": sum(r.verdict == "pass" for r in ext), "po_extended_n": len(ext),
        "po_compact_pass": sum(r.verdict == "pass" for r in cmp_), "po_compact_n": len(cmp_),
        "po_max_drift_err_mps": float(max(abs(r.drift_err_mps) for r in po)),
        "po_max_angle_err_deg": float(max(abs(r.angle_err_deg) for r in po)),
        "po_max_hdrift_err_mps": float(max(abs(r.hdrift_err_mps) for r in cmp_)),
        "lag_4p4_drift_rms_mps": lag100[4.4]["drift_rms_mps"],
        "lag_4p4_tilt_rms_deg": lag100[4.4]["along_tilt_rms_deg"],
        "same_plane_cross_offset_km_range": [float(min(abs(x["cross_track_offset_km"]) for x in sp)),
                                             float(max(abs(x["cross_track_offset_km"]) for x in sp))],
        "same_plane_drift_err_per_deg_cross_tilt_mps_range": [
            float(min(abs(x["drift_err_1deg_cross_tilt_mps"]) for x in sp)),
            float(max(abs(x["drift_err_1deg_cross_tilt_mps"]) for x in sp))],
    }
    summary = f"""
Route two, two nodes.  Forward-model checks (physical optics, free space and
plasma): extended {metrics["po_extended_pass"]}/{metrics["po_extended_n"]}, compact {metrics["po_compact_pass"]}/{metrics["po_compact_n"]} pass;
worst drift error {metrics["po_max_drift_err_mps"]:.3g} m/s, worst angle error {metrics["po_max_angle_err_deg"]:.2g} deg,
worst horizontal-drift error {metrics["po_max_hdrift_err_mps"]:.3g} m/s.
Same ground track, 4.4 min, range 100 m, rate 10 m/s: drift {metrics["lag_4p4_drift_rms_mps"]:.2f} m/s rms,
along-track tilt {metrics["lag_4p4_tilt_rms_deg"]:.3f} deg rms.
Same-plane follower at 4.4 min sits {metrics["same_plane_cross_offset_km_range"][0]:.0f} to {metrics["same_plane_cross_offset_km_range"][1]:.0f} km
off track; a 1 deg cross-track tilt then costs {metrics["same_plane_drift_err_per_deg_cross_tilt_mps_range"][0]:.1f} to
{metrics["same_plane_drift_err_per_deg_cross_tilt_mps_range"][1]:.1f} m/s of drift.  See the CSVs for evolution and crossing tracks.
""".strip()
    return StudyResult(outputs=outputs, summary=summary, metrics=metrics)
