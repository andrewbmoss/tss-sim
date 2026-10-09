"""C3 verification: the noise-free route-one and extended-plane tests.

Layer:        study
WP:           WP5
Claims:       C3
Status:       implemented
Requirements: SCI-15, SCI-23

Technical Volume Section 8.15 asks for two tests before anything else in C3.

*  A compact tilted facet with a known drift, processed twice.  With the
   overflight time left free the fit must show the degeneracy of TV
   Figure 8.x.  With the overflight time taken from the echo power it must
   recover the drift and the tilt against truth.
*  The same drift on an extended tilted plane.  The estimator must report that
   it cannot separate drift from tilt rather than return a biased drift.

Each is run in free space and inside a smooth plasma profile.  Committee
Review 5 made the pair condition 1 for signing the proposal.

How the test is kept honest
---------------------------
The forward model (:mod:`tss.l1_propagation.plate_echo`) is a physical-optics
sum over a flat plate with exact geometry, and inside the plasma an exact
Snell's-law ray to every point of the plate.  It shares no equation with the
estimator (:mod:`tss.l3_signal.novel.route_one_estimator`), which assumes the
linearised Doppler history of TV equation 8.x and the ``sinc^2`` envelope of
equation 8.y.  Whatever the estimator recovers it has to find in the field.

Cases
-----
A   Untilted compact facet, at nadir and displaced along track, with drift.
    Tests the degeneracy with the overflight time free, and drift and
    along-track position with it taken from the power.
B   Compact facet whose normal is tilted, with drift.  This is the case the
    TV's test names.  The power peaks when the platform crosses the facet
    normal, which for a tilted facet is not the time of closest approach, so
    the estimator is expected to return drift plus the tilt's line-of-sight
    rate, ``v_s mu_f sin(a) / mu_sat``, with nothing to flag it.
C   Extended plate, tapered, tilted, with drift.  The estimator must decline.
D   Plate length from 0.5 to 32 km, untilted, at 25 m/s: where the
    identifiability rule switches from compact to extended.

Pass criteria
-------------
These were set before the run (session of 2026-10-08) and are recorded in the
CSV and the summary: drift within 1 m/s of truth, along-track arrival angle
within 0.01 deg, and for the extended plate a refusal to separate.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
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
from tss.l3_signal.novel.route_one_estimator import (
    doppler_history_from_field,
    estimate_route_one,
    free_t0_degeneracy,
)
from tss.studies import StudyResult, register_study
from tss.studies.route_one import PLATFORM_SPEED_MPS, envelope_width_s

__all__ = [
    "Environment",
    "TrialResult",
    "ENVIRONMENTS",
    "run_trial",
    "normal_crossing_offset_m",
    "predicted_tilt_bias_mps",
]

DRIFT_TOL_MPS = 1.0
ANGLE_TOL_DEG = 0.01
REPETITION_HZ = 45.0
PROFILE = ExponentialTopside(f_n_sat_hz=1.5e6, scale_height_m=200e3)


@dataclass(frozen=True)
class Environment:
    """Where a trial is run.  ``mu_facet`` is ignored in free space."""

    name: str
    f_hz: float
    plasma: bool
    depth_m: float = 500e3
    mu_facet: float = 1.0

    def plate_depth(self) -> float:
        if not self.plasma:
            return self.depth_m
        return PROFILE.depth_for_mu(self.f_hz, self.mu_facet)


ENVIRONMENTS: tuple[Environment, ...] = (
    Environment("free_3MHz_500km", 3.0e6, False, 500e3),
    Environment("free_5MHz_500km", 5.0e6, False, 500e3),
    Environment("free_10MHz_500km", 10.0e6, False, 500e3),
    Environment("free_5MHz_1000km", 5.0e6, False, 1000e3),
    Environment("plasma_3MHz_mu0.5", 3.0e6, True, mu_facet=0.5),
    Environment("plasma_5MHz_mu0.5", 5.0e6, True, mu_facet=0.5),
    Environment("plasma_10MHz_mu0.8", 10.0e6, True, mu_facet=0.8),
)


@dataclass(frozen=True)
class TrialResult:
    case: str
    environment: str
    f_mhz: float
    plate_km: float
    x_true_km: float
    tilt_deg: float
    drift_true_mps: float
    mu_sat: float
    mu_facet: float
    depth_km: float
    group_path_km: float
    separable: bool
    reason: str
    slope_ratio: float
    peak_contrast: float
    drift_est_mps: float
    drift_err_mps: float
    predicted_bias_mps: float
    angle_err_deg: float
    los_rate_mps: float
    los_rate_predicted_mps: float
    degeneracy_spread_hz2: float
    degeneracy_reference_hz2: float
    verdict: str


def _mu_at_depth(f_hz: float, depth_m: float) -> float:
    return float(np.sqrt(PROFILE.mu_squared(f_hz, depth_m)))


def normal_crossing_offset_m(env: Environment, tilt_rad: float) -> float:
    """Along-track offset of the platform from the plate when it is on the normal.

    In free space ``h tan(a)``.  In the plasma the ray must arrive at the plate
    along its normal, so ``q = mu_f sin(a)``, and the offset is ``X(q)``.
    """
    depth = env.plate_depth()
    if not env.plasma:
        return float(depth * np.tan(tilt_rad))
    nodes, weights = np.polynomial.legendre.leggauss(64)
    d = 0.5 * depth * (1.0 + nodes)
    mu2 = PROFILE.mu_squared(env.f_hz, d)
    q = _mu_at_depth(env.f_hz, depth) * np.sin(tilt_rad)
    return float(np.sum(0.5 * depth * weights * q / np.sqrt(mu2 - q * q)))


def predicted_tilt_bias_mps(env: Environment, tilt_rad: float) -> float:
    """Route one's drift error for a facet whose normal is tilted.

    At the power peak the ray meets the facet along its normal, so its Snell
    invariant is ``q = mu_f sin(a)`` and the platform's contribution to the
    line-of-sight rate is ``v_s q``.  Route one converts the Doppler at the
    power peak into drift by dividing by ``mu_sat``.  The facet normal leans
    toward ``+x``, so the platform crosses it after closest approach, while
    receding, and the bias is negative.
    """
    mu_f = _mu_at_depth(env.f_hz, env.plate_depth()) if env.plasma else 1.0
    mu_s = _mu_at_depth(env.f_hz, 0.0) if env.plasma else 1.0
    return float(-PLATFORM_SPEED_MPS * mu_f * np.sin(tilt_rad) / mu_s)


def _echo(env: Environment, plate: Plate, t: np.ndarray) -> PlateEcho:
    if env.plasma:
        return echo_stratified(t, env.f_hz, plate, PROFILE, v_s_mps=PLATFORM_SPEED_MPS)
    return echo_free_space(t, env.f_hz, plate, v_s_mps=PLATFORM_SPEED_MPS)


def _window(env: Environment, plate: Plate, t_centre: float, half_s: float, rng) -> np.ndarray:
    dt = 1.0 / REPETITION_HZ
    n = int(np.ceil(half_s / dt))
    # Offset the grid by a fraction of a sample so no event sits on it.
    return t_centre + float(rng.uniform(-0.5, 0.5)) * dt + np.arange(-n, n + 1) * dt


def run_trial(
    case: str,
    env: Environment,
    *,
    plate_m: float = 2.0e3,
    x_m: float = 0.0,
    tilt_deg: float = 0.0,
    drift_mps: float = 0.0,
    extended: bool = False,
    half_window_s: float | None = None,
    rng: np.random.Generator | None = None,
    degeneracy: bool = False,
) -> TrialResult:
    """Run one forward model and the estimator on it."""
    rng = np.random.default_rng(0) if rng is None else rng
    a = np.radians(tilt_deg)
    depth = env.plate_depth()
    plate = Plate(
        x_m=x_m, depth_m=depth, tilt_rad=a, length_m=plate_m, drift_mps=drift_mps,
        taper="tukey" if extended else "uniform",
    )
    # Window: centred where the platform is on the normal, wide enough for the
    # envelope of the plate or of a 2 km facet, whichever is wider.
    probe = _echo(env, Plate(x_m=0.0, depth_m=depth, length_m=1.0), np.array([0.0]))
    p_vertical = float(probe.group_path_m[0])
    mu_s = probe.mu_sat
    w_env = max(
        float(envelope_width_s(env.f_hz * mu_s, p_vertical, min(plate_m, 2.0e3) if extended else plate_m)),
        float(envelope_width_s(env.f_hz * mu_s, p_vertical, 2.0e3)),
    )
    half = 2.5 * w_env if half_window_s is None else float(half_window_s)
    t_centre = (x_m + normal_crossing_offset_m(env, a)) / PLATFORM_SPEED_MPS
    t = _window(env, plate, t_centre, half, rng)
    echo = _echo(env, plate, t)

    i_pk = int(np.argmax(np.abs(echo.field)))
    p_meas = float(echo.group_path_m[i_pk])
    est = estimate_route_one(
        t, echo.field, f_hz=env.f_hz, group_path_m=p_meas, mu_sat=echo.mu_sat,
        v_s_mps=PLATFORM_SPEED_MPS,
    )

    mu_f = _mu_at_depth(env.f_hz, depth) if env.plasma else 1.0
    bias = predicted_tilt_bias_mps(env, a)
    los_pred = drift_mps * np.cos(a) + bias  # the line-of-sight rate at the normal crossing
    if est.separable:
        drift_err = est.drift_mps - drift_mps
        # Along-track position: route one reads it as v_s t_pk.
        angle_err = float(np.degrees((est.along_track_m - x_m) / p_meas))
    else:
        drift_err = float("nan")
        angle_err = float("nan")

    spread = ref = float("nan")
    if degeneracy:
        hist = doppler_history_from_field(t, echo.field)
        slope = 2.0 * env.f_hz * PLATFORM_SPEED_MPS**2 / (C_LIGHT * p_meas)
        trial = np.linspace(-100.0, 100.0, 41)
        best = free_t0_degeneracy(hist, env.f_hz, echo.mu_sat, slope, trial)
        spread = float(best.max() - best.min())
        # Reference: misfit with t0 held at truth and drift 25 m/s wrong.
        c_err = 2.0 * env.f_hz * echo.mu_sat * 25.0 / C_LIGHT
        ref = float(c_err**2)

    if extended:
        verdict = "pass" if not est.separable else "fail: separated an extended plate"
    elif case.startswith("D"):
        verdict = "compact" if est.separable else "declined"
    elif not est.separable:
        verdict = "fail: declined a compact facet"
    else:
        ok_v = abs(drift_err) <= DRIFT_TOL_MPS
        ok_a = abs(angle_err) <= ANGLE_TOL_DEG or tilt_deg != 0.0
        if tilt_deg == 0.0:
            verdict = "pass" if (ok_v and ok_a) else "fail"
        else:
            verdict = "fail: tilt read as drift" if not ok_v else "pass"
        if degeneracy and not (spread < 1e-6 * ref):
            verdict += "; degeneracy not reproduced"

    return TrialResult(
        case=case, environment=env.name, f_mhz=env.f_hz / 1e6, plate_km=plate_m / 1e3,
        x_true_km=x_m / 1e3, tilt_deg=tilt_deg, drift_true_mps=drift_mps, mu_sat=echo.mu_sat,
        mu_facet=mu_f, depth_km=depth / 1e3, group_path_km=p_meas / 1e3,
        separable=est.separable, reason=est.reason, slope_ratio=est.slope_ratio,
        peak_contrast=est.peak_contrast, drift_est_mps=est.drift_mps, drift_err_mps=drift_err,
        predicted_bias_mps=bias if tilt_deg else 0.0, angle_err_deg=angle_err,
        los_rate_mps=est.los_rate_mps, los_rate_predicted_mps=float(los_pred),
        degeneracy_spread_hz2=spread, degeneracy_reference_hz2=ref, verdict=verdict,
    )


def _all_trials(seed: int) -> list[TrialResult]:
    rng = np.random.default_rng(seed)
    out: list[TrialResult] = []
    env = {e.name: e for e in ENVIRONMENTS}

    # A: untilted compact facet.
    for e in ENVIRONMENTS:
        for x in (0.0, 5.0e3):
            for v in (-50.0, 0.0, 25.0, 100.0):
                out.append(run_trial("A_untilted", e, x_m=x, drift_mps=v, rng=rng,
                                     degeneracy=(v == 25.0)))
    # B: tilted compact facet.
    for name in ("free_5MHz_500km", "free_3MHz_500km", "plasma_5MHz_mu0.5", "plasma_3MHz_mu0.5"):
        for tilt in (0.1, 0.25, 0.5, 1.0):
            for v in (0.0, 25.0):
                out.append(run_trial("B_tilted", env[name], tilt_deg=tilt, drift_mps=v, rng=rng))
    # C: extended plate.
    for name in ("free_5MHz_500km", "plasma_5MHz_mu0.5"):
        for tilt in (0.0, 0.5, 1.0):
            for v in (0.0, 25.0):
                out.append(run_trial("C_extended", env[name], plate_m=80e3, tilt_deg=tilt,
                                     drift_mps=v, extended=True, rng=rng))
    # D: plate length.
    e5 = env["free_5MHz_500km"]
    w2 = float(envelope_width_s(e5.f_hz, e5.depth_m, 2.0e3))
    for d_km in (0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0):
        half = max(2.5 * float(envelope_width_s(e5.f_hz, e5.depth_m, d_km * 1e3)), 2.5 * w2)
        out.append(run_trial("D_length", e5, plate_m=d_km * 1e3, drift_mps=25.0,
                             half_window_s=half, rng=rng))
    return out


def _figure(rows: list[TrialResult], out_dir: Path) -> list[Path]:
    try:
        import matplotlib  # noqa: F401
    except ImportError:  # pragma: no cover - optional dependency
        return []
    from tss.figures import style

    style.use()
    fig, axes = style.figure(1, 3, figsize=(11.0, 3.4))

    # (a) Case A drift error.
    ax = axes[0]
    a_rows = [r for r in rows if r.case == "A_untilted"]
    names = [e.name for e in ENVIRONMENTS]
    for r in a_rows:
        ax.plot(names.index(r.environment), abs(r.drift_err_mps), "o", color=style.SERIES[0], ms=3)
    ax.axhline(DRIFT_TOL_MPS, color=style.INK_MUTED, ls="--", lw=0.8)
    ax.set_yscale("log")
    ax.set_xticks(range(len(names)))
    short = {"free_3MHz_500km": "F 3/500", "free_5MHz_500km": "F 5/500",
             "free_10MHz_500km": "F 10/500", "free_5MHz_1000km": "F 5/1000",
             "plasma_3MHz_mu0.5": "P 3", "plasma_5MHz_mu0.5": "P 5", "plasma_10MHz_mu0.8": "P 10"}
    ax.set_xticklabels([short.get(n, n) for n in names], fontsize=7, rotation=45)
    ax.set_ylabel("|drift error| (m/s)")
    ax.set_title("(a) Untilted facet", fontsize=9)

    # (b) Case B returned drift against tilt.
    ax = axes[1]
    for k, name in enumerate(("free_5MHz_500km", "plasma_5MHz_mu0.5", "plasma_3MHz_mu0.5")):
        b = [r for r in rows if r.case == "B_tilted" and r.environment == name and r.drift_true_mps == 25.0]
        b.sort(key=lambda r: r.tilt_deg)
        tilt = [r.tilt_deg for r in b]
        ax.plot(tilt, [r.drift_est_mps for r in b], "o", color=style.SERIES[k], ms=4, label=name)
        ax.plot(tilt, [25.0 + r.predicted_bias_mps for r in b], "-", color=style.SERIES[k], lw=0.8)
    ax.axhline(25.0, color=style.INK, lw=0.8)
    for lim, lab in ((-50.0, "SCI-15 goal"), (-100.0, "SCI-15 threshold")):
        ax.axhline(25.0 + lim, color=style.INK_MUTED, ls=":", lw=0.8)
        ax.text(1.0, 25.0 + lim + 2.0, lab, fontsize=6, color=style.INK_MUTED, ha="right")
    ax.set_xlabel("facet normal tilt (deg)")
    ax.set_ylabel("route-one drift (m/s), truth 25")
    ax.set_title("(b) Tilted facet", fontsize=9)
    ax.legend(fontsize=6, frameon=False)

    # (c) Case D slope ratio against plate length.
    ax = axes[2]
    d = sorted((r for r in rows if r.case == "D_length"), key=lambda r: r.plate_km)
    ax.semilogx([r.plate_km for r in d], [r.slope_ratio for r in d], "o-", color=style.SERIES[0], ms=4)
    ax.axhspan(0.8, 1.2, color=style.SERIES[2], alpha=0.15, lw=0)
    ax.axvline(np.sqrt(C_LIGHT / 5e6 * 500e3 / 2.0) / 1e3, color=style.INK_MUTED, ls="--", lw=0.8)
    ax.set_xlabel("plate length (km)")
    ax.set_ylabel("fitted / ephemeris Doppler slope")
    ax.set_title("(c) Compact or extended, 5 MHz", fontsize=9)

    fig.tight_layout()
    return style.savefig(fig, out_dir / "c3_verify", formats=("png", "pdf"))


@register_study(
    "c3_verify",
    description="C3 verification: noise-free route-one and extended-plane tests (TV 8.15)",
    work_packages=["WP5"],
    claims=["C3"],
    status="implemented",
    requirements=["SCI-15", "SCI-23"],
)
def run(out_dir: Path, *, config: dict[str, Any], manifest: Any) -> StudyResult:
    opts = dict(config.get("c3_verify", {}))
    seed = int(opts.get("seed", 20261008))
    rows = _all_trials(seed)

    csv_path = out_dir / "c3_verify_trials.csv"
    fields = list(asdict(rows[0]).keys())
    with open(csv_path, "w", encoding="utf-8") as fh:
        fh.write(f"# pass criteria: drift within {DRIFT_TOL_MPS} m/s, angle within {ANGLE_TOL_DEG} deg\n")
        fh.write(",".join(fields) + "\n")
        for r in rows:
            vals = []
            for k in fields:
                v = getattr(r, k)
                vals.append(f'"{v}"' if isinstance(v, str) else (f"{v:.6g}" if isinstance(v, float) else str(v)))
            fh.write(",".join(vals) + "\n")

    def sel(case: str) -> list[TrialResult]:
        return [r for r in rows if r.case == case]

    a, b, c, d = sel("A_untilted"), sel("B_tilted"), sel("C_extended"), sel("D_length")
    b_ratio = [
        (r.drift_err_mps / r.predicted_bias_mps) for r in b if r.separable and r.predicted_bias_mps
    ]
    metrics = {
        "seed": seed,
        "n_trials": len(rows),
        "A_pass": sum(r.verdict == "pass" for r in a),
        "A_n": len(a),
        "A_max_abs_drift_err_mps": float(np.nanmax([abs(r.drift_err_mps) for r in a])),
        "A_max_abs_angle_err_deg": float(np.nanmax([abs(r.angle_err_deg) for r in a])),
        "A_degeneracy_max_spread_over_reference": float(
            np.nanmax([r.degeneracy_spread_hz2 / r.degeneracy_reference_hz2 for r in a
                       if np.isfinite(r.degeneracy_spread_hz2)])
        ),
        "B_pass": sum(r.verdict == "pass" for r in b),
        "B_n": len(b),
        "B_flagged": sum(not r.separable for r in b),
        "B_err_over_predicted_min": float(min(b_ratio)) if b_ratio else float("nan"),
        "B_err_over_predicted_max": float(max(b_ratio)) if b_ratio else float("nan"),
        "B_err_mps_at_0p1deg_free5": float(next(
            r.drift_err_mps for r in b if r.environment == "free_5MHz_500km" and r.tilt_deg == 0.1)),
        "B_err_mps_at_1deg_free5": float(next(
            r.drift_err_mps for r in b if r.environment == "free_5MHz_500km" and r.tilt_deg == 1.0)),
        "C_pass": sum(r.verdict == "pass" for r in c),
        "C_n": len(c),
        "C_max_abs_slope_ratio": float(max(abs(r.slope_ratio) for r in c)),
        "D": {f"{r.plate_km:g}km": [round(r.slope_ratio, 3), r.verdict] for r in d},
    }

    summary = f"""
Noise-free tests of TV Section 8.15 with an independent physical-optics forward
model (exact geometry; exact Snell rays in an exponential topside with
f_N,sat = {PROFILE.f_n_sat_hz / 1e6:.1f} MHz and H = {PROFILE.scale_height_m / 1e3:.0f} km).
Criteria fixed before the run: drift within {DRIFT_TOL_MPS} m/s, along-track
arrival angle within {ANGLE_TOL_DEG} deg, extended plate declined.

A  Untilted compact facet: {metrics["A_pass"]}/{metrics["A_n"]} pass.  Worst drift error
   {metrics["A_max_abs_drift_err_mps"]:.3g} m/s, worst angle error {metrics["A_max_abs_angle_err_deg"]:.2g} deg.
   With t0 free the best misfit is flat in drift to {metrics["A_degeneracy_max_spread_over_reference"]:.1e}
   of the misfit a 25 m/s error gives with t0 known: the degeneracy is reproduced.
B  Tilted compact facet: {metrics["B_pass"]}/{metrics["B_n"]} pass, {metrics["B_flagged"]} flagged.
   The estimator returns drift plus v_s mu_f sin(a)/mu_sat: the error is
   {metrics["B_err_over_predicted_min"]:.3f} to {metrics["B_err_over_predicted_max"]:.3f} of that prediction,
   {metrics["B_err_mps_at_0p1deg_free5"]:.1f} m/s at 0.1 deg and {metrics["B_err_mps_at_1deg_free5"]:.1f} m/s at 1 deg (5 MHz, 500 km).
   The power peaks where the platform crosses the facet normal, not at closest
   approach, so a single pass separates drift from the facet's position but not
   from its tilt.
C  Extended plate: {metrics["C_pass"]}/{metrics["C_n"]} declined (largest slope ratio
   {metrics["C_max_abs_slope_ratio"]:.3f}).
D  Plate length at 5 MHz, 500 km (slope ratio, decision): {metrics["D"]}.
""".strip()

    outputs = [csv_path, *_figure(rows, out_dir)]
    return StudyResult(outputs=outputs, summary=summary, metrics=metrics)
