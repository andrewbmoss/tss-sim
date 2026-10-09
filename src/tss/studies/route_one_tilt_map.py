"""Route one's drift error against facet tilt, frequency, depth and range.

Layer:        study
WP:           WP5
Claims:       C3
Status:       implemented
Requirements: SCI-15

The C3 verification study (:mod:`tss.studies.c3_verify`) showed that route one
reads the tilt of a compact facet as drift, with an error of

    b = v_s mu_f sin(a) / mu_sat = v_s sin(theta_sat),

where ``mu_f`` is the index at the facet, ``mu_sat`` at the platform and
``theta_sat`` the along-track arrival angle at the platform when the echo
power peaks.  This study maps that error over the instrument's operating space
and turns it into the facet tilt route one can tolerate against SCI-15.

What depends on what
--------------------
*  The systematic term depends on frequency and depth only through
   ``mu_f / mu_sat``.  Inside the plasma ``mu_f < mu_sat``, so a facet deep in
   the plasma, near its reflection level, gives a smaller error per degree of
   tilt than one near the platform.
*  Apparent range does not enter the systematic term.  It enters the random
   budget of :func:`tss.studies.route_one.route_one_budget`, through the width
   of the power envelope and the ephemeris term.
*  The tilt allowance combines the two.  The allowance is the largest tilt
   for which the systematic error and the random budget together stay within
   the SCI-15 accuracy, 50 m/s at goal and 100 m/s at threshold.  Two
   combinations are reported: root-sum-square, which is how the TV budget
   combines its terms, and the linear sum ``|b| + sigma``, which is the honest
   one for a bias that does not average down.

Spot checks
-----------
The map uses the closed form.  A subset of grid points is run through the
physical-optics forward model and the estimator to confirm it, including the
deepest facets and the highest frequencies.  Points where the Doppler history
would exceed the 22.5 Hz Nyquist limit of the 45 Hz drift mode are marked
aliased and not run: there the tilt cannot be read even in principle without
removing an unknown constant.

Profiles
--------
Two field-free exponential topsides: the reference case of TV Chapter 8,
``f_N,sat = 1.5 MHz``, and a low-density case, ``f_N,sat = 0.5 MHz``, both with
``H = 200 km``.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np

from tss.constants import C_LIGHT
from tss.l1_propagation.plate_echo import ExponentialTopside, Plate, echo_stratified
from tss.l3_signal.novel.route_one_estimator import estimate_route_one
from tss.studies import StudyResult, register_study
from tss.studies.route_one import PLATFORM_SPEED_MPS, envelope_width_s, route_one_budget

__all__ = ["vertical_group_path_m", "tilt_bias_mps", "tilt_allowance_deg", "GridPoint"]

SCI15 = {"goal": 50.0, "threshold": 100.0}
REPETITION_HZ = 45.0
PROFILES = {
    "ref_1p5MHz": ExponentialTopside(f_n_sat_hz=1.5e6, scale_height_m=200e3),
    "low_0p5MHz": ExponentialTopside(f_n_sat_hz=0.5e6, scale_height_m=200e3),
}
FREQUENCIES_MHZ = (3.0, 5.0, 7.5, 10.0, 15.0, 20.0)
MU_FACET = (0.95, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2)
TILTS_DEG = (0.1, 0.25, 0.5, 1.0, 2.0)


def vertical_group_path_m(profile: ExponentialTopside, f_hz: float, depth_m: float) -> float:
    """Apparent range of a facet straight below, ``int dz / mu``."""
    nodes, weights = np.polynomial.legendre.leggauss(96)
    d = 0.5 * depth_m * (1.0 + nodes)
    return float(np.sum(0.5 * depth_m * weights / np.sqrt(profile.mu_squared(f_hz, d))))


def tilt_bias_mps(tilt_deg: float, mu_f: float, mu_sat: float) -> float:
    """Magnitude of route one's drift error from facet tilt."""
    return float(PLATFORM_SPEED_MPS * mu_f * np.sin(np.radians(tilt_deg)) / mu_sat)


def tilt_allowance_deg(limit_mps: float, sigma_mps: float, mu_f: float, mu_sat: float,
                       *, linear: bool) -> float:
    """Largest tilt whose error, with ``sigma``, stays within ``limit``.  0 if none."""
    if not np.isfinite(sigma_mps):
        return float("nan")
    room = limit_mps - sigma_mps if linear else np.sqrt(max(limit_mps**2 - sigma_mps**2, 0.0))
    if room <= 0.0:
        return 0.0
    s = room * mu_sat / (PLATFORM_SPEED_MPS * mu_f)
    return float(np.degrees(np.arcsin(min(s, 1.0))))


def _max_doppler_hz(f_hz: float, mu_sat: float, tilt_deg: float, p_group_m: float,
                    facet_m: float = 2.0e3) -> float:
    """Largest |f_D| over the estimator's window: tilt term plus the sweep."""
    sin_sat = tilt_bias_mps(tilt_deg, 1.0, 1.0) / PLATFORM_SPEED_MPS  # upper bound, mu_f <= mu_sat
    w = float(envelope_width_s(f_hz * mu_sat, p_group_m, facet_m))
    slope = 2.0 * f_hz * PLATFORM_SPEED_MPS**2 / (C_LIGHT * p_group_m)
    return 2.0 * f_hz * mu_sat * PLATFORM_SPEED_MPS * sin_sat / C_LIGHT + slope * 2.5 * w


@dataclass(frozen=True)
class GridPoint:
    profile: str
    f_mhz: float
    mu_sat: float
    mu_facet: float
    ratio: float
    depth_km: float
    group_path_km: float
    sigma_none_mps: float
    sigma_fading_mps: float
    bias_per_0p1deg_mps: float
    allow_goal_rss_none_deg: float
    allow_goal_rss_fading_deg: float
    allow_goal_lin_none_deg: float
    allow_goal_lin_fading_deg: float
    allow_thr_rss_none_deg: float
    allow_thr_rss_fading_deg: float
    allow_thr_lin_none_deg: float
    allow_thr_lin_fading_deg: float
    nyquist_tilt_deg: float
    envelope_width_s: float
    envelope_resolved: bool


#: The power envelope must span at least this many sweeps at half power for its
#: peak to be timed at all.  Below it the fit has nothing to fit.
MIN_SWEEPS_PER_WIDTH = 4.0


@dataclass(frozen=True)
class SpotCheck:
    profile: str
    f_mhz: float
    mu_facet: float
    tilt_deg: float
    group_path_km: float
    predicted_err_mps: float
    simulated_err_mps: float
    separable: bool
    aliased: bool


def _nyquist_tilt_deg(f_hz: float, mu_sat: float, p_group_m: float) -> float:
    """Largest tilt for which the Doppler history stays inside +-22.5 Hz."""
    lo, hi = 0.0, 10.0
    if _max_doppler_hz(f_hz, mu_sat, 0.0, p_group_m) >= REPETITION_HZ / 2:
        return 0.0
    for _ in range(50):
        mid = 0.5 * (lo + hi)
        if _max_doppler_hz(f_hz, mu_sat, mid, p_group_m) < REPETITION_HZ / 2:
            lo = mid
        else:
            hi = mid
    return lo


def _grid(n_trials: int, seed: int) -> list[GridPoint]:
    rows: list[GridPoint] = []
    for pname, prof in PROFILES.items():
        for f_mhz in FREQUENCIES_MHZ:
            f = f_mhz * 1e6
            if f < 2.0 * prof.f_n_sat_hz:
                continue  # below the Doppler-mode floor, 2 f_N,sat
            mu_s = float(np.sqrt(prof.mu_squared(f, 0.0)))
            for k, mu_f in enumerate(MU_FACET):
                if mu_f >= mu_s:
                    continue
                depth = prof.depth_for_mu(f, mu_f)
                if depth < 10e3:
                    continue
                pg = vertical_group_path_m(prof, f, depth)
                width = float(envelope_width_s(f, pg, 2.0e3))
                resolved = width * REPETITION_HZ >= MIN_SWEEPS_PER_WIDTH
                sig = []
                for j, tau in enumerate((None, 0.5)):
                    if not resolved:
                        sig.append(float("nan"))
                        continue
                    b = route_one_budget(
                        f, pg, fading_tau_s=tau, f_n_sat_hz=prof.f_n_sat_hz, n_trials=n_trials,
                        seed=seed + 1000 * j + 10 * k + int(f_mhz),
                    )
                    sig.append(b.total_mps)
                allow = {}
                for lname, lim in SCI15.items():
                    for comb, lin in (("rss", False), ("lin", True)):
                        for case, s in (("none", sig[0]), ("fading", sig[1])):
                            allow[f"{lname}_{comb}_{case}"] = tilt_allowance_deg(lim, s, mu_f, mu_s, linear=lin)
                rows.append(GridPoint(
                    profile=pname, f_mhz=f_mhz, mu_sat=mu_s, mu_facet=mu_f, ratio=mu_f / mu_s,
                    depth_km=depth / 1e3, group_path_km=pg / 1e3, sigma_none_mps=sig[0],
                    sigma_fading_mps=sig[1], bias_per_0p1deg_mps=tilt_bias_mps(0.1, mu_f, mu_s),
                    allow_goal_rss_none_deg=allow["goal_rss_none"],
                    allow_goal_rss_fading_deg=allow["goal_rss_fading"],
                    allow_goal_lin_none_deg=allow["goal_lin_none"],
                    allow_goal_lin_fading_deg=allow["goal_lin_fading"],
                    allow_thr_rss_none_deg=allow["threshold_rss_none"],
                    allow_thr_rss_fading_deg=allow["threshold_rss_fading"],
                    allow_thr_lin_none_deg=allow["threshold_lin_none"],
                    allow_thr_lin_fading_deg=allow["threshold_lin_fading"],
                    nyquist_tilt_deg=_nyquist_tilt_deg(f, mu_s, pg),
                    envelope_width_s=width, envelope_resolved=resolved,
                ))
    return rows


def _spot_checks(seed: int) -> list[SpotCheck]:
    rng = np.random.default_rng(seed)
    out: list[SpotCheck] = []
    cases = [
        ("ref_1p5MHz", 3.0, mu_f, t) for mu_f in (0.8, 0.5, 0.25) for t in (0.25, 1.0, 2.0)
    ] + [
        ("ref_1p5MHz", 10.0, mu_f, t) for mu_f in (0.9, 0.5) for t in (0.25, 1.0)
    ] + [
        ("ref_1p5MHz", 20.0, 0.5, t) for t in (0.25, 0.5, 1.0)
    ] + [
        ("low_0p5MHz", 5.0, mu_f, t) for mu_f in (0.9, 0.3) for t in (0.5, 2.0)
    ]
    for pname, f_mhz, mu_f, tilt in cases:
        prof = PROFILES[pname]
        f = f_mhz * 1e6
        mu_s = float(np.sqrt(prof.mu_squared(f, 0.0)))
        depth = prof.depth_for_mu(f, mu_f)
        pg = vertical_group_path_m(prof, f, depth)
        pred = -tilt_bias_mps(tilt, mu_f, mu_s)
        aliased = _max_doppler_hz(f, mu_s, tilt, pg) >= REPETITION_HZ / 2
        if aliased:
            out.append(SpotCheck(pname, f_mhz, mu_f, tilt, pg / 1e3, pred, float("nan"), False, True))
            continue
        a = np.radians(tilt)
        # Platform on the facet normal: offset X(q) with q = mu_f sin(a).
        nodes, weights = np.polynomial.legendre.leggauss(64)
        d = 0.5 * depth * (1.0 + nodes)
        q = mu_f * np.sin(a)
        x_off = float(np.sum(0.5 * depth * weights * q / np.sqrt(prof.mu_squared(f, d) - q * q)))
        w = float(envelope_width_s(f * mu_s, pg, 2.0e3))
        dt = 1.0 / REPETITION_HZ
        n = int(np.ceil(2.5 * w / dt))
        t = x_off / PLATFORM_SPEED_MPS + float(rng.uniform(-0.5, 0.5)) * dt + np.arange(-n, n + 1) * dt
        plate = Plate(x_m=0.0, depth_m=depth, tilt_rad=a, length_m=2.0e3, drift_mps=25.0)
        echo = echo_stratified(t, f, plate, prof, v_s_mps=PLATFORM_SPEED_MPS)
        i = int(np.argmax(np.abs(echo.field)))
        est = estimate_route_one(t, echo.field, f_hz=f, group_path_m=float(echo.group_path_m[i]),
                                 mu_sat=echo.mu_sat, v_s_mps=PLATFORM_SPEED_MPS)
        out.append(SpotCheck(pname, f_mhz, mu_f, tilt, pg / 1e3, pred,
                             est.drift_mps - 25.0 if est.separable else float("nan"),
                             est.separable, False))
    return out


def _rng(v: list[float]) -> list[float]:
    a = np.array(v, dtype=float)
    return [float(np.nanmin(a)), float(np.nanmax(a))]


def _write_csv(path: Path, rows: list, header_note: str) -> None:
    fields = list(asdict(rows[0]).keys())
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(f"# {header_note}\n")
        fh.write(",".join(fields) + "\n")
        for r in rows:
            vals = []
            for k in fields:
                v = getattr(r, k)
                vals.append(f'"{v}"' if isinstance(v, str) else (f"{v:.6g}" if isinstance(v, float) else str(v)))
            fh.write(",".join(vals) + "\n")


def _figure(grid: list[GridPoint], spots: list[SpotCheck], out_dir: Path) -> list[Path]:
    try:
        import matplotlib  # noqa: F401
    except ImportError:  # pragma: no cover
        return []
    from tss.figures import style

    style.use()
    fig, axes = style.figure(1, 3, figsize=(11.5, 3.6))

    # (a) Error against tilt for several index ratios, with spot checks.
    ax = axes[0]
    tilt = np.linspace(0.0, 2.0, 101)
    for k, ratio in enumerate((1.0, 0.8, 0.5, 0.25)):
        ax.plot(tilt, [tilt_bias_mps(t, ratio, 1.0) for t in tilt], color=style.SERIES[k], lw=1.0,
                label=f"mu_f/mu_sat = {ratio:g}")
    ok = [s for s in spots if not s.aliased and s.separable]
    for s in ok:
        ax.plot(s.tilt_deg, abs(s.simulated_err_mps), "o", ms=3.5, mfc="none",
                color=style.INK, zorder=5)
    for lim, lab in ((50.0, "SCI-15 goal"), (100.0, "SCI-15 threshold")):
        ax.axhline(lim, color=style.INK_MUTED, ls=":", lw=0.8)
        ax.text(0.02, lim + 4, lab, fontsize=6, color=style.INK_MUTED)
    ax.set_xlabel("facet normal tilt (deg)")
    ax.set_ylabel("|drift error| (m/s)")
    ax.set_title("(a) Error per degree; circles: forward model", fontsize=9)
    ax.legend(fontsize=6, frameon=False, loc="upper left", bbox_to_anchor=(0.0, 0.92))

    # (b) Goal allowance (linear, fading) against frequency, by depth, reference profile.
    ax = axes[1]
    ref = [g for g in grid if g.profile == "ref_1p5MHz"]
    for k, mu_f in enumerate((0.8, 0.5, 0.2)):
        r = sorted((g for g in ref if g.mu_facet == mu_f), key=lambda g: g.f_mhz)
        if not r:
            continue
        ax.plot([g.f_mhz for g in r], [g.allow_goal_lin_none_deg for g in r], "o-", ms=3,
                color=style.SERIES[k], label=f"mu_f = {mu_f:g}, no fading")
        ax.plot([g.f_mhz for g in r], [g.allow_goal_lin_fading_deg for g in r], "s--", ms=3,
                color=style.SERIES[k], label=f"mu_f = {mu_f:g}, fading 0.5 s")
    ax.set_xlabel("frequency (MHz)")
    ax.set_ylabel("tolerable tilt, SCI-15 goal (deg)")
    ax.set_title("(b) Goal, linear sum, f_N,sat 1.5 MHz", fontsize=9)
    ax.legend(fontsize=5.5, frameon=False)

    # (c) Same at threshold.
    ax = axes[2]
    for k, mu_f in enumerate((0.8, 0.5, 0.2)):
        r = sorted((g for g in ref if g.mu_facet == mu_f), key=lambda g: g.f_mhz)
        if not r:
            continue
        ax.plot([g.f_mhz for g in r], [g.allow_thr_lin_none_deg for g in r], "o-", ms=3,
                color=style.SERIES[k], label=f"mu_f = {mu_f:g}")
        ax.plot([g.f_mhz for g in r], [g.allow_thr_lin_fading_deg for g in r], "s--", ms=3,
                color=style.SERIES[k])
        ax.plot([g.f_mhz for g in r], [g.nyquist_tilt_deg for g in r], ":", lw=0.8,
                color=style.SERIES[k])
    ax.set_xlabel("frequency (MHz)")
    ax.set_ylabel("tolerable tilt, SCI-15 threshold (deg)")
    ax.set_title("(c) Threshold; dotted: 45 Hz Nyquist limit", fontsize=9)
    ax.legend(fontsize=6, frameon=False)

    fig.tight_layout()
    return style.savefig(fig, out_dir / "route_one_tilt_map", formats=("png", "pdf"))


@register_study(
    "route_one_tilt_map",
    description="Route one drift error against facet tilt, frequency, depth and range",
    work_packages=["WP5"],
    claims=["C3"],
    status="implemented",
    requirements=["SCI-15"],
)
def run(out_dir: Path, *, config: dict[str, Any], manifest: Any) -> StudyResult:
    opts = dict(config.get("route_one_tilt_map", {}))
    seed = int(opts.get("seed", 20261009))
    n_trials = int(opts.get("n_trials", 300))
    grid = _grid(n_trials, seed)
    spots = _spot_checks(seed + 7)

    grid_path = out_dir / "route_one_tilt_map_grid.csv"
    spot_path = out_dir / "route_one_tilt_map_spotchecks.csv"
    _write_csv(grid_path, grid, f"allowances in deg; SCI-15 goal 50, threshold 100 m/s; {n_trials} trials")
    _write_csv(spot_path, spots, "drift 25 m/s, 2 km facet; errors in m/s")

    run_spots = [s for s in spots if not s.aliased and s.separable]
    rel = [abs(s.simulated_err_mps / s.predicted_err_mps - 1.0) for s in run_spots]
    ref = [g for g in grid if g.profile == "ref_1p5MHz"]
    metrics = {
        "seed": seed,
        "n_trials": n_trials,
        "n_grid": len(grid),
        "spot_checks_run": len(run_spots),
        "spot_checks_aliased": sum(s.aliased for s in spots),
        "spot_checks_declined": sum((not s.aliased) and (not s.separable) for s in spots),
        "spot_check_max_rel_dev": float(max(rel)) if rel else float("nan"),
        "ratio_range_ref": [float(min(g.ratio for g in ref)), float(max(g.ratio for g in ref))],
        "n_unresolved": sum(not g.envelope_resolved for g in grid),
        "goal_lin_none_range_deg": _rng([g.allow_goal_lin_none_deg for g in grid]),
        "goal_lin_fading_range_deg": _rng([g.allow_goal_lin_fading_deg for g in grid]),
        "thr_lin_none_range_deg": _rng([g.allow_thr_lin_none_deg for g in grid]),
        "thr_lin_fading_range_deg": _rng([g.allow_thr_lin_fading_deg for g in grid]),
        "free_space_goal_rss_deg_bound": tilt_allowance_deg(50.0, 0.0, 1.0, 1.0, linear=False),
        "free_space_threshold_deg_bound": tilt_allowance_deg(100.0, 0.0, 1.0, 1.0, linear=False),
    }
    summary = f"""
Route one's tilt error is v_s mu_f sin(a)/mu_sat, confirmed by the forward model at
{metrics["spot_checks_run"]} spot checks to within {100 * metrics["spot_check_max_rel_dev"]:.2f} %
({metrics["spot_checks_aliased"]} aliased at the 45 Hz drift-mode Nyquist limit,
{metrics["spot_checks_declined"]} declined).  Tolerable facet tilt over {len(grid)} grid points
(two profiles, 3 to 20 MHz, facet index 0.2 to 0.95; {metrics["n_unresolved"]} with the power
envelope narrower than {MIN_SWEEPS_PER_WIDTH:g} sweeps, where the peak cannot be timed), linear sum with the random budget:
goal {metrics["goal_lin_none_range_deg"][0]:.2f} to {metrics["goal_lin_none_range_deg"][1]:.2f} deg without fading,
{metrics["goal_lin_fading_range_deg"][0]:.2f} to {metrics["goal_lin_fading_range_deg"][1]:.2f} deg with fading at 0.5 s;
threshold {metrics["thr_lin_none_range_deg"][0]:.2f} to {metrics["thr_lin_none_range_deg"][1]:.2f} and
{metrics["thr_lin_fading_range_deg"][0]:.2f} to {metrics["thr_lin_fading_range_deg"][1]:.2f} deg.
""".strip()
    outputs = [grid_path, spot_path, *_figure(grid, spots, out_dir)]
    return StudyResult(outputs=outputs, summary=summary, metrics=metrics)
