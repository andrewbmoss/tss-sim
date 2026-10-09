"""Drift of a field-aligned sheet with its orientation from the field model.

Layer:        study
WP:           WP5
Claims:       C3
Status:       implemented
Requirements: SCI-15

Question
--------
Topside irregularities are field aligned (Dyson 1969), so the normal of a
reflecting sheet is perpendicular to B and the field model fixes its tilt along
B.  Does that let one pass separate drift from tilt, where route one could not
(:mod:`tss.studies.c3_verify`)?

Forward model
-------------
:mod:`tss.l1_propagation.sheet_echo`: a three-dimensional physical-optics sum
over a sheet 1.3 km across B and tens of km along it, curved with the dipole
field line, moving rigidly at ``V = V_n n + V_u u`` (perpendicular to B).  Free
space for the full matrix; exact Snell rays in a stratified plasma frozen to
the sheet for a subset.

Estimators
----------
All are given the ephemeris and the field line the echo comes from.  The sheet's
rotation ``phi`` about B is what each treats differently
(:mod:`tss.l3_signal.novel.field_aligned_estimator`):

``assume``  ``phi = 0`` (normal in the magnetic meridian).
``fit``     ``phi`` chosen to flatten the compensated Doppler history.
``given``   ``phi`` true, as an ideal cross-track arrival angle would give it.
``route1``  route one as it stands, with no field model.

Truth is the line-of-sight velocity the radar measures, ``s . V`` with ``s``
the ray direction at the platform, power weighted over the window.

Grid
----
Free space, 5 MHz, platform 500 km above the sheet: dip 0, 10, 30 deg; track
heading from magnetic north 0, 2, 5, 15, 45, 90 deg; rotation -5, 0, 2, 5, 10
deg; ``V_n = 25``, ``V_u = 30`` m/s.  Plasma subset: 5 MHz, ``f_N,sat`` 1.5 MHz,
``H`` 200 km, sheet where ``mu = 0.5``; dip 0 and 10, heading 0 and 15, rotation
0 and 5.

The arrival-angle accuracy route three would need is found from the
sensitivity of ``given`` to a 0.5 deg error in ``phi``.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np

from tss.l1_propagation import sheet_echo as se
from tss.l1_propagation.plate_echo import ExponentialTopside
from tss.l3_signal.novel import field_aligned_estimator as fe
from tss.l3_signal.novel.route_one_estimator import estimate_route_one
from tss.studies import StudyResult, register_study

__all__ = ["SheetTrial", "run_sheet_trial"]

F_HZ = 5.0e6
V_S = 7600.0
HEIGHT_M = 500e3
SHEET_ALT_M = 350e3
REPETITION_HZ = 45.0
PROFILE = ExponentialTopside(f_n_sat_hz=1.5e6, scale_height_m=200e3)
SCI15 = {"goal": 50.0, "threshold": 100.0}


@dataclass(frozen=True)
class SheetTrial:
    medium: str
    dip_deg: float
    heading_deg: float
    rotation_deg: float
    v_n_mps: float
    v_u_mps: float
    truth_los_mps: float
    err_assume_mps: float
    err_given_mps: float
    err_fit_mps: float
    fit_rotation_deg: float
    err_route1_mps: float
    route1_reason: str
    dlos_dphi_mps_per_deg: float
    aoa_for_goal_deg: float
    aoa_for_threshold_deg: float
    predicted_assume_err_mps: float
    assume_aliased: bool


def _sheet(dip, rot, vn, vu, length):
    r = se.EARTH_RADIUS_M + SHEET_ALT_M
    return se.FieldAlignedSheet(
        dip_rad=np.radians(dip), rotation_rad=np.radians(rot),
        curvature_radius_m=se.dipole_curvature_radius_m(np.radians(dip), r),
        length_m=length, velocity_n_mps=vn, velocity_u_mps=vu,
    )


def _predicted_assume_error(dip, heading, rot) -> float:
    """Free-space closed form: v_s (t . n(phi) - t . n(0)), flat sheet."""
    t = se.heading_vector(np.radians(heading))
    n_true = se.sheet_frame(_sheet(dip, rot, 0, 0, 1e3))["n"]
    n_zero = se.sheet_frame(_sheet(dip, 0.0, 0, 0, 1e3))["n"]
    # The platform term is -v_s t.n (receding when t.n > 0); a wrong n leaves
    # its difference in the drift.
    return float(-V_S * (np.dot(t, n_true) - np.dot(t, n_zero)))


def run_sheet_trial(
    dip: float, heading: float, rot: float, *, v_n: float = 25.0, v_u: float = 30.0,
    plasma: bool = False, fit: bool = False, half_window_s: float | None = None,
    rng: np.random.Generator | None = None,
) -> SheetTrial:
    rng = np.random.default_rng(0) if rng is None else rng
    half = (1.0 if plasma else 2.0) if half_window_s is None else half_window_s
    th = se.heading_vector(np.radians(heading))
    b = se.sheet_frame(_sheet(dip, rot, 0, 0, 1e3))["b"]
    along = abs(float(np.dot(th, b))) * V_S * half
    length = max(30e3, 2.0 * (along + 12e3) / 0.7)
    sheet = _sheet(dip, rot, v_n, v_u, length)

    dt = 1.0 / REPETITION_HZ
    n = int(np.ceil(half / dt))
    t = float(rng.uniform(-0.5, 0.5)) * dt + np.arange(-n, n + 1) * dt

    if plasma:
        depth = PROFILE.depth_for_mu(F_HZ, 0.5)
        echo = se.echo_stratified(t, F_HZ, sheet, PROFILE, height_m=depth, heading_rad=np.radians(heading))
        p0 = se.platform_start_stratified(sheet, depth, PROFILE, F_HZ)
        prof, h = PROFILE, depth
    else:
        echo = se.echo_free_space(t, F_HZ, sheet, height_m=HEIGHT_M, heading_rad=np.radians(heading))
        p0 = se.platform_start(sheet, HEIGHT_M)
        prof, h = None, HEIGHT_M
    pos = p0 + V_S * t[:, None] * th
    vel = se.sheet_velocity(sheet)

    truth_rays = fe.predict_path(t, pos - vel * t[:, None], sheet, f_hz=F_HZ, profile=prof, height_m=h)
    pw = np.abs(echo.field) ** 2
    w = pw / pw.sum()
    truth = float(np.sum(w * (truth_rays.ray_at_platform @ vel)))

    def est(method, phi=0.0):
        return fe.estimate_field_aligned(
            t, echo.field, pos, sheet, f_hz=F_HZ, mu_sat=echo.mu_sat, method=method,
            rotation_rad=np.radians(phi), profile=prof, height_m=h,
        )

    e_assume = est("assume").los_mps - truth
    e_given = est("given", rot).los_mps - truth
    e_off = est("given", rot + 0.5).los_mps - truth
    slope = (e_off - e_given) / 0.5
    if fit:
        fres = est("fit")
        e_fit, phi_fit = fres.los_mps - truth, float(np.degrees(fres.rotation_rad))
    else:
        e_fit, phi_fit = float("nan"), float("nan")

    i = int(np.argmax(pw))
    r1 = estimate_route_one(t, echo.field, f_hz=F_HZ, group_path_m=float(echo.group_path_m[i]),
                            mu_sat=echo.mu_sat, v_s_mps=V_S)
    e_r1 = (r1.drift_mps - truth) if r1.separable else float("nan")

    def aoa(limit):
        room = limit - abs(e_given)
        return float(room / abs(slope)) if abs(slope) > 1e-6 and room > 0 else float("inf")

    return SheetTrial(
        medium="plasma_mu0.5" if plasma else "free", dip_deg=dip, heading_deg=heading,
        rotation_deg=rot, v_n_mps=v_n, v_u_mps=v_u, truth_los_mps=truth, err_assume_mps=e_assume,
        err_given_mps=e_given, err_fit_mps=e_fit, fit_rotation_deg=phi_fit, err_route1_mps=e_r1,
        route1_reason=r1.reason, dlos_dphi_mps_per_deg=slope, aoa_for_goal_deg=aoa(SCI15["goal"]),
        aoa_for_threshold_deg=aoa(SCI15["threshold"]),
        predicted_assume_err_mps=_predicted_assume_error(dip, heading, rot) if not plasma else float("nan"),
        assume_aliased=bool(
            2.0 * F_HZ * abs(_predicted_assume_error(dip, heading, rot)) / 299792458.0 > 0.8 * REPETITION_HZ / 2.0
        ),
    )


def _figure(rows: list[SheetTrial], out_dir: Path) -> list[Path]:
    try:
        import matplotlib  # noqa: F401
    except ImportError:  # pragma: no cover
        return []
    from tss.figures import style

    style.use()
    fig, axes = style.figure(1, 3, figsize=(11.5, 3.6))
    free = [r for r in rows if r.medium == "free"]
    headings = sorted({r.heading_deg for r in free})

    ax = axes[0]
    for k, rot in enumerate((2.0, 5.0, 10.0)):
        r = sorted((x for x in free if x.dip_deg == 10.0 and x.rotation_deg == rot), key=lambda x: x.heading_deg)
        ax.semilogy([x.heading_deg for x in r], [max(abs(x.predicted_assume_err_mps), 1e-3) for x in r], "-",
                    lw=0.8, color=style.SERIES[k])
        ok = [x for x in r if not x.assume_aliased]
        ax.semilogy([x.heading_deg for x in ok], [max(abs(x.err_assume_mps), 1e-3) for x in ok], "o",
                    ms=3, color=style.SERIES[k], label=f"phi = {rot:g} deg, phi assumed 0")
    r = sorted((x for x in free if x.dip_deg == 10.0 and x.rotation_deg == 5.0), key=lambda x: x.heading_deg)
    ax.semilogy([x.heading_deg for x in r], [max(abs(x.err_given_mps), 1e-3) for x in r], "s--", ms=3,
                color=style.INK, label="phi = 5 deg, phi given")
    for lim in SCI15.values():
        ax.axhline(lim, color=style.INK_MUTED, ls=":", lw=0.8)
    ax.set_xlabel("track heading from magnetic north (deg)")
    ax.set_ylabel("|drift error| (m/s)")
    ax.set_title("(a) Dip 10 deg, free space", fontsize=9)
    ax.legend(fontsize=6, frameon=False)

    ax = axes[1]
    for k, dip in enumerate((0.0, 10.0, 30.0)):
        r = sorted((x for x in free if x.dip_deg == dip and x.rotation_deg == 5.0), key=lambda x: x.heading_deg)
        ax.semilogy([x.heading_deg for x in r], [min(x.aoa_for_goal_deg, 30.0) for x in r], "o-", ms=3,
                    color=style.SERIES[k], label=f"dip {dip:g} deg")
    ax.set_xlabel("track heading from magnetic north (deg)")
    ax.set_ylabel("cross-track AoA accuracy for SCI-15 goal (deg)")
    ax.set_title("(b) What route three must supply", fontsize=9)
    ax.legend(fontsize=6, frameon=False)
    del headings

    ax = axes[2]
    for k, dip in enumerate((0.0, 10.0, 30.0)):
        r = sorted((x for x in free if x.dip_deg == dip and x.heading_deg == 0.0), key=lambda x: x.rotation_deg)
        ax.plot([x.rotation_deg for x in r], [x.err_assume_mps for x in r], "o-", ms=3,
                color=style.SERIES[k], label=f"dip {dip:g} deg")
    ax.axhline(0.0, color=style.INK, lw=0.6)
    ax.set_xlabel("sheet rotation about B (deg)")
    ax.set_ylabel("drift error, phi assumed 0 (m/s)")
    ax.set_title("(c) Meridional track (heading 0)", fontsize=9)
    ax.legend(fontsize=6, frameon=False)

    fig.tight_layout()
    return style.savefig(fig, out_dir / "field_aligned_verify", formats=("png", "pdf"))


@register_study(
    "field_aligned_verify",
    description="Drift of a field-aligned sheet with orientation from the field model",
    work_packages=["WP5"],
    claims=["C3"],
    status="implemented",
    requirements=["SCI-15"],
)
def run(out_dir: Path, *, config: dict[str, Any], manifest: Any) -> StudyResult:
    opts = dict(config.get("field_aligned_verify", {}))
    seed = int(opts.get("seed", 20261009))
    include_plasma = bool(opts.get("plasma", True))
    rng = np.random.default_rng(seed)
    rows: list[SheetTrial] = []
    for dip in (0.0, 10.0, 30.0):
        for heading in (0.0, 2.0, 5.0, 15.0, 45.0, 90.0):
            for rot in (-5.0, 0.0, 2.0, 5.0, 10.0):
                fit = dip == 10.0 and rot == 5.0 and heading in (0.0, 15.0, 45.0)
                rows.append(run_sheet_trial(dip, heading, rot, fit=fit, rng=rng))
    if include_plasma:
        for dip in (0.0, 10.0):
            for heading in (0.0, 15.0):
                for rot in (0.0, 5.0):
                    rows.append(run_sheet_trial(dip, heading, rot, plasma=True, rng=rng))

    path = out_dir / "field_aligned_verify.csv"
    fields = list(asdict(rows[0]).keys())
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(",".join(fields) + "\n")
        for r in rows:
            fh.write(",".join(f'"{v}"' if isinstance(v, str) else (f"{v:.6g}" if isinstance(v, float) else str(v))
                              for v in asdict(r).values()) + "\n")

    free = [r for r in rows if r.medium == "free"]
    pl = [r for r in rows if r.medium != "free"]
    meridional = [r for r in free if r.heading_deg == 0.0]
    fit_rows = [r for r in free if np.isfinite(r.err_fit_mps)]
    metrics = {
        "seed": seed,
        "n_free": len(free),
        "n_plasma": len(pl),
        "given_max_abs_err_free": float(max(abs(r.err_given_mps) for r in free)),
        "given_max_abs_err_plasma": float(max((abs(r.err_given_mps) for r in pl), default=float("nan"))),
        "assume_max_abs_err_meridional": float(max(abs(r.err_assume_mps) for r in meridional)),
        "assume_vs_closed_form_max_dev": float(max(abs(r.err_assume_mps - r.predicted_assume_err_mps)
                                                   for r in free if not r.assume_aliased)),
        "assume_aliased": sum(r.assume_aliased for r in free),
        "fit": [(r.heading_deg, round(r.fit_rotation_deg, 2), round(r.err_fit_mps, 1)) for r in fit_rows],
        "route1_separated": sum(np.isfinite(r.err_route1_mps) for r in rows),
        "aoa_goal_deg_dip10_rot5": {
            f"{r.heading_deg:g}": round(r.aoa_for_goal_deg, 3) for r in free
            if r.dip_deg == 10.0 and r.rotation_deg == 5.0
        },
        "plasma": [(r.dip_deg, r.heading_deg, r.rotation_deg, round(r.err_assume_mps, 2), round(r.err_given_mps, 2))
                   for r in pl],
    }
    summary = f"""
Field-aligned sheet (1.3 km across B, curved with the dipole field line), drift
V_n 25 m/s and V_u 30 m/s, 5 MHz.  {len(free)} free-space and {len(pl)} plasma trials.

With the rotation about B given (an ideal cross-track arrival angle) the drift is
recovered to {metrics["given_max_abs_err_free"]:.2f} m/s in free space and
{metrics["given_max_abs_err_plasma"]:.2f} m/s in the plasma.  With it assumed zero, the
error follows the closed form v_s (t.n(phi) - t.n(0)) to {metrics["assume_vs_closed_form_max_dev"]:.2f} m/s
where it does not alias ({metrics["assume_aliased"]} cases put the residual past the 45 Hz Nyquist limit):
first order in v_s sin(heading) sin(phi), and second order on a meridional track
(worst {metrics["assume_max_abs_err_meridional"]:.1f} m/s over dip 0-30 deg and phi -5 to 10 deg).
Fitting phi from the Doppler history does not work: {metrics["fit"]}.
Route one separated {metrics["route1_separated"]} of {len(rows)} echoes.
""".strip()
    outputs = [path, *_figure(rows, out_dir)]
    return StudyResult(outputs=outputs, summary=summary, metrics=metrics)
