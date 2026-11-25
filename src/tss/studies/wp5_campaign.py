"""Layer: studies
WP: WP5
Claims: C3
Status: implemented
Requirements: SCI-15, INS-15, INS-26

WP5 Phase 0: the evaluation harness for claim C3.

Every WP5 phase is scored against the same fixed list of scenes.  In Phase 0
each scene runs through the existing parameter-level code (the physical-optics
forward models and the route-one and field-aligned estimators), so this study
is a harness and its first run, ``data/reference/wp5_baseline``, is the
baseline: it must reproduce the verified results before anything is built on
it.  Phase 3 runs the same scenes through the signal-level chain and is scored
against the same file.

THE SCENES
----------
The scenes are the grids that have already been verified, tagged by echo
class, so that every scene has a published number to reproduce:

* ``A_untilted``  c3_verify case A: an untilted 2 km facet at nadir and 5 km
  along track, drift -50, 0, 25 and 100 m/s, in the seven c3_verify
  environments (free space 3, 5, 10 MHz at 500 km and 5 MHz at 1000 km;
  plasma 3 and 5 MHz at mu_f 0.5 and 10 MHz at mu_f 0.8).  56 scenes.
* ``B_tilted``    c3_verify case B: facet normals tilted 0.1, 0.25, 0.5, 1
  degree, drift 0 and 25 m/s, four environments.  32 scenes.
* ``C_extended``  c3_verify case C: an 80 km tapered plate tilted 0, 0.5, 1
  degree, drift 0 and 25 m/s, free space and plasma at 5 MHz.  12 scenes.
* ``D_length``    c3_verify case D: plate length 0.5 to 32 km at 5 MHz and
  500 km, 25 m/s.  7 scenes.
* ``S_sheet``     field_aligned_verify: a field-aligned sheet at a dip of 10
  degrees, track heading psi 0, 15 and 45 degrees from magnetic north,
  rotation phi about B 0 and 5 degrees, in free space; and psi 0 and 15 at
  phi 5 inside the plasma.  8 scenes.
* ``E_budget``    added in Phase 1 (2026-11-24): an untilted 2 km facet at
  nadir, 25 m/s, inside the reference plasma at an apparent range of a
  budget point, so the budget has a forward-model check.  2 scenes.

The sheet's forward model fixes 5 MHz and 500 km, so the sheets are at 5 MHz
only.  Plasma facets in groups A to C are placed by their refractive index
mu_f, so their ``P'`` is recorded rather than set.  The ``E_budget`` facets are
placed by apparent range instead: mu_f is solved from
``P' = int_0^d dz / mu(z)`` in the c3_verify profile (f_N,sat 1.5 MHz,
exponential, H = 200 km).  Only two of the six budget points can be reached.
A facet below the reflection level is seen at no more than 527 km at 3 MHz,
750 km at 5 MHz and 1034 km at 10 MHz, so 3 MHz and 5 MHz at 1000 km do not
exist in this profile, and 3 MHz at 500 km and 10 MHz at 1000 km fall at the
reflection level (mu_f 0.07 and 0.08), where the plate model refuses to run.
Those four are listed in :data:`BUDGET_POINTS_UNREACHABLE` and not added
(decision of 2026-11-24 evening).  The analytic budget below still covers all
six, because it treats ``P'`` and ``mu_sat`` as independent.

Each scene draws from its own generator, seeded from the campaign seed and a
checksum of the scene's id, so a scene gives the same row alone, in a subset or
in the full list, whatever the order (the verified studies thread one
generator through their whole grid instead).  The generator sets only the
sub-sample offset of the sweep grid; the forward models are deterministic.

THE BUDGET
----------
The route-one error budget of TV Table ``c3r1budget`` at 3, 5 and 10 MHz and
apparent ranges of 500 and 1000 km, without fading and with fading at
tau_c 0.5 s, recomputed with :func:`tss.studies.route_one.route_one_budget`
with the seeds and trial count of the ``route_one`` study (seed 20261006 +
1000 + 10 i + j, 1000 trials), so the reproduction is exact rather than within
Monte Carlo sampling error.

WHAT IS SCORED
--------------
Rows carry the names of the single-node output contract (``docs/drift-mode.md``).
Truth is the plasma's own velocity along the wave normal (the SCI-15 rule of
2026-11-24): for a facet the plate's vertical drift ``w``, since
``cos(theta_sat)`` differs from 1 by under 1.6e-4 at the tilts here; for a
sheet the line-of-sight velocity ``s . V`` of field_aligned_verify.  Products
that have no parameter-level logic yet (``tau_c_s``, ``facet_size_m``,
``drift_mf_mps``, ``doppler_spread_hz``, the field-aligned classification
rule) are listed in each row's ``not_evaluated`` and are not computed here.

The expected values in :data:`EXPECTED` are copied from the Technical Volume
and the project notes before the run.  A run that does not match them reports
the mismatch; the expected values are not adjusted to fit.

CONFIGURATION
-------------
``wp5_campaign`` block, all optional::

    seed: 20261124
    groups: [A_untilted, B_tilted, C_extended, D_length, S_sheet]
    sheets: true            # include the S_sheet scenes
    plasma_sheets: true     # include the two plasma sheets (about 70 s each)
    budget_cases: [[3.0, 500.0], ...]     # default: all six
    budget_trials: 1000
"""

from __future__ import annotations

import json
import math
import zlib
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Any

import numpy as np

from tss.studies import StudyResult, register_study

__all__ = [
    "Scene",
    "SCENES",
    "GROUPS",
    "EXPECTED",
    "scene_seed",
    "run_scene",
    "run_budget",
    "check_budget",
    "check_reproduction",
]

DEFAULT_SEED = 20261124
ROUTE_ONE_SEED = 20261006
ROUTE_ONE_TRIALS = 1000
GROUPS = ("A_untilted", "B_tilted", "C_extended", "D_length", "S_sheet", "E_budget")

#: Budget points a facet can reach in the reference profile, and those it cannot.
BUDGET_POINTS_REACHABLE = ((5.0, 500.0), (10.0, 500.0))
BUDGET_POINTS_UNREACHABLE = ((3.0, 500.0), (3.0, 1000.0), (5.0, 1000.0), (10.0, 1000.0))

#: Smallest facet refractive index accepted for a budget-point scene.  Below
#: it the facet is at the reflection level and the plate model refuses it.
MIN_MU_FACET = 0.2

#: Products of the contract with no parameter-level logic in Phase 0.
NOT_EVALUATED_FACET = "tau_c_s;facet_size_m;drift_mf_mps;doppler_spread_hz"
NOT_EVALUATED_SHEET = "echo_class;phi_error_flags;tau_c_s;doppler_spread_hz"

#: Values the baseline must reproduce, copied before the run.  Sources:
#: c3_verify, TV Table c3verify and note c3-verify-route-one-tilt-finding
#: (2026-10-08); TV Table c3r1budget and note route-one-study-and-ch8-plasma-values
#: (2026-10-06); TV Table c3sheet (field_aligned_verify, seed 20261009).
EXPECTED: dict[str, Any] = {
    "A_pass": 56,
    "A_n": 56,
    "B_pass": 0,
    "B_flagged": 0,
    "B_n": 32,
    "B_err_over_predicted": (0.99, 1.01),
    "C_declined": 12,
    "C_n": 12,
    "D_compact_km": (0.5, 1.0, 2.0),
    "D_declined_km": (4.0, 8.0, 16.0, 32.0),
    "budget_total_mps": {
        (3.0, 500.0): (12.2, 33.7),
        (3.0, 1000.0): (11.3, 30.9),
        (5.0, 500.0): (10.8, 23.1),
        (5.0, 1000.0): (10.4, 20.1),
        (10.0, 500.0): (10.6, 14.8),
        (10.0, 1000.0): (10.5, 14.5),
    },
    "budget_tolerance": 0.05,
    "sheet_given_max_abs_err_free_mps": 0.38,
    "sheet_given_max_abs_err_plasma_mps": 0.08,
    "sheet_assume_err_psi15_phi5_mps": -167.0,
    "sheet_assume_err_psi45_phi5_mps": -465.0,
    "sheet_assume_abs_err_meridian_dip10_phi5_mps": 5.0,
    "sheet_assume_err_plasma_psi15_mps": (-90.0, -87.0),
    "sheet_assume_abs_err_plasma_meridian_mps": 2.6,
    # E_budget: new points (2026-11-24), not in the TV.  Expected from the
    # route-one acceptance rule (drift within 1 m/s) and from a screening
    # trial of the same evening (5 and 10 MHz at 500 km: pass, 0.001 m/s).
    "E_pass": 2,
    "E_n": 2,
    "E_group_path_tol_km": 1.0,
}


@dataclass(frozen=True)
class Scene:
    """One scene of the campaign.  Facet fields are ignored for sheets and back."""

    scene_id: str
    group: str
    echo_class_truth: str
    medium: str
    f_mhz: float
    environment: str = ""
    x_m: float = 0.0
    tilt_deg: float = 0.0
    drift_mps: float = 0.0
    plate_m: float = 2.0e3
    half_window_s: float | None = None
    degeneracy: bool = False
    dip_deg: float = 0.0
    heading_deg: float = 0.0
    rotation_deg: float = 0.0
    mu_facet: float | None = None
    p_target_km: float | None = None


def _group_path_m(f_hz: float, depth_m: float) -> float:
    """Vertical apparent range to ``depth_m`` in the reference profile, m."""
    from scipy.integrate import quad

    from tss.studies.c3_verify import PROFILE

    def inv_mu(z: float) -> float:
        return 1.0 / math.sqrt(max(float(PROFILE.mu_squared(f_hz, np.array([z]))[0]), 1e-30))

    return float(quad(inv_mu, 0.0, float(depth_m), limit=200)[0])


def max_facet_group_path_m(f_hz: float) -> float:
    """Apparent range of the reflection level in the reference profile, m."""
    from tss.studies.c3_verify import PROFILE

    return _group_path_m(f_hz, PROFILE.depth_for_mu(f_hz, 1e-9) * (1.0 - 1e-12))


def budget_point_mu(f_hz: float, p_group_m: float) -> float | None:
    """Facet refractive index that puts a facet at apparent range ``p_group_m``.

    ``None`` when no facet below the reflection level is that far, or when
    the facet would sit below :data:`MIN_MU_FACET`.
    """
    from scipy.optimize import brentq

    from tss.studies.c3_verify import PROFILE

    if float(p_group_m) >= max_facet_group_path_m(f_hz):
        return None
    d_ref = PROFILE.depth_for_mu(f_hz, 1e-9) * (1.0 - 1e-12)
    depth = brentq(lambda d: _group_path_m(f_hz, d) - float(p_group_m), 0.0, d_ref, xtol=1e-3)
    mu = math.sqrt(float(PROFILE.mu_squared(f_hz, np.array([depth]))[0]))
    return mu if mu >= MIN_MU_FACET else None


def _build_scenes() -> tuple[Scene, ...]:
    from tss.studies.c3_verify import ENVIRONMENTS
    from tss.studies.route_one import envelope_width_s

    env = {e.name: e for e in ENVIRONMENTS}
    med = {e.name: ("plasma" if e.plasma else "free") for e in ENVIRONMENTS}
    out: list[Scene] = []
    for e in ENVIRONMENTS:
        for x in (0.0, 5.0e3):
            for v in (-50.0, 0.0, 25.0, 100.0):
                out.append(Scene(f"A_{e.name}_x{x / 1e3:g}_v{v:g}", "A_untilted", "compact_facet", med[e.name],
                                 e.f_hz / 1e6, e.name, x_m=x, drift_mps=v, degeneracy=(v == 25.0)))
    for name in ("free_5MHz_500km", "free_3MHz_500km", "plasma_5MHz_mu0.5", "plasma_3MHz_mu0.5"):
        for tilt in (0.1, 0.25, 0.5, 1.0):
            for v in (0.0, 25.0):
                out.append(Scene(f"B_{name}_t{tilt:g}_v{v:g}", "B_tilted", "compact_facet", med[name],
                                 env[name].f_hz / 1e6, name, tilt_deg=tilt, drift_mps=v))
    for name in ("free_5MHz_500km", "plasma_5MHz_mu0.5"):
        for tilt in (0.0, 0.5, 1.0):
            for v in (0.0, 25.0):
                out.append(Scene(f"C_{name}_t{tilt:g}_v{v:g}", "C_extended", "extended", med[name],
                                 env[name].f_hz / 1e6, name, tilt_deg=tilt, drift_mps=v, plate_m=80e3))
    e5 = env["free_5MHz_500km"]
    w2 = float(envelope_width_s(e5.f_hz, e5.depth_m, 2.0e3))
    for d_km in (0.5, 1.0, 2.0, 4.0, 8.0, 16.0, 32.0):
        half = max(2.5 * float(envelope_width_s(e5.f_hz, e5.depth_m, d_km * 1e3)), 2.5 * w2)
        truth = "compact_facet" if d_km in EXPECTED["D_compact_km"] else "extended"
        out.append(Scene(f"D_free_5MHz_500km_L{d_km:g}", "D_length", truth, "free", 5.0,
                         "free_5MHz_500km", drift_mps=25.0, plate_m=d_km * 1e3, half_window_s=half))
    for heading in (0.0, 15.0, 45.0):
        for rot in (0.0, 5.0):
            out.append(Scene(f"S_free_dip10_psi{heading:g}_phi{rot:g}", "S_sheet", "field_aligned", "free",
                             5.0, dip_deg=10.0, heading_deg=heading, rotation_deg=rot))
    for heading in (0.0, 15.0):
        out.append(Scene(f"S_plasma_dip10_psi{heading:g}_phi5", "S_sheet", "field_aligned", "plasma",
                         5.0, dip_deg=10.0, heading_deg=heading, rotation_deg=5.0))
    for f_mhz, p_km in BUDGET_POINTS_REACHABLE:
        mu = budget_point_mu(f_mhz * 1e6, p_km * 1e3)
        if mu is None:  # pragma: no cover - guarded by the tests
            raise RuntimeError(f"budget point {f_mhz} MHz {p_km} km is out of reach")
        name = f"plasma_{f_mhz:g}MHz_P{p_km:g}km"
        out.append(Scene(f"E_{name}_v25", "E_budget", "compact_facet", "plasma", f_mhz, name,
                         drift_mps=25.0, mu_facet=round(mu, 6), p_target_km=p_km))
    return tuple(out)


SCENES: tuple[Scene, ...] = _build_scenes()


def scene_seed(scene: Scene, campaign_seed: int) -> int:
    """The scene's seed: the campaign seed combined with a checksum of its id."""
    return int((int(campaign_seed) * 1_000_003 + zlib.crc32(scene.scene_id.encode())) % (2**63))


def _facet_row(scene: Scene, seed: int) -> dict[str, Any]:
    from tss.studies.c3_verify import ENVIRONMENTS, Environment, run_trial

    known = {e.name: e for e in ENVIRONMENTS}
    if scene.environment in known:
        env = known[scene.environment]
    else:
        env = Environment(scene.environment, scene.f_mhz * 1e6, scene.medium == "plasma",
                          mu_facet=float(scene.mu_facet))
    case = {"A_untilted": "A_untilted", "B_tilted": "B_tilted", "C_extended": "C_extended",
            "D_length": "D_length", "E_budget": "E_budget"}[scene.group]
    r = run_trial(
        case, env, plate_m=scene.plate_m, x_m=scene.x_m, tilt_deg=scene.tilt_deg,
        drift_mps=scene.drift_mps, extended=(scene.group == "C_extended"),
        half_window_s=scene.half_window_s, rng=np.random.default_rng(seed), degeneracy=scene.degeneracy,
    )
    flags = ["TILT_UNFLAGGED"] if r.separable else ["TILT_DEGENERATE", "DRIFT_NOT_SEPARABLE"]
    ratio = (r.degeneracy_spread_hz2 / r.degeneracy_reference_hz2
             if math.isfinite(r.degeneracy_spread_hz2) else float("nan"))
    return {
        "scene_id": scene.scene_id, "group": scene.group, "seed": int(seed),
        "echo_class_truth": scene.echo_class_truth,
        "echo_class": "compact_facet" if r.separable else "extended",
        "medium": scene.medium, "environment": scene.environment, "f_mhz": scene.f_mhz,
        "mu_sat": float(r.mu_sat), "mu_facet": float(r.mu_facet), "depth_km": float(r.depth_km),
        "group_path_km": float(r.group_path_km), "x_true_km": float(r.x_true_km),
        "tilt_deg": scene.tilt_deg, "plate_km": scene.plate_m / 1e3, "heading_deg": float("nan"),
        "rotation_deg": float("nan"),
        "truth_mps": float(scene.drift_mps),
        "drift_mps": float(r.drift_est_mps) if r.separable else float("nan"),
        "drift_err_mps": float(r.drift_err_mps),
        "drift_given_err_mps": float("nan"),
        "predicted_bias_mps": float(r.predicted_bias_mps),
        "angle_err_deg": float(r.angle_err_deg),
        "los_rate_mps": float(r.los_rate_mps),
        "los_rate_predicted_mps": float(r.los_rate_predicted_mps),
        "slope_ratio": float(r.slope_ratio), "peak_contrast": float(r.peak_contrast),
        "degeneracy_spread_over_reference": float(ratio),
        "reason": r.reason, "flags": ";".join(flags), "not_evaluated": NOT_EVALUATED_FACET,
        "verdict": r.verdict,
    }


def _sheet_row(scene: Scene, seed: int) -> dict[str, Any]:
    from tss.studies.field_aligned_verify import run_sheet_trial

    r = run_sheet_trial(scene.dip_deg, scene.heading_deg, scene.rotation_deg,
                        plasma=(scene.medium == "plasma"), rng=np.random.default_rng(seed))
    return {
        "scene_id": scene.scene_id, "group": scene.group, "seed": int(seed),
        "echo_class_truth": scene.echo_class_truth, "echo_class": "not_evaluated",
        "medium": scene.medium, "environment": r.medium, "f_mhz": scene.f_mhz,
        "mu_sat": float("nan"), "mu_facet": float("nan"), "depth_km": float("nan"),
        "group_path_km": float("nan"), "x_true_km": float("nan"),
        "tilt_deg": float("nan"), "plate_km": float("nan"), "heading_deg": scene.heading_deg,
        "rotation_deg": scene.rotation_deg,
        "truth_mps": float(r.truth_los_mps),
        "drift_mps": float(r.truth_los_mps + r.err_assume_mps),
        "drift_err_mps": float(r.err_assume_mps),
        "drift_given_err_mps": float(r.err_given_mps),
        "predicted_bias_mps": float(r.predicted_assume_err_mps),
        "angle_err_deg": float("nan"), "los_rate_mps": float("nan"), "los_rate_predicted_mps": float("nan"),
        "slope_ratio": float("nan"), "peak_contrast": float("nan"),
        "degeneracy_spread_over_reference": float("nan"),
        "reason": "phi assumed 0 (drift_fa_mps); phi given in drift_given_err_mps",
        "flags": "ALIAS_RISK" if r.assume_aliased else "",
        "not_evaluated": NOT_EVALUATED_SHEET, "verdict": "",
    }


def run_scene(scene: Scene, *, campaign_seed: int = DEFAULT_SEED) -> dict[str, Any]:
    """Run one scene at parameter level and return its row."""
    seed = scene_seed(scene, campaign_seed)
    if scene.group == "S_sheet":
        return _sheet_row(scene, seed)
    return _facet_row(scene, seed)


def budget_seed(index: int, fading: int) -> int:
    """The route_one study's seed for budget row ``index`` and fading case ``fading``."""
    return ROUTE_ONE_SEED + 1000 + 10 * int(index) + int(fading)


def run_budget(cases: list[tuple[float, float]], *, n_trials: int = ROUTE_ONE_TRIALS) -> list[dict[str, Any]]:
    """TV Table c3r1budget, recomputed with the route_one study's seeds."""
    from tss.studies.route_one import BUDGET_CASES, route_one_budget

    rows = []
    for f_mhz, p_km in cases:
        i = [tuple(map(float, c)) for c in BUDGET_CASES].index((float(f_mhz), float(p_km)))
        for j, tau in enumerate((None, 0.5)):
            b = route_one_budget(f_mhz * 1e6, p_km * 1e3, fading_tau_s=tau, n_trials=n_trials,
                                 seed=budget_seed(i, j))
            rows.append({"f_mhz": float(f_mhz), "p_km": float(p_km), "fading": j,
                         "tau_c_s": float("nan") if tau is None else tau, "seed": budget_seed(i, j),
                         "n_trials": int(n_trials), "mu_sat": b.mu_sat, "width_s": b.width_s,
                         "sigma_t0_s": b.sigma_t0_s, "overflight_mps": b.overflight_mps,
                         "doppler_mps": b.doppler_mps, "oscillator_mps": b.oscillator_mps,
                         "ephemeris_mps": b.ephemeris_mps, "range_scale_mps": b.range_scale_mps,
                         "total_mps": b.total_mps})
    return rows


def check_budget(rows: list[dict[str, Any]]) -> list[dict[str, Any]]:
    """Each row against TV Table c3r1budget, within the stated 5 %."""
    exp = EXPECTED["budget_total_mps"]
    tol = float(EXPECTED["budget_tolerance"])
    out = []
    for r in rows:
        e = exp[(float(r["f_mhz"]), float(r["p_km"]))][int(r["fading"])]
        dev = float(r["total_mps"]) / e - 1.0
        out.append({**r, "expected_mps": e, "deviation": dev, "within_5pct": abs(dev) <= tol})
    return out


def _sheet_lookup(rows, scene_id):
    return next((r for r in rows if r["scene_id"] == scene_id), None)


def check_reproduction(rows: list[dict[str, Any]], budget: list[dict[str, Any]]) -> dict[str, Any]:
    """Compare the scene rows and the budget with :data:`EXPECTED`; list every mismatch."""
    e = EXPECTED
    g = {k: [r for r in rows if r["group"] == k] for k in GROUPS}
    m: dict[str, Any] = {}
    bad: list[str] = []

    if g["A_untilted"]:
        a = g["A_untilted"]
        m["A_pass"], m["A_n"] = sum(r["verdict"] == "pass" for r in a), len(a)
        m["A_max_abs_drift_err_mps"] = float(np.nanmax([abs(r["drift_err_mps"]) for r in a]))
        m["A_max_abs_angle_err_deg"] = float(np.nanmax([abs(r["angle_err_deg"]) for r in a]))
        deg = [r["degeneracy_spread_over_reference"] for r in a
               if math.isfinite(r["degeneracy_spread_over_reference"])]
        m["A_degeneracy_max_ratio"] = float(max(deg)) if deg else float("nan")
        if m["A_n"] == e["A_n"] and m["A_pass"] != e["A_pass"]:
            bad.append(f"A: {m['A_pass']}/{m['A_n']} pass, expected {e['A_pass']}")
        if m["A_n"] < e["A_n"] and m["A_pass"] != m["A_n"]:
            bad.append(f"A: {m['A_pass']}/{m['A_n']} pass")
    if g["B_tilted"]:
        b = g["B_tilted"]
        ratios = [r["drift_err_mps"] / r["predicted_bias_mps"] for r in b
                  if r["echo_class"] == "compact_facet" and r["predicted_bias_mps"]]
        m["B_pass"] = sum(r["verdict"] == "pass" for r in b)
        m["B_flagged"] = sum(r["echo_class"] != "compact_facet" for r in b)
        m["B_n"] = len(b)
        m["B_err_over_predicted_min"] = float(min(ratios)) if ratios else float("nan")
        m["B_err_over_predicted_max"] = float(max(ratios)) if ratios else float("nan")
        lo, hi = e["B_err_over_predicted"]
        if m["B_flagged"] != e["B_flagged"]:
            bad.append(f"B: {m['B_flagged']} flagged, expected {e['B_flagged']}")
        if m["B_pass"] != e["B_pass"]:
            bad.append(f"B: {m['B_pass']} pass, expected {e['B_pass']}")
        if ratios and not (lo <= m["B_err_over_predicted_min"] and m["B_err_over_predicted_max"] <= hi):
            bad.append(f"B: error/predicted {m['B_err_over_predicted_min']:.4f}-{m['B_err_over_predicted_max']:.4f}")
    if g["C_extended"]:
        c = g["C_extended"]
        m["C_declined"], m["C_n"] = sum(r["echo_class"] == "extended" for r in c), len(c)
        m["C_max_abs_slope_ratio"] = float(max(abs(r["slope_ratio"]) for r in c))
        if m["C_declined"] != m["C_n"]:
            bad.append(f"C: {m['C_declined']}/{m['C_n']} declined")
    if g["D_length"]:
        d = g["D_length"]
        m["D"] = {f"{r['plate_km']:g}km": [round(r["slope_ratio"], 3), r["echo_class"]] for r in d}
        for r in d:
            if r["echo_class"] != r["echo_class_truth"]:
                bad.append(f"D: {r['plate_km']:g} km classed {r['echo_class']}")
    if g["E_budget"]:
        eb = g["E_budget"]
        target = {s.scene_id: s.p_target_km for s in SCENES if s.group == "E_budget"}
        m["E_pass"], m["E_n"] = sum(r["verdict"] == "pass" for r in eb), len(eb)
        m["E_max_abs_drift_err_mps"] = float(np.nanmax([abs(r["drift_err_mps"]) for r in eb]))
        m["E_max_abs_group_path_dev_km"] = float(max(abs(r["group_path_km"] - target[r["scene_id"]])
                                                     for r in eb))
        if m["E_pass"] != m["E_n"]:
            bad.append(f"E: {m['E_pass']}/{m['E_n']} pass")
        if m["E_max_abs_group_path_dev_km"] > e["E_group_path_tol_km"]:
            bad.append(f"E: apparent range off by {m['E_max_abs_group_path_dev_km']:.2f} km")
    if g["S_sheet"]:
        s = g["S_sheet"]
        free = [r for r in s if r["medium"] == "free"]
        pl = [r for r in s if r["medium"] == "plasma"]
        if free:
            m["sheet_given_max_abs_err_free_mps"] = float(max(abs(r["drift_given_err_mps"]) for r in free))
            if m["sheet_given_max_abs_err_free_mps"] > e["sheet_given_max_abs_err_free_mps"]:
                bad.append(f"S: phi given {m['sheet_given_max_abs_err_free_mps']:.3f} m/s in free space")
            for sid, key, tol in (
                ("S_free_dip10_psi15_phi5", "sheet_assume_err_psi15_phi5_mps", 1.0),
                ("S_free_dip10_psi45_phi5", "sheet_assume_err_psi45_phi5_mps", 1.0),
            ):
                r = _sheet_lookup(s, sid)
                if r is not None:
                    m[key] = r["drift_err_mps"]
                    if abs(r["drift_err_mps"] - e[key]) > tol:
                        bad.append(f"S: {sid} error {r['drift_err_mps']:.1f} m/s, expected {e[key]}")
            r = _sheet_lookup(s, "S_free_dip10_psi0_phi5")
            if r is not None:
                key = "sheet_assume_abs_err_meridian_dip10_phi5_mps"
                m[key] = abs(r["drift_err_mps"])
                if abs(m[key] - e[key]) > 0.4:
                    bad.append(f"S: meridian error {m[key]:.2f} m/s, expected {e[key]}")
        if pl:
            m["sheet_given_max_abs_err_plasma_mps"] = float(max(abs(r["drift_given_err_mps"]) for r in pl))
            if m["sheet_given_max_abs_err_plasma_mps"] > e["sheet_given_max_abs_err_plasma_mps"]:
                bad.append(f"S: phi given {m['sheet_given_max_abs_err_plasma_mps']:.3f} m/s in the plasma")
            r = _sheet_lookup(s, "S_plasma_dip10_psi15_phi5")
            if r is not None:
                lo, hi = e["sheet_assume_err_plasma_psi15_mps"]
                m["sheet_assume_err_plasma_psi15_mps"] = r["drift_err_mps"]
                if not lo <= r["drift_err_mps"] <= hi:
                    bad.append(f"S: plasma psi 15 error {r['drift_err_mps']:.1f} m/s, expected {lo} to {hi}")
            r = _sheet_lookup(s, "S_plasma_dip10_psi0_phi5")
            if r is not None:
                key = "sheet_assume_abs_err_plasma_meridian_mps"
                m[key] = abs(r["drift_err_mps"])
                if abs(m[key] - e[key]) > 0.3:
                    bad.append(f"S: plasma meridian error {m[key]:.2f} m/s, expected {e[key]}")
    if budget:
        chk = check_budget(budget)
        m["budget_max_abs_deviation"] = float(max(abs(r["deviation"]) for r in chk))
        for r in chk:
            if not r["within_5pct"]:
                bad.append(f"budget {r['f_mhz']:g} MHz {r['p_km']:g} km fading {r['fading']}: "
                           f"{r['total_mps']:.2f} m/s against {r['expected_mps']}")
    m["mismatches"] = bad
    m["reproduced"] = not bad
    return m


def _scenes_for(opts: dict[str, Any]) -> list[Scene]:
    groups = list(opts.get("groups") or GROUPS)
    unknown = [x for x in groups if x not in GROUPS]
    if unknown:
        raise ValueError(f"unknown group(s) {unknown}; known: {list(GROUPS)}")
    sheets = bool(opts.get("sheets", True))
    plasma_sheets = bool(opts.get("plasma_sheets", True))
    out = []
    for s in SCENES:
        if s.group not in groups:
            continue
        if s.group == "S_sheet" and (not sheets or (s.medium == "plasma" and not plasma_sheets)):
            continue
        out.append(s)
    return out


def _fmt(v: Any) -> str:
    if isinstance(v, bool):
        return str(v)
    if isinstance(v, float):
        return "nan" if not math.isfinite(v) else f"{v:.6g}"
    if isinstance(v, str):
        return f'"{v}"'
    return str(v)


def _write_csv(path: Path, rows: list[dict[str, Any]]) -> Path:
    cols = list(rows[0]) if rows else []
    with path.open("w", encoding="utf-8") as fh:
        fh.write(",".join(cols) + "\n")
        for r in rows:
            fh.write(",".join(_fmt(r[c]) for c in cols) + "\n")
    return path


@register_study(
    "wp5_campaign",
    description="WP5 evaluation harness: C3 scenes by echo class and the route-one budget (Phase 0 baseline)",
    work_packages=["WP5"],
    claims=["C3"],
    status="implemented",
    requirements=["SCI-15", "INS-15", "INS-26"],
)
def run(out_dir: Path, *, config: dict[str, Any], manifest: Any) -> StudyResult:
    opts = dict(config.get("wp5_campaign", {}) or {})
    seed = int(opts.get("seed", DEFAULT_SEED))
    scenes = _scenes_for(opts)
    from tss.studies.route_one import BUDGET_CASES

    cases = [tuple(map(float, c)) for c in (opts.get("budget_cases") or BUDGET_CASES)]
    n_trials = int(opts.get("budget_trials", ROUTE_ONE_TRIALS))

    rows = [run_scene(s, campaign_seed=seed) for s in scenes]
    budget = run_budget(cases, n_trials=n_trials) if cases else []
    repro = check_reproduction(rows, budget)

    manifest.seeds.update({
        "campaign": seed,
        "scenes": {r["scene_id"]: r["seed"] for r in rows},
        "budget": {f"{r['f_mhz']:g}MHz_{r['p_km']:g}km_fading{r['fading']}": r["seed"] for r in budget},
        "budget_trials": n_trials,
    })

    out_dir = Path(out_dir)
    outputs = [_write_csv(out_dir / "wp5_campaign_scenes.csv", rows)]
    outputs.append(_write_csv(out_dir / "wp5_campaign_budget.csv", check_budget(budget)))
    jpath = out_dir / "wp5_campaign.json"
    jpath.write_text(json.dumps({
        "seed": seed, "n_scenes": len(rows), "groups": sorted({s.group for s in scenes}),
        "budget_cases": [list(c) for c in cases], "budget_trials": n_trials,
        "expected": {k: (list(v) if isinstance(v, tuple) else v) for k, v in EXPECTED.items()
                     if k != "budget_total_mps"} | {
            "budget_total_mps": {f"{f:g}MHz_{p:g}km": list(v) for (f, p), v in EXPECTED["budget_total_mps"].items()}},
        "reproduction": repro,
        "scenes": [asdict(s) for s in scenes],
    }, indent=2, default=str), encoding="utf-8")
    outputs.append(jpath)

    lines = [f"# wp5_campaign: WP5 Phase 0 baseline (parameter level), seed {seed}", "",
             f"{len(rows)} scenes; budget {len(budget)} rows at {n_trials} trials.",
             f"Reproduced: **{repro['reproduced']}**.", ""]
    for k in ("A_pass", "B_flagged", "C_declined", "E_pass"):
        if k in repro:
            n = repro[k.split("_")[0] + "_n"]
            lines.append(f"- {k}: {repro[k]}/{n}")
    for k in ("A_max_abs_drift_err_mps", "A_degeneracy_max_ratio", "B_err_over_predicted_min",
              "B_err_over_predicted_max", "C_max_abs_slope_ratio", "sheet_given_max_abs_err_free_mps",
              "sheet_given_max_abs_err_plasma_mps", "sheet_assume_err_psi15_phi5_mps",
              "sheet_assume_err_psi45_phi5_mps", "sheet_assume_abs_err_meridian_dip10_phi5_mps",
              "sheet_assume_err_plasma_psi15_mps", "sheet_assume_abs_err_plasma_meridian_mps",
              "E_max_abs_drift_err_mps", "E_max_abs_group_path_dev_km", "budget_max_abs_deviation"):
        if k in repro:
            lines.append(f"- {k}: {repro[k]:.4g}")
    if "D" in repro:
        lines.append(f"- D (slope ratio, class): {repro['D']}")
    if budget:
        lines += ["", "| f MHz | P' km | fading | total m/s | TV m/s | deviation |", "|---|---|---|---|---|---|"]
        for r in check_budget(budget):
            lines.append(f"| {r['f_mhz']:g} | {r['p_km']:g} | {'tau 0.5 s' if r['fading'] else 'none'} | "
                         f"{r['total_mps']:.2f} | {r['expected_mps']} | {100 * r['deviation']:+.1f} % |")
    if repro["mismatches"]:
        lines += ["", "Mismatches:"] + [f"- {x}" for x in repro["mismatches"]]
    metrics = {k: v for k, v in repro.items() if not isinstance(v, dict)} | {"n_scenes": len(rows), "seed": seed}
    return StudyResult(outputs=outputs, summary="\n".join(lines), metrics=metrics)
