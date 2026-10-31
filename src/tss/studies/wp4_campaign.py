"""Layer: studies
WP: WP4
Claims: C1, C2, C4
Status: implemented
Requirements: SCI-05, SCI-07, INS-05, INS-10

WP4 Phase 0: the evaluation harness of the baseline signal-processing chain.

Every WP4 phase is scored against the same campaign: a set of ionospheres
(latitude, local time, solar activity) run end to end through L0 to L5 over
several seeds.  The first run of this study, on the chain as WP3 left it, is the
WP4 baseline (``data/reference/wp4_baseline``); later phases are better or worse
than that, and the number says which.

WHAT IS RUN
-----------
Each profile is an alpha-Chapman ionosphere (:mod:`tss.l0_ionosphere.chapman`)
with stated foF2, hmF2 and plasma scale height, placed at a stated latitude and
local time.  The Chapman tier does not vary with location, so the place sets the
geomagnetic field (IGRF, which sets the O-X split) and the local time a later
climatology can read; the profile parameters carry the latitude, local time and
solar activity.  The values are representative, not fitted to any epoch:

* ``mid_day_max``    40 N, 12 LT, foF2 11 MHz, hmF2 320 km, Hp 120 km
* ``mid_day_min``    40 N, 12 LT, foF2 6.5 MHz, hmF2 260 km, Hp 80 km
* ``mid_night_min``  40 N, 00 LT, foF2 3.5 MHz, hmF2 300 km, Hp 150 km
  (Meehan 2017 measures 150 +/- 30 km at night)
* ``eq_evening_max``  5 N, 21 LT, foF2 12 MHz, hmF2 420 km, Hp 180 km
* ``high_night``     70 N, 00 LT, foF2 3 MHz, hmF2 300 km, Hp 110 km
* ``smoke``           0 N, 12 LT, the Chapman defaults (foF2 6.35 MHz,
  hmF2 300 km, Hp 100 km), which ties the campaign to the smoke numbers

The instrument is the campaign's (``config/instrument/baseline_6u.yaml``) as
the named configuration flies it, A1 stepped sweep, scalar receive.

WHAT IS SCORED
--------------
Per run: SCI-05 and SCI-07 (achieved and verdict), INS-10 (the minimum SNR at
the truth echoes, cutoff steps excluded per G5c), INS-05 (worst band), and the
trace measures below.  Since WP4 Phase 2 INS-10 is measured against the noise
floor and ``scr_min_db`` reports the same minimum against the column, which
clutter raises (signal to noise plus clutter).  Per profile and overall: medians and worst cases over
seeds.

Trace measures, per ionogram column that has an O truth echo (cutoff steps
excluded).  The extracted range is

* ``o``      nearer the O echo's virtual range than the X echo's, and within
  ``tol`` of it;
* ``x``      nearer the X echo's, and within ``tol`` of it;
* ``false``  finite but within ``tol`` of neither;
* ``missed`` not extracted.

``tol`` is twice the column's range resolution plus 15 km.  The 15 km is
there because the chain as WP3 left it reports ranges 5 to 15 km short of the
truth (its leading edge, measured by this study), and a tighter tolerance would
call that bias a false trace.  The bias itself is reported: the mean and RMS
error of the matched O and X samples.  Columns with no truth
echo at all that carry an extracted range are counted as ``spurious``.  The
completeness of the O trace is ``o / (columns with an O echo)``.

The clutter generator of Phase 1 is applied as the configuration says
(``clutter`` block) unless ``wp4_campaign.clutter`` is ``on`` or ``off``.

CONFIGURATION
-------------
``wp4_campaign`` block, all optional::

    seeds: [1, 2, 3, 4, 5]
    profiles: [mid_day_max, ...]      # default: all
    configuration: leo_500km_single
    clutter: config | on | off
"""

from __future__ import annotations

import copy
import io
import json
import math
from dataclasses import replace
from pathlib import Path
from typing import Any

import numpy as np

from tss.studies import StudyResult, register_study

C_LIGHT = 299_792_458.0
EPOCH = "2026-01-01T12:00:00Z"
DEFAULT_SEEDS = (1, 2, 3, 4, 5)
TOL_EXTRA_M = 15.0e3


def _nm_from_fo(fo_hz: float) -> float:
    """NmF2 (m^-3) from foF2 (Hz): f_N^2 = 80.6 N."""
    return float(fo_hz) ** 2 / 80.6


#: name: (lat_deg, local_time_h, foF2_Hz, hmF2_m, Hp_m, description)
PROFILES: dict[str, tuple[float, float, float, float, float, str]] = {
    "mid_day_max": (40.0, 12.0, 11.0e6, 320.0e3, 120.0e3, "mid-latitude noon, solar maximum"),
    "mid_day_min": (40.0, 12.0, 6.5e6, 260.0e3, 80.0e3, "mid-latitude noon, solar minimum"),
    "mid_night_min": (40.0, 0.0, 3.5e6, 300.0e3, 150.0e3, "mid-latitude midnight, solar minimum"),
    "eq_evening_max": (5.0, 21.0, 12.0e6, 420.0e3, 180.0e3, "equatorial evening, solar maximum"),
    "high_night": (70.0, 0.0, 3.0e6, 300.0e3, 110.0e3, "high-latitude midnight"),
    "smoke": (0.0, 12.0, None, 300.0e3, 100.0e3, "the smoke campaign's ionosphere (Chapman defaults)"),
}


def _lon_for_local_time(lt_h: float, epoch: str = EPOCH) -> float:
    ut = float(epoch[11:13]) + float(epoch[14:16]) / 60.0
    lon = (lt_h - ut) * 15.0
    return float((lon + 180.0) % 360.0 - 180.0) if lon != 180.0 else 180.0


def _l0_for(name: str) -> Any:
    from tss.registry import build

    _lat, _lt, fo, hm, hp, _d = PROFILES[name]
    kw: dict[str, Any] = {"hmf2_m": hm, "plasma_scale_height_m": hp}
    if fo is not None:
        kw["nmf2"] = _nm_from_fo(fo)
    return build("l0", "chapman", **kw)


def _clutter_for(config: dict[str, Any], option: str) -> Any:
    """The clutter generator the run uses, or ``None`` (Phase 1 hook)."""
    try:
        from tss.l1_propagation.clutter import ClutterGenerator
    except ImportError:  # before Phase 1
        return None
    block = copy.deepcopy(dict(config.get("clutter", {}) or {}))
    if option == "off":
        return None
    if option == "on":
        block["enabled"] = True
    return ClutterGenerator.from_config(block)


def _as_float(v: Any) -> float:
    try:
        x = float(v)
    except (TypeError, ValueError):
        return float("nan")
    return x


def trace_measures(ionogram: Any, echoes: Any) -> dict[str, Any]:
    """The O-trace completeness and the X, false and missed counts (module docstring)."""
    from tss.l5_assessment.scoring import truth_echo_cells

    cells = truth_echo_cells(echoes)
    f_c = np.asarray(cells["freq_hz"], dtype=float)
    r_c = np.asarray(cells["range_m"], dtype=float)
    m_c = np.asarray([str(m) for m in cells["mode"]], dtype=object)
    cut = np.asarray(cells["cutoff_step"], dtype=bool)
    f_i = np.asarray(ionogram.freq_hz, dtype=float)
    r_i = np.asarray(ionogram.virtual_range_m, dtype=float)
    r_i = r_i[:, 0] if r_i.ndim == 2 else r_i
    res = (ionogram.flags or {}).get("range_resolution_m_by_column")
    res = np.asarray(res, dtype=float) if res is not None and np.size(res) == f_i.size else np.full(
        f_i.size, _as_float((ionogram.flags or {}).get("range_resolution_m")))
    counts = {"o": 0, "x": 0, "false": 0, "missed": 0, "spurious": 0}
    err_o: list[float] = []
    err_x: list[float] = []
    n_truth_o = 0
    per_col: list[str] = []
    for j, f in enumerate(f_i):
        here = np.isclose(f_c, f, rtol=1e-6, atol=1.0)
        tol = 2.0 * res[j] + TOL_EXTRA_M if np.isfinite(res[j]) else 20.0e3 + TOL_EXTRA_M
        o = here & (m_c == "O") & ~cut
        x = here & (m_c == "X")
        r = r_i[j]
        if not np.any(o):
            if np.isfinite(r) and not np.any(here):
                counts["spurious"] += 1
                per_col.append("spurious")
            else:
                per_col.append("")
            continue
        n_truth_o += 1
        ro = float(r_c[o][0])
        d_o = abs(r - ro) if np.isfinite(r) else np.inf
        rx = float(r_c[x][np.argmin(np.abs(r - r_c[x]))]) if np.any(x) and np.isfinite(r) else np.nan
        d_x = abs(r - rx) if np.isfinite(rx) else np.inf
        if not np.isfinite(r):
            counts["missed"] += 1
            per_col.append("missed")
        elif d_o <= tol and d_o <= d_x:
            counts["o"] += 1
            err_o.append(r - ro)
            per_col.append("o")
        elif d_x <= tol:
            counts["x"] += 1
            err_x.append(r - rx)
            per_col.append("x")
        else:
            counts["false"] += 1
            per_col.append("false")
    e = np.asarray(err_o, dtype=float)
    ex = np.asarray(err_x, dtype=float)
    return {
        "n_truth_o_columns": n_truth_o,
        **{f"n_{k}": v for k, v in counts.items()},
        "completeness_o": counts["o"] / n_truth_o if n_truth_o else float("nan"),
        "fraction_on_x": counts["x"] / n_truth_o if n_truth_o else float("nan"),
        "fraction_false": counts["false"] / n_truth_o if n_truth_o else float("nan"),
        "range_error_o_rms_m": float(np.sqrt(np.mean(e**2))) if e.size else float("nan"),
        "range_error_o_mean_m": float(np.mean(e)) if e.size else float("nan"),
        "range_error_x_rms_m": float(np.sqrt(np.mean(ex**2))) if ex.size else float("nan"),
        "range_error_x_mean_m": float(np.mean(ex)) if ex.size else float("nan"),
        "per_column": per_col,
    }


def run_one(profile: str, seed: int, *, config: dict[str, Any], configuration: str,
            clutter: Any, shared: dict[str, Any], keep: bool = False) -> dict[str, Any]:
    """One profile, one seed, end to end.  Returns the scored row."""
    from tss.config import load_yaml, resolve_path
    from tss.pipeline import Pipeline
    from tss.studies import _grids, _platform_for, _sweep_from_config, build_instrument

    lat, lt, *_ = PROFILES[profile]
    conf = load_yaml(resolve_path(configuration, subdir="configurations"))
    sweep, _ = _sweep_from_config(config, io.StringIO())
    alt_m, freq_hz = _grids(config)
    platform = replace(_platform_for(conf, EPOCH), lat_deg=float(lat), lon_deg=_lon_for_local_time(lt))
    tiers = {**shared, "l0": _l0_for(profile),
             "l2": build_instrument(config["tiers"]["l2"], config, conf)}
    kwargs: dict[str, Any] = {"seed": int(seed), "keep_intermediates": True}
    if clutter is not None:
        kwargs["clutter"] = clutter
    pipe = Pipeline(tiers, **kwargs)
    res = pipe.run_configuration(configuration, platform=platform, sweep=sweep, freq_hz=freq_hz,
                                 alt_m=alt_m, epoch=EPOCH)
    a = res.assessment
    req = a.requirements if a is not None else {}

    def ach(rid: str) -> float:
        return _as_float((req.get(rid) or {}).get("achieved"))

    def ver(rid: str) -> str:
        return str((req.get(rid) or {}).get("verdict", "not_evaluated"))

    tm = trace_measures(res.ionogram, res.echoes)
    fo = float(res.truth.extra.get("fof2_hz", float("nan"))) if res.truth is not None else float("nan")
    row: dict[str, Any] = {
        "profile": profile, "seed": int(seed), "configuration": configuration,
        "lat_deg": lat, "local_time_h": lt, "fof2_mhz": fo / 1e6,
        "clutter": "on" if clutter is not None else "off",
        "sci05": ach("SCI-05"), "sci05_verdict": ver("SCI-05"),
        "sci07_m": ach("SCI-07"), "sci07_verdict": ver("SCI-07"),
        "ins10_db": ach("INS-10"), "ins10_verdict": ver("INS-10"),
        "scr_min_db": _scr_min(res.truth_echo_snr),
        "ins05_m": ach("INS-05"), "ins05_verdict": ver("INS-05"),
        "n_detected": int((res.ionogram.flags or {}).get("n_detected", 0)),
        **{k: v for k, v in tm.items() if k != "per_column"},
        "per_column": "".join({"o": "o", "x": "x", "false": "f", "missed": "m", "spurious": "s", "": "."}[c]
                              for c in tm["per_column"]),
        "clutter_rows": int(sum(1 for fl in (res.echoes.flags if res.echoes.flags is not None else [])
                                if str(fl).startswith("clutter"))),
        "run_s": float(sum(res.timing_s.values())),
    }
    if keep:
        row["_result"] = res
    return row


def _scr_min(tes: Any) -> float:
    """Worst SNR at the scored truth echoes against the column (noise plus clutter), dB."""
    if not tes or "snr_db_column" not in tes:
        return float("nan")
    col = np.asarray(tes["snr_db_column"], dtype=float)
    cut = np.asarray(tes.get("cutoff_step", np.zeros(col.size, dtype=bool)), dtype=bool)
    v = col[~cut] if cut.size == col.size else col
    return float(np.nanmin(v)) if np.any(np.isfinite(v)) else float("nan")


def _nanmed(v: list[float]) -> float:
    a = np.asarray(v, dtype=float)
    return float(np.nanmedian(a)) if np.any(np.isfinite(a)) else float("nan")


def _nanmax(v: list[float]) -> float:
    a = np.asarray(v, dtype=float)
    return float(np.nanmax(a)) if np.any(np.isfinite(a)) else float("nan")


def _nanmin(v: list[float]) -> float:
    a = np.asarray(v, dtype=float)
    return float(np.nanmin(a)) if np.any(np.isfinite(a)) else float("nan")


def aggregate(rows: list[dict[str, Any]]) -> dict[str, Any]:
    """Medians and worst cases per profile and overall."""
    def block(rs: list[dict[str, Any]]) -> dict[str, Any]:
        out = {
            "n_runs": len(rs),
            "sci05_median": _nanmed([r["sci05"] for r in rs]),
            "sci05_worst": _nanmax([r["sci05"] for r in rs]),
            "sci07_median_m": _nanmed([r["sci07_m"] for r in rs]),
            "sci07_worst_m": _nanmax([r["sci07_m"] for r in rs]),
            "ins10_median_db": _nanmed([r["ins10_db"] for r in rs]),
            "ins10_worst_db": _nanmin([r["ins10_db"] for r in rs]),
            "scr_min_median_db": _nanmed([r["scr_min_db"] for r in rs]),
            "completeness_o_median": _nanmed([r["completeness_o"] for r in rs]),
            "completeness_o_worst": _nanmin([r["completeness_o"] for r in rs]),
            "fraction_on_x_median": _nanmed([r["fraction_on_x"] for r in rs]),
            "fraction_false_median": _nanmed([r["fraction_false"] for r in rs]),
            "spurious_median": _nanmed([r["n_spurious"] for r in rs]),
            "range_error_o_rms_median_m": _nanmed([r["range_error_o_rms_m"] for r in rs]),
            "range_error_o_mean_median_m": _nanmed([r["range_error_o_mean_m"] for r in rs]),
            "range_error_x_mean_median_m": _nanmed([r["range_error_x_mean_m"] for r in rs]),
        }
        for rid in ("sci05", "sci07", "ins10"):
            vs = [r[f"{rid}_verdict"] for r in rs]
            out[f"{rid}_verdicts"] = {v: vs.count(v) for v in sorted(set(vs))}
        return out

    profiles = sorted({r["profile"] for r in rows}, key=list(PROFILES).index)
    return {"overall": block(rows), "by_profile": {p: block([r for r in rows if r["profile"] == p])
                                                   for p in profiles}}


def compute(config: dict[str, Any]) -> dict[str, Any]:
    from tss.pipeline import build_tiers

    opts = dict(config.get("wp4_campaign", {}) or {})
    seeds = [int(s) for s in (opts.get("seeds") or DEFAULT_SEEDS)]
    profiles = list(opts.get("profiles") or PROFILES)
    configuration = str(opts.get("configuration", "leo_500km_single"))
    option = str(opts.get("clutter", "config")).lower()
    clutter = _clutter_for(config, option)
    tiers = dict(config["tiers"])
    shared = build_tiers({k: v for k, v in tiers.items() if k not in ("l0", "l2")})
    rows = [run_one(p, s, config=config, configuration=configuration, clutter=clutter, shared=shared)
            for p in profiles for s in seeds]
    agg = aggregate(rows)
    return {
        "seeds": seeds, "profiles": profiles, "configuration": configuration,
        "clutter": "on" if clutter is not None else "off",
        "clutter_description": clutter.describe() if clutter is not None else None,
        "profile_definitions": {p: dict(zip(("lat_deg", "local_time_h", "fof2_hz", "hmf2_m", "hp_m",
                                              "description"), PROFILES[p], strict=True)) for p in profiles},
        "tolerance": "2 x column range resolution + 15 km",
        "rows": rows,
        "aggregate": agg,
    }


def _f(x: float, nd: int = 3) -> str:
    return "nan" if x is None or (isinstance(x, float) and not math.isfinite(x)) else f"{x:.{nd}f}"


def _summary(res: dict[str, Any]) -> tuple[str, dict[str, Any]]:
    agg = res["aggregate"]
    o = agg["overall"]
    lines = [
        f"# wp4_campaign: clutter {res['clutter']}",
        "",
        f"{len(res['profiles'])} profiles x {len(res['seeds'])} seeds on {res['configuration']}"
        f" ({o['n_runs']} runs).  Tolerance for a trace match: {res['tolerance']}.",
        "",
        "| Profile | SCI-05 med (worst) | SCI-07 km med (worst) | INS-10 dB med (worst) | S/(N+C) dB med"
        " | O complete med (worst) | on X | false | spurious | O bias km | X bias km |",
        "|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for p, b in list(agg["by_profile"].items()) + [("**all**", o)]:
        lines.append(
            f"| {p} | {_f(b['sci05_median'])} ({_f(b['sci05_worst'])}) | {_f(b['sci07_median_m'] / 1e3, 1)}"
            f" ({_f(b['sci07_worst_m'] / 1e3, 1)}) | {_f(b['ins10_median_db'], 1)} ({_f(b['ins10_worst_db'], 1)})"
            f" | {_f(b['scr_min_median_db'], 1)}"
            f" | {_f(b['completeness_o_median'], 2)} ({_f(b['completeness_o_worst'], 2)})"
            f" | {_f(b['fraction_on_x_median'], 2)} | {_f(b['fraction_false_median'], 2)}"
            f" | {_f(b['spurious_median'], 0)} | {_f(b['range_error_o_mean_median_m'] / 1e3, 1)}"
            f" | {_f(b['range_error_x_mean_median_m'] / 1e3, 1)} |")
    lines += ["", "Verdicts over all runs: "
              + "; ".join(f"{k.upper()} {o[k + '_verdicts']}" for k in ("sci05", "sci07", "ins10")), ""]
    metrics = {
        "clutter": res["clutter"],
        "n_runs": o["n_runs"],
        "sci05_median": o["sci05_median"], "sci05_worst": o["sci05_worst"],
        "sci07_median_m": o["sci07_median_m"], "sci07_worst_m": o["sci07_worst_m"],
        "ins10_median_db": o["ins10_median_db"], "ins10_worst_db": o["ins10_worst_db"],
        "completeness_o_median": o["completeness_o_median"], "completeness_o_worst": o["completeness_o_worst"],
        "fraction_on_x_median": o["fraction_on_x_median"],
    }
    return "\n".join(lines), metrics


@register_study(
    "wp4_campaign",
    description="WP4 evaluation harness: profiles x seeds end to end; SCI-05/07, INS-05/10 and O-trace completeness",
    work_packages=["WP4"],
    claims=["C1", "C2", "C4"],
    status="implemented",
    requirements=["SCI-05", "SCI-07", "INS-05", "INS-10"],
)
def run(out_dir: Path, *, config: dict[str, Any], manifest: Any) -> StudyResult:
    res = compute(dict(config))
    rows = res["rows"]
    cols = [c for c in rows[0] if not c.startswith("_")]
    path = out_dir / "wp4_campaign_runs.csv"
    with path.open("w", encoding="utf-8") as fh:
        fh.write(",".join(cols) + "\n")
        for r in rows:
            fh.write(",".join(_f(r[c], 6) if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")
    text, metrics = _summary(res)
    jpath = out_dir / "wp4_campaign.json"
    jpath.write_text(json.dumps({k: v for k, v in res.items() if k != "rows"} | {"metrics": metrics},
                                indent=2, default=str), encoding="utf-8")
    return StudyResult(outputs=[path, jpath], summary=text, metrics=metrics)
