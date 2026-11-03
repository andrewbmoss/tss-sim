"""Layer: studies
WP: WP4, WP8
Claims: C1, C3
Status: implemented
Requirements: SCI-05, SCI-07, INS-10

WP4 Phase 5: the measured-ionosphere check.  Every usable ISIS-1 topside
profile is flown through the P4 chain, with clutter on and off, and the
retrieval is scored against the profile.

WHAT IS RUN
-----------
The ionosphere is the ISIS-1 profile above its lowest scaled height, completed
below by an alpha-Chapman peak (:class:`tss.studies.wp4_dispersion.MeasuredTopside`);
the spacecraft is put where ISIS-1 was (altitude, latitude, longitude, epoch),
so the field model and the local plasma frequency are the record's own.  The
instrument and chain are the campaign's (``leo_500km_single`` as the baseline
6U flies it, the P4 tracker), one seed per record.

Records: ``iquality`` in the chosen set, 4 to 9 by default (the archive's
index, 0 best to 10 worst, 4 to 10 in these files; Andrew, 2026-11-03), at least 6 points, spacecraft at 500-1100 km.

WHAT IS SCORED
--------------
Per record and clutter setting: SCI-05 and SCI-07 as the assessment tier scores
them (against the L0 truth, which below the lowest scaled height is the
Chapman cap), the coverage gap, INS-10, the trace measures of
:func:`tss.studies.wp4_campaign.trace_measures`, and two measures restricted to
the heights ISIS-1 actually measured: the median absolute relative N_e error
(``span_ne_err``) and the fraction of the measured span the retrieval reaches
(``span_reached``).  Aggregates are given overall, by stratum (latitude band,
day or night, solar activity, as in ``wp4_dispersion``) and by quality index.

Cost: a record takes several seconds per clutter setting, so the run is spread
over a process pool (``workers``, default the CPU count less two).  Reported,
not gated.
"""

from __future__ import annotations

import csv
import json
import math
import os
from concurrent.futures import ProcessPoolExecutor
from dataclasses import replace
from pathlib import Path
from typing import Any

import numpy as np

from tss.studies import StudyResult, register_study

DEFAULTS: dict[str, Any] = {"qualities": [4, 5, 6, 7, 8, 9], "clutter": ["on", "off"], "limit": None,
                            "workers": None, "configuration": "leo_500km_single",
                            "checkpoint": None}
_STATE: dict[str, Any] = {}


def select_records(qualities: list[int], limit: int | None = None) -> list[Any]:
    from tss.data import spdf_isis

    path = spdf_isis.local_ne_path()
    if not path:
        return []
    out = []
    qs = {int(q) for q in qualities}
    for p in spdf_isis.iter_profiles(path):
        alt = np.asarray(p.alt_km, dtype=float)
        if p.quality is None or int(p.quality) not in qs or alt.size < 6 or not 500.0 <= p.sat_alt_km <= 1100.0:
            continue
        out.append(p)
        if limit is not None and len(out) >= int(limit):
            break
    return out


def stratum(p: Any) -> str:
    ex = dict(getattr(p, "extra", {}) or {})
    lat = abs(float(p.lat_deg))
    band = "low" if lat < 20 else ("mid" if lat < 55 else "high")
    day = "day" if 6.0 <= float(ex.get("lt_hour", 12.0)) < 18.0 else "night"
    sun = "high_sun" if float(ex.get("rz12", 0.0)) >= 60 else "low_sun"
    return f"{band}/{day}/{sun}"


def span_measures(retrieved: Any, prof: Any) -> tuple[float, float]:
    """Median |N_ret/N_meas - 1| over the measured heights reached, and the fraction of the span reached."""
    if retrieved is None:
        return float("nan"), 0.0
    a_m, n_m = prof.ascending()
    ra = np.asarray(retrieved.alt_m, dtype=float)
    rn = np.asarray(retrieved.ne, dtype=float)
    ok = np.isfinite(ra) & np.isfinite(rn) & (rn > 0)
    if np.count_nonzero(ok) < 2:
        return float("nan"), 0.0
    ra, rn = ra[ok], rn[ok]
    order = np.argsort(ra)
    ra, rn = ra[order], rn[order]
    inside = (a_m >= ra[0]) & (a_m <= ra[-1])
    span = float(a_m[-1] - a_m[0])
    reached = float(min(ra[-1], a_m[-1]) - max(ra[0], a_m[0])) / span if span > 0 else 0.0
    if not np.any(inside):
        return float("nan"), max(reached, 0.0)
    est = np.exp(np.interp(a_m[inside], ra, np.log(rn)))
    return float(np.median(np.abs(est / n_m[inside] - 1.0))), max(min(reached, 1.0), 0.0)


def _init(config: dict[str, Any], qualities: list[int], limit: int | None) -> None:
    from tss.config import load_yaml, resolve_path
    from tss.pipeline import build_tiers
    from tss.studies import wp4_campaign as w

    _STATE["config"] = config
    _STATE["records"] = select_records(qualities, limit)
    tiers = dict(config["tiers"])
    _STATE["shared"] = build_tiers({k: v for k, v in tiers.items() if k not in ("l0", "l2")})
    _STATE["clutter"] = {"on": w._clutter_for(config, "on"), "off": None}
    opts = {**DEFAULTS, **dict(config.get("isis1_retrieval", {}) or {})}
    _STATE["configuration"] = str(opts["configuration"])
    _STATE["conf"] = load_yaml(resolve_path(_STATE["configuration"], subdir="configurations"))


def run_record(index: int, clutter: str) -> dict[str, Any]:
    """One record, one clutter setting, end to end (in a worker)."""
    import io
    import time

    from tss.pipeline import Pipeline
    from tss.studies import _grids, _platform_for, _sweep_from_config, build_instrument
    from tss.studies import wp4_campaign as w
    from tss.studies.wp4_dispersion import MeasuredTopside

    config, prof = _STATE["config"], _STATE["records"][index]
    row: dict[str, Any] = {"index": index, "clutter": clutter, "quality": int(prof.quality),
                           "epoch": prof.epoch.strftime("%Y-%m-%dT%H:%M:%SZ"), "lat_deg": float(prof.lat_deg),
                           "lon_deg": float(prof.lon_deg), "sat_alt_km": float(prof.sat_alt_km),
                           "stratum": stratum(prof), "n_points": int(len(prof.alt_km)),
                           "h_low_km": float(np.min(prof.alt_km))}
    t0 = time.perf_counter()
    try:
        l0 = MeasuredTopside(prof)
        ep = row["epoch"]
        conf = _STATE["conf"]
        pl = replace(_platform_for(conf, ep), alt_m=float(prof.sat_alt_km) * 1e3, lat_deg=float(prof.lat_deg),
                     lon_deg=float(prof.lon_deg), epoch=ep)
        sweep, _ = _sweep_from_config(config, io.StringIO())
        alt_m, freq_hz = _grids(config)
        tiers = {**_STATE["shared"], "l0": l0, "l2": build_instrument(config["tiers"]["l2"], config, conf)}
        kwargs: dict[str, Any] = {"seed": int(index) + 1, "keep_intermediates": True}
        if _STATE["clutter"][clutter] is not None:
            kwargs["clutter"] = _STATE["clutter"][clutter]
        res = Pipeline(tiers, **kwargs).run_configuration(_STATE["configuration"], platform=pl, sweep=sweep,
                                                          freq_hz=freq_hz, alt_m=alt_m, epoch=ep)
        req = res.assessment.requirements if res.assessment is not None else {}

        def ach(rid: str) -> float:
            return w._as_float((req.get(rid) or {}).get("achieved"))

        def ver(rid: str) -> str:
            return str((req.get(rid) or {}).get("verdict", "not_evaluated"))

        tm = w.trace_measures(res.ionogram, res.echoes)
        err, reach = span_measures(res.retrieved, prof)
        fo = float(res.truth.extra.get("fof2_hz", float("nan"))) if res.truth is not None else float("nan")
        row.update({
            "fof2_mhz": fo / 1e6, "fn_sat_mhz": float(prof.plasma_freq_sat_mhz),
            "sci05": ach("SCI-05"), "sci05_verdict": ver("SCI-05"),
            "sci07_m": ach("SCI-07"), "sci07_verdict": ver("SCI-07"),
            "ins10_db": ach("INS-10"), "ins10_verdict": ver("INS-10"),
            "coverage_gap_m": w._coverage_gap(res.retrieved, res.truth),
            "span_ne_err": err, "span_reached": reach,
            "completeness_o": tm.get("completeness_o", float("nan")),
            "labelling_error": tm.get("labelling_error", float("nan")),
            "n_o": tm.get("n_o", 0), "n_x": tm.get("n_x", 0), "n_false": tm.get("n_false", 0),
            "error": "",
        })
    except Exception as exc:  # noqa: BLE001  (a record that fails is reported, not dropped)
        row["error"] = f"{type(exc).__name__}: {str(exc)[:200]}"
    row["run_s"] = time.perf_counter() - t0
    return row


def _nan(f: Any, v: list[float]) -> float:
    a = np.asarray([x for x in v if x is not None], dtype=float)
    a = a[np.isfinite(a)]
    return float(f(a)) if a.size else float("nan")


def block(rs: list[dict[str, Any]]) -> dict[str, Any]:
    ok = [r for r in rs if not r.get("error")]
    out: dict[str, Any] = {"n": len(rs), "n_failed": len(rs) - len(ok)}
    for key, worst in (("sci05", np.max), ("sci07_m", np.max), ("coverage_gap_m", np.max),
                       ("span_ne_err", np.max), ("span_reached", np.min), ("ins10_db", np.min),
                       ("completeness_o", np.min), ("labelling_error", np.max)):
        vals = [r.get(key, float("nan")) for r in ok]
        out[f"{key}_median"] = _nan(np.median, vals)
        q = 10 if worst is np.min else 90
        out[f"{key}_p90"] = _nan(lambda a, q=q: np.percentile(a, q), vals)
        out[f"{key}_worst"] = _nan(worst, vals)
    for rid in ("sci05", "sci07", "ins10"):
        vs = [str(r.get(f"{rid}_verdict")) for r in ok]
        out[f"{rid}_verdicts"] = {v: vs.count(v) for v in sorted(set(vs))}
    out["run_s_median"] = _nan(np.median, [r["run_s"] for r in rs])
    return out


def aggregate(rows: list[dict[str, Any]]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for c in sorted({r["clutter"] for r in rows}):
        rc = [r for r in rows if r["clutter"] == c]
        out[c] = {"overall": block(rc),
                  "by_stratum": {s: block([r for r in rc if r["stratum"] == s])
                                 for s in sorted({r["stratum"] for r in rc})},
                  "by_quality": {str(q): block([r for r in rc if r["quality"] == q])
                                 for q in sorted({r["quality"] for r in rc})}}
    return out


def _load_checkpoint(path: str | None) -> dict[tuple[int, str], dict[str, Any]]:
    done: dict[tuple[int, str], dict[str, Any]] = {}
    if path and os.path.exists(path):
        with open(path, encoding="utf-8") as fh:
            for line in fh:
                try:
                    r = json.loads(line)
                except json.JSONDecodeError:   # a line cut short when the run was stopped
                    continue
                done[(int(r["index"]), str(r["clutter"]))] = r
    return done


def compute(config: dict[str, Any]) -> dict[str, Any]:
    """Run every record and clutter setting.

    ``checkpoint`` (a JSON-lines path) makes the run resumable: each finished
    row is appended as it arrives, and rows already there are not rerun.  The
    full run takes hours, and a stopped run otherwise loses everything.
    """
    import time
    from concurrent.futures import as_completed

    opts = {**DEFAULTS, **dict(config.get("isis1_retrieval", {}) or {})}
    quals = [int(q) for q in opts["qualities"]]
    limit = opts.get("limit")
    limit = int(limit) if limit not in (None, "", "none") else None
    clutters = [str(c) for c in opts["clutter"]]
    ckpt = opts.get("checkpoint") or None
    n = len(select_records(quals, limit))
    done = _load_checkpoint(ckpt)
    jobs = [(i, c) for i in range(n) for c in clutters if (i, c) not in done]
    workers = int(opts["workers"] or max((os.cpu_count() or 2) - 2, 1))
    t0 = time.perf_counter()
    fh = open(ckpt, "a", encoding="utf-8") if ckpt else None
    new_rows: list[dict[str, Any]] = []
    try:
        def keep(row: dict[str, Any]) -> None:
            new_rows.append(row)
            if fh is not None:
                fh.write(json.dumps(row) + "\n")
                fh.flush()

        if workers <= 1 or not jobs:
            if jobs:
                _init(config, quals, limit)
            for i, c in jobs:
                keep(run_record(i, c))
        else:
            with ProcessPoolExecutor(max_workers=workers, initializer=_init, initargs=(config, quals, limit)) as ex:
                for fut in as_completed([ex.submit(run_record, i, c) for i, c in jobs]):
                    keep(fut.result())
    finally:
        if fh is not None:
            fh.close()
    wall = time.perf_counter() - t0
    allrows = {(int(r["index"]), str(r["clutter"])): r for r in [*done.values(), *new_rows]}
    rows = [allrows[k] for k in sorted(allrows) if k[0] < n and k[1] in clutters]
    return {"options": {**opts, "workers": workers}, "n_records": n, "rows": rows, "aggregate": aggregate(rows),
            "wall_s": wall, "resumed_rows": len(done), "cpu_s": float(sum(r["run_s"] for r in rows))}


def _g(x: float, nd: int = 3) -> str:
    return "nan" if x is None or not math.isfinite(x) else f"{x:.{nd}f}"


def _summary(res: dict[str, Any]) -> tuple[str, dict[str, Any]]:
    agg = res["aggregate"]
    lines = [f"{res['n_records']} ISIS-1 records (iquality {res['options']['qualities']}), clutter "
             f"{', '.join(agg)}; {res['wall_s'] / 3600:.2f} h wall for this session ({res['resumed_rows']} rows resumed from the checkpoint), {res['cpu_s'] / 3600:.2f} CPU h in all, on "
             f"{res['options']['workers']} workers.", "",
             "| Clutter | N (failed) | SCI-05 median (p90, worst) | SCI-05 verdicts | SCI-07 km median (p90, worst) "
             "| SCI-07 verdicts | Span N_e error median (p90) | Span reached median (p10) | Coverage gap km median "
             "(p90) | O completeness median |",
             "|---|---|---|---|---|---|---|---|---|---|"]
    metrics: dict[str, Any] = {}
    for c, a in agg.items():
        o = a["overall"]
        lines.append(
            f"| {c} | {o['n']} ({o['n_failed']}) | {_g(o['sci05_median'])} ({_g(o['sci05_p90'])}, "
            f"{_g(o['sci05_worst'])}) | {o['sci05_verdicts']} | {_g(o['sci07_m_median'] / 1e3, 1)} "
            f"({_g(o['sci07_m_p90'] / 1e3, 1)}, {_g(o['sci07_m_worst'] / 1e3, 1)}) | {o['sci07_verdicts']} "
            f"| {_g(o['span_ne_err_median'])} ({_g(o['span_ne_err_p90'])}) | {_g(o['span_reached_median'], 2)} "
            f"({_g(o['span_reached_p90'], 2)}) | {_g(o['coverage_gap_m_median'] / 1e3, 0)} "
            f"({_g(o['coverage_gap_m_p90'] / 1e3, 0)}) | {_g(o['completeness_o_median'], 2)} |")
        metrics[c] = {k: o[k] for k in ("n", "n_failed", "sci05_median", "sci05_p90", "sci07_m_median",
                                        "sci07_m_p90", "span_ne_err_median", "span_reached_median",
                                        "coverage_gap_m_median", "sci05_verdicts", "sci07_verdicts")}
    for title, key in (("By stratum", "by_stratum"), ("By quality index", "by_quality")):
        lines += ["", f"{title}:", "", "| Clutter | Group | N | SCI-05 median (p90) | SCI-07 km median (p90) | "
                  "Span N_e error median | Coverage gap km median | SCI-05 fails |", "|---|---|---|---|---|---|---|---|"]
        for c, a in agg.items():
            for g, o in a[key].items():
                lines.append(f"| {c} | {g} | {o['n']} | {_g(o['sci05_median'])} ({_g(o['sci05_p90'])}) | "
                             f"{_g(o['sci07_m_median'] / 1e3, 1)} ({_g(o['sci07_m_p90'] / 1e3, 1)}) | "
                             f"{_g(o['span_ne_err_median'])} | {_g(o['coverage_gap_m_median'] / 1e3, 0)} | "
                             f"{o['sci05_verdicts'].get('fail', 0)} |")
    fails = [r for r in res["rows"] if r.get("error")]
    if fails:
        lines += ["", f"Failed runs: {len(fails)}; first: {fails[0]['error']}"]
    return "\n".join(lines) + "\n", metrics


@register_study(
    "isis1_retrieval",
    description="WP4 Phase 5: every usable ISIS-1 topside profile through the P4 chain, clutter on and off",
    work_packages=["WP4", "WP8"],
    claims=["C1", "C3"],
    status="implemented",
    requirements=["SCI-05", "SCI-07", "INS-10"],
)
def run(out_dir: Path, *, config: dict[str, Any], manifest: Any) -> StudyResult:
    res = compute(dict(config))
    path = out_dir / "isis1_retrieval_records.csv"
    cols = sorted({k for r in res["rows"] for k in r})
    with path.open("w", encoding="utf-8", newline="") as fh:
        wr = csv.DictWriter(fh, fieldnames=cols)
        wr.writeheader()
        for r in res["rows"]:
            wr.writerow(r)
    text, metrics = _summary(res)
    jpath = out_dir / "isis1_retrieval.json"
    jpath.write_text(json.dumps({k: v for k, v in res.items() if k != "rows"} | {"metrics": metrics}, indent=2,
                                default=str), encoding="utf-8")
    return StudyResult(outputs=[path, jpath], summary=text, metrics=metrics)
