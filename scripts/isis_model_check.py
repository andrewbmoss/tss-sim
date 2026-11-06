#!/usr/bin/env python3
"""WP4 Phase 6: where the data and the models fall short, on the ISIS-1 records.

    python scripts/isis_model_check.py run      OUT_DIR [--matlab matlab] [--limit N]
    python scripts/isis_model_check.py analyse  OUT_DIR

Decisions (Andrew, 2026-11-05): M2 includes measured profiles, both to score the
chain as close to truth as possible and to show where current data and models
are lacking, as the case for the dissertation and for follow-on missions.  This
script provides the second part, on the same records as the study
``isis1_retrieval`` (iquality 4-9, spacecraft at 500-1100 km, at least 6
points, in that study's order):

``run`` writes the records' times and places and evaluates PHaRLAP's IRI-2020
at each from 100 to 1200 km (``scripts/isis_model_profiles.m``) on the epoch's
own indices, with each of its three topsides (NeQuick, IRI-2001, IRI-2001
corrected) and the foF2 storm model on: the epoch configurations of WP0's
``topside_isis_epoch_iri.m``, over the whole profile.

``analyse`` writes ``summary.md``, ``summary.json``, ``model_records.csv`` and
``coverage.png``:

* **The models against ISIS-1.**  Each model profile is scored against the
  measured one over the measured span exactly as the chain's retrieval is in
  ``isis1_retrieval`` (rms fractional N_e error, SCI-05's metric, and rms height
  error at matched density, SCI-07's), with the same Rev G limits, and the
  plasma frequency at the spacecraft as a ratio in dB.  By plasma frequency at
  the spacecraft, stratum, solar activity and quality index.
* **The coverage of the data.**  Counts of every ISIS-1 record in the archive
  and of the usable set, by latitude, local time, solar activity and
  spacecraft height: where the historical data are thin or absent.
* **The records' own quality**, as the limit of the historical data.

A per-record comparison of the chain's error with the models' was not chosen
(2026-11-05); the two are reported on the same metric, separately.
"""

from __future__ import annotations

import argparse
import json
import math
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

PHARLAP = r"C:\Users\Andrew\PhD\Chartier Verification\pharlap_4.7.4"
QUALITIES = (4, 5, 6, 7, 8, 9)
SOLAR = (("solar_min", 0.0, 30.0), ("intermediate", 30.0, 100.0), ("solar_max", 100.0, 1e9))


def records(limit: int | None = None):
    from tss.studies.isis1_retrieval import select_records

    return select_records(list(QUALITIES), limit)


def solar_class(rz: float) -> str:
    for name, lo, hi in SOLAR:
        if lo <= rz < hi or (name == "solar_min" and rz <= 30.0):
            return name
    return "unknown"


def cmd_run(out: Path, matlab: str, limit: int | None) -> None:
    out.mkdir(parents=True, exist_ok=True)
    recs = records(limit)
    pts = np.array([[p.epoch.year, p.epoch.month, p.epoch.day, p.epoch.hour, p.epoch.minute,
                     float(p.lat_deg), (float(p.lon_deg) + 180.0) % 360.0 - 180.0] for p in recs], dtype=float)
    pcsv = out / "points.csv"
    np.savetxt(pcsv, pts, delimiter=",", fmt="%.6f")
    sd, p, o, root = (str(x).replace("\\", "/") for x in (ROOT / "scripts", pcsv, out / "iri_profiles.mat", PHARLAP))
    call = f"addpath('{sd}'); isis_model_profiles('{p}','{o}','{root}')"
    t0 = time.time()
    subprocess.run([matlab, "-batch", call], check=True)
    (out / "run.json").write_text(json.dumps({"n_records": len(recs), "matlab_s": time.time() - t0,
                                              "qualities": QUALITIES, "pharlap": PHARLAP}, indent=2))
    print(f"{len(recs)} records, IRI-2020 profiles in {time.time() - t0:.0f} s")


def _nan(f, v):
    a = np.asarray(v, dtype=float)
    a = a[np.isfinite(a)]
    return float(f(a)) if a.size else float("nan")


def model_rows(out: Path) -> list[dict]:
    from scipy.io import loadmat

    from tss.studies.isis1_retrieval import SCI05, SCI07_M, fn_band, span_scores, stratum, verdict

    mat = loadmat(out / "iri_profiles.mat")
    ne = np.asarray(mat["ne"], dtype=float)
    alt = np.asarray(mat["alt_km"], dtype=float).ravel() * 1e3
    names = [str(np.asarray(x).ravel()[0]) for x in np.asarray(mat["conf_names"]).ravel()]
    n = ne.shape[0]
    recs = records(n)
    rows = []
    from types import SimpleNamespace

    for k, p in enumerate(recs):
        a_m, n_m = p.ascending()
        ex = dict(p.extra or {})
        row = {"index": k, "quality": int(p.quality), "stratum": stratum(p), "lat_deg": float(p.lat_deg),
               "lt_hour": float(ex.get("lt_hour", float("nan"))), "rz12": float(ex.get("rz12", float("nan"))),
               "sat_alt_km": float(p.sat_alt_km), "h_low_km": float(a_m[0] / 1e3),
               "fn_sat_mhz": float(p.plasma_freq_sat_mhz), "fn_band": fn_band(float(p.plasma_freq_sat_mhz))}
        row["solar"] = solar_class(row["rz12"])
        # the measured profile as the truth, on a 5 km grid in log density
        g = np.arange(a_m[0], a_m[-1] + 1.0, 5.0e3)
        truth = SimpleNamespace(alt_m=g, ne=np.exp(np.interp(g, a_m, np.log(n_m))))
        for c, name in enumerate(names):
            v = ne[k, :, c]
            ok = np.isfinite(v) & (v > 0)
            if np.count_nonzero(ok) < 2:
                s05 = s07 = ratio_db = float("nan")
            else:
                s05, s07 = span_scores(SimpleNamespace(alt_m=alt[ok], ne=v[ok], flags={}), truth, a_m[0], a_m[-1])
                n_sat = math.exp(float(np.interp(a_m[-1], alt[ok], np.log(v[ok]))))
                ratio_db = 10.0 * math.log10(n_sat / float(n_m[-1]))       # = 20 log10 of the f_N ratio
            row[f"{name}_sci05"] = s05
            row[f"{name}_sci05_verdict"] = verdict(s05, SCI05)
            row[f"{name}_sci07_m"] = s07
            row[f"{name}_sci07_verdict"] = verdict(s07, SCI07_M)
            row[f"{name}_fnsat_db"] = ratio_db
        rows.append(row)
    return rows, names


def block(rs: list[dict], names: list[str]) -> dict:
    out = {"n": len(rs)}
    for name in names:
        s05 = [r[f"{name}_sci05"] for r in rs]
        s07 = [r[f"{name}_sci07_m"] for r in rs]
        db = [r[f"{name}_fnsat_db"] for r in rs]
        out[name] = {
            "sci05_median": _nan(np.median, s05), "sci05_p90": _nan(lambda a: np.percentile(a, 90), s05),
            "sci07_median_m": _nan(np.median, s07), "sci07_p90_m": _nan(lambda a: np.percentile(a, 90), s07),
            "sci05_at_threshold_or_better": float(np.mean([r[f"{name}_sci05_verdict"] in ("goal", "threshold")
                                                           for r in rs])) if rs else float("nan"),
            "sci07_at_threshold_or_better": float(np.mean([r[f"{name}_sci07_verdict"] in ("goal", "threshold")
                                                           for r in rs])) if rs else float("nan"),
            "fnsat_db_median": _nan(np.median, db),
            "fnsat_db_p05": _nan(lambda a: np.percentile(a, 5), db),
            "fnsat_db_p95": _nan(lambda a: np.percentile(a, 95), db),
        }
    return out


def coverage() -> dict:
    """Counts of every ISIS-1 record in the archive and of the usable set."""
    from tss.data import spdf_isis

    allp = list(spdf_isis.iter_profiles(spdf_isis.local_ne_path()))
    usable = records()
    lat_edges = np.arange(-90, 91, 15)
    lt_edges = np.arange(0, 25, 3)
    alt_edges = np.array([0, 500, 700, 900, 1100, 1400, 2000, 4000])

    def cube(ps):
        lat = np.array([float(p.lat_deg) for p in ps])
        lt = np.array([float(dict(p.extra or {}).get("lt_hour", np.nan)) for p in ps])
        rz = np.array([float(dict(p.extra or {}).get("rz12", np.nan)) for p in ps])
        alt = np.array([float(p.sat_alt_km) for p in ps])
        sol = np.array([solar_class(x) for x in rz])
        h_lat_lt = {s: np.histogram2d(lat[sol == s], lt[sol == s], bins=[lat_edges, lt_edges])[0].astype(int).tolist()
                    for s, _, _ in SOLAR}
        return {"n": len(ps),
                "by_solar": {s: int(np.count_nonzero(sol == s)) for s, _, _ in SOLAR},
                "by_alt": {f"{alt_edges[i]}-{alt_edges[i + 1]} km": int(np.count_nonzero((alt >= alt_edges[i])
                                                                                       & (alt < alt_edges[i + 1])))
                           for i in range(alt_edges.size - 1)},
                "by_quality": {str(q): int(np.count_nonzero(np.array([p.quality for p in ps]) == q))
                               for q in sorted({p.quality for p in ps if p.quality is not None})},
                "years": [int(min(p.epoch.year for p in ps)), int(max(p.epoch.year for p in ps))],
                "max_sat_alt_km": float(alt.max()),
                "lat_lt": h_lat_lt, "lat_edges": lat_edges.tolist(), "lt_edges": lt_edges.tolist()}

    return {"archive": cube(allp), "usable": cube(usable)}


def figure(cov: dict, path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    u = cov["usable"]
    fig, axs = plt.subplots(1, 3, figsize=(12, 4.2), sharey=True, constrained_layout=True)
    vmax = max(max(max(r) for r in u["lat_lt"][s]) for s, _, _ in SOLAR)
    for ax, (s, _, _) in zip(axs, SOLAR, strict=True):
        h = np.array(u["lat_lt"][s], dtype=float)
        h[h == 0] = np.nan
        im = ax.pcolormesh(u["lt_edges"], u["lat_edges"], h, cmap="viridis", vmin=1, vmax=vmax, shading="flat")
        ax.set_title(f"{s.replace('_', ' ')} ({u['by_solar'][s]} records)")
        ax.set_xlabel("local time (h)")
        ax.set_xticks(u["lt_edges"])
        ax.set_facecolor("#f0f0f0")
    axs[0].set_ylabel("latitude (deg)")
    fig.colorbar(im, ax=axs, label="usable ISIS-1 records (grey: none)")
    fig.suptitle("ISIS-1 usable records (quality 4-9, 500-1100 km): coverage by latitude and local time")
    fig.savefig(path, dpi=150)
    plt.close(fig)


def _g(x, nd=3):
    return "nan" if x is None or not math.isfinite(x) else f"{x:.{nd}f}"


def cmd_analyse(out: Path) -> None:
    rows, names = model_rows(out)
    groups = {"overall": {"all": rows},
              "by_fn_band": {b: [r for r in rows if r["fn_band"] == b] for b in sorted({r["fn_band"] for r in rows})},
              "by_solar": {s: [r for r in rows if r["solar"] == s] for s, _, _ in SOLAR},
              "by_stratum": {s: [r for r in rows if r["stratum"] == s] for s in sorted({r["stratum"] for r in rows})},
              "by_quality": {str(q): [r for r in rows if r["quality"] == q] for q in QUALITIES}}
    agg = {g: {k: block(v, names) for k, v in d.items() if v} for g, d in groups.items()}
    cov = coverage()
    figure(cov, out / "coverage.png")
    cols = sorted({k for r in rows for k in r})
    import csv

    with (out / "model_records.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=cols)
        w.writeheader()
        w.writerows(rows)
    (out / "summary.json").write_text(json.dumps({"models": names, "aggregate": agg, "coverage": cov}, indent=2),
                                      encoding="utf-8")
    L = ["# IRI-2020 against ISIS-1, and the coverage of the data (WP4 P6)", "",
         f"{len(rows)} usable ISIS-1 records (quality 4-9, 500-1100 km), the records of `isis1_retrieval`.  "
         "Each IRI-2020 profile (epoch indices, foF2 storm model on) is scored against the measured profile over "
         "the measured span on SCI-05's metric (rms fractional N_e) and SCI-07's (rms height at matched density), "
         "with the Rev G limits (0.15 / 0.10; 20 / 10 km).  The plasma frequency at the spacecraft is a ratio "
         "model/measured in dB (20 log10).", ""]
    for g, d in agg.items():
        L += [f"## {g.replace('_', ' ')}", "",
              "| Group | N | Topside | SCI-05 metric median (p90) | at threshold or better | SCI-07 metric km median "
              "(p90) | at threshold or better | f_N at spacecraft dB median [p5, p95] |",
              "|---|---|---|---|---|---|---|---|"]
        for k, b in d.items():
            for name in names:
                m = b[name]
                L.append(f"| {k} | {b['n']} | {name} | {_g(m['sci05_median'])} ({_g(m['sci05_p90'])}) | "
                         f"{m['sci05_at_threshold_or_better']:.0%} | {_g(m['sci07_median_m'] / 1e3, 1)} "
                         f"({_g(m['sci07_p90_m'] / 1e3, 1)}) | {m['sci07_at_threshold_or_better']:.0%} | "
                         f"{_g(m['fnsat_db_median'], 2)} [{_g(m['fnsat_db_p05'], 2)}, {_g(m['fnsat_db_p95'], 2)}] |")
        L.append("")
    for which in ("archive", "usable"):
        c = cov[which]
        L += [f"## Coverage: {which} ({c['n']} records, {c['years'][0]}-{c['years'][1]}, spacecraft to "
              f"{c['max_sat_alt_km']:.0f} km)", "",
              "| Solar activity | Records |", "|---|---|"] + [f"| {k} | {v} |" for k, v in c["by_solar"].items()]
        L += ["", "| Spacecraft height | Records |", "|---|---|"] + [f"| {k} | {v} |" for k, v in c["by_alt"].items()]
        L += ["", "| Quality index | Records |", "|---|---|"] + [f"| {k} | {v} |" for k, v in c["by_quality"].items()]
        empty = {s: sum(1 for row in c["lat_lt"][s] for v in row if v == 0) for s, _, _ in SOLAR}
        cells = (len(c["lat_edges"]) - 1) * (len(c["lt_edges"]) - 1)
        L += ["", "Empty latitude x local-time cells (15 deg x 3 h, of " + str(cells) + "): "
              + "; ".join(f"{k} {v}" for k, v in empty.items()), ""]
    L += ["`coverage.png`: the usable records by latitude and local time, per solar-activity class.", ""]
    (out / "summary.md").write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L))


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("out", type=Path)
    r.add_argument("--matlab", default="matlab")
    r.add_argument("--limit", type=int, default=None)
    a = sub.add_parser("analyse")
    a.add_argument("out", type=Path)
    args = ap.parse_args()
    if args.cmd == "run":
        cmd_run(args.out, args.matlab, args.limit)
    else:
        cmd_analyse(args.out)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
