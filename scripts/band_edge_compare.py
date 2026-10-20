#!/usr/bin/env python3
"""WP0 Phase 2: compare an epoch-driven solar-minimum band-edge run with Rev G.

    python scripts/band_edge_compare.py REFERENCE_RUN EPOCH_RUN [OUT_DIR]

REFERENCE_RUN is a ``tss study min_frequency`` output of the plan
``band_edge_reference`` (or ``data/reference/band_edge``); EPOCH_RUN is one of
``band_edge_solarmin_epoch``.  Writes ``phase2_comparison.md`` and
``phase2_comparison.json`` to OUT_DIR (default EPOCH_RUN): per cell and window
the ensemble band edge (minimum across topsides) at 95 and 99 per cent against
the Rev G solar-minimum value, the indices IRI actually ran on, the storm split
for the window run both ways, and the INS-01a decision rule of the WP0 plan
(>= 0.33 MHz stands; 0.264-0.33 within the -20 per cent tolerance; below 0.264
moves).
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

import numpy as np

THRESHOLD_MHZ, TOLERANCE_MHZ = 0.33, 0.264


def _rows(run: Path) -> list[dict]:
    with open(run / "minfreq_summary.csv", newline="") as fh:
        return list(csv.DictReader(fh))


def _ensemble(rows: list[dict], solar: str) -> dict[tuple[str, float], dict]:
    cells: dict[tuple[str, float], list[dict]] = {}
    for r in rows:
        if r["solar"] == solar:
            cells.setdefault((r["orbit_class"], float(r["altitude_km"])), []).append(r)
    out = {}
    for key, rs in cells.items():
        f95 = np.array([float(r["f_min_for_95pct_coverage_MHz"]) for r in rs])
        f99 = np.array([float(r["f_min_for_99pct_coverage_MHz"]) for r in rs])
        out[key] = {"f95": float(f95.min()), "f99": float(f99.min()),
                    "model95": rs[int(f95.argmin())]["topside"],
                    "per_model95": {r["topside"]: float(r["f_min_for_95pct_coverage_MHz"]) for r in rs}}
    return out


def _indices(run: Path, rows: list[dict], solar: str) -> dict:
    """Medians of the indices IRI reported, over every sample of a window."""
    want = ("rz12", "ig12", "f107", "f107_81", "ap_daily", "kp", "fof2_storm_ratio")
    acc: dict[str, list[np.ndarray]] = {k: [] for k in want}
    for r in rows:
        if r["solar"] != solar or not r.get("samples_file"):
            continue
        z = np.load(run / r["samples_file"])
        cols = [str(c) for c in z["extra_columns"]]
        ex = z["extra"]
        for k in want:
            if k in cols and ex.shape[1] > cols.index(k):
                acc[k].append(ex[:, cols.index(k)])
    out = {}
    for k, parts in acc.items():
        if parts:
            v = np.concatenate(parts)
            out[k] = {"p05": float(np.nanpercentile(v, 5)), "median": float(np.nanmedian(v)),
                      "p95": float(np.nanpercentile(v, 95))}
    return out


def verdict(f95: float) -> str:
    if f95 >= THRESHOLD_MHZ - 5e-4:
        return "stands"
    if f95 >= TOLERANCE_MHZ:
        return "within tolerance"
    return "moves (G1)"


def main(argv: list[str]) -> int:
    ref_run, epoch_run = Path(argv[1]), Path(argv[2])
    out = Path(argv[3]) if len(argv) > 3 else epoch_run
    ref = _ensemble(_rows(ref_run), "min")
    erows = _rows(epoch_run)
    windows = list(dict.fromkeys(r["solar"] for r in erows))
    res: dict = {"reference_run": str(ref_run), "epoch_run": str(epoch_run), "windows": {}}
    md = ["# WP0 Phase 2: epoch-driven solar minimum against Rev G", "",
          f"Reference: `{ref_run}` (user R12 = 10, storm off).  Epoch run: `{epoch_run}`.", "",
          "Ensemble band edge = minimum across the three topsides of the 5th (99 per cent: 1st) "
          "percentile of fN at the spacecraft, MHz.", ""]
    ref_gov = min(ref.items(), key=lambda kv: kv[1]["f95"])
    for w in windows:
        ens = _ensemble(erows, w)
        idx = _indices(epoch_run, erows, w)
        cells = {}
        for key, e in sorted(ens.items()):
            r = ref.get(key)
            cells[f"{key[0]}_{int(key[1])}"] = {
                "f95": e["f95"], "f99": e["f99"], "model95": e["model95"],
                "ref_f95": r["f95"] if r else None, "ref_f99": r["f99"] if r else None,
                "shift95_pct": 100.0 * (e["f95"] / r["f95"] - 1.0) if r else None,
                "shift99_pct": 100.0 * (e["f99"] / r["f99"] - 1.0) if r else None,
                "per_model95": e["per_model95"]}
        gov_key, gov = min(ens.items(), key=lambda kv: kv[1]["f95"])
        res["windows"][w] = {"indices": idx, "cells": cells,
                             "governing": {"cell": f"{gov_key[0]}_{int(gov_key[1])}",
                                           "f95": gov["f95"], "f99": gov["f99"],
                                           "model95": gov["model95"], "verdict": verdict(gov["f95"])}}
        md += [f"## {w}", ""]
        if idx:
            md.append("Indices IRI ran on (5th / median / 95th percentile over the samples): " + "; ".join(
                f"{k} {v['p05']:.1f} / {v['median']:.1f} / {v['p95']:.1f}" for k, v in idx.items()) + ".")
            md.append("")
        md += ["| cell | ensemble f95 | Rev G f95 | shift | ensemble f99 | Rev G f99 | shift | model |",
               "|---|---|---|---|---|---|---|---|"]
        for c, v in cells.items():
            md.append(f"| {c} | {v['f95']:.3f} | {v['ref_f95']:.3f} | {v['shift95_pct']:+.1f}% | "
                      f"{v['f99']:.3f} | {v['ref_f99']:.3f} | {v['shift99_pct']:+.1f}% | {v['model95']} |")
        md += ["", f"Governing: {gov_key[0]} {gov_key[1]:.0f} km, {gov['f95']:.3f} MHz at 95 per cent "
               f"({gov['f99']:.3f} at 99), {gov['model95']}: INS-01a **{verdict(gov['f95'])}**.", ""]
    if "min2008" in windows and "min2008q" in windows:
        a, b = res["windows"]["min2008"]["cells"], res["windows"]["min2008q"]["cells"]
        split = {c: 100.0 * (a[c]["f95"] / b[c]["f95"] - 1.0) for c in a if c in b}
        res["storm_split_2008_pct"] = split
        md += ["## Storm split, 2008-09", "",
               "Ensemble f95 with the foF2 storm model on, relative to every storm model off: " +
               ", ".join(f"{c} {v:+.1f}%" for c, v in split.items()) + ".", ""]
    allgov = min(((w, v["governing"]) for w, v in res["windows"].items()), key=lambda t: t[1]["f95"])
    res["overall"] = {"window": allgov[0], **allgov[1],
                      "rev_g_governing": {"cell": f"{ref_gov[0][0]}_{int(ref_gov[0][1])}",
                                          "f95": ref_gov[1]["f95"], "f99": ref_gov[1]["f99"]}}
    md += ["## Decision rule (WP0 plan, Phase 2)", "",
           f"Rev G governing value {ref_gov[1]['f95']:.3f} MHz ({ref_gov[0][0]} {ref_gov[0][1]:.0f} km).  "
           f"Lowest epoch-driven value {allgov[1]['f95']:.3f} MHz ({allgov[0]}, {allgov[1]['cell']}): "
           f"INS-01a **{allgov[1]['verdict']}** (>= 0.33 stands; 0.264-0.33 within tolerance; "
           "below 0.264 moves, gate G1).", ""]
    out.mkdir(parents=True, exist_ok=True)
    (out / "phase2_comparison.json").write_text(json.dumps(res, indent=1), encoding="utf-8")
    (out / "phase2_comparison.md").write_text("\n".join(md), encoding="utf-8")
    print("\n".join(md))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv))
