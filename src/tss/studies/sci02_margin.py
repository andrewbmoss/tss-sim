"""SCI-02 threshold: the spacecraft's height above the F2 peak, from the band-edge samples.

Layer:        study
WP:           WP0
Claims:       C1
Status:       implemented
Requirements: SCI-02, MIS-01, SCI-08

SCI-02's threshold asks for continuous altitude coverage from the spacecraft
down to and including the F2 peak, for h_sat 500-900 km and hmF2 anywhere in
200-500 km.  MIS-01 adds that the spacecraft must be above the F2 peak at all
times sampled, and SCI-08 (hence INS-04's 1050 km) assumes hmF2 >= 200 km.
Three things are therefore checked, per orbit class, altitude, solar level,
local time and geomagnetic activity:

* the margin h_sat - hmF2, as the fraction of samples below 0, 50 and 100 km.
  Below 0 the spacecraft is under the peak and there is no topside trace to the
  peak at all; within 100 km the trace is short and starts close to foF2;
* the hmF2 distribution's tails (0.1 and 99.9 percentiles) against the
  200-500 km envelope the threshold names;
* the same by Kp bin, from runs on the epoch's own indices with the storm
  models on.

The inputs are the per-sample ``.npz`` files ``tss study min_frequency``
writes.  hmF2 is a bottomside quantity and is identical across the three
topside models (checked by ``topside_independence``), so one topside per cell
is enough.

Nothing here runs IRI; ``scripts/sci02_threshold.py`` is the driver.
"""

from __future__ import annotations

import re
from collections.abc import Iterable
from pathlib import Path
from typing import Any

import numpy as np

MARGINS_KM = (0.0, 50.0, 100.0)
ENVELOPE_KM = (200.0, 500.0)
#: MIS-01 as Andrew restated it on 2026-10-20 (G3): the spacecraft above the
#: F2 peak in at least this fraction of the band-edge samples, in place of
#: Rev G's "at all times sampled".  99.9 per cent matches the 0.1/99.9
#: percentile reading of SCI-02's hmF2 envelope.
MIS01_ABOVE_PEAK_FRACTION = 0.999
LT_EDGES_H = tuple(range(0, 25, 3))
# Kp bins: quiet, unsettled/active, minor-moderate storm, strong-severe storm.
KP_EDGES = (0.0, 3.0, 5.0, 7.0, 9.01)
_NAME = re.compile(r"^(?P<cls>O\d)_(?P<alt>\d+)km_(?P<solar>.+)_(?P<topside>nequick|iri-2001|"
                   r"iri-2001-corrected)\.npz$")


def parse_sample_name(name: str) -> dict[str, Any] | None:
    """Orbit class, altitude, solar level and topside from a samples file name."""
    m = _NAME.match(Path(name).name)
    if not m:
        return None
    return dict(orbit_class=m["cls"], altitude_km=float(m["alt"]), solar=m["solar"],
                topside=m["topside"])


def sample_files(samples_dir: Path, topside: str = "nequick") -> list[Path]:
    """The samples files of one topside, sorted."""
    out = []
    for p in sorted(Path(samples_dir).glob("*.npz")):
        info = parse_sample_name(p.name)
        if info and info["topside"] == topside:
            out.append(p)
    return out


def margin_stats(h_sat_m: float, hmf2_m: np.ndarray,
                 margins_km: Iterable[float] = MARGINS_KM,
                 envelope_km: tuple[float, float] = ENVELOPE_KM) -> dict[str, float]:
    """Margin fractions and hmF2 tails for one set of samples.

    The fractions are of samples with a valid hmF2; NaN hmF2 (a failed IRI
    point) is excluded and counted.
    """
    hm = np.asarray(hmf2_m, dtype=float)
    ok = np.isfinite(hm)
    hm = hm[ok]
    out: dict[str, float] = {"n": float(hm.size), "n_invalid": float((~ok).sum())}
    if hm.size == 0:
        return out
    margin = (float(h_sat_m) - hm) / 1e3
    for m in margins_km:
        out[f"frac_margin_lt_{int(m)}km"] = float(np.mean(margin < m))
    out["min_margin_km"] = float(margin.min())
    hk = hm / 1e3
    for q, tag in ((0.1, "p001"), (1.0, "p01"), (50.0, "p50"), (99.0, "p99"), (99.9, "p999")):
        out[f"hmf2_{tag}_km"] = float(np.percentile(hk, q))
    out["hmf2_min_km"] = float(hk.min())
    out["hmf2_max_km"] = float(hk.max())
    lo, hi = envelope_km
    out["frac_hmf2_below_envelope"] = float(np.mean(hk < lo))
    out["frac_hmf2_above_envelope"] = float(np.mean(hk > hi))
    return out


def binned(values: np.ndarray, edges: Iterable[float], h_sat_m: float,
           hmf2_m: np.ndarray, label: str) -> list[dict[str, Any]]:
    """``margin_stats`` per bin of ``values`` (local time or Kp), half-open bins."""
    edges = list(edges)
    v = np.asarray(values, dtype=float)
    rows = []
    for lo, hi in zip(edges[:-1], edges[1:], strict=True):
        sel = (v >= lo) & (v < hi)
        st = margin_stats(h_sat_m, np.asarray(hmf2_m)[sel])
        rows.append({f"{label}_lo": lo, f"{label}_hi": hi, **st})
    return rows


def load_samples(path: Path) -> dict[str, Any]:
    """The arrays one samples file holds, with the extras as named columns."""
    d = np.load(path, allow_pickle=False)
    out = {k: d[k] for k in d.files}
    cols = [str(c) for c in out.get("extra_columns", [])]
    ex = out.get("extra")
    if ex is not None and ex.ndim == 2:
        for i, c in enumerate(cols):
            out[f"x_{c}"] = ex[:, i].astype(float)
    return out


def topside_independence(samples_dir: Path) -> float:
    """Largest |hmF2 difference| in metres between topsides of the same cell."""
    groups: dict[tuple, list[Path]] = {}
    for p in sorted(Path(samples_dir).glob("*.npz")):
        info = parse_sample_name(p.name)
        if info:
            key = (info["orbit_class"], info["altitude_km"], info["solar"])
            groups.setdefault(key, []).append(p)
    worst = 0.0
    for paths in groups.values():
        if len(paths) < 2:
            continue
        ref = np.load(paths[0])["hmf2_m"].astype(float)
        for p in paths[1:]:
            h = np.load(p)["hmf2_m"].astype(float)
            worst = max(worst, float(np.nanmax(np.abs(h - ref))))
    return worst


def cell_tables(samples_dir: Path, run: str, topside: str = "nequick"):
    """Per-cell, per-local-time and per-Kp tables for one run's samples."""
    cells, lt_rows, kp_rows = [], [], []
    for p in sample_files(samples_dir, topside):
        info = parse_sample_name(p.name)
        s = load_samples(p)
        h_sat = float(s["h_sat_m"])
        hm = s["hmf2_m"].astype(float)
        key = dict(run=run, **info)
        cells.append({**key, **margin_stats(h_sat, hm)})
        for row in binned(s["local_time_hours"], LT_EDGES_H, h_sat, hm, "lt_h"):
            lt_rows.append({**key, **row})
        if "x_kp" in s:
            for row in binned(s["x_kp"], KP_EDGES, h_sat, hm, "kp"):
                kp_rows.append({**key, **row})
    return cells, lt_rows, kp_rows


def verdict(cells: list[dict[str, Any]], altitudes_km=(500.0, 900.0),
            min_fraction_above: float = MIS01_ABOVE_PEAK_FRACTION) -> dict[str, Any]:
    """The SCI-02 threshold reading of a set of cell rows.

    * ``above_peak``: no sample in the h_sat range has hmF2 above the
      spacecraft (Rev G's MIS-01, "above the F2 peak at all times sampled").
    * ``above_peak_fraction_met``: in every cell at least
      ``min_fraction_above`` of the samples have the spacecraft above the peak
      (MIS-01 as restated on 2026-10-20).
    * ``in_envelope``: the 0.1 and 99.9 percentiles of hmF2 lie in 200-500 km.
    The worst cells are named so the reading can be checked.
    """
    lo, hi = altitudes_km
    rows = [c for c in cells if lo <= c["altitude_km"] <= hi and c.get("n", 0) > 0]
    if not rows:
        return {"n_cells": 0}
    worst_below = max(rows, key=lambda c: c["frac_margin_lt_0km"])
    worst_100 = max(rows, key=lambda c: c["frac_margin_lt_100km"])
    lowest = min(rows, key=lambda c: c["hmf2_p001_km"])
    highest = max(rows, key=lambda c: c["hmf2_p999_km"])

    def name(c):
        return f"{c['run']}:{c['orbit_class']}/{c['altitude_km']:.0f} km/{c['solar']}"

    return {
        "n_cells": len(rows),
        "above_peak": worst_below["frac_margin_lt_0km"] == 0.0,
        "above_peak_fraction_required": min_fraction_above,
        "min_fraction_above_peak": 1.0 - worst_below["frac_margin_lt_0km"],
        "above_peak_fraction_met": (1.0 - worst_below["frac_margin_lt_0km"]) >= min_fraction_above,
        "worst_below_peak": {"cell": name(worst_below),
                             "fraction": worst_below["frac_margin_lt_0km"],
                             "min_margin_km": worst_below["min_margin_km"]},
        "worst_within_100km": {"cell": name(worst_100),
                               "fraction": worst_100["frac_margin_lt_100km"]},
        "in_envelope": (lowest["hmf2_p001_km"] >= ENVELOPE_KM[0]
                        and highest["hmf2_p999_km"] <= ENVELOPE_KM[1]),
        "lowest_hmf2_p001": {"cell": name(lowest), "km": lowest["hmf2_p001_km"],
                             "min_km": lowest["hmf2_min_km"]},
        "highest_hmf2_p999": {"cell": name(highest), "km": highest["hmf2_p999_km"],
                              "max_km": highest["hmf2_max_km"]},
    }
