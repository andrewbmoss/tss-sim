"""WP0 Phase 3: the SCI-02 threshold from band-edge samples.

    python scripts/sci02_threshold.py OUT NAME=RUN_DIR [NAME=RUN_DIR ...]

Each RUN_DIR is the output of ``tss study min_frequency`` with
``keep_samples: true``.  Writes to OUT:

* ``sci02_cells.csv``     one row per run, orbit class, altitude and solar level
* ``sci02_local_time.csv``  the same in 3-hour local-time bins
* ``sci02_kp.csv``        the same in Kp bins (0-3, 3-5, 5-7, 7-9)
* ``sci02_threshold.json``  the verdict per run and over all runs
* ``sci02_threshold.md``  the tables a reader needs

Only the NeQuick files are read: hmF2 does not depend on the topside model,
and the script checks that on every run it reads.
"""

from __future__ import annotations

import csv
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from tss.studies.sci02_margin import cell_tables, topside_independence, verdict  # noqa: E402


def _write(path: Path, rows: list[dict]) -> None:
    if not rows:
        return
    keys: list[str] = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)


def _pct(x: float) -> str:
    return f"{100 * x:.2f}"


def main(argv: list[str]) -> int:
    if len(argv) < 2:
        print(__doc__)
        return 2
    out = Path(argv[0])
    out.mkdir(parents=True, exist_ok=True)
    cells, lts, kps, verdicts, indep = [], [], [], {}, {}
    for spec in argv[1:]:
        name, _, path = spec.partition("=")
        samples = Path(path) / "samples"
        indep[name] = topside_independence(samples)
        c, lt, kp = cell_tables(samples, name)
        cells += c
        lts += lt
        kps += kp
        verdicts[name] = verdict(c)
    verdicts["all"] = verdict(cells)
    _write(out / "sci02_cells.csv", cells)
    _write(out / "sci02_local_time.csv", lts)
    _write(out / "sci02_kp.csv", kps)
    (out / "sci02_threshold.json").write_text(
        json.dumps({"verdicts": verdicts, "topside_hmf2_max_abs_diff_m": indep}, indent=1),
        encoding="utf-8")

    lines = ["# SCI-02 threshold: spacecraft height above the F2 peak", "",
             "Fractions in per cent of samples.  margin = h_sat - hmF2.", "",
             "## Per cell", "",
             "| run | class | h_sat km | solar | n | margin<0 | <50 km | <100 km | min margin km "
             "| hmF2 p0.1 | p99.9 | min | max |",
             "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for c in cells:
        lines.append(
            f"| {c['run']} | {c['orbit_class']} | {c['altitude_km']:.0f} | {c['solar']} | "
            f"{int(c['n'])} | {_pct(c['frac_margin_lt_0km'])} | {_pct(c['frac_margin_lt_50km'])} | "
            f"{_pct(c['frac_margin_lt_100km'])} | {c['min_margin_km']:.0f} | "
            f"{c['hmf2_p001_km']:.0f} | {c['hmf2_p999_km']:.0f} | {c['hmf2_min_km']:.0f} | "
            f"{c['hmf2_max_km']:.0f} |")
    lines += ["", "## Local time, cells with any sample within 100 km of the peak", "",
              "| run | class | h_sat km | solar | LT h | n | margin<0 | <100 km | hmF2 p99.9 | max |",
              "|---|---|---|---|---|---|---|---|---|---|"]
    for r in lts:
        if r.get("n", 0) and r["frac_margin_lt_100km"] > 0:
            lines.append(
                f"| {r['run']} | {r['orbit_class']} | {r['altitude_km']:.0f} | {r['solar']} | "
                f"{r['lt_h_lo']:.0f}-{r['lt_h_hi']:.0f} | {int(r['n'])} | "
                f"{_pct(r['frac_margin_lt_0km'])} | {_pct(r['frac_margin_lt_100km'])} | "
                f"{r['hmf2_p999_km']:.0f} | {r['hmf2_max_km']:.0f} |")
    lines += ["", "## Kp bins, 500 km", "",
              "| run | class | solar | Kp | n | margin<0 | <100 km | hmF2 p50 | p99.9 | max |",
              "|---|---|---|---|---|---|---|---|---|---|"]
    for r in kps:
        if r["altitude_km"] == 500 and r.get("n", 0):
            lines.append(
                f"| {r['run']} | {r['orbit_class']} | {r['solar']} | {r['kp_lo']:.0f}-"
                f"{min(r['kp_hi'], 9):.0f} | {int(r['n'])} | {_pct(r['frac_margin_lt_0km'])} | "
                f"{_pct(r['frac_margin_lt_100km'])} | {r['hmf2_p50_km']:.0f} | "
                f"{r['hmf2_p999_km']:.0f} | {r['hmf2_max_km']:.0f} |")
    lines += ["", "## Verdicts", "", "```", json.dumps(verdicts, indent=1), "```", "",
              f"Largest hmF2 difference between topsides, m: {json.dumps(indep)}", ""]
    (out / "sci02_threshold.md").write_text("\n".join(lines), encoding="utf-8")
    print(json.dumps(verdicts["all"], indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
