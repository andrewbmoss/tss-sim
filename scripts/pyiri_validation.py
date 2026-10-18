#!/usr/bin/env python3
"""Compare PyIRI with Fortran IRI-2020 (PHaRLAP) at a matched solar level.

Reproduces ``docs/pyiri-validation.md``.

    python scripts/pyiri_validation.py OUT_DIR [--f107 120]

The ``iri2020`` tier is PyIRI, an independent Python implementation of IRI.
This script asks how far it sits from the Fortran IRI-2020 that PHaRLAP's mex
runs, at the same points, the same times and the same solar level, on the
quantities this work rests on: NmF2, hmF2, and the topside density at the
orbit altitudes, where the plasma frequency at the spacecraft sets the lower
band edge (INS-01a).

Both models are run with IRI-2020's default options and the NeQuick topside.
PHaRLAP takes a user R12, which is set from the F10.7 by inverting the Davies
relation :func:`tss.l0_ionosphere.base.r12_to_f107` already uses; PyIRI takes
the F10.7.  That pairing is not exact, because the two codes interpolate their
coefficient maps in different solar indices, and the comparison carries that
difference with it.

Needs PyIRI and MATLAB with PHaRLAP.  Writes ``pyiri_vs_iri2020.csv`` (one row
per point) and ``pyiri_vs_iri2020.json`` (the summary) into ``OUT_DIR``.
"""

from __future__ import annotations

import argparse
import datetime as dt
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from tss.constants import plasma_frequency  # noqa: E402
from tss.l0_ionosphere.backends import MatlabIRIBackend, PyIRIBackend  # noqa: E402
from tss.l0_ionosphere.base import r12_to_f107  # noqa: E402

ALT_M = np.arange(1.0e5, 2.0e6 + 1.0, 1.0e4)
ORBIT_ALTITUDES_M = (5.0e5, 7.0e5, 9.0e5)
LATS = np.arange(-80.0, 81.0, 20.0)
LONS = np.arange(-180.0, 180.0, 45.0)
UTS = (0, 6, 12, 18)
DAYS = (dt.date(2024, 3, 15), dt.date(2024, 6, 15))


def f107_to_r12(f107: float) -> float:
    """Invert 63.7 + 0.728 R + 0.00089 R^2 = F10.7 for R."""
    a, b, c = 0.00089, 0.728, 63.7 - float(f107)
    return float((-b + np.sqrt(b * b - 4.0 * a * c)) / (2.0 * a))


def _points():
    lat, lon, when = [], [], []
    for day in DAYS:
        for ut in UTS:
            for la in LATS:
                for lo in LONS:
                    lat.append(float(la))
                    lon.append(float(lo))
                    when.append(dt.datetime(day.year, day.month, day.day, ut))
    return lat, lon, when


def _stats(x: np.ndarray) -> dict[str, float]:
    x = np.asarray(x, dtype=float)
    x = x[np.isfinite(x)]
    return {
        "median": float(np.median(x)),
        "mean": float(np.mean(x)),
        "p10": float(np.percentile(x, 10)),
        "p90": float(np.percentile(x, 90)),
        "abs_p90": float(np.percentile(np.abs(x), 90)),
        "n": int(x.size),
    }


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("out_dir", type=Path)
    parser.add_argument("--f107", type=float, default=120.0)
    args = parser.parse_args(argv)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    r12 = f107_to_r12(args.f107)
    lat, lon, when = _points()

    pyiri = PyIRIBackend(f107=args.f107)
    ne_py = pyiri.profiles(lat, lon, when, ALT_M)
    extra_py = list(pyiri._last_extra)

    fortran = MatlabIRIBackend(topside="NeQuick", r12=r12, study="pyiri_validation")
    ne_ft = fortran.profiles(lat, lon, when, ALT_M)

    def peak(ne: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
        i = np.argmax(ne, axis=1)
        return ne[np.arange(ne.shape[0]), i], ALT_M[i]

    nm_ft, hm_ft = peak(ne_ft)
    nm_py = np.array([e["nmf2"] for e in extra_py])
    hm_py = np.array([e["hmf2_km"] * 1e3 for e in extra_py])

    rows = []
    header = ["lat", "lon", "month", "ut", "nmf2_pyiri", "nmf2_iri2020", "hmf2_pyiri_km",
              "hmf2_iri2020_km"]
    for h in ORBIT_ALTITUDES_M:
        header += [f"fp_pyiri_{h/1e3:.0f}km_mhz", f"fp_iri2020_{h/1e3:.0f}km_mhz"]
    ix_orbit = [int(np.argmin(np.abs(ALT_M - h))) for h in ORBIT_ALTITUDES_M]
    for k in range(len(lat)):
        row = [lat[k], lon[k], when[k].month, when[k].hour, nm_py[k], nm_ft[k],
               hm_py[k] / 1e3, hm_ft[k] / 1e3]
        for j in ix_orbit:
            row += [float(plasma_frequency(ne_py[k, j])) / 1e6,
                    float(plasma_frequency(ne_ft[k, j])) / 1e6]
        rows.append(row)
    np.savetxt(args.out_dir / "pyiri_vs_iri2020.csv", np.array(rows), delimiter=",",
               header=",".join(header), comments="")

    summary: dict[str, object] = {
        "f107": args.f107,
        "r12_for_pharlap": r12,
        "r12_back_to_f107": float(r12_to_f107(r12)),
        "n_points": len(lat),
        "days": [d.isoformat() for d in DAYS],
        "uts_hours": list(UTS),
        "pyiri": pyiri.describe(),
        "iri2020": {"label": fortran.label, "provenance": fortran.provenance},
        "nmf2_ratio_pyiri_over_iri2020": _stats(nm_py / nm_ft),
        "hmf2_difference_km": _stats((hm_py - hm_ft) / 1e3),
        "hmf2_grid_step_km": float((ALT_M[1] - ALT_M[0]) / 1e3),
    }
    for h, j in zip(ORBIT_ALTITUDES_M, ix_orbit, strict=True):
        ratio = ne_py[:, j] / ne_ft[:, j]
        fp_py = plasma_frequency(ne_py[:, j]) / 1e6
        fp_ft = plasma_frequency(ne_ft[:, j]) / 1e6
        summary[f"ne_ratio_at_{h/1e3:.0f}km"] = _stats(ratio)
        summary[f"fp_difference_mhz_at_{h/1e3:.0f}km"] = _stats(fp_py - fp_ft)
        summary[f"fp_5th_percentile_mhz_at_{h/1e3:.0f}km"] = {
            "pyiri": float(np.percentile(fp_py, 5)),
            "iri2020": float(np.percentile(fp_ft, 5)),
        }
    (args.out_dir / "pyiri_vs_iri2020.json").write_text(json.dumps(summary, indent=2))
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
