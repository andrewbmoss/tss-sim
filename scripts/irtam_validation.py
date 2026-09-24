#!/usr/bin/env python3
"""Compare IRTAM maps with GIRO Digisonde measurements, and IRI-CCIR with both.

Reproduces ``docs/irtam-validation.md``.  Reads cached files only; it never
fetches, so the numbers in the note come from exactly the bytes on disk.

    python scripts/irtam_validation.py DATA_DIR [--f107 175.3] [--figure OUT.png]

``DATA_DIR`` holds ``irtam/IRTAM_{foF2,hmF2}_COEFFS_<YYYYMMDD>_<HH>0000.ASC``
(hourly archive files for one day) and ``giro/<URSI>.txt`` (DIDBase responses
for the same day with foF2 and hmF2).  Each hourly coefficient set is evaluated
at its time of validity at every station and compared with the station's
record nearest that hour, within 7.5 minutes, keeping only records with an
autoscaling confidence score of at least 70 or manual scaling.

The CCIR baseline needs PyIRI (``pip install PyIRI``) and the 12-month
smoothed F10.7 for the month.  Without PyIRI the script reports IRTAM only.
DIDBase's own ``foF2p`` column is not used as a baseline: on 2024-03-01 it gave
about 4.8 MHz for the Millstone Hill afternoon near solar maximum, where
IRI-CCIR with the observed smoothed F10.7 gives about 11 MHz.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))

from tss.data.giro import parse_scaled  # noqa: E402
from tss.data.irtam import parse_coefficients  # noqa: E402
from tss.l0_ionosphere.geomag import modip_deg  # noqa: E402
from tss.l0_ionosphere.irtam_maps import evaluate  # noqa: E402

MATCH_S = 450.0
MIN_CS = 70


def load_irtam(data: Path) -> tuple[dict, dict]:
    sets: dict[str, dict[int, object]] = {"foF2": {}, "hmF2": {}}
    for ch in sets:
        for f in sorted((data / "irtam").glob(f"IRTAM_{ch}_COEFFS_*.ASC")):
            c = parse_coefficients(f.read_bytes())
            sets[ch][c.valid_time.hour] = c
    if not sets["foF2"]:
        raise SystemExit(f"no IRTAM files under {data / 'irtam'}")
    assimilated = set().union(*(set(c.stations) for c in sets["foF2"].values()))
    return sets, {"assimilated": assimilated}


def ccir_baseline(stations, day, f107):
    try:
        import PyIRI
        import PyIRI.edp_update as eu
    except ImportError:
        return None
    lat = np.array([s.lat_deg for s in stations])
    lon = np.array([s.lon_deg for s in stations])
    lon = np.where(lon > 180, lon - 360, lon)
    f2, *_ = eu.IRI_density_1day(day.year, day.month, day.day, np.arange(24.0), lon, lat,
                                 np.array([300.0]), f107, PyIRI.coeff_dir, 0)
    return np.sqrt(f2["Nm"] / 1.2404e10), np.asarray(f2["hm"])


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("data", type=Path)
    ap.add_argument("--f107", type=float, default=None, help="12-month smoothed F10.7 for the CCIR baseline")
    ap.add_argument("--figure", type=Path, default=None)
    ap.add_argument("--json", type=Path, default=None)
    ap.add_argument("--note", action="append", default=[], metavar="URSI=TEXT",
                    help="replace a station's subtitle; '|' starts a new line")
    args = ap.parse_args(argv)

    sets, meta = load_irtam(args.data)
    any_set = next(iter(sets["foF2"].values()))
    day = any_set.valid_time.date()
    stations = []
    for f in sorted((args.data / "giro").glob("*.txt")):
        r = parse_scaled(f.read_bytes())
        if len(r) and "foF2" in r.values and "hmF2" in r.values:
            stations.append(r)
    if not stations:
        raise SystemExit("no GIRO station files with data")

    base = ccir_baseline(stations, day, args.f107) if args.f107 else None
    mid = f"{day.year}-{day.month:02d}-15"
    lat = np.array([s.lat_deg for s in stations])
    lon = np.array([s.lon_deg for s in stations])
    mod, mod_src = modip_deg(lat, lon, epoch=mid)

    report = {"day": str(day), "f107_smoothed": args.f107, "modip": mod_src,
              "n_assimilated": len(meta["assimilated"]), "stations": {}}
    series = {}
    for j, s in enumerate(stations):
        rows = []
        for h in sorted(sets["foF2"]):
            t = np.datetime64(f"{day}T{h:02d}:00:00", "ms")
            dts = np.abs((s.times - t).astype("timedelta64[ms]").astype(float)) / 1e3
            ok = s.good(MIN_CS) & (dts <= MATCH_S)
            if not ok.any():
                continue
            i = np.where(ok)[0][np.argmin(dts[ok])]
            fo = evaluate(sets["foF2"][h], lat[j:j + 1], lon[j:j + 1], modip_deg=mod[j:j + 1])[0]
            hm = evaluate(sets["hmF2"][h], lat[j:j + 1], lon[j:j + 1], modip_deg=mod[j:j + 1])[0]
            rows.append({
                "hour": h,
                "foF2_meas": float(s.values["foF2"][i]), "hmF2_meas": float(s.values["hmF2"][i]),
                "foF2_irtam": float(fo), "hmF2_irtam": float(hm),
                "foF2_ccir": float(base[0][h, j]) if base else None,
                "hmF2_ccir": float(base[1][h, j]) if base else None,
            })

        def stat(model, q, rows=rows):
            d = np.array([r[f"{q}_{model}"] - r[f"{q}_meas"] for r in rows
                          if r[f"{q}_{model}"] is not None and np.isfinite(r[f"{q}_meas"])])
            if not d.size:
                return None
            return {"n": int(d.size), "bias": float(d.mean()), "rms": float(np.sqrt((d ** 2).mean()))}

        rec = {
            "name": s.station, "lat_deg": s.lat_deg, "lon_deg": s.lon_deg,
            "assimilated": s.ursi_code in meta["assimilated"],
            "foF2_irtam": stat("irtam", "foF2"), "hmF2_irtam": stat("irtam", "hmF2"),
            "foF2_ccir": stat("ccir", "foF2") if base else None,
            "hmF2_ccir": stat("ccir", "hmF2") if base else None,
        }
        report["stations"][s.ursi_code] = rec
        # Model curves at every hour, so a plot never draws a model line straight
        # across hours with no measurement.
        full = []
        for h in sorted(sets["foF2"]):
            full.append({
                "hour": h,
                "foF2_irtam": float(evaluate(sets["foF2"][h], lat[j:j + 1], lon[j:j + 1], modip_deg=mod[j:j + 1])[0]),
                "hmF2_irtam": float(evaluate(sets["hmF2"][h], lat[j:j + 1], lon[j:j + 1], modip_deg=mod[j:j + 1])[0]),
                "foF2_ccir": float(base[0][h, j]) if base else None,
                "hmF2_ccir": float(base[1][h, j]) if base else None,
            })
        series[s.ursi_code] = {"matched": rows, "model": full}
        line = f"{s.ursi_code} {s.station:16s} {'assimilated' if rec['assimilated'] else 'not assimilated':16s}"
        for q, unit in (("foF2", "MHz"), ("hmF2", "km")):
            a, c = rec[f"{q}_irtam"], rec[f"{q}_ccir"]
            line += f" | {q} rms IRTAM {a['rms']:.2f}" if a else f" | {q} -"
            line += f" CCIR {c['rms']:.2f} {unit}" if c else f" {unit}"
        print(line + f" (n={rec['foF2_irtam']['n'] if rec['foF2_irtam'] else 0})")

    if args.json:
        args.json.write_text(json.dumps({"report": report, "series": series}, indent=1))
    if args.figure:
        notes = {k: v.replace("|", "\n") for k, v in (n.split("=", 1) for n in args.note)}
        _figure(series, report, args.figure, notes)
    return 0


def _figure(series, report, out: Path, notes: dict | None = None) -> None:
    from tss.figures import style

    order = sorted(series, key=lambda k: (not report["stations"][k]["assimilated"], k))
    notes = notes or {}
    fig, axes = style.figure(2, len(order), figsize=(3.2 * len(order), 5.4), sharex=True)
    axes = np.atleast_2d(axes).reshape(2, len(order))
    for col, code in enumerate(order):
        rows = series[code]["matched"]
        model = series[code]["model"]
        rec = report["stations"][code]
        hrs = np.array([r["hour"] for r in rows])
        mh = np.array([r["hour"] for r in model])
        for row, (q, unit) in enumerate((("foF2", "MHz"), ("hmF2", "km"))):
            ax = axes[row, col]
            style.grid(ax, "y")
            ax.plot(hrs, [r[f"{q}_meas"] for r in rows], "o", ms=4.5, color=style.INK_2,
                    label="GIRO measured", zorder=3)
            ax.plot(mh, [r[f"{q}_irtam"] for r in model], "-", lw=2, color=style.SERIES[0],
                    label="IRTAM", zorder=2)
            if model and model[0][f"{q}_ccir"] is not None:
                ax.plot(mh, [r[f"{q}_ccir"] for r in model], "-", lw=2, color=style.SERIES[1],
                        label="IRI-CCIR", zorder=1)
            if row == 0:
                tag = notes.get(code) or ("assimilated" if rec["assimilated"] else "not assimilated")
                ax.set_title(f"{rec['name'].title()} ({code})\n{tag}", fontsize=9)
            if col == 0:
                ax.set_ylabel(f"{q} [{unit}]")
            if row == 1:
                ax.set_xlabel("UT [h]")
                ax.set_xticks(range(0, 25, 6))
    axes[0, 0].legend(loc="upper left", fontsize=7, frameon=False)
    style.caption(fig, f"{report['day']}: hourly IRTAM coefficient sets evaluated at each "
                  "station vs the nearest DIDBase record (CS >= 70) within 7.5 min. "
                  f"IRI-CCIR with 12-month smoothed F10.7 = {report['f107_smoothed']}.")
    fig.tight_layout()
    style.savefig(fig, out.with_suffix(""), formats=("png",))


if __name__ == "__main__":
    raise SystemExit(main())
