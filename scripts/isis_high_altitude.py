#!/usr/bin/env python3
"""WP0 Phase 5: ISIS-1 against epoch-driven IRI-2020 at 900-1400 km.

    python scripts/isis_high_altitude.py run     OUT_DIR --crc CRC.asc [--matlab matlab]
    python scripts/isis_high_altitude.py analyse OUT_DIR [--below PREVIOUS_DIR]

The Rev G check of the band-edge method against ISIS-1 (T1, and the epoch
rerun of 2026-10-18 in ``docs/topside-model-campaign.md``) covered spacecraft at
500-900 km.  SCI-02's goal runs to 1400 km, where the Phase 4 band edges rest
on the topside models' low tail, so this repeats the comparison for the CRC
hand-scaled profiles taken from 900 to 1400 km: 5,233 profiles with quality
index 0-7 and a recorded Rz12, 65 of them at solar minimum (Rz12 <= 30, all
1975-77).

``run`` evaluates PHaRLAP's iri2020 at each profile's own time, place and
spacecraft height through ``scripts/topside_isis_epoch_iri.m``, in the six
configurations that script defines: user R12 at the record's 5-unit Rz12 bucket
with NeQuick (the Rev G route), and on the epoch's own indices (R12 = -1) with
NeQuick storm-off, NeQuick with the foF2 storm model (IRI's epoch default),
NeQuick with every storm model, IRI-2001 and IRI-2001 corrected.  It also fits
the NeQuick topside H0 to each measured profile, through the profile's deepest
point, as ``topside_campaign.py isis-epoch`` does.

``analyse`` writes ``summary.json`` and ``summary.md``.  Per 100 km band and
solar level:

* the median model/measured ratio of the plasma frequency at the spacecraft,
  in dB (20 log10);
* the 5th percentile of the measured fN_sat against each configuration's, and
  against the ensemble edge INS-01a's construction uses (the minimum of the
  three epoch topsides' 5th percentiles);
* the realised coverage of that ensemble edge: the fraction of measured fN_sat
  at or above it, 95 per cent if the models' low tail is right;
* the fitted H0 against IRI-2020's H0 law on the measured peak.

``analyse`` also splits the solar-maximum tail by latitude, local time and CRC
quality index, and with ``--goal-run`` bounds Phase 4's coverage at 0.33 MHz by
dividing the model's fN_sat by the measured tail ratio of the matching band.

``--below`` adds the 500-900 km rows of an earlier ``isis-epoch`` run
(``summary_epoch.json``) to the markdown for comparison.
"""

from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from topside_campaign import (  # noqa: E402
    EPOCH_CONFS,
    EXI,
    PHARLAP,
    _jsonable,
    fit_h0,
    fp_mhz,
    h0_iri,
    isis_points,
)

H_LO, H_HI = 900.0, 1400.0
BANDS = tuple((float(b), float(b + 100)) for b in range(900, 1400, 100))
EPOCH_TOPSIDES = ("epoch_NeQuick", "epoch_IRI2001", "epoch_IRI2001corr")
REGIMES = {"solar_max": lambda rz: rz >= 100, "intermediate": lambda rz: (rz > 30) & (rz < 100),
           "solar_min": lambda rz: rz <= 30, "all": lambda rz: np.ones_like(rz, bool)}
MIN_N = 5


def select(crc: Path, lo: float = H_LO, hi: float = H_HI):
    """CRC records with lo < h_sat <= hi, quality 0-7, a finite fN_sat and an Rz12."""
    sys.path.insert(0, str(crc.parent))
    import crc_read  # Dissertation/ISISValidation/crc_read.py

    _, recs = crc_read.read_crc(str(crc), want_profiles=True)
    sel = [x for x in recs if lo < x["h_sat_km"] <= hi and 0 <= x["iquality"] <= 7
           and np.isfinite(x["fn_sat_mhz"]) and x["rz12"] is not None]
    for x in sel:
        x["rz_bucket"] = int(round(x["rz12"] / 5.0) * 5)
    return sel


def cmd_run(out: Path, crc: Path, matlab: str = "matlab", reuse_iri: bool = False) -> None:
    import scipy.io as sio

    out.mkdir(parents=True, exist_ok=True)
    order = select(crc)
    pts = isis_points(order)
    n = len(order)
    print(f"{n} CRC profiles at {H_LO:.0f}-{H_HI:.0f} km")
    mat = out / "iri_epoch.mat"
    if not (reuse_iri and mat.is_file()):
        np.savetxt(out / "isis_points.csv", pts, delimiter=",", fmt="%.6g")
        sd = str(ROOT / "scripts").replace("\\", "/")
        pc, oc = (str(x).replace("\\", "/") for x in (out / "isis_points.csv", mat))
        call = (f"addpath('{sd}'); topside_isis_epoch_iri('{pc}','{oc}',"
                f"'{PHARLAP.replace(chr(92), '/')}')")
        t0 = time.time()
        subprocess.run([matlab, "-batch", call], check=True)
        print(f"IRI-2020 at the ISIS-1 points, {len(EPOCH_CONFS)} configurations, "
              f"{time.time() - t0:.0f} s")
    M = sio.loadmat(mat)
    names = [str(c[0]) for c in M["conf_names"].ravel()]
    if tuple(names) != EPOCH_CONFS:
        raise RuntimeError(f"unexpected configurations in {mat}: {names}")
    ne_iri = np.asarray(M["ne_sat"], dtype=float)
    ex = np.asarray(M["extra"], dtype=float)[:, :, EPOCH_CONFS.index("epoch_NeQuick")]
    idx = {k: ex[:, c] for k, c in EXI.items()}
    n_failed = np.asarray(M["n_failed"]).ravel().astype(int)

    obs = {k: np.full(n, np.nan) for k in ("nmf2_obs", "hmf2_obs", "h0_fit", "fit_rms", "depth_km")}
    for i, x in enumerate(order):
        h, ne = np.asarray(x["h_km"], float), np.asarray(x["ne_cm3"], float) * 1.0e6
        j = int(np.argmax(ne))
        obs["nmf2_obs"][i], obs["hmf2_obs"][i] = ne[j], h[j]
        obs["depth_km"][i] = h[0] - h[j]
        use = (h >= h[j] + 30.0) & np.isfinite(ne) & (ne > 0)
        if use.sum() >= 4 and obs["depth_km"][i] >= 150.0:
            obs["h0_fit"][i], obs["fit_rms"][i] = fit_h0(h[use], ne[use], h[j], ne[j])

    np.savez_compressed(
        out / "isis_high.npz", conf_names=np.array(EPOCH_CONFS), n_failed=n_failed,
        year=pts[:, 0], month=pts[:, 1], day=pts[:, 2], ut=pts[:, 3] + pts[:, 4] / 60.0,
        lat=pts[:, 5], lon=pts[:, 6], h_sat=pts[:, 7], rz_bucket=pts[:, 8],
        lt=np.array([x["lt_hour"] for x in order], float),
        rz12_crc=np.array([x["rz12"] for x in order], float),
        quality=np.array([x["iquality"] for x in order], float),
        fn_obs=np.array([x["fn_sat_mhz"] for x in order], float),
        ne_iri=ne_iri, **{f"idx_{k}": v for k, v in idx.items()}, **obs)
    print(f"wrote {out / 'isis_high.npz'}; IRI failures per configuration {n_failed.tolist()}")


def band_rows(fn_obs, models, hs, rz, h0_ratio, bands=BANDS, regimes=REGIMES, min_n=MIN_N):
    """The per-band, per-regime table.  ``models`` maps a name to fN_sat in MHz."""
    out: dict = {}
    db = {k: 20.0 * np.log10(v / fn_obs) for k, v in models.items()}
    for rname, rfun in regimes.items():
        reg = rfun(rz)
        out[rname] = {}
        for lo, hi in list(bands) + [(bands[0][0], bands[-1][1])]:
            m = reg & (hs > lo) & (hs <= hi)
            if m.sum() < min_n:
                continue
            row = {"n": int(m.sum()), "median_db": {}, "p05_mhz": {
                "observed": float(np.percentile(fn_obs[m], 5))}}
            for k in models:
                row["median_db"][k] = float(np.nanmedian(db[k][m]))
                row["p05_mhz"][k] = float(np.nanpercentile(models[k][m], 5))
            ens = min(row["p05_mhz"][k] for k in EPOCH_TOPSIDES if k in row["p05_mhz"])
            row["ensemble_p05_mhz"] = ens
            row["ensemble_over_observed_p05"] = ens / row["p05_mhz"]["observed"]
            row["realised_coverage_of_ensemble_edge"] = float(np.mean(fn_obs[m] >= ens))
            hr = h0_ratio[m]
            hr = hr[np.isfinite(hr)]
            if hr.size >= min_n:
                row["h0_fit_over_iri_law"] = [float(np.percentile(hr, q)) for q in (10, 50, 90)]
                row["n_h0_fit"] = int(hr.size)
            out[rname][f"{int(lo)}-{int(hi)}"] = row
    return out


def _tail(fn_obs, models, m):
    ens = min(float(np.nanpercentile(models[k][m], 5)) for k in EPOCH_TOPSIDES)
    po = float(np.percentile(fn_obs[m], 5))
    return {"n": int(m.sum()), "p05_observed_mhz": po, "ensemble_p05_mhz": ens,
            "ensemble_over_observed_p05": ens / po,
            "realised_coverage_of_ensemble_edge": float(np.mean(fn_obs[m] >= ens)),
            "median_db_nequick": float(np.nanmedian(
                20 * np.log10(models["epoch_NeQuick"][m] / fn_obs[m])))}


def splits(fn_obs, models, hs, rz, lat, lt, quality, min_n=20):
    """Where the low tail is overstated: latitude, local time, CRC quality index.

    Polar here is |lat| >= 60 deg.  The CRC quality index of the polar records
    is mostly 6-7, so latitude and quality are confounded in this archive.
    """
    sm, smin = rz >= 100, rz <= 30
    pol = np.abs(lat) >= 60.0
    night = (lt >= 18.0) | (lt < 6.0)
    groups = {
        "solar_max_abs_lat_0_30": sm & (np.abs(lat) < 30),
        "solar_max_abs_lat_30_60": sm & (np.abs(lat) >= 30) & ~pol,
        "solar_max_abs_lat_60_70": sm & (np.abs(lat) >= 60) & (np.abs(lat) < 70),
        "solar_max_abs_lat_70_80": sm & (np.abs(lat) >= 70) & (np.abs(lat) < 80),
        "solar_max_abs_lat_80_90": sm & (np.abs(lat) >= 80),
        "solar_max_day_06_18": sm & ~night, "solar_max_night_18_06": sm & night,
        "solar_max_quality_le5": sm & (quality <= 5), "solar_max_quality_ge6": sm & (quality >= 6),
        "solar_max_nonpolar_quality_le5": sm & ~pol & (quality <= 5),
        "solar_max_nonpolar_quality_ge6": sm & ~pol & (quality >= 6),
        "solar_max_polar_quality_ge6": sm & pol & (quality >= 6),
        "solar_min_polar": smin & pol,
    }
    for lo, hi in BANDS:
        b = (hs > lo) & (hs <= hi)
        groups[f"solar_max_nonpolar_{int(lo)}_{int(hi)}"] = sm & ~pol & b
        groups[f"solar_max_polar_{int(lo)}_{int(hi)}"] = sm & pol & b
    out = {k: _tail(fn_obs, models, m) for k, m in groups.items() if m.sum() >= min_n}
    out["n_solar_min_polar"] = int((smin & pol).sum())
    out["n_solar_min_nonpolar"] = int((smin & ~pol).sum())
    return out


def goal_sensitivity(goal_dir: Path, solar_max_bands: dict) -> dict:
    """Phase 4's coverage at 0.33 MHz if the model's fN_sat is divided by the tail ratio.

    A uniform division by the solar-maximum ensemble/observed 5th-percentile
    ratio of the matching ISIS-1 band is a pessimistic bound: the medians agree
    to about 0.5 dB, so the overstatement is in the tail, and 0.33 MHz sits
    near the 35th percentile at 1400 km.  NeQuick, the governing topside.
    """
    res = {}
    for alt in (1000, 1100, 1200, 1300, 1400):
        band = f"{alt - 100}-{alt}"
        if band not in solar_max_bands:
            continue
        r = solar_max_bands[band]["ensemble_over_observed_p05"]
        for cls in ("O1", "O2"):
            for sol in ("min", "min2008", "min2019", "max"):
                f = Path(goal_dir) / "samples" / f"{cls}_{alt}km_{sol}_nequick.npz"
                if not f.is_file():
                    continue
                fn = np.load(f)["fn_sat_hz"].astype(float) / 1e6
                fn = fn[np.isfinite(fn)]
                res[f"{cls}_{alt}_{sol}"] = {
                    "tail_ratio": r, "coverage_at_033": float(np.mean(fn >= 0.33)),
                    "coverage_at_033_bound": float(np.mean(fn / r >= 0.33)),
                    "p05_mhz": float(np.percentile(fn, 5)),
                    "p05_mhz_bound": float(np.percentile(fn, 5)) / r}
    return res


def cmd_analyse(out: Path, below: Path | None = None, goal_samples: Path | None = None) -> dict:
    J = dict(np.load(out / "isis_high.npz", allow_pickle=False))
    hs, fn_obs, rz = J["h_sat"], J["fn_obs"], J["rz12_crc"]
    names = [str(x) for x in J["conf_names"]]
    models = {k: fp_mhz(J["ne_iri"][:, c]) for c, k in enumerate(names)}
    r12, m3 = J["idx_r12"], J["idx_m3000"]
    fo = np.sqrt(J["nmf2_obs"] / 1.24e10)
    h0_law = h0_iri(fo, m3, J["hmf2_obs"], r12)
    ratio = J["h0_fit"] / h0_law
    s = {"n": int(hs.size), "n_failed": J["n_failed"].tolist(),
         "n_by_regime": {k: int(f(rz).sum()) for k, f in REGIMES.items()},
         "years_solar_min": sorted({int(y) for y in J["year"][rz <= 30]}),
         "bands": band_rows(fn_obs, models, hs, rz, ratio)}
    kp = J["idx_kp"]
    sm = rz >= 100
    s["storm_solar_max"] = {}
    for cname, cm in {"kp_lt3": kp < 3, "kp_3_4": (kp >= 3) & (kp < 4), "kp_ge4": kp >= 4}.items():
        m = sm & cm
        if m.sum() >= 20:
            s["storm_solar_max"][cname] = {"n": int(m.sum()), "median_db": {
                k: float(np.nanmedian(20 * np.log10(v[m] / fn_obs[m]))) for k, v in models.items()}}
    s["splits"] = splits(fn_obs, models, hs, rz, J["lat"], J["lt"], J["quality"])
    if goal_samples is not None:
        s["goal_sensitivity"] = goal_sensitivity(goal_samples, s["bands"]["solar_max"])
    (out / "summary.json").write_text(json.dumps(_jsonable(s), indent=1), encoding="utf-8")

    L = ["# ISIS-1 against epoch-driven IRI-2020, 900-1400 km", "",
         f"{s['n']} CRC profiles (quality 0-7): " + ", ".join(
             f"{k} {v}" for k, v in s["n_by_regime"].items()) + ".", "",
         "Median model/measured fN_sat in dB; 5th percentiles in MHz; the ensemble edge is the "
         "minimum of the three epoch topsides' 5th percentiles; realised coverage is the fraction "
         "of measured fN_sat at or above it.", ""]
    short = {"binned_NeQuick": "binned NQ", "epoch_NeQuick_quiet": "NQ quiet",
             "epoch_NeQuick": "NQ", "epoch_NeQuick_allstorm": "NQ all storms",
             "epoch_IRI2001": "IRI-2001", "epoch_IRI2001corr": "IRI-2001c"}
    for rname, rows in s["bands"].items():
        if not rows:
            continue
        L += [f"## {rname}", "",
              "| band | n | " + " | ".join(short[k] for k in names) +
              " | p05 obs | p05 NQ / 2001 / 2001c | ensemble | ens/obs | realised coverage | "
              "H0 fit / IRI law (p10, p50, p90) |",
              "|---" * (len(names) + 8) + "|"]
        for b, row in rows.items():
            p = row["p05_mhz"]
            h0 = row.get("h0_fit_over_iri_law")
            L.append(f"| {b} | {row['n']} | " + " | ".join(
                f"{row['median_db'][k]:+.2f}" for k in names) +
                f" | {p['observed']:.3f} | {p['epoch_NeQuick']:.3f} / {p['epoch_IRI2001']:.3f} / "
                f"{p['epoch_IRI2001corr']:.3f} | {row['ensemble_p05_mhz']:.3f} | "
                f"{row['ensemble_over_observed_p05']:.3f} | "
                f"{100 * row['realised_coverage_of_ensemble_edge']:.1f} % | " +
                (", ".join(f"{v:.2f}" for v in h0) if h0 else "-") + " |")
        L.append("")
    if below is not None and (below / "summary_epoch.json").is_file():
        B = json.loads((below / "summary_epoch.json").read_text())
        L += ["## For comparison: 500-900 km (2026-10-18 isis-epoch run)", "",
              "| regime | band | n | NQ epoch dB | p05 obs | p05 NQ epoch |", "|---|---|---|---|---|---|"]
        for rname in ("solar_max", "solar_min"):
            for b, row in B["bands"].get(rname, {}).items():
                L.append(f"| {rname} | {b} | {row['n']} | "
                         f"{row['median_db']['iri_epoch_NeQuick']:+.2f} | "
                         f"{row['p05_mhz']['observed']:.3f} | "
                         f"{row['p05_mhz']['iri_epoch_NeQuick']:.3f} |")
        L.append("")
    L += ["## Where the tail is overstated", "",
          "| group | n | p05 obs | ensemble p05 | ens/obs | realised coverage | NQ median dB |",
          "|---|---|---|---|---|---|---|"]
    for k, v in s["splits"].items():
        if isinstance(v, dict):
            L.append(f"| {k} | {v['n']} | {v['p05_observed_mhz']:.3f} | {v['ensemble_p05_mhz']:.3f} | "
                     f"{v['ensemble_over_observed_p05']:.3f} | "
                     f"{100 * v['realised_coverage_of_ensemble_edge']:.1f} % | "
                     f"{v['median_db_nequick']:+.2f} |")
    L.append("")
    if "goal_sensitivity" in s:
        L += ["## Phase 4 sensitivity: coverage at 0.33 MHz with the tail ratio applied (bound)", "",
              "| cell | tail ratio | coverage | bound | p05 MHz | bound |", "|---|---|---|---|---|---|"]
        for k, v in s["goal_sensitivity"].items():
            L.append(f"| {k} | {v['tail_ratio']:.3f} | {100 * v['coverage_at_033']:.1f} % | "
                     f"{100 * v['coverage_at_033_bound']:.1f} % | {v['p05_mhz']:.3f} | "
                     f"{v['p05_mhz_bound']:.3f} |")
        L.append("")
    L += ["## Storm activity, solar maximum, 900-1400 km", "", "```",
          json.dumps(_jsonable(s["storm_solar_max"]), indent=1), "```", ""]
    (out / "summary.md").write_text("\n".join(L), encoding="utf-8")
    print("\n".join(L))
    return s


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    rp = sub.add_parser("run")
    rp.add_argument("out_dir", type=Path)
    rp.add_argument("--crc", type=Path, required=True)
    rp.add_argument("--matlab", default="matlab")
    rp.add_argument("--reuse-iri", action="store_true")
    an = sub.add_parser("analyse")
    an.add_argument("out_dir", type=Path)
    an.add_argument("--below", type=Path, default=None)
    an.add_argument("--goal-run", type=Path, default=None,
                    help="band_edge_goal_1400 output, for the Phase 4 sensitivity bound")
    args = ap.parse_args(argv)
    if args.cmd == "run":
        cmd_run(args.out_dir, args.crc, args.matlab, args.reuse_iri)
    else:
        cmd_analyse(args.out_dir, args.below, args.goal_run)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
