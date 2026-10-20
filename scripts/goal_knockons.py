"""WP0 Phase 4: SCI-02's goal at 1000-1400 km, and what it costs the instrument.

    python scripts/goal_knockons.py GOAL_RUN OUT

GOAL_RUN is the output of ``tss study min_frequency`` with the
``band_edge_goal_1400`` plan.  Writes ``goal_altitudes.csv`` (one row per
altitude, class and solar level), ``goal_knockons.json`` and ``goal_knockons.md``
to OUT.

Per cell, from the run's summary and ensemble:

* the 95 and 99 per cent band edge, combined across the three topsides at the
  minimum of the per-model percentiles (as INS-01a);
* the profile a sweep starting at INS-01a's 0.33 MHz retrieves: coverage and
  altitude completeness, taking the worst topside;
* IRI's O+/H+ transition height (5th, 50th, 95th percentile) and the fraction
  of samples in which it lies below the spacecraft, which is what the goal is
  for.

Then the knock-ons of flying at the goal's top altitude, from the Rev G
derivations with h_sat in place of 900 km:

* INS-04 / SCI-08: 1.5 x (h_sat - 200 km), the hmF2 floor of SCI-02;
* INS-06: f_s >= 4 K_f R_max / c, at A1's peak rate and A2's 37.3 MHz/s;
* INS-03 and INS-13: A1 steps from the band edge to 20 MHz, 10 km range bins
  over R_max, 8 bit, at the A1 sweep time and INS-07's 6.7 s cadence;
* the two-way spreading term against 900 km, 20 log10 of the true-range ratio;
* INS-20: the half-wave length at the goal's band edge;
* INS-10, provisional: the link-budget SNR per range cell over the band with
  the true range to a 300 km peak, at 900 and at h_sat (the model's 7 dB drop
  on PHaRLAP is not yet traced, so this is a difference, not a verdict).
"""

from __future__ import annotations

import csv
import json
import math
import sys
from collections import defaultdict
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

C = 299_792_458.0
INS01A_HZ = 0.33e6
HMF2_FLOOR_KM = 200.0       # SCI-02's hmF2 envelope floor, the INS-04 basis
REF_ALT_KM = 900.0          # top of MIS-01, the Rev G basis
A2_RATE = 37.3e6            # Hz/s, configuration A2 (INS-06)
INS07_CADENCE_S = 6.7
BIN_KM = 10.0               # INS-13 product range bins


def _f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def cell_rows(run: Path) -> list[dict]:
    summ = list(csv.DictReader(open(run / "minfreq_summary.csv")))
    ens = {(r["orbit_class"], _f(r["altitude_km"]), r["solar"]): r
           for r in csv.DictReader(open(run / "ensemble.csv"))}
    groups = defaultdict(list)
    for r in summ:
        groups[(r["orbit_class"], _f(r["altitude_km"]), r["solar"])].append(r)
    out = []
    for key in sorted(groups, key=lambda k: (k[1], k[0], k[2])):
        rs = groups[key]
        e = ens[key]
        neq = next((r for r in rs if r["topside"] == "NeQuick"), rs[0])
        oh50 = [_f(r["oh_transition_p50_km"]) for r in rs]
        out.append(dict(
            orbit_class=key[0], altitude_km=key[1], solar=key[2], n_models=len(rs),
            f_min_95_MHz=_f(e["f_min_95_ensemble_MHz"]),
            f_min_99_MHz=_f(e["f_min_99_ensemble_MHz"]),
            governing_model_95=e["governing_model_95"],
            coverage_at_033=min(_f(r["coverage_at_ins01a"]) for r in rs),
            mean_alt_completeness_at_033=min(_f(r["mean_alt_completeness_at_ins01a"]) for r in rs),
            p05_alt_completeness_at_033=min(_f(r["p05_alt_completeness_at_ins01a"]) for r in rs),
            hmf2_median_km=_f(neq["hmf2_median_km"]),
            hmf2_p001_km=min(_f(r["hmf2_p001_km"]) for r in rs),
            oh_p05_km=_f(neq["oh_transition_p05_km"]),
            oh_p50_km=_f(neq["oh_transition_p50_km"]),
            oh_p95_km=_f(neq["oh_transition_p95_km"]),
            oh_valid=_f(neq["frac_oh_transition_valid"]),
            oh_below_hsat=_f(neq["frac_oh_transition_below_hsat"]),
            oh_p50_spread_km=float(np.nanmax(oh50) - np.nanmin(oh50)) if np.isfinite(oh50).any()
            else float("nan"),
        ))
    return out


def knockons(h_sat_km: float, f_edge_hz: float, hmf2_median_km: float) -> dict:
    from tss.sweep_law import active_duration_s, step_count, stepped_schedule

    r_true = h_sat_km - HMF2_FLOOR_KM
    r_max = 1.5 * r_true
    r_max_ref = 1.5 * (REF_ALT_KM - HMF2_FLOOR_KM)
    out = {"h_sat_km": h_sat_km, "ins04_r_max_km": r_max, "ins04_rev_g_km": r_max_ref}
    a1 = stepped_schedule(INS01A_HZ, 20e6, 26.8e-3)
    k_a1 = a1.k_f_hz_per_s
    for name, k in (("a1", k_a1), ("a2", A2_RATE)):
        out[f"ins06_fs_{name}_kHz"] = 4 * k * r_max * 1e3 / C / 1e3
        out[f"ins06_fs_{name}_rev_g_kHz"] = 4 * k * r_max_ref * 1e3 / C / 1e3
    bins = math.ceil(r_max / BIN_KM)
    bins_ref = math.ceil(r_max_ref / BIN_KM)
    for tag, f0 in (("at_033", INS01A_HZ), ("at_goal_edge", f_edge_hz)):
        sch = stepped_schedule(f0, 20e6, 26.8e-3)
        n = step_count(sch)
        t_sweep = active_duration_s(sch)
        size = n * bins
        out[f"ins03_steps_{tag}"] = n
        out[f"sweep_s_{tag}"] = t_sweep
        out[f"ins13_ionogram_kB_{tag}"] = size / 1e3
        out[f"ins13_MB_per_day_sweep_{tag}"] = size * 86400 / t_sweep / 1e6
        out[f"ins13_MB_per_day_6p7s_{tag}"] = size * 86400 / INS07_CADENCE_S / 1e6
    n_ref = step_count(a1)
    out["ins13_rev_g_ionogram_kB"] = n_ref * bins_ref / 1e3
    out["ins13_rev_g_MB_per_day_sweep"] = n_ref * bins_ref * 86400 / active_duration_s(a1) / 1e6
    out["range_bins"] = bins
    out["spreading_db_vs_900_at_floor"] = 20 * math.log10(r_true / (REF_ALT_KM - HMF2_FLOOR_KM))
    out["spreading_db_vs_900_at_median_hmf2"] = 20 * math.log10(
        (h_sat_km - hmf2_median_km) / (REF_ALT_KM - hmf2_median_km))
    out["ins20_half_wave_m_at_goal_edge"] = C / (2 * f_edge_hz)
    out["ins20_half_wave_m_at_033"] = C / (2 * INS01A_HZ)
    out["ins10_provisional"] = ins10(h_sat_km)
    return out


def ins10(h_sat_km: float, hmf2_km: float = 300.0) -> dict:
    from dataclasses import replace

    from tss.l2_instrument import linkbudget as lb

    f = np.concatenate([np.linspace(0.33, 1.0, 40), np.linspace(1.0, 20.0, 200)[1:]])
    res = {}
    for tag, h in (("h900", REF_ALT_KM), ("hsat", h_sat_km)):
        cfg = replace(lb.SounderConfig(), h_sat=h * 1e3, true_range=(h - hmf2_km) * 1e3)
        snr = np.asarray(lb.snr_vs_frequency(f, cfg), dtype=float)
        ok = np.isfinite(snr)
        res[tag] = {"snr_median_db": float(np.median(snr[ok])), "snr_min_db": float(np.min(snr[ok])),
                    "f_at_min_MHz": float(f[ok][np.argmin(snr[ok])]), "n_valid": int(ok.sum())}
    res["median_change_db"] = res["hsat"]["snr_median_db"] - res["h900"]["snr_median_db"]
    res["min_change_db"] = res["hsat"]["snr_min_db"] - res["h900"]["snr_min_db"]
    res["note"] = ("tss.l2_instrument.linkbudget with its analytic ionosphere; true range to a "
                   f"{hmf2_km:.0f} km peak.  Provisional until the INS-10 drop on PHaRLAP is traced.")
    return res


def main(argv: list[str]) -> int:
    if len(argv) != 2:
        print(__doc__)
        return 2
    run, out = Path(argv[0]), Path(argv[1])
    out.mkdir(parents=True, exist_ok=True)
    cells = cell_rows(run)
    with open(out / "goal_altitudes.csv", "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(cells[0].keys()))
        w.writeheader()
        w.writerows(cells)
    top = max(c["altitude_km"] for c in cells)
    by_alt = {}
    for a in sorted({c["altitude_km"] for c in cells}):
        cs = [c for c in cells if c["altitude_km"] == a]
        gov = min(cs, key=lambda c: c["f_min_95_MHz"])
        worst_cov = min(cs, key=lambda c: c["coverage_at_033"])
        by_alt[a] = dict(
            governing=f"{gov['orbit_class']}/{gov['solar']}", f_min_95_MHz=gov["f_min_95_MHz"],
            f_min_99_MHz=min(c["f_min_99_MHz"] for c in cs),
            worst_coverage_at_033=worst_cov["coverage_at_033"],
            worst_coverage_cell=f"{worst_cov['orbit_class']}/{worst_cov['solar']}",
            worst_mean_alt_completeness_at_033=min(c["mean_alt_completeness_at_033"] for c in cs),
            worst_p05_alt_completeness_at_033=min(c["p05_alt_completeness_at_033"] for c in cs),
            oh_below_hsat_min=min(c["oh_below_hsat"] for c in cs),
            oh_below_hsat_max=max(c["oh_below_hsat"] for c in cs),
        )
    cs_top = [c for c in cells if c["altitude_km"] == top]
    edge = min(c["f_min_95_MHz"] for c in cs_top) * 1e6
    hm_med = float(np.median([c["hmf2_median_km"] for c in cs_top]))
    k = knockons(top, edge, hm_med)
    result = {"by_altitude": by_alt, "knockons_at_top": k,
              "oh_p50_topside_spread_km_max": float(np.nanmax([c["oh_p50_spread_km"] for c in cells]))}
    (out / "goal_knockons.json").write_text(json.dumps(result, indent=1), encoding="utf-8")

    L = ["# SCI-02 goal: 1000-1400 km", "",
         "Band edge combined at the minimum of the per-topside percentiles; coverage and "
         "completeness at a 0.33 MHz start take the worst topside; O+/H+ from NeQuick runs.", "",
         "| h_sat | class | solar | 95 % edge MHz | 99 % edge | coverage at 0.33 | mean alt. "
         "completeness | p05 | hmF2 median | O+/H+ p05 | p50 | p95 km | below h_sat |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for c in cells:
        L.append(f"| {c['altitude_km']:.0f} | {c['orbit_class']} | {c['solar']} | "
                 f"{c['f_min_95_MHz']:.3f} | {c['f_min_99_MHz']:.3f} | {100 * c['coverage_at_033']:.1f} % | "
                 f"{c['mean_alt_completeness_at_033']:.3f} | {c['p05_alt_completeness_at_033']:.3f} | "
                 f"{c['hmf2_median_km']:.0f} | {c['oh_p05_km']:.0f} | {c['oh_p50_km']:.0f} | "
                 f"{c['oh_p95_km']:.0f} | {100 * c['oh_below_hsat']:.0f} % |")
    L += ["", "## By altitude", "", "```", json.dumps(by_alt, indent=1), "```", "",
          f"## Knock-ons at {top:.0f} km", "", "```", json.dumps(k, indent=1), "```", ""]
    (out / "goal_knockons.md").write_text("\n".join(L), encoding="utf-8")
    print(json.dumps({"by_altitude": by_alt}, indent=1))
    print(json.dumps(k, indent=1))
    return 0


if __name__ == "__main__":
    raise SystemExit(main(sys.argv[1:]))
