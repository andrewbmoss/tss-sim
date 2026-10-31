"""Layer: studies
WP: WP4
Claims: C1, C2
Status: implemented
Requirements: INS-05, INS-07, INS-10, SCI-05, SCI-07, SCI-09

WP4 Phase 3: the comparison of the dispersion and step-transition treatments.

Decisions of 2026-10-30: score at the echo level, independent of the trace
extractor, and through the chain with clutter off; option 3 (narrower steps)
only where the smear needs it; one treatment chosen band-wide for each of
dispersion (gate G2) and the step transition (G3); the matched option is the
hybrid matched filter.

THE TREATMENTS
--------------
Dispersion, each with the transition accepted:

* ``D0`` one transform of the whole dwell (the baseline);
* ``D1x2``, ``D1x4`` the dwell cut into 2 or 4 sub-blocks added in power;
* ``D2`` the rate search, beat chirp rates of ``j / T^2``, ``j`` = -2..32, the
  rate with the highest peak kept for each column;
* ``D3`` narrower steps: a 50 kHz cap on the steps that start between 2 and
  8 MHz, sized so that the B2 reference (H = 150 km, f_N at the spacecraft
  2 MHz) crosses about two cells at 3 MHz; the 150 kHz cap elsewhere.

Step transition, each with the D0 transform:

* ``T0`` accept the skirt (the 2026-09-30 default; ``D0`` itself);
* ``T1`` gate the first ``2 R_max / c`` (7 ms) of every dwell;
* ``T2`` the hybrid matched filter in the steps whose rate changes.

ECHO-LEVEL MEASURES
-------------------
Each scene is run once per magnetoionic mode, so every column holds one echo
and nothing else (clutter off).  Each treatment is applied to the record with
in-dwell dispersion modelled, and the reference is the same record with
dispersion switched off in L2 (same noise draw) under ``D0``: the loss of a
treatment therefore includes what it costs an undispersed echo (the coherent
gain a sub-block gives up) as well as what it recovers.  Per truth echo:

* ``snr_db`` the peak within +-150 km of the echo's range at the step centre,
  over the detector's noise reference;
* ``loss_db`` reference SNR minus ``snr_db``;
* ``width_km`` the contiguous extent within 3 dB of that peak;
* ``edge_err_km`` the first detected cell in the window minus the echo's
  range (negative: early);
* ``detected`` whether the detector declared any cell in the window;
* ``cell_km`` the column's range resolution (INS-05);
* for the rate search, ``dpdf_err`` the relative error of the dP'/df it
  returns against the truth's.

THE SCENES
----------
The six Chapman ionospheres of :mod:`tss.studies.wp4_campaign` and a sample of
ISIS-1 hand-scaled topside profiles (Jackson's method, the CRC archive), each
completed below its lowest point by an alpha-Chapman peak placed at the
record's IRI hmF2 (or 30 km below the lowest point if that is lower), with
NmF2 the larger of 1.15 times the lowest density and the record's IRI NmF2, and
the scale height that makes the two meet.  The spacecraft is at the record's
altitude, place and time.

The chain-level part runs :mod:`tss.studies.wp4_campaign` (clutter off) once
per treatment.
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
WINDOW_M = 150.0e3
BANDS_MHZ = ((0.0, 0.5), (0.5, 2.0), (2.0, 5.0), (5.0, 10.0), (10.0, 20.0))
DEFAULT_SEEDS = (1, 2, 3)

#: name: (L3 keyword arguments, sweep overrides)
NARROW = {"narrow_lo_hz": 2.0e6, "narrow_hi_hz": 8.0e6, "narrow_cap_hz": 50.0e3}
TREATMENTS: dict[str, tuple[dict[str, Any], dict[str, Any]]] = {
    "D0": ({}, {}),
    "D1x2": ({"dispersion": "blocks", "dispersion_blocks": 2}, {}),
    "D1x4": ({"dispersion": "blocks", "dispersion_blocks": 4}, {}),
    "D2": ({"dispersion": "rate_search"}, {}),
    "D3": ({}, NARROW),
    "T1": ({"transition": "gate"}, {}),
    "T2": ({"transition": "matched"}, {}),
}
GROUP = {"D0": "both", "D1x2": "G2", "D1x4": "G2", "D2": "G2", "D3": "G2", "T1": "G3", "T2": "G3"}


# -- ISIS-1 scenes -----------------------------------------------------------------


class MeasuredTopside:
    """An L0 tier from one ISIS-1 profile, completed below by an alpha-Chapman peak."""

    def __init__(self, profile: Any) -> None:
        from tss.l0_ionosphere.chapman import ChapmanIonosphere

        alt = np.asarray(profile.alt_km, dtype=float) * 1e3
        ne = np.asarray(profile.ne_m3, dtype=float)
        order = np.argsort(alt)
        self.alt, self.ln = alt[order], np.log(ne[order])
        h_low, n_low = float(self.alt[0]), float(ne[order][0])
        ex = dict(getattr(profile, "extra", {}) or {})
        hm = min(float(ex.get("hmf2_iri_km", 300.0)) * 1e3, h_low - 30.0e3)
        fo = float(ex.get("fof2_iri_mhz", 0.0) or 0.0) * 1e6
        nm = max(1.15 * n_low, fo**2 / 80.6)
        # Solve n_low = nm exp(0.5 (1 - z - exp(-z))) for z > 0.
        target = 2.0 * math.log(n_low / nm)  # = 1 - z - exp(-z), negative
        lo, hi = 1e-6, 50.0
        for _ in range(200):
            mid = 0.5 * (lo + hi)
            if 1.0 - mid - math.exp(-mid) > target:
                lo = mid
            else:
                hi = mid
        z = 0.5 * (lo + hi)
        h_ch = (h_low - hm) / z
        self.cap = ChapmanIonosphere(nmf2=nm, hmf2_m=hm, plasma_scale_height_m=2.0 * h_ch)
        self.h_low = h_low
        top = slice(-2, None) if self.alt.size >= 2 else slice(None)
        da = float(np.diff(self.alt[top])[0]) if self.alt.size >= 2 else 1.0
        self.top_slope = float(np.diff(self.ln[top])[0]) / da if self.alt.size >= 2 else -1.0 / 300e3
        self.label = (f"ISIS-1 {profile.epoch:%Y-%m-%d %H:%M} {profile.lat_deg:.1f},{profile.lon_deg:.1f} "
                      f"sat {profile.sat_alt_km:.0f} km + Chapman cap (hmF2 {hm / 1e3:.0f} km)")

    def describe(self) -> dict[str, Any]:
        return {"name": "measured_topside", "label": self.label}

    def state_at(self, lat_deg: float, lon_deg: float, epoch: str, alt_m: np.ndarray) -> Any:
        from tss.l0_ionosphere.chapman import chapman_density

        st = self.cap.state_at(lat_deg, lon_deg, epoch, alt_m)
        a = np.asarray(st.alt_m, dtype=float)
        ln = np.interp(a, self.alt, self.ln)
        above = a > self.alt[-1]
        ln[above] = self.ln[-1] + self.top_slope * (a[above] - self.alt[-1])
        below = a < self.h_low
        cap_ne = chapman_density(a, self.cap.nmf2, self.cap.hmf2_m, self.cap.plasma_scale_height_m)
        ne = np.where(below, cap_ne, np.exp(ln))
        from tss.l0_ionosphere.chapman import collision_frequency

        return replace(st, ne=ne, nu=collision_frequency(ne, np.asarray(st.te, dtype=float)), source=self.label,
                       extra={**dict(st.extra), "fof2_hz": float(np.sqrt(80.6 * np.max(ne))), "measured_topside": True})


def isis_scenes(n: int, seed: int = 0) -> list[Any]:
    """A stratified sample of ISIS-1 profiles: latitude band x day/night x solar activity."""
    from tss.data import spdf_isis

    path = spdf_isis.local_ne_path()
    if not path:
        return []
    cells: dict[tuple[int, int, int], list[Any]] = {}
    for p in spdf_isis.iter_profiles(path):
        alt = np.asarray(p.alt_km, dtype=float)
        if getattr(p, "quality", 0) < 6 or alt.size < 6 or not 500.0 <= p.sat_alt_km <= 1100.0:
            continue
        ex = dict(getattr(p, "extra", {}) or {})
        lat = abs(float(p.lat_deg))
        band = 0 if lat < 20 else (1 if lat < 55 else 2)
        day = 1 if 6.0 <= float(ex.get("lt_hour", 12.0)) < 18.0 else 0
        hi_sun = 1 if float(ex.get("rz12", 0.0)) >= 60 else 0
        cells.setdefault((band, day, hi_sun), []).append(p)
    rng = np.random.default_rng(seed)
    keys = sorted(cells)
    out: list[Any] = []
    while len(out) < n and any(cells[k] for k in keys):
        for k in keys:
            if cells[k] and len(out) < n:
                out.append(cells[k].pop(int(rng.integers(len(cells[k])))))
    return out


# -- one scene ---------------------------------------------------------------------


def _sweep(config: dict[str, Any], over: dict[str, Any]) -> tuple[Any, np.ndarray, dict[str, Any]]:
    from tss.studies import _grids, _sweep_from_config

    cfg = copy.deepcopy(config)
    cfg["sweep"] = {**dict(cfg.get("sweep", {})), **over}
    sweep, notes = _sweep_from_config(cfg, io.StringIO())
    _, freq = _grids(cfg)
    return sweep, freq, notes


def _instrument(config: dict[str, Any], dispersion: bool) -> Any:
    from tss.l2_instrument.baseline import BaselineInstrument
    from tss.studies import instrument_config

    icfg = instrument_config(config)
    icfg["in_dwell_dispersion"] = bool(dispersion)
    inst = BaselineInstrument.from_config(icfg)
    assert inst.in_dwell_dispersion == bool(dispersion)
    return inst


def _column_for(spectra: list[Any], f: float) -> int | None:
    c = np.array([s.centre_frequency_hz for s in spectra])
    j = int(np.argmin(np.abs(c - f)))
    return j if abs(c[j] - f) < 1.0 else None


def echo_measures(spectra: list[Any], echoes: Any, l3: Any, ref_snr: dict[float, float] | None,
                  dpdf_truth: dict[float, float] | None) -> list[dict[str, Any]]:
    """Per-truth-echo measures for one single-mode record (module docstring)."""
    from tss.l3_signal.detect import cfar_detect, noise_reference
    from tss.l5_assessment.scoring import truth_echo_cells

    cells = truth_echo_cells(echoes)
    rows = []
    for f, r, cut in zip(cells["freq_hz"], cells["range_m"], cells["cutoff_step"], strict=True):
        j = _column_for(spectra, float(f))
        if j is None:
            continue
        sp = spectra[j]
        if r > sp.unambiguous_range_m:
            continue
        rng = np.asarray(sp.range_m, dtype=float)
        p = np.asarray(sp.power_db, dtype=float)[0]
        win = np.abs(rng - r) <= WINDOW_M
        if not np.any(win):
            continue
        floor = float(noise_reference(sp)["mean_db"])
        idx = np.flatnonzero(win)
        k = int(idx[np.argmax(p[idx])])
        snr = float(p[k] - floor)
        lo = k
        while lo > 0 and p[lo - 1] >= p[k] - 3.0:
            lo -= 1
        hi = k
        while hi < p.size - 1 and p[hi + 1] >= p[k] - 3.0:
            hi += 1
        det = cfar_detect(sp, **l3.detector_settings())
        hits = np.flatnonzero(det.mask & win)
        row = {
            "f_mhz": float(f) / 1e6, "range_km": float(r) / 1e3, "cutoff_step": bool(cut),
            "snr_db": snr, "width_km": float(rng[hi] - rng[lo]) / 1e3 + float(sp.bin_spacing_m) / 1e3,
            "detected": bool(hits.size), "edge_err_km": float(rng[hits[0]] - r) / 1e3 if hits.size else float("nan"),
            "cell_km": float(sp.range_resolution_m) / 1e3,
            "loss_db": (ref_snr.get(_key(float(f)), float("nan")) - snr) if ref_snr else float("nan"),
        }
        rate = getattr(sp, "rate_hz_per_s2", None)
        if rate is not None and dpdf_truth:
            # K at this column from the range scale: R = c f_b / (2 K).
            kk = C_LIGHT * float(sp.beat_hz[1] - sp.beat_hz[0]) / (2.0 * float(sp.range_m[1] - sp.range_m[0]))
            est = float(rate[k]) * C_LIGHT / (2.0 * kk**2)
            tru = dpdf_truth.get(_key(float(f)), float("nan"))
            row["dpdf_err"] = (est - tru) / tru if tru and np.isfinite(tru) else float("nan")
        rows.append(row)
    return rows


def _dpdf_truth(echoes: Any) -> dict[float, float]:
    from tss.l5_assessment.scoring import truth_echo_cells

    c = truth_echo_cells(echoes)
    f = np.asarray(c["freq_hz"], dtype=float)
    r = np.asarray(c["range_m"], dtype=float)
    if f.size < 3:
        return {}
    o = np.argsort(f)
    g = np.gradient(r[o], f[o])
    return {_key(float(fi)): float(gi) for fi, gi in zip(f[o], g, strict=True)}


def _key(f_hz: float) -> int:
    return int(round(f_hz / 100.0))


def run_scene(name: str, l0: Any, platform: Any, epoch: str, *, config: dict[str, Any], seeds: list[int],
              treatments: dict[str, tuple[dict[str, Any], dict[str, Any]]]) -> list[dict[str, Any]]:
    from tss.l3_signal.ionogram import BaselineSignalProcessing
    from tss.registry import build
    from tss.studies import _grids

    alt_m, _ = _grids(config)
    state = l0.state_at(platform.lat_deg, platform.lon_deg, epoch, alt_m)
    inst_on, inst_off = _instrument(config, True), _instrument(config, False)
    rows: list[dict[str, Any]] = []
    by_sweep: dict[str, list[str]] = {}
    for t, (_, sw) in treatments.items():
        by_sweep.setdefault(json.dumps(sw, sort_keys=True), []).append(t)
    for sw_key, names in by_sweep.items():
        sweep, freq, notes = _sweep(config, json.loads(sw_key))
        for mode in ("O", "X"):
            l1 = build("l1", "analytic", modes=[mode])
            echoes = l1.echoes(state, platform, freq, np.random.default_rng(0))
            truth_dpdf = _dpdf_truth(echoes)
            for seed in seeds:
                rec_on = inst_on.receive(echoes, sweep, np.random.default_rng(seed))
                rec_off = inst_off.receive(echoes, sweep, np.random.default_rng(seed))
                l3_ref = BaselineSignalProcessing()
                ref_rows = echo_measures(l3_ref.spectra(rec_off), echoes, l3_ref, None, None)
                ref = {_key(r["f_mhz"] * 1e6): r["snr_db"] for r in ref_rows}
                for tname in names:
                    kw, _ = treatments[tname]
                    l3 = BaselineSignalProcessing(**kw)
                    for r in echo_measures(l3.spectra(rec_on), echoes, l3, ref, truth_dpdf):
                        rows.append({"scene": name, "mode": mode, "seed": seed, "treatment": tname,
                                     "n_steps": int(notes.get("n_steps", 0)), **r})
    return rows


# -- aggregation -----------------------------------------------------------------


def _band(f: float) -> str:
    for lo, hi in BANDS_MHZ:
        if lo <= f < hi:
            return f"{lo:g}-{hi:g}"
    return "other"


def _stats(rows: list[dict[str, Any]]) -> dict[str, Any]:
    def med(k: str) -> float:
        v = np.array([r[k] for r in rows if k in r], dtype=float)
        return float(np.nanmedian(v)) if np.any(np.isfinite(v)) else float("nan")

    def pct(k: str, q: float) -> float:
        v = np.array([r[k] for r in rows if k in r], dtype=float)
        return float(np.nanpercentile(v, q)) if np.any(np.isfinite(v)) else float("nan")

    det = np.array([r["detected"] for r in rows], dtype=float)
    return {
        "n": len(rows),
        "loss_db_median": med("loss_db"), "loss_db_p90": pct("loss_db", 90),
        "snr_db_median": med("snr_db"), "snr_db_p10": pct("snr_db", 10),
        "width_km_median": med("width_km"), "width_km_p90": pct("width_km", 90),
        "edge_err_km_median": med("edge_err_km"), "edge_abs_km_p90": float(np.nanpercentile(
            np.abs([r["edge_err_km"] for r in rows]), 90)) if rows else float("nan"),
        "detected": float(np.mean(det)) if det.size else float("nan"),
        "cell_km_max": float(np.nanmax([r["cell_km"] for r in rows])) if rows else float("nan"),
        "dpdf_err_median_abs": float(np.nanmedian(np.abs([r.get("dpdf_err", np.nan) for r in rows])))
        if any("dpdf_err" in r for r in rows) else float("nan"),
    }


def aggregate(rows: list[dict[str, Any]]) -> dict[str, Any]:
    out: dict[str, Any] = {}
    for t in TREATMENTS:
        mine = [r for r in rows if r["treatment"] == t]
        if not mine:
            continue
        blk = {"all": _stats(mine)}
        for lo, hi in BANDS_MHZ:
            sel = [r for r in mine if lo <= r["f_mhz"] < hi]
            if sel:
                blk[f"{lo:g}-{hi:g}"] = _stats(sel)
        # Near the cusp: each (scene, mode, seed) record's three highest echoes
        # (cutoff steps included); near f_N: its three lowest.
        keyf = {}
        for r in mine:
            keyf.setdefault((r["scene"], r["mode"], r["seed"]), []).append(r)
        cusp, low = [], []
        for v in keyf.values():
            v = sorted(v, key=lambda x: x["f_mhz"])
            low += v[:3]
            cusp += v[-3:]
        blk["near_fN"] = _stats(low)
        blk["near_cusp"] = _stats(cusp)
        out[t] = blk
    return out


def sweep_costs(config: dict[str, Any]) -> dict[str, Any]:
    from tss import sweep_law

    out = {}
    for t, (_, sw) in TREATMENTS.items():
        sweep, _, _ = _sweep(config, sw)
        dur = float(sweep_law.active_duration_s(sweep))
        out[t] = {"n_steps": int(sweep_law.step_count(sweep)), "sweep_s": dur,
                  "along_track_km_500": dur * 7.06, "along_track_km_900": dur * 6.42}
    return out


def chain(config: dict[str, Any], seeds: list[int]) -> dict[str, Any]:
    from tss.studies import wp4_campaign as w

    out = {}
    for t, (kw, sw) in TREATMENTS.items():
        cfg = copy.deepcopy(config)
        cfg["tiers"] = {**cfg["tiers"], "l3": {"name": "baseline", **kw}}
        cfg["sweep"] = {**dict(cfg.get("sweep", {})), **sw}
        cfg["wp4_campaign"] = {"seeds": seeds, "clutter": "off"}
        res = w.compute(cfg)
        a = res["aggregate"]["overall"]
        ins05 = [r["ins05_m"] for r in res["rows"]]
        out[t] = {k: a[k] for k in ("sci05_median", "sci05_worst", "sci07_median_m", "sci07_worst_m",
                                     "ins10_median_db", "ins10_worst_db", "completeness_o_median",
                                     "sci05_verdicts", "sci07_verdicts", "ins10_verdicts")}
        out[t]["ins05_worst_m"] = float(np.nanmax(ins05))
        out[t]["ins05_verdicts"] = {v: sum(1 for r in res["rows"] if r["ins05_verdict"] == v)
                                    for v in sorted({r["ins05_verdict"] for r in res["rows"]})}
        out[t]["run_s_median"] = float(np.median([r["run_s"] for r in res["rows"]]))
    return out


def compute(config: dict[str, Any]) -> dict[str, Any]:
    from tss.config import load_yaml, resolve_path
    from tss.studies import _platform_for
    from tss.studies import wp4_campaign as w

    opts = dict(config.get("wp4_dispersion", {}) or {})
    seeds = [int(s) for s in (opts.get("seeds") or DEFAULT_SEEDS)]
    n_isis = int(opts.get("n_isis", 12))
    do_chain = bool(opts.get("chain", True))
    conf = load_yaml(resolve_path("leo_500km_single", subdir="configurations"))
    rows: list[dict[str, Any]] = []
    scenes: list[dict[str, Any]] = []
    for name in (opts.get("profiles") or list(w.PROFILES)):
        lat, lt, *_ = w.PROFILES[name]
        pl = replace(_platform_for(conf, w.EPOCH), lat_deg=float(lat), lon_deg=w._lon_for_local_time(lt))
        rows += run_scene(name, w._l0_for(name), pl, w.EPOCH, config=config, seeds=seeds, treatments=TREATMENTS)
        scenes.append({"scene": name, "source": "chapman"})
    for i, prof in enumerate(isis_scenes(n_isis)):
        l0 = MeasuredTopside(prof)
        ep = prof.epoch.strftime("%Y-%m-%dT%H:%M:%SZ")
        pl = replace(_platform_for(conf, ep), alt_m=float(prof.sat_alt_km) * 1e3, lat_deg=float(prof.lat_deg),
                     lon_deg=float(prof.lon_deg), epoch=ep)
        name = f"isis{i:02d}"
        try:
            rows += run_scene(name, l0, pl, ep, config=config, seeds=seeds, treatments=TREATMENTS)
            scenes.append({"scene": name, "source": l0.label})
        except Exception as exc:  # noqa: BLE001
            scenes.append({"scene": name, "source": l0.label, "error": str(exc)[:200]})
    out = {"seeds": seeds, "scenes": scenes, "treatments": {k: {"l3": v[0], "sweep": v[1]} for k, v in TREATMENTS.items()},
           "rows": rows, "echo": aggregate(rows), "sweep_costs": sweep_costs(config)}
    if do_chain:
        out["chain"] = chain(config, seeds)
    return out


def _f(x: Any, nd: int = 1) -> str:
    return "nan" if x is None or (isinstance(x, float) and not math.isfinite(x)) else f"{x:.{nd}f}"


def _summary(res: dict[str, Any]) -> tuple[str, dict[str, Any]]:
    e = res["echo"]
    lines = ["# wp4_dispersion", "",
             f"{len(res['scenes'])} scenes x O, X x {len(res['seeds'])} seeds; echo-level, clutter off.", "",
             "| Treatment | Loss dB med (p90) | Width km med (p90) | Edge err km med | Detected | Cell km max"
             " | Loss near f_N | Loss near cusp | Steps | Sweep s |",
             "|---|---|---|---|---|---|---|---|---|---|"]
    for t, b in e.items():
        a = b["all"]
        c = res["sweep_costs"][t]
        lines.append(f"| {t} | {_f(a['loss_db_median'])} ({_f(a['loss_db_p90'])}) | {_f(a['width_km_median'])}"
                     f" ({_f(a['width_km_p90'])}) | {_f(a['edge_err_km_median'])} | {_f(a['detected'], 3)}"
                     f" | {_f(a['cell_km_max'])} | {_f(b['near_fN']['loss_db_median'])}"
                     f" | {_f(b['near_cusp']['loss_db_median'])} | {c['n_steps']} | {_f(c['sweep_s'], 2)} |")
    lines += ["", "Loss by band (median dB):", "",
              "| Treatment | " + " | ".join(f"{lo:g}-{hi:g} MHz" for lo, hi in BANDS_MHZ) + " |",
              "|---|" + "---|" * len(BANDS_MHZ)]
    for t, b in e.items():
        lines.append(f"| {t} | " + " | ".join(_f(b.get(f"{lo:g}-{hi:g}", {}).get("loss_db_median", float("nan")))
                                              for lo, hi in BANDS_MHZ) + " |")
    if "chain" in res:
        lines += ["", "Chain (wp4_campaign, clutter off):", "",
                  "| Treatment | SCI-05 med | SCI-07 km med | INS-10 dB med (worst) | INS-05 worst km | O complete"
                  " | SCI-05 verdicts | run s |", "|---|---|---|---|---|---|---|---|"]
        for t, c in res["chain"].items():
            lines.append(f"| {t} | {_f(c['sci05_median'], 3)} | {_f(c['sci07_median_m'] / 1e3)} | "
                         f"{_f(c['ins10_median_db'])} ({_f(c['ins10_worst_db'])}) | {_f(c['ins05_worst_m'] / 1e3)} | "
                         f"{_f(c['completeness_o_median'], 2)} | {c['sci05_verdicts']} | {_f(c['run_s_median'], 1)} |")
    metrics = {f"loss_median_{t}": b["all"]["loss_db_median"] for t, b in e.items()}
    return "\n".join(lines), metrics


@register_study(
    "wp4_dispersion",
    description="WP4 Phase 3: dispersion and step-transition treatments compared at echo level and through the chain",
    work_packages=["WP4"],
    claims=["C1", "C2"],
    status="implemented",
    requirements=["INS-05", "INS-07", "INS-10", "SCI-05", "SCI-07", "SCI-09"],
)
def run(out_dir: Path, *, config: dict[str, Any], manifest: Any) -> StudyResult:
    res = compute(dict(config))
    rows = res["rows"]
    cols = sorted({k for r in rows for k in r}, key=lambda k: list(rows[0]).index(k) if k in rows[0] else 99)
    path = out_dir / "wp4_dispersion_echoes.csv"
    with path.open("w", encoding="utf-8") as fh:
        fh.write(",".join(cols) + "\n")
        for r in rows:
            fh.write(",".join(_f(r.get(c), 4) if isinstance(r.get(c), float) else str(r.get(c, "")) for c in cols)
                     + "\n")
    text, metrics = _summary(res)
    j = out_dir / "wp4_dispersion.json"
    j.write_text(json.dumps({k: v for k, v in res.items() if k != "rows"} | {"metrics": metrics}, indent=2,
                            default=str), encoding="utf-8")
    return StudyResult(outputs=[path, j], summary=text, metrics=metrics)
