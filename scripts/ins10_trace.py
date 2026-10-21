"""WP3 Phase 0: where the INS-10 minimum at the truth echoes comes from.

The baseline campaign scored INS-10 at 7.0 dB at 800 km on PHaRLAP's IRI-2020
and 18.7 dB on PyIRI (RESUME, 2026-10-18).  Two commands:

``echoes MODEL OUT``
    Runs the baseline configurations through the pipeline exactly as
    ``tss run baseline`` does (same seed, same generator per configuration)
    and writes every truth echo with its L1 path gain, absorption, reflection
    height and the SNR L3 measures there, plus the profile scalars and the
    distance of the highest echo of each mode below that mode's cutoff.

``cusp MODEL OUT [CONFIGURATION]``
    Scales the truth density by 0.97 to 1.00 (foF2 by 1.5 %, one 150 kHz
    step at 10.6 MHz) with one fixed noise draw and records the INS-10
    minimum each time.  If the minimum were a property of the topside it would
    not move; it moves with where the last step lands below the cutoff.

``MODEL`` is ``pharlap`` (needs MATLAB and PHaRLAP) or ``pyiri``.
"""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np

from tss.config import load_campaign, load_yaml, resolve_path
from tss.l5_assessment.scoring import _NO_ECHO_FLAGS, truth_echo_cells, truth_echo_metrics
from tss.pipeline import Pipeline, build_tiers
from tss.studies import _grids, _platform_for, _sweep_from_config, build_instrument

GYRO_HZ_PER_T = 2.799249e10


def _setup(model: str):
    loaded = load_campaign("baseline", [f"tiers.l0.model={model}"])
    config = loaded.config
    shared = build_tiers({k: v for k, v in loaded.tiers.items() if k != "l2"})
    sweep, _ = _sweep_from_config(config, sys.stdout)
    alt_m, freq_hz = _grids(config)
    epoch = str(config.get("epoch", "2026-01-01T12:00:00Z"))
    return loaded, config, shared, sweep, alt_m, freq_hz, epoch


def _cutoffs(truth) -> dict[str, float]:
    """foF2, fxF2 (from |B| at hmF2) and hmF2 of a truth profile."""
    ne = np.asarray(truth.ne, dtype=float)
    i = int(np.nanargmax(ne))
    fo = 8.98 * float(np.sqrt(ne[i]))
    fh = float("nan")
    if truth.b_field_t is not None:
        b = np.asarray(truth.b_field_t, dtype=float)
        fh = GYRO_HZ_PER_T * float(np.linalg.norm(b[:, i]))
    fx = fh / 2.0 + float(np.sqrt(fh**2 / 4.0 + fo**2)) if np.isfinite(fh) else float("nan")
    return {"fof2_hz": fo, "fxf2_hz": fx, "fh_at_hmf2_hz": fh, "hmf2_m": float(truth.alt_m[i])}


def _last_below_cutoff(cells: dict, snr: np.ndarray, cut: dict[str, float]) -> dict[str, dict]:
    out: dict[str, dict] = {}
    mode = np.array(cells["mode"])
    f = np.asarray(cells["freq_hz"], dtype=float)
    for lab, key in (("O", "fof2_hz"), ("X", "fxf2_hz")):
        sel = np.flatnonzero(mode == lab)
        if sel.size:
            i = sel[np.argmax(f[sel])]
            out[lab] = {"f_hz": float(f[i]), "below_cutoff_hz": float(cut[key] - f[i]),
                        "range_m": float(cells["range_m"][i]), "snr_db": float(snr[i])}
    return out


def cmd_echoes(model: str, out: Path) -> None:
    loaded, config, shared, sweep, alt_m, freq_hz, epoch = _setup(model)
    out.mkdir(parents=True, exist_ok=True)
    for name in loaded.configurations:
        configuration = load_yaml(resolve_path(name, subdir="configurations"))
        pipe = Pipeline({**shared, "l2": build_instrument(loaded.tiers["l2"], config, configuration)},
                        seed=loaded.seed, keep_intermediates=True)
        platform = _platform_for(configuration, epoch)
        res = pipe.run_configuration(name, platform=platform, sweep=sweep, freq_hz=freq_hz,
                                     alt_m=alt_m, epoch=epoch)
        e = res.echoes
        cells = truth_echo_cells(e)
        freq = np.asarray(e.freq_hz, float)
        delay = np.asarray(e.group_delay_s, float)
        flags = ([str(x) for x in np.asarray(e.flags, dtype=object).ravel()[: freq.size]]
                 if e.flags is not None else [""] * freq.size)
        ok = (np.isfinite(freq) & np.isfinite(delay) & (delay > 0)
              & np.array([fl not in _NO_ECHO_FLAGS for fl in flags]))
        snr = np.asarray(res.truth_echo_snr["snr_db"], float)
        amp = np.asarray(e.amplitude_db, float)[ok]
        absn = np.asarray(e.absorption_db, float)[ok]
        refl = (np.asarray(e.reflection_alt_m, float)[ok] if e.reflection_alt_m is not None
                else np.full(snr.size, np.nan))
        with (out / f"{model}_{name}.csv").open("w", encoding="utf-8") as fh:
            fh.write("freq_hz,mode,range_m,amplitude_db,absorption_db,reflection_alt_m,snr_db\n")
            for k in range(snr.size):
                fh.write(f"{cells['freq_hz'][k]:.1f},{cells['mode'][k]},{cells['range_m'][k]:.1f},"
                         f"{amp[k]:.3f},{absn[k]:.3f},{refl[k]:.1f},{snr[k]:.3f}\n")
        t = res.truth
        ne_sat = float(np.interp(platform.alt_m, np.asarray(t.alt_m, float), np.asarray(t.ne, float)))
        cut = _cutoffs(t)
        record = {
            "model": model, "configuration": name, "alt_sat_m": float(platform.alt_m),
            "lat_deg": platform.lat_deg, "lon_deg": platform.lon_deg, "epoch": epoch,
            "fn_sat_hz": 8.98 * float(np.sqrt(ne_sat)), **cut,
            "last_echo_below_cutoff": _last_below_cutoff(cells, snr, cut),
            "metrics": truth_echo_metrics({**cells, "snr_db": snr}),
        }
        (out / f"{model}_{name}.json").write_text(json.dumps(record, indent=2), encoding="utf-8")
        m = record["metrics"]
        print(f"{model} {name}: worst {m['snr_db_min_at_truth_echoes']:.1f} dB at "
              f"{m['snr_db_min_at_truth_echoes_freq_hz'] / 1e6:.3f} MHz, median "
              f"{m['snr_db_median_at_truth_echoes']:.1f} dB; last X "
              f"{record['last_echo_below_cutoff'].get('X', {}).get('below_cutoff_hz', float('nan')) / 1e3:.0f}"
              " kHz below fxF2")


def cmd_cusp(model: str, out: Path, cfgname: str) -> None:
    loaded, config, shared, sweep, alt_m, freq_hz, epoch = _setup(model)
    out.mkdir(parents=True, exist_ok=True)
    configuration = load_yaml(resolve_path(cfgname, subdir="configurations"))
    l2 = build_instrument(loaded.tiers["l2"], config, configuration)
    platform = _platform_for(configuration, epoch)
    truth0 = shared["l0"].state_at(platform.lat_deg, platform.lon_deg, epoch, alt_m)
    rows = []
    for s in np.linspace(0.97, 1.0, 13):
        truth = replace(truth0, ne=np.asarray(truth0.ne, float) * s)
        rng = np.random.default_rng(20260921)
        echoes = shared["l1"].echoes(truth, platform, freq_hz, rng)
        signal = l2.receive(echoes, sweep, rng)
        cells = truth_echo_cells(echoes)
        snr = np.asarray(shared["l3"].snr_at(signal, cells["freq_hz"], cells["range_m"]))
        m = truth_echo_metrics({**cells, "snr_db": snr})
        cut = _cutoffs(truth)
        f = np.asarray(cells["freq_hz"], float)
        mode = np.array(cells["mode"])
        # Leave out each mode's last step below its own cutoff.
        last = np.zeros(f.size, dtype=bool)
        for lab in ("O", "X"):
            sel = np.flatnonzero(mode == lab)
            if sel.size:
                last[sel[np.argmax(f[sel])]] = True
        rows.append({
            "scale": float(s), **cut,
            "worst_db": m["snr_db_min_at_truth_echoes"],
            "worst_f_hz": m["snr_db_min_at_truth_echoes_freq_hz"],
            "worst_o_db": m["snr_db_min_at_truth_echoes_o"],
            "worst_x_db": m["snr_db_min_at_truth_echoes_x"],
            "median_db": m["snr_db_median_at_truth_echoes"],
            "worst_without_last_step_db": float(np.nanmin(snr[~last])),
            "last_echo_below_cutoff": _last_below_cutoff(cells, snr, cut),
        })
        r = rows[-1]
        print(f"s={s:.4f} foF2={r['fof2_hz'] / 1e6:.4f} worst={r['worst_db']:.1f} at "
              f"{r['worst_f_hz'] / 1e6:.3f}  without last step {r['worst_without_last_step_db']:.1f}  "
              f"median {r['median_db']:.1f}")
    (out / f"cusp_{model}_{cfgname}.json").write_text(json.dumps(rows, indent=2), encoding="utf-8")


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = p.add_subparsers(dest="cmd", required=True)
    a = sub.add_parser("echoes")
    a.add_argument("model", choices=["pharlap", "pyiri"])
    a.add_argument("out", type=Path)
    b = sub.add_parser("cusp")
    b.add_argument("model", choices=["pharlap", "pyiri"])
    b.add_argument("out", type=Path)
    b.add_argument("configuration", nargs="?", default="leo_800km_single")
    args = p.parse_args(argv)
    if args.cmd == "echoes":
        cmd_echoes(args.model, args.out)
    else:
        cmd_cusp(args.model, args.out, args.configuration)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
