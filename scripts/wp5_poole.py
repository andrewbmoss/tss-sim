#!/usr/bin/env python3
"""Poole's velocity as the comparator on the drift sequences (WP5 Phase 3 part 3).

    python scripts/wp5_poole.py OUT_DIR [--seeds 48] [--jobs N]

Technical Volume section 8.21, item 10, velocity half only (Andrew,
2026-11-27, option A).  Poole's (1985) phase velocity is formed between
successive sweeps of the same drift sequences the tier processes
(:func:`tss.l3_signal.poole_range.velocity_series`).  **The range half cannot
be formed from this data**: it needs a cell at a second start frequency, and a
drift sequence has one.

Poole's phase velocity ``V*`` is the rate of the phase range; ``-V* / mu_sat``
is the line-of-sight rate in route one's sense (positive closing), which holds
the platform's motion, the drift and any tilt together.  For each sequence:
the RMS difference of Poole's pair rates from the line-of-sight rate L1 put in
(over the main lobe, pairs above 0.1 of the peak power), and from the tier's
own phase steps; and, as drift estimates, Poole's rate at the strongest pair and
its power-weighted mean over the lobe, beside ``pinned``, ``joint`` and
``joint_fading``.

Scenes (``docs/drift-tier.md`` sections 6 and 14): the 2 km facet overhead, the
same 2 km along track, the 80 km tapered plate, and the facet with a halo of 2
lengths at K = 3.  48 noise draws (halo: draws) at 20 and 10 dB per sweep, and
none.  About 3 minutes on 14 processes.  Writes ``sequences.csv`` and
``summary.txt``.
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import math
import os
from concurrent.futures import ProcessPoolExecutor
from functools import cache, lru_cache
from pathlib import Path
from typing import Any

import numpy as np

SCENES: dict[str, dict[str, Any]] = {
    "nadir": {},
    "2 km along": {"x_m": 2.0e3},
    "extended": {"length_m": 80e3, "taper": "tukey"},
    "halo r2 K3": {"halo": {"ratio": 2.0, "k": 3.0}},
}
LEVELS: tuple[float | None, ...] = (None, 20.0, 10.0)
ESTIMATORS = ("pinned", "joint", "joint_fading")
F_HZ = 5.0e6
DRIFT_MPS = 25.0
N_SWEEPS = 209
#: Two range cells of beat frequency, Hz.
_NEAR_HZ = 2 * 66.0


def _run(scene: str, seed: int):
    from tss import sweep_law
    from tss.interfaces import PlatformState
    from tss.l1_propagation.plate_echo import ExponentialTopside
    from tss.l2_instrument.baseline import BaselineInstrument
    from tss.pipeline import Pipeline
    from tss.registry import build

    facet: dict[str, Any] = {"id": "f1", "depth_m": ExponentialTopside().depth_for_mu(F_HZ, 0.5),
                             "drift_mps": DRIFT_MPS, **SCENES[scene]}
    if "halo" in facet:
        facet["halo"] = {**facet["halo"], "seed": seed}
    tiers = {"l0": build("l0", "chapman"), "l1": build("l1", "facet", facets=[facet]),
             "l2": BaselineInstrument(quantise=False, thermal_noise=False)}
    t = sweep_law.drift_mode_times(N_SWEEPS, t0_s=-(N_SWEEPS // 2) / 45 - 1 / 90)
    plats = [PlatformState(alt_m=600e3, lat_deg=0.0, lon_deg=0.0) for _ in t]
    return Pipeline(tiers, seed=11).run_sequence(
        scene, sweep=sweep_law.drift_mode(F_HZ), t_s=t, platforms=plats,
        alt_m=np.linspace(100e3, 600e3, 101), epoch="2026-01-01T12:00:00Z")


@cache
def _fixed_run(scene: str):
    return _run(scene, 0)


@lru_cache(maxsize=1)
def _noise_reference() -> float:
    """The raw-sample noise variance that is 0 dB per sweep for the overhead facet."""
    from tss.l3_signal.drift import DriftProcessing

    seq = _fixed_run("nadir").sequence
    x = np.stack([r.samples[0] for r in seq.records])
    rec = seq.records[0]
    fb = DriftProcessing(detector="summed")._detect(x, rec.fs_hz, rec.sweep.k_f_hz_per_s)[0]
    tt = np.arange(x.shape[1]) / rec.fs_hz
    i0 = int(math.ceil(fb / rec.sweep.k_f_hz_per_s * rec.fs_hz)) + 3
    z = x[:, i0:] @ np.exp(-2j * math.pi * fb * (tt[i0:] - tt[i0:].mean()))
    return float(np.max(np.abs(z) ** 2)) / (x.shape[1] - i0)


def _noisy(seq, sig2: float, seed: tuple[int, ...]):
    rng = np.random.default_rng(seed)
    records = [dataclasses.replace(r, samples=r.samples + math.sqrt(sig2 / 2) * (
        rng.standard_normal(r.samples.shape) + 1j * rng.standard_normal(r.samples.shape))) for r in seq.records]
    return dataclasses.replace(seq, records=records)


def _truth_rate(run, f_hz: float, mu_sat: float) -> np.ndarray:
    """The line-of-sight rate L1 put in, between successive sweeps, from the facet row's Doppler."""
    from tss.constants import C_LIGHT

    fd = np.array([float(e.doppler_hz[[str(f).startswith("facet:") for f in e.flags]][0]) for e in run.echoes])
    return C_LIGHT * 0.5 * (fd[1:] + fd[:-1]) / (2.0 * f_hz * mu_sat)


def _tier_rate(seq, fb: float, f_hz: float, mu_sat: float) -> tuple[np.ndarray, np.ndarray]:
    """The tier's own phase steps (its gate and demodulation), as a rate and a pair power."""
    from tss.constants import C_LIGHT

    x = np.stack([r.samples[0] for r in seq.records])
    rec = seq.records[0]
    tt = np.arange(x.shape[1]) / rec.fs_hz
    i0 = int(math.ceil(fb / rec.sweep.k_f_hz_per_s * rec.fs_hz)) + 3
    z = x[:, i0:] @ np.exp(-2j * math.pi * fb * (tt[i0:] - tt[i0:].mean()))
    step = np.conj(z[1:]) * z[:-1]  # the tier's field is conj(z)
    period = float(np.mean(np.diff(seq.t_s)))
    rate = C_LIGHT * np.angle(step) / (2.0 * np.pi * period) / (2.0 * f_hz * mu_sat)
    return rate, np.abs(step)


def trial(job: tuple[str, int], estimators: tuple[str, ...] = ESTIMATORS,
          levels: tuple[float | None, ...] = LEVELS) -> list[dict[str, Any]]:
    """``(scene, seed)``: one sequence per noise level, Poole beside the tier."""
    from tss.constants import C_LIGHT
    from tss.l3_signal.drift import DriftProcessing, OnboardPrior
    from tss.l3_signal.poole_range import velocity_series

    scene, seed = job
    run = _run(scene, seed) if "halo" in SCENES[scene] else _fixed_run(scene)
    sig2_0 = _noise_reference()
    prior = OnboardPrior.exponential(f_n_sat_hz=1.5e6, scale_height_m=200e3, ref_alt_m=600e3, source="truth")
    rows = []
    for i, level in enumerate(levels):
        seq = run.sequence if level is None else _noisy(run.sequence, sig2_0 / 10 ** (level / 10), (seed, i))
        rec = seq.records[0]
        sweep = rec.sweep
        k = sweep.k_f_hz_per_s
        row: dict[str, Any] = {"scene": scene, "seed": seed, "level_db": level, "detected": False,
                               "poole_pair_rms_mps": math.nan, "tier_pair_rms_mps": math.nan,
                               "poole_minus_tier_rms_mps": math.nan, "err_poole_peak": math.nan,
                               "err_poole_mean": math.nan, **{f"err_{e}": math.nan for e in estimators}}
        echo = None
        for est in estimators:
            prod = DriftProcessing(prior=prior, estimator=est).process_sequence(seq)
            if not prod.echo_id:
                continue
            j = int(np.argmax(prod.flags["detection"]["excess_db"]))
            row[f"err_{est}"] = float(prod.products["drift_mps"][j]) - DRIFT_MPS
            echo = {n: prod.products[n][j] for n in ("range_beat_m", "f_eff_hz", "mu_sat")}
        if echo is None:
            rows.append(row)
            continue
        row["detected"] = True
        fb = 2.0 * k * float(echo["range_beat_m"]) / C_LIGHT
        f_eff, mu = float(echo["f_eff_hz"]), float(echo["mu_sat"])
        cells = [r.samples[0] for r in seq.records]
        # Poole picks a spectral component: the strongest bin near the echo's beat.
        freqs = np.arange(cells[0].size // 2 + 1) * rec.fs_hz / cells[0].size
        near = np.nonzero(np.abs(freqs - fb) < _NEAR_HZ)[0]
        tries = [velocity_series(cells, rec.fs_hz, k_hz_per_s=k, f0_hz=sweep.f_start_hz,
                                 period_s=sweep.sweep_period_s, bin_index=int(b)) for b in near]
        got = max(tries, key=lambda q: float(q["power"].sum()))
        rate = -got["velocity_mps"] / mu
        power = got["power"]
        lobe = power >= 0.1 * power.max()
        truth = _truth_rate(run, f_eff, mu)
        tier, _ = _tier_rate(seq, fb, f_eff, mu)
        row.update({
            "poole_pair_rms_mps": float(np.sqrt(np.mean((rate - truth)[lobe] ** 2))),
            "tier_pair_rms_mps": float(np.sqrt(np.mean((tier - truth)[lobe] ** 2))),
            "poole_minus_tier_rms_mps": float(np.sqrt(np.mean((rate - tier)[lobe] ** 2))),
            "err_poole_peak": float(rate[int(np.argmax(power))]) - DRIFT_MPS,
            "err_poole_mean": float(np.sum(power[lobe] * rate[lobe]) / power[lobe].sum()) - DRIFT_MPS,
        })
        rows.append(row)
    return rows


def _stats(values: list[float]) -> dict[str, float]:
    e = np.array([v for v in values if math.isfinite(v)], dtype=float)
    if not e.size:
        return {"n": 0, "rms": math.nan, "robust": math.nan, "bias": math.nan}
    return {"n": int(e.size), "rms": float(np.sqrt(np.mean(e**2))),
            "robust": float(1.4826 * np.median(np.abs(e - np.median(e)))), "bias": float(e.mean())}


def summarise(rows: list[dict[str, Any]]) -> dict[tuple[str, float | None], dict[str, Any]]:
    out: dict[tuple[str, float | None], dict[str, Any]] = {}
    for key in dict.fromkeys((r["scene"], r["level_db"]) for r in rows):
        sel = [r for r in rows if (r["scene"], r["level_db"]) == key]
        case: dict[str, Any] = {"n": len(sel), "detected": sum(r["detected"] for r in sel)}
        for name in sel[0]:
            if name.startswith("err_"):
                case[name] = _stats([r[name] for r in sel])
            elif name.endswith("_rms_mps"):
                v = [r[name] for r in sel if math.isfinite(r[name])]
                case[name] = float(np.median(v)) if v else math.nan
        out[key] = case
    return out


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("out", type=Path)
    ap.add_argument("--seeds", type=int, default=48)
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 6))
    a = ap.parse_args(argv)
    a.out.mkdir(parents=True, exist_ok=True)
    # The noise-free level is one sequence per fixed scene; the halo draws anew each seed.
    jobs = [(scene, seed) for scene in SCENES for seed in range(a.seeds)]
    with ProcessPoolExecutor(a.jobs) as pool:
        rows = [r for got in pool.map(trial, jobs) for r in got]
    with open(a.out / "sequences.csv", "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    lines = []
    for (scene, level), c in summarise(rows).items():
        name = "no noise" if level is None else f"{level:g} dB"
        lines.append(f"== {scene}, {name}: detected {c['detected']}/{c['n']}; pair rate against the injected "
                     f"line of sight, median RMS: Poole {c['poole_pair_rms_mps']:.2f}, tier "
                     f"{c['tier_pair_rms_mps']:.2f}; Poole against the tier {c['poole_minus_tier_rms_mps']:.2f} m/s")
        for key in ("err_poole_peak", "err_poole_mean", *(f"err_{e}" for e in ESTIMATORS)):
            v = c[key]
            lines.append(f"   drift error {key[4:]:13s} formed {v['n']:3d}/{c['n']:<3d} rms {v['rms']:8.2f} "
                         f"robust {v['robust']:8.2f} bias {v['bias']:+8.2f}")
    text = "\n".join(lines) + "\n"
    (a.out / "summary.txt").write_text(text, encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
