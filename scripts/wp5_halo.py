#!/usr/bin/env python3
"""Route one with a diffuse halo around the facet (WP5 Phase 3 part 3, decision (e)).

    python scripts/wp5_halo.py OUT_DIR [--seeds 48] [--jobs N]

The straight scene of ``docs/drift-tier.md`` section 6 (a 2 km facet 425 km
below the platform, 5 MHz, drift 25 m/s, 209 sweeps) through L1, L2 and the
drift tier, with the facet's halo (:class:`tss.l1_propagation.facet.HaloSpec`)
in the cases agreed on 2026-11-27: 2 facet lengths at K = 10, 3 and 1, 1 and 4
lengths at K = 3, K = 0 (no facet, which should be declined), and no halo.
Each realisation draws a new halo.  White noise is added to the records at 20
and 10 dB per sweep of the facet's own echo (the same noise in every case), and
one level without noise.  Availability and drift error for ``pinned``, ``joint``
and ``joint_fading``.  L2 quantisation and thermal noise are off; no clutter.
About 4 minutes on 14 processes at 48 realisations.  Writes ``echoes.csv`` and
``summary.txt``.
"""

from __future__ import annotations

import argparse
import csv
import dataclasses
import math
import os
from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache
from pathlib import Path
from typing import Any

import numpy as np

CASES: dict[str, tuple[float, float] | None] = {
    "no halo": None, "r2 K10": (2.0, 10.0), "r2 K3": (2.0, 3.0), "r2 K1": (2.0, 1.0),
    "r1 K3": (1.0, 3.0), "r4 K3": (4.0, 3.0), "r2 K0": (2.0, 0.0),
}
#: Facet SNR per sweep, dB; None is no noise.
LEVELS: tuple[float | None, ...] = (None, 20.0, 10.0)
ESTIMATORS = ("pinned", "joint", "joint_fading")
F_HZ = 5.0e6
DRIFT_MPS = 25.0
N_SWEEPS = 209


def _sequence(halo: dict[str, Any] | None):
    from tss import sweep_law
    from tss.interfaces import PlatformState
    from tss.l1_propagation.plate_echo import ExponentialTopside
    from tss.l2_instrument.baseline import BaselineInstrument
    from tss.pipeline import Pipeline
    from tss.registry import build

    facet: dict[str, Any] = {"id": "f1", "depth_m": ExponentialTopside().depth_for_mu(F_HZ, 0.5),
                             "drift_mps": DRIFT_MPS}
    if halo is not None:
        facet["halo"] = halo
    tiers = {"l0": build("l0", "chapman"), "l1": build("l1", "facet", facets=[facet]),
             "l2": BaselineInstrument(quantise=False, thermal_noise=False)}
    t = sweep_law.drift_mode_times(N_SWEEPS, t0_s=-(N_SWEEPS // 2) / 45 - 1 / 90)
    plats = [PlatformState(alt_m=600e3, lat_deg=0.0, lon_deg=0.0) for _ in t]
    return Pipeline(tiers, seed=11).run_sequence(
        "halo", sweep=sweep_law.drift_mode(F_HZ), t_s=t, platforms=plats,
        alt_m=np.linspace(100e3, 600e3, 101), epoch="2026-01-01T12:00:00Z").sequence


def _prior():
    from tss.l3_signal.drift import OnboardPrior

    return OnboardPrior.exponential(f_n_sat_hz=1.5e6, scale_height_m=200e3, ref_alt_m=600e3, source="truth")


@lru_cache(maxsize=1)
def _reference() -> tuple[Any, float, float]:
    """The sequence without a halo, its beat, and the raw-sample noise variance at 0 dB.

    SNR per sweep is the budget's definition: the per-sweep complex value at the
    envelope's peak over its noise (the facet's own echo).
    """
    from tss.l3_signal.drift import DriftProcessing

    seq = _sequence(None)
    x = np.stack([r.samples[0] for r in seq.records])
    rec = seq.records[0]
    fb = DriftProcessing(detector="summed")._detect(x, rec.fs_hz, rec.sweep.k_f_hz_per_s)[0]
    tt = np.arange(x.shape[1]) / rec.fs_hz
    i0 = int(math.ceil(fb / rec.sweep.k_f_hz_per_s * rec.fs_hz)) + 3
    z = x[:, i0:] @ np.exp(-2j * math.pi * fb * (tt[i0:] - tt[i0:].mean()))
    return seq, float(fb), float(np.max(np.abs(z) ** 2)) / (x.shape[1] - i0)


def _noisy(seq, sig2: float, seed: tuple[int, ...]):
    rng = np.random.default_rng(seed)
    records = [dataclasses.replace(r, samples=r.samples + math.sqrt(sig2 / 2) * (
        rng.standard_normal(r.samples.shape) + 1j * rng.standard_normal(r.samples.shape))) for r in seq.records]
    return dataclasses.replace(seq, records=records)


def trial(job: tuple[str, int], estimators: tuple[str, ...] = ESTIMATORS) -> list[dict[str, Any]]:
    """``(case, seed)``: one halo draw, every noise level and estimator."""
    from tss.constants import C_LIGHT
    from tss.l3_signal.drift import DriftProcessing

    case, seed = job
    plain, fb, sig2_0 = _reference()
    k = plain.records[0].sweep.k_f_hz_per_s
    halo = CASES[case]
    seq = plain if halo is None else _sequence({"ratio": halo[0], "k": halo[1], "seed": seed})
    prior = _prior()
    rows = []
    for i, level in enumerate(LEVELS):
        data = seq if level is None else _noisy(seq, sig2_0 / 10 ** (level / 10), (seed, i))
        for est in estimators:
            prod = DriftProcessing(prior=prior, estimator=est).process_sequence(data)
            beats = 2.0 * k * np.asarray(prod.products.get("range_beat_m", []), dtype=float) / C_LIGHT
            # The echo within two range cells (about 66 Hz each) of the facet's beat.
            near = np.nonzero(np.abs(beats - fb) < 2 * 66.0)[0]
            err, formed = math.nan, False
            if near.size:
                j = int(near[np.argmin(np.abs(beats[near] - fb))])
                v = float(prod.products["drift_mps"][j])
                formed = math.isfinite(v)
                err = v - DRIFT_MPS if formed else math.nan
            rows.append({"case": case, "seed": seed, "level_db": level, "estimator": est,
                         "detected": bool(near.size), "formed": formed, "err_mps": err})
    return rows


def summarise(rows: list[dict[str, Any]]) -> dict[tuple[str, float | None, str], dict[str, Any]]:
    out: dict[tuple[str, float | None, str], dict[str, Any]] = {}
    for key in dict.fromkeys((r["case"], r["level_db"], r["estimator"]) for r in rows):
        sel = [r for r in rows if (r["case"], r["level_db"], r["estimator"]) == key]
        e = np.array([r["err_mps"] for r in sel if r["formed"]], dtype=float)
        out[key] = {"n": len(sel), "detected": sum(r["detected"] for r in sel), "formed": int(e.size),
                    "rms": float(np.sqrt(np.mean(e**2))) if e.size else math.nan,
                    "robust": float(1.4826 * np.median(np.abs(e - np.median(e)))) if e.size else math.nan,
                    "bias": float(e.mean()) if e.size else math.nan}
    return out


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("out", type=Path)
    ap.add_argument("--seeds", type=int, default=48)
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 6))
    a = ap.parse_args(argv)
    a.out.mkdir(parents=True, exist_ok=True)
    jobs = [(case, seed) for case in CASES for seed in range(a.seeds)]
    with ProcessPoolExecutor(a.jobs) as pool:
        rows = [r for got in pool.map(trial, jobs) for r in got]
    with open(a.out / "echoes.csv", "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    lines = [f"{'case':8s} {'level':>8s} {'estimator':13s}  detected  formed     rms  robust    bias  (m/s)"]
    for (case, level, est), v in summarise(rows).items():
        name = "no noise" if level is None else f"{level:g} dB"
        lines.append(f"{case:8s} {name:>8s} {est:13s}  {v['detected']:3d}/{v['n']:<3d}  {v['formed']:3d}/{v['n']:<3d} "
                     f"{v['rms']:7.2f} {v['robust']:7.2f} {v['bias']:+7.2f}")
    text = "\n".join(lines) + "\n"
    (a.out / "summary.txt").write_text(text, encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
