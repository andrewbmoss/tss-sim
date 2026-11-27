#!/usr/bin/env python3
"""Can a statistic of the echo spot a rough patch? (WP5 Phase 3 part 3, 2026-11-27).

    python scripts/wp5_rough_test.py OUT_DIR [--seeds 100] [--jobs N]

**Report only: nothing here is used by route one or the drift tier.**  Route
one's compact test forms a drift on a randomly phased patch about half the time
(``docs/drift-tier.md`` section 14).  This tries statistics on every echo the
test accepts, from the deterministic joint fit and the echo itself:

``residual_over_noise``
    the fit's residual variance over the measured noise variance;
``residual_over_echo``
    in the fitted lobe (model power above 0.1), residual power over fitted power;
``outside_over_noise``
    mean residual power outside the lobe (model power below 0.05, the fitted
    sidelobes taken out) over the noise;
``contrast``
    route one's own peak contrast (low is suspect).

Each statistic's threshold is set to reject a stated share (5 and 10 %) of a
calibration set: true reflectors without halo or fading, at the same SNR, from
noise seeds of their own.  Then the share rejected is counted, of the formed
echoes, for true reflectors (held-out seeds; also with Rayleigh specular fading
at 0.5 and 2 s through L2), for facets with a halo of 2 lengths at K = 10 and 3,
and for the patch alone (K = 0).  The scene is section 14's, through L1, L2 and
the tier, at 20 and 10 dB per sweep.  About 5 minutes on 14 processes at 100
seeds.  Writes ``echoes.csv`` and ``summary.txt``.
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

#: name: (halo ratio and K, or None; fading time in s, or None)
CASES: dict[str, tuple[tuple[float, float] | None, float | None]] = {
    "no halo": (None, None), "fading 0.5 s": (None, 0.5), "fading 2 s": (None, 2.0),
    "r2 K10": ((2.0, 10.0), None), "r2 K3": ((2.0, 3.0), None), "r2 K0": ((2.0, 0.0), None),
}
#: Statistic: True when a high value is suspect.
STATISTICS = {"residual_over_noise": True, "residual_over_echo": True, "outside_over_noise": True,
              "contrast": False}
LEVELS = (20.0, 10.0)
RATES = (0.05, 0.10)
#: Noise seeds of the calibration set, apart from the held-out ones (from 0).
CALIBRATION_SEED0 = 5000
F_HZ = 5.0e6
DRIFT_MPS = 25.0
N_SWEEPS = 209


def statistics(t: np.ndarray, field: np.ndarray, fit: Any, *, noise_var: float, **kw: float) -> dict[str, float]:
    """The fit's statistics for one echo (module docstring); ``contrast`` is route one's."""
    from tss.l3_signal.novel.route_one_joint import joint_model

    s = joint_model(t, t0_s=fit.t0_s, drift_mps=fit.drift_mps, width_s=fit.width_s, **kw)
    model = fit.amplitude * s
    r2 = np.abs(field - model) ** 2
    s2 = np.abs(s) ** 2
    lobe, outside = s2 > 0.1, s2 < 0.05
    return {
        "residual_over_noise": float(fit.noise_var / noise_var),
        "residual_over_echo": float(r2[lobe].sum() / max(np.sum(np.abs(model[lobe]) ** 2), 1e-300)),
        "outside_over_noise": float(np.mean(r2[outside]) / noise_var) if outside.any() else math.nan,
    }


def threshold(calibration: list[float], rate: float, high_is_suspect: bool) -> float:
    """The value that rejects ``rate`` of the calibration set."""
    return float(np.quantile(calibration, 1.0 - rate if high_is_suspect else rate))


def rejected(values: list[float], limit: float, high_is_suspect: bool) -> np.ndarray:
    v = np.asarray(values, dtype=float)
    return v > limit if high_is_suspect else v < limit


def _pipeline(halo: dict[str, Any] | None):
    from tss.l1_propagation.plate_echo import ExponentialTopside
    from tss.l2_instrument.baseline import BaselineInstrument
    from tss.pipeline import Pipeline
    from tss.registry import build

    facet: dict[str, Any] = {"id": "f1", "depth_m": ExponentialTopside().depth_for_mu(F_HZ, 0.5),
                             "drift_mps": DRIFT_MPS}
    if halo is not None:
        facet["halo"] = halo
    return Pipeline({"l0": build("l0", "chapman"), "l1": build("l1", "facet", facets=[facet]),
                     "l2": BaselineInstrument(quantise=False, thermal_noise=False)}, seed=11)


def _run(halo: dict[str, Any] | None):
    from tss import sweep_law
    from tss.interfaces import PlatformState

    t = sweep_law.drift_mode_times(N_SWEEPS, t0_s=-(N_SWEEPS // 2) / 45 - 1 / 90)
    plats = [PlatformState(alt_m=600e3, lat_deg=0.0, lon_deg=0.0) for _ in t]
    return _pipeline(halo).run_sequence(
        "rough", sweep=sweep_law.drift_mode(F_HZ), t_s=t, platforms=plats,
        alt_m=np.linspace(100e3, 600e3, 101), epoch="2026-01-01T12:00:00Z")


@lru_cache(maxsize=1)
def _plain():
    """The run without halo or fading, and the raw noise variance that is 0 dB per sweep."""
    from tss.l3_signal.drift import DriftProcessing

    run = _run(None)
    seq = run.sequence
    x = np.stack([r.samples[0] for r in seq.records])
    rec = seq.records[0]
    fb = DriftProcessing(detector="summed")._detect(x, rec.fs_hz, rec.sweep.k_f_hz_per_s)[0]
    tt = np.arange(x.shape[1]) / rec.fs_hz
    i0 = int(math.ceil(fb / rec.sweep.k_f_hz_per_s * rec.fs_hz)) + 3
    z = x[:, i0:] @ np.exp(-2j * math.pi * fb * (tt[i0:] - tt[i0:].mean()))
    return run, float(np.max(np.abs(z) ** 2)) / (x.shape[1] - i0)


def _sequence(case: str, seed: int):
    from tss.l2_instrument.baseline import BaselineInstrument
    from tss.l2_instrument.fading import FadingModel

    halo, tau = CASES[case]
    if halo is not None:
        return _run({"ratio": halo[0], "k": halo[1], "seed": seed}).sequence
    run = _plain()[0]
    if tau is None:
        return run.sequence
    l2 = BaselineInstrument(quantise=False, thermal_noise=False, fading=FadingModel(
        k_factor=math.inf, specular="rayleigh", tau_i_s=tau, tau_log_sigma=0.0))
    records = l2.receive_sequence(list(run.echoes), run.sequence.records[0].sweep, run.sequence.t_s,
                                  np.random.default_rng((seed, 77)))
    return dataclasses.replace(run.sequence, records=records)


def trial(job: tuple[str, int], levels: tuple[float, ...] = LEVELS) -> list[dict[str, Any]]:
    """``(case, seed)``: one sequence per level, the drift and the statistics of its echo."""
    from tss.l3_signal import drift as drift_module
    from tss.l3_signal.drift import DriftProcessing, OnboardPrior
    from tss.l3_signal.novel.route_one_joint import sample_noise_variance

    case, seed = job
    clean = _sequence(case, seed)
    sig2_0 = _plain()[1]
    prior = OnboardPrior.exponential(f_n_sat_hz=1.5e6, scale_height_m=200e3, ref_alt_m=600e3, source="truth")
    real = drift_module.estimate_route_one_joint
    rows = []
    for i, level in enumerate(levels):
        rng = np.random.default_rng((seed, i))
        sig2 = sig2_0 / 10 ** (level / 10)
        records = [dataclasses.replace(r, samples=r.samples + math.sqrt(sig2 / 2) * (
            rng.standard_normal(r.samples.shape) + 1j * rng.standard_normal(r.samples.shape)))
            for r in clean.records]
        seq = dataclasses.replace(clean, records=records)
        seen: list[tuple[Any, ...]] = []

        def capture(t, field, seen=seen, **kw):
            fit = real(t, field, **kw)
            seen.append((np.array(t), np.array(field), fit, kw))
            return fit

        drift_module.estimate_route_one_joint = capture
        try:
            # The strongest echo only.
            prod = DriftProcessing(prior=prior, estimator="joint", max_echoes=1).process_sequence(seq)
        finally:
            drift_module.estimate_route_one_joint = real
        row: dict[str, Any] = {"case": case, "seed": seed, "level_db": level, "detected": bool(prod.echo_id),
                               "formed": False, "err_pinned_mps": math.nan, "err_joint_mps": math.nan,
                               **{name: math.nan for name in STATISTICS}}
        if seen and seen[-1][2].ok:
            t, field, fit, kw = seen[-1]
            p = {n: v[0] for n, v in prod.products.items()}
            det = prod.flags["detection"]
            n = records[0].samples.shape[1]
            noise = ((1.0 - float(p["blanked_fraction"])) * n
                     * sample_noise_variance(det["noise_power_per_bin"], det["segment_samples"])
                     / 10.0 ** (float(p["power_peak_db"]) / 10.0))
            model_kw = {name: kw[name] for name in ("f_hz", "group_path_m", "mu_sat", "v_s_mps")}
            row.update(statistics(t, field, fit, noise_var=noise, **model_kw))
            row.update({"formed": True, "contrast": float(p["peak_contrast"]),
                        "err_pinned_mps": float(p["drift_pinned_mps"]) - DRIFT_MPS,
                        "err_joint_mps": float(p["drift_mps"]) - DRIFT_MPS})
        rows.append(row)
    return rows


def _report(rows: list[dict[str, Any]]) -> str:
    lines = []
    for level in LEVELS:
        at = [r for r in rows if r["level_db"] == level and r["formed"]]
        calib = [r for r in at if r["case"] == "no halo" and r["seed"] >= CALIBRATION_SEED0]
        lines.append(f"== {level:g} dB: calibration {len(calib)} formed true reflectors")
        for name, high in STATISTICS.items():
            for rate in RATES:
                limit = threshold([r[name] for r in calib], rate, high)
                lines.append(f"   {name} {'>' if high else '<'} {limit:.4g} (set for {rate:.0%})")
                for case in CASES:
                    sel = [r for r in at if r["case"] == case and r["seed"] < CALIBRATION_SEED0]
                    n_all = sum(r["case"] == case and r["level_db"] == level and r["seed"] < CALIBRATION_SEED0
                                for r in rows)
                    out = rejected([r[name] for r in sel], limit, high)
                    kept = np.array([r["err_joint_mps"] for r, o in zip(sel, out, strict=True) if not o])
                    before = np.array([r["err_joint_mps"] for r in sel])
                    rms = lambda e: math.sqrt(np.mean(e**2)) if e.size else math.nan  # noqa: E731
                    lines.append(f"      {case:13s} formed {len(sel):3d}/{n_all:<3d} rejected {out.sum():3d} "
                                 f"({out.mean() if out.size else math.nan:6.1%})  joint rms {rms(before):7.2f} -> "
                                 f"{rms(kept):7.2f} of {kept.size}")
    return "\n".join(lines) + "\n"


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("out", type=Path)
    ap.add_argument("--seeds", type=int, default=100)
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 6))
    a = ap.parse_args(argv)
    a.out.mkdir(parents=True, exist_ok=True)
    jobs = [(case, seed) for case in CASES for seed in range(a.seeds)]
    jobs += [("no halo", CALIBRATION_SEED0 + i) for i in range(2 * a.seeds)]
    with ProcessPoolExecutor(a.jobs) as pool:
        rows = [r for got in pool.map(trial, jobs) for r in got]
    with open(a.out / "echoes.csv", "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    text = _report(rows)
    (a.out / "summary.txt").write_text(text, encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
