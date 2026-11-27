#!/usr/bin/env python3
"""The power cut of route one's Doppler history (WP5 Phase 3 part 3, decision (n)).

    python scripts/wp5_power_cut.py OUT_DIR [--seeds 200] [--seed0 40000] [--jobs N]

Route one forms its Doppler history from sweep pairs whose power is at least
``main_lobe`` of the peak (0.1 in use).  The history feeds the free-slope test,
which decides whether a drift is formed at all, and the pinned zero crossing.
This compares 0.1 with two higher cuts, 0.25 and 0.5: availability (the
fraction formed) and the drift spread, for the pinned estimator and for
``joint_fading`` behind the same test.  **An experiment: nothing is adopted.**

Parameter level, the straight scene of ``docs/drift-tier.md`` section 9, SNR 15
and 10 dB per sweep, without fading and with Rayleigh fading at 0.5 and 2 s.
``joint_fading`` is fitted once per echo, started from the pinned estimate of
the lowest cut that formed one; its search covers the grid whatever the start.
Also runs ``c3_verify``'s noise-free scenes at each cut, which show what a cut
does to the compact-facet test itself.  About 4 minutes on 14 processes.
Writes ``echoes.csv`` and ``summary.txt``.
"""

from __future__ import annotations

import argparse
import csv
import math
import os
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from typing import Any

import numpy as np

CUTS = (0.1, 0.25, 0.5)
SNRS_DB = (15.0, 10.0)
TAUS_S = (None, 0.5, 2.0)
F_HZ, P_M, MU, WIDTH_S, V_S = 4.988e6, 529.8e3, 0.9537, 0.926, 7600.0
DRIFT_MPS = 25.0
N_SWEEPS = 209
KW = {"f_hz": F_HZ, "group_path_m": P_M, "mu_sat": MU, "v_s_mps": V_S}


def _times() -> np.ndarray:
    return (np.arange(N_SWEEPS) - N_SWEEPS // 2) / 45.0


def _noise(rng: np.random.Generator, n: int) -> np.ndarray:
    return (rng.standard_normal(n) + 1j * rng.standard_normal(n)) / math.sqrt(2)


def _echo(snr_db: float, tau_c_s: float | None, seed: int) -> np.ndarray:
    """One echo of the scene: unit noise, the facet at ``snr_db`` per sweep."""
    from tss.l3_signal.novel.route_one_joint import joint_model

    t = _times()
    rng = np.random.default_rng(seed)
    g: Any = 1.0
    if tau_c_s is not None:
        lag = np.abs(t[:, None] - t[None, :])
        g = np.linalg.cholesky(np.exp(-lag / tau_c_s) + 1e-12 * np.eye(t.size)) @ _noise(rng, t.size)
    s = joint_model(t, t0_s=0.0, drift_mps=DRIFT_MPS, width_s=WIDTH_S, **KW)
    return 10 ** (snr_db / 20) * s * g + _noise(rng, t.size)


def trial(job: tuple[float, float | None, int], fit_fading: bool = True) -> dict[str, Any]:
    """``(snr_db, tau_c_s, seed)``: route one at every cut, and the fading fit."""
    from tss.l3_signal.novel.route_one_estimator import estimate_route_one
    from tss.l3_signal.novel.route_one_joint import estimate_route_one_joint_fading

    snr, tau, seed = job
    t = _times()
    y = _echo(snr, tau, seed)
    row: dict[str, Any] = {"snr_db": snr, "tau_c_s": tau, "seed": seed}
    first = None
    for c in CUTS:
        e = estimate_route_one(t, y, main_lobe=c, **KW)
        row[f"formed_{c:g}"] = bool(e.separable)
        row[f"pinned_err_{c:g}"] = e.drift_mps - DRIFT_MPS if e.separable else math.nan
        row[f"slope_ratio_{c:g}"] = e.slope_ratio
        if e.separable and first is None:
            first = e
    row["fading_err"] = math.nan
    if fit_fading and first is not None:
        fit = estimate_route_one_joint_fading(t, y, noise_var=1.0, t0_start_s=first.t_power_peak_s,
                                              drift_start_mps=first.drift_mps, **KW)
        if fit.ok:
            row["fading_err"] = fit.drift_mps - DRIFT_MPS
    return row


def _spread(err: np.ndarray) -> dict[str, float]:
    e = err[np.isfinite(err)]
    if e.size < 2:
        return {"n": int(e.size), "rms": math.nan, "robust": math.nan, "bias": math.nan}
    return {"n": int(e.size), "rms": float(np.sqrt(np.mean(e**2))),
            "robust": float(1.4826 * np.median(np.abs(e - np.median(e)))), "bias": float(e.mean())}


def summarise(rows: list[dict[str, Any]]) -> dict[tuple[float, float | None], dict[str, Any]]:
    """Availability and spread per case and cut; ``common`` is formed at every cut."""
    out: dict[tuple[float, float | None], dict[str, Any]] = {}
    for key in dict.fromkeys((r["snr_db"], r["tau_c_s"]) for r in rows):
        sel = [r for r in rows if (r["snr_db"], r["tau_c_s"]) == key]
        formed = np.array([[r[f"formed_{c:g}"] for c in CUTS] for r in sel], dtype=bool)
        common = formed.all(axis=1)
        fad = np.array([r["fading_err"] for r in sel], dtype=float)
        cuts = {}
        for i, c in enumerate(CUTS):
            err = np.array([r[f"pinned_err_{c:g}"] for r in sel], dtype=float)
            s = formed[:, i]
            cuts[c] = {"formed": int(s.sum()), "gained": int((s & ~formed[:, 0]).sum()),
                       "lost": int((~s & formed[:, 0]).sum()),
                       "pinned": _spread(err[s]), "pinned_common": _spread(err[common]),
                       "fading": _spread(fad[s]), "fading_common": _spread(fad[common])}
        out[key] = {"n": len(sel), "common": int(common.sum()), "cuts": cuts}
    return out


def noise_free(cut: float) -> str:
    """``c3_verify``'s scenes with route one at ``cut``."""
    from tss.l3_signal.novel import route_one_estimator as r1
    from tss.studies import c3_verify as cv

    def at_cut(t, field, **kw):
        return r1.estimate_route_one(t, field, main_lobe=cut, **kw)

    keep, cv.estimate_route_one = cv.estimate_route_one, at_cut
    try:
        rows = cv._all_trials(0)
    finally:
        cv.estimate_route_one = keep
    a = [r for r in rows if r.case.startswith("A")]
    b = [r for r in rows if r.case.startswith("B")]
    c = [r for r in rows if r.case.startswith("C")]
    d = [r for r in rows if r.case.startswith("D")]
    return (f"cut {cut:g}: A formed {sum(r.separable for r in a)}/{len(a)}, max |err| "
            f"{max(abs(r.drift_err_mps) for r in a if r.separable):.3f} m/s, slope ratio "
            f"{min(r.slope_ratio for r in a):.3f}-{max(r.slope_ratio for r in a):.3f}; B formed "
            f"{sum(r.separable for r in b)}/{len(b)}; C declined {sum(not r.separable for r in c)}/{len(c)}; "
            f"D formed at plate lengths (km) {[r.plate_km for r in d if r.separable]}")


def _report(summary: dict, checks: list[str]) -> str:
    def s(x: dict[str, float]) -> str:
        return f"rms {x['rms']:6.2f} robust {x['robust']:6.2f} bias {x['bias']:+6.2f} (n {x['n']})"

    lines = []
    for (snr, tau), case in summary.items():
        lines.append(f"== {snr:g} dB, tau_c {tau}: {case['n']} echoes, formed at all cuts {case['common']}")
        for c, v in case["cuts"].items():
            lines.append(f"   cut {c:g}: formed {v['formed'] / case['n']:6.1%} (+{v['gained']} -{v['lost']} "
                         f"against {CUTS[0]:g})")
            lines.append(f"      pinned        formed {s(v['pinned'])} | common {s(v['pinned_common'])}")
            lines.append(f"      joint_fading  formed {s(v['fading'])} | common {s(v['fading_common'])}")
    return "\n".join([*lines, "noise-free c3_verify", *(f"   {c}" for c in checks)]) + "\n"


def main(argv: list[str] | None = None) -> None:
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("out", type=Path)
    ap.add_argument("--seeds", type=int, default=200)
    ap.add_argument("--seed0", type=int, default=40000)
    ap.add_argument("--jobs", type=int, default=max(1, (os.cpu_count() or 2) - 6))
    a = ap.parse_args(argv)
    a.out.mkdir(parents=True, exist_ok=True)
    jobs = [(snr, tau, a.seed0 + i) for tau in TAUS_S for snr in SNRS_DB for i in range(a.seeds)]
    with ProcessPoolExecutor(a.jobs) as pool:
        checks_running = pool.map(noise_free, CUTS)
        rows = list(pool.map(trial, jobs, chunksize=2))
        checks = list(checks_running)
    with open(a.out / "echoes.csv", "w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    text = _report(summarise(rows), checks)
    (a.out / "summary.txt").write_text(text, encoding="utf-8")
    print(text)


if __name__ == "__main__":
    main()
