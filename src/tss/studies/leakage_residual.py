"""Layer: studies
WP: WP3
Claims: C4
Status: implemented
Requirements: INS-10, INS-11, INS-23, INS-24

WP3 Phase 4: the residual leakage in the echo band, step by step.

The default instrument since Phase 4 carries the leakage (the moment-method
isolation of the crossed, offset dipoles on the 6U bus), the transmitter's
allocations at threshold (the baseline oscillator, a -123 dBc/Hz phase-noise
floor, -133 dBc/Hz of AM noise and two PM spurs at INS-22's -82 dBc), and runs
primed: the transmitter runs between sweeps (decisions of 2026-10-28).  For
each step of the A1 and A2 sweeps this study records what the L2 baseline
scores INS-24 on (:func:`tss.l2_instrument.leakage.dechirped_residual`): the
residual's pedestal, any line, and its margin below the weakest echo INS-10
admits, with the R_max reading and the ideal-canceller bound beside it.

Cases: the default (``alloc_primed``); a cold start (``alloc_cold``); the
default without AM noise (``am_off``), which shows what sets the pedestal; and
the default on a flat 40 dB isolation (``flat40``), INS-23's threshold.

Two quantities follow per step, each read as "what would meet INS-24 here":

* ``am_needed_dbc_hz``: the AM density at which the pedestal would sit 10 dB
  below the weakest echo, ``am - (10 - margin)``.  It holds where AM sets the
  pedestal, which the ``am_off`` case shows is everywhere for the default;
* ``canceller_needed_db``: the suppression ahead of the receiver that would
  close the margin, ``max(0, 10 - margin)``.
"""

from __future__ import annotations

import copy
import json
from pathlib import Path
from typing import Any

import numpy as np

from tss.studies import StudyResult, register_study

INS24_MARGIN_DB = 10.0
SUMMARY_MHZ = (0.33, 0.5, 1.0, 1.26, 2.0, 3.0, 5.0, 10.0, 20.0)
CASES = ("alloc_primed", "alloc_cold", "am_off", "flat40")


def _cfg(case: str) -> dict[str, Any]:
    from tss.config import load_yaml, resolve_path

    cfg = copy.deepcopy(load_yaml(resolve_path("instrument/baseline_6u.yaml")))
    cfg["leakage"]["enabled"] = True
    # Pinned to the Phase 4 default, which this study records: the masks' AM
    # and no canceller (Phase 5 moved the default to INS-30 and INS-31).
    cfg["transmit_impairments"]["am_noise_dbc_hz"] = -133.0
    cfg["canceller"] = {"enabled": False}
    if case == "alloc_cold":
        cfg["leakage"]["primed"] = False
    elif case == "am_off":
        cfg["transmit_impairments"]["am_noise_dbc_hz"] = None
    elif case == "flat40":
        cfg["leakage"]["profile"] = "flat"
    return cfg


def _no_echoes() -> Any:
    from tss.interfaces import EchoSet

    e = np.array([], dtype=float)
    return EchoSet(freq_hz=e, mode=np.array([], dtype=object), group_delay_s=e, amplitude_db=e,
                   absorption_db=e, flags=np.array([], dtype=object), source="none")


def per_step(case: str, sweep: Any, seed: int = 0) -> tuple[dict[str, np.ndarray], dict[str, Any]]:
    """Per-step residual and margins, and the record's L2 metrics, for one case."""
    from tss.l2_instrument.baseline import BaselineInstrument

    cfg = _cfg(case)
    inst = BaselineInstrument.from_config(cfg)
    rec = inst.receive(_no_echoes(), sweep, np.random.default_rng(seed))
    ps = {k: np.asarray(v, dtype=float) for k, v in rec.flags["leakage_per_step"].items()}
    m = dict(rec.flags["l2_metrics"])
    am = inst.impairments.am_noise_dbc_hz
    margin = ps["ins24_margin_db"]
    ps["noise_cell_dbw"] = ps["weakest_echo_dbw"] - 10.0
    ps["pedestal_over_noise_db"] = ps["residual_pedestal_dbw"] - ps["noise_cell_dbw"]
    ps["canceller_needed_db"] = np.maximum(0.0, INS24_MARGIN_DB - margin)
    ps["am_needed_dbc_hz"] = (am - (INS24_MARGIN_DB - margin)) if am is not None else np.full(margin.shape, np.nan)
    return ps, m


def compute(seed: int = 0) -> dict[str, Any]:
    from tss.studies.data_volume import A1, A2

    out: dict[str, Any] = {}
    for name, cfgs in (("A1", A1), ("A2", A2)):
        sweep = cfgs.schedule()
        out[name] = {case: per_step(case, sweep, seed) for case in CASES}
    return out


def _metrics(ps: dict[str, np.ndarray], m: dict[str, Any]) -> dict[str, Any]:
    keep = (
        "isolation_db_min", "isolation_db_min_freq_hz", "leakage_residual_margin_db",
        "leakage_residual_margin_freq_hz", "leakage_residual_margin_db_rmax",
        "leakage_residual_steps_failing", "leakage_residual_steps_scored",
        "leakage_residual_margin_db_canceller", "leakage_residual_steps_failing_canceller",
        "dynamic_range_db_at_digitiser", "ins11_rejection_margin_db",
        "ins11_residual_peak_over_clip_db",
    )
    out = {k: m.get(k) for k in keep}
    fin = np.isfinite(ps["am_needed_dbc_hz"])
    out["am_needed_dbc_hz_min"] = float(np.min(ps["am_needed_dbc_hz"][fin])) if np.any(fin) else None
    out["canceller_needed_db_max"] = float(np.nanmax(ps["canceller_needed_db"]))
    out["pedestal_over_noise_db_max"] = float(np.nanmax(ps["pedestal_over_noise_db"]))
    out["steps_with_a_line"] = int(np.count_nonzero(np.isfinite(ps["residual_line_dbw"])))
    return out


def _table(ps: dict[str, np.ndarray]) -> str:
    rows = [
        "| f (MHz) | iso (dB) | pedestal (dBW) | over noise (dB) | weakest echo (dBW) | margin (dB) "
        "| margin, R_max (dB) | margin, canceller (dB) | AM needed (dBc/Hz) |",
        "|---|---|---|---|---|---|---|---|---|",
    ]
    f = ps["f_hz"] / 1e6
    canc = ps["weakest_echo_dbw"] - ps.get("residual_level_dbw_canceller", np.full(f.shape, np.nan))
    for target in SUMMARY_MHZ:
        i = int(np.argmin(np.abs(f - target)))
        rows.append(
            f"| {f[i]:.3f} | {ps['isolation_db'][i]:.1f} | {ps['residual_pedestal_dbw'][i]:.1f} | "
            f"{ps['pedestal_over_noise_db'][i]:.1f} | {ps['weakest_echo_dbw'][i]:.1f} | "
            f"{ps['ins24_margin_db'][i]:.1f} | {ps['ins24_margin_db_rmax'][i]:.1f} | {canc[i]:.1f} | "
            f"{ps['am_needed_dbc_hz'][i]:.1f} |"
        )
    return "\n".join(rows)


@register_study(
    "leakage_residual",
    description="WP3 Phase 4: residual leakage in the echo band and the INS-24 margin, per step, A1 and A2",
    work_packages=["WP3"],
    claims=["C4"],
    status="implemented",
    requirements=["INS-10", "INS-11", "INS-23", "INS-24"],
)
def run(out_dir: Path, *, config: dict[str, Any], manifest: Any) -> StudyResult:
    res = compute(int(config.get("seed", 0)))
    outputs: list[Path] = []
    metrics: dict[str, Any] = {}
    parts: list[str] = []
    for name in ("A1", "A2"):
        for case, (ps, m) in res[name].items():
            cols = [k for k in ps if ps[k].size == ps["f_hz"].size]
            path = out_dir / f"{name}_{case}.csv"
            np.savetxt(path, np.column_stack([ps[c] for c in cols]), delimiter=",", header=",".join(cols),
                       comments="", fmt="%.6g")
            outputs.append(path)
            metrics[f"{name}_{case}"] = _metrics(ps, m)
        parts.append(f"## {name}, default (allocations at threshold, primed)\n\n{_table(res[name]['alloc_primed'][0])}\n")
    a = metrics["A1_alloc_primed"]
    off = metrics["A1_am_off"]
    flat = metrics["A1_flat40"]
    summary = "\n".join([
        "Residual leakage in each step's echo band, scored against the weakest echo INS-10 admits "
        "(INS-24, 10 dB), with no canceller.",
        "",
        f"* Default, A1: margin {a['leakage_residual_margin_db']:.1f} dB at "
        f"{a['leakage_residual_margin_freq_hz'] / 1e6:.3f} MHz; {a['leakage_residual_steps_failing']} of "
        f"{a['leakage_residual_steps_scored']} steps fail.  R_max reading {a['leakage_residual_margin_db_rmax']:.1f} dB.  "
        f"Pedestal up to {a['pedestal_over_noise_db_max']:.1f} dB above the noise in a cell.",
        f"* Behind an ideal 30 dB canceller: margin {a['leakage_residual_margin_db_canceller']:.1f} dB, "
        f"{a['leakage_residual_steps_failing_canceller']} steps failing.  Canceller needed: "
        f"{a['canceller_needed_db_max']:.1f} dB.",
        f"* Without AM noise: margin {off['leakage_residual_margin_db']:.1f} dB.  AM sets the pedestal; "
        f"the AM density that would meet INS-24 at every step is {a['am_needed_dbc_hz_min']:.1f} dBc/Hz.",
        f"* On a flat 40 dB isolation: margin {flat['leakage_residual_margin_db']:.1f} dB.",
        f"* INS-11 verified at {a['dynamic_range_db_at_digitiser']:.0f} dB; least rejection margin on the noise "
        f"rule {a['ins11_rejection_margin_db']:.1f} dB.  INS-23: {a['isolation_db_min']:.1f} dB at "
        f"{a['isolation_db_min_freq_hz'] / 1e6:.2f} MHz.",
        "",
        *parts,
    ])
    path = out_dir / "leakage_residual.json"
    path.write_text(json.dumps({"metrics": metrics}, indent=2, default=str), encoding="utf-8")
    outputs.append(path)
    return StudyResult(outputs=outputs, summary=summary, metrics=metrics)
