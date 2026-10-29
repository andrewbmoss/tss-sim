"""Layer: studies
WP: WP3
Claims: C4
Status: implemented
Requirements: INS-11, INS-23, INS-24

WP3 Phase 3: the leakage through the receive chain, step by step.

For each step of the A1 and A2 sweeps, with the leakage on the moment-method
isolation of the dipoles on the 6U bus (the default when leakage is on,
decision of 2026-10-28) and, for comparison, flat at INS-23's 40 dB:

* the leakage at the receiver input, and its ratio to the noise in the
  digitiser band and to the free-space echo from 1050 km;
* INS-11's required high-pass rejection on the noise-referenced rule
  ``L/N_fs - 60 + 10`` (decision G5b, 2026-10-28) and on Rev G's rule
  ``L/E - 60 + 10``, beside the rejection achieved: the leakage power at the
  filter output in the last quarter of the step over its input power;
* the settling time: from the start of the step until the filtered leakage
  stays below the noise in the digitiser band, against the dwell and the
  echo's own transition time ``2 R_max / c``;
* the front end: the leakage's margin below the input 1 dB point (WP9-2 asks
  10 dB), the desensitisation of an echo beside it, and the image the
  compression puts at the negative beat;
* the high-pass corner and the echo loss at 100 km and 1050 km.

The chain is run cold, with the transmitter switching on at the start of the
sweep and the filter at rest, and primed, with the transmitter already on and
the filter settled, as it would be if the instrument transmitted continuously
between sweeps.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from tss.studies import StudyResult, register_study

INS11_DR_DB = 60.0
HEADROOM_DB = 10.0
WP9_2_MARGIN_DB = 10.0
SUMMARY_MHZ = (0.33, 0.5, 1.0, 1.26, 2.0, 3.0, 5.0, 10.0, 20.0)


def _instrument(**overrides: Any) -> Any:
    from tss.config import load_yaml, resolve_path
    from tss.l2_instrument.baseline import BaselineInstrument

    cfg = load_yaml(resolve_path("instrument/baseline_6u.yaml"))
    cfg = json.loads(json.dumps(cfg))
    cfg.setdefault("leakage", {})["enabled"] = True
    # The receive chain is studied on an ideal transmitter, as it was when this
    # study was written (Phase 3); the allocations the default carries since
    # Phase 4 are the business of the leakage_residual study.
    cfg["transmit_impairments"] = {}
    cfg["canceller"] = {"enabled": False}
    return BaselineInstrument.from_config(cfg, **overrides)


def _steady_state(hp: Any, x: np.ndarray, sweep: Any, fs: float) -> np.ndarray:
    """The filter's state after a long run, at the first step's corner, on the
    leakage as it is at the start of the sweep."""
    import dataclasses

    from tss import sweep_law

    first = dataclasses.replace(hp, fixed_corner_hz=float(hp.corner_hz(sweep_law.rate_hz_per_s(sweep, 0.0))))
    pre = np.full(int(0.5 * fs), x[np.flatnonzero(np.abs(x) > 0)[0]])
    _, st = first.apply(pre, sweep, fs, t0_s=0.0, return_state=True)
    return st


def per_step(inst: Any, sweep: Any, *, model: Any | None = None, primed: bool = False) -> dict[str, np.ndarray]:
    """Per-step quantities for one instrument, sweep and leakage model."""
    from tss import sweep_law
    from tss.constants import C_LIGHT, K_BOLTZ
    from tss.l2_instrument import leakage as lk
    from tss.l2_instrument import receive_chain as rc

    model = model if model is not None else inst.leakage
    fs = inst.sample_rate_for(sweep)
    dwell = float(sweep.step_dwell_s)
    n = int(math.floor(sweep_law.active_duration_s(sweep) * fs))
    rng = np.random.default_rng(0)
    leak = lk.leakage_signal(sweep, fs, model, rng, n_samples=n, transmit_power_w=inst.tx_power_w,
                             tx_chain_efficiency_db=inst.tx_chain_efficiency_db,
                             impairments=inst.impairments)
    if primed:
        # The transmitter was already on: the leakage does not switch on at t = tau.
        first = int(np.flatnonzero(np.abs(leak) > 0)[0])
        leak = leak.copy()
        leak[:first] = leak[first]
    x = inst.front_end.apply(leak) if not inst.front_end.linear else leak
    hp = inst.highpass
    state = _steady_state(hp, x, sweep, fs) if primed else None
    y = hp.apply(x, sweep, fs, state=state)

    centres = sweep_law.step_centres_hz(sweep)
    nsteps = centres.size
    t = np.arange(n) / fs
    idx = sweep_law.step_index(sweep, t)
    f_tx = sweep_law.frequency_hz(sweep, t)
    t_sys = inst.antenna_temperature_k(np.clip(f_tx, 1.0, None) / 1e6) + inst.receiver.noise_temperature_k
    n_fs = K_BOLTZ * t_sys * fs                                    # noise in the digitiser band, W

    lev = lk.leakage_over_echo_db(inst, sweep, model)
    out: dict[str, list[float]] = {k: [] for k in (
        "f_hz", "rate_hz_per_s", "corner_hz", "isolation_db", "leak_in_dbm", "leak_over_noise_fs_db",
        "leak_over_rmax_echo_db", "required_noise_rule_db", "required_revg_rule_db", "achieved_db",
        "settle_s", "settle_over_dwell", "echo_transition_s", "residual_peak_over_clip_db",
        "fe_margin_db", "desens_db", "image_dbc", "echo_loss_100km_db", "echo_loss_1050km_db")}
    for i in range(nsteps):
        sel = idx == i
        if not np.any(sel):
            continue
        f_c = float(sweep_law.frequency_hz(sweep, t[sel][len(t[sel]) // 2]))
        rate = float(sweep_law.rate_hz_per_s(sweep, t[sel][0]))
        corner = float(hp.corner_hz(rate))
        p_in = np.abs(x[sel]) ** 2
        p_out = np.abs(y[sel]) ** 2
        tail = slice(int(0.75 * p_in.size), p_in.size)
        nfs = float(np.median(n_fs[sel]))
        j = int(np.argmin(np.abs(lev["f_hz"] - f_c)))
        lin = float(np.mean(p_in[tail]))
        lout = float(np.mean(p_out[tail]))
        above = np.flatnonzero(p_out > nfs)
        settle = (above[-1] + 1) / fs if above.size else 0.0
        clip_w = nfs * 10 ** (INS11_DR_DB / 10)
        r_fe = math.sqrt(float(np.mean(np.abs(leak[sel]) ** 2)))
        fe = inst.front_end
        out["f_hz"].append(f_c)
        out["rate_hz_per_s"].append(rate)
        out["corner_hz"].append(corner)
        out["isolation_db"].append(float(lev["isolation_db"][j]))
        out["leak_in_dbm"].append(10 * math.log10(lin) + 30)
        out["leak_over_noise_fs_db"].append(10 * math.log10(lin / nfs))
        out["leak_over_rmax_echo_db"].append(float(lev["leakage_over_rmax_echo_db"][j]))
        out["required_noise_rule_db"].append(10 * math.log10(lin / nfs) - INS11_DR_DB + HEADROOM_DB)
        out["required_revg_rule_db"].append(float(lev["leakage_over_rmax_echo_db"][j]) - INS11_DR_DB + HEADROOM_DB)
        out["achieved_db"].append(10 * math.log10(lin / max(lout, 1e-300)))
        out["settle_s"].append(settle)
        out["settle_over_dwell"].append(settle / dwell)
        out["echo_transition_s"].append(2.0 * float(inst.design_range_m) / C_LIGHT)
        out["residual_peak_over_clip_db"].append(10 * math.log10(max(float(p_out.max()), 1e-300) / clip_w))
        if fe.linear:
            out["fe_margin_db"].append(float("inf"))
            out["desens_db"].append(0.0)
            out["image_dbc"].append(float("-inf"))
        else:
            out["fe_margin_db"].append(float(fe.input_p1db_dbm) - (20 * math.log10(r_fe) + 30))
            out["desens_db"].append(rc.small_signal_gain_db(fe, r_fe**2))
            out["image_dbc"].append(rc.compression_image_dbc(fe, r_fe**2))
        out["echo_loss_100km_db"].append(-float(hp.response_db(2 * rate * 1.0e5 / C_LIGHT, corner)[0]))
        out["echo_loss_1050km_db"].append(-float(hp.response_db(2 * rate * 1.05e6 / C_LIGHT, corner)[0]))
    return {k: np.asarray(v, dtype=float) for k, v in out.items()}


def compute() -> dict[str, Any]:
    from tss.l2_instrument import leakage as lk
    from tss.studies.data_volume import A1, A2

    inst = _instrument()
    flat = lk.LeakageModel.from_config(inst.isolation)
    out: dict[str, Any] = {}
    for name, cfgs in (("A1", A1), ("A2", A2)):
        sweep = cfgs.schedule()
        out[name] = {
            "mom_cold": per_step(inst, sweep),
            "mom_primed": per_step(inst, sweep, primed=True),
            "flat40_cold": per_step(inst, sweep, model=flat),
        }
    out["front_end"] = dict(inst.front_end.describe())
    out["highpass"] = dict(inst.highpass.describe())
    out["leakage_model"] = dict(inst.leakage.describe())
    return out


def _table(r: dict[str, np.ndarray]) -> str:
    rows = [
        "| f (MHz) | corner (Hz) | iso (dB) | L in (dBm) | L/N_fs (dB) | need, noise rule (dB) | need, Rev G rule (dB) "
        "| achieved (dB) | settle (ms) | settle / dwell | FE margin (dB) | desens (dB) | image (dBc) | loss at 100 km (dB) |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    f = r["f_hz"] / 1e6
    for target in SUMMARY_MHZ:
        i = int(np.argmin(np.abs(f - target)))
        rows.append(
            f"| {f[i]:.3f} | {r['corner_hz'][i]:.0f} | {r['isolation_db'][i]:.1f} | {r['leak_in_dbm'][i]:.1f} | "
            f"{r['leak_over_noise_fs_db'][i]:.1f} | {r['required_noise_rule_db'][i]:.1f} | "
            f"{r['required_revg_rule_db'][i]:.1f} | {min(r['achieved_db'][i], 300):.0f} | "
            f"{1e3 * r['settle_s'][i]:.2f} | {r['settle_over_dwell'][i]:.3f} | {r['fe_margin_db'][i]:.1f} | "
            f"{r['desens_db'][i]:.3f} | {r['image_dbc'][i]:.1f} | {r['echo_loss_100km_db'][i]:.3f} |"
        )
    return "\n".join(rows)


def _metrics(r: dict[str, np.ndarray]) -> dict[str, Any]:
    margin = r["achieved_db"] - r["required_noise_rule_db"]
    return {
        "steps": int(r["f_hz"].size),
        "ins11_noise_rule_met_all_steps": bool(np.all(margin >= 0)),
        "ins11_noise_rule_min_margin_db": float(np.min(margin)),
        "ins11_revg_rule_met_all_steps": bool(np.all(r["achieved_db"] >= r["required_revg_rule_db"])),
        "required_noise_rule_max_db": float(np.max(r["required_noise_rule_db"])),
        "required_revg_rule_max_db": float(np.max(r["required_revg_rule_db"])),
        "settle_max_ms": float(1e3 * np.max(r["settle_s"])),
        "settle_max_over_dwell": float(np.max(r["settle_over_dwell"])),
        "steps_settling_longer_than_echo_transition": int(np.count_nonzero(r["settle_s"] > r["echo_transition_s"])),
        "first_step_settle_ms": float(1e3 * r["settle_s"][0]),
        "fe_margin_min_db": float(np.min(r["fe_margin_db"])),
        "steps_below_wp9_2_margin": int(np.count_nonzero(r["fe_margin_db"] < WP9_2_MARGIN_DB)),
        "desens_worst_db": float(np.min(r["desens_db"])),
        "image_worst_dbc": float(np.max(r["image_dbc"])),
        "residual_peak_over_clip_max_db": float(np.max(r["residual_peak_over_clip_db"])),
        "echo_loss_100km_max_db": float(np.max(r["echo_loss_100km_db"])),
        "corner_range_hz": [float(np.min(r["corner_hz"])), float(np.max(r["corner_hz"]))],
    }


@register_study(
    "leakage_receive_chain",
    description="WP3 Phase 3: leakage through the front end and the leakage high-pass, per step, A1 and A2",
    work_packages=["WP3"],
    claims=["C4"],
    status="implemented",
    requirements=["INS-11", "INS-23", "INS-24"],
)
def run(out_dir: Path, *, config: dict[str, Any], manifest: Any) -> StudyResult:
    res = compute()
    outputs: list[Path] = []
    metrics: dict[str, Any] = {}
    parts: list[str] = []
    for name in ("A1", "A2"):
        for case, r in res[name].items():
            cols = [k for k in r if r[k].size == r["f_hz"].size]
            path = out_dir / f"{name}_{case}.csv"
            np.savetxt(path, np.column_stack([r[c] for c in cols]), delimiter=",", header=",".join(cols),
                       comments="", fmt="%.6g")
            outputs.append(path)
            metrics[f"{name}_{case}"] = _metrics(r)
        parts.append(f"## {name}, MoM bus isolation, filter cold\n\n{_table(res[name]['mom_cold'])}\n")
    m = metrics["A1_mom_cold"]
    summary = "\n".join([
        "Leakage on the moment-method isolation of the dipoles on the 6U bus, front end "
        f"{res['front_end']['input_p1db_dbm']} dBm P1dB (Rapp), high-pass order {res['highpass']['order']} "
        f"at {res['highpass']['corner_fraction']} of each step's 100 km beat.",
        "",
        f"* INS-11, noise-referenced rule (G5b): {'met' if m['ins11_noise_rule_met_all_steps'] else 'NOT met'} "
        f"at every step of A1, least margin {m['ins11_noise_rule_min_margin_db']:.0f} dB; the rule asks at most "
        f"{m['required_noise_rule_max_db']:.1f} dB.  Rev G's rule asks at most {m['required_revg_rule_max_db']:.1f} dB.",
        f"* Settling, cold: longest {m['settle_max_ms']:.2f} ms ({m['settle_max_over_dwell']:.2f} of a dwell); "
        f"first step {m['first_step_settle_ms']:.2f} ms; steps settling longer than the 7.0 ms echo transition: "
        f"{m['steps_settling_longer_than_echo_transition']}.",
        f"* Front end: least margin below P1dB {m['fe_margin_min_db']:.1f} dB; steps under WP9-2's 10 dB: "
        f"{m['steps_below_wp9_2_margin']}; worst echo desensitisation {m['desens_worst_db']:.3f} dB; worst image "
        f"{m['image_worst_dbc']:.1f} dBc.",
        f"* Echo loss at 100 km through the high-pass: at most {m['echo_loss_100km_max_db']:.3f} dB.",
        "",
        *parts,
    ])
    path = out_dir / "leakage_receive_chain.json"
    path.write_text(json.dumps({"metrics": metrics, "front_end": res["front_end"], "highpass": res["highpass"],
                                "leakage_model": res["leakage_model"]}, indent=2, default=str), encoding="utf-8")
    outputs.append(path)
    return StudyResult(outputs=outputs, summary=summary, metrics=metrics)
