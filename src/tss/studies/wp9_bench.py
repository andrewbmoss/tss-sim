"""Layer: studies
WP: WP9
Claims: C4
Status: implemented
Requirements: INS-11, INS-24, INS-29, INS-30, INS-31

WP3 Phase 6: predictions and pass limits for the WP9 conducted bench.

WP9's conducted bench (Technical Volume Ch. 10; BOM Rev 1; the WP9 Bench Test
Specification) replaces the bus coupling with a coax delay line and calibrated
attenuators.  In Phase 1a it emulates INS-31's canceller by attenuation.  An
ideal canceller of depth D at isolation I leaves the leakage, with the AM it
carries, I + D below the transmitter, and a flat path of I + D does the same.
This study gives, per step, what the bench should see, so that every bench
result has a prediction to be compared with (WP9-9):

* the INS-29 limit k T_sys at the receiver input, in dBm/Hz, that a bench
  result is referred to;
* the leakage-borne noise density at the receiver input, for each bench path
  and each length of the coax delay set, at INS-30's threshold and goal AM;
* how far that noise lifts the bench's own floor, which is what the bench
  resolves;
* the front-end drive, margin and compression image for each path;
* the limits on stray coupling (WP9-1), and the weakest echo a channel of the
  synthesiser can provide before its crosstalk, converted by the leakage
  carrier, exceeds a tenth of it.

The bench paths run from the PA output at 1 W to the receiver input:

* E0, 39 dB: no canceller.  This is the most leakage a level-7 mixer takes
  with WP9-2's 10 dB margin.
* E1, 61 dB: INS-23's threshold (35 dB) plus INS-31's threshold (26 dB).
* E2, 69 dB: INS-23's threshold plus INS-31's goal (34 dB).

METHOD
------
The method is the leakage_budget study's.  The L2 chain is run with the
leakage on a flat path at two AM levels with matched draws.  That separates
exactly the AM pedestal ``X`` and the rest ``Y``, per unit leakage, for each
coupling delay.  The pedestal is INS-29's (``residual_pedestal_dbw``).  The
leakage-borne noise density at the receiver input is then
``P_leak (X 10^((AM + 133)/10) + Y) / B_cell``.

The bench's own noise floor is an estimate until WP9's floor measurement
replaces it.  It assumes a ZAD-6+ mixer (4.65 dB conversion loss) followed by
a low-noise baseband stage, about 8.6 dB of noise figure in all, so
-165 dBm/Hz.  The transfer from the AM to the beat is the same in the bench
receiver as in the simulated one, since the de-chirp passes AM unchanged.
The bench does not reproduce the flight noise level; it measures the
leakage-borne density in dBm/Hz and refers it to the flight k T_sys.
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from tss.studies import StudyResult, register_study

#: Transmitter output at the PA, dBm (INS-14's 1 W).
P_TX_DBM = 30.0
#: Velocity factor of the RG-58 delay set.
RG58_VF = 0.66
#: The coax delay set (BOM A16): RG-58 lengths for 1, 3, 10 and 30 m of free space.
DELAY_SET_M = (0.66, 1.98, 6.6, 19.8)
#: The delay used for the headline tables: 3 m of free space.
HEADLINE_DELAY_M = 1.98
#: Bench paths, PA output to receiver input, dB.
BENCH_PATHS_DB = {"E0": 39.0, "E1": 61.0, "E2": 69.0}
#: What each path stands for.
BENCH_PATH_MEANING = {
    "E0": "no canceller; level-7 mixer at WP9-2's margin",
    "E1": "INS-23 threshold 35 dB + INS-31 threshold 26 dB",
    "E2": "INS-23 threshold 35 dB + INS-31 goal 34 dB",
}
#: INS-30 threshold and goal, dBc/Hz single sideband.
AM_CASES_DBC_HZ = (-140.0, -150.0)
#: Estimated bench receiver floor at the receiver input, dBm/Hz.
BENCH_FLOOR_DBM_HZ = -165.0
#: The same with the ZFL-500LN+ ahead of the mixer (NF 2.9 dB), allowed by
#: WP9-2 only below about -26 dBm of leakage, so only behind a canceller.
BENCH_FLOOR_LNA_DBM_HZ = -171.0
#: Mixer input P1dB: level-7 (ZAD-6+) and level-17 (BOM B02).
MIXER_P1DB_DBM = (1.0, 11.0)
WP9_2_MARGIN_DB = 10.0
#: Echo SNR on the bench: INS-10's weakest admitted echo, and a strong echo.
WEAKEST_SNR_DB = 10.0
STRONG_SNR_DB = 60.0
#: AD9959 channel-to-channel isolation (data sheet: "> 65 dB").
AD9959_ISOLATION_DB = 65.0
#: DDS channel output level assumed for the de-chirp LO and an echo channel, dBm.
DDS_CHANNEL_DBM = 0.0
STRAY_MARGIN_DB = 10.0
#: Averages assumed for a two-chain cross-spectrum AM measurement.
CROSS_SPECTRUM_AVERAGES = 1000
#: INS-31 weight-step limit as a fraction of the leakage: 0.001 dB or 0.0066 degrees.
INS31_JUMP_DBC = 20.0 * math.log10(10.0 ** (0.001 / 20.0) - 1.0)
#: INS-31 added-noise limit at the receiver input, dBm/Hz.
INS31_ADDED_NOISE_DBM_HZ = -170.0
SEEDS = (0, 1, 2)


def coax_delay_s(length_m: float) -> float:
    """One-way delay of an RG-58 length."""
    from tss.constants import C_LIGHT

    return float(length_m) / (RG58_VF * C_LIGHT)


def _lin(db: Any) -> Any:
    return 10.0 ** (np.asarray(db, dtype=float) / 10.0)


def _db(x: Any) -> Any:
    return 10.0 * np.log10(np.maximum(np.asarray(x, dtype=float), 1e-300))


def characterise(sweep: Any, delay_s: float, seed: int) -> dict[str, np.ndarray]:
    """``X`` and ``Y`` of INS-29's pedestal, per unit leakage, on a flat path of one delay."""
    from tss.l2_instrument.reference_echo import reference_levels
    from tss.studies import leakage_budget as lb

    def cfg(am: float) -> dict[str, Any]:
        c = lb._cfg(profile="flat", am=am, canceller="none")
        c["isolation"]["coupling_delay_s"] = float(delay_s)
        return c

    on, _, inst = lb._run(cfg(lb.AM_REF_DBC_HZ), sweep, seed)
    low, _, _ = lb._run(cfg(lb.AM_LOW_DBC_HZ), sweep, seed)
    p_leak = _lin(-on["isolation_db"]) * float(inst.tx_power_w) * 10.0 ** (inst.tx_chain_efficiency_db / 10.0)
    a, b = _lin(on["residual_pedestal_dbw"]), _lin(low["residual_pedestal_dbw"])
    r = 10.0 ** ((lb.AM_LOW_DBC_HZ - lb.AM_REF_DBC_HZ) / 10.0)
    x = np.maximum(a - b, 0.0) / (1.0 - r)
    lev = reference_levels(inst, sweep)
    return {
        "f_hz": np.asarray(lev["f_hz"], dtype=float),
        "X": x / p_leak,
        "Y": np.maximum(a - x, 0.0) / p_leak,
        "cell_bandwidth_hz": np.asarray(lev["cell_bandwidth_hz"], dtype=float),
        "noise_cell_dbw": np.asarray(lev["noise_cell_dbw"], dtype=float),
        "beat_min_hz": np.asarray(lev["beat_min_hz"], dtype=float),
        "beat_max_hz": np.asarray(lev["beat_max_hz"], dtype=float),
        "t_sys_k": np.asarray(lev["t_sys_k"], dtype=float),
    }


def kt_sys_dbm_hz(ch: dict[str, np.ndarray]) -> np.ndarray:
    """INS-29's limit, k T_sys at the receiver input, dBm/Hz."""
    return ch["noise_cell_dbw"] - _db(ch["cell_bandwidth_hz"]) + 30.0


def leakage_noise_dbm_hz(ch: dict[str, np.ndarray], path_db: float, am_dbc_hz: float,
                         p_tx_dbm: float = P_TX_DBM) -> np.ndarray:
    """Leakage-borne noise density at the receiver input behind a flat path, dBm/Hz."""
    p_leak_w = 10.0 ** ((p_tx_dbm - float(path_db) - 30.0) / 10.0)
    ped = p_leak_w * (ch["X"] * 10.0 ** ((float(am_dbc_hz) - (-133.0)) / 10.0) + ch["Y"])
    return _db(ped / ch["cell_bandwidth_hz"]) + 30.0


def floor_rise_db(noise_dbm_hz: Any, floor_dbm_hz: float = BENCH_FLOOR_DBM_HZ) -> Any:
    """How far a noise density lifts a floor: ``10 log10(1 + N/F)``."""
    return 10.0 * np.log10(1.0 + _lin(np.asarray(noise_dbm_hz, dtype=float) - float(floor_dbm_hz)))


def am_allowed_dbc_hz(ch: dict[str, np.ndarray], path_db: float, p_tx_dbm: float = P_TX_DBM) -> np.ndarray:
    """Per step, the AM at which the leakage-borne noise behind ``path_db`` equals k T_sys."""
    p_leak_w = 10.0 ** ((p_tx_dbm - float(path_db) - 30.0) / 10.0)
    allow = _lin(kt_sys_dbm_hz(ch) - 30.0) * ch["cell_bandwidth_hz"] - p_leak_w * ch["Y"]
    with np.errstate(divide="ignore", invalid="ignore"):
        need = -133.0 + _db(allow / np.maximum(p_leak_w * ch["X"], 1e-300))
    return np.where(allow > 0.0, need, np.nan)


def am_measurement_floor_dbc_hz(floor_dbm_hz: float = BENCH_FLOOR_DBM_HZ,
                                carrier_dbm: float = MIXER_P1DB_DBM[0] - WP9_2_MARGIN_DB,
                                averages: int | None = None) -> float:
    """Single-sideband AM floor of in-phase homodyne detection through the bench receiver.

    The AM sidebands of a carrier ``P_c`` add coherently in the in-phase
    channel (power ``2 L P_c``), and half of the receiver noise from the two
    sidebands lands there (``N0``), so the floor is ``N0 / (2 P_c)``.  A
    cross-spectrum of two independent chains, each with the same carrier,
    lowers it by ``5 log10(m)`` for ``m`` averages.
    """
    floor = float(floor_dbm_hz) - float(carrier_dbm) - 10.0 * math.log10(2.0)
    if averages:
        floor -= 5.0 * math.log10(float(averages))
    return floor


def stray_limits(cell_bandwidth_hz: float) -> dict[str, float]:
    """WP9-1: stray-coupling limits on the bench."""
    weakest = BENCH_FLOOR_DBM_HZ + 10.0 * math.log10(cell_bandwidth_hz) + WEAKEST_SNR_DB
    out = {
        "weakest_echo_dbm": weakest,
        "strong_echo_dbm": weakest + STRONG_SNR_DB - WEAKEST_SNR_DB,
        # The transmit chain's stray is more leakage: it must sit 10 dB below the weakest emulated path.
        "tx_stray_carrier_db": max(BENCH_PATHS_DB.values()) + STRAY_MARGIN_DB,
        # ... and the noise it carries 10 dB below the bench floor (AM at INS-30's threshold).
        "tx_stray_noise_db": P_TX_DBM + AM_CASES_DBC_HZ[0] - BENCH_FLOOR_DBM_HZ + STRAY_MARGIN_DB,
    }
    for name, path in BENCH_PATHS_DB.items():
        p_leak = P_TX_DBM - path
        # Crosstalk from an echo channel into the TX or LO channel, converted by
        # the leakage carrier: spurious/echo = P_leak - iso - P_LO + L_echo <= -10 dB.
        l_echo_max = DDS_CHANNEL_DBM - p_leak + AD9959_ISOLATION_DB - STRAY_MARGIN_DB
        out[f"dds_echo_path_max_db_{name}"] = l_echo_max
        out[f"dds_echo_min_dbm_{name}"] = DDS_CHANNEL_DBM - l_echo_max
        out[f"dds_echo_min_snr_db_{name}"] = DDS_CHANNEL_DBM - l_echo_max - (weakest - WEAKEST_SNR_DB)
    return out


def front_end(path_db: float) -> dict[str, float]:
    from tss.studies.leakage_budget import front_end_metrics

    p_w = np.array([10.0 ** ((P_TX_DBM - float(path_db) - 30.0) / 10.0)])
    out: dict[str, float] = {"leakage_dbm": P_TX_DBM - float(path_db)}
    for p1 in MIXER_P1DB_DBM:
        m = front_end_metrics(p_w, p1)
        out[f"margin_db_p1db_{p1:+.0f}"] = m["fe_margin_db"]
        out[f"image_dbc_p1db_{p1:+.0f}"] = m["fe_image_dbc"]
    return out


def compute(*, seeds: tuple[int, ...] = SEEDS) -> dict[str, Any]:
    from tss.studies.data_volume import A1, A2

    sweeps = {"A1": A1.schedule(), "A2": A2.schedule()}
    rows: list[dict[str, Any]] = []
    limits: dict[str, list[dict[str, Any]]] = {}
    metrics: dict[str, Any] = {}
    for sname, sweep in sweeps.items():
        chars = {(d, s): characterise(sweep, coax_delay_s(d), s) for d in DELAY_SET_M for s in seeds}
        ref = chars[(DELAY_SET_M[0], seeds[0])]
        kt = kt_sys_dbm_hz(ref)
        n_steps = kt.size
        worst: dict[tuple[str, float, float], np.ndarray] = {}
        for d in DELAY_SET_M:
            for path, pdb in BENCH_PATHS_DB.items():
                for am in AM_CASES_DBC_HZ:
                    n = np.nanmax(np.vstack([leakage_noise_dbm_hz(chars[(d, s)], pdb, am) for s in seeds]), axis=0)
                    worst[(path, am, d)] = n
                    for i in range(n_steps):
                        rows.append({
                            "sweep": sname, "step": i + 1, "f_hz": float(ref["f_hz"][i]),
                            "path": path, "path_db": pdb, "am_dbc_hz": am, "delay_m": d,
                            "delay_ns": coax_delay_s(d) * 1e9,
                            "leak_noise_dbm_hz": float(n[i]), "kt_sys_dbm_hz": float(kt[i]),
                            "ins29_db": float(n[i] - kt[i]),
                            "floor_rise_db": float(floor_rise_db(n[i])),
                        })
        allowed = np.nanmin(np.vstack([am_allowed_dbc_hz(chars[(HEADLINE_DELAY_M, s)], BENCH_PATHS_DB["E1"])
                                       for s in seeds]), axis=0)
        h = HEADLINE_DELAY_M
        lim = []
        for i in range(n_steps):
            lim.append({
                "step": i + 1, "f_hz": float(ref["f_hz"][i]),
                "beat_min_hz": float(ref["beat_min_hz"][i]), "beat_max_hz": float(ref["beat_max_hz"][i]),
                "hpf_corner_hz": 0.5 * float(ref["beat_min_hz"][i]),
                "t_sys_k": float(ref["t_sys_k"][i]),
                "kt_sys_dbm_hz": float(kt[i]),
                "e0_noise_m140_dbm_hz": float(worst[("E0", -140.0, h)][i]),
                "e0_floor_rise_m140_db": float(floor_rise_db(worst[("E0", -140.0, h)][i])),
                "e1_noise_m140_dbm_hz": float(worst[("E1", -140.0, h)][i]),
                "e1_ins29_m140_db": float(worst[("E1", -140.0, h)][i] - kt[i]),
                "e1_floor_rise_m140_db": float(floor_rise_db(worst[("E1", -140.0, h)][i])),
                "e1_floor_rise_at_limit_db": float(floor_rise_db(kt[i])),
                "e1_am_allowed_dbc_hz": float(allowed[i]),
                "e2_noise_m140_dbm_hz": float(worst[("E2", -140.0, h)][i]),
            })
        limits[sname] = lim
        # Headline metrics for this sweep.
        for path in BENCH_PATHS_DB:
            for am in AM_CASES_DBC_HZ:
                ins29 = worst[(path, am, h)] - kt
                k = int(np.nanargmax(ins29))
                tag = f"{sname}_{path}_{int(-am)}"
                metrics[f"ins29_worst_db_{tag}"] = float(ins29[k])
                metrics[f"ins29_worst_freq_mhz_{tag}"] = float(ref["f_hz"][k] / 1e6)
                metrics[f"floor_rise_min_db_{tag}"] = float(np.nanmin(floor_rise_db(worst[(path, am, h)])))
        for am in AM_CASES_DBC_HZ:
            spread = np.nanmax(np.abs(worst[("E0", am, DELAY_SET_M[-1])] - worst[("E0", am, DELAY_SET_M[0])]))
            metrics[f"delay_sensitivity_db_{sname}_{int(-am)}"] = float(spread)
        lin_err = worst[("E1", -140.0, h)] - (worst[("E0", -140.0, h)] - (BENCH_PATHS_DB["E1"] - BENCH_PATHS_DB["E0"]))
        metrics[f"e0_to_e1_scaling_error_db_{sname}"] = float(np.nanmax(np.abs(lin_err)))
        metrics[f"e1_am_allowed_worst_dbc_hz_{sname}"] = float(np.nanmin(allowed))
        metrics[f"e1_am_allowed_worst_freq_mhz_{sname}"] = float(ref["f_hz"][int(np.nanargmin(allowed))] / 1e6)
        metrics[f"floor_rise_at_ins29_limit_min_db_{sname}"] = float(np.nanmin(floor_rise_db(kt)))
        metrics[f"kt_sys_min_dbm_hz_{sname}"] = float(np.nanmin(kt))
        metrics[f"cell_bandwidth_hz_{sname}"] = float(np.nanmedian(ref["cell_bandwidth_hz"]))
        metrics[f"beat_band_hz_{sname}"] = [float(np.nanmin(ref["beat_min_hz"])), float(np.nanmax(ref["beat_max_hz"]))]
        metrics[f"hpf_corner_hz_{sname}"] = [0.5 * float(np.nanmin(ref["beat_min_hz"])),
                                             0.5 * float(np.nanmax(ref["beat_min_hz"]))]
    fe = {p: front_end(v) for p, v in BENCH_PATHS_DB.items()}
    fe["INS-23 35 dB, no canceller"] = front_end(35.0)
    stray = stray_limits(metrics["cell_bandwidth_hz_A1"])
    am_floor = {
        "single_chain_dbc_hz": am_measurement_floor_dbc_hz(),
        "cross_spectrum_dbc_hz": am_measurement_floor_dbc_hz(averages=CROSS_SPECTRUM_AVERAGES),
        "single_chain_level17_dbc_hz": am_measurement_floor_dbc_hz(carrier_dbm=MIXER_P1DB_DBM[1] - WP9_2_MARGIN_DB),
        "carrier_dbm": MIXER_P1DB_DBM[0] - WP9_2_MARGIN_DB,
        "averages": CROSS_SPECTRUM_AVERAGES,
    }
    canceller = {
        "jump_limit_dbc": INS31_JUMP_DBC,
        "added_noise_limit_dbm_hz": INS31_ADDED_NOISE_DBM_HZ,
        # The weight's own noise, under drive, must stay inside the added-noise limit.
        "weight_noise_limit_dbc_hz_35db": INS31_ADDED_NOISE_DBM_HZ - (P_TX_DBM - 35.0),
        "residual_carrier_dbm_35_26": P_TX_DBM - 35.0 - 26.0,
    }
    return {"rows": rows, "limits": limits, "metrics": metrics, "front_end": fe, "stray": stray,
            "am_floor": am_floor, "canceller": canceller,
            "assumptions": {"p_tx_dbm": P_TX_DBM, "bench_floor_dbm_hz": BENCH_FLOOR_DBM_HZ,
                            "bench_floor_lna_dbm_hz": BENCH_FLOOR_LNA_DBM_HZ,
                            "paths_db": BENCH_PATHS_DB, "path_meaning": BENCH_PATH_MEANING,
                            "delay_set_m": DELAY_SET_M, "rg58_vf": RG58_VF, "seeds": list(seeds),
                            "ad9959_isolation_db": AD9959_ISOLATION_DB, "dds_channel_dbm": DDS_CHANNEL_DBM}}


def _summary(res: dict[str, Any]) -> str:
    m = res["metrics"]
    lines = ["WP9 conducted bench: predictions and pass limits (WP3 Phase 6).", "",
             f"Paths from the PA output at 1 W; delay {HEADLINE_DELAY_M} m RG-58 unless stated; "
             f"worst of {len(res['assumptions']['seeds'])} seeds; bench floor assumed "
             f"{BENCH_FLOOR_DBM_HZ:.0f} dBm/Hz.", "",
             "## INS-29 behind each bench path (worst step)", "",
             "| sweep | path | AM (dBc/Hz) | INS-29 (dB re k T_sys) | at (MHz) | least floor rise (dB) |",
             "|---|---|---|---|---|---|"]
    for s in ("A1", "A2"):
        for p in BENCH_PATHS_DB:
            for am in AM_CASES_DBC_HZ:
                t = f"{s}_{p}_{int(-am)}"
                lines.append(f"| {s} | {p} ({BENCH_PATHS_DB[p]:.0f} dB) | {am:.0f} | {m['ins29_worst_db_' + t]:+.1f} | "
                             f"{m['ins29_worst_freq_mhz_' + t]:.3f} | {m['floor_rise_min_db_' + t]:.2f} |")
    lines += ["", "## Checks the bench makes on the model", "",
              f"- Delay sensitivity, E0, 0.66 m to 19.8 m: A1 {m['delay_sensitivity_db_A1_140']:.2f} dB at -140, "
              f"{m['delay_sensitivity_db_A1_150']:.2f} dB at -150 (AM binds, not phase noise).",
              f"- Linear scaling E0 to E1 (22 dB): largest departure {m['e0_to_e1_scaling_error_db_A1']:.3f} dB (A1).",
              f"- AM at which E1 just meets INS-29: A1 {m['e1_am_allowed_worst_dbc_hz_A1']:.1f} dBc/Hz at "
              f"{m['e1_am_allowed_worst_freq_mhz_A1']:.3f} MHz; A2 {m['e1_am_allowed_worst_dbc_hz_A2']:.1f}.",
              f"- Floor rise at the INS-29 limit, least over the band: {m['floor_rise_at_ins29_limit_min_db_A1']:.2f} dB "
              f"(k T_sys down to {m['kt_sys_min_dbm_hz_A1']:.1f} dBm/Hz).", "",
              "## Front end", "", "| path | leakage (dBm) | level-7: margin, image | level-17: margin, image |",
              "|---|---|---|---|"]
    for k, v in res["front_end"].items():
        lines.append(f"| {k} | {v['leakage_dbm']:.1f} | {v['margin_db_p1db_+1']:.1f} dB, {v['image_dbc_p1db_+1']:.1f} dBc | "
                     f"{v['margin_db_p1db_+11']:.1f} dB, {v['image_dbc_p1db_+11']:.1f} dBc |")
    st, af, cn = res["stray"], res["am_floor"], res["canceller"]
    lines += ["", "## Stray coupling and echo levels (WP9-1)", "",
              f"- Weakest admitted echo on the bench {st['weakest_echo_dbm']:.1f} dBm; strong echo {st['strong_echo_dbm']:.1f} dBm.",
              f"- Transmit-chain stray: >= {st['tx_stray_carrier_db']:.0f} dB from the PA output (10 dB under E2); "
              f"its noise needs >= {st['tx_stray_noise_db']:.0f} dB.",
              f"- A synthesiser-channel echo (AD9959, > {AD9959_ISOLATION_DB:.0f} dB channel isolation): weakest usable "
              f"echo {st['dds_echo_min_dbm_E0']:.0f} / {st['dds_echo_min_dbm_E1']:.0f} / {st['dds_echo_min_dbm_E2']:.0f} dBm "
              f"behind E0 / E1 / E2, i.e. SNR {st['dds_echo_min_snr_db_E0']:.0f} / {st['dds_echo_min_snr_db_E1']:.0f} / "
              f"{st['dds_echo_min_snr_db_E2']:.0f} dB on the bench.  Weak echoes with leakage present need the DRFM.", "",
              "## AM measurement (WP9-5) and canceller (WP9-8)", "",
              f"- AM floor, one chain, level-7 mixer at {af['carrier_dbm']:.0f} dBm: {af['single_chain_dbc_hz']:.1f} dBc/Hz; "
              f"two-chain cross-spectrum, {af['averages']} averages: {af['cross_spectrum_dbc_hz']:.1f}; level-17 mixer: "
              f"{af['single_chain_level17_dbc_hz']:.1f}.",
              f"- INS-31 weight step limit {cn['jump_limit_dbc']:.1f} dB of the leakage; added noise "
              f"<= {cn['added_noise_limit_dbm_hz']:.0f} dBm/Hz, which bounds the weight's own noise under drive at "
              f"{cn['weight_noise_limit_dbc_hz_35db']:.0f} dBc/Hz of the leakage at 35 dB isolation.", ""]
    return "\n".join(lines)


@register_study(
    "wp9_bench",
    description="WP3 Phase 6: WP9 conducted-bench predictions and pass limits (INS-29 behind emulated canceller paths)",
    work_packages=["WP3", "WP9"],
    claims=["C4"],
    status="implemented",
    requirements=["INS-11", "INS-24", "INS-29", "INS-30", "INS-31"],
)
def run(out_dir: Path, *, config: dict[str, Any], manifest: Any) -> StudyResult:
    res = compute()
    rows = res["rows"]
    cols = list(rows[0].keys())
    pred = out_dir / "wp9_bench_predictions.csv"
    with pred.open("w", encoding="utf-8") as fh:
        fh.write(",".join(cols) + "\n")
        for r in rows:
            fh.write(",".join(str(r[c]) for c in cols) + "\n")
    outs = [pred]
    for sname, lim in res["limits"].items():
        p = out_dir / f"wp9_bench_limits_{sname}.csv"
        lc = list(lim[0].keys())
        with p.open("w", encoding="utf-8") as fh:
            fh.write(",".join(lc) + "\n")
            for r in lim:
                fh.write(",".join(f"{r[c]:.6g}" if isinstance(r[c], float) else str(r[c]) for c in lc) + "\n")
        outs.append(p)
    j = out_dir / "wp9_bench.json"
    j.write_text(json.dumps({k: v for k, v in res.items() if k not in ("rows", "limits")}, indent=2, default=str),
                 encoding="utf-8")
    outs.append(j)
    flat = {k: v for k, v in res["metrics"].items() if not isinstance(v, list)}
    return StudyResult(outputs=outs, summary=_summary(res), metrics=flat)
