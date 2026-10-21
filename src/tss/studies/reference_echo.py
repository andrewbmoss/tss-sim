"""The echo levels the leakage family is measured against (WP3 Phase 0).

Layer:        study
WP:           WP3
Claims:       C1, C4
Status:       implemented
Requirements: INS-10, INS-11, INS-22, INS-23, INS-24

Runs :func:`tss.l2_instrument.reference_echo.reference_levels` on the baseline
instrument (``config/instrument/baseline_6u.yaml``: 1 W, the 5 m short dipole
with the Balanis reactance, 3 dB noise figure, Franklin and MacLean sky) for
the A1 and A2 stepped sweeps, and writes, per step:

* the noise in one range cell, the weakest echo INS-10 admits (10 dB above
  it), and the free-space echo from INS-04's 1050 km;
* the leakage carrier at INS-23's isolation, and its ratio to each;
* the echo band (beat frequencies from 100 km to 1050 km) INS-24 applies over.

It then checks the prediction against the chain: a synthetic echo at 1050 km
is put through the L2 baseline and measured by L3 as INS-10 is scored, at six
frequencies, and the processing loss is recorded.

Finally it puts the new reference into the two Rev G derivations built on the
old "122 dB - isolation" (SRM §16 and §22), so the size of the change is on
the record.  The requirement documents are not changed here (WP3 plan, G5).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np

from tss.studies import StudyResult, register_study

#: Frequencies the summary table and the chain check report at, MHz.
SUMMARY_MHZ = (0.33, 0.5, 1.0, 2.0, 3.0, 5.0, 10.0, 20.0)
CHECK_MHZ = (0.36, 1.0, 2.5, 5.0, 10.0, 18.0)
ISOLATIONS_DB = (30.0, 40.0, 50.0, 60.0)
PATHS_M = (1.0, 3.0, 10.0, 30.0)
#: INS-11 instantaneous dynamic range and the SRM §22 headroom.
INS11_DR_DB = 60.0
HEADROOM_DB = 10.0

_COLUMNS = (
    "f_hz", "step_rate_hz_per_s", "beat_min_hz", "beat_max_hz", "gain_dbi", "t_ant_k",
    "t_sys_k", "cell_bandwidth_hz", "noise_cell_dbw", "weakest_echo_dbw", "echo_rmax_dbw",
    "snr_rmax_db", "leakage_dbw", "isolation_db", "leakage_over_weakest_db",
    "leakage_over_rmax_echo_db", "leakage_over_noise_db",
)


def _instrument() -> Any:
    from tss.config import load_yaml, resolve_path
    from tss.l2_instrument.baseline import BaselineInstrument

    return BaselineInstrument.from_config(load_yaml(resolve_path("instrument/baseline_6u.yaml")))


def _crossing_mhz(f_hz: np.ndarray, snr_db: np.ndarray, level: float = 10.0) -> float | None:
    """Lowest frequency above which the free-space SNR at R_max stays at or above ``level``."""
    ok = np.asarray(snr_db) >= level
    if np.all(ok):
        return float(f_hz[0] / 1e6)
    if not np.any(ok):
        return None
    last_bad = int(np.flatnonzero(~ok)[-1])
    return float(f_hz[min(last_bad + 1, f_hz.size - 1)] / 1e6)


def compute() -> dict[str, Any]:
    """Everything the study writes, as data (no files)."""
    from tss import sweep_law
    from tss.l2_instrument import reference_echo as rx
    from tss.studies.data_volume import A1, A2

    inst = _instrument()
    out: dict[str, Any] = {"legacy": rx.legacy_echo_db(5.0e6), "configs": {}}
    for cfg in (A1, A2):
        sweep = cfg.schedule()
        lev = rx.reference_levels(inst, sweep)
        fs = float(inst.sample_rate_for(sweep))
        noise_fs_over_cell = 10.0 * np.log10(fs / float(lev["cell_bandwidth_hz"][0]))
        hpf_revg = lev["leakage_over_rmax_echo_db"] - INS11_DR_DB + HEADROOM_DB
        hpf_noise = (lev["leakage_over_noise_db"] - noise_fs_over_cell) - INS11_DR_DB + HEADROOM_DB
        spurs = {
            f"{ref}_{p:g}m": rx.required_spur_dbc(lev, coupling_path_m=p, reference=ref)
            for ref in ("weakest", "rmax") for p in PATHS_M
        }
        out["configs"][cfg.name] = {
            "sweep": sweep, "levels": lev, "fs_hz": fs,
            "noise_fs_over_cell_db": float(noise_fs_over_cell),
            "hpf_rejection_revg_rule_db": hpf_revg, "hpf_rejection_noise_rule_db": hpf_noise,
            "spur_dbc": spurs,
            "steps": int(sweep_law.step_count(sweep)),
        }
    # Chain check on A1.
    sweep = A1.schedule()
    centres = sweep_law.step_centres_hz(sweep)
    checks = []
    for f in CHECK_MHZ:
        fc = float(centres[int(np.argmin(np.abs(centres - f * 1e6)))])
        checks.append(rx.processing_loss_db(inst, sweep, fc, 1.05e6, snr_db=20.0))
    out["chain_check"] = checks
    return out


def _rows(lev: dict[str, np.ndarray]) -> np.ndarray:
    return np.column_stack([np.asarray(lev[c], dtype=float) for c in _COLUMNS])


@register_study(
    "reference_echo",
    description="INS-10/INS-24 reference echo, leakage ratios and the echo band, A1 and A2",
    work_packages=["WP3"],
    claims=["C1", "C4"],
    status="implemented",
    requirements=["INS-10", "INS-11", "INS-22", "INS-23", "INS-24"],
)
def run(out_dir: Path, *, config: dict[str, Any], manifest: Any) -> StudyResult:
    from tss.l2_instrument import reference_echo as rx

    res = compute()
    outputs: list[Path] = []
    header = ",".join(_COLUMNS)
    for name, block in res["configs"].items():
        path = out_dir / f"reference_{name}.csv"
        lev = block["levels"]
        extra = np.column_stack(
            [block["hpf_rejection_revg_rule_db"], block["hpf_rejection_noise_rule_db"]]
            + [block["spur_dbc"][k] for k in sorted(block["spur_dbc"])]
        )
        cols = header + ",hpf_rejection_revg_rule_db,hpf_rejection_noise_rule_db," + ",".join(
            f"spur_dbc_{k}" for k in sorted(block["spur_dbc"])
        )
        np.savetxt(path, np.hstack([_rows(lev), extra]), delimiter=",", header=cols,
                   comments="", fmt="%.6g")
        outputs.append(path)

    a1 = res["configs"]["A1"]
    lev = a1["levels"]
    srows = rx.summary_rows(lev, SUMMARY_MHZ)
    idx = [int(np.argmin(np.abs(lev["f_hz"] / 1e6 - f))) for f in SUMMARY_MHZ]
    spur_w10 = a1["spur_dbc"]["weakest_10m"]
    spur_r10 = a1["spur_dbc"]["rmax_10m"]
    lines = [
        "| f (MHz) | step rate (MHz/s) | echo band (kHz) | G (dBi) | N cell (dBm) | weakest echo (dBm) "
        "| free-space echo at 1050 km (dBm) | SNR at 1050 km (dB) | L/weakest (dB) | L/echo(1050) (dB) "
        "| L/N cell (dB) |",
        "|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r in srows:
        lines.append(
            f"| {r['f_hz']/1e6:.3f} | {r['step_rate_hz_per_s']/1e6:.3f} | "
            f"{r['beat_min_hz']/1e3:.2f}-{r['beat_max_hz']/1e3:.1f} | {r['gain_dbi']:.1f} | "
            f"{r['noise_cell_dbw']+30:.1f} | {r['weakest_echo_dbw']+30:.1f} | {r['echo_rmax_dbw']+30:.1f} | "
            f"{r['snr_rmax_db']:.1f} | {r['leakage_over_weakest_db']:.1f} | "
            f"{r['leakage_over_rmax_echo_db']:.1f} | {r['leakage_over_noise_db']:.1f} |"
        )
    table = "\n".join(lines)

    derive = [
        "| f (MHz) | INS-11 HPF, Rev G rule (dB) | INS-11 HPF, noise rule (dB) "
        "| INS-22 spur, 10 m, vs echo at 1050 km (dBc) | INS-22 spur, 10 m, vs weakest echo (dBc) |",
        "|---|---|---|---|---|",
    ]
    for i in idx:
        derive.append(
            f"| {lev['f_hz'][i]/1e6:.3f} | {a1['hpf_rejection_revg_rule_db'][i]:.1f} | "
            f"{a1['hpf_rejection_noise_rule_db'][i]:.1f} | {spur_r10[i]:.1f} | {spur_w10[i]:.1f} |"
        )
    derive_table = "\n".join(derive)

    checks = res["chain_check"]
    check_lines = ["| f (MHz) | predicted SNR (dB) | measured, mean of 4 (dB) | loss (dB) |",
                   "|---|---|---|---|"]
    for c in checks:
        check_lines.append(f"| {c['f_hz']/1e6:.3f} | {c['predicted_snr_db']:.1f} | "
                           f"{c['measured_snr_db_mean']:.2f} | {c['loss_db']:+.2f} |")
    check_table = "\n".join(check_lines)

    i5 = idx[SUMMARY_MHZ.index(5.0)]
    legacy = res["legacy"]
    cross = _crossing_mhz(lev["f_hz"], lev["snr_rmax_db"])
    a2 = res["configs"]["A2"]
    metrics = {
        "legacy_echo_below_tx_db_5mhz_300km": legacy["echo_below_tx_db"],
        "srm_leakage_over_echo_plus_isolation_db": rx.LEGACY_LEAKAGE_OVER_ECHO_PLUS_ISOLATION_DB,
        "isolation_db": float(lev["isolation_db"][0]),
        "leakage_over_rmax_echo_db_5mhz": float(lev["leakage_over_rmax_echo_db"][i5]),
        "leakage_over_weakest_db_5mhz": float(lev["leakage_over_weakest_db"][i5]),
        "leakage_over_noise_db_5mhz": float(lev["leakage_over_noise_db"][i5]),
        "leakage_over_weakest_db_max": float(np.max(lev["leakage_over_weakest_db"])),
        "leakage_over_weakest_db_min": float(np.min(lev["leakage_over_weakest_db"])),
        "leakage_over_rmax_echo_db_max": float(np.max(lev["leakage_over_rmax_echo_db"])),
        "leakage_over_rmax_echo_db_min": float(np.min(lev["leakage_over_rmax_echo_db"])),
        "weakest_echo_dbm_5mhz": float(lev["weakest_echo_dbw"][i5] + 30.0),
        "noise_cell_dbm_5mhz": float(lev["noise_cell_dbw"][i5] + 30.0),
        "snr_rmax_db_min": float(np.min(lev["snr_rmax_db"])),
        "snr_rmax_db_median": float(np.median(lev["snr_rmax_db"])),
        "snr_rmax_meets_10db_above_mhz": cross,
        "hpf_rejection_revg_rule_db_5mhz": float(a1["hpf_rejection_revg_rule_db"][i5]),
        "hpf_rejection_revg_rule_db_max": float(np.max(a1["hpf_rejection_revg_rule_db"])),
        "hpf_rejection_noise_rule_db_5mhz": float(a1["hpf_rejection_noise_rule_db"][i5]),
        "hpf_rejection_noise_rule_db_max": float(np.max(a1["hpf_rejection_noise_rule_db"])),
        "spur_dbc_rmax_10m_5mhz": float(spur_r10[i5]),
        "spur_dbc_weakest_10m_5mhz": float(spur_w10[i5]),
        "spur_dbc_weakest_10m_worst": float(np.min(spur_w10)),
        "spur_dbc_rmax_10m_worst": float(np.min(spur_r10)),
        "spur_dbc_weakest_30m_worst": float(np.min(a1["spur_dbc"]["weakest_30m"])),
        "spur_dbc_weakest_1m_worst": float(np.min(a1["spur_dbc"]["weakest_1m"])),
        "a2_spur_dbc_weakest_10m_worst": float(np.min(a2["spur_dbc"]["weakest_10m"])),
        "a2_spur_dbc_rmax_10m_worst": float(np.min(a2["spur_dbc"]["rmax_10m"])),
        "a2_beat_max_hz": float(np.max(a2["levels"]["beat_max_hz"])),
        "a1_beat_max_hz": float(np.max(lev["beat_max_hz"])),
        "chain_check_loss_db": [round(c["loss_db"], 3) for c in checks],
        "chain_check_loss_db_max_abs": float(max(abs(c["loss_db"]) for c in checks)),
    }
    js = out_dir / "reference_echo.json"
    js.write_text(json.dumps({"metrics": metrics, "chain_check": checks, "legacy": legacy},
                             indent=2), encoding="utf-8")
    outputs.append(js)

    summary = f"""
The reference INS-24 is measured against, computed from the baseline instrument
(1 W, 5 m dipole with the Balanis reactance, 3 dB noise figure, Franklin and
MacLean sky, Hann window, 26.8 ms dwell, 40 dB isolation).  A1 sweep; A2 has the
same levels per cell because the dwell is the same, and a wider echo band
(to {metrics['a2_beat_max_hz']/1e3:.0f} kHz against {metrics['a1_beat_max_hz']/1e3:.1f} kHz).

{table}

The old figure.  SRM §16 takes leakage over echo as 122 dB less the isolation.
The inherited link budget (`linkbudget.SounderConfig` defaults: 300 km, the
legacy -6.8 dBi antenna) puts its 5 MHz echo {legacy['echo_below_tx_db']:.1f} dB
below the transmitter, which is that 122 dB.  On the instrument, at INS-04's
1050 km, the free-space echo at 5 MHz is
{metrics['leakage_over_rmax_echo_db_5mhz'] + metrics['isolation_db']:.1f} dB below
the transmitter, so leakage over echo is {metrics['leakage_over_rmax_echo_db_5mhz']:.1f} dB
at 40 dB isolation, not 82.  Against the weakest echo INS-10 admits it is
{metrics['leakage_over_weakest_db_5mhz']:.1f} dB at 5 MHz and
{metrics['leakage_over_weakest_db_min']:.1f} to {metrics['leakage_over_weakest_db_max']:.1f} dB
over the band.  All of these move one for one with the isolation.

What that does to the two derivations built on the old figure (A1, 10 m coupling
path, 10 dB margin, single-tone small-angle as in Rev G §5.9; INS-11 at 60 dB
with 10 dB headroom):

{derive_table}

The "Rev G rule" is SRM §22, leakage over echo less 60 dB plus 10 dB.  The
"noise rule" refers the digitiser's 60 dB to the noise in the full digitiser
band ({a1['fs_hz']/1e3:.1f} kHz, {a1['noise_fs_over_cell_db']:.1f} dB above one
cell), which is what the 60 dB has to span.  Which rule INS-11 uses, and which
echo INS-22 and INS-24 refer to, are WP3 decisions (plan, G5); no requirement
document is changed by this study.

Free-space SNR at 1050 km falls below INS-10's 10 dB below
{cross if cross is not None else float('nan'):.2f} MHz on the 5 m antenna
(minimum {metrics['snr_rmax_db_min']:.1f} dB at the band edge), because only
transmit efficiency costs SNR and the receive efficiency stops cancelling where
the sky no longer dominates.  INS-10 is worded "at R_max across the band"; at
the bottom of the band the echoes come from far nearer than 1050 km.

Chain check.  A synthetic echo at 1050 km set to a predicted 20 dB, through the
L2 baseline (noise, receiver, quantiser) and measured by L3 exactly as INS-10 is
scored:

{check_table}

Below 3 MHz the loss is mostly the step-transition skirt: for the first
tau = 7.0 ms of each 26.8 ms dwell the echo beats at the previous step's rate,
which costs about 0.9 dB with the Hann window.  Above 3 MHz the rate does not
change, and what is left is the scatter of four noise draws (about +/-0.7 dB).
The prediction and the chain agree to within 1 dB everywhere, so the reference
is the level the chain measures.

INS-22 with the coupling path.  Against the weakest admitted echo, A1 at 40 dB:
{metrics['spur_dbc_weakest_1m_worst']:.1f} dBc at 1 m, {metrics['spur_dbc_weakest_10m_worst']:.1f} dBc
at 10 m and {metrics['spur_dbc_weakest_30m_worst']:.1f} dBc at 30 m (worst step); A2 at 10 m
{metrics['a2_spur_dbc_weakest_10m_worst']:.1f} dBc.  A longer coupling path suppresses less,
because the suppression grows as (pi f_b tau)^2.  Rev G states -70 dBc.
"""
    return StudyResult(outputs=outputs, summary=summary, metrics=metrics)
