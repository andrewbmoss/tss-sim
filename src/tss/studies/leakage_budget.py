"""Layer: studies
WP: WP3
Claims: C4
Status: implemented
Requirements: INS-10, INS-11, INS-14, INS-22, INS-23, INS-24

WP3 Phase 5: the leakage budget.

Phase 4 found that the transmitter's AM noise, carried on the leakage, sets the
residual in the echo band, and that INS-24 fails by 27 dB at the bottom of the
band under the threshold allocations.  This study trades what can close it:

* the transmitter's AM noise density (dBc/Hz, single sideband);
* an RF canceller ahead of the front end, fed by a sample of the PA output
  (:class:`tss.l2_instrument.leakage.RFCanceller`; decision of 2026-10-28),
  with amplitude, phase and delay errors, three retuning policies, and passive
  or active (added noise);
* the isolation (the moment-method bus profile and flat brackets);
* the transmit power below INS-14's 1 W (decision of 2026-10-28), whose INS-24
  gain is paid for in INS-10;
* the front end's 1 dB compression point (level-7, -13 and -17 mixers);
* the start-up: primed (the operating assumption), cold, and one settling step.

METHOD
------
The residual in a step's echo band is the sum of terms that scale with the
leakage power and terms that do not:

    ped_i = P_leak,i [ X_i 10^((AM + 133)/10) + Y_i ] + N_c B_cell

``X_i`` is the AM pedestal per unit leakage at the -133 dBc/Hz reference after
whatever the canceller does to it; ``Y_i`` is everything else per unit leakage
(the range-correlated phase noise, the spurs, and the canceller's own
dynamics: the carrier it leaves changes whenever its weight moves, and the
high-pass passes that change into the echo band); ``N_c`` is an active
canceller's noise.  ``X_i`` and ``Y_i`` are measured by simulation, the full L2
chain with the AM on and off, per sweep, isolation shape, canceller and seed.
The grid over AM, power and isolation level then follows by scaling, and a
sample of grid points is re-simulated directly to check it.

INS-10's cost is taken from the SNR that L3 measured at every truth echo of
the smoke campaign on an ideal instrument (leakage off), lowered by the power
reduction, by ``10 log10(1 + ped/N)`` in the echo's step, and by the insertion
loss of the canceller's coupler where one is fitted.  The smoke campaign
itself is re-run at sample points to check that too.
"""

from __future__ import annotations

import copy
import itertools
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from tss.studies import StudyResult, register_study

INS24_MARGIN_DB = 10.0
INS10_SNR_DB = 10.0
WP9_2_FE_MARGIN_DB = 10.0
AM_REF_DBC_HZ = -133.0
AM_LEVELS_DBC_HZ = (-133.0, -136.0, -140.0, -145.0, -150.0, -155.0, -160.0, -165.0)
#: The second AM level of the characterisation (see :func:`characterise`).
AM_LOW_DBC_HZ = -153.0
POWERS_W = (1.0, 0.5, 0.25, 0.1)
#: Front-end input 1 dB compression points: level-7, level-13 and level-17 diode mixers.
FRONT_END_P1DB_DBM = (1.0, 7.0, 11.0)
#: Isolations: the moment-method bus profile, and flat brackets.  20 dB stands
#: for the bus mount with about 1 degree of deployment error (14 to 21 dB).
ISOLATIONS = ("mom_bus", "flat20", "flat30", "flat40", "flat50", "flat60")
#: Insertion loss the canceller's injection coupler adds to the receive path, dB.
COUPLER_LOSS_DB = 0.5
#: An active vector modulator's output noise (-148 dBm/Hz) behind a 10 dB coupler.
ACTIVE_NOISE_DBM_HZ = -158.0
SEEDS = (0, 1, 2)

_GOOD = {"amplitude_error_db": 0.1, "phase_error_deg": 1.0, "delay_mismatch_s": 1.0e-9}
#: Canceller cases: name -> RFCanceller keywords (None for no canceller).
CANCELLERS: dict[str, dict[str, Any] | None] = {
    "none": None,
    "ideal20": {"ideal_depth_db": 20.0},
    "ideal26_ins31": {"ideal_depth_db": 26.0},
    "ideal30": {"ideal_depth_db": 30.0},
    "ideal40": {"ideal_depth_db": 40.0},
    "rf_good": {**_GOOD, "retune": "track", "error_correlation_steps": 10.0},
    "rf_moderate": {"amplitude_error_db": 0.3, "phase_error_deg": 2.0, "delay_mismatch_s": 1.0e-9,
                    "retune": "track", "error_correlation_steps": 10.0},
    "rf_poor": {"amplitude_error_db": 0.5, "phase_error_deg": 3.0, "delay_mismatch_s": 1.0e-9,
                "retune": "track", "error_correlation_steps": 10.0},
    "rf_good_bias_only": {**_GOOD, "error_draw": "worst"},
    "rf_good_uncorrelated": {**_GOOD, "retune": "track", "error_correlation_steps": 0.0},
    "rf_good_stepped": {**_GOOD, "retune": "step", "error_correlation_steps": 10.0},
    "rf_good_stepped_uncorrelated": {**_GOOD, "retune": "step", "error_correlation_steps": 0.0},
    "rf_one_setting": {**_GOOD, "retune": "sweep", "calibration_hz": 1.0e6, "error_draw": "worst"},
}
#: How large a jump of the weight at a step boundary can be: an otherwise
#: perfect canceller switched at each boundary by an independent amplitude
#: error of this RMS (a jump of about 1.4 times it).
JUMP_RMS_DB = (0.001, 0.003, 0.01, 0.03)
for _j in JUMP_RMS_DB:
    CANCELLERS[f"jump_{_j:g}dB"] = {"amplitude_error_db": _j, "phase_error_deg": 0.0, "retune": "step",
                                    "error_correlation_steps": 0.0}
del _j
#: Cases computed from another by adding the active modulator's noise.
ACTIVE_FROM = {"rf_good_active": "rf_good"}


# ----------------------------------------------------------------------------
# Building the instrument
# ----------------------------------------------------------------------------


def _base_cfg() -> dict[str, Any]:
    from tss.config import load_yaml, resolve_path

    cfg = copy.deepcopy(load_yaml(resolve_path("instrument/baseline_6u.yaml")))
    cfg["leakage"]["enabled"] = True
    cfg["leakage"]["primed"] = True
    return cfg


def _cfg(*, profile: str = "mom_bus", am: float | None = AM_REF_DBC_HZ, canceller: str = "none",
         power_w: float = 1.0, flat_db: float = 40.0) -> dict[str, Any]:
    cfg = _base_cfg()
    cfg["transmit"]["power_w"] = float(power_w)
    cfg["transmit_impairments"]["am_noise_dbc_hz"] = am
    if profile == "flat":
        cfg["leakage"]["profile"] = "flat"
        cfg["isolation"]["transmit_receive_db"] = float(flat_db)
    spec = CANCELLERS.get(canceller)
    if canceller in ACTIVE_FROM:
        spec = {**CANCELLERS[ACTIVE_FROM[canceller]], "added_noise_dbm_hz": ACTIVE_NOISE_DBM_HZ}
    cfg["canceller"] = {"enabled": spec is not None, **(spec or {})}
    return cfg


def _no_echoes() -> Any:
    from tss.studies.leakage_residual import _no_echoes as ne

    return ne()


def _run(cfg: dict[str, Any], sweep: Any, seed: int) -> tuple[dict[str, np.ndarray], dict[str, Any], Any]:
    from tss.l2_instrument.baseline import BaselineInstrument

    inst = BaselineInstrument.from_config(cfg)
    rec = inst.receive(_no_echoes(), sweep, np.random.default_rng(seed))
    ps = {k: np.asarray(v, dtype=float) for k, v in rec.flags["leakage_per_step"].items()}
    return ps, dict(rec.flags["l2_metrics"]), inst


def _lin(db: np.ndarray) -> np.ndarray:
    return 10.0 ** (np.asarray(db, dtype=float) / 10.0)


def _db(x: np.ndarray) -> np.ndarray:
    return 10.0 * np.log10(np.maximum(np.asarray(x, dtype=float), 1e-300))


# ----------------------------------------------------------------------------
# Characterising the residual by simulation
# ----------------------------------------------------------------------------


def characterise(sweep: Any, profile: str, canceller: str, seed: int) -> dict[str, np.ndarray]:
    """``X``, ``Y`` and the references for one sweep, isolation shape, canceller and seed.

    The leakage at 1 W on the profile's own levels (40 dB for the flat shape).
    ``X`` and ``Y`` are per unit leakage power: watts in a cell per watt of
    leakage at the receiver input.
    """
    from tss.l2_instrument.reference_echo import reference_levels

    # Two AM levels rather than AM on and off: the AM draw is then the same
    # draw scaled, and every later draw from the record's generator (the
    # canceller's weight errors among them) is the same draw too, so the two
    # records differ only in the AM.
    on, _, inst = _run(_cfg(profile=profile, canceller=canceller), sweep, seed)
    low, _, _ = _run(_cfg(profile=profile, canceller=canceller, am=AM_LOW_DBC_HZ), sweep, seed)
    p_leak = _lin(-on["isolation_db"]) * float(inst.tx_power_w) * 10.0 ** (inst.tx_chain_efficiency_db / 10.0)
    level_on, level_low = _lin(on["residual_level_dbw"]), _lin(low["residual_level_dbw"])
    r = 10.0 ** ((AM_LOW_DBC_HZ - AM_REF_DBC_HZ) / 10.0)
    x = np.maximum(level_on - level_low, 0.0) / (1.0 - r)
    lev = reference_levels(inst, sweep)
    return {
        "f_hz": on["f_hz"],
        "isolation_db": on["isolation_db"],
        "X": x / p_leak,
        "Y": np.maximum(level_on - x, 0.0) / p_leak,
        "weakest_echo_dbw": on["weakest_echo_dbw"],
        "noise_cell_dbw": np.asarray(lev["noise_cell_dbw"], dtype=float),
        "cell_bandwidth_hz": np.asarray(lev["cell_bandwidth_hz"], dtype=float),
        "level_on_dbw": on["residual_level_dbw"],
    }


def canceller_carrier_factor(sweep: Any, profile: str, canceller: str, seed: int, fs_hz: float) -> np.ndarray:
    """Mean ``|g|^2`` per step: the carrier the canceller leaves, per unit leakage.

    It sets the front end's drive.  Drawn from its own generator, so it is the
    same statistics as the simulated record, not the same draw.
    """
    from tss import sweep_law
    from tss.l2_instrument.baseline import BaselineInstrument

    n_steps = sweep_law.step_count(sweep)
    inst = BaselineInstrument.from_config(_cfg(profile=profile, canceller=canceller))
    if inst.canceller is None:
        return np.ones(n_steps)
    dwell = float(sweep.step_dwell_s)
    t = np.arange(int(n_steps * dwell * 2000.0)) / 2000.0
    g = inst.canceller.residual_factor(sweep, t, inst.leakage, np.random.default_rng(seed))
    idx = np.clip((t / dwell).astype(int), 0, n_steps - 1)
    p = np.bincount(idx, weights=np.abs(g) ** 2, minlength=n_steps) / np.maximum(np.bincount(idx, minlength=n_steps), 1)
    return p


# ----------------------------------------------------------------------------
# The grid
# ----------------------------------------------------------------------------


def _iso_case(name: str) -> tuple[str, float]:
    return ("mom_bus", 0.0) if name == "mom_bus" else ("flat", float(name[4:]))


def predict(ch: dict[str, np.ndarray], *, iso_name: str, am: float, power_w: float,
            added_noise_dbm_hz: float | None = None) -> dict[str, np.ndarray]:
    """Per-step residual and INS-24 margin at one grid point, from a characterisation."""
    profile, flat_db = _iso_case(iso_name)
    iso = ch["isolation_db"] if profile == "mom_bus" else np.full(ch["isolation_db"].shape, flat_db)
    p_leak = float(power_w) * _lin(-iso)
    ped = p_leak * (ch["X"] * 10.0 ** ((am - AM_REF_DBC_HZ) / 10.0) + ch["Y"])
    if added_noise_dbm_hz is not None:
        ped = ped + 10.0 ** ((added_noise_dbm_hz - 30.0) / 10.0) * ch["cell_bandwidth_hz"]
    level = _db(ped)
    return {"level_dbw": level, "margin_db": ch["weakest_echo_dbw"] - level,
            "ped_over_noise_db": level - ch["noise_cell_dbw"], "p_leak_w": p_leak}


def am_needed_dbc_hz(ch: dict[str, np.ndarray], *, iso_name: str, power_w: float,
                     added_noise_dbm_hz: float | None = None) -> float | None:
    """The AM density at which the worst step just meets INS-24, under the model.

    ``None`` when no AM level would: the terms that do not scale with AM
    already exceed the allowance at some step.
    """
    profile, flat_db = _iso_case(iso_name)
    iso = ch["isolation_db"] if profile == "mom_bus" else np.full(ch["isolation_db"].shape, flat_db)
    p_leak = float(power_w) * _lin(-iso)
    allow = _lin(ch["weakest_echo_dbw"] - INS24_MARGIN_DB) - p_leak * ch["Y"]
    if added_noise_dbm_hz is not None:
        allow = allow - 10.0 ** ((added_noise_dbm_hz - 30.0) / 10.0) * ch["cell_bandwidth_hz"]
    ok = np.isfinite(allow) & np.isfinite(ch["X"])     # steps the record scored
    if np.any(allow[ok] <= 0.0):
        return None
    need = AM_REF_DBC_HZ + _db(allow[ok] / np.maximum(p_leak[ok] * ch["X"][ok], 1e-300))
    return float(np.min(need))


def front_end_metrics(p_fe_w: np.ndarray, p1db_dbm: float) -> dict[str, float]:
    from tss.l2_instrument import receive_chain as rc

    fe = rc.FrontEnd(input_p1db_dbm=float(p1db_dbm))
    worst = float(np.max(p_fe_w))
    return {
        "fe_margin_db": float(p1db_dbm) - (10.0 * math.log10(worst) + 30.0),
        "fe_image_dbc": float(rc.compression_image_dbc(fe, worst)),
        "fe_desens_db": float(rc.small_signal_gain_db(fe, worst)),
    }


# ----------------------------------------------------------------------------
# INS-10: truth echoes
# ----------------------------------------------------------------------------


def smoke_truth_echoes(instrument_cfg: dict[str, Any] | None = None, campaign: str = "smoke") -> dict[str, Any]:
    """Run the campaign's first configuration once and return the truth-echo SNR.

    ``instrument_cfg`` replaces the instrument file (``None``: an ideal one,
    leakage and transmitter impairments off).  Returns the frequencies, ranges
    and L3-measured SNR of every truth echo, and the INS-10 minimum.
    """
    from tss.config import load_campaign, load_yaml, resolve_path
    from tss.l2_instrument.baseline import BaselineInstrument
    from tss.pipeline import Pipeline, build_tiers
    from tss.studies import _grids, _platform_for, _sweep_from_config, instrument_config

    loaded = load_campaign(campaign, ())
    config = loaded.config
    shared = build_tiers({k: v for k, v in loaded.tiers.items() if k != "l2"})
    import io

    sweep, _ = _sweep_from_config(config, io.StringIO())
    alt_m, freq_hz = _grids(config)
    epoch = str(config.get("epoch", "2026-01-01T12:00:00Z"))
    name = loaded.configurations[0]
    configuration = load_yaml(resolve_path(name, subdir="configurations"))
    icfg = instrument_config(config, configuration)
    if instrument_cfg is None:
        icfg["leakage"]["enabled"] = False
        icfg["transmit_impairments"] = {}
    else:
        # The configuration's own overrides (antenna length, channels) stay;
        # the blocks this study trades come from ``instrument_cfg``.
        for key in ("transmit", "transmit_impairments", "leakage", "isolation", "canceller",
                    "front_end", "leakage_highpass", "oscillator"):
            if key in instrument_cfg:
                icfg[key] = copy.deepcopy(instrument_cfg[key])
    pipe = Pipeline({**shared, "l2": BaselineInstrument.from_config(icfg)}, seed=loaded.seed,
                    keep_intermediates=False)
    res = pipe.run_configuration(name, platform=_platform_for(configuration, epoch), sweep=sweep,
                                 freq_hz=freq_hz, alt_m=alt_m, epoch=epoch)
    tes = res.truth_echo_snr or {}
    snr = np.asarray(tes.get("snr_db", []), dtype=float)
    return {
        "freq_hz": np.asarray(tes.get("freq_hz", []), dtype=float),
        "range_m": np.asarray(tes.get("range_m", []), dtype=float),
        "snr_db": snr,
        "ins10_min_db": float(np.nanmin(snr)) if np.any(np.isfinite(snr)) else float("nan"),
        "assessment": res.assessment,
        "sweep": sweep,
    }


def ins10_predicted(truth: dict[str, Any], f_steps: np.ndarray, ped_over_noise_db: np.ndarray, *,
                    power_w: float, coupler_loss_db: float = 0.0) -> dict[str, float]:
    """INS-10 at a grid point from the ideal instrument's truth-echo SNR."""
    f = truth["freq_hz"]
    snr0 = truth["snr_db"]
    ok = np.isfinite(snr0)
    j = np.clip(np.searchsorted(f_steps, f[ok]) , 0, f_steps.size - 1)
    j = np.where((j > 0) & (np.abs(f_steps[np.maximum(j - 1, 0)] - f[ok]) < np.abs(f_steps[j] - f[ok])), j - 1, j)
    loss = 10.0 * np.log10(1.0 + _lin(ped_over_noise_db[j]))
    snr = snr0[ok] + 10.0 * math.log10(float(power_w)) - loss - float(coupler_loss_db)
    return {
        "ins10_min_db": float(np.min(snr)) if snr.size else float("nan"),
        "ins10_echoes_below": int(np.count_nonzero(snr < INS10_SNR_DB)),
        "ins10_echoes": int(snr.size),
        "ins10_median_db": float(np.median(snr)) if snr.size else float("nan"),
    }


# ----------------------------------------------------------------------------
# Start-up
# ----------------------------------------------------------------------------


def start_up(sweep: Any, canceller: str = "none", seed: int = 0) -> dict[str, dict[str, float]]:
    """Primed, cold, and one settling step, on the default (MoM bus, -133 dBc/Hz AM)."""
    from tss import sweep_law
    from tss.l2_instrument import leakage as lk
    from tss.l2_instrument.baseline import BaselineInstrument
    from tss.l2_instrument.reference_echo import reference_levels

    out: dict[str, dict[str, float]] = {}
    for case in ("primed", "cold", "settling_step"):
        cfg = _cfg(canceller=canceller)
        cfg["leakage"]["primed"] = case == "primed"
        inst = BaselineInstrument.from_config(cfg)
        if case != "settling_step":
            ps, m, _ = _run(cfg, sweep, seed)
            out[case] = {
                "ins11_db": m["dynamic_range_db_at_digitiser"],
                "ins11_peak_over_clip_db": m["ins11_residual_peak_over_clip_db"],
                "ins24_margin_db": m["leakage_residual_margin_db"],
                "ins24_step1_margin_db": float(ps["ins24_margin_db"][0]),
            }
            continue
        # One dwell of transmission before the sweep, the filter at rest when it
        # starts: the transmitter switches on one step early.
        fs = inst.sample_rate_for(sweep)
        dwell = float(sweep.step_dwell_s)
        n_pre = int(round(dwell * fs))
        n = int(math.floor(sweep_law.active_duration_s(sweep) * fs))
        rng = np.random.default_rng(seed)
        leak = lk.leakage_signal(sweep, fs, inst.leakage, rng, n_samples=n + n_pre, t0_s=-n_pre / fs,
                                 transmit_power_w=inst.tx_power_w,
                                 tx_chain_efficiency_db=inst.tx_chain_efficiency_db,
                                 impairments=inst.impairments, primed=True)
        if inst.canceller is not None:
            t = -n_pre / fs + np.arange(n + n_pre) / fs
            leak = leak * inst.canceller.residual_factor(sweep, t, inst.leakage, rng)
        x = inst._front(leak)
        y = inst.highpass.apply(x, sweep, fs, t0_s=-n_pre / fs)[n_pre:]
        lev = reference_levels(inst, sweep)
        r = lk.dechirped_residual(y, sweep, fs, levels=lev)
        margin = np.asarray(lev["weakest_echo_dbw"]) - r["level_dbw"]
        from tss.constants import K_BOLTZ

        f_mid = sweep_law.step_centres_hz(sweep)
        t_sys = inst.antenna_temperature_k(f_mid / 1e6) + inst.receiver.noise_temperature_k
        clip_in = float(np.median(K_BOLTZ * t_sys * fs)) * 10.0 ** (inst.receiver.instantaneous_dynamic_range_db / 10.0)
        over = 10.0 * math.log10(float(np.max(np.abs(y) ** 2)) / clip_in)
        out[case] = {
            "ins11_db": inst.receiver.instantaneous_dynamic_range_db - max(0.0, over),
            "ins11_peak_over_clip_db": over,
            "ins24_margin_db": float(np.nanmin(margin)),
            "ins24_step1_margin_db": float(margin[0]),
        }
    return out


# ----------------------------------------------------------------------------
# The study
# ----------------------------------------------------------------------------


def compute(*, seeds: tuple[int, ...] = SEEDS, validate: int = 16) -> dict[str, Any]:
    from tss import sweep_law
    from tss.studies.data_volume import A1, A2

    truth = smoke_truth_echoes()
    rows: list[dict[str, Any]] = []
    chars: dict[tuple[str, str, str, int], dict[str, np.ndarray]] = {}
    carriers: dict[tuple[str, str, str, int], np.ndarray] = {}
    sweeps = {"A1": A1.schedule(), "A2": A2.schedule()}
    sim_cancellers = list(CANCELLERS)
    for (sname, sweep), shape, canc, seed in itertools.product(
        sweeps.items(), ("mom_bus", "flat"), sim_cancellers, seeds
    ):
        chars[(sname, shape, canc, seed)] = characterise(sweep, shape, canc, seed)
        carriers[(sname, shape, canc, seed)] = canceller_carrier_factor(sweep, shape, canc, seed, 0.0)

    all_cancellers = sim_cancellers + list(ACTIVE_FROM)
    for sname, iso_name, canc, am, pw in itertools.product(
        sweeps, ISOLATIONS, all_cancellers, AM_LEVELS_DBC_HZ, POWERS_W
    ):
        shape = "mom_bus" if iso_name == "mom_bus" else "flat"
        base = ACTIVE_FROM.get(canc, canc)
        noise = ACTIVE_NOISE_DBM_HZ if canc in ACTIVE_FROM else None
        per_seed = []
        for seed in seeds:
            ch = chars[(sname, shape, base, seed)]
            pr = predict(ch, iso_name=iso_name, am=am, power_w=pw, added_noise_dbm_hz=noise)
            per_seed.append(pr)
        margins = np.array([float(np.nanmin(p["margin_db"])) for p in per_seed])
        needs = [am_needed_dbc_hz(chars[(sname, shape, base, sd)], iso_name=iso_name, power_w=pw,
                                  added_noise_dbm_hz=noise) for sd in seeds]
        worst = int(np.argmin(margins))
        pr = per_seed[worst]
        ch = chars[(sname, shape, base, seeds[worst])]
        fails = int(np.count_nonzero(pr["margin_db"] < INS24_MARGIN_DB))
        p_fe = pr["p_leak_w"] * carriers[(sname, shape, base, seeds[worst])]
        row: dict[str, Any] = {
            "sweep": sname, "isolation": iso_name, "canceller": canc, "am_dbc_hz": am, "power_w": pw,
            "ins24_margin_db": float(margins.min()), "ins24_margin_db_median_seed": float(np.median(margins)),
            "ins24_worst_freq_hz": float(ch["f_hz"][int(np.nanargmin(pr["margin_db"]))]),
            "ins24_steps_failing": fails, "ins24_met": bool(margins.min() >= INS24_MARGIN_DB),
            "ped_over_noise_max_db": float(np.nanmax(pr["ped_over_noise_db"])),
            "am_needed_dbc_hz": None if any(n is None for n in needs) else float(min(needs)),
        }
        coupler = COUPLER_LOSS_DB if CANCELLERS.get(base) is not None and "ideal" not in canc else 0.0
        if sname == "A1":
            row.update(ins10_predicted(truth, ch["f_hz"], pr["ped_over_noise_db"], power_w=pw,
                                       coupler_loss_db=coupler))
        for p1 in FRONT_END_P1DB_DBM:
            fe = front_end_metrics(p_fe, p1)
            row[f"fe_margin_db_p1db_{p1:+.0f}"] = fe["fe_margin_db"]
            row[f"fe_image_dbc_p1db_{p1:+.0f}"] = fe["fe_image_dbc"]
        rows.append(row)

    # Re-simulate a sample of grid points directly.
    rng = np.random.default_rng(12345)
    checks = []
    pick = [r for r in rows if r["canceller"] not in ACTIVE_FROM]
    for k in rng.choice(len(pick), size=min(validate, len(pick)), replace=False):
        r = pick[int(k)]
        shape, flat = _iso_case(r["isolation"])
        direct = []
        for seed in seeds:
            _, m, _ = _run(_cfg(profile=shape, flat_db=flat, am=r["am_dbc_hz"], canceller=r["canceller"],
                                power_w=r["power_w"]), sweeps[r["sweep"]], seed)
            direct.append(m["leakage_residual_margin_db"])
        checks.append({k2: r[k2] for k2 in ("sweep", "isolation", "canceller", "am_dbc_hz", "power_w")}
                      | {"predicted_db": r["ins24_margin_db"], "direct_db": float(min(direct)),
                         "error_db": float(min(direct) - r["ins24_margin_db"])})

    # INS-10 checked by running the campaign at two grid points.
    ins10_checks = []
    for am, pw, canc in ((AM_REF_DBC_HZ, 1.0, "none"), (-150.0, 0.25, "none"), (-145.0, 1.0, "rf_good")):
        cfg = _cfg(am=am, power_w=pw, canceller=canc)
        direct = smoke_truth_echoes(cfg)["ins10_min_db"]
        row = next(r for r in rows if r["sweep"] == "A1" and r["isolation"] == "mom_bus" and r["canceller"] == canc
                   and r["am_dbc_hz"] == am and r["power_w"] == pw)
        ins10_checks.append({"am_dbc_hz": am, "power_w": pw, "canceller": canc,
                             "predicted_db": row["ins10_min_db"], "direct_db": direct,
                             "error_db": direct - row["ins10_min_db"]})

    startup = {c: start_up(sweeps["A1"], canceller=c) for c in ("none", "rf_good")}
    depths = {}
    from tss.l2_instrument.leakage import RFCanceller

    for name, spec in CANCELLERS.items():
        if spec is not None:
            depths[name] = RFCanceller(**spec).nominal_depth_db()
    return {"rows": rows, "checks": checks, "ins10_checks": ins10_checks, "startup": startup,
            "truth_echoes": {"n": int(truth["snr_db"].size), "ins10_min_ideal_db": truth["ins10_min_db"]},
            "nominal_depth_db": depths,
            "steps": {k: int(sweep_law.step_count(v)) for k, v in sweeps.items()}}


def _row(rows: list[dict[str, Any]], **where: Any) -> dict[str, Any]:
    return next(r for r in rows if all(r.get(k) == v for k, v in where.items()))


def _fmt(v: float | None, spec: str = ".1f", none: str = "none") -> str:
    return none if v is None else format(v, spec)


def _summary(res: dict[str, Any]) -> tuple[str, dict[str, Any]]:
    rows = res["rows"]
    metrics: dict[str, Any] = {}
    base = {"sweep": "A1", "isolation": "mom_bus", "power_w": 1.0}

    lines = ["| canceller | nominal depth (dB) | AM needed (dBc/Hz) | INS-24 at -133 (dB) | at -140 | at -150 | INS-10 at -140 (dB) |",
             "|---|---|---|---|---|---|---|"]
    for canc in list(CANCELLERS) + list(ACTIVE_FROM):
        r133 = _row(rows, **base, canceller=canc, am_dbc_hz=-133.0)
        r140 = _row(rows, **base, canceller=canc, am_dbc_hz=-140.0)
        r150 = _row(rows, **base, canceller=canc, am_dbc_hz=-150.0)
        depth = res["nominal_depth_db"].get(ACTIVE_FROM.get(canc, canc))
        lines.append(f"| {canc} | {_fmt(depth, '.1f', '')} | {_fmt(r133['am_needed_dbc_hz'], '.1f', 'none suffices')} | "
                     f"{r133['ins24_margin_db']:.1f} | {r140['ins24_margin_db']:.1f} | {r150['ins24_margin_db']:.1f} | "
                     f"{r140['ins10_min_db']:.1f} |")
        metrics[f"am_needed_mom_1w_{canc}"] = r133["am_needed_dbc_hz"]
    table_canc = "\n".join(lines)

    lines = ["| canceller | A2: AM needed (dBc/Hz) | A2: INS-24 at -140 (dB) |", "|---|---|---|"]
    for canc in ("none", "ideal26_ins31", "rf_good", "rf_moderate", "rf_one_setting", "rf_good_stepped", "rf_good_active"):
        r = _row(rows, sweep="A2", isolation="mom_bus", power_w=1.0, canceller=canc, am_dbc_hz=-140.0)
        lines.append(f"| {canc} | {_fmt(r['am_needed_dbc_hz'], '.1f', 'none suffices')} | {r['ins24_margin_db']:.1f} |")
        metrics[f"am_needed_a2_mom_1w_{canc}"] = r["am_needed_dbc_hz"]
    table_a2 = "\n".join(lines)

    lines = ["| isolation | AM needed, no canceller | AM needed, rf_good | AM needed, rf_one_setting |", "|---|---|---|---|"]
    for iso in ISOLATIONS:
        vals = [_row(rows, sweep="A1", isolation=iso, power_w=1.0, canceller=c, am_dbc_hz=-133.0)["am_needed_dbc_hz"]
                for c in ("none", "rf_good", "rf_one_setting")]
        lines.append(f"| {iso} | " + " | ".join(_fmt(v, '.1f', 'none suffices') for v in vals) + " |")
    table_iso = "\n".join(lines)

    lines = ["| power (W) | INS-24 at -133, no canceller (dB) | AM needed, no canceller | AM needed, rf_good | "
             "INS-10 min at -133, no canceller (dB) | INS-10 min, rf_good at its AM limit (dB) | echoes below 10 dB |",
             "|---|---|---|---|---|---|---|"]
    for pw in POWERS_W:
        r = _row(rows, sweep="A1", isolation="mom_bus", canceller="none", am_dbc_hz=-133.0, power_w=pw)
        g = _row(rows, sweep="A1", isolation="mom_bus", canceller="rf_good", am_dbc_hz=-133.0, power_w=pw)
        am_g = g["am_needed_dbc_hz"]
        g_at = None
        if am_g is not None:
            ams = [a for a in AM_LEVELS_DBC_HZ if a <= am_g]
            if ams:
                g_at = _row(rows, sweep="A1", isolation="mom_bus", canceller="rf_good", am_dbc_hz=max(ams), power_w=pw)
        lines.append(f"| {pw:g} | {r['ins24_margin_db']:.1f} | {_fmt(r['am_needed_dbc_hz'], '.1f', 'none suffices')} | "
                     f"{_fmt(am_g, '.1f', 'none suffices')} | {r['ins10_min_db']:.1f} | "
                     f"{'' if g_at is None else format(g_at['ins10_min_db'], '.1f')} | {r['ins10_echoes_below']} of {r['ins10_echoes']} |")
    table_pw = "\n".join(lines)

    lines = ["| canceller | AM (dBc/Hz) | " + " | ".join(f"P1dB {p:+.0f} dBm: margin (dB), image (dBc)" for p in FRONT_END_P1DB_DBM) + " |",
             "|---|---|" + "---|" * len(FRONT_END_P1DB_DBM)]
    for canc in ("none", "rf_good", "rf_one_setting"):
        r = _row(rows, **base, canceller=canc, am_dbc_hz=-140.0)
        cells = [f"{r[f'fe_margin_db_p1db_{p:+.0f}']:.1f}, {r[f'fe_image_dbc_p1db_{p:+.0f}']:.1f}" for p in FRONT_END_P1DB_DBM]
        lines.append(f"| {canc} | -140 | " + " | ".join(cells) + " |")
    table_fe = "\n".join(lines)

    err = np.array([c["error_db"] for c in res["checks"]])
    e10 = np.array([c["error_db"] for c in res["ins10_checks"]])
    metrics["validation_max_abs_error_db"] = float(np.max(np.abs(err))) if err.size else None
    metrics["ins10_validation_max_abs_error_db"] = float(np.max(np.abs(e10))) if e10.size else None
    su = res["startup"]
    lines = ["| start-up | canceller | INS-11 (dB) | residual peak over clip (dB) | INS-24 (dB) | INS-24 step 1 (dB) |",
             "|---|---|---|---|---|---|"]
    for canc, cases in su.items():
        for case, v in cases.items():
            lines.append(f"| {case} | {canc} | {v['ins11_db']:.1f} | {v['ins11_peak_over_clip_db']:.1f} | "
                         f"{v['ins24_margin_db']:.1f} | {v['ins24_step1_margin_db']:.1f} |")
    table_su = "\n".join(lines)
    text = "\n".join([
        "Leakage budget, A1 unless stated.  INS-24 margin is the worst step of the worst of "
        f"{len(SEEDS)} seeds against the weakest admitted echo; it must be at least 10 dB.  'AM needed' is the "
        "AM density at which the worst step just meets it, from the characterisation.",
        "",
        f"Model checked against direct simulation at {len(res['checks'])} grid points: largest error "
        f"{metrics['validation_max_abs_error_db']:.2f} dB.  INS-10 prediction checked against the smoke "
        f"campaign at {len(res['ins10_checks'])} points: largest error {metrics['ins10_validation_max_abs_error_db']:.2f} dB.",
        "",
        "## Canceller against AM (MoM bus isolation, 1 W)", "", table_canc, "",
        "## Configuration A2 (MoM bus isolation, 1 W)", "", table_a2, "",
        "## Isolation (1 W)", "", table_iso, "",
        "## Transmit power (MoM bus)", "", table_pw, "",
        "## Front end (MoM bus, 1 W)", "", table_fe, "",
        "## Start-up (MoM bus, -133 dBc/Hz)", "", table_su, "",
    ])
    return text, metrics


@register_study(
    "leakage_budget",
    description="WP3 Phase 5: leakage budget; AM noise x RF canceller x isolation x power x front end x start-up",
    work_packages=["WP3"],
    claims=["C4"],
    status="implemented",
    requirements=["INS-10", "INS-11", "INS-14", "INS-22", "INS-23", "INS-24"],
)
def run(out_dir: Path, *, config: dict[str, Any], manifest: Any) -> StudyResult:
    res = compute()
    rows = res["rows"]
    cols = list(rows[0].keys())
    for r in rows:
        for c in r:
            if c not in cols:
                cols.append(c)
    path = out_dir / "leakage_budget_grid.csv"
    with path.open("w", encoding="utf-8") as fh:
        fh.write(",".join(cols) + "\n")
        for r in rows:
            fh.write(",".join("" if r.get(c) is None else str(r.get(c)) for c in cols) + "\n")
    text, metrics = _summary(res)
    jpath = out_dir / "leakage_budget.json"
    jpath.write_text(json.dumps({k: v for k, v in res.items() if k != "rows"} | {"metrics": metrics},
                                indent=2, default=str), encoding="utf-8")
    return StudyResult(outputs=[path, jpath], summary=text, metrics=metrics)
