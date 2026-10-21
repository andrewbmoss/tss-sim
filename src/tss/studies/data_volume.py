"""What the instrument produces, product by product and mode by mode, per day.

Layer:        study
WP:           WP0
Claims:       -
Status:       implemented
Requirements: INS-13, INS-19, INS-03, INS-04, INS-06, INS-07, INS-27, SCI-09, SCI-18

WP0 closure task 2 (Phase 7.1-7.3).  INS-13 states the data product and its
volume; Rev G gave only the raw ionogram and trace-extraction arithmetic.  This
module defines every product the four observing modes of INS-27 make, sizes it
from the requirements that set its axes, and sums it over a day for a duty
cycle per orbit class.  It stops there: the downlink, the contact plan and the
on-board processor are mission designs (Andrew, 2026-10-20), and what a
mission must size them for is the daily volume of the product set it chooses.

THE PRODUCTS
------------
* **L0 raw**: digitiser samples, 16 bit, at INS-06's minimum rate
  ``f_s = 4 K_f R_max / c`` for the configuration's peak sweep rate, per receive
  channel.  Kept only as calibration bursts (one sweep per orbit by default);
  the continuous figure is given for reference.
* **L1a ionogram**: one amplitude per frequency step and 10 km range bin over
  R_max, 8 bit (INS-13), plus a header (time, position, attitude, sweep
  settings) and one AGC byte per step.  In the vector configuration one matrix
  per receive channel.
* **L1b traces**: up to four echoes per step (O, X, Z and one off-vertical or
  ground echo), each a 4-byte record (virtual range 2 B, amplitude 1 B, mode and
  quality flags 1 B, the INS-16 flag among them).  The vector configuration adds
  the arrival angle (azimuth 2 B, zenith 1 B) and polarization sense (1 B),
  8 bytes per echo.
* **L2 profile**: Ne(h) and its uncertainty (SCI-18) at SCI-06's height sampling
  from h_sat down to hmF2, float32; the goal product adds the scale height H(h)
  and its uncertainty.  A 32-byte block carries NmF2, hmF2, foF2, fxF2 and flags.
* **Doppler product** (Doppler mode, INS-09, INS-15): per sub-band dwell, a
  range-Doppler map, 10 km bins by N Doppler bins at 8 bit, N = 2 v / dv
  sub-sweeps (20 at threshold, 60 at goal), plus up to eight extracted echoes of
  8 bytes (range, line-of-sight velocity, its sigma, amplitude, flag).
* **Drift product** (drift mode, route one of TV Chapter 8): per dwell, the
  phase, amplitude and range of up to four echoes at the drift repetition rate
  (45 Hz), 5 bytes per echo per sub-sweep.
* **Housekeeping**: a fixed rate, 50 B/s by default.

The drift estimate SCI-23 asks for (route two) is made from the displacement of
the same depletion between two nodes' L1b and L2 products, so it adds no
product of its own.

THE PRODUCT SETS
----------------
* ``threshold``: L1b and L2 for every ionogram, full L1a for one ionogram in
  ``l1a_subset`` (10 by default) as a calibration and validation subset, the
  Doppler and drift products, housekeeping and the L0 calibration bursts.
* ``goal``: the same with L1a for every ionogram.
* ``raw_reference``: L0 for every sweep, which is what the instrument makes
  before any on-board processing; not a delivery option.

DUTY CYCLES
-----------
* O1 (about 80 deg, global): ionogram mode at the configuration's sweep all the
  time.  ``O1_6p7s`` holds INS-07's 6.7 s threshold cadence instead.
* O2 (30 deg or below, bubble science): ionogram mode outside 18-03 LT; inside
  it (9 h in 24, SCI-26, MIS-07), a bubble mode in which Doppler dwells take 20
  per cent of the time and drift dwells 10 per cent, and ionograms the rest.
These are reference duty cycles for sizing, not requirements; a mission
declares its own (MIS-08).
"""

from __future__ import annotations

import csv
import json
import math
from dataclasses import asdict, dataclass, replace
from pathlib import Path
from typing import Any

from tss.studies import StudyResult, register_study

C = 299_792_458.0
DAY_S = 86_400.0
MU_EARTH = 3.986004418e14
R_EARTH = 6_378_137.0


@dataclass(frozen=True)
class SweepConfig:
    """A sweep configuration: the INS-02 step law over the INS-01 band."""

    name: str
    f_start_hz: float = 0.33e6
    f_stop_hz: float = 20.0e6
    step_fraction: float = 0.05
    step_cap_hz: float = 150.0e3
    dwell_s: float = 26.8e-3
    n_rx: int = 1

    def schedule(self):
        from tss.sweep_law import stepped_schedule

        return stepped_schedule(self.f_start_hz, self.f_stop_hz, self.dwell_s,
                                step_fraction=self.step_fraction, step_cap_hz=self.step_cap_hz)

    @property
    def steps(self) -> int:
        from tss.sweep_law import step_count

        return step_count(self.schedule())

    @property
    def sweep_s(self) -> float:
        from tss.sweep_law import active_duration_s

        return float(active_duration_s(self.schedule()))

    @property
    def peak_rate_hz_per_s(self) -> float:
        return float(self.schedule().k_f_hz_per_s)


#: Configuration A1: one receive channel, INS-02's 150 kHz cap above f_x.
A1 = SweepConfig("A1")
#: Configuration A2: the vector configuration (INS-18), O and X separated by
#: polarization, so no 150 kHz cap; three receive channels.
A2 = SweepConfig("A2", step_cap_hz=0.0, n_rx=3)


@dataclass(frozen=True)
class Envelope:
    """Spacecraft altitude and the hmF2 floor that set R_max (SCI-08, INS-04)."""

    h_sat_km: float = 900.0
    hmf2_floor_km: float = 200.0

    @property
    def r_max_km(self) -> float:
        return 1.5 * (self.h_sat_km - self.hmf2_floor_km)

    @property
    def orbit_period_s(self) -> float:
        a = R_EARTH + self.h_sat_km * 1e3
        return 2.0 * math.pi * math.sqrt(a ** 3 / MU_EARTH)


THRESHOLD_ENVELOPE = Envelope(900.0)
GOAL_ENVELOPE = Envelope(1400.0)


@dataclass(frozen=True)
class ProductParams:
    """Byte-level choices.  Each default is stated in the module docstring."""

    range_bin_km: float = 10.0
    amplitude_bits: int = 8
    header_bytes: int = 64
    agc_bytes_per_step: int = 1
    echoes_per_step: int = 4
    echo_bytes: int = 4
    echo_bytes_vector_extra: int = 4
    raw_bits: int = 16
    l2_block_bytes: int = 32
    l2_value_bytes: int = 4
    height_step_threshold_km: float = 20.0
    height_step_goal_km: float = 10.0
    doppler_bins_threshold: int = 20
    doppler_bins_goal: int = 60
    doppler_subsweep_s: float = 18e-3
    doppler_echoes: int = 8
    doppler_echo_bytes: int = 8
    drift_rate_hz: float = 45.0
    drift_dwell_s: float = 2.0
    drift_echoes: int = 4
    drift_echo_bytes: int = 5
    housekeeping_bytes_per_s: float = 50.0
    l1a_subset: int = 10
    raw_bursts_per_orbit: float = 1.0


DEFAULT_PARAMS = ProductParams()


@dataclass(frozen=True)
class DutyCycle:
    """Fractions of the day per mode, and the ionogram cadence."""

    name: str
    ionogram_fraction: float = 1.0
    doppler_fraction: float = 0.0
    drift_fraction: float = 0.0
    cadence_s: float | None = None     # None: back to back at the sweep time
    description: str = ""


BUBBLE_FRACTION = 9.0 / 24.0          # 18-03 LT, SCI-26 / MIS-07
DUTY_CYCLES = {
    "O1": DutyCycle("O1", description="ionogram mode continuously at the sweep time"),
    "O1_6p7s": DutyCycle("O1_6p7s", cadence_s=6.7,
                         description="ionogram mode at INS-07's 6.7 s threshold cadence"),
    "O2": DutyCycle("O2", ionogram_fraction=1.0 - BUBBLE_FRACTION * 0.3,
                    doppler_fraction=BUBBLE_FRACTION * 0.2, drift_fraction=BUBBLE_FRACTION * 0.1,
                    description="ionograms outside 18-03 LT; inside it, Doppler 20 %, drift 10 %, "
                                "ionograms the rest"),
}


# ---------------------------------------------------------------------------
# Product sizes, bytes
# ---------------------------------------------------------------------------
def range_bins(env: Envelope, p: ProductParams = DEFAULT_PARAMS) -> int:
    return math.ceil(env.r_max_km / p.range_bin_km - 1e-9)


def l1a_matrix_bytes(cfg: SweepConfig, env: Envelope, p: ProductParams = DEFAULT_PARAMS) -> int:
    """The INS-13 raw ionogram: steps x range bins x amplitude bits, one channel."""
    return cfg.steps * range_bins(env, p) * p.amplitude_bits // 8


def l1a_bytes(cfg: SweepConfig, env: Envelope, p: ProductParams = DEFAULT_PARAMS) -> int:
    """A delivered L1a ionogram: a matrix per channel, the header and the AGC record."""
    return (cfg.n_rx * l1a_matrix_bytes(cfg, env, p) + p.header_bytes
            + cfg.steps * p.agc_bytes_per_step)


def l1b_bytes(cfg: SweepConfig, p: ProductParams = DEFAULT_PARAMS) -> int:
    per_echo = p.echo_bytes + (p.echo_bytes_vector_extra if cfg.n_rx > 1 else 0)
    return cfg.steps * p.echoes_per_step * per_echo + p.header_bytes


def l2_bytes(env: Envelope, p: ProductParams = DEFAULT_PARAMS, goal: bool = False) -> int:
    dh = p.height_step_goal_km if goal else p.height_step_threshold_km
    n_h = math.ceil((env.h_sat_km - env.hmf2_floor_km) / dh - 1e-9) + 1
    per_h = (4 if goal else 2) * p.l2_value_bytes     # Ne, sigma (+ H, sigma_H)
    return n_h * per_h + p.l2_block_bytes + p.header_bytes


def raw_rate_bytes_per_s(cfg: SweepConfig, env: Envelope, p: ProductParams = DEFAULT_PARAMS) -> float:
    fs = 4.0 * cfg.peak_rate_hz_per_s * env.r_max_km * 1e3 / C
    return fs * p.raw_bits / 8 * cfg.n_rx


def raw_sweep_bytes(cfg: SweepConfig, env: Envelope, p: ProductParams = DEFAULT_PARAMS) -> float:
    return raw_rate_bytes_per_s(cfg, env, p) * cfg.sweep_s


def doppler_dwell(env: Envelope, p: ProductParams = DEFAULT_PARAMS, goal: bool = False):
    """(bytes, seconds) of one Doppler-mode sub-band dwell."""
    n = p.doppler_bins_goal if goal else p.doppler_bins_threshold
    size = range_bins(env, p) * n * p.amplitude_bits // 8 + p.doppler_echoes * p.doppler_echo_bytes \
        + p.header_bytes
    return size, n * p.doppler_subsweep_s


def drift_dwell(p: ProductParams = DEFAULT_PARAMS):
    n = round(p.drift_rate_hz * p.drift_dwell_s)
    return n * p.drift_echoes * p.drift_echo_bytes + p.header_bytes, p.drift_dwell_s


# ---------------------------------------------------------------------------
# Daily volumes
# ---------------------------------------------------------------------------
def daily_volume(cfg: SweepConfig, env: Envelope, duty: DutyCycle,
                 p: ProductParams = DEFAULT_PARAMS) -> dict[str, Any]:
    """Bytes per day per product, and the three product sets, for one case.

    The threshold set carries the threshold-fidelity products (L2 at 20 km
    with Ne and sigma; 20-bin Doppler maps) and an L1a subset; the goal set the
    goal-fidelity ones (L2 at 10 km with H(h); 60-bin Doppler maps) and every
    L1a.  Doppler and drift dwells displace ionograms, they are not added on
    top of a full ionogram day.
    """
    period = duty.cadence_s if duty.cadence_s else cfg.sweep_s
    if period < cfg.sweep_s - 1e-9:
        raise ValueError(f"cadence {period} s is shorter than the {cfg.sweep_s:.2f} s sweep")
    if duty.ionogram_fraction + duty.doppler_fraction + duty.drift_fraction > 1.0 + 1e-9:
        raise ValueError(f"duty cycle {duty.name} adds to more than the whole day")
    n_iono = duty.ionogram_fraction * DAY_S / period
    dt_size, dt_t = doppler_dwell(env, p, goal=False)
    dg_size, dg_t = doppler_dwell(env, p, goal=True)
    r_size, r_t = drift_dwell(p)
    n_dr = duty.drift_fraction * DAY_S / r_t
    orbits = DAY_S / env.orbit_period_s
    v = {
        "l1a_all": n_iono * l1a_bytes(cfg, env, p),
        "l1a_subset": n_iono * l1a_bytes(cfg, env, p) / p.l1a_subset,
        "l1b": n_iono * l1b_bytes(cfg, p),
        "l2_threshold": n_iono * l2_bytes(env, p, goal=False),
        "l2_goal": n_iono * l2_bytes(env, p, goal=True),
        "doppler_threshold": duty.doppler_fraction * DAY_S / dt_t * dt_size,
        "doppler_goal": duty.doppler_fraction * DAY_S / dg_t * dg_size,
        "drift": n_dr * r_size,
        "housekeeping": DAY_S * p.housekeeping_bytes_per_s,
        "l0_calibration": orbits * p.raw_bursts_per_orbit * raw_sweep_bytes(cfg, env, p),
        "l0_continuous": DAY_S * raw_rate_bytes_per_s(cfg, env, p),
    }
    base = v["l1b"] + v["drift"] + v["housekeeping"] + v["l0_calibration"]
    sets = {"threshold": base + v["l2_threshold"] + v["doppler_threshold"] + v["l1a_subset"],
            "goal": base + v["l2_goal"] + v["doppler_goal"] + v["l1a_all"],
            "raw_reference": v["l0_continuous"]}
    return {
        "configuration": cfg.name, "n_rx": cfg.n_rx, "h_sat_km": env.h_sat_km,
        "duty_cycle": duty.name, "steps": cfg.steps, "sweep_s": cfg.sweep_s,
        "ionogram_period_s": period, "r_max_km": env.r_max_km, "range_bins": range_bins(env, p),
        "ionograms_per_day": n_iono,
        "doppler_dwells_per_day": {"threshold": duty.doppler_fraction * DAY_S / dt_t,
                                   "goal": duty.doppler_fraction * DAY_S / dg_t},
        "drift_dwells_per_day": n_dr,
        "bytes_per_item": {"l1a": l1a_bytes(cfg, env, p), "l1a_matrix": l1a_matrix_bytes(cfg, env, p),
                           "l1b": l1b_bytes(cfg, p), "l2_threshold": l2_bytes(env, p, False),
                           "l2_goal": l2_bytes(env, p, True), "doppler_dwell_threshold": dt_size,
                           "doppler_dwell_goal": dg_size, "drift_dwell": r_size,
                           "l0_sweep": raw_sweep_bytes(cfg, env, p)},
        "raw_rate_bytes_per_s": raw_rate_bytes_per_s(cfg, env, p),
        "bytes_per_day": v, "sets_bytes_per_day": sets,
    }


def vector_factor(env: Envelope = THRESHOLD_ENVELOPE, p: ProductParams = DEFAULT_PARAMS) -> dict:
    """INS-19: what three receive channels do to each product, same sweep."""
    one = A1
    three = replace(A1, name="A1x3", n_rx=3)
    return {
        "raw": raw_rate_bytes_per_s(three, env, p) / raw_rate_bytes_per_s(one, env, p),
        "l1a": l1a_bytes(three, env, p) / l1a_bytes(one, env, p),
        "l1b": l1b_bytes(three, p) / l1b_bytes(one, p),
        "threshold_set_O1": (daily_volume(three, env, DUTY_CYCLES["O1"], p)["sets_bytes_per_day"]["threshold"]
                             / daily_volume(one, env, DUTY_CYCLES["O1"], p)["sets_bytes_per_day"]["threshold"]),
        "goal_set_O1": (daily_volume(three, env, DUTY_CYCLES["O1"], p)["sets_bytes_per_day"]["goal"]
                        / daily_volume(one, env, DUTY_CYCLES["O1"], p)["sets_bytes_per_day"]["goal"]),
    }


def reference_table(p: ProductParams = DEFAULT_PARAMS) -> list[dict[str, Any]]:
    rows = []
    for env in (THRESHOLD_ENVELOPE, GOAL_ENVELOPE):
        for cfg in (A1, A2):
            for duty in DUTY_CYCLES.values():
                rows.append(daily_volume(cfg, env, duty, p))
    return rows


# ---------------------------------------------------------------------------
# The study
# ---------------------------------------------------------------------------
@register_study(
    "data_volume",
    description="Data products per observing mode and their daily volume per orbit class (INS-13)",
    work_packages=["WP0"],
    claims=[],
    status="implemented",
    requirements=["INS-13", "INS-19"],
)
def run(out_dir: Path, *, config: dict[str, Any], manifest: Any) -> StudyResult:
    opts = dict(config.get("data_volume", {}) or {})
    p = replace(ProductParams(), **{k: v for k, v in opts.items() if k in ProductParams.__dataclass_fields__})
    rows = reference_table(p)
    mb = 1e6
    out_csv = out_dir / "data_volume.csv"
    keys = ["configuration", "n_rx", "h_sat_km", "duty_cycle", "steps", "sweep_s",
            "ionogram_period_s", "r_max_km", "range_bins", "ionograms_per_day"]
    prods = list(rows[0]["bytes_per_day"])
    with open(out_csv, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(keys + [f"{k}_MB_per_day" for k in prods] +
                   ["threshold_set_MB_per_day", "goal_set_MB_per_day", "raw_reference_MB_per_day"])
        for r in rows:
            w.writerow([r[k] for k in keys] + [f"{r['bytes_per_day'][k] / mb:.3f}" for k in prods] +
                       [f"{r['sets_bytes_per_day'][k] / mb:.3f}" for k in ("threshold", "goal", "raw_reference")])
    vf = vector_factor(THRESHOLD_ENVELOPE, p)
    out_md = out_dir / "data_volume.md"
    L = ["# Daily data volume per product set (MB/day)", "",
         "| config | h_sat | duty | steps | sweep s | period s | R_max km | ionograms/day | L1a | L1b | "
         "L2 thr | Doppler thr | drift | L0 cal | threshold set | goal set | raw reference |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        v, s = r["bytes_per_day"], r["sets_bytes_per_day"]
        L.append(f"| {r['configuration']} | {r['h_sat_km']:.0f} | {r['duty_cycle']} | {r['steps']} | "
                 f"{r['sweep_s']:.2f} | {r['ionogram_period_s']:.2f} | {r['r_max_km']:.0f} | "
                 f"{r['ionograms_per_day']:.0f} | {v['l1a_all'] / mb:.0f} | {v['l1b'] / mb:.1f} | "
                 f"{v['l2_threshold'] / mb:.1f} | {v['doppler_threshold'] / mb:.1f} | {v['drift'] / mb:.1f} | "
                 f"{v['l0_calibration'] / mb:.1f} | **{s['threshold'] / mb:.0f}** | **{s['goal'] / mb:.0f}** | "
                 f"{s['raw_reference'] / 1e9:.1f} GB |")
    L += ["", "Vector factor (three channels against one, A1 sweep, 900 km):", "", "```",
          json.dumps(vf, indent=1), "```", ""]
    out_md.write_text("\n".join(L), encoding="utf-8")
    out_json = out_dir / "data_volume.json"
    out_json.write_text(json.dumps({"params": asdict(p), "duty_cycles": {k: asdict(v) for k, v in DUTY_CYCLES.items()},
                                    "rows": rows, "vector_factor": vf}, indent=1), encoding="utf-8")
    ref = daily_volume(A1, THRESHOLD_ENVELOPE, DUTY_CYCLES["O1"], p)
    metrics = {
        "l1a_matrix_bytes_A1_900km": ref["bytes_per_item"]["l1a_matrix"],
        "l1a_matrix_MB_per_day_A1_4p26s": ref["bytes_per_item"]["l1a_matrix"] * DAY_S / A1.sweep_s / mb,
        "l1a_matrix_MB_per_day_A1_6p7s": ref["bytes_per_item"]["l1a_matrix"] * DAY_S / 6.7 / mb,
        "threshold_set_MB_per_day_O1_A1_900km": ref["sets_bytes_per_day"]["threshold"] / mb,
        "goal_set_MB_per_day_O1_A1_900km": ref["sets_bytes_per_day"]["goal"] / mb,
        "vector_factor": vf,
    }
    summary = (f"{len(rows)} cases.  A1 at 900 km: raw ionogram {metrics['l1a_matrix_bytes_A1_900km']} B, "
               f"{metrics['l1a_matrix_MB_per_day_A1_4p26s']:.0f} MB/day back to back; threshold product set "
               f"{metrics['threshold_set_MB_per_day_O1_A1_900km']:.0f} MB/day, goal "
               f"{metrics['goal_set_MB_per_day_O1_A1_900km']:.0f} MB/day (O1).  Three channels: "
               f"threshold set x{vf['threshold_set_O1']:.2f}, goal set x{vf['goal_set_O1']:.2f}.")
    return StudyResult(outputs=[out_csv, out_json, out_md], summary=summary, metrics=metrics)


__all__ = ["SweepConfig", "A1", "A2", "Envelope", "THRESHOLD_ENVELOPE", "GOAL_ENVELOPE",
           "ProductParams", "DutyCycle", "DUTY_CYCLES", "range_bins", "l1a_matrix_bytes", "l1a_bytes",
           "l1b_bytes", "l2_bytes", "raw_rate_bytes_per_s", "raw_sweep_bytes", "doppler_dwell",
           "drift_dwell", "daily_volume", "vector_factor", "reference_table"]
