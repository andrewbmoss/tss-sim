"""Layer: L2
WP: WP3
Claims: C1, C4
Status: implemented
Requirements: INS-10, INS-11, INS-22, INS-23, INS-24

The echo levels the leakage requirements are measured against.

Why this module exists
----------------------
INS-24 asks that the de-chirped leakage pedestal lie at least 10 dB below "the
weakest echo admitted by INS-10".  INS-22 (spur level) and INS-11 (the
high-pass rejection ahead of the digitiser) are derived from the same
comparison.  Until WP3 Phase 0 (2026-10-20) that comparison was a single
number, "leakage over echo = 122 dB - isolation" (SRM Derivations §16, row
179, from the M&S plan), with no frequency, range or antenna attached to it.
Phase 0 traced it: it is the inherited link budget's echo at 5 MHz from a
reflection 300 km away, with the legacy antenna gain of -6.8 dBi at each end
(:class:`tss.l2_instrument.linkbudget.SounderConfig` defaults).  That echo is
neither at INS-04's 1050 km nor on the antenna the instrument flies.

This module computes the reference from the instrument itself, with the
same terms the L2 baseline puts into the samples, so the number INS-24 is
scored against and the samples it is scored on cannot drift apart:

``noise_cell_dbw``
    Noise power in one range cell at the receiver input,
    ``k (T_ant + T_rx) ENBW / T_d``.  ``T_ant = eta T_sky + (1 - eta) T_phys``
    follows the sounding frequency (the SCI-20 sky reference the instrument
    is configured with), ``T_d`` is the step dwell and ``ENBW`` the window's
    equivalent noise bandwidth in bins (1.5 for the Hann window L3 uses).

``weakest_echo_dbw``
    The weakest echo INS-10 admits: ``noise_cell_dbw + 10 dB``.  This is the
    reading of INS-24 this module implements.  Processing losses (window
    scalloping, the step-transition skirt, the estimator) are not folded in;
    :func:`processing_loss_db` measures them against the L2/L3 chain and
    ``tests/unit/test_reference_echo.py`` holds the two together.

``echo_rmax_dbw``
    A specular echo from INS-04's range in free space (no absorption), with
    the same transmit power, antenna gain at both ends and polarization loss
    the L2 baseline applies.  This is the echo the old 122 dB figure stood in
    for, now at the right range and on the right antenna.

``leakage_dbw``
    The leakage carrier at the receiver input, transmit power less INS-23's
    isolation.

and the three ratios that the requirements use: leakage over the weakest
admitted echo (INS-24), leakage over the free-space echo at R_max (INS-22 as
Rev G derives it), and leakage over the noise in one cell (INS-11).

The echo band of each step is ``2 K_i R / c`` for R between ``r_min_m`` and
``r_max_m``, with ``K_i`` that step's sweep rate, which is where INS-24 applies.

Nothing here models leakage.  It says what the leakage has to be measured
against; :mod:`tss.l2_instrument.leakage` (WP3 Phase 1) produces it.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any

import numpy as np

from tss import sweep_law
from tss.constants import C_LIGHT, K_BOLTZ
from tss.interfaces import SweepSchedule
from tss.units import db

__all__ = [
    "INS10_SNR_DB",
    "DEFAULT_R_MIN_M",
    "DEFAULT_ISOLATION_DB",
    "LEGACY_LEAKAGE_OVER_ECHO_PLUS_ISOLATION_DB",
    "window_enbw_bins",
    "reference_levels",
    "summary_rows",
    "required_spur_dbc",
    "legacy_echo_db",
    "processing_loss_db",
]

#: INS-10 threshold, SNR per range cell.
INS10_SNR_DB = 10.0
#: Nearest echo range used for the bottom of the echo band, m.  SRM §22 takes
#: 100 km (the 3.7 kHz beat at 5.6 MHz/s); it is a convention, not a requirement.
DEFAULT_R_MIN_M = 1.0e5
#: INS-23 threshold, used when the instrument does not state an isolation.
DEFAULT_ISOLATION_DB = 40.0
#: SRM Derivations §16 row 179: "leakage over echo = 122 - isolation".
LEGACY_LEAKAGE_OVER_ECHO_PLUS_ISOLATION_DB = 122.0


def window_enbw_bins(window: str = "hann", n: int = 4096) -> float:
    """Equivalent noise bandwidth of an L3 window, in bins."""
    from tss.l3_signal.transform import window_gains, window_taper

    return float(window_gains(window_taper(window, int(n)))[1])


def _isolation_db(instrument: Any, isolation_db: float | None) -> float:
    if isolation_db is not None:
        return float(isolation_db)
    iso: Mapping[str, Any] = getattr(instrument, "isolation", {}) or {}
    return float(iso.get("transmit_receive_db", DEFAULT_ISOLATION_DB))


def reference_levels(
    instrument: Any,
    sweep: SweepSchedule,
    *,
    f_hz: np.ndarray | None = None,
    r_max_m: float | None = None,
    r_min_m: float = DEFAULT_R_MIN_M,
    snr_db: float = INS10_SNR_DB,
    isolation_db: float | None = None,
    window: str = "hann",
) -> dict[str, np.ndarray]:
    """The INS-10/INS-24 reference levels at each frequency, as parallel arrays.

    ``instrument`` is an L2 :class:`~tss.l2_instrument.baseline.BaselineInstrument`
    (its antenna gain, antenna temperature, receiver noise temperature, transmit
    power and losses are used, so the numbers are the ones its samples carry).
    ``sweep`` must be stepped; ``f_hz`` defaults to its step centres.
    ``r_max_m`` defaults to the instrument's design range (INS-04).
    Powers are in dBW at the receiver input.  Nothing is modified.
    """
    if not sweep_law.is_stepped(sweep):
        raise ValueError("reference_levels needs a stepped sweep (one dwell per step)")
    dwell = float(sweep.step_dwell_s)
    centres = sweep_law.step_centres_hz(sweep)
    f = np.asarray(centres if f_hz is None else f_hz, dtype=float)
    if np.any(f <= 0.0):
        raise ValueError("frequencies must be positive")
    r_max = float(getattr(instrument, "design_range_m", 1.05e6) if r_max_m is None else r_max_m)
    iso = _isolation_db(instrument, isolation_db)

    edges = sweep_law.step_edges_hz(sweep)
    rates = sweep_law.step_rates_hz_per_s(sweep)
    idx = np.clip(np.searchsorted(edges, f, side="right") - 1, 0, rates.size - 1)
    k_i = np.asarray(rates, dtype=float)[idx]

    f_mhz = f / 1.0e6
    gain = np.asarray(instrument.gain_dbi(f_mhz), dtype=float)
    t_ant = np.asarray(instrument.antenna_temperature_k(f_mhz), dtype=float)
    t_rx = float(instrument.receiver.noise_temperature_k)
    t_sys = t_ant + t_rx
    enbw = window_enbw_bins(window)
    b_cell = enbw / dwell
    noise = db(K_BOLTZ * t_sys * b_cell)

    p_tx = db(float(instrument.tx_power_w)) + float(instrument.tx_chain_efficiency_db)
    lam = C_LIGHT / f
    spreading = 20.0 * np.log10(lam / (4.0 * math.pi * 2.0 * r_max))
    losses = float(instrument.polarization_loss_db) + float(instrument.extra_loss_db)
    echo_rmax = p_tx + 2.0 * gain + spreading + losses
    weakest = noise + float(snr_db)
    leak = np.full(f.shape, p_tx - iso)

    return {
        "f_hz": f,
        "step_rate_hz_per_s": k_i,
        "beat_min_hz": 2.0 * k_i * float(r_min_m) / C_LIGHT,
        "beat_max_hz": 2.0 * k_i * r_max / C_LIGHT,
        "gain_dbi": gain,
        "t_ant_k": t_ant,
        "t_sys_k": t_sys,
        "cell_bandwidth_hz": np.full(f.shape, b_cell),
        "noise_cell_dbw": noise,
        "weakest_echo_dbw": weakest,
        "echo_rmax_dbw": echo_rmax,
        "snr_rmax_db": echo_rmax - noise,
        "leakage_dbw": leak,
        "isolation_db": np.full(f.shape, iso),
        "leakage_over_weakest_db": leak - weakest,
        "leakage_over_rmax_echo_db": leak - echo_rmax,
        "leakage_over_noise_db": leak - noise,
    }


def summary_rows(
    levels: Mapping[str, np.ndarray], at_mhz: tuple[float, ...] = (0.33, 0.5, 1, 2, 3, 5, 10, 20)
) -> list[dict[str, float]]:
    """The row of ``levels`` nearest each frequency in ``at_mhz``."""
    f = np.asarray(levels["f_hz"], dtype=float) / 1.0e6
    out: list[dict[str, float]] = []
    for target in at_mhz:
        i = int(np.argmin(np.abs(f - float(target))))
        out.append({k: float(np.asarray(v)[i]) for k, v in levels.items()})
    return out


def required_spur_dbc(
    levels: Mapping[str, np.ndarray],
    *,
    coupling_path_m: float = 10.0,
    margin_db: float = 10.0,
    reference: str = "weakest",
) -> np.ndarray:
    """Transmit spur level that keeps the spread leakage ``margin_db`` below the reference.

    The SRM §16 rule, ``-(L/E + margin) - S``, with ``S = 4 sin^2(pi f_b tau)``
    the range-correlation suppression at the top of each step's echo band
    (the longest range, where it is least) and ``tau`` the coupling delay.
    ``reference`` is ``"weakest"`` (the INS-10-admitted echo, INS-24's
    reading) or ``"rmax"`` (the free-space echo at R_max, Rev G's reading).
    Single-tone, small-angle, as in Rev G §5.9.
    """
    key = {"weakest": "leakage_over_weakest_db", "rmax": "leakage_over_rmax_echo_db"}[reference]
    tau = float(coupling_path_m) / C_LIGHT
    fb = np.asarray(levels["beat_max_hz"], dtype=float)
    suppression = 10.0 * np.log10(4.0 * np.sin(math.pi * fb * tau) ** 2)
    return -(np.asarray(levels[key], dtype=float) + float(margin_db)) - suppression


def legacy_echo_db(f_hz: float = 5.0e6) -> dict[str, float]:
    """The inherited budget's echo, relative to transmit power, and what it implies.

    :class:`~tss.l2_instrument.linkbudget.SounderConfig` defaults: 300 km
    reflection, the legacy antenna (-6.8 dBi at 5 MHz), -6 dB polarization,
    integrated absorption on its Chapman layer.  This reproduces the 122 dB of
    SRM Derivations §16 to within its rounding, which is the evidence that the
    figure is this budget and not the instrument.
    """
    from tss.l2_instrument import linkbudget as lb

    cfg = lb.SounderConfig()
    r = lb.link_budget(float(f_hz), cfg)
    return {
        "f_hz": float(f_hz),
        "true_range_m": float(cfg.true_range),
        "tx_gain_dbi": float(r.tx_gain_dbi),
        "echo_below_tx_db": float(r.tx_power_dbw - r.rx_power_dbw),
        "absorption_db": float(r.absorption_db),
    }


def processing_loss_db(
    instrument: Any,
    sweep: SweepSchedule,
    f_hz: float,
    range_m: float,
    *,
    snr_db: float = 20.0,
    seeds: tuple[int, ...] = (1, 2, 3, 4),
    l3: Any = None,
) -> dict[str, float]:
    """Predicted against measured SNR for one synthetic echo, through L2 and L3.

    One ordinary-mode echo at ``f_hz`` (moved to its step centre by the
    caller if wanted) and ``range_m``, with the path gain set so that
    :func:`reference_levels` predicts ``snr_db``.  The L2 instrument makes the
    samples (noise, receiver, quantiser) and the L3 baseline measures the SNR
    at the echo's cell exactly as INS-10 is scored.  Returns the mean
    measured SNR over ``seeds`` and the loss, predicted less measured.
    """
    from tss.interfaces import EchoSet

    if l3 is None:
        from tss.registry import build

        l3 = build("l3", "baseline")
    lev = reference_levels(instrument, sweep, f_hz=np.array([float(f_hz)]), r_max_m=float(range_m))
    # path gain that puts the echo at the requested SNR
    lam = C_LIGHT / float(f_hz)
    spreading = 20.0 * math.log10(lam / (4.0 * math.pi * 2.0 * float(range_m)))
    path_gain = spreading + (float(snr_db) - float(lev["snr_rmax_db"][0]))
    echoes = EchoSet(
        freq_hz=np.array([float(f_hz)]),
        mode=np.array(["O"], dtype=object),
        group_delay_s=np.array([2.0 * float(range_m) / C_LIGHT]),
        amplitude_db=np.array([path_gain]),
        absorption_db=np.array([0.0]),
        flags=np.array([""], dtype=object),
        source="reference_echo.processing_loss_db",
    )
    measured = []
    for s in seeds:
        rng = np.random.default_rng(int(s))
        signal = instrument.receive(echoes, sweep, rng)
        measured.append(float(np.asarray(l3.snr_at(signal, [float(f_hz)], [float(range_m)]))[0]))
    m = np.asarray(measured, dtype=float)
    mean_lin = float(db(np.mean(10.0 ** (m / 10.0))))
    return {
        "f_hz": float(f_hz),
        "range_m": float(range_m),
        "predicted_snr_db": float(snr_db),
        "measured_snr_db_mean": mean_lin,
        "measured_snr_db_spread": float(np.max(m) - np.min(m)),
        "loss_db": float(snr_db) - mean_lin,
    }
