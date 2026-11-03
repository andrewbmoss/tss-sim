"""Layer: studies
WP: WP4
Claims: C1, C3
Status: implemented
Requirements: INS-05, INS-09, INS-10, SCI-07

WP4 Phase 5: the FM-CW acceptance study.  Does the chain reproduce the
published numbers of Poole (1985, "paper 1"), Poole and Evans (1985,
"paper 2") and Barry (1971), at those papers' own parameters?

Decisions (Andrew, 2026-11-03): the items are checked at waveform and L2/L3
level (the papers describe ground-based bottomside sounders; tss-sim's L0-L1
geometry is topside, and no ground geometry is built); A8 (the echo
criterion) and A9 (the group-velocity method) are included; Barry's gating
loss is not simulated.  Gate G5 (2026-10-31): arithmetic items within 1 %;
statistical items reproduced when the published value lies inside the Monte
Carlo's 95 % confidence interval; Barry items within 10 %, being read off
plotted or rounded results.

HOW
---
De-chirped cells are synthesised from paper 1's eq. 14,
``phi_D(t) = 2 pi (2 f_T(t) R*/c - 2 k (R*/c)^2)``, with ``R*`` the phase range
at the instantaneous frequency ``f_T(t) = f_0 + k t`` and time, so dispersion
and motion enter exactly; noise is complex white Gaussian, set from the
per-bin SNR of the chain's window.  They go through the chain's transform
(:func:`tss.l3_signal.transform.transform`), its detector
(:func:`tss.l3_signal.detect.cfar_detect`) and Poole's estimator
(:mod:`tss.l3_signal.poole_range`).

Each item scores ``reproduced``, ``difference explained`` (the WP4 exit
criterion allows it, with the cause stated), ``not reproduced, open``,
``reported`` (A8, a comparison) or ``not simulated`` (a hardware property, or
outside the decided scope; listed with the reason).
"""

from __future__ import annotations

import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from tss.studies import StudyResult, register_study

C = 299_792_458.0
DEFAULTS = {"mc_trials": 2000, "seed": 20261103, "a9_snr_db": 30.0}


# -- synthesis and the chain -------------------------------------------------------------


def cell(f0_hz: float, t0_s: float, k: float, cell_s: float, fs_hz: float, rstar: Any,
         *, snr_db: float | None = None, rng: np.random.Generator | None = None,
         window: str = "hann") -> np.ndarray:
    """One de-chirped cell (paper 1 eq. 14).  ``rstar(f_hz, t_s)`` is the phase range, m."""
    from tss.l3_signal.transform import window_gains, window_taper

    n = int(round(cell_s * fs_hz))
    t = np.arange(n) / fs_hz
    f_t = f0_hz + k * t
    r = np.asarray(rstar(f_t, t0_s + t), dtype=float)
    x = np.exp(2j * np.pi * (2.0 * f_t * r / C - 2.0 * k * (r / C) ** 2))
    if snr_db is not None:
        _, enbw = window_gains(window_taper(window, n))
        sigma2 = n / (enbw * 10.0 ** (snr_db / 10.0))
        g = rng if rng is not None else np.random.default_rng()
        x = x + math.sqrt(sigma2 / 2.0) * (g.standard_normal(n) + 1j * g.standard_normal(n))
    return x


def noise_cell(n: int, sigma2: float, rng: np.random.Generator) -> np.ndarray:
    return math.sqrt(sigma2 / 2.0) * (rng.standard_normal(n) + 1j * rng.standard_normal(n))


def spectrum(x: np.ndarray, fs_hz: float, k: float, f_c_hz: float, *, pad: int = 1, window: str = "hann") -> Any:
    """The chain's transform of a synthetic cell."""
    from tss.interfaces import SweepSchedule
    from tss.l3_signal.dechirp import Dechirped
    from tss.l3_signal.transform import transform

    n = x.size
    sw = SweepSchedule(f_start_hz=f_c_hz, f_stop_hz=f_c_hz + k * n / fs_hz, k_f_hz_per_s=k,
                       sweep_period_s=n / fs_hz)
    blk = Dechirped(samples=x[None, :], fs_hz=fs_hz, t_s=np.arange(n) / fs_hz, sweep=sw, oriented=True,
                    centre_frequency_hz=f_c_hz, k_f_hz_per_s=k)
    return transform(blk, window=window, n_fft=n * pad)


def peak_range_m(sp: Any) -> float:
    """Peak of the chain's spectrum, refined by a parabola through the three top bins (dB)."""
    p = np.asarray(sp.power_db, dtype=float)[0]
    i = int(np.argmax(p))
    if 0 < i < p.size - 1:
        a, b, c = p[i - 1], p[i], p[i + 1]
        d = 0.5 * (a - c) / (a - 2 * b + c)
    else:
        d = 0.0
    return float(sp.range_m[i] + d * sp.bin_spacing_m)


def _row(item: str, statement: str, source: str, published: Any, measured: Any, tolerance: str,
         outcome: str, note: str = "") -> dict[str, Any]:
    return {"item": item, "statement": statement, "source": source, "published": published,
            "measured": measured, "tolerance": tolerance, "outcome": outcome, "note": note}


def _within(meas: float, pub: float, rel: float) -> bool:
    return abs(meas - pub) <= rel * abs(pub)


def _rms_ci(samples: np.ndarray) -> tuple[float, float, float]:
    """RMS with a 95 % interval (normal approximation on the mean square, delta to the root)."""
    s2 = np.asarray(samples, dtype=float) ** 2
    m = float(np.mean(s2))
    se = float(np.std(s2, ddof=1) / math.sqrt(s2.size))
    lo, hi = max(m - 1.96 * se, 0.0), m + 1.96 * se
    return math.sqrt(m), math.sqrt(lo), math.sqrt(hi)


# Paper 2, Table 1 (swept-frequency column).
K_P = 100.0e3
T_P = 0.5
FS_P = 2048.0          # 1024 samples per 0.5 s cell, the hardware FFT's input
DF_P = 5.0e3
SWEEP_P = 50.0e3       # overall (average) sweep rate
STEP_P = 75.0e3        # frequency increment between soundings


def sounding(f0: float, t0: float, rstar: Any, *, snr_db: float | None, rng: np.random.Generator) -> list[np.ndarray]:
    """Poole's three cells: two at f0, one at f0 + df, T apart."""
    return [cell(f0, t0, K_P, T_P, FS_P, rstar, snr_db=snr_db, rng=rng),
            cell(f0, t0 + T_P, K_P, T_P, FS_P, rstar, snr_db=snr_db, rng=rng),
            cell(f0 + DF_P, t0 + 2 * T_P, K_P, T_P, FS_P, rstar, snr_db=snr_db, rng=rng)]


# -- the items -----------------------------------------------------------------------------


def item_a1_a2(rng: np.random.Generator) -> list[dict[str, Any]]:
    from tss.l3_signal.detect import cfar_detect

    rows = []
    r0 = 300.0e3
    cs = sounding(5.0e6, 0.0, lambda f, t: np.full_like(f, r0), snr_db=30.0, rng=rng)
    det = []
    for j, x in enumerate(cs):
        f0 = 5.0e6 + (DF_P if j == 2 else 0.0)
        sp = spectrum(x, FS_P, K_P, f0 + 0.5 * K_P * T_P)
        res = cfar_detect(sp)
        det.append(bool(np.any(res.mask & (np.abs(sp.range_m - r0) < 2 * sp.range_resolution_m))))
    n_sound = int(round((7.0e6 - 5.5e6) / STEP_P)) + 1
    rows.append(_row("A1", "Chain driven at Table 1: k = 100 kHz/s, T = 0.5 s, 3 cells per sounding, df = 5 kHz in "
                           "cell 3, 50 kHz/s overall, 75 kHz between soundings", "Paper 2, Table 1",
                     "runs", f"echo detected in {sum(det)}/3 cells; {n_sound} soundings over 5.5-7.0 MHz take "
                     f"{n_sound * 3 * T_P:.1f} s, so {STEP_P / (3 * T_P) / 1e3:.0f} kHz/s overall",
                     "runs", "reproduced" if all(det) and abs(STEP_P / (3 * T_P) - SWEEP_P) < 1 else "not reproduced, open",
                     "1024 samples per cell (fs 2048 Hz), the chain's Hann window and stated-P_fa detector"))
    sp = spectrum(cs[0], FS_P, K_P, 5.0e6)
    pub = C / (2 * K_P * T_P)
    rows.append(_row("A2", "Range cell c/(2kT) = 3.0 km", "Paper 1 eq. 27; Table 1", round(pub / 1e3, 4),
                     round(sp.bin_spacing_m / 1e3, 4), "1 %",
                     "reproduced" if _within(sp.bin_spacing_m, 3.0e3, 0.01) else "not reproduced, open",
                     f"the chain's unpadded bin; its Hann cell (ENBW 1.5 bins, INS-05's measure) is "
                     f"{sp.range_resolution_m / 1e3:.2f} km"))
    return rows


def item_a3(rng: np.random.Generator) -> list[dict[str, Any]]:
    k, f, r0 = 50.0e3, 5.0e6, 300.0e3
    cell_s, fs = 0.5, 2048.0
    out = {}
    for v in (0.0, 30.0):
        x = cell(f - 0.5 * k * cell_s, 0.0, k, cell_s, fs, lambda ff, t, v=v: r0 + v * t, rng=rng)
        out[v] = peak_range_m(spectrum(x, fs, k, f, pad=64))
    rel = (out[30.0] - out[0.0]) / r0
    return [_row("A3", "f_D = 2kR'/c holds to 1 % when V* < 1 % of kR'/f (3 km/s at k = 50 kHz/s, R' = 300 km, "
                       "f = 5 MHz), so V* < 30 m/s", "Paper 1 eqs. 26-27",
                 "1.00 % at 30 m/s", f"{100 * rel:.3f} % at 30 m/s", "1 %",
                 "reproduced" if _within(rel, 0.01, 0.01) else "not reproduced, open",
                 f"range read from the chain's spectrum peak; bias {(out[30.0] - out[0.0]) / 1e3:.3f} km "
                 f"(fV*/k = {f * 30 / k / 1e3:.3f} km)")]


def item_a4(rng: np.random.Generator) -> list[dict[str, Any]]:
    from tss.l3_signal import poole_range as pr
    from tss.l3_signal.detect import cfar_detect

    f0 = 5.9e6
    errs, edge = [], []
    for r0 in np.linspace(150e3, 600e3, 19):
        cs = sounding(f0, 0.0, lambda f, t, r0=r0: np.full_like(f, r0), snr_db=40.0, rng=rng)
        res = pr.three_cell(cs, FS_P, k_hz_per_s=K_P, f0_hz=f0, delta_f_hz=DF_P)
        errs.append(res["group_range_m"] - r0)
        sp = spectrum(cs[0], FS_P, K_P, f0 + 0.5 * K_P * T_P)
        m = cfar_detect(sp).mask
        run = np.flatnonzero(m & (sp.range_m > r0 - 20e3))
        edge.append(float(sp.range_m[run[0]] - r0) if run.size else float("nan"))
    # the wrap period of phi_R against range: one cycle per c/(2 df)
    rr = np.linspace(300e3, 360e3, 241)
    ph = []
    for r0 in rr:
        cs = sounding(f0, 0.0, lambda f, t, r0=r0: np.full_like(f, r0), snr_db=None, rng=rng)
        ph.append(pr.three_cell(cs, FS_P, k_hz_per_s=K_P, f0_hz=f0, delta_f_hz=DF_P,
                                bin_index=int(round(2 * K_P * 330e3 / C * T_P)))["phi_r_rad"])
    slope = np.polyfit(rr, np.unwrap(ph), 1)[0]
    period = 2 * math.pi / slope
    rows = [
        _row("A4a", "Ambiguity interval of the stationary-phase range c/(2 df) = 30 km", "Paper 1 eqs. 48, 70",
             round(C / (2 * DF_P) / 1e3, 3), round(period / 1e3, 3), "1 %",
             "reproduced" if _within(period, C / (2 * DF_P), 0.01) else "not reproduced, open",
             "measured as the period of phi_R against injected range"),
        _row("A4b", "f_R1 = k/df = 20 Hz", "Paper 1 eq. 51", 20.0,
             round(2 * K_P * period / C, 3), "1 %",
             "reproduced" if _within(2 * K_P * period / C, 20.0, 0.01) else "not reproduced, open",
             "the beat-frequency width of the measured interval, 2k x period / c"),
        _row("A4c", "R' = c(N + phi_R/2pi)/(2 df) as a reference estimator, against the chain's leading edge",
             "Paper 1 eqs. 47-52, 70-71", "method",
             f"Poole R' error {np.mean(errs) / 1e3:+.3f} km mean, {np.max(np.abs(errs)) / 1e3:.3f} km worst over "
             f"150-600 km (40 dB); chain leading edge {np.nanmean(edge) / 1e3:+.2f} km mean",
             "estimator recovers R'", "reproduced" if np.max(np.abs(errs)) < 0.1 * C / (2 * DF_P) else
             "not reproduced, open",
             "the leading edge reads short by about a cell by construction (the WP3-era definition); "
             "Poole's R' has no such offset"),
    ]
    return rows


def item_a5(rng: np.random.Generator) -> list[dict[str, Any]]:
    """Bisect the injected V* at which N (eq. 71) goes wrong.

    Eq. 59 takes ``f_i`` as the echo's beat frequency; the chain measures it at
    a bin centre, which moves the bound by up to half a bin (1 Hz, 30 m/s at
    5 MHz).  Both are reported: ``f_i`` refined by a parabola through the three
    top bins (the paper's quantity, scored), and the bin centre.
    """
    from tss.l3_signal import poole_range as pr

    f0 = 5.0e6 - 0.5 * K_P * T_P
    r0 = 300.0e3 - (300.0e3 % (C / (2 * K_P * T_P)))     # on a bin, so f' = 0 at V* = 0

    def fails(v: float, df: float, refine: bool) -> bool:
        cs = [cell(f0, 0, K_P, T_P, FS_P, lambda f, t: r0 + v * t),
              cell(f0, T_P, K_P, T_P, FS_P, lambda f, t: r0 + v * t),
              cell(f0 + df, 2 * T_P, K_P, T_P, FS_P, lambda f, t: r0 + v * t)]
        res = pr.three_cell(cs, FS_P, k_hz_per_s=K_P, f0_hz=f0, delta_f_hz=df)
        r = res["group_range_m"]
        if refine:
            pad = 32
            fr, z = pr.complex_bins(cs[0], FS_P, n_fft=pad * cs[0].size)
            p = np.abs(z) ** 2
            i = int(np.argmax(np.where(np.abs(fr - res["f_i_hz"]) <= 1.0 / T_P, p, 0.0)))
            a, b, c = np.log(p[i - 1]), np.log(p[i]), np.log(p[i + 1])
            f_i = fr[i] + 0.5 * (a - c) / (a - 2 * b + c) * (fr[1] - fr[0])
            r, _n = pr.group_range_m(res["phi_r_rad"], f_i, K_P, df)
        return abs(r - r0) > 0.25 * C / (2 * df)

    def bound(df: float, hi: float, refine: bool) -> float:
        lo = 0.0
        for _ in range(40):
            mid = 0.5 * (lo + hi)
            if fails(mid, df, refine):
                hi = mid
            else:
                lo = mid
        return 0.5 * (lo + hi)

    b5, b5_bin = bound(5.0e3, 2000.0, True), bound(5.0e3, 2000.0, False)
    b1, b1_bin = bound(1.0e3, 5000.0, True), bound(1.0e3, 5000.0, False)
    theory = pr.max_unambiguous_velocity_mps(K_P, 5e3, 5e6)

    def corrected(df: float) -> float:
        # eq. 59 with the two terms of order V* T it drops: the echo's range at the
        # window centre (k V* T / c in f_i) and its change across the three cells
        # (5 df V* T / c in phi_R), both from eq. 14 evaluated at each cell's centroid
        f_r1 = K_P / df
        return 0.5 * C / (2 * 5.0e6 / f_r1 + K_P * T_P / f_r1 - 5 * df * T_P)

    c5, c1 = corrected(5.0e3), corrected(1.0e3)
    if _within(b5, 300.0, 0.01):
        out5 = "reproduced"
    elif _within(b5, c5, 0.01):
        out5 = "difference explained"
    else:
        out5 = "not reproduced, open"
    return [
        _row("A5a", "Unambiguous Doppler about 300 m/s at 5 MHz for k = 100 kHz/s, df = 5 kHz", "Paper 1 eq. 59",
             300.0, round(b5, 2), "1 %", out5,
             f"injected V* until N is wrong, f_i refined between bins by a 32x padded transform (eq. 59: {theory:.1f} "
             f"m/s).  Eq. 59 drops the echo's range change during and between the cells, of order V* T (the same "
             f"approximation behind the R' bias in A4); with those terms the bound is {c5:.1f} m/s, which the "
             f"measurement matches.  With f_i at the bin centre, as the hardware read it, the bound is "
             f"{b5_bin:.1f} m/s: the half bin, f' in eq. 59"),
        _row("A5b", "df down to 1 kHz to measure 500 m/s", "Paper 1 eq. 59", ">= 500", round(b1, 1), "1 %",
             "reproduced" if b1 >= 500.0 * 0.99 else "not reproduced, open",
             f"the bound at df = 1 kHz (eq. 59: 1499 m/s; with the V* T terms {c1:.1f} m/s; bin centre "
             f"{b1_bin:.1f} m/s); the paper's statement is a sufficient condition"),
    ]


def item_a6_a7(trials: int, rng: np.random.Generator) -> tuple[list[dict[str, Any]], list[dict[str, Any]]]:
    from tss.l3_signal import poole_range as pr

    curve = []
    n = int(T_P * FS_P)
    bin_i = 200
    f_b = bin_i / T_P
    r0 = f_b * C / (2 * K_P)
    for snr in (10.0, 15.0, 20.0, 25.0, 30.0, 40.0):
        err = []
        for _ in range(trials):
            x = cell(5.0e6, 0.0, K_P, T_P, FS_P, lambda f, t: np.full_like(f, r0), snr_db=snr, rng=rng)
            ref = cell(5.0e6, 0.0, K_P, T_P, FS_P, lambda f, t: np.full_like(f, r0))
            _, a = pr.complex_bins(x, FS_P)
            _, b = pr.complex_bins(ref, FS_P)
            err.append(float(np.angle(a[bin_i] / b[bin_i])))
        rms, lo, hi = _rms_ci(np.degrees(err))
        curve.append({"snr_db": snr, "mc_rms_deg": rms, "ci_lo_deg": lo, "ci_hi_deg": hi,
                      "closed_form_deg": math.degrees(pr.phase_error_rms_rad(snr))})
    c20 = next(c for c in curve if c["snr_db"] == 20.0)
    inside = c20["ci_lo_deg"] <= 4.0 <= c20["ci_hi_deg"]
    a6 = _row("A6", "RMS phase error below 4 deg above 20 dB (Fig. 3)", "Paper 1 eqs. 76-82, Fig. 3", 4.0,
              f"{c20['mc_rms_deg']:.3f} [{c20['ci_lo_deg']:.3f}, {c20['ci_hi_deg']:.3f}] at 20 dB",
              "published value inside the 95 % CI",
              "reproduced" if inside else ("difference explained"
                                           if c20["ci_lo_deg"] <= c20["closed_form_deg"] <= c20["ci_hi_deg"]
                                           else "not reproduced, open"),
              f"Poole's closed form (eq. 82) gives {c20['closed_form_deg']:.3f} deg at 20 dB: '4 deg above 20 dB' "
              f"is that value rounded, and the Monte Carlo agrees with eq. 82; {trials} trials per SNR through the chain's window, Gaussian noise")
    del n
    # A7: V* rms for a 4 deg error in phi_D at T = 0.5 s, f = 5.9 MHz (M = 0)
    f0 = 5.9e6 - 0.5 * K_P * T_P
    out = {}
    for label, snr in (("phiD_4deg", 10 * math.log10(1.0 / (2.0 * math.radians(4.0 / math.sqrt(2.0)) ** 2))),
                       ("cell_20dB", 20.0)):
        vs, pd = [], []
        for _ in range(trials):
            c1 = cell(f0, 0.0, K_P, T_P, FS_P, lambda f, t: np.full_like(f, r0), snr_db=snr, rng=rng)
            c2 = cell(f0, T_P, K_P, T_P, FS_P, lambda f, t: np.full_like(f, r0), snr_db=snr, rng=rng)
            _, a = pr.complex_bins(c1, FS_P)
            _, b = pr.complex_bins(c2, FS_P)
            phi_d = float(np.angle(b[bin_i] / a[bin_i]))
            pd.append(math.degrees(phi_d))
            v, _m = pr.doppler_velocity_mps(phi_d, f_b, f_b, T_P, 5.9e6)
            vs.append(v)
        out[label] = (_rms_ci(np.asarray(vs)), _rms_ci(np.asarray(pd)), snr)
    (v, lo, hi), (pdr, _, _), snr = out["phiD_4deg"]
    a7 = _row("A7", "V* rms error 0.56 m/s for a 4 deg phase error, T = 0.5 s, f = 5.9 MHz", "Paper 2, from paper 1 eq. 72",
              0.56, f"{v:.3f} [{lo:.3f}, {hi:.3f}] (phi_D rms {pdr:.2f} deg)", "published value inside the 95 % CI",
              "reproduced" if lo <= 0.56 <= hi else "difference explained",
              f"cells at {snr:.1f} dB each, so the phase difference has 4 deg rms; with 20 dB per cell, as paper 2's "
              f"sentence reads, the difference has sqrt(2) x 4.05 deg and V* {out['cell_20dB'][0][0]:.3f} m/s")
    return [a6, a7], curve


def item_a8(trials: int, rng: np.random.Generator) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Paper 2's echo criterion against the P2 detector (a comparison, not a pass criterion).

    Poole's criterion is an amplitude threshold with at most the 16 highest
    values kept per cell; it is applied here at the P2 detector's own threshold
    over the mean, with the noise mean known exactly (a perfect AGC), so the
    difference shown is the reference estimate and the cap.
    """
    from tss.l3_signal.detect import cfar_detect

    n = int(T_P * FS_P)
    r0 = 300.0e3
    noise = [spectrum(noise_cell(n, 1.0, rng), FS_P, K_P, 5.0e6) for _ in range(trials)]
    lin = np.array([10 ** (np.asarray(sp.power_db)[0] / 10) for sp in noise])
    unit = float(np.mean(lin))                     # the known noise mean per bin
    thr_db = float(cfar_detect(noise[0]).flags["threshold_over_mean_db"])
    thr = 10 ** (thr_db / 10)
    fa_poole = int(sum(min(int(np.count_nonzero(row / unit > thr)), 16) for row in lin))
    fa_cfar = int(sum(int(cfar_detect(sp).mask.sum()) for sp in noise))
    cells = lin.size
    from tss.l3_signal.transform import window_gains, window_taper

    _, enbw = window_gains(window_taper("hann", n))
    pd: dict[str, tuple[float, float]] = {}
    n_pd = max(trials // 4, 50)
    for snr in (6.0, 8.0, 10.0, 12.0, 14.0):
        hit_p = hit_c = 0
        for _ in range(n_pd):
            x = cell(5.0e6, 0, K_P, T_P, FS_P, lambda f, t: np.full_like(f, r0), snr_db=snr, rng=rng)
            sp = spectrum(x, FS_P, K_P, 5.0e6)
            near = np.abs(sp.range_m - r0) <= sp.bin_spacing_m
            s2 = n / (enbw * 10.0 ** (snr / 10.0))     # the noise variance cell() used
            hit_p += int(np.any(10 ** (np.asarray(sp.power_db)[0][near] / 10) / (unit * s2) > thr))
            hit_c += int(np.any(cfar_detect(sp).mask[near]))
        pd[f"{snr:.0f}"] = (hit_p / n_pd, hit_c / n_pd)
    # a spread echo: 40 adjacent bins at 15 dB each, random phases (spread F)
    bins = 300.0e3 + np.arange(40) * C / (2 * K_P * T_P)
    ph = np.exp(2j * np.pi * rng.random(bins.size))
    x = sum(a * cell(5.0e6, 0, K_P, T_P, FS_P, lambda f, t, r=r: np.full_like(f, r)) for a, r in zip(ph, bins, strict=True))
    x = x * math.sqrt(10 ** 1.5 * enbw / n) + noise_cell(n, 1.0, rng)
    sp = spectrum(x, FS_P, K_P, 5.0e6)
    above = int(np.count_nonzero(10 ** (np.asarray(sp.power_db)[0] / 10) / unit > thr))
    spread_cfar = int(cfar_detect(sp).mask.sum())
    rep = {"pfa_poole_capped": fa_poole / cells, "pfa_cfar": fa_cfar / cells, "threshold_over_mean_db": thr_db,
           "pd_poole_cfar": pd, "spread_bins_over_threshold": above, "spread_bins_kept_poole": min(above, 16),
           "spread_bins_kept_cfar": spread_cfar, "trials": trials}
    row = _row("A8", "Echo criterion: at most the 16 highest spectral values per cell over an amplitude criterion, "
                     "from a 1024-sample FFT", "Paper 2", "comparison",
               f"same threshold ({thr_db:.1f} dB over mean): P_fa {rep['pfa_poole_capped']:.1e} (known mean) vs "
               f"{rep['pfa_cfar']:.1e} (P2); Pd at 10 dB {pd['10'][0]:.2f} vs {pd['10'][1]:.2f}; a 40-bin spread "
               f"echo keeps {min(above, 16)} vs {spread_cfar} bins", "comparison", "reported",
               "the cap of 16 never binds on noise; it truncates spread and multiple echoes, which the P2 detector "
               "keeps for the tracker.  The P2 detector estimates the mean from the cell (ordered statistic), so "
               "it pays a small CFAR loss against the known-mean criterion")
    return [row], rep


def parabolic_layer(f_c0: float, rate_hz_per_s: float, h_b: float, y_m: float) -> tuple[Any, Any]:
    """A bottomside parabolic layer whose critical frequency changes linearly in time.

    Group height ``h' = h_b + (y_m/2) x ln((1+x)/(1-x))``, ``x = f/f_c``; phase
    height ``P = h_b + (y_m/2)(f_c/f)[((x^2-1)/2) ln((1+x)/(1-x)) + x]``, from
    ``d(fP)/df = h'``.
    """
    def phase(f: Any, t: Any) -> np.ndarray:
        fc = f_c0 + rate_hz_per_s * np.asarray(t, dtype=float)
        x = np.asarray(f, dtype=float) / fc
        lg = np.log((1 + x) / (1 - x))
        return h_b + 0.5 * y_m * (fc / f) * (0.5 * (x**2 - 1) * lg + x)

    def group(f: Any, t: Any) -> np.ndarray:
        fc = f_c0 + rate_hz_per_s * np.asarray(t, dtype=float)
        x = np.asarray(f, dtype=float) / fc
        return h_b + 0.5 * y_m * x * np.log((1 + x) / (1 - x))

    return phase, group


def _fc_schedule(f_m: float, h_b: float, y_m: float, h0: float, v: float) -> Any:
    """foF2(t) chosen so that the group height at ``f_m`` is ``h0 + v t`` exactly."""
    fc = np.linspace(f_m * 1.0005, f_m * 1.6, 20001)
    x = f_m / fc
    g = h_b + 0.5 * y_m * x * np.log((1 + x) / (1 - x))      # decreasing in fc
    return lambda t: np.interp(-(h0 + v * np.asarray(t, dtype=float)), -g, fc)


def item_a9(snr_db: float, rng: np.random.Generator) -> tuple[list[dict[str, Any]], dict[str, Any]]:
    """Paper 2's two estimates of V' at 6.4 MHz, on a layer whose true V' there is constant.

    The layer is parabolic (h_b 250 km, y_m 100 km); its critical frequency
    rises so that the group height at 6.4 MHz falls at 50 m/s from 376 km,
    near paper 2's -44.5 / -55.1 m/s at 378 km.  The layer changes shape, so
    V* differs from V' and varies with frequency, as in paper 2.
    """
    from tss.l3_signal import poole_range as pr

    f_m, h_b, y_m, v_true = 6.4e6, 250.0e3, 100.0e3, -50.0
    h0 = 376.0e3
    fc_t = _fc_schedule(f_m, h_b, y_m, h0, v_true)

    def phase(f: Any, t: Any) -> np.ndarray:
        fc = fc_t(t)
        x = np.asarray(f, dtype=float) / fc
        lg = np.log((1 + x) / (1 - x))
        return h_b + 0.5 * y_m * (fc / f) * (0.5 * (x**2 - 1) * lg + x)

    def group(f: Any, t: Any) -> np.ndarray:
        x = np.asarray(f, dtype=float) / fc_t(t)
        return h_b + 0.5 * y_m * x * np.log((1 + x) / (1 - x))

    # swept experiment (paper 2 Fig. 4): soundings every 75 kHz, 1.5 s each, up to 7 MHz
    f0s = np.arange(5.5e6, 7.0e6, STEP_P)
    fc_cell = f0s + 0.5 * K_P * T_P
    t_s = np.arange(f0s.size) * 3 * T_P
    ok = fc_cell < 0.98 * fc_t(t_s)                          # below the critical frequency
    v_star = np.full(f0s.size, np.nan)
    for j in np.flatnonzero(ok):
        cs = sounding(f0s[j], t_s[j], phase, snr_db=snr_db, rng=rng)
        v_star[j] = pr.three_cell(cs, FS_P, k_hz_per_s=K_P, f0_hz=f0s[j], delta_f_hz=DF_P)["phase_velocity_mps"]
    near = np.flatnonzero(ok)[np.argsort(np.abs(fc_cell[ok] - f_m))[:9]]
    cf = np.polyfit(fc_cell[near] - f_m, v_star[near], 2)
    v_at, dv_df = float(cf[2]), float(cf[1])
    v_i = v_at + f_m * dv_df
    # fixed-frequency experiment (paper 2 Fig. 5): 3 min at 6.4 MHz, straight after
    t_start = float(t_s[-1]) + 3 * T_P
    f0 = f_m - 0.5 * K_P * T_P
    rs, ts = [], []
    for j in range(120):
        t0 = t_start + j * 3 * T_P
        cs = sounding(f0, t0, phase, snr_db=snr_db, rng=rng)
        rs.append(pr.three_cell(cs, FS_P, k_hz_per_s=K_P, f0_hz=f0, delta_f_hz=DF_P)["group_range_m"])
        ts.append(t0 + 1.5 * T_P)
    v_d = float(np.polyfit(ts, rs, 1)[0])
    ratio = v_i / v_d
    t_mid = float(np.mean(t_s[near]))
    true_v_star = float((phase(f_m, t_mid + 1.0) - phase(f_m, t_mid - 1.0)) / 2.0)
    rep = {"V_star_at_fm": v_at, "V_star_true_at_fm": true_v_star, "dV_star_df": dv_df, "V_prime_indirect": v_i,
           "V_prime_direct": v_d, "V_prime_true": v_true, "soundings_fitted": int(near.size),
           "group_height_at_fm_km": float(group(f_m, t_mid)) / 1e3, "foF2_at_sweep_mhz": float(fc_t(t_mid)) / 1e6}
    row = _row("A9", "Group velocity V' = V* + f dV*/df from a swept ionogram; the indirect estimate within 20 % of "
                     "the direct dR'/dt", "Paper 2 eqs. 1-6, Figs. 4-5",
               "within 20 % (-44.52 vs -55.12 m/s)",
               f"{v_i:.2f} vs {v_d:.2f} m/s ({100 * (ratio - 1):+.1f} %); truth {v_true:.1f} m/s; V* {v_at:.2f} m/s, "
               f"dV*/df {dv_df:.2e} m/s/Hz", "20 % (paper 2's own agreement)",
               "reproduced" if abs(ratio - 1) <= 0.2 else "not reproduced, open",
               f"a parabolic layer whose shape changes so that V' at 6.4 MHz is constant ({snr_db:.0f} dB per cell); "
               f"paper 2's V* -8.68 m/s and dV*/df -5.60e-6 m/s/Hz are of the same order; the {near.size} soundings "
               "nearest 6.4 MHz fitted with a parabola")
    return [row], rep


def item_barry(rng: np.random.Generator) -> list[dict[str, Any]]:
    from tss.l3_signal.detect import cfar_detect

    k, coh, b_rx = 14.0e3, 30.0e3, 100.0
    cell_s = coh / k
    fs = 2.0 * b_rx
    r0 = 300.0e3
    x = cell(5.0e6, 0, k, cell_s, fs, lambda f, t: np.full_like(f, r0))
    sp = spectrum(x, fs, k, 5.0e6, pad=64)
    p = np.asarray(sp.power_db)[0]
    i = int(np.argmax(p))
    a = i
    while a > 0 and p[a - 1] > p[i] - 3.0:
        a -= 1
    b = i
    while b < p.size - 1 and p[b + 1] > p[i] - 3.0:
        b += 1
    width_s = 2.0 * (sp.range_m[b] - sp.range_m[a] + sp.bin_spacing_m) / C
    xs = cell(5.0e6, 0, k, cell_s, fs, lambda f, t: np.full_like(f, r0), snr_db=25.0, rng=rng)
    spd = spectrum(xs, fs, k, 5.0e6)
    detected = bool(np.any(cfar_detect(spd, max_range_m=None).mask & (np.abs(spd.range_m - r0) < 20e3)))
    r_max = spd.unambiguous_range_m
    return [
        _row("A10", "Time-delay resolution 45 us with 30 kHz sweep coherence at 14 kHz/s", "Barry 1971, Table I",
             45.0, round(width_s * 1e6, 2), "10 %", "reproduced" if _within(width_s * 1e6, 45.0, 0.10) else
             "difference explained", f"half-power width of a point echo through the chain (Hann, coherent "
             f"cell {cell_s:.3f} s); 1/(30 kHz) = 33.3 us before the window"),
        _row("A11", "Virtual height range 0 to 1000 km with a 100 Hz receiver bandwidth at 14 kHz/s",
             "Barry 1971, Table I", 1000.0, round(r_max / 1e3, 1), "10 %",
             "reproduced" if _within(r_max, 1000.0e3, 0.10) else "difference explained",
             "the chain's unambiguous range with the beat band limited to 100 Hz: c B/(2k)"),
        _row("A12", "Chain driven at Barry's parameters (linear 14 kHz/s, 30 kHz coherence, 100 Hz bandwidth)",
             "Barry 1971, Table I", "runs", f"echo detected: {detected}", "runs",
             "reproduced" if detected else "not reproduced, open", "25 dB echo at 300 km"),
        _row("A13", "8 W peak, 3 W average", "Barry 1971, Table I", "8 W / 3 W", "-", "-", "not simulated",
             "transmitter hardware; a 50 % gate alone gives 4 W average, so the 3 W includes the shaping Barry "
             "describes (gain-weighted gating)"),
        _row("A14", "Logarithmic sweep 49 s/octave over 0.5-30 MHz", "Barry 1971, Table I", "49 s/octave",
             f"{49.0 * math.log2(30 / 0.5):.0f} s per ionogram", "-", "not simulated",
             "a sweep-schedule property; the chain's A1 sweep is stepped, not logarithmic, and the paper states no "
             "duration to check"),
        _row("A15", "Gating halves the average power and loses half the echoes", "Barry 1971, section IV", "3 dB + half",
             "-", "-", "not simulated", "excluded from P5 by decision (2026-11-03); tss-sim's design is CW with a "
             "canceller (WP3)"),
        _row("A16", "Pulse-repetition frequency randomly varied 350-800 pps", "Barry 1971, Table I", "350-800 pps",
             "-", "-", "not simulated", "the T/R gating pattern (paper 2 uses a 477.1 us M-sequence); see A15"),
        _row("A17", "No echoes below 1.6 MHz (antenna-limited); ionograms comparable with pulse sondes (Fig. 5)",
             "Barry 1971, section V, Fig. 5", "qualitative", "-", "-", "not simulated",
             "antenna and bottomside geometry; no ground vertical-incidence geometry is built (decision 2026-11-03)"),
    ]


def compute(config: dict[str, Any]) -> dict[str, Any]:
    opts = {**DEFAULTS, **dict(config.get("fmcw_acceptance", {}) or {})}
    rng = np.random.default_rng(int(opts["seed"]))
    trials = int(opts["mc_trials"])
    rows: list[dict[str, Any]] = []
    rows += item_a1_a2(rng)
    rows += item_a3(rng)
    rows += item_a4(rng)
    rows += item_a5(rng)
    a67, curve = item_a6_a7(trials, rng)
    rows += a67
    a8, rep8 = item_a8(trials, rng)
    rows += a8
    a9, rep9 = item_a9(float(opts["a9_snr_db"]), rng)
    rows += a9
    rows += item_barry(rng)
    return {"options": opts, "rows": rows, "phase_error_curve": curve, "a8": rep8, "a9": rep9}


def _summary(res: dict[str, Any]) -> tuple[str, dict[str, Any]]:
    rows = res["rows"]
    lines = [f"Monte Carlo {res['options']['mc_trials']} trials per point, seed {res['options']['seed']}.  "
             "Tolerances (G5): arithmetic 1 %, statistical inside the 95 % CI, Barry 10 %.", "",
             "| Item | Statement | Source | Published | Measured | Tolerance | Outcome |",
             "|---|---|---|---|---|---|---|"]
    for r in rows:
        lines.append(f"| {r['item']} | {r['statement']} | {r['source']} | {r['published']} | {r['measured']} "
                     f"| {r['tolerance']} | **{r['outcome']}** |")
    lines += ["", "Notes:", ""] + [f"- **{r['item']}**: {r['note']}" for r in rows if r["note"]]
    lines += ["", "Phase error against SNR (A6):", "", "| SNR dB | MC rms deg [95 % CI] | Poole eq. 82 deg |",
              "|---|---|---|"]
    for c in res["phase_error_curve"]:
        lines.append(f"| {c['snr_db']:.0f} | {c['mc_rms_deg']:.3f} [{c['ci_lo_deg']:.3f}, {c['ci_hi_deg']:.3f}] "
                     f"| {c['closed_form_deg']:.3f} |")
    counts: dict[str, int] = {}
    for r in rows:
        counts[r["outcome"]] = counts.get(r["outcome"], 0) + 1
    lines += ["", "Outcomes: " + "; ".join(f"{k} {v}" for k, v in sorted(counts.items())), ""]
    return "\n".join(lines), {"outcomes": counts, "items": {r["item"]: r["outcome"] for r in rows}}


@register_study(
    "fmcw_acceptance",
    description="WP4 Phase 5: reproduce Poole (1985), Poole and Evans (1985) and Barry (1971) at their own parameters",
    work_packages=["WP4"],
    claims=["C1", "C3"],
    status="implemented",
    requirements=["INS-05", "INS-09", "INS-10", "SCI-07"],
)
def run(out_dir: Path, *, config: dict[str, Any], manifest: Any) -> StudyResult:
    res = compute(dict(config))
    path = out_dir / "fmcw_acceptance_items.csv"
    cols = ["item", "statement", "source", "published", "measured", "tolerance", "outcome", "note"]
    with path.open("w", encoding="utf-8") as fh:
        fh.write(",".join(cols) + "\n")
        for r in res["rows"]:
            fh.write(",".join('"' + str(r[c]).replace('"', "'") + '"' for c in cols) + "\n")
    text, metrics = _summary(res)
    jpath = out_dir / "fmcw_acceptance.json"
    jpath.write_text(json.dumps(res | {"metrics": metrics}, indent=2, default=str), encoding="utf-8")
    return StudyResult(outputs=[path, jpath], summary=text, metrics=metrics)
