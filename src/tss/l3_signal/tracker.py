"""Layer: L3
WP: WP4
Claims: C1, C2
Status: implemented
Requirements: SCI-05, SCI-07, INS-05, INS-08, INS-10

The trace tracker (WP4 Phase 4): labelled O and X traces from the column stack.

Decisions of 2026-10-31 (Andrew): track-before-detect, with the O and X traces
labelled by self-consistency; only what is available on board may be used (the
spacecraft position, so the IGRF field, and the ionogram itself); the policy for
blended O and X samples (gate G4) is built as options and chosen on data; the
rate search keeps a rate per detection (:func:`tss.l3_signal.ionogram.assign_rates`).

WHAT IT DOES
------------
1. **Onboard priors.**  The gyrofrequency at the spacecraft comes from the
   geomagnetic field at its position (:func:`tss.l1_propagation.lamination.geomagnetic_field_model`,
   IGRF with a dipole fallback).  The local plasma frequency is estimated from
   the ionogram, as the Alouette and ISIS analysts did: the plasma resonances
   ring through a whole column, so the lowest streak column that is not a
   gyroharmonic is ``f_N`` (``f_UH`` shares its step when ``f_N >> f_H``).
   With no resonance (clutter off) the lowest column with a near-range echo
   is used.  The X trace starts at ``f_x = f_H/2 + sqrt(f_H^2/4 + f_N^2)``.
2. **Masks.**  Streak columns (resonances, broadcast carriers) are left out of
   the tracking.  The ground echo, which every topside ionogram above foF2
   carries at about the spacecraft height and which took the WP3-era extractor
   (P1), is tracked from the top of the band down, where its range rises
   towards the cusp, and every cell at or beyond it is masked.
3. **Track-before-detect.**  A Viterbi search on a common range grid finds the
   monotone path (range non-decreasing with frequency, up to a slope limit,
   with gaps) that maximises the summed log-likelihood ratio of a 10 dB echo
   against noise.  Weak samples on a strong trace are kept; isolated false
   alarms (about 12 per ionogram at 1e-4) are not on any trace.
4. **O and X by self-consistency.**  At a given frequency the X echo is always
   nearer than the O echo, so the two paths are found with that ordering, each
   started near its own cut-off at the spacecraft.  The O path is inverted
   (magnetoionic lamination, the field from step 1) and the X trace that
   profile implies is predicted (:func:`predict_x_range_m`); the X path is then
   gated on the prediction and the O path kept out of it, and the two are
   iterated.  The X residual against the prediction is the consistency
   statistic, and each sample's distance from the predicted X trace gives its
   labelling confidence.
5. **Samples.**  Each path point is refined to the leading edge of the stated
   P_fa detector's run at that range (or kept as a track-only sample, flagged,
   when the echo is under the cell threshold).  Each gets a range uncertainty
   and quality flags: track-only, blended with the other mode, spread, beyond
   the rate-search span (``N > 64``, P3), cusp outlier, rejected.
6. **Cusp model.**  The top of the O trace is fitted with the virtual range of
   a parabolic layer near its critical frequency, ``a + b x ln((1+x)/(1-x))/2``,
   ``x = f / foF2``; foF2 is reported and samples far off the fit are flagged
   (``cusp_policy="flag"``) or replaced by it (``"replace"``).
7. **Gate G4** (``blend_policy``): where the O and X samples of a column are
   within ``blend_cells`` range cells, ``"reject"`` drops the O sample,
   ``"inflate"`` keeps it with its uncertainty raised to the separation, and
   ``"separate"`` fits two echoes to the column (after Paul 1983).

The tracker never reads the truth.  Everything it uses is in the received
record, the sweep and the platform state.
"""

from __future__ import annotations

import math
from collections.abc import Callable, Sequence
from dataclasses import dataclass, field, replace
from typing import Any

import numpy as np

__all__ = [
    "TrackerSettings",
    "TrackResult",
    "Q_TRACK_ONLY", "Q_BLENDED", "Q_SPREAD", "Q_CUSP", "Q_OVER_SPAN", "Q_LOW_CONFIDENCE", "Q_REJECTED",
    "Q_STOPPED", "stop_at_loss",
    "track_traces",
    "predict_x_range_m",
    "viterbi_path",
    "column_image",
    "streak_columns",
    "onboard_priors",
    "cusp_fit",
]

C_LIGHT = 299_792_458.0

#: Quality bits per sample.
Q_TRACK_ONLY = 1      # on the track, but no cell over the detector threshold
Q_BLENDED = 2         # within blend_cells of the other mode's sample
Q_SPREAD = 4          # the echo spreads over more than spread_cells
Q_CUSP = 8            # far off the cusp model
Q_OVER_SPAN = 16      # compressed at the edge of the rate-search span (N > span)
Q_LOW_CONFIDENCE = 32  # labelling confidence under 0.5
Q_REJECTED = 64       # removed by the G4 policy
Q_STOPPED = 128       # beyond the point where the O trace was lost (stop_at_loss)

FieldModel = Callable[[np.ndarray], tuple[np.ndarray, np.ndarray]]


@dataclass(frozen=True)
class TrackerSettings:
    """Tracker settings.  Ranges in metres."""

    grid_m: float = 2.0e3
    max_range_m: float = 1.05e6
    design_snr_db: float = 10.0
    score_clip_db: float = 30.0
    tbd_bias: float = 6.0
    max_gap: int = 3
    gap_cost: float = 1.0
    late_start_cost: float = 8.0
    down_tol_m: float = 4.0e3
    up_max_m: float = 60.0e3
    up_frac: float = 0.35
    up_cap_m: float = 150.0e3
    start_max_m: float = 80.0e3
    start_steps: int = 4
    min_track_score: float = 12.0
    min_track_columns: int = 3
    streak_snr_db: float = 6.0
    streak_near_m: float = 200.0e3
    streak_near_fraction: float = 0.8
    streak_fraction: float = 0.3
    ground_min_frac_hsat: float = 0.6
    ground_margin_m: float = 30.0e3
    ground_min_columns: int = 4
    order_sep_m: float = 2.0e3
    x_gate_m: float = 5.0e3
    x_gate_frac: float = 0.03
    x_gate_slope: float = 0.5
    iterations: int = 3
    blend_cells: float = 2.0
    blend_policy: str = "reject"
    cusp_policy: str = "flag"
    cusp_min_samples: int = 4
    cusp_sigma_factor: float = 3.0
    spread_cells: float = 8.0
    shadow_m: float = 100.0e3
    shadow_db: float = 12.0
    edge_window_m: float = 6.0e3
    edge_drop_db: float = 15.0
    # WP4 Phase 6 (2026-11-05; defaults set from the development sample).
    # ``fn_prior``: ``near_start`` is P4's rule (f_N
    # from the resonances or the lowest near-range echo a trace follows);
    # ``trace_floor`` also takes the lowest frequency at which any persistent
    # trace starts, at any range, as an upper bound on f_N (no echo returns
    # below the local plasma frequency).  Under a shallow topside gradient the
    # O trace starts hundreds of kilometres out while the X trace starts near
    # range, so P4's rule took f_x for f_N.  ``far_start_cost`` is the start
    # cost of a far-range start in the first columns above the prior.
    fn_prior: str = "trace_floor"
    far_start_cost: float = 2.0
    floor_persist: int = 2
    # ``stop_at_loss`` (off: it raised SCI-05 on ISIS-1 by 1.5-3.5 points but
    # cut O completeness on the Chapman campaign from 0.95 to 0.69 with clutter
    # off): end the O trace where it is lost near the top of the band (see
    # :func:`stop_at_loss`).  ``cusp_up_cap_m``: the largest rise per
    # column allowed in the steep parts of the O trace, the top ``cusp_band``
    # of the O search and the first columns above f_N (cusp tracking); equal
    # to ``up_cap_m`` it is off (400 km changed nothing measurable).
    stop_at_loss: bool = False
    stop_run: int = 2
    cusp_up_cap_m: float = 150.0e3
    cusp_band: float = 0.25

    def __post_init__(self) -> None:
        if self.blend_policy not in ("reject", "inflate", "separate"):
            raise ValueError(f"unknown blend_policy {self.blend_policy!r}; expected reject, inflate or separate")
        if self.cusp_policy not in ("flag", "replace"):
            raise ValueError(f"unknown cusp_policy {self.cusp_policy!r}; expected flag or replace")
        if self.fn_prior not in ("near_start", "trace_floor"):
            raise ValueError(f"unknown fn_prior {self.fn_prior!r}; expected near_start or trace_floor")


@dataclass(frozen=True)
class TrackResult:
    """Labelled traces.  Arrays are (n_freq, 2), column 0 the O mode, 1 the X mode."""

    freq_hz: np.ndarray
    range_m: np.ndarray
    snr_db: np.ndarray
    sigma_m: np.ndarray
    confidence: np.ndarray
    quality: np.ndarray
    info: dict[str, Any] = field(default_factory=dict)


# -- the image ------------------------------------------------------------------------


def column_image(spectra: Sequence[Any], floors_db: Sequence[float], grid_m: float,
                 max_range_m: float) -> np.ndarray:
    """SNR in dB over each column's noise floor, max-pooled onto a common range grid."""
    n_r = int(math.floor(max_range_m / grid_m)) + 1
    img = np.full((len(spectra), n_r), -np.inf)
    for i, (sp, fl) in enumerate(zip(spectra, floors_db, strict=True)):
        rng = np.asarray(sp.range_m, dtype=float)
        p = np.asarray(sp.power_db, dtype=float)[0] - float(fl)
        ok = (rng >= 0.0) & (rng <= max_range_m) & np.isfinite(p)
        idx = np.floor(rng[ok] / grid_m).astype(int)
        np.maximum.at(img[i], idx, p[ok])
    return img


def _flatten_background(img_db: np.ndarray, streak: np.ndarray) -> np.ndarray:
    """Divide out the range dependence of the noise background.

    Each column is already referred to its own noise floor (the negative-beat
    half), but the floor is not flat in range: the receiver's anti-alias skirt
    and the leakage pedestal shape it, by a decibel or two at the far end, and
    the track search, which sums evidence over many columns, turns a small bias
    into a path through noise.  Echoes occupy few columns at any one range, so
    the median over columns at each range cell is the background there; for a
    unit exponential the median is ``ln 2``, which sets the reference.
    """
    keep = ~np.asarray(streak, dtype=bool)
    if np.count_nonzero(keep) < 8:
        return img_db
    z = 10.0 ** (img_db[keep] / 10.0)
    med = np.nanmedian(np.where(np.isfinite(z), z, np.nan), axis=0)
    good = np.isfinite(med) & (med > 0)
    corr = np.zeros(med.shape)
    corr[good] = 10.0 * np.log10(med[good] / math.log(2.0))
    corr = np.clip(corr, -6.0, 6.0)
    return img_db - corr[None, :]


def _shadow_suppress(img_db: np.ndarray, grid_m: float, s: TrackerSettings) -> np.ndarray:
    """Cells far below a stronger cell just nearer in the same column are noise to the tracker.

    Spread F trails an echo by tens of kilometres at well under its level (in
    the clutter model, 12 rows 17 dB below it), and a window's sidelobes do the
    same at 31 dB.  A cell more than ``shadow_db`` below the strongest cell in
    the ``shadow_m`` nearer than it (beyond the main lobe) is set to the mean
    noise level.  An O echo beside an X echo of comparable strength is kept.
    """
    from scipy.ndimage import maximum_filter1d

    w = max(int(round(s.shadow_m / grid_m)), 2)
    guard = 2
    filled = np.where(np.isfinite(img_db), img_db, -np.inf)
    # running maximum over [r - w, r - guard]
    prev = maximum_filter1d(filled, size=w - guard + 1, axis=1, origin=(w - guard) // 2, mode="constant",
                            cval=-np.inf)
    shifted = np.full_like(prev, -np.inf)
    shifted[:, guard:] = prev[:, :-guard]
    out = img_db.copy()
    out[np.isfinite(img_db) & (img_db < shifted - s.shadow_db)] = 0.0
    return out


def streak_columns(img_db: np.ndarray, grid_m: float, s: TrackerSettings) -> np.ndarray:
    """Columns rung through by a resonance or a broadcast carrier.

    A resonance rings from zero range outwards, so its column is hot over
    nearly all of the first ``streak_near_m``; a carrier, de-chirped, fills the
    whole column.  A trace column with spread F behind its echoes is hot in
    patches that do not start at zero range.
    """
    hot = img_db > s.streak_snr_db
    near = max(int(s.streak_near_m / grid_m), 1)
    return (np.mean(hot[:, :near], axis=1) >= s.streak_near_fraction) | (np.mean(hot, axis=1) > s.streak_fraction)


def _step_halfwidths(freq_hz: np.ndarray) -> np.ndarray:
    f = np.asarray(freq_hz, dtype=float)
    if f.size < 2:
        return np.full(f.size, 0.025 * float(f[0]) if f.size else 0.0)
    d = np.diff(f)
    return 0.5 * np.concatenate([[d[0]], 0.5 * (d[:-1] + d[1:]), [d[-1]]])


def onboard_priors(freq_hz: np.ndarray, img_db: np.ndarray, streak: np.ndarray, f_h_hz: float,
                   grid_m: float, s: TrackerSettings) -> dict[str, Any]:
    """Local plasma frequency from the resonances or the trace start; f_x from it and f_H."""
    f = np.asarray(freq_hz, dtype=float)
    half = _step_halfwidths(f)
    gyro = np.zeros(f.size, dtype=bool)
    gyro_wide = np.zeros(f.size, dtype=bool)
    for n in (1, 2, 3, 4):
        gyro |= np.abs(f - n * f_h_hz) <= half + 1.0
        # a harmonic can ring in the next column too: the field model's f_H
        # and the column centres differ by up to a step
        gyro_wide |= np.abs(f - n * f_h_hz) <= 2.0 * half + 1.0
    floor = float("nan")
    if s.fn_prior == "trace_floor":
        # The lowest column at which a persistent trace starts, at any range:
        # a strong cell, and strong cells in at least ``floor_persist`` of the
        # next three columns.  No echo returns below the local plasma
        # frequency, so f_N lies within a step below it (the first sample is
        # at most one step above the cut-off); its midpoint is the estimate.
        # Columns at or next to a gyroharmonic are left out: the resonance
        # there can ring in two or three adjacent columns and look like a
        # trace start.
        strong = np.any(img_db > s.design_snr_db + 2.0, axis=1) & ~streak & ~gyro_wide
        for i in range(f.size - 1):
            if strong[i] and np.count_nonzero(strong[i + 1: i + 4]) >= s.floor_persist:
                floor = float(f[i] - half[i])
                break
    cand = np.flatnonzero(streak & ~gyro)
    if np.isfinite(floor):
        # f_N is at or below the trace floor, and the O trace starts within a
        # step or two above it: the resonance is the highest candidate not
        # above the floor (a lower one is a harmonic missed by the tolerance,
        # or a carrier), stepped down to its partner when it is f_UH.
        k = int(np.searchsorted(f, floor))
        cand = cand[f[cand] <= floor + 2.0 * half[min(k, f.size - 1)]]
        cand = cand[::-1]
    method = None
    f_n = float("nan")
    if cand.size:
        c1 = f[cand[0]]
        f_n, method = float(c1), "resonance"
        if np.isfinite(floor):
            for c in cand[1:]:
                if abs(c1 - math.sqrt(f[c] ** 2 + f_h_hz**2)) <= 2.0 * half[cand[0]]:
                    f_n, method = float(f[c]), "resonance_pair"
                    break
        elif cand.size > 1:
            c2 = f[cand[1]]
            if abs(c2 - math.sqrt(c1**2 + f_h_hz**2)) <= 2.0 * half[cand[1]]:
                method = "resonance_pair"
    if not np.isfinite(f_n):
        # The lowest column with a strong near-range echo that a trace follows
        # (any echo in the next three columns): near f_N the trace rises
        # steeply, so the next point may already be far out.
        near = max(int(s.start_max_m / grid_m), 1)
        hit_near = np.any(img_db[:, :near] > s.design_snr_db + 2.0, axis=1) & ~streak
        hit_any = np.any(img_db > s.design_snr_db, axis=1) & ~streak
        for i in range(f.size - 1):
            if hit_near[i] and np.any(hit_any[i + 1: i + 4]):
                f_n, method = float(f[i] - half[i]), "trace_start"
                break
    near_rule = f_n
    if np.isfinite(floor) and (not np.isfinite(f_n) or floor < f_n - half[min(np.searchsorted(f, f_n), f.size - 1)]):
        f_n, method = floor, "trace_floor"
    f_x = float(0.5 * f_h_hz + math.sqrt(0.25 * f_h_hz**2 + f_n**2)) if np.isfinite(f_n) else float("nan")
    return {"f_n_hz": f_n, "f_h_hz": float(f_h_hz), "f_x_hz": f_x, "f_n_method": method,
            "f_n_near_rule_hz": near_rule, "trace_floor_hz": floor,
            "gyro_columns": [float(v) for v in f[gyro & streak]],
            "streak_columns": [float(v) for v in f[streak]]}


# -- the Viterbi search ---------------------------------------------------------------


def _llr(snr_db: np.ndarray, s: TrackerSettings) -> np.ndarray:
    """Per-cell track evidence ``min(z, z_clip) - tbd_bias``, ``z`` the SNR as a power ratio.

    The search takes, in effect, the best of the few tens of cells its slope
    limit allows in each column, and the largest of K unit exponentials is
    about ``ln K + 0.58``; the bias of 6 (7.8 dB) makes a path through noise
    lose score and one along a 10 dB echo gain 4 per column.  The clip keeps a
    single strong cell from carrying a path.
    """
    z = 10.0 ** (np.minimum(snr_db, s.score_clip_db) / 10.0)
    out = z - s.tbd_bias
    return np.where(np.isfinite(snr_db), out, -np.inf)


def viterbi_path(score: np.ndarray, allowed: np.ndarray, start_cost: np.ndarray, grid_m: float,
                 s: TrackerSettings, up_cap_cols: np.ndarray | None = None) -> tuple[np.ndarray, float]:
    """Best monotone path through ``score`` (columns in frequency order).

    The range of consecutive path points may fall by at most ``down_tol_m`` and
    rise by at most ``max(up_max_m, up_frac * r)`` (capped at ``up_cap_m``);
    up to ``max_gap`` columns may be skipped at ``gap_cost`` each.  A path may
    start at a cell at ``start_cost`` (``inf`` where it may not).  Returns the
    range index per column (-1 where the path has no point) and the score.
    ``up_cap_cols`` replaces ``up_cap_m`` column by column (cusp tracking).
    """
    n_c, n_r = score.shape
    sc = np.where(allowed, score, -np.inf)
    d = int(round(s.down_tol_m / grid_m))
    rr = np.arange(n_r)
    up_base = np.maximum(s.up_max_m, s.up_frac * rr * grid_m)
    caps = (np.full(n_c, s.up_cap_m) if up_cap_cols is None
            else np.asarray(up_cap_cols, dtype=float))
    up_cells_c = np.round(np.minimum(up_base[None, :], caps[:, None]) / grid_m).astype(int)
    u_max = int(up_cells_c.max()) if n_r and n_c else 0
    v = np.full((n_c, n_r), -np.inf)
    back_c = np.full((n_c, n_r), -1, dtype=np.int32)
    back_r = np.full((n_c, n_r), -1, dtype=np.int32)
    carry = np.full(n_r, -np.inf)
    carry_c = np.full(n_r, -1, dtype=np.int32)
    carry_gap = np.zeros(n_r, dtype=np.int32)
    for i in range(n_c):
        up_cells = up_cells_c[i]
        best = -np.inf * np.ones(n_r)
        bc = np.full(n_r, -1, dtype=np.int32)
        br = np.full(n_r, -1, dtype=np.int32)
        if np.any(np.isfinite(carry)):
            for k in range(-d, u_max + 1):
                # previous range rp = r - k: k > 0 is a rise, k < 0 a fall
                if k >= 0:
                    src = np.concatenate([np.full(k, -np.inf), carry[: n_r - k]]) if k else carry
                    src_r = rr - k
                else:
                    src = np.concatenate([carry[-k:], np.full(-k, -np.inf)])
                    src_r = rr - k
                ok = (k <= up_cells) & (src > best)
                if np.any(ok):
                    best = np.where(ok, src, best)
                    srcc = np.clip(src_r, 0, n_r - 1)
                    bc = np.where(ok, carry_c[srcc], bc)
                    br = np.where(ok, srcc, br)
        st = -np.asarray(start_cost[i], dtype=float)
        use_start = st > best
        base = np.where(use_start, st, best)
        v[i] = sc[i] + base
        back_c[i] = np.where(use_start, -1, bc)
        back_r[i] = np.where(use_start, -1, br)
        # carry: the best path ending at each range in the last max_gap + 1 columns
        decayed = carry - s.gap_cost
        gap = carry_gap + 1
        take = (v[i] >= decayed) | (gap > s.max_gap)
        carry = np.where(take, v[i], decayed)
        carry_c = np.where(take, i, carry_c).astype(np.int32)
        carry_gap = np.where(take, 0, gap).astype(np.int32)
    path = np.full(n_c, -1, dtype=int)
    if not np.any(np.isfinite(v)):
        return path, float("-inf")
    i, r = np.unravel_index(int(np.argmax(np.where(np.isfinite(v), v, -np.inf))), v.shape)
    while i >= 0 and r >= 0:
        path[i] = r
        i, r = int(back_c[i, r]), int(back_r[i, r])
    return _trim(path, sc, s)


def _trim(path: np.ndarray, sc: np.ndarray, s: TrackerSettings) -> tuple[np.ndarray, float]:
    """Keep the run of path points with the largest summed evidence (gaps cost ``gap_cost``).

    The search's end points are free, so a path can carry a tail through noise
    that a few strong cells pay for; the maximum-sum run drops it.
    """
    cols = np.flatnonzero(path >= 0)
    if cols.size == 0:
        return path, float("-inf")
    vals = np.array([sc[c, path[c]] for c in cols], dtype=float)
    gaps = np.concatenate([[0], np.diff(cols) - 1]) * s.gap_cost
    best, best_lo, best_hi = -np.inf, 0, 0
    run, lo = 0.0, 0
    for k in range(cols.size):
        add = vals[k] - (gaps[k] if k > lo else 0.0)
        if run <= 0.0 or k == lo:
            if run <= 0.0:
                run, lo = vals[k], k
            else:
                run += add
        else:
            run += add
        if run > best:
            best, best_lo, best_hi = run, lo, k
    out = np.full(path.size, -1, dtype=int)
    keep = cols[best_lo: best_hi + 1]
    out[keep] = path[keep]
    return out, float(best)


# -- the magnetoionic forward model for X from an O inversion ----------------------------


def _profile_from_o(freq_hz: np.ndarray, range_m: np.ndarray, f_start_hz: float, h_sat_m: float,
                    fld: FieldModel) -> tuple[np.ndarray, np.ndarray] | None:
    from tss.l1_propagation.lamination import lamina_thicknesses_magnetoionic_m

    f = np.asarray(freq_hz, dtype=float)
    r = np.asarray(range_m, dtype=float)
    ok = np.isfinite(f) & np.isfinite(r) & (f > f_start_hz)
    if np.count_nonzero(ok) < 2:
        return None
    order = np.argsort(f[ok])
    f, r = f[ok][order], np.maximum.accumulate(np.maximum(r[ok][order], 0.0))
    try:
        t = lamina_thicknesses_magnetoionic_m(r, f, float(f_start_hz), float(h_sat_m), fld)
    except (ValueError, FloatingPointError):
        return None
    t = np.where(np.isfinite(t), np.maximum(t, 0.0), 0.0)
    h = float(h_sat_m) - np.concatenate([[0.0], np.cumsum(t)])
    fp = np.concatenate([[float(f_start_hz)], f])
    return h, fp


def predict_x_range_m(freq_o_hz: np.ndarray, range_o_m: np.ndarray, f_start_hz: float, h_sat_m: float,
                      fld: FieldModel, f_query_hz: np.ndarray, *, n_quad: int = 48) -> np.ndarray:
    """X-mode virtual range implied by an O trace, at ``f_query_hz``.

    The O trace is inverted by the magnetoionic lamination (plasma frequency
    linear in height in each lamina).  For each query frequency the X
    reflection height solves ``f_p(h)^2 = f (f - f_H(h))`` and the virtual
    range is the integral of the X group index (Budden 12.8, through
    :func:`tss.l1_propagation.magnetoionic.group_index`) from there to the
    spacecraft, with ``h = h_r + s^2`` taking out the square-root singularity
    at reflection.  NaN where the X wave reflects outside the inverted part
    of the profile.
    """
    from tss.l1_propagation.magnetoionic import group_index

    fq = np.atleast_1d(np.asarray(f_query_hz, dtype=float))
    out = np.full(fq.shape, np.nan)
    prof = _profile_from_o(freq_o_hz, range_o_m, f_start_hz, h_sat_m, fld)
    if prof is None:
        return out
    h_e, fp_e = prof
    h_asc, fp_asc = h_e[::-1], fp_e[::-1]
    lo, hi = float(h_e[-1]), float(h_sat_m)
    xg, wg = np.polynomial.legendre.leggauss(int(n_quad))
    for k, f in enumerate(fq):
        if not np.isfinite(f) or f <= 0:
            continue

        def g(h: np.ndarray, f: float = float(f)) -> np.ndarray:
            fp = np.interp(h, h_asc, fp_asc)
            return fp**2 - f * (f - fld(np.atleast_1d(h))[0])

        if g(np.array([hi]))[0] >= 0.0 or g(np.array([lo]))[0] < 0.0:
            continue
        a, b = lo, hi
        for _ in range(48):
            m = 0.5 * (a + b)
            if g(np.array([m]))[0] >= 0.0:
                a = m
            else:
                b = m
        h_r = b
        smax = math.sqrt(max(hi - h_r, 0.0))
        if smax <= 0.0:
            out[k] = 0.0
            continue
        sv = 0.5 * smax * (xg + 1.0)
        hv = h_r + sv**2
        fp = np.interp(hv, h_asc, fp_asc)
        f_h, th = fld(hv)
        mu = group_index((fp / f) ** 2, np.asarray(f_h) / f, np.asarray(th), mode="X")
        val = np.where(np.isfinite(mu), mu, 0.0) * 2.0 * sv
        out[k] = float(0.5 * smax * np.sum(wg * val))
    return out


# -- the cusp model ---------------------------------------------------------------------


def _cusp_shape(x: np.ndarray) -> np.ndarray:
    x = np.clip(x, 0.0, 1.0 - 1e-9)
    return 0.5 * x * np.log((1.0 + x) / (1.0 - x))


def cusp_fit(freq_hz: np.ndarray, range_m: np.ndarray, *, fo_hint_hz: float | None = None,
             n_top: int = 8, max_trim: int = 3) -> dict[str, Any] | None:
    """Fit the cusp, dropping up to ``max_trim`` trailing points that do not belong to it.

    An O trace rises without bound into foF2; points at the top end that
    flatten instead (an X sample, a ground-echo cell, spread) cannot be on it.
    The fit with ``k`` trailing points left out is preferred to the one with
    ``k - 1`` when its residual is under half as large.  ``n_dropped`` says
    how many were left out.
    """
    f = np.asarray(freq_hz, dtype=float)
    r = np.asarray(range_m, dtype=float)
    ok = np.isfinite(f) & np.isfinite(r)
    order = np.argsort(f[ok])
    f, r = f[ok][order], r[ok][order]
    best = None
    for k in range(0, max_trim + 1):
        if f.size - k < 4:
            break
        fit = _cusp_fit_once(f[: f.size - k], r[: r.size - k], fo_hint_hz=fo_hint_hz, n_top=n_top)
        if fit is None:
            continue
        fit["n_dropped"] = k
        # a trailing point is dropped only when the residual it causes is
        # material (over a kilometre) and halves without it
        if best is None or (best["rms_m"] > 1.0e3 and fit["rms_m"] < 0.5 * best["rms_m"]):
            best = fit
    return best


def _cusp_fit_once(freq_hz: np.ndarray, range_m: np.ndarray, *, fo_hint_hz: float | None,
                   n_top: int) -> dict[str, Any] | None:
    f = np.asarray(freq_hz, dtype=float)
    r = np.asarray(range_m, dtype=float)
    ok = np.isfinite(f) & np.isfinite(r)
    f, r = f[ok], r[ok]
    if f.size < 4:
        return None
    order = np.argsort(f)
    f, r = f[order][-n_top:], r[order][-n_top:]
    fmax = float(f[-1])
    lo = fmax * 1.0005
    hi = fmax * 1.25 if fo_hint_hz is None or not np.isfinite(fo_hint_hz) else max(float(fo_hint_hz) * 1.02, lo * 1.001)
    best = None
    for fc in np.linspace(lo, hi, 200):
        gx = _cusp_shape(f / fc)
        a_mat = np.column_stack([np.ones_like(gx), gx])
        coef, *_ = np.linalg.lstsq(a_mat, r, rcond=None)
        if coef[1] <= 0:
            continue
        res = r - a_mat @ coef
        sse = float(np.sum(res**2))
        if best is None or sse < best[0]:
            best = (sse, fc, coef, res)
    if best is None:
        return None
    sse, fc, coef, res = best
    dof = max(f.size - 3, 1)
    return {"fo_hz": float(fc), "a_m": float(coef[0]), "b_m": float(coef[1]),
            "rms_m": float(math.sqrt(sse / dof)), "n": int(f.size)}


def _cusp_eval(fit: dict[str, Any], f: np.ndarray) -> np.ndarray:
    return fit["a_m"] + fit["b_m"] * _cusp_shape(np.asarray(f, dtype=float) / fit["fo_hz"])


# -- samples ---------------------------------------------------------------------------------


def _runs(mask: np.ndarray) -> list[tuple[int, int]]:
    m = np.concatenate([[False], np.asarray(mask, dtype=bool), [False]])
    d = np.diff(m.astype(int))
    return list(zip(np.flatnonzero(d == 1), np.flatnonzero(d == -1), strict=True))


def _refine(sp: Any, cfar: Any, r_dp: float, s: TrackerSettings, span_rate: float | None,
            exclude: tuple[float, float] | None = None) -> dict[str, Any]:
    """Leading edge of the detector run at the path range, or a track-only sample."""
    rng = np.asarray(sp.range_m, dtype=float)
    pw = np.asarray(sp.power_db, dtype=float)[0]
    floor = float(cfar.noise_floor_db)
    cell = float(sp.range_resolution_m)
    lo, hi = r_dp - s.edge_window_m, r_dp + s.grid_m + s.edge_window_m
    q = 0
    chosen = None
    for a, b in _runs(np.asarray(cfar.mask, dtype=bool)):
        if b - a < 2:
            continue
        if rng[b - 1] >= lo and rng[a] <= hi:
            if exclude is not None and exclude[0] <= rng[a] <= exclude[1]:
                continue
            chosen = (a, b)
            break
    if chosen is not None:
        a, b = chosen
        # the peak nearest the path range within the run, not the run's maximum:
        # a run can hold the other mode's echo and the spread between them
        near = np.flatnonzero((rng[a:b] >= lo) & (rng[a:b] <= hi)) + a
        k = int(near[np.argmax(pw[near])]) if near.size else a + int(np.argmax(pw[a:b]))
        # If the run is held below threshold from a nearer echo (the other
        # mode's, or the spread trailing it), the edge is just after the valley
        # in front of this echo's peak, when there is a real one (6 dB).
        valley = int(a + np.argmin(pw[a:k + 1]))
        start = a
        if valley > a and pw[k] - pw[valley] >= 6.0:
            start = valley + 1
            q |= Q_BLENDED
        # The edge is the first cell from there that is within edge_drop_db of
        # this echo's peak: for an echo under about 25 dB SNR that is the
        # detector's first cell, as before; spread trailing a nearer echo, 17 dB
        # down, no longer extends it.
        e = start
        while e < k and pw[e] < pw[k] - s.edge_drop_db:
            e += 1
        edge = float(rng[e])
        snr = float(pw[k] - floor)
        n_cells = (rng[b - 1] - rng[a]) / max(cell, 1.0) + 1.0
        if n_cells > s.spread_cells:
            q |= Q_SPREAD
    else:
        win = (rng >= lo) & (rng <= hi)
        if not np.any(win):
            return {"range_m": float(r_dp), "snr_db": float("nan"), "q": Q_TRACK_ONLY, "cell_m": cell, "bin": -1}
        idx = np.flatnonzero(win)
        k = int(idx[np.argmax(pw[idx])])
        half = pw[k] - 6.0
        a = k
        while a - 1 >= idx[0] and pw[a - 1] > half:
            a -= 1
        edge = float(rng[a])
        snr = float(pw[k] - floor)
        q |= Q_TRACK_ONLY
    rate = getattr(sp, "rate_hz_per_s2", None)
    if rate is not None and span_rate is not None and np.size(rate) == rng.size:
        if abs(float(np.asarray(rate)[k])) >= span_rate * (1.0 - 1e-9):
            q |= Q_OVER_SPAN
    return {"range_m": edge, "snr_db": snr, "q": q, "cell_m": cell, "bin": k}


def _sigma(snr_db: float, cell_m: float) -> float:
    z = 10.0 ** (snr_db / 10.0) if np.isfinite(snr_db) else 1.0
    return float(math.sqrt((cell_m / (2.0 * math.sqrt(max(z, 1.0)))) ** 2 + cell_m**2 / 12.0))


def _separate(sp: Any, r_near: float, r_far: float) -> tuple[float, float] | None:
    """Two-echo least-squares fit to the column's power (after Paul 1983)."""
    rng = np.asarray(sp.range_m, dtype=float)
    p = 10.0 ** (np.asarray(sp.power_db, dtype=float)[0] / 10.0)
    cell = float(sp.range_resolution_m)
    win = (rng >= r_near - 3 * cell) & (rng <= r_far + 3 * cell)
    if np.count_nonzero(win) < 6:
        return None
    x, y = rng[win], p[win]

    def kern(c: float) -> np.ndarray:
        u = (x - c) / cell
        return np.where(np.abs(u) < 1.0, np.cos(0.5 * np.pi * u) ** 4, 0.0)

    best = None
    for c1 in np.linspace(r_near - cell, r_near + cell, 21):
        k1 = kern(c1)
        for c2 in np.linspace(r_far - cell, r_far + cell, 21):
            if c2 <= c1:
                continue
            a = np.column_stack([k1, kern(c2)])
            coef, *_ = np.linalg.lstsq(a, y, rcond=None)
            if np.any(coef <= 0):
                continue
            sse = float(np.sum((y - a @ coef) ** 2))
            if best is None or sse < best[0]:
                best = (sse, c1, c2)
    if best is None:
        return None
    half = 0.5 * cell
    return best[1] - half, best[2] - half


# -- the tracker -----------------------------------------------------------------------------


def stop_at_loss(freq_hz: np.ndarray, range_m: np.ndarray, quality: np.ndarray, confidence: np.ndarray,
                 s: TrackerSettings) -> float:
    """The frequency from which the O trace is treated as lost; ``nan`` if it is not.

    Two signs, each judged on the samples in frequency order after the first
    three: a run of ``stop_run`` samples (or the rest of the trace, if shorter)
    whose labelling confidence is under 0.5, meaning they sit on the X trace
    the inverted O trace predicts, so the O path has moved onto X; and a tail
    of track-only samples (no cell over the detector threshold).  Samples
    from that frequency up are not inverted: the profile ends higher, which
    shows as coverage, rather than continuing on the wrong echo.
    """
    ok = np.isfinite(range_m) & ((quality & Q_REJECTED) == 0)
    idx = np.flatnonzero(ok)
    idx = idx[np.argsort(freq_hz[idx])]
    if idx.size < 4:
        return float("nan")
    low = np.asarray([np.isfinite(confidence[j]) and confidence[j] < 0.5 for j in idx])
    stop = None
    for k in range(3, idx.size):
        run = low[k: k + s.stop_run]
        if low[k] and np.all(run):
            stop = k
            break
    track_only = (quality[idx] & Q_TRACK_ONLY) != 0
    k = idx.size
    while k > 3 and track_only[k - 1]:
        k -= 1
    if k < idx.size and (stop is None or k < stop):
        stop = k
    return float(freq_hz[idx[stop]]) if stop is not None else float("nan")


def track_traces(spectra: Sequence[Any], cfars: Sequence[Any], *, h_sat_m: float, fld: FieldModel,
                 settings: TrackerSettings | None = None, span_rate_hz_per_s2: Sequence[float] | None = None,
                 ) -> TrackResult:
    """Labelled O and X traces from a stack of range spectra (module docstring)."""
    s = settings or TrackerSettings()
    n = len(spectra)
    freq = np.array([float(sp.centre_frequency_hz) for sp in spectra])
    floors = [float(c.noise_floor_db) for c in cfars]
    img = column_image(spectra, floors, s.grid_m, s.max_range_m)
    n_r = img.shape[1]
    r_grid = (np.arange(n_r) + 0.5) * s.grid_m
    streak = streak_columns(img, s.grid_m, s)
    img = _flatten_background(img, streak)
    f_h_sat = float(np.asarray(fld(np.array([h_sat_m]))[0])[0])
    pri = onboard_priors(freq, img, streak, f_h_sat, s.grid_m, s)
    score = _llr(_shadow_suppress(img, s.grid_m, s), s)
    half = _step_halfwidths(freq)
    base = np.ones_like(img, dtype=bool)
    base[streak] = False

    # The ground echo: from the top of the band down its range rises towards foF2.
    ground_path = np.full(n, -1, dtype=int)
    ground_mask = np.zeros_like(base)
    ground_branches = []
    g_allowed = base & (r_grid[None, :] >= s.ground_min_frac_hsat * h_sat_m)
    top = freq >= np.quantile(freq, 0.75) if n else np.zeros(0, dtype=bool)
    for _ in range(2):
        st = np.where(top[:, None], 0.0, np.inf) * np.ones((1, n_r))
        p_rev, tot = viterbi_path(score[::-1], g_allowed[::-1], st[::-1], s.grid_m, s)
        p = p_rev[::-1]
        cols = np.flatnonzero(p >= 0)
        if cols.size < s.ground_min_columns or tot < s.min_track_score:
            break
        # A ground echo runs to the top of the band, mostly over the detector
        # threshold; a path that does not is not one.
        det = img[cols, p[cols]] > s.design_snr_db
        if freq[cols[-1]] < np.quantile(freq, 0.85) or np.mean(det) < 0.5:
            break
        # The branch starts where its contiguous run from the top of the band
        # starts (gaps of one column allowed, streak columns not counted); a
        # gappy tail below that is not
        # the ground echo, and would put foF2 too low.
        k = cols.size - 1
        while k > 0 and np.count_nonzero(~streak[cols[k - 1] + 1: cols[k]]) <= s.max_gap:
            k -= 1
        cols = cols[k:]
        strong = img[cols, p[cols]] > s.design_snr_db
        k = 0
        while k < cols.size and not (strong[k] and np.mean(strong[k: k + 5]) >= 0.8):
            k += 1
        cols = cols[k:]
        if cols.size < s.ground_min_columns:
            break
        ground_branches.append({"f_start_hz": float(freq[cols[0]]), "n": int(cols.size), "score": tot,
                                "range_top_m": float(r_grid[p[cols[-1]]])})
        marg = int(round(s.ground_margin_m / s.grid_m))
        # fill the gaps and mask at and beyond the echo in every column of the run
        span_cols = np.arange(cols[0], n)
        r_fill = np.interp(span_cols, cols, p[cols].astype(float))
        for i, rf in zip(span_cols, r_fill, strict=True):
            lo = max(int(round(rf)) - marg, 0)
            g_allowed[i, lo: int(round(rf)) + marg + 1] = False
            ground_mask[i, lo:] = True
            if ground_path[i] < 0 or rf < ground_path[i]:
                ground_path[i] = int(round(rf))
    gcols = np.flatnonzero(ground_path >= 0)
    fo_ground = float(min(b["f_start_hz"] for b in ground_branches)) if ground_branches else float("nan")
    allowed = base & ~ground_mask

    f_n, f_x = pri["f_n_hz"], pri["f_x_hz"]
    if not np.isfinite(f_n):
        f_n = float(freq[0]) if n else float("nan")
        f_x = float(0.5 * f_h_sat + math.sqrt(0.25 * f_h_sat**2 + f_n**2)) if np.isfinite(f_n) else f_n
    near_start = r_grid <= s.start_max_m

    def start_cost(f0: float, far: bool = False, lo_hz: float | None = None) -> np.ndarray:
        # ``far``: a far-range start in the first columns above the prior costs
        # ``far_start_cost`` (the O trace under a shallow gradient; X starts
        # near range).  ``lo_hz``: the lowest column a path may start in.
        out = np.full((n, n_r), np.inf)
        ok_col = freq >= (f0 if lo_hz is None else lo_hz) - 1.5 * half
        early = (freq >= f0 - 1.5 * half) & (freq <= f0 * (1.0 + 0.05 * s.start_steps) + half)
        out[ok_col] = s.late_start_cost
        if far:
            out[early] = min(s.far_start_cost, s.late_start_cost)
        out[np.ix_(early, near_start)] = 0.0
        return out

    # Above its penetration frequency a mode has no trace: the ground echo
    # branches give foF2 (the lower) and fxF2 (the higher) when they are seen.
    starts = sorted(b["f_start_hz"] for b in ground_branches)
    cap_o = starts[0] + 1.5 * float(np.max(half)) if starts else np.inf
    # fxF2 - foF2 is about half the gyrofrequency at the peak; with one branch
    # seen, the X cap is set that far above foF2 (f_H at the spacecraft is a
    # lower bound for the peak's, so 0.75 f_H leaves room).
    if len(starts) > 1:
        cap_x = starts[1] + 1.5 * float(np.max(half))
    elif starts:
        cap_x = starts[0] + 0.75 * f_h_sat + 1.5 * float(np.max(half))
    else:
        cap_x = np.inf
    # The trace floor is an upper bound on f_N (a weak start below it may have
    # been missed), so with that prior the O search is not cut off below it.
    o_lo = float(freq[0]) if pri["f_n_method"] == "trace_floor" and n else f_n
    region_o = ((freq >= o_lo - 1.5 * half) & (freq <= cap_o))[:, None] & np.ones((1, n_r), dtype=bool)
    region_x = ((freq >= f_x - 1.5 * half) & (freq <= cap_x))[:, None] & np.ones((1, n_r), dtype=bool)
    sc_o = start_cost(f_n, far=pri["f_n_method"] == "trace_floor", lo_hz=o_lo)
    sc_x = start_cost(f_x)
    sep = int(round(s.order_sep_m / s.grid_m))

    def order_mask(other: np.ndarray, nearer: bool) -> np.ndarray:
        m = np.ones((n, n_r), dtype=bool)
        for i in np.flatnonzero(other >= 0):
            if nearer:
                m[i, : other[i] + sep + 1] = False
            else:
                m[i, max(other[i] - sep, 0):] = False
        return m

    def gates(pred: np.ndarray) -> np.ndarray:
        # The prediction is good to a few kilometres along the trace and lags
        # near the cusp, where the trace is steep: the gate is a floor, a
        # fraction of the range, and half the predicted change per column.
        g = np.maximum(s.x_gate_m, s.x_gate_frac * np.nan_to_num(pred))
        slope = np.full(n, 0.0)
        fin = np.isfinite(pred)
        if np.count_nonzero(fin) >= 2:
            dp = np.abs(np.gradient(np.where(fin, pred, np.interp(np.arange(n), np.flatnonzero(fin), pred[fin]))))
            slope = np.where(fin, dp, 0.0)
        return g + s.x_gate_slope * slope

    def gate_mask(pred: np.ndarray, inside: bool) -> np.ndarray:
        m = np.ones((n, n_r), dtype=bool)
        g = gates(pred)
        for i in np.flatnonzero(np.isfinite(pred)):
            near = np.abs(r_grid - pred[i]) <= g[i]
            m[i] = near if inside else ~near
        return m

    p_x, tot_x = viterbi_path(score, allowed & region_x, sc_x, s.grid_m, s)
    p_o = np.full(n, -1, dtype=int)
    tot_o = float("-inf")
    x_pred = np.full(n, np.nan)
    history = []
    for it in range(max(int(s.iterations), 1)):
        am = allowed & region_o & order_mask(p_x, nearer=True)
        if np.any(np.isfinite(x_pred)):
            # Where the X path has a point the ordering already keeps O beyond
            # it; the predicted X trace is kept out of O only where X was not
            # found.  Excluding it everywhere fed a prediction error back: a
            # prediction a few km long excluded the true O, O moved onto the
            # spread behind it, and the inversion moved the prediction further.
            no_x = np.where(p_x >= 0, np.nan, x_pred)
            am &= gate_mask(no_x, inside=False)
        caps = None
        if s.cusp_up_cap_m > s.up_cap_m:
            # The O trace rises steeply at both ends: into the cusp below foF2,
            # and, under a shallow topside gradient, just above f_N, where the
            # group retardation grows from zero within a step or two.
            prev = np.flatnonzero(p_o >= 0)
            f_top = cap_o if np.isfinite(cap_o) else (1.1 * float(freq[prev[-1]]) if prev.size else np.nan)
            steep = np.zeros(n, dtype=bool)
            if np.isfinite(f_top) and np.isfinite(f_n):
                steep |= freq >= f_top - s.cusp_band * (f_top - f_n)
            if np.isfinite(f_n):
                steep |= (freq >= o_lo - 1.5 * half) & (freq <= f_n * (1.0 + 0.1 * s.start_steps) + half)
            caps = np.where(steep, s.cusp_up_cap_m, s.up_cap_m)
        p_o, tot_o = viterbi_path(score, am, sc_o, s.grid_m, s, up_cap_cols=caps)
        cols_o = np.flatnonzero(p_o >= 0)
        if cols_o.size >= 2:
            x_pred = predict_x_range_m(freq[cols_o], r_grid[p_o[cols_o]], min(f_n, float(freq[cols_o[0]]) * 0.999),
                                       h_sat_m, fld, freq)
        ax = allowed & region_x & order_mask(p_o, nearer=False)
        if np.any(np.isfinite(x_pred)):
            ax &= gate_mask(x_pred, inside=True) | ~np.isfinite(x_pred)[:, None]
        p_x_new, tot_x = viterbi_path(score, ax, sc_x, s.grid_m, s)
        history.append({"iteration": it, "o_points": int(cols_o.size),
                        "x_points": int(np.count_nonzero(p_x_new >= 0)), "o_score": tot_o, "x_score": tot_x})
        if np.array_equal(p_x_new, p_x):
            p_x = p_x_new
            break
        p_x = p_x_new

    def accept(p: np.ndarray, tot: float) -> np.ndarray:
        if np.count_nonzero(p >= 0) < s.min_track_columns or not np.isfinite(tot) or tot < s.min_track_score:
            return np.full(n, -1, dtype=int)
        return p

    p_o, p_x = accept(p_o, tot_o), accept(p_x, tot_x)

    # samples
    rng_out = np.full((n, 2), np.nan)
    snr_out = np.full((n, 2), np.nan)
    sig_out = np.full((n, 2), np.nan)
    conf = np.full((n, 2), np.nan)
    qual = np.zeros((n, 2), dtype=int)
    spans = list(span_rate_hz_per_s2) if span_rate_hz_per_s2 is not None else [None] * n
    for m, p in ((0, p_o), (1, p_x)):
        for i in np.flatnonzero(p >= 0):
            smp = _refine(spectra[i], cfars[i], float(r_grid[p[i]] - 0.5 * s.grid_m), s, spans[i])
            rng_out[i, m] = smp["range_m"]
            snr_out[i, m] = smp["snr_db"]
            sig_out[i, m] = _sigma(smp["snr_db"], smp["cell_m"])
            qual[i, m] = smp["q"]
    # labelling confidence from the predicted X trace
    g_all = gates(x_pred)
    for i in range(n):
        pred = x_pred[i]
        g = g_all[i]
        if np.isfinite(rng_out[i, 0]):
            if freq[i] < f_x - half[i]:
                conf[i, 0] = 1.0
            elif np.isfinite(pred):
                conf[i, 0] = float(1.0 - math.exp(-0.5 * ((rng_out[i, 0] - pred) / g) ** 2))
            else:
                conf[i, 0] = 0.75
        if np.isfinite(rng_out[i, 1]):
            if np.isfinite(pred):
                conf[i, 1] = float(math.exp(-0.5 * ((rng_out[i, 1] - pred) / g) ** 2))
            else:
                conf[i, 1] = 0.5
        for m in (0, 1):
            if np.isfinite(conf[i, m]) and conf[i, m] < 0.5:
                qual[i, m] |= Q_LOW_CONFIDENCE
    # consistency statistic
    both = np.isfinite(rng_out[:, 1]) & np.isfinite(x_pred)
    resid = rng_out[both, 1] - x_pred[both]
    consistency = {"n": int(both.sum()),
                   "median_abs_m": float(np.median(np.abs(resid))) if resid.size else float("nan"),
                   "rms_m": float(np.sqrt(np.mean(resid**2))) if resid.size else float("nan")}

    # G4: blended samples
    for i in range(n):
        if not (np.isfinite(rng_out[i, 0]) and np.isfinite(rng_out[i, 1])):
            continue
        cell = float(spectra[i].range_resolution_m)
        sep_m = rng_out[i, 0] - rng_out[i, 1]
        if sep_m >= s.blend_cells * cell:
            continue
        qual[i, 0] |= Q_BLENDED
        qual[i, 1] |= Q_BLENDED
        if s.blend_policy == "reject":
            qual[i, 0] |= Q_REJECTED
        elif s.blend_policy == "inflate":
            sig_out[i, 0] = float(math.hypot(sig_out[i, 0], max(cell - max(sep_m, 0.0), 0.5 * cell)))
        else:
            fit = _separate(spectra[i], rng_out[i, 1], max(rng_out[i, 0], rng_out[i, 1] + 0.5 * cell))
            if fit is None:
                qual[i, 0] |= Q_REJECTED
            else:
                rng_out[i, 1], rng_out[i, 0] = fit

    # cusp model on the O trace
    cusp = None
    ok_o = np.isfinite(rng_out[:, 0]) & ((qual[:, 0] & Q_REJECTED) == 0)
    if np.count_nonzero(ok_o) >= s.cusp_min_samples:
        cusp = cusp_fit(freq[ok_o], rng_out[ok_o, 0], fo_hint_hz=fo_ground)
        if cusp is not None:
            idx = np.flatnonzero(ok_o)
            idx = idx[np.argsort(freq[idx])]
            # O cannot reflect at or above foF2: points the fit left out, or
            # beyond its critical frequency, are removed.
            for j in idx:
                if freq[j] >= cusp["fo_hz"] or j in idx[idx.size - cusp["n_dropped"]:]:
                    qual[j, 0] |= Q_CUSP | Q_REJECTED
            idx = idx[: idx.size - cusp["n_dropped"]]
            top_idx = idx[-cusp["n"]:]
            model = _cusp_eval(cusp, freq[top_idx])
            tol = s.cusp_sigma_factor * np.maximum(sig_out[top_idx, 0], cusp["rms_m"])
            bad = np.abs(rng_out[top_idx, 0] - model) > tol
            over = (qual[top_idx, 0] & Q_OVER_SPAN) != 0
            for j, b, o, mv in zip(top_idx, bad, over, model, strict=True):
                if b or o:
                    qual[j, 0] |= Q_CUSP
                    if s.cusp_policy == "replace":
                        rng_out[j, 0] = float(mv)
                        sig_out[j, 0] = float(max(sig_out[j, 0], cusp["rms_m"]))

    f_stop = float("nan")
    if s.stop_at_loss:
        f_stop = stop_at_loss(freq, rng_out[:, 0], qual[:, 0], conf[:, 0], s)
        if np.isfinite(f_stop):
            beyond = np.isfinite(rng_out[:, 0]) & (freq >= f_stop)
            qual[beyond, 0] |= Q_STOPPED | Q_REJECTED

    out_rng = rng_out.copy()
    rej = (qual & Q_REJECTED) != 0
    out_rng[rej] = np.nan
    info = {
        "priors": pri,
        "ground": {"branches": ground_branches, "fo_from_ground_hz": fo_ground,
                   "n_columns": int(gcols.size)},
        "x_consistency": consistency,
        "x_pred_km": [round(float(v) / 1e3, 1) if np.isfinite(v) else None for v in x_pred],
        "iterations": history,
        "cusp": cusp,
        "stop_hz": f_stop,
        "n_o": int(np.count_nonzero(np.isfinite(out_rng[:, 0]))),
        "n_x": int(np.count_nonzero(np.isfinite(out_rng[:, 1]))),
        "settings": {k: getattr(s, k) for k in s.__dataclass_fields__},
    }
    return TrackResult(freq_hz=freq, range_m=out_rng, snr_db=snr_out, sigma_m=sig_out, confidence=conf,
                       quality=qual, info=info)


def settings_from(**kw: Any) -> TrackerSettings:
    """TrackerSettings with unknown keys refused."""
    return replace(TrackerSettings(), **kw)
