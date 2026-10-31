"""Layer: L3
WP: WP4
Claims: C1, C4
Status: implemented
Requirements: INS-10

The detector (WP4 Phase 2): a cell test with a stated probability of false alarm.

Decisions of 2026-10-30 (gate G1): a false-alarm probability of 1e-4 per range
cell; the noise reference taken from the echo-free negative-beat half of the
spectrum, with a local check so that a strong return's window sidelobes are not
called detections; single-cell detection here, with integration along the trace
left to the tracker of Phase 4; and INS-10 measured against this noise
reference rather than against noise plus clutter.

The local check as built is sidelobe blanking against the window's known
envelope (below).  The local ordered-statistic veto first built for it is kept
as an option and is off by default: an echo smeared over many range cells by
dispersion within the dwell (50 to 100 km near f_N and the cusp) fills its own
reference window and was vetoed, which lost whole traces on the equatorial
evening and mid-latitude night profiles of ``wp4_campaign``.  That is the
self-masking of an extended target that an ordered-statistic CFAR is known for.

THE NOISE REFERENCE
-------------------
The receiver samples complex baseband, and after the de-chirp's orientation
every echo beats at a positive frequency (``+K tau`` on an up-sweep).  The
negative half of each block's transform therefore holds the receiver noise, the
sky noise and any carrier, and no echo and no clutter return.  It is about as
long as the positive half, a thousand cells or more, so its noise level is
nearly a known quantity: the detector is close to the Neyman-Pearson test for
a known noise level, without the loss a short CFAR window carries and without
the bias that a strong echo, a resonance continuum or spread F puts into an
estimate taken where they are.

Three things can still reach the negative half, and the estimator is chosen to
be indifferent to them:

* the image of the leakage near zero beat, which the guard ``near_guard_cells``
  removes;
* returns beyond the unambiguous range, which alias there if the anti-alias
  filter passes them (a ground echo just above penetration);
* in hardware, the image of strong positive-beat returns through I/Q
  imbalance, typically 30 to 40 dB down.

The estimator is the sample median, an ordered statistic that a few
contaminated cells do not move.  The threshold is ``alpha * X_(k)``, with ``k``
the middle of ``N`` and ``alpha`` from the ordered-statistic false-alarm
probability for independent exponential cells (Rohling 1983):

    P_fa = prod_{i=0}^{k-1} (N - i) / (N - i + alpha)

with ``N`` the number of independent cells, which is the number of bins used
over the zero-padding factor and the window's noise bandwidth.  As ``N`` grows
this tends to the known-noise threshold ``-ln(P_fa)`` times the mean.

THE LOCAL CHECKS
----------------
A Hann window's first sidelobe is 31.5 dB down, so a return 40 dB over the
noise, such as the ground echo or a resonance, has sidelobes that pass the
global test.

*Ordered-statistic veto (optional, ``local_pfa``; off by default).*  Each cell
must also exceed ``beta`` times the ``local_order`` ordered statistic of
``local_ref_cells`` cells on each side, beyond ``guard_cells``.  ``beta`` is set for ``local_pfa`` (1e-2) under the same
formula.  In homogeneous noise the global test is far the stricter of the two,
so the veto costs nothing.  Beside a strong return its window is filled by the
sidelobes or the continuum, and a cell must stand out of that to count.  An
ordered statistic at three quarters keeps a weaker echo beside a stronger one
(O beside X) from being vetoed by it.

It masks range-extended echoes (above), and it cannot remove the first one or
two sidelobes of a very strong return either: just outside the guard cells, the cell's own neighbourhood holds
the main lobe on one side and lower sidelobes on the other, and the cell stands
out of it.  On the near side of a strong echo those sidelobes would move the
leading edge early.

*Sidelobe blanking (the default local check, ``sidelobe_margin_db``).*  The
window's response is known exactly, so the detector blanks every cell that lies under a stronger detected peak's sidelobe
envelope plus ``sidelobe_margin_db`` (6 dB): the envelope is the upper bound
of the window's transform over all fractional bin offsets, so it holds for an
echo anywhere in its bin.  A genuine return under the skirt that is more than
the margin above it survives.

WHAT IS STATED
--------------
``P_fa`` is per range cell: the marginal probability that a noise-only cell is
declared.  Zero padding makes adjacent bins correlated without changing that
marginal probability; it changes how many distinct false alarms an ionogram
holds.  The number of false alarms an ionogram holds is reported, not
assumed.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Any

import numpy as np

__all__ = [
    "DEFAULT_PFA",
    "DEFAULT_LOCAL_PFA",
    "CfarResult",
    "os_alpha",
    "spectrum_pad_factor",
    "noise_reference",
    "sidelobe_envelope",
    "cfar_detect",
]

DEFAULT_PFA = 1.0e-4
DEFAULT_LOCAL_PFA = 1.0e-2
LN2 = math.log(2.0)


def os_alpha(pfa: float, n: int, k: int) -> float:
    """Threshold multiplier on the ``k``-th smallest of ``n`` exponential cells for ``pfa``.

    Solves ``sum_{i<k} ln((n - i) / (n - i + alpha)) = ln(pfa)`` by bisection.
    """
    if not 0.0 < pfa < 1.0:
        raise ValueError("pfa must be in (0, 1)")
    n, k = int(n), int(k)
    if n < 1 or not 1 <= k <= n:
        raise ValueError(f"need 1 <= k <= n; got n={n}, k={k}")
    i = np.arange(k, dtype=float)
    target = math.log(pfa)

    def f(a: float) -> float:
        return float(np.sum(np.log((n - i) / (n - i + a)))) - target

    lo, hi = 0.0, 1.0
    while f(hi) > 0.0:
        hi *= 2.0
        if hi > 1e12:
            raise ValueError("no threshold reaches that pfa with this reference")
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if f(mid) > 0.0:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def spectrum_pad_factor(spectrum: Any) -> float:
    """Zero-padding factor of a :class:`~tss.l3_signal.transform.RangeSpectrum`."""
    cell_bin = float(spectrum.range_resolution_m) / float(spectrum.enbw_bins)
    return max(cell_bin / float(spectrum.bin_spacing_m), 1.0)


def _independent(n_bins: int, spectrum: Any) -> int:
    return max(int(round(n_bins / (spectrum_pad_factor(spectrum) * float(spectrum.enbw_bins)))), 1)


def noise_reference(spectrum: Any, *, channel: int = 0, near_guard_cells: float = 4.0
                    ) -> dict[str, Any]:
    """The noise reference of one spectrum: median, mean and cell counts.

    From the negative-beat half when the spectrum carries it; otherwise from the
    positive half (the WP3-era estimate), which is flagged.
    """
    neg = getattr(spectrum, "negative_power_db", None)
    pad = spectrum_pad_factor(spectrum)
    if neg is not None and np.size(neg):
        row = np.asarray(neg, dtype=float)[int(channel)]
        g0 = int(math.ceil(near_guard_cells * pad * float(spectrum.enbw_bins)))
        use = row[g0:]
        source = "negative_half"
    else:
        use = np.asarray(spectrum.power_db, dtype=float)[int(channel)]
        source = "positive_half_median"
    use = use[np.isfinite(use)]
    if use.size == 0:
        return {"source": source, "median_w": float("nan"), "mean_db": float("-inf"), "n_bins": 0, "n_eff": 0}
    med = float(np.median(10.0 ** (use / 10.0)))
    m = int(getattr(spectrum, "n_incoherent", 1) or 1)
    return {
        "source": source,
        "median_w": med,
        "mean_db": 10.0 * math.log10(med / _median_over_mean(m)) if med > 0 else float("-inf"),
        "n_bins": int(use.size),
        "n_eff": _independent(int(use.size), spectrum),
    }


def _median_over_mean(m: int) -> float:
    """Median over mean of the average of ``m`` exponential cells (ln 2 for one)."""
    if m <= 1:
        return LN2
    from scipy.special import gammaincinv

    return float(gammaincinv(m, 0.5)) / m


def _quantile_over_mean(m: int, p: float) -> float:
    """Upper ``p`` point over the mean of the average of ``m`` exponential cells."""
    if m <= 1:
        return -math.log(p)
    from scipy.special import gammaincinv

    return float(gammaincinv(m, 1.0 - p)) / m


_GAMMA_ALPHA: dict[tuple[int, int, float], float] = {}


def gamma_median_alpha(pfa: float, m: int, n_eff: int, draws: int = 4000) -> float:
    """Threshold over the reference median for averages of ``m`` exponential cells.

    The reference median is itself estimated from ``n_eff`` cells, and its
    scatter raises the false-alarm probability of a known-noise threshold
    (by a factor of 1.5 at ``m = 4`` and ``n_eff`` of about 170, measured
    2026-10-30).  So the threshold solves ``E[Q(m, m alpha med)] = pfa`` over the
    sampling distribution of the median, drawn once with a fixed seed and
    cached.  ``Q`` is the regularised upper incomplete gamma function.
    """
    from scipy.special import gammaincc

    key = (int(m), int(n_eff), float(pfa))
    if key not in _GAMMA_ALPHA:
        rng = np.random.default_rng(20261030)
        med = np.median(rng.gamma(m, 1.0 / m, size=(draws, int(n_eff))), axis=1)

        def f(a: float) -> float:
            return float(np.mean(gammaincc(m, m * a * med))) - pfa

        lo, hi = 0.0, 1.0
        while f(hi) > 0.0:
            hi *= 2.0
        for _ in range(100):
            mid = 0.5 * (lo + hi)
            if f(mid) > 0.0:
                lo = mid
            else:
                hi = mid
        _GAMMA_ALPHA[key] = 0.5 * (lo + hi)
    return _GAMMA_ALPHA[key]


@dataclass(frozen=True)
class CfarResult:
    """The detector's decision for one spectrum."""

    mask: np.ndarray
    noise_floor_db: float
    threshold_db: float
    local_threshold_db: np.ndarray | None
    pfa: float
    alpha: float
    n_reference: int
    n_eff: int
    source: str
    flags: dict[str, Any] = field(default_factory=dict)


def _local_os(p: np.ndarray, ref: int, guard: int, order: float) -> np.ndarray:
    """Per-bin ``order`` quantile (ordered statistic) of the cells beside it, NaN-aware."""
    n = p.size
    pad = guard + ref
    pp = np.concatenate([np.full(pad, np.nan), p, np.full(pad, np.nan)])
    from numpy.lib.stride_tricks import sliding_window_view

    win = sliding_window_view(pp, ref)
    left = win[np.arange(n)]  # bins i - guard - ref .. i - guard - 1
    right = win[np.arange(n) + 2 * guard + ref + 1]  # bins i + guard + 1 .. i + guard + ref
    both = np.sort(np.concatenate([left, right], axis=1), axis=1)  # NaN sort last
    m = np.sum(np.isfinite(both), axis=1)
    idx = np.clip(np.floor(order * m).astype(int), 0, both.shape[1] - 1)
    out = both[np.arange(n), idx]
    out[m == 0] = np.nan
    return out


_ENVELOPES: dict[tuple[str, int], tuple[np.ndarray, float]] = {}


def sidelobe_envelope(window: str, n: int, oversample: int = 32) -> tuple[np.ndarray, float]:
    """Upper envelope of a window's normalised transform, and its first null.

    Returns ``(env, null)``: ``env[j]`` is the largest ``|W(u)| / W(0)`` over
    every ``u >= j / oversample`` (in unpadded bins), and ``null`` the first
    null in bins.  Cached per (window, length).
    """
    from tss.l3_signal.transform import window_taper

    key = (str(window).lower(), int(n))
    if key not in _ENVELOPES:
        taper = window_taper(key[0], key[1])
        w = np.abs(np.fft.fft(taper, oversample * key[1]))[: oversample * key[1] // 2]
        w = w / w[0]
        dif = np.diff(w)
        mins = np.flatnonzero((dif[:-1] < 0) & (dif[1:] >= 0)) + 1
        null = float(mins[0]) / oversample if mins.size else 1.0
        env = np.maximum.accumulate(w[::-1])[::-1]
        _ENVELOPES[key] = (env, null)
    return _ENVELOPES[key]


def _sidelobe_blank(p: np.ndarray, candidates: np.ndarray, spectrum: Any, margin_db: float) -> np.ndarray:
    """Cells under a stronger candidate peak's sidelobe envelope plus the margin."""
    pad = spectrum_pad_factor(spectrum)
    n_fft = 2 * (p.size - 1)
    n = max(int(round(n_fft / pad)), 2)
    env, null = sidelobe_envelope(str(spectrum.window), n)
    over = env.size / (n / 2.0)  # envelope samples per unpadded bin
    blank = np.zeros(p.size, dtype=bool)
    idx = np.flatnonzero(candidates)
    if idx.size == 0:
        return blank
    # Peaks: candidate cells that are local maxima.
    left = np.concatenate([[-np.inf], p[:-1]])
    right = np.concatenate([p[1:], [-np.inf]])
    peaks = idx[(p[idx] >= left[idx]) & (p[idx] >= right[idx])]
    margin = 10.0 ** (margin_db / 10.0)
    cells = np.arange(p.size)
    for c in peaks[np.argsort(p[peaks])[::-1]]:
        u = np.abs(cells - c) / pad
        side = u >= null
        j = np.minimum((u * over).astype(int), env.size - 1)
        under = side & (p <= p[c] * env[j] ** 2 * margin)
        blank |= under
    return blank


def cfar_detect(spectrum: Any, *, channel: int = 0, pfa: float = DEFAULT_PFA,
                local_pfa: float | None = None, local_ref_cells: int = 16,
                guard_cells: int = 2, local_order: float = 0.75, near_guard_cells: float = 4.0,
                sidelobe_margin_db: float | None = 6.0,
                max_range_m: float | None = None, exclude_below_m: float = 0.0) -> CfarResult:
    """Cell decisions for one range spectrum (module docstring)."""
    power_db = np.asarray(spectrum.power_db, dtype=float)[int(channel)]
    rng = np.asarray(spectrum.range_m, dtype=float)
    p = np.where(np.isfinite(power_db), 10.0 ** (power_db / 10.0), np.nan)
    ref = noise_reference(spectrum, channel=channel, near_guard_cells=near_guard_cells)
    n_eff = int(ref["n_eff"])
    if n_eff < 2 or not np.isfinite(ref["median_w"]):
        return CfarResult(np.zeros(p.size, dtype=bool), float("-inf"), float("inf"), None, pfa,
                          float("nan"), int(ref["n_bins"]), n_eff, str(ref["source"]),
                          {"no_reference": True})
    looks = max(int(getattr(spectrum, "n_looks", 1) or 1), 1)
    m_inc = max(int(getattr(spectrum, "n_incoherent", 1) or 1), 1)
    # The stated P_fa holds over all looks (WP4 Phase 3 rate search): each
    # look is tested at 1 - (1 - P_fa)^(1/L), a bound since the looks are
    # correlated.
    pfa_look = -math.expm1(math.log1p(-pfa) / looks)
    if m_inc == 1:
        k = max(int(round(0.5 * n_eff)), 1)
        alpha = os_alpha(pfa_look, n_eff, k)
    else:
        # Average of m_inc sub-blocks (option 1): a gamma variate, against a
        # median estimated from fewer independent cells (gamma_median_alpha).
        alpha = gamma_median_alpha(pfa_look, m_inc, n_eff)
    thr = alpha * float(ref["median_w"])
    mask = np.isfinite(p) & (p > thr)
    local_db = None
    if local_pfa is not None and local_ref_cells > 0:
        padf = spectrum_pad_factor(spectrum)
        r_bins = max(int(round(local_ref_cells * padf)), 1)
        g_bins = int(round(guard_cells * padf))
        n_loc = _independent(2 * r_bins, spectrum)
        beta = os_alpha(float(local_pfa), n_loc, max(min(int(round(local_order * n_loc)), n_loc), 1))
        loc = beta * _local_os(np.where(np.isfinite(p), p, np.nan), r_bins, g_bins, local_order)
        mask &= ~(np.isfinite(loc) & (p <= loc))
        with np.errstate(divide="ignore"):
            local_db = 10.0 * np.log10(loc)
    if sidelobe_margin_db is not None:
        mask &= ~_sidelobe_blank(np.where(np.isfinite(p), p, 0.0), mask, spectrum, float(sidelobe_margin_db))
    if max_range_m is not None:
        mask &= rng <= float(max_range_m)
    if exclude_below_m > 0.0:
        mask &= rng >= float(exclude_below_m)
    return CfarResult(
        mask=mask,
        noise_floor_db=float(ref["mean_db"]),
        threshold_db=10.0 * math.log10(thr),
        local_threshold_db=local_db,
        pfa=float(pfa),
        alpha=float(alpha),
        n_reference=int(ref["n_bins"]),
        n_eff=n_eff,
        source=str(ref["source"]),
        flags={"threshold_over_mean_db": 10.0 * math.log10(alpha * _median_over_mean(m_inc)),
               "n_looks": looks, "n_incoherent": m_inc, "pfa_per_look": pfa_look},
    )
