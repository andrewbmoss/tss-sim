"""Layer: L3
WP: WP5
Claims: C3
Status: partial
Requirements: SCI-15, INS-15

The drift tier (``l3/drift``): route one at signal level on a drift-mode
sequence (WP5 Phase 3; ``docs/drift-mode.md`` sections 2 to 5).

WHAT IT DOES
------------
:meth:`DriftProcessing.process_sequence` takes the :class:`ReceivedSequence` of
``Pipeline.run_sequence`` and returns a :class:`DriftProduct` with the
contract's names.  One receive channel (channel 0).  Steps:

1.  **Detect** (decision (j)).  Each record is transformed over the samples
    after the longest delay searched (``max_range_m``, or a range gate from
    the ionogram's trace, :func:`range_gate_from_ionogram`); the power is
    summed over ``window_sweeps`` consecutive sweeps (21, sliding), and a
    (window, range bin) cell is declared against a threshold set from the
    noise in the negative-beat half.  The false-alarm rate is stated per
    sequence (``pfa_per_sequence``, 1e-2; decisions (m) and (q), 2026-11-27)
    and divided by the number of cells in the sequence; ``pfa``, per cell as
    WP4 states it, overrides it when given.
    The strongest ``max_echoes`` range peaks become echoes, sorted by range;
    the beat of each is refined on the window that peaked.  The summed
    detector of part 1 (the whole sequence, a 13 dB rule) is kept as
    ``detector="summed"``.
2.  **Gate.**  For each echo the first group delay of every record holds the
    previous sweep's echo at another beat; L2 leaves it empty and this tier
    drops it (decision (f), 2026-11-25): the gate starts ``guard_samples``
    after the echo's delay.  The fraction dropped is reported
    (``blanked_fraction``) and is a loss of SNR per sweep.
3.  **One complex value per sweep**, the gated record demodulated at the echo's
    beat with time measured from the gate's centre ``t_ref``, stamped with the
    time the echo was reflected, ``t_n + t_ref - tau / 2``.  Its phase is the
    echo's carrier phase then, at the frequency transmitted for it,
    ``f_eff = f_start + K (t_ref - tau)`` (D8 condition 2): velocities are
    converted at ``f_eff``, not at the sub-band's centre.
4.  **Platform Doppler.**  ``f_D,plat = (2 f_eff / c) mu_sat v_down`` for a
    local-down wave normal is removed by a phase ramp (section 2, step 1).
5.  **Route one** (:func:`estimate_route_one`, reused): Doppler from phase
    steps, the power envelope's peak ``t_pk``, the zero crossing ``t_zc`` with
    the ephemeris slope, and ``drift = v_eff^2 (t_zc - t_pk) / (mu_sat P')``,
    with its refusals (extended, no peak).
6.  **Joint estimators** (``estimator``; part 3, decisions (o), (r) and (s)).
    ``"joint"`` fits amplitude, ``t0``, drift and width to the complex values
    (deterministic likelihood); ``"joint_fading"`` fits power, ``t0``, drift,
    width and ``tau_c`` to a Rayleigh-fading likelihood with the noise fixed at
    the detector's (:mod:`tss.l3_signal.novel.route_one_joint`).  Step 5 stays:
    it is the gate, so all three form a drift on the same echoes, and its
    answer is the start and is kept as ``drift_pinned_mps``.  A joint fit that
    fails leaves the pinned drift, flagged ``JOINT_FAILED_PINNED_USED``.
    **Ambiguity check** (``ambiguity_margin``, ``"joint_fading"`` only, off by
    default; decision 2B, 2026-11-27): when the likelihood has a second
    maximum, well apart in drift, within the margin of the best in log
    likelihood, the echo is declined: ``drift_mps``, its sigma, ``t_pk_s``,
    ``t_zc_s`` and ``along_track_m`` are NaN, flagged ``DECLINED: ambiguous``.
    ``ambiguity_nll_margin`` and ``drift_alt_mps`` are reported either way.
    The default stays ``"pinned"`` for now; a later switch is to
    ``"joint_fading"``, not ``"joint"``, whose sigma is far too small under
    fading (decision (t) as amended by Andrew, 2026-11-27).

THE EFFECTIVE SPEED ON AN ORBIT
-------------------------------
Near overflight the phase path's second derivative is ``v^2 / P'_flat`` with
``P'_flat = int (r_p / r)^2 dz / mu`` (the group-path integral of the
Earth-flattened medium about the platform's radius), not ``v^2 / P'``.  So
route one uses ``v_eff^2 = v^2 P' / P'_flat``, which is ``v^2 r_f / r_p`` in
free space (4.3 % under ``v^2`` for a facet 300 km below a 600 km orbit).
``P'`` is measured; the facet's true depth, which ``P'_flat`` needs, comes from
the N(h) prior by inverting ``P' = int_0^D dz / mu`` (decision (d)).  The prior
is an explicit input (:class:`OnboardPrior`), and its source is in
``flags["prior_source"]`` (``"truth"`` when the truth profile was used).  The
onboard prior (Option A, 2026-11-25) is built by
:func:`onboard_prior_from_ionogram` from an ionogram of the same scene
(``Pipeline.run_prior_ionogram``): the WP4 tracker's ``f_N`` at the
spacecraft and an assumed scale height, reported in ``flags["prior"]``.  No
valid ``f_N`` means no prior, and no drift (``NO_PRIOR``).
On the straight track (platform states without velocity) the geometry is flat
and ``v_eff`` is ``v_s_mps``.

WHAT IS NOT DONE YET
--------------------
``tau_c_s``, ``doppler_spread_hz``, the field-aligned and
several-frequency products, the per-sweep products and the rest listed in
``flags["not_done"]``.  ``drift_sigma_mps`` is formed by the joint estimators
only (NaN with ``"pinned"``).  ``"joint_fading"`` also reports
``tau_c_fit_s``, the maximum-likelihood time of an exponential correlation;
``tau_c_s`` is left for the lag-one definition of ``docs/drift-mode.md``
section 4 (Andrew, 2026-11-27).  Every facet product carries
``TILT_UNFLAGGED``.

**Noise** (``docs/drift-tier.md`` sections 6 to 8).  With the beat given,
20 dB per sweep matches the parameter-level budget; at 10 dB route one's
slope test declines about half the compact facets, on an ideal field too.
The field is scaled to a unit peak before the estimator (its fit has absolute
tolerances).  The windowed detector finds the facet at 10 dB in every trial.
At 1e-4 per cell it also declared about one false echo per sequence, so the
rate is now stated per sequence (section 10).
"""

from __future__ import annotations

import math
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from typing import Any

import numpy as np

from tss.constants import C_LIGHT
from tss.interfaces import DriftProduct, ReceivedSequence
from tss.l3_signal.novel.route_one_estimator import estimate_route_one
from tss.l3_signal.novel.route_one_joint import (
    JointEstimate,
    estimate_route_one_joint,
    estimate_route_one_joint_fading,
    sample_noise_variance,
)
from tss.platform.orbits import RE
from tss.registry import register

__all__ = ["OnboardPrior", "DriftProcessing", "effective_speed_ratio", "onboard_prior_from_ionogram",
           "range_gate_from_ionogram"]

#: ``f_N^2 / N_e``, Hz^2 m^3.
_FN2_PER_NE = 80.6163

ESTIMATORS = ("pinned", "joint", "joint_fading")


@dataclass(frozen=True)
class OnboardPrior:
    """The N(h) profile the drift tier works from, and where it came from.

    ``alt_m``, ``ne``: the profile (m, m^-3), interpolated in ``log N_e``.
    ``source``: ``"truth"``, ``"onboard"`` (the ionogram chain's retrieval)
    or ``"configured"``; reported in the product's flags.
    """

    alt_m: np.ndarray
    ne: np.ndarray
    source: str = "configured"
    extra: Mapping[str, Any] = field(default_factory=dict)

    @classmethod
    def exponential(cls, *, f_n_sat_hz: float = 1.5e6, scale_height_m: float = 200e3,
                    ref_alt_m: float = 600e3, source: str = "configured") -> OnboardPrior:
        """The reference topside of TV Ch. 8 about ``ref_alt_m``, on a fine grid."""
        alt = np.linspace(ref_alt_m - 900e3, ref_alt_m + 200e3, 4401)
        fn2 = float(f_n_sat_hz) ** 2 * np.exp((ref_alt_m - alt) / float(scale_height_m))
        return cls(alt_m=alt, ne=fn2 / _FN2_PER_NE, source=source,
                   extra={"f_n_sat_hz": float(f_n_sat_hz), "scale_height_m": float(scale_height_m),
                          "ref_alt_m": float(ref_alt_m)})

    def mu_squared(self, f_hz: float, alt_m: np.ndarray) -> np.ndarray:
        alt = np.asarray(alt_m, dtype=float)
        order = np.argsort(self.alt_m)
        ln = np.interp(alt, np.asarray(self.alt_m)[order], np.log(np.clip(np.asarray(self.ne)[order], 1.0, None)))
        return 1.0 - _FN2_PER_NE * np.exp(ln) / float(f_hz) ** 2

    def mu_at(self, f_hz: float, alt_m: float) -> float:
        m2 = float(self.mu_squared(f_hz, np.array([alt_m]))[0])
        return math.sqrt(m2) if m2 > 0.0 else math.nan


def _group_integrals(prior: OnboardPrior, f_hz: float, alt_p: float, depth: float, r_p: float,
                     n: int = 400) -> tuple[float, float]:
    """``int_0^D dz / mu`` and ``int_0^D (r_p / r)^2 dz / mu`` below ``alt_p``."""
    x, w = np.polynomial.legendre.leggauss(n)
    z = 0.5 * depth * (1.0 + x)
    hw = 0.5 * depth * w
    m2 = prior.mu_squared(f_hz, alt_p - z)
    if np.any(m2 <= 0.0):
        return math.nan, math.nan
    inv = 1.0 / np.sqrt(m2)
    ratio = r_p / (r_p - z)
    return float(np.sum(hw * inv)), float(np.sum(hw * ratio * ratio * inv))


def effective_speed_ratio(prior: OnboardPrior, f_hz: float, alt_p: float, group_path_m: float,
                          r_p: float) -> tuple[float, float]:
    """``(v_eff / v, depth)`` for an echo of group path ``P'`` (module docstring).

    The depth is found from ``P' = int_0^D dz / mu`` by bisection.  When
    ``D = P'`` lies beyond the prior's reflection level the bracket is first
    cut back to that level, found by bisection between a depth where the wave
    propagates and one where it does not.  (Until 2026-11-25 it was halved
    until it propagated, which could cut it below the answer: on a Chapman
    topside with a facet 170 km down it returned 108 km.)  If ``P'`` exceeds
    the group path to the reflection level the depth returned is that level.
    """
    target = float(group_path_m)
    lo, hi = 0.0, target  # group path >= geometric depth
    if not math.isfinite(_group_integrals(prior, f_hz, alt_p, hi, r_p)[0]):
        ok, bad = 0.0, hi
        while bad - ok > 1.0:
            mid = 0.5 * (ok + bad)
            if math.isfinite(_group_integrals(prior, f_hz, alt_p, mid, r_p)[0]):
                ok = mid
            else:
                bad = mid
        hi = ok
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        g, _ = _group_integrals(prior, f_hz, alt_p, mid, r_p)
        if g < target:
            lo = mid
        else:
            hi = mid
    depth = 0.5 * (lo + hi)
    g, g_flat = _group_integrals(prior, f_hz, alt_p, depth, r_p)
    return math.sqrt(g / g_flat), depth


def onboard_prior_from_ionogram(
    ionogram: Any,
    platform: Any,
    *,
    scale_height_m: float = 200e3,
    mode_index: int = 0,
) -> tuple[OnboardPrior | None, dict[str, Any]]:
    """The drift tier's prior from an ionogram's tracker (WP5 Phase 3, Option A).

    The plasma frequency at the spacecraft is the WP4 tracker's estimate,
    ``ionogram.flags["tracker"]["priors"]["f_n_hz"]``, accepted only when it is
    finite, positive and below the lowest frequency detected on the trace
    (column ``mode_index``), as a plasma frequency at the spacecraft must be
    (the rule of :func:`tss.l4_retrieval.base.starting_plasma_frequency_hz`).
    The profile below the spacecraft is an exponential in ``N_e`` of the stated
    ``scale_height_m`` about the platform's altitude: **assumed, not
    measured**.  Only ``f_N,sat`` comes from the ionogram; the retrieved N(h)
    is later work (decision (i), 2026-11-25).  ``f_N,sat`` sets ``mu_sat``,
    which scales the drift; the scale height sets the facet's inferred depth
    and ``v_eff``, which move the slope ratio and the drift only at second
    order (``docs/drift-tier.md`` section 5).

    Returns ``(prior, report)``; ``prior`` is None when no valid ``f_N,sat``
    is found, and the drift tier then forms no drift and flags ``NO_PRIOR``
    (decision (k)).  ``report`` says why, and what the tracker gave.
    """
    flags = dict(getattr(ionogram, "flags", None) or {})
    pri = dict((flags.get("tracker") or {}).get("priors") or {})
    f_n = float(pri.get("f_n_hz", math.nan)) if pri.get("f_n_hz") is not None else math.nan
    report: dict[str, Any] = {
        "f_n_hz": f_n, "f_n_method": pri.get("f_n_method"),
        "trace_floor_hz": float(pri.get("trace_floor_hz", math.nan) or math.nan),
        "scale_height_m": float(scale_height_m), "scale_height_source": "assumed",
    }
    freq = np.asarray(getattr(ionogram, "freq_hz", []), dtype=float)
    vr = np.asarray(getattr(ionogram, "virtual_range_m", np.empty((0, 1))), dtype=float)
    if vr.ndim == 2 and vr.shape[0] == freq.size and 0 <= int(mode_index) < vr.shape[1]:
        det = freq[np.isfinite(vr[:, int(mode_index)]) & np.isfinite(freq)]
    else:
        det = np.empty(0)
    lowest = float(det.min()) if det.size else math.nan
    report["lowest_trace_hz"] = lowest
    if not pri:
        report["reason"] = "the ionogram carries no tracker priors (extractor not 'tracker')"
        return None, report
    if not (math.isfinite(f_n) and f_n > 0.0):
        report["reason"] = "the tracker found no plasma frequency at the spacecraft"
        return None, report
    if not math.isfinite(lowest):
        report["reason"] = "no trace was detected to check f_N against"
        return None, report
    if not f_n < lowest:
        report["reason"] = (f"f_N {f_n / 1e6:.4f} MHz is not below the lowest detected trace "
                            f"frequency {lowest / 1e6:.4f} MHz")
        return None, report
    report["reason"] = "accepted"
    alt = float(platform.alt_m)
    prior = OnboardPrior.exponential(f_n_sat_hz=f_n, scale_height_m=float(scale_height_m), ref_alt_m=alt,
                                     source="onboard")
    prior = OnboardPrior(alt_m=prior.alt_m, ne=prior.ne, source="onboard",
                         extra={**dict(prior.extra), **report})
    return prior, report


def range_gate_from_ionogram(
    ionogram: Any,
    f_hz: float,
    *,
    margin_m: float = 20e3,
    mode_index: int | None = None,
) -> tuple[float | None, dict[str, Any]]:
    """An upper bound on a facet echo's range at ``f_hz``, from the ionogram's trace.

    A facet the wave reaches lies above the reflection level of ``f_hz``, so
    its group path is shorter than the trace's virtual range there:
    ``int_0^D dz / mu`` with ``D`` above the reflection level is less than
    the same integral to it.  For each trace column (or only ``mode_index``)
    the virtual range at ``f_hz`` is interpolated linearly between the two
    detected columns either side, when they are at most two ionogram steps
    apart; the gate is the largest of these plus ``margin_m``, and None (no
    gate) when no column brackets ``f_hz``.  The largest, because the X trace
    reflects higher (a shorter range) and the tracker's O/X labels can swap
    near 5 MHz (seen on the Chapman scene, 2026-11-25): the O trace's range is
    the bound, and taking the larger of the two keeps a mislabel from
    tightening the gate.  It is a bound only when the O trace is among the
    bracketing columns: if only a true X trace brackets ``f_hz``, the gate is
    tighter than the O bound by the O-X separation less ``margin_m`` and can
    cut a facet near the reflection level.  ``report["columns_used"]`` says
    which columns gave it.
    """
    freq = np.asarray(getattr(ionogram, "freq_hz", []), dtype=float)
    vr = np.asarray(getattr(ionogram, "virtual_range_m", np.empty((0, 1))), dtype=float)
    rep: dict[str, Any] = {"f_hz": float(f_hz), "margin_m": float(margin_m), "mode_index": mode_index}
    if vr.ndim != 2 or vr.shape[0] != freq.size or freq.size < 2:
        rep["reason"] = "no usable trace column"
        return None, rep
    cols = range(vr.shape[1]) if mode_index is None else [int(mode_index)]
    grid = np.sort(freq[np.isfinite(freq)])
    kk = int(np.clip(np.searchsorted(grid, float(f_hz)), 1, grid.size - 1))
    step = float(grid[kk] - grid[kk - 1])  # the ionogram's step at f
    ranges: dict[int, float] = {}
    for c in cols:
        if not 0 <= c < vr.shape[1]:
            continue
        ok = np.isfinite(vr[:, c]) & np.isfinite(freq)
        f_ok, r_ok = freq[ok], vr[ok, c]
        order = np.argsort(f_ok)
        f_ok, r_ok = f_ok[order], r_ok[order]
        j = int(np.searchsorted(f_ok, float(f_hz)))
        if 0 < j < f_ok.size and f_ok[j] - f_ok[j - 1] <= 2.0 * step + 1.0:
            ranges[c] = float(np.interp(float(f_hz), f_ok[j - 1: j + 1], r_ok[j - 1: j + 1]))
    rep["virtual_range_by_column_m"] = {str(c): r for c, r in ranges.items()}
    if not ranges:
        rep["reason"] = "no trace column brackets the frequency without a gap"
        return None, rep
    r = max(ranges.values())
    rep.update({"virtual_range_m": r, "reason": "gated", "columns_used": sorted(ranges)})
    return r + float(margin_m), rep


@register("l3", "drift")
class DriftProcessing:
    """Route one at signal level on a drift-mode sequence.  See the module docstring."""

    layer = "l3"
    name = "drift"

    def __init__(
        self,
        *,
        prior: OnboardPrior | None = None,
        max_range_m: float = 1050e3,
        min_range_m: float = 20e3,
        max_echoes: int = 2,
        detect_db: float = 13.0,
        dynamic_db: float = 25.0,
        guard_samples: int = 3,
        zero_pad: int = 8,
        v_s_mps: float = 7600.0,
        earth_radius_m: float = RE,
        slope_tolerance: float = 0.2,
        min_contrast: float = 4.0,
        use_effective_speed: bool = True,
        detector: str = "windowed",
        window_sweeps: int = 21,
        pfa: float | None = None,
        pfa_per_sequence: float = 1.0e-2,
        estimator: str = "pinned",
        ambiguity_margin: float | None = None,
    ) -> None:
        if estimator not in ESTIMATORS:
            raise ValueError(f"estimator must be one of {ESTIMATORS}")
        self.estimator = estimator
        if ambiguity_margin is not None:
            if estimator != "joint_fading":
                raise ValueError("ambiguity_margin needs estimator='joint_fading'")
            if not float(ambiguity_margin) >= 0.0:
                raise ValueError("ambiguity_margin must not be negative")
        self.ambiguity_margin = None if ambiguity_margin is None else float(ambiguity_margin)
        if detector not in ("windowed", "summed"):
            raise ValueError("detector must be 'windowed' or 'summed'")
        if pfa is not None and not 0.0 < float(pfa) < 1.0:
            raise ValueError("pfa must be in (0, 1)")
        if not 0.0 < float(pfa_per_sequence) < 1.0:
            raise ValueError("pfa_per_sequence must be in (0, 1)")
        if int(window_sweeps) < 1:
            raise ValueError("window_sweeps must be at least 1")
        self.detector = detector
        self.window_sweeps = int(window_sweeps)
        self.pfa = None if pfa is None else float(pfa)
        self.pfa_per_sequence = float(pfa_per_sequence)
        self.last_detection: dict[str, Any] = {}
        self.prior = prior
        self.max_range_m = float(max_range_m)
        self.min_range_m = float(min_range_m)
        self.max_echoes = int(max_echoes)
        self.detect_db = float(detect_db)
        self.dynamic_db = float(dynamic_db)
        self.guard_samples = int(guard_samples)
        self.zero_pad = int(zero_pad)
        self.v_s_mps = float(v_s_mps)
        self.earth_radius_m = float(earth_radius_m)
        self.slope_tolerance = float(slope_tolerance)
        self.min_contrast = float(min_contrast)
        self.use_effective_speed = bool(use_effective_speed)

    def describe(self) -> Mapping[str, Any]:
        return {"name": self.name, "layer": self.layer, "max_range_m": self.max_range_m,
                "min_range_m": self.min_range_m, "max_echoes": self.max_echoes, "detect_db": self.detect_db, "dynamic_db": self.dynamic_db,
                "guard_samples": self.guard_samples, "zero_pad": self.zero_pad, "v_s_mps": self.v_s_mps,
                "earth_radius_m": self.earth_radius_m, "slope_tolerance": self.slope_tolerance,
                "min_contrast": self.min_contrast, "use_effective_speed": self.use_effective_speed,
                "detector": self.detector, "window_sweeps": self.window_sweeps, "pfa": self.pfa,
                "pfa_per_sequence": self.pfa_per_sequence, "estimator": self.estimator,
                "ambiguity_margin": self.ambiguity_margin,
                "prior_source": None if self.prior is None else self.prior.source}

    # -- steps 1 to 3 ----------------------------------------------------------------

    def _detect(self, x: np.ndarray, fs: float, k: float, max_range_m: float | None = None) -> list[float]:
        """Beats of the echoes, Hz, nearest first (step 1)."""
        if self.detector == "windowed":
            return self._detect_windowed(x, fs, k, self.max_range_m if max_range_m is None else max_range_m)
        self.last_detection = {"method": "summed"}
        n = x.shape[1]
        i0 = min(n - 8, int(math.ceil(2.0 * self.max_range_m / C_LIGHT * fs)))
        seg = x[:, i0:]
        m = seg.shape[1]
        nfft = int(2 ** math.ceil(math.log2(m * self.zero_pad)))
        win = np.hanning(m)
        p = np.sum(np.abs(np.fft.fft(seg * win, n=nfft, axis=1)) ** 2, axis=0)
        fr = np.fft.fftfreq(nfft, 1.0 / fs)
        lo, hi = k * 2.0 * self.min_range_m / C_LIGHT, k * 2.0 * self.max_range_m / C_LIGHT
        ok = np.nonzero((fr > lo) & (fr < hi))[0]
        if ok.size < 3:
            return []
        floor = float(np.median(p[ok]))
        thresh = floor * 10.0 ** (self.detect_db / 10.0)
        cell = max(3, int(round(4 * nfft / m)))  # four resolution cells either side
        found: list[float] = []
        p_top = float(p[ok].max())
        pp = p.copy()
        for _ in range(self.max_echoes):
            j = ok[int(np.argmax(pp[ok]))]
            if pp[j] < thresh or (found and pp[j] < p_top * 10.0 ** (-self.dynamic_db / 10.0)):
                break
            a, b, c = np.log(p[j - 1]), np.log(p[j]), np.log(p[j + 1])
            d = 0.5 * (a - c) / (a - 2.0 * b + c) if (a - 2.0 * b + c) != 0 else 0.0
            found.append(float(fr[j] + d * fs / nfft))
            pp[max(0, j - cell): j + cell + 1] = 0.0
        return sorted(found)

    def _detect_windowed(self, x: np.ndarray, fs: float, k: float, max_range_m: float) -> list[float]:
        """The windowed detector with a stated false-alarm probability (decision (j)).

        Each sweep's samples after ``2 max_range / c`` (so every echo nearer
        than ``max_range`` is present in all of them) are Hann windowed and
        transformed without padding; one bin is one range cell.  The power is
        summed over ``window_sweeps`` consecutive sweeps, sliding by one; a
        cell is a (window position, range bin) pair between ``min_range`` and
        ``max_range``.  The noise reference is the negative-beat half, where no
        echo falls (WP4's practice, :mod:`tss.l3_signal.detect`): its median
        over all sweeps and bins, over ``ln 2``, is the mean noise power per
        bin.  Under noise alone a cell is the sum of ``L`` independent
        exponential powers, so the threshold is the mean times the Gamma(L)
        quantile at ``1 - pfa_cell``.  ``pfa_cell`` is ``pfa`` when that is
        given (per cell, as WP4 states it); otherwise it is
        ``pfa_per_sequence`` over the number of cells in this sequence
        (decision (q)): a union bound on the chance that any cell of a
        noise-only sequence crosses, conservative because the sliding windows
        share sweeps.  A range gate leaves fewer cells and so a higher rate per
        cell.  Echoes are the strongest ``max_echoes`` range peaks above it, four bins apart
        and within ``dynamic_db`` of the strongest (the Hann sidelobe rule of
        the summed detector), then sorted by range.  The beat is refined by a
        parabola through the zero-padded power summed over the window that
        peaked.

        A noise-free record has no noise reference (its negative half holds
        only leakage of the window); its threshold is then far below the peak
        and the dynamic-range rule is the one that binds (``noise_free`` in the
        report: the threshold more than 60 dB under the strongest cell).
        """
        from scipy import ndimage
        from scipy.stats import gamma

        n_sw, n = x.shape
        i0 = min(n - 8, int(math.ceil(2.0 * float(max_range_m) / C_LIGHT * fs)))
        seg = x[:, i0:]
        m = seg.shape[1]
        win = np.hanning(m)
        p = np.abs(np.fft.fft(seg * win, axis=1)) ** 2
        fr = np.fft.fftfreq(m, 1.0 / fs)
        lo, hi = k * 2.0 * self.min_range_m / C_LIGHT, k * 2.0 * float(max_range_m) / C_LIGHT
        cells = np.nonzero((fr > lo) & (fr < hi))[0]
        neg = np.nonzero((fr < -lo) & (fr > -k * 2.0 * self.max_range_m / C_LIGHT))[0]
        L = min(self.window_sweeps, n_sw)
        report: dict[str, Any] = {"method": "windowed", "window_sweeps": L, "pfa_per_cell": self.pfa,
                                  "pfa_per_sequence": None if self.pfa is not None else self.pfa_per_sequence,
                                  "segment_samples": int(m), "max_range_m": float(max_range_m),
                                  "range_bins": int(cells.size)}
        self.last_detection = report
        if cells.size < 3 or neg.size < 3:
            report["n_cells"] = 0
            return []
        csum = np.cumsum(np.vstack([np.zeros((1, cells.size)), p[:, cells]]), axis=0)
        s = csum[L:] - csum[:-L]  # (window position, range bin)
        noise = float(np.median(p[:, neg])) / math.log(2.0)
        peak = float(s.max())
        pfa_cell = self.pfa if self.pfa is not None else self.pfa_per_sequence / float(s.size)
        report["pfa_per_cell"] = pfa_cell
        factor = float(gamma.isf(pfa_cell, L))
        thresh = noise * factor
        noise_free = bool(thresh * 1e6 < peak)
        report.update({"noise_power_per_bin": noise, "threshold_factor": factor, "noise_free": noise_free,
                       "n_cells": int(s.size), "noise_ref_bins": int(neg.size * n_sw)})
        over = s > thresh
        report["n_cells_over"] = int(over.sum())
        report["n_clusters_over"] = int(ndimage.label(over)[1])
        best = s.max(axis=0)
        bestj = s.argmax(axis=0)
        floor_dyn = peak * 10.0 ** (-self.dynamic_db / 10.0)
        cand = best.copy()
        picked: list[int] = []
        for _ in range(self.max_echoes):
            b = int(np.argmax(cand))
            v = float(cand[b])
            if v <= thresh or v < floor_dyn or v <= 0.0:
                break
            picked.append(b)
            cand[max(0, b - 4): b + 5] = 0.0
        report["range_bin_m"] = float(C_LIGHT * fs / m / (2.0 * k))
        beats = []
        pad = max(self.zero_pad, 1)
        nfft = int(2 ** math.ceil(math.log2(m * pad)))
        frp = np.fft.fftfreq(nfft, 1.0 / fs)
        for b in picked:
            j0 = int(bestj[b])
            pw = np.sum(np.abs(np.fft.fft(seg[j0: j0 + L] * win, n=nfft, axis=1)) ** 2, axis=0)
            centre = fr[cells[b]]
            near = np.nonzero(np.abs(frp - centre) <= 1.5 * fs / m)[0]
            jj = int(near[np.argmax(pw[near])])
            a, c0, c = np.log(pw[jj - 1]), np.log(pw[jj]), np.log(pw[jj + 1])
            den = a - 2.0 * c0 + c
            d = 0.5 * (a - c) / den if den != 0 else 0.0
            beats.append(float(frp[jj] + d * fs / nfft))
        # The per-echo entries in the order of the echoes (by range), not of picking.
        order = sorted(range(len(beats)), key=lambda i: beats[i])
        report["excess_db"] = [10.0 * math.log10(float(best[picked[i]]) / thresh) if thresh > 0 else math.inf
                               for i in order]
        report["window_start_sweep"] = [int(bestj[picked[i]]) for i in order]
        return [beats[i] for i in order]

    # -- the tier ---------------------------------------------------------------------

    def process_sequence(self, sequence: ReceivedSequence, *, prior: OnboardPrior | None = None,
                         max_range_m: float | None = None) -> DriftProduct:
        """Steps 1 to 5 on one sequence.

        ``max_range_m`` overrides the search's far edge for this sequence, as a
        range gate from the ionogram does (:func:`range_gate_from_ionogram`):
        the detector's segment then starts earlier and is longer.
        """
        prior = prior if prior is not None else self.prior
        recs = list(sequence.records)
        if not recs:
            raise ValueError("an empty sequence")
        sweep = recs[0].sweep
        fs = float(recs[0].fs_hz)
        k = float(sweep.k_f_hz_per_s)
        f0 = float(sweep.f_start_hz)
        n_min = min(r.samples.shape[1] for r in recs)
        x = np.stack([np.asarray(r.samples[0][:n_min], dtype=complex) for r in recs])
        t_n = np.asarray(sequence.t_s, dtype=float)
        tt = np.arange(n_min) / fs
        gate = self.max_range_m if max_range_m is None else min(float(max_range_m), self.max_range_m)
        self.last_detection = {}
        beats = self._detect(x, fs, k, gate) if self.detector == "windowed" else self._detect(x, fs, k)
        # The raw samples' noise variance, from the windowed detector's reference.
        sample_var = math.nan
        det = self.last_detection
        if det.get("noise_power_per_bin") is not None:
            sample_var = sample_noise_variance(det["noise_power_per_bin"], det["segment_samples"])
        not_done = ["tau_c_s", "doppler_spread_hz", "drift_sigma_mps",
                    "field_aligned", "several_frequency", "per_sweep_range_m", "per_sweep_power_db",
                    "per_sweep_doppler_hz", "snr_db", "resolution_mps", "v_unambiguous_mps",
                    "range_debiased_m", "los_rate_sigma_mps"]
        formed = () if self.estimator == "pinned" else ("drift_sigma_mps",)
        flags: dict[str, Any] = {
            "detection": dict(self.last_detection), "max_range_m": gate,
            "range_gate": max_range_m is not None,
            "prior_source": None if prior is None else prior.source,
            # What the prior was built from (onboard_prior_from_ionogram's report).
            "prior": {} if prior is None else {kk: vv for kk, vv in dict(prior.extra).items()
                                               if isinstance(vv, (str, int, float, type(None)))},
            "channel": 0, "n_sweeps": len(recs), "detected": len(beats),
            "estimator": self.estimator,
            "not_done": [n for n in not_done if n not in formed],
        }
        ids, rows = [], []
        for i, fb in enumerate(beats):
            ids.append(f"e{i + 1}")
            rows.append(self._echo(x, tt, t_n, fs, k, f0, fb, sequence.platforms, prior, sample_var))
        names = sorted({key for r in rows for key in r if key != "flags"})
        products: dict[str, np.ndarray] = {}
        for key in names:
            vals = [r.get(key, math.nan) for r in rows]
            products[key] = np.array(vals, dtype=object if any(isinstance(v, str) for v in vals) else float)
        flags["echo_flags"] = [r["flags"] for r in rows]
        return DriftProduct(t_s=t_n, echo_id=ids, products=products, flags=flags,
                            source=f"l3/drift ({len(ids)} echoes)")

    def _joint(self, times, field_e, est, field_var: float, **kw):
        """Step 6 on an echo the pinned test accepted: ``(products, flags)``."""
        start = {"t0_start_s": est.t_power_peak_s, "drift_start_mps": est.drift_mps}
        extra: dict[str, Any] = {"drift_pinned_mps": est.drift_mps, "t_pk_pinned_s": est.t_power_peak_s}
        flags: list[str] = []
        if self.estimator == "joint":
            fit = estimate_route_one_joint(times, field_e, **start, **kw)
        else:
            extra["field_noise_var"] = field_var
            if math.isfinite(field_var):
                # A noise-free record's reference is tiny; the fit floors it.
                fit = estimate_route_one_joint_fading(
                    times, field_e, noise_var=max(field_var, 1e-30), **start, **kw)
            else:
                fit = JointEstimate(False, "no noise reference (the windowed detector gives it)")
        extra["joint_ok"] = float(fit.ok)
        if not fit.ok:
            flags.append("JOINT_FAILED_PINNED_USED: " + fit.reason)
            return extra, flags
        extra.update({"drift_mps": fit.drift_mps, "drift_sigma_mps": fit.drift_sigma_mps,
                      "t_pk_s": fit.t0_s, "t_zc_s": fit.t_zc_s, "t0_sigma_s": fit.t0_sigma_s,
                      "width_s": fit.width_s, "corr_t0_v": fit.corr_t0_v})
        if self.estimator == "joint_fading":
            extra.update({"tau_c_fit_s": fit.tau_c_fit_s, "ambiguity_nll_margin": fit.nll_margin,
                          "drift_alt_mps": fit.drift_alt_mps})
            # The rule of docs/drift-mode.md section 4 for tau_c_s, on the fitted time.
            if fit.tau_c_fit_s > 0.5 * float(times[-1] - times[0]):
                flags.append("TAU_FIT_LOWER_BOUND")
            if fit.at_bound:
                flags.append("FIT_AT_BOUND: " + ", ".join(fit.at_bound))
            # A margin that could not be formed (NaN) is not a pass.
            if self.ambiguity_margin is not None and not fit.nll_margin >= self.ambiguity_margin:
                flags.append(f"DECLINED: ambiguous (a second maximum within {self.ambiguity_margin:g} "
                             f"in log likelihood: margin {fit.nll_margin:.2f}, drift {fit.drift_alt_mps:.1f} m/s)")
                for name in ("drift_mps", "drift_sigma_mps", "t_pk_s", "t_zc_s"):
                    extra[name] = math.nan
        return extra, flags

    def _echo(self, x, tt, t_n, fs, k, f0, fb, platforms: Sequence[Any], prior, sample_var: float):
        n = tt.size
        tau0 = fb / k
        i0 = min(n - 8, int(math.ceil(tau0 * fs)) + self.guard_samples)
        g = slice(i0, n)
        t_ref = float(np.mean(tt[g]))
        osc = np.exp(-2j * math.pi * fb * (tt[g] - t_ref))
        z = x[:, g] @ osc  # one complex value per sweep, phase +2 pi f tau
        f_eff = f0 + k * (t_ref - tau0)
        # The geometry time: received at t_ref, reflected half a delay earlier
        # (L2 module docstring, geometry time and receive time).
        times = t_n + t_ref - 0.5 * tau0
        echo_flags: list[str] = ["TILT_UNFLAGGED"]
        p0 = platforms[0] if len(platforms) else None
        alt_p = float(getattr(p0, "alt_m", math.nan))
        mu_sat = prior.mu_at(f_eff, alt_p) if prior is not None and math.isfinite(alt_p) else math.nan
        if not math.isfinite(mu_sat):
            echo_flags.append("NO_PRIOR" if prior is None else "EVANESCENT")
            return {"range_m": C_LIGHT * tau0 / 2.0, "f_eff_hz": f_eff, "blanked_fraction": i0 / n,
                    "echo_class": "unresolved", "flags": echo_flags}
        # Step 4: the platform term for a local-down wave normal.
        v_down = np.array([float(np.asarray(p.velocity_ned_mps)[2]) if getattr(p, "velocity_ned_mps", None)
                           is not None else 0.0 for p in platforms])
        fd_plat = 2.0 * f_eff * mu_sat * v_down / C_LIGHT
        ramp = 2.0 * math.pi * np.concatenate([[0.0], np.cumsum(0.5 * (fd_plat[1:] + fd_plat[:-1]) * np.diff(times))])
        field_e = np.conj(z) * np.exp(-1j * ramp)  # plate convention exp(-j 2 k0 P): closing -> +Doppler
        # Scale to a unit peak.  The envelope fit (scipy least_squares, absolute
        # gradient tolerance) stops early on physical-scale powers (a peak of
        # about -66 dB in the tier's units on the test scenes) and returns the
        # width it started from: with noise that put t_pk's spread at 66 ms
        # against the parameter-level 10 ms at 20 dB (2026-11-25).  A scale
        # changes nothing else route one forms.
        peak_power = float(np.max(np.abs(field_e) ** 2))
        if peak_power > 0.0:
            field_e = field_e / math.sqrt(peak_power)
        power = np.abs(field_e) ** 2
        ipk = int(np.argmax(power))
        # Apparent range at the peak, debiased for the Doppler (beat = K tau - f_D).
        step = field_e[1:] * np.conj(field_e[:-1])
        fd_local = float(np.angle(step[min(ipk, step.size - 1)]) / (2.0 * math.pi * np.diff(times).mean()))
        p_group = C_LIGHT * (fb + fd_local) / k / 2.0
        geometry = "flat"
        v = self.v_s_mps
        ratio = 1.0
        depth = math.nan
        p_pk = platforms[ipk] if len(platforms) > ipk else p0
        vel = getattr(p_pk, "velocity_ned_mps", None)
        if vel is not None:
            vv = np.asarray(vel, dtype=float)
            v = float(math.hypot(vv[0], vv[1]))
            geometry = "orbit"
            if self.use_effective_speed:
                r_p = self.earth_radius_m + float(p_pk.alt_m)
                ratio, depth = effective_speed_ratio(prior, f_eff, float(p_pk.alt_m), p_group, r_p)
                # The prior's group path cannot reach P': its reflection level
                # is above the echo, and the depth is that level, not the facet's.
                g_reach, _ = _group_integrals(prior, f_eff, float(p_pk.alt_m), depth, r_p)
                if not g_reach >= 0.99 * p_group:
                    echo_flags.append("PRIOR_REFLECTS_ABOVE")
            else:
                echo_flags.append("FLAT_SPEED")
        v_eff = v * ratio
        est = estimate_route_one(times, field_e, f_hz=f_eff, group_path_m=p_group, mu_sat=mu_sat,
                                 v_s_mps=v_eff, slope_tolerance=self.slope_tolerance,
                                 min_contrast=self.min_contrast)
        joint: dict[str, Any] = {}
        if not est.separable:
            echo_flags.append("DECLINED: " + est.reason)
        elif self.estimator != "pinned":
            # The per-sweep value sums the gated samples; then the unit-peak scaling.
            field_var = (n - i0) * sample_var / peak_power if peak_power > 0.0 else math.nan
            joint, joint_flags = self._joint(times, field_e, est, field_var, f_hz=f_eff,
                                             group_path_m=p_group, mu_sat=mu_sat, v_s_mps=v_eff)
            echo_flags.extend(joint_flags)
        t0 = joint.get("t_pk_s", est.t_power_peak_s)
        echo_class = "compact_facet" if est.separable else (
            "extended" if abs(est.slope_ratio - 1.0) > self.slope_tolerance else "unresolved")
        row = {
            "range_m": p_group, "range_beat_m": C_LIGHT * tau0 / 2.0, "f_eff_hz": f_eff, "mu_sat": mu_sat,
            "blanked_fraction": i0 / n, "t_ref_s": t_ref,
            "platform_doppler_hz": float(np.mean(fd_plat)), "los_rate_mps": est.los_rate_mps,
            "drift_mps": est.drift_mps, "drift_sigma_mps": math.nan, "t_pk_s": est.t_power_peak_s,
            "t_zc_s": est.t_zero_crossing_s, "slope_ratio": est.slope_ratio, "peak_contrast": est.peak_contrast,
            "along_track_m": v_eff * t0 if est.separable else math.nan,
            "v_mps": v, "v_eff_mps": v_eff, "facet_depth_m": depth, "geometry": geometry,
            "power_peak_db": 10.0 * math.log10(peak_power) if peak_power > 0 else -math.inf,
            "echo_class": echo_class, "flags": echo_flags,
        }
        row.update(joint)
        return row
