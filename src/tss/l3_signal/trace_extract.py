"""Layer: L3
WP: WP4
Claims: C1, C2
Status: partial
Requirements: INS-05, INS-10

Detect the echo in a range spectrum and extract the leading edge of the trace.

What is delivered here is the easy half.  On a clean synthetic ionogram, where
each sounding frequency has one echo well above a flat noise floor, finding the
trace is a threshold and a connectivity rule and nothing more.  That is what
:func:`extract_trace` does, and it is enough to close the L3 to L4 to L5 chain
and to make the lamination inversion testable end to end.

What is not delivered is the half that a real ionogram needs, and the difference
between the two is worth being precise about, because a leading-edge extractor
that works on clean data and is quietly used on dirty data produces a profile
that is wrong and looks fine.

Why the leading edge and not the peak
--------------------------------------
The quantity the inversion of :mod:`tss.l4_retrieval.lamination` needs is the
group range at which the echo *starts*, because that is the range corresponding
to the reflection at the sounding frequency.  The peak of the return sits behind
the leading edge by an amount that depends on the shape of the profile near the
reflection height, on the window main lobe, and on how much the echo is spread
by the ionosphere.  Taking the peak therefore builds a profile-dependent bias
into the retrieval.  Taking the first cell over threshold builds in a
noise-dependent bias instead, which is smaller and, more importantly, can be
bounded from the signal to noise ratio.  Neither is unbiased.  The choice made
here is the leading edge, and the bias it carries is recorded in the returned
record rather than being left for the reader to discover.

The threshold
-------------
The noise floor is estimated as the median of the power in the cells beyond a
stated exclusion range, converted from median to mean with the factor for an
exponentially distributed power variate.  The median is used rather than the
mean because a single strong echo drags the mean up and raises the threshold,
which is the failure mode that loses the weakest traces.  The conversion factor
is ``1 / ln 2`` in power, or 1.59 dB, and it is applied so that the threshold is
stated relative to a mean noise power and can be compared against the link
budget prediction, which is also a mean.

A detection is declared where ``min_run`` consecutive cells exceed the floor by
``threshold_db``.  The run requirement is what stops a single noise spike from
starting a trace at the wrong range.  With a Hann window an echo occupies about
two cells at the -3 dB points, so ``min_run = 2`` is the smallest value that is
not simply the single-cell test again.

What the robust version has to do, and why it is a stub
--------------------------------------------------------
:func:`robust_trace_extract` is the WP4 deliverable and it is not written.  The
problems it has to solve, none of which the threshold-and-connect version
addresses:

* **The leakage skirt.**  The transmit-to-receive leakage sits at essentially
  zero range with tens of decibels more power than any echo, and its window
  skirt extends across the near range cells.  The exclusion range used here is a
  blunt instrument that throws away the near range along with the leakage.  The
  real answer is to cancel the leakage first, which is
  :mod:`tss.l3_signal.novel.leakage_cancel`, and until that exists the
  near-range behaviour of this extractor is not meaningful.
* **Multiple traces.**  O and X are two traces, the X displaced upward in
  frequency by roughly half the gyrofrequency, and a scalar receiver sees them
  superposed.  Labelling them needs either polarization (WP11, C5) or a model of
  their separation, and picking the first edge without labelling silently mixes
  them near the cusp.
* **Spread F and multiple hop.**  A spread trace has no well defined leading
  edge at all, and a second-hop return looks like a legitimate echo at twice the
  range.
* **The trace cusp.**  Near the critical frequency the group range rises
  steeply and the leading-edge estimate becomes very sensitive to the threshold.
  This is where most of the retrieval error in a real inversion comes from, and
  it is exactly where a fixed threshold performs worst.
* **A stated detection statistic.**  A fixed threshold in decibels has no stated
  false-alarm rate.  A CFAR detector with a cell-averaging or ordered-statistic
  estimator does, and INS-10 should be scored against a detector whose
  probability of false alarm is a number rather than a habit.

Until those are addressed, results from this module are valid for the clean
synthetic case and for nothing else, and a study that uses it should say so.
"""

from __future__ import annotations

import math
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

import numpy as np

from tss.l3_signal.transform import RangeSpectrum

__all__ = [
    "TraceSample",
    "MEDIAN_TO_MEAN_DB",
    "noise_floor_db",
    "detect_mask",
    "leading_edge",
    "extract_trace",
    "connect_trace",
    "robust_trace_extract",
]

#: Decibels to add to a median power estimate to get the mean, for an
#: exponentially distributed power variate: ``10 log10(1 / ln 2)`` = 1.592 dB.
#: Applied so the estimated floor is a mean and is comparable with the L2 link
#: budget, which predicts a mean.
MEDIAN_TO_MEAN_DB = 10.0 * math.log10(1.0 / math.log(2.0))

_ROBUST_MESSAGE = (
    "WP4: the robust trace extractor is not implemented.  The threshold-and-connect "
    "extractor in tss.l3_signal.trace_extract.extract_trace works on clean synthetic "
    "data only.  A robust extractor needs leakage cancellation "
    "(tss.l3_signal.novel.leakage_cancel, WP5), a stated detection statistic with a "
    "false-alarm rate rather than a fixed threshold, O and X labelling (WP11, C5), "
    "and a cusp treatment.  See the module docstring and docs/work-packages.md row WP4."
)


@dataclass(frozen=True)
class TraceSample:
    """One sounding frequency's extracted trace point.

    Attributes
    ----------
    freq_hz
        Sounding frequency this sample belongs to.
    range_m
        Apparent range of the leading edge, metres, or NaN when nothing was
        detected.
    snr_db
        Peak power in the detected run above the estimated noise floor, or NaN.
    peak_range_m
        Range of the peak of the detected run, for readers who want to see how
        far behind the leading edge it sits.  NaN when nothing was detected.
    n_cells
        Length of the detected run, in cells.  Zero when nothing was detected.
    noise_floor_db
        The estimated mean noise floor the threshold was taken from.
    flags
        Free-form record of why a sample is what it is, notably
        ``"no_detection"`` and ``"at_exclusion_edge"``.
    """

    freq_hz: float
    range_m: float
    snr_db: float
    peak_range_m: float
    n_cells: int
    noise_floor_db: float
    flags: Mapping[str, Any]


def noise_floor_db(
    power_db: np.ndarray,
    range_m: np.ndarray,
    *,
    exclude_below_m: float = 0.0,
    exclude_above_m: float | None = None,
) -> float:
    """Estimate the mean noise floor of one spectrum, dB.

    Median of the power in the cells outside the exclusion window, plus
    :data:`MEDIAN_TO_MEAN_DB`.  The median is taken in decibels, which is the
    median in power as well because the conversion is monotonic.

    ``exclude_below_m`` removes the near range where the leakage and its window
    skirt live.  Excluding it is a blunt substitute for cancelling it; see the
    module docstring.

    Returns ``-inf`` when every usable cell is empty, which is a result and not
    an error: a block that contains no signal and no noise has no floor.
    """
    p = np.asarray(power_db, dtype=float)
    r = np.asarray(range_m, dtype=float)
    usable = r >= float(exclude_below_m)
    if exclude_above_m is not None:
        usable &= r <= float(exclude_above_m)
    usable &= np.isfinite(p)
    if not np.any(usable):
        return -np.inf
    return float(np.median(p[usable])) + MEDIAN_TO_MEAN_DB


def detect_mask(
    power_db: np.ndarray,
    range_m: np.ndarray,
    floor_db: float,
    *,
    threshold_db: float = 10.0,
    exclude_below_m: float = 0.0,
) -> np.ndarray:
    """Boolean mask of cells exceeding the floor by ``threshold_db``.

    Cells inside the near-range exclusion are never marked, whatever their
    power, because the leakage is there and it is not an echo.
    """
    p = np.asarray(power_db, dtype=float)
    r = np.asarray(range_m, dtype=float)
    mask = np.isfinite(p) & (p >= floor_db + float(threshold_db))
    mask &= r >= float(exclude_below_m)
    return mask


def _runs(mask: np.ndarray) -> list[tuple[int, int]]:
    """Contiguous ``True`` runs as ``(start, stop)`` half-open index pairs."""
    m = np.asarray(mask, dtype=bool)
    if m.size == 0:
        return []
    edges = np.diff(m.astype(np.int8))
    starts = list(np.flatnonzero(edges == 1) + 1)
    stops = list(np.flatnonzero(edges == -1) + 1)
    if m[0]:
        starts.insert(0, 0)
    if m[-1]:
        stops.append(m.size)
    # strict=False: the two lists are built to the same length by
    # construction, and pairing them is the whole point.
    return list(zip(starts, stops, strict=False))


def leading_edge(
    spectrum: RangeSpectrum,
    *,
    channel: int = 0,
    threshold_db: float = 10.0,
    min_run: int = 2,
    exclude_below_m: float = 0.0,
    exclude_above_m: float | None = None,
) -> TraceSample:
    """Leading-edge range and signal to noise ratio for one range spectrum.

    Finds the first contiguous run of at least ``min_run`` cells that exceeds
    the estimated noise floor by ``threshold_db``, and reports the range of the
    first cell of that run.  Returns a sample with NaN range and
    ``flags["no_detection"]`` when there is no such run.  A non-detection is a
    result: the wave was evanescent, or it penetrated, or there was not enough
    signal, and contract rule 3 in :mod:`tss.interfaces` says a tier reports
    that rather than raising.

    The range reported is the range of the first cell over threshold, not an
    interpolated edge.  Interpolating the edge would need a model of the echo
    shape, and on clean data the quantisation error is half a bin, which is
    already recorded in :attr:`RangeSpectrum.bin_spacing_m`.
    """
    power = np.asarray(spectrum.power_db, dtype=float)
    if power.ndim != 2:
        raise ValueError(f"RangeSpectrum.power_db must be 2-D; got {power.shape}")
    if not 0 <= int(channel) < power.shape[0]:
        raise ValueError(f"channel {channel} outside 0..{power.shape[0] - 1}")
    row = power[int(channel)]
    rng = np.asarray(spectrum.range_m, dtype=float)

    floor = noise_floor_db(
        row, rng, exclude_below_m=exclude_below_m, exclude_above_m=exclude_above_m
    )
    mask = detect_mask(
        row, rng, floor, threshold_db=threshold_db, exclude_below_m=exclude_below_m
    )
    if exclude_above_m is not None:
        mask &= rng <= float(exclude_above_m)

    for start, stop in _runs(mask):
        if stop - start < int(min_run):
            continue
        segment = row[start:stop]
        peak_offset = int(np.argmax(segment))
        flags: dict[str, Any] = {}
        if exclude_below_m > 0.0 and start > 0 and rng[start - 1] < float(exclude_below_m):
            flags["at_exclusion_edge"] = True
        return TraceSample(
            freq_hz=float(spectrum.centre_frequency_hz),
            range_m=float(rng[start]),
            snr_db=float(segment[peak_offset] - floor),
            peak_range_m=float(rng[start + peak_offset]),
            n_cells=int(stop - start),
            noise_floor_db=float(floor),
            flags=flags,
        )

    return TraceSample(
        freq_hz=float(spectrum.centre_frequency_hz),
        range_m=float("nan"),
        snr_db=float("nan"),
        peak_range_m=float("nan"),
        n_cells=0,
        noise_floor_db=float(floor),
        flags={"no_detection": True},
    )


def connect_trace(
    samples: list[TraceSample],
    *,
    max_jump_m: float | None = None,
) -> list[TraceSample]:
    """Reject samples that jump too far from their neighbours.

    A topside trace is continuous in range against frequency except at the cusp,
    so a point that sits far from both of its neighbours is almost always a
    noise spike or a second-hop return rather than the trace.  This drops such
    points by setting their range to NaN and flagging them, rather than removing
    them from the list, so that the frequency axis of the ionogram stays the
    frequency axis the sweep visited.

    ``max_jump_m`` of ``None`` disables the rule and returns equivalent samples
    unchanged.  The default is off because a sensible value depends on the
    frequency step and on how close to the cusp the sweep runs, and a guessed
    value would quietly delete real cusp points.  New objects are returned; the
    input list and its members are not modified.
    """
    out = list(samples)
    if max_jump_m is None or len(out) < 3:
        return out

    ranges = np.array([s.range_m for s in out], dtype=float)
    limit = float(max_jump_m)
    keep = np.ones(ranges.size, dtype=bool)
    for i in range(1, ranges.size - 1):
        here, before, after = ranges[i], ranges[i - 1], ranges[i + 1]
        if not np.isfinite(here) or not np.isfinite(before) or not np.isfinite(after):
            continue
        if abs(here - before) > limit and abs(here - after) > limit:
            keep[i] = False

    for i, ok in enumerate(keep):
        if ok:
            continue
        s = out[i]
        flags = dict(s.flags)
        flags["rejected_jump"] = True
        out[i] = TraceSample(
            freq_hz=s.freq_hz,
            range_m=float("nan"),
            snr_db=float("nan"),
            peak_range_m=float("nan"),
            n_cells=0,
            noise_floor_db=s.noise_floor_db,
            flags=flags,
        )
    return out


def extract_trace(
    spectra: list[RangeSpectrum],
    *,
    channel: int = 0,
    threshold_db: float = 10.0,
    min_run: int = 2,
    exclude_below_m: float = 0.0,
    exclude_above_m: float | None = None,
    max_jump_m: float | None = None,
) -> list[TraceSample]:
    """Leading-edge extraction across a stack of range spectra.

    One spectrum per sounding frequency, in sweep order.  Returns one
    :class:`TraceSample` per spectrum, in the same order, with NaN ranges where
    nothing was detected.

    Valid for clean synthetic data.  See the module docstring for what it does
    not handle and :func:`robust_trace_extract` for where that work lives.
    """
    samples = [
        leading_edge(
            s,
            channel=channel,
            threshold_db=threshold_db,
            min_run=min_run,
            exclude_below_m=exclude_below_m,
            exclude_above_m=exclude_above_m,
        )
        for s in spectra
    ]
    return connect_trace(samples, max_jump_m=max_jump_m)


def robust_trace_extract(
    spectra: list[RangeSpectrum],
    *,
    channel: int = 0,
    pfa: float = 1.0e-4,
    gyrofrequency_hz: float | None = None,
) -> list[TraceSample]:
    """Trace extraction for real ionograms.  Not implemented.

    Intended approach: an ordered-statistic CFAR detector with the stated
    probability of false alarm ``pfa`` replacing the fixed threshold, run after
    leakage cancellation rather than after a near-range exclusion, with O and X
    hypotheses tracked separately using ``gyrofrequency_hz`` to predict their
    separation, and with an explicit cusp model so that the steep part of the
    trace is fitted rather than thresholded.

    The open question is what to do when the hypotheses cannot be separated.
    A scalar receiver near the cusp sees one blended return, and the choice is
    between rejecting the frequency, which loses exactly the part of the trace
    that constrains the peak, and accepting a blended range, which biases the
    inversion.  That choice has to be made against a measured error budget and
    the budget does not exist yet.

    Raises
    ------
    NotImplementedError
        Always.  WP4.
    """
    del spectra, channel, pfa, gyrofrequency_hz
    raise NotImplementedError(_ROBUST_MESSAGE)
