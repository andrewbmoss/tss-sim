"""Layer: L3
WP: WP4
Claims: C1
Status: implemented
Requirements: INS-04, INS-05, INS-10

Windowed transform of a de-chirped block, and the bin-index to apparent-range
mapping.

This is the module where the sub-watt argument of C1 becomes arithmetic, so the
argument is written out before the code.

Why the range cell bandwidth is the number that matters
-------------------------------------------------------
A pulse sounder that wants a range resolution ``dR`` must transmit a pulse of
bandwidth ``c / (2 dR)`` and its receiver must admit that whole bandwidth while
it listens.  For a 1 km cell that is 150 kHz of noise bandwidth, and at the
external noise levels of the topside HF band that sets the sensitivity.

An FMCW sounder gets the same range resolution from a sweep that is slow.  After
de-chirping, an echo at two-way range ``R`` sits at beat frequency
``f_b = 2 K_f R / c``, so a range cell of width ``dR`` occupies a beat-frequency
interval of

    B = 2 * K_f * dR / c

which is :func:`tss.constants.range_cell_noise_bandwidth`.  At the baseline
sweep rate ``K_f = 5.6 MHz/s`` a 1 km cell is about 37 Hz wide.  That is a
factor of about 4000, or 36 dB, less noise admitted per cell than the pulse
case at the same resolution, and it is the entire reason a transmitter under
1 W can close the link.  The swept bandwidth of about 20 MHz never appears in
the sensitivity calculation; it appears only in the unambiguous range and in
how long the sweep takes.

Two consistency statements are worth holding on to, because they are what the
unit tests check.

1.  The per-cell noise bandwidth equals the reciprocal of the coherent block
    duration.  A block of ``n`` samples at ``fs`` lasts ``T = n / fs`` and its
    Fourier bins are ``1 / T`` apart, which converts to a range spacing of
    ``c / (2 K_f T)``; feeding that spacing back through
    :func:`tss.constants.range_cell_noise_bandwidth` returns ``1 / T``.  The
    two ways of saying it agree exactly, and the FMCW sensitivity credit is
    just the statement that ``T`` is long.
2.  The mapping from bin index to range is the exact inverse of
    :func:`tss.constants.beat_frequency`.  No approximation enters between the
    two, so any range error in a result comes from the beat frequency being
    wrong (Doppler, sweep rate error, nonlinearity) and never from the
    transform.

Windows, resolution and bin spacing
-----------------------------------
Windowing is not optional here.  A topside ionogram has a strong leading-edge
echo and much weaker structure behind it, and the sidelobes of a rectangular
window would bury the weak structure under the strong echo's skirts.  The cost
is a wider main lobe.

The module keeps the two widths separate because conflating them is the usual
way a resolution claim gets overstated:

``bin_spacing_m``
    ``c fs / (2 K_f n_fft)``.  How finely the spectrum is sampled.  Zero
    padding reduces it without adding information.
``range_resolution_m``
    ``c / (2 K_f T) * enbw_bins``, where ``enbw_bins`` is the equivalent noise
    bandwidth of the window in bins and ``T`` is the unpadded block duration.
    This is the number INS-05 is scored against, and it is the one the noise
    bandwidth is computed from.

For a Hann window ``enbw_bins`` is 1.5, so a Hann-windowed block gives half
again the range cell width, and half again the noise bandwidth, that the
rectangular case would.  That is a real cost and it is reported rather than
absorbed.

Unambiguous range
-----------------
The highest beat frequency the sampled record can represent without folding is
``fs / 2``, which maps to ``c fs / (4 K_f)``.  INS-04 is scored against that
number.  It is a property of the sample rate and the sweep rate and nothing
else, which is why :func:`unambiguous_range_m` takes only those two.
"""

from __future__ import annotations

from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

import numpy as np

from tss.constants import C_LIGHT, range_cell_noise_bandwidth
from tss.l3_signal.base import range_from_beat_hz
from tss.l3_signal.dechirp import Dechirped

__all__ = [
    "RangeSpectrum",
    "WINDOWS",
    "window_taper",
    "window_gains",
    "transform",
    "bin_spacing_m",
    "unambiguous_range_m",
    "cell_noise_bandwidth_hz",
    "range_bin_index",
]

#: Window tapers implemented here.  They are written out in numpy rather than
#: taken from ``scipy.signal`` because this module may not import scipy at
#: module scope and a lazy import for four cosine coefficients is not worth the
#: indirection.  The periodic (``n / N``) convention is used throughout, which
#: is the correct one for spectral analysis; the symmetric (``n / (N-1)``)
#: convention that :func:`numpy.hanning` uses belongs to filter design and
#: leaks a little more.
WINDOWS: tuple[str, ...] = ("rect", "hann", "hamming", "blackman", "blackmanharris")

_COSINE_COEFFS: dict[str, tuple[float, ...]] = {
    "rect": (1.0,),
    "hann": (0.5, -0.5),
    "hamming": (0.54, -0.46),
    "blackman": (0.42, -0.5, 0.08),
    "blackmanharris": (0.35875, -0.48829, 0.14128, -0.01168),
}


@dataclass(frozen=True)
class RangeSpectrum:
    """One de-chirped block transformed into apparent range.

    Attributes
    ----------
    power_db
        Power spectrum in decibels, shape ``(n_channels, n_bins)``, scaled so
        that a unit-amplitude tone filling the block reads 0 dB whatever window
        was used.  The window's coherent gain is divided out, which is what
        makes a level comparable across window choices.
    range_m
        Apparent one-way range of each bin, metres.  Parallel to the last axis
        of ``power_db``.
    beat_hz
        Beat frequency of each bin, Hz.  Kept because a diagnostic that wants
        to look at the spectrum as a spectrum should not have to undo the range
        mapping.
    centre_frequency_hz
        Instantaneous transmitted frequency at the centre of the block.  This
        is the sounding frequency the resulting ionogram column is filed under.
    bin_spacing_m
        Range spacing between adjacent bins.  Zero padding shrinks this.
    range_resolution_m
        Width of one range cell from the window's equivalent noise bandwidth
        over the unpadded block.  This is what INS-05 is scored against.
    noise_bandwidth_hz
        Noise bandwidth admitted by one range cell, from
        :func:`tss.constants.range_cell_noise_bandwidth`.  Equals the
        reciprocal of the coherent block duration times ``enbw_bins``.
    unambiguous_range_m
        Range at the folding frequency ``fs / 2``.  INS-04.
    window
        Name of the taper applied.
    enbw_bins
        Equivalent noise bandwidth of that taper, in bins.
    negative_power_db
        The negative-beat half, same scaling, ordered by increasing ``|beat|``
        (index ``j`` is beat ``-(j + 1) fs / n_fft``).  No echo beats there, so
        it is the detector's noise reference (WP4 Phase 2,
        :mod:`tss.l3_signal.detect`).  ``None`` when not kept.
    """

    power_db: np.ndarray
    range_m: np.ndarray
    beat_hz: np.ndarray
    centre_frequency_hz: float
    bin_spacing_m: float
    range_resolution_m: float
    noise_bandwidth_hz: float
    unambiguous_range_m: float
    window: str
    enbw_bins: float
    negative_power_db: np.ndarray | None = None

    def describe(self) -> Mapping[str, Any]:
        """Small JSON-serialisable record of the transform, for diagnostics."""
        return {
            "n_channels": int(self.power_db.shape[0]),
            "n_bins": int(self.power_db.shape[-1]),
            "centre_frequency_hz": float(self.centre_frequency_hz),
            "bin_spacing_m": float(self.bin_spacing_m),
            "range_resolution_m": float(self.range_resolution_m),
            "noise_bandwidth_hz": float(self.noise_bandwidth_hz),
            "unambiguous_range_m": float(self.unambiguous_range_m),
            "window": str(self.window),
            "enbw_bins": float(self.enbw_bins),
        }


def window_taper(name: str, n: int) -> np.ndarray:
    """Return the named taper of length ``n``, periodic convention.

    Raises
    ------
    ValueError
        For an unknown name.  A silent fall back to a rectangular window would
        change the sidelobe floor by tens of decibels without saying so, which
        is exactly the sort of quiet substitution that makes a dynamic-range
        result untrustworthy.
    """
    key = str(name).lower()
    if key not in _COSINE_COEFFS:
        raise ValueError(f"unknown window {name!r}; expected one of {WINDOWS}")
    if int(n) < 1:
        raise ValueError(f"window length must be at least 1; got {n}")
    n = int(n)
    coeffs = _COSINE_COEFFS[key]
    k = np.arange(n, dtype=float)
    taper = np.zeros(n, dtype=float)
    for order, a in enumerate(coeffs):
        taper += a * np.cos(2.0 * np.pi * order * k / n)
    return taper


def window_gains(taper: np.ndarray) -> tuple[float, float]:
    """Coherent gain and equivalent noise bandwidth in bins for a taper.

    Coherent gain is ``sum(w) / n``: what the taper does to the amplitude of a
    tone sitting exactly on a bin.  Equivalent noise bandwidth is
    ``n sum(w^2) / (sum w)^2``: how many bins' worth of white noise one bin
    collects.  For a rectangular window both are 1.0 and 1.0; for Hann they are
    0.5 and 1.5.

    They are returned together because using one without the other is how a
    signal-to-noise figure acquires a factor of the window gain.
    """
    w = np.asarray(taper, dtype=float)
    n = w.size
    total = float(np.sum(w))
    if total == 0.0:
        raise ValueError("taper sums to zero, so its gains are undefined")
    coherent = total / n
    enbw = n * float(np.sum(w**2)) / (total**2)
    return coherent, enbw


def bin_spacing_m(fs_hz: float, n_fft: int, k_f_hz_per_s: float) -> float:
    """Range spacing between adjacent Fourier bins, metres.

    ``c fs / (2 K_f n_fft)``.  This is sampling, not resolution: padding the
    transform makes it smaller and adds no information.  Use
    :attr:`RangeSpectrum.range_resolution_m` for the number INS-05 wants.
    """
    return C_LIGHT * float(fs_hz) / (2.0 * abs(float(k_f_hz_per_s)) * int(n_fft))


def unambiguous_range_m(fs_hz: float, k_f_hz_per_s: float) -> float:
    """Apparent range at the folding beat frequency ``fs / 2``, metres.

    ``c fs / (4 K_f)``.  An echo beyond this folds back into the record and is
    reported at the wrong range with no indication that it did, which is why
    INS-04 is a requirement on the instrument and not a note in a docstring.
    """
    return C_LIGHT * float(fs_hz) / (4.0 * abs(float(k_f_hz_per_s)))


def cell_noise_bandwidth_hz(range_resolution_m: float, k_f_hz_per_s: float) -> float:
    """Noise bandwidth of one range cell, Hz.

    A thin wrapper on :func:`tss.constants.range_cell_noise_bandwidth`, kept
    here so that the number can be found from the transform that produces it.
    At the baseline ``K_f`` of 5.6 MHz/s a 1 km cell is about 37 Hz.  That
    number, rather than the roughly 20 MHz of swept bandwidth, is what enters
    the sensitivity budget, and it is why INS-10 closes at a sub-watt transmit
    power.
    """
    return float(range_cell_noise_bandwidth(float(range_resolution_m), abs(float(k_f_hz_per_s))))


def range_bin_index(
    range_m: float | np.ndarray,
    fs_hz: float,
    n_fft: int,
    k_f_hz_per_s: float,
) -> np.ndarray:
    """Fourier bin a given apparent range lands in, as a real number.

    ``k = 2 K_f R n_fft / (c fs)``.  Not rounded, because the fractional part
    is the bin straddling loss and a caller that rounds it away should do so
    deliberately.  This is the inverse of the ``range_m`` axis of
    :func:`transform` and the unit tests check the round trip.
    """
    r = np.asarray(range_m, dtype=float)
    return 2.0 * abs(float(k_f_hz_per_s)) * r * int(n_fft) / (C_LIGHT * float(fs_hz))


def transform(
    block: Dechirped,
    *,
    window: str = "hann",
    n_fft: int | None = None,
) -> RangeSpectrum:
    """Windowed FFT of a de-chirped block, mapped to apparent range.

    Parameters
    ----------
    block
        Output of :func:`tss.l3_signal.dechirp.dechirp`.  It must have been
        oriented, or the echo sits on the negative frequency axis and the
        positive half of the spectrum this function keeps is empty.  That is
        not enforced, because a caller looking at the raw mixer product has a
        legitimate reason to be here, but it is the first thing to check when a
        spectrum comes back empty.
    window
        One of :data:`WINDOWS`.
    n_fft
        Transform length.  The default is the block length.  A larger value
        zero pads, which reduces :attr:`RangeSpectrum.bin_spacing_m` and leaves
        :attr:`RangeSpectrum.range_resolution_m` alone.

    Returns
    -------
    RangeSpectrum
        A new frozen object.  ``block`` is not modified.

    Notes
    -----
    The non-negative beat frequencies, indices ``0`` to ``n_fft // 2``
    inclusive, are the range spectrum.  After the orientation applied by
    :func:`tss.l3_signal.dechirp.dechirp` an echo always has positive beat
    frequency, so the negative half carries only the image of the leakage, the
    noise and any carrier.  Since WP4 Phase 2 it is kept, as
    ``negative_power_db``, because that makes it the cleanest noise reference
    the detector has (:mod:`tss.l3_signal.detect`).
    """
    samples = np.atleast_2d(np.asarray(block.samples))
    n = int(samples.shape[1])
    if n < 2:
        raise ValueError(f"need at least two samples to transform; got {n}")
    length = n if n_fft is None else int(n_fft)
    if length < n:
        raise ValueError(
            f"n_fft={length} is shorter than the block ({n} samples); truncating a "
            "block silently changes the range resolution"
        )

    fs = float(block.fs_hz)
    # The block's own rate: for a stepped sweep the rate of its step, which
    # differs by a factor of about nine across the A1 band.  Falls back to the
    # schedule's rate for a block that does not carry one.
    local = getattr(block, "k_f_hz_per_s", None)
    k_f = abs(float(local if local is not None else block.sweep.k_f_hz_per_s))
    if k_f <= 0.0:
        raise ValueError("SweepSchedule.k_f_hz_per_s must be non-zero")

    taper = window_taper(window, n)
    coherent_gain, enbw = window_gains(taper)

    full = np.fft.fft(samples * taper[None, :], n=length, axis=-1)
    keep = length // 2 + 1
    spectrum = full[:, :keep]
    negative = full[:, keep:][:, ::-1]

    # Scale so a unit-amplitude tone filling the block reads 0 dB whatever the
    # window.  Dividing by n * coherent_gain is dividing by sum(w).
    amplitude = np.abs(spectrum) / (n * coherent_gain)
    with np.errstate(divide="ignore"):
        power_db = 20.0 * np.log10(np.where(amplitude > 0.0, amplitude, np.nan))
    power_db = np.where(np.isnan(power_db), -np.inf, power_db)
    neg_amp = np.abs(negative) / (n * coherent_gain)
    with np.errstate(divide="ignore"):
        neg_db = 20.0 * np.log10(np.where(neg_amp > 0.0, neg_amp, np.nan))
    neg_db = np.where(np.isnan(neg_db), -np.inf, neg_db)

    beat_hz = np.arange(keep, dtype=float) * fs / float(length)
    range_m = range_from_beat_hz(beat_hz, k_f)

    block_duration_s = n / fs
    resolution_m = enbw * C_LIGHT / (2.0 * k_f * block_duration_s)

    return RangeSpectrum(
        power_db=power_db,
        range_m=range_m,
        beat_hz=beat_hz,
        centre_frequency_hz=float(block.centre_frequency_hz),
        bin_spacing_m=bin_spacing_m(fs, length, k_f),
        range_resolution_m=float(resolution_m),
        noise_bandwidth_hz=cell_noise_bandwidth_hz(resolution_m, k_f),
        unambiguous_range_m=unambiguous_range_m(fs, k_f),
        window=str(window).lower(),
        enbw_bins=float(enbw),
        negative_power_db=neg_db,
    )
