"""Layer: L5
WP: WP6, WP7
Claims: C1, C2
Status: implemented
Requirements: INS-01, INS-02, INS-04, INS-05, INS-07, INS-08, INS-09, INS-10, SCI-09

Truth-versus-retrieval error metrics, and the achieved instrument numbers read
off the ionogram.

This module computes; it does not judge.  It returns a plain dictionary of
metrics keyed by name, and :mod:`tss.l5_assessment.requirements_check` maps
those names onto requirement IDs and decides verdicts.  The separation is
deliberate.  Changing what is measured and changing what counts as passing are
different changes, and a module that did both would make it impossible to see
which of the two moved when a verdict changed.

Every metric is either a number, a boolean, or ``None``.  ``None`` means the
metric could not be computed from what was supplied, and the requirements check
turns it into ``not_evaluated``.  There is no metric that defaults to a
plausible value, because a default would be scored as though it had been
measured.

The metrics
-----------

*Density accuracy.*  Fractional and absolute, both as an RMS and as a bias.
Both forms are reported because they answer different questions: the fractional
error is what the science requirement is written in, and the absolute error is
what shows whether the fractional error is being driven by the low-density top
of the profile where a small absolute error is a large fraction.

*Height error at matched density.*  For each retrieved point, the altitude at
which the truth profile has the same electron density, subtracted from the
altitude the retrieval put it at.  This is the metric that a topside inversion
should be judged on, because the inversion's job is to place a density at a
height, and comparing densities at a common height instead converts a height
error into a density error through whatever the local gradient happens to be.
Where the profile is steep a small height error looks like a large density
error, and where it is flat a large height error looks like none at all.  Both
are reported; the height error is the more honest of the two.

*Achieved instrument numbers.*  Range resolution, unambiguous range, ionogram
period, the realised frequency step law, the signal to noise ratio per range
cell, the lowest frequency at which an echo was obtained, and along-track
sampling.  These are read from the ionogram and its flags rather than from the
configuration, so that what is scored is what the processing achieved and not
what the campaign asked for.

A note on the along-track number
---------------------------------
``along_track_sampling_m`` is the ionogram period times the horizontal speed of
the platform.  Two things about it.  The horizontal speed is taken as the
magnitude of the first two components of ``PlatformState.velocity_ecef_mps``,
which assumes the caller has supplied local north-east-down components; that is
the same convention :func:`tss.l1_propagation.doppler.doppler_geometry` assumes
and it is stated rather than enforced, for the same reason.  And the spacing
that matters for SCI-09 is the spacing of the reflection points, which are
below the spacecraft and therefore move slightly slower, by the factor
``(R_E + h_r) / (R_E + h_s)``.  Using the spacecraft speed overstates the
spacing by a few per cent, which is the conservative direction, and the
correction is not applied because the reflection altitude varies across the
sweep and picking one would be arbitrary.

It is also sampling and not resolution.  The along-track cell is the antenna
footprint, which for an electrically small antenna with no directivity is very
much larger than the sampling interval.  Making the resolution match the
sampling is what :mod:`tss.l3_signal.novel.aperture_synthesis` would be for, and
it is exploratory.  SCI-09 as written asks for sampling, so sampling is what is
reported, and the gap between the two is recorded here rather than left
implicit.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import numpy as np

from tss.interfaces import DensityProfile, Ionogram, IonosphereState, PlatformState

__all__ = [
    "INS02_STEP_CROSSOVER_HZ",
    "compute_metrics",
    "profile_metrics",
    "ionogram_metrics",
    "height_error_at_matched_density_m",
    "along_track_sampling_m",
]

#: Frequency at which the INS-02 step law changes from the fractional rule to
#: the absolute cap: ``150 kHz / 0.05``.  Below it the 5 per cent fraction is
#: the binding constraint; above it the 150 kHz cap is.  The requirement's own
#: note names the 3 MHz crossover, and this is where that number comes from.
INS02_STEP_CROSSOVER_HZ = 150.0e3 / 0.05


def _finite_pairs(a: np.ndarray, b: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Rows where both arrays are finite."""
    a = np.asarray(a, dtype=float)
    b = np.asarray(b, dtype=float)
    keep = np.isfinite(a) & np.isfinite(b)
    return a[keep], b[keep]


def _topside_segment(truth: IonosphereState) -> tuple[np.ndarray, np.ndarray]:
    """The part of the truth profile above the density peak, altitude ascending.

    Restricting to the topside is what makes density a monotonic function of
    altitude, which is what lets a density be inverted back to a height.  Below
    the peak the same density occurs twice and the question "at what height does
    the truth have this density" has two answers.
    """
    alt = np.asarray(truth.alt_m, dtype=float)
    ne = np.asarray(truth.ne, dtype=float)
    if alt.ndim != 1 or ne.shape != alt.shape:
        raise ValueError("IonosphereState.alt_m and .ne must be parallel 1-D arrays")
    peak = int(np.argmax(ne))
    return alt[peak:], ne[peak:]


def height_error_at_matched_density_m(
    retrieved: DensityProfile,
    truth: IonosphereState,
) -> np.ndarray:
    """Retrieved altitude minus the truth altitude at the same density, metres.

    One entry per retrieved point whose density falls inside the range the
    truth topside covers.  Points outside that range are dropped rather than
    extrapolated, because the height at a density the truth never reaches is
    not defined.

    Returns an empty array when there is no overlap, which is a result and not
    an error.
    """
    alt_t, ne_t = _topside_segment(truth)
    if alt_t.size < 2:
        return np.empty(0, dtype=float)

    # Density falls with altitude on the topside, so reverse to get an
    # increasing x axis for the interpolation and work in log density, which is
    # close to linear in altitude for an exponential topside.
    good = np.isfinite(ne_t) & (ne_t > 0.0)
    alt_t, ne_t = alt_t[good], ne_t[good]
    if alt_t.size < 2:
        return np.empty(0, dtype=float)
    log_ne = np.log(ne_t)[::-1]
    alt_rev = alt_t[::-1]
    if not np.all(np.diff(log_ne) > 0.0):
        # Not strictly monotonic, usually because the grid undersamples the
        # peak.  Keep the longest strictly increasing run from the top.
        keep = np.concatenate(([True], np.diff(log_ne) > 0.0))
        log_ne, alt_rev = log_ne[keep], alt_rev[keep]
        if log_ne.size < 2:
            return np.empty(0, dtype=float)

    alt_r = np.asarray(retrieved.alt_m, dtype=float)
    ne_r = np.asarray(retrieved.ne, dtype=float)
    inside = (
        np.isfinite(alt_r)
        & np.isfinite(ne_r)
        & (ne_r > 0.0)
        & (np.log(np.where(ne_r > 0.0, ne_r, 1.0)) >= log_ne[0])
        & (np.log(np.where(ne_r > 0.0, ne_r, 1.0)) <= log_ne[-1])
    )
    if not np.any(inside):
        return np.empty(0, dtype=float)

    alt_truth_at = np.interp(np.log(ne_r[inside]), log_ne, alt_rev)
    return alt_r[inside] - alt_truth_at


def profile_metrics(
    retrieved: DensityProfile,
    truth: IonosphereState,
) -> dict[str, Any]:
    """Density and height error metrics for one retrieval against its truth.

    Densities are compared on the retrieved altitude grid, restricted to the
    topside the sounder can observe: from the truth density peak up to the
    highest altitude both profiles cover, which is the scope SW-01 of
    Requirements Rev G states.  Retrieved points below the peak are not
    scored, because a topside sounder has no echo from there and any such
    point is an artefact of the retrieval rather than a measurement.  Heights
    are compared at matched density.  Every entry is ``None`` when it could not
    be computed.
    """
    out: dict[str, Any] = {
        "ne_rms_fractional": None,
        "ne_bias_fractional": None,
        "ne_rms_absolute": None,
        "ne_bias_absolute": None,
        "height_error_rms_m": None,
        "height_error_bias_m": None,
        "n_points_compared": 0,
        "n_height_points_compared": 0,
        "overlap_min_alt_m": None,
        "overlap_max_alt_m": None,
        "n_points_below_peak": 0,
    }

    alt_r = np.asarray(retrieved.alt_m, dtype=float)
    ne_r = np.asarray(retrieved.ne, dtype=float)
    alt_t = np.asarray(truth.alt_m, dtype=float)
    ne_t = np.asarray(truth.ne, dtype=float)
    if alt_r.size == 0 or alt_t.size < 2:
        return out

    finite_t = np.isfinite(ne_t)
    if not np.any(finite_t):
        return out
    hm_t = float(alt_t[finite_t][int(np.argmax(ne_t[finite_t]))])
    lo = max(float(np.nanmin(alt_r)), float(np.min(alt_t)), hm_t)
    hi = min(float(np.nanmax(alt_r)), float(np.max(alt_t)))
    if not np.isfinite(lo) or not np.isfinite(hi) or hi <= lo:
        return out

    out["n_points_below_peak"] = int(np.count_nonzero(np.isfinite(alt_r) & (alt_r < hm_t)))
    inside = np.isfinite(alt_r) & np.isfinite(ne_r) & (alt_r >= lo) & (alt_r <= hi)
    a, r = alt_r[inside], ne_r[inside]
    if a.size == 0:
        return out

    t = np.interp(a, alt_t, ne_t)
    a, r = _finite_pairs(a, r)
    t = np.interp(a, alt_t, ne_t)
    usable = np.isfinite(t) & (t > 0.0)
    if not np.any(usable):
        return out
    r, t = r[usable], t[usable]

    absolute = r - t
    fractional = absolute / t
    out["ne_rms_fractional"] = float(np.sqrt(np.mean(fractional**2)))
    out["ne_bias_fractional"] = float(np.mean(fractional))
    out["ne_rms_absolute"] = float(np.sqrt(np.mean(absolute**2)))
    out["ne_bias_absolute"] = float(np.mean(absolute))
    out["n_points_compared"] = int(r.size)
    out["overlap_min_alt_m"] = float(lo)
    out["overlap_max_alt_m"] = float(hi)

    dh = height_error_at_matched_density_m(retrieved, truth)
    if dh.size:
        out["height_error_rms_m"] = float(np.sqrt(np.mean(dh**2)))
        out["height_error_bias_m"] = float(np.mean(dh))
        out["n_height_points_compared"] = int(dh.size)
    return out


def along_track_sampling_m(
    ionogram_period_s: float | None,
    platform: PlatformState | None,
) -> float | None:
    """Along-track spacing between successive ionograms, metres.

    ``None`` when either the period or the platform velocity is missing.  See
    the module docstring for the frame convention assumed and for why this is
    sampling and not resolution.
    """
    if ionogram_period_s is None or platform is None:
        return None
    v = platform.velocity_ecef_mps
    if v is None:
        return None
    v = np.asarray(v, dtype=float)
    if v.shape != (3,):
        return None
    horizontal = float(np.hypot(v[0], v[1]))
    if not np.isfinite(horizontal) or horizontal <= 0.0:
        return None
    return float(ionogram_period_s) * horizontal


def ionogram_metrics(
    ionogram: Ionogram,
    *,
    platform: PlatformState | None = None,
) -> dict[str, Any]:
    """Achieved instrument numbers read off an ionogram and its flags.

    Anything the L3 tier did not record is ``None`` here and
    ``not_evaluated`` downstream.

    Two pairs of band edges are reported and they mean different things.
    ``sweep_min_frequency_hz`` and ``sweep_max_frequency_hz`` are the band the
    transmitted sweep covered, which is the instrument capability INS-01a and
    INS-01b ask for.  ``min_sounding_frequency_hz`` and
    ``max_sounding_frequency_hz`` are the lowest and highest frequencies at
    which an echo was detected, which the ionosphere sets: the plasma frequency
    at the spacecraft and the X-mode penetration frequency.  The detected pair
    is a diagnostic and must not be scored against INS-01a/b.  Nothing in this function invents a value
    from the campaign configuration; what is scored is what the processing
    achieved.
    """
    flags: Mapping[str, Any] = ionogram.flags or {}
    freq = np.asarray(ionogram.freq_hz, dtype=float)
    vr = np.asarray(ionogram.virtual_range_m, dtype=float)
    snr = np.asarray(ionogram.snr_db, dtype=float)

    out: dict[str, Any] = {
        "virtual_range_resolution_m": _as_float(flags.get("range_resolution_m")),
        "unambiguous_range_m": _as_float(flags.get("unambiguous_range_m")),
        "ionogram_period_s": _as_float(flags.get("ionogram_period_s")),
        "coherent_interval_s": _as_float(flags.get("coherent_interval_s")),
        "block_duration_s": _as_float(flags.get("block_duration_s")),
        "noise_bandwidth_hz": _as_float(flags.get("noise_bandwidth_hz")),
        "min_sounding_frequency_hz": None,
        "max_sounding_frequency_hz": None,
        "sweep_min_frequency_hz": _as_float(flags.get("sweep_min_frequency_hz")),
        "sweep_max_frequency_hz": _as_float(flags.get("sweep_max_frequency_hz")),
        "frequency_step_fraction_max": None,
        "snr_db_min": None,
        "snr_db_median": None,
        "n_detected": None,
        "mode_discrimination": None,
        "doppler_resolved": bool(ionogram.doppler_mps is not None),
        "residual_range_bias_cells": None,
        "drift_resolution_mps": _as_float(flags.get("drift_resolution_mps")),
    }

    if "mode_discrimination" in flags:
        out["mode_discrimination"] = bool(flags["mode_discrimination"])

    detected = np.isfinite(vr).any(axis=1) if vr.ndim == 2 else np.isfinite(vr)
    detected &= np.isfinite(freq)
    out["n_detected"] = int(np.count_nonzero(detected))
    if np.any(detected):
        f_det = freq[detected]
        out["min_sounding_frequency_hz"] = float(np.min(f_det))
        out["max_sounding_frequency_hz"] = float(np.max(f_det))

    # INS-02 is scored as df <= min(0.05 f, 150 kHz).  Dividing the step by
    # min(f, crossover) turns both branches into one fraction, so the single
    # number can be compared against the 0.05 threshold directly.
    f_sorted = np.sort(freq[np.isfinite(freq)])
    if f_sorted.size >= 2:
        steps = np.diff(f_sorted)
        reference = np.minimum(f_sorted[:-1], INS02_STEP_CROSSOVER_HZ)
        with np.errstate(divide="ignore", invalid="ignore"):
            fraction = steps / reference
        fraction = fraction[np.isfinite(fraction)]
        if fraction.size:
            out["frequency_step_fraction_max"] = float(np.max(fraction))

    snr_values = snr[np.isfinite(snr)] if snr.size else np.empty(0)
    if snr_values.size:
        out["snr_db_min"] = float(np.min(snr_values))
        out["snr_db_median"] = float(np.median(snr_values))

    resolution = out["virtual_range_resolution_m"]
    if ionogram.residual_range_bias_m is not None and resolution:
        residual = np.asarray(ionogram.residual_range_bias_m, dtype=float)
        residual = residual[np.isfinite(residual)]
        if residual.size:
            out["residual_range_bias_cells"] = float(
                np.max(np.abs(residual)) / float(resolution)
            )

    out["along_track_sampling_m"] = along_track_sampling_m(
        out["ionogram_period_s"], platform
    )
    return out


def _as_float(value: Any) -> float | None:
    """Coerce a flag value to float, or ``None`` if it is absent or not a number."""
    if value is None or isinstance(value, bool):
        return None
    if type(value).__module__ == "numpy" and getattr(value, "dtype", None) is not None:
        if value.dtype == bool:
            # A numpy bool is not a Python bool and would otherwise coerce to
            # 1.0, turning a capability flag into a scored number.
            return None
    try:
        out = float(value)
    except (TypeError, ValueError):
        return None
    return out if np.isfinite(out) else None


def compute_metrics(
    retrieved: DensityProfile,
    truth: IonosphereState,
    ionogram: Ionogram,
    *,
    platform: PlatformState | None = None,
    extra: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """All metrics for one run, as a flat dictionary.

    ``extra`` merges in metrics computed elsewhere, which is how an L2-owned
    number such as the transmit power or the achieved leakage isolation reaches
    the requirements check without this module pretending to have measured it.
    Keys in ``extra`` override the computed ones, and a study that uses that to
    override a measured metric should say why in its own record.

    No input is modified.
    """
    metrics: dict[str, Any] = {}
    metrics.update(profile_metrics(retrieved, truth))
    metrics.update(ionogram_metrics(ionogram, platform=platform))
    if extra:
        metrics.update(dict(extra))
    return metrics
