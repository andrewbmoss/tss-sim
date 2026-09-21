"""Unit tests for L3, signal processing.

Two kinds of test live here and the difference matters.

The first kind checks a closed-form case: the beat-frequency round trip, the
noise bandwidth identity, and a synthetic tone landing in the bin the mapping
predicts.  These are where a sign error or a factor of two is actually caught.

The second kind checks that a stub is still a stub.  Every unimplemented
function in this layer has a test asserting it raises
:class:`NotImplementedError` with its work package named in the message.  That
is not ceremony.  A stub that quietly acquired a plausible-looking return value
would let a campaign produce a number that looks like evidence for C3 or C4,
and nothing else in the suite would notice.  These tests are what stops that.
"""

from __future__ import annotations

import numpy as np
import pytest

from tss.constants import C_LIGHT, beat_frequency, range_cell_noise_bandwidth
from tss.interfaces import Ionogram, ReceivedSignal, SweepSchedule
from tss.l3_signal import (
    BaselineSignalProcessing,
    MotionSeparatedSignalProcessing,
    dechirp,
    leading_edge,
    range_bin_index,
    range_from_beat_hz,
    robust_trace_extract,
    sweep_phase_rad,
    transform,
    vector,
    window_gains,
    window_taper,
)
from tss.l3_signal.novel import (
    aperture_synthesis,
    leakage_cancel,
    motion_compensation,
    sweep_pairing,
)
from tss.registry import available

K_F = 5.6e6
FS = 2.0e5


def _sweep(direction: int = 1, nonlinearity_ppm: float = 0.0) -> SweepSchedule:
    return SweepSchedule(
        f_start_hz=1.0e6,
        f_stop_hz=5.0e6,
        k_f_hz_per_s=K_F,
        sweep_period_s=1.0,
        direction=direction,
        nonlinearity_ppm=nonlinearity_ppm,
    )


def _point_target(range_m: float, n: int = 4096, noise: float = 0.0,
                  rng: np.random.Generator | None = None,
                  sweep: SweepSchedule | None = None) -> ReceivedSignal:
    """A single stationary echo at a known range, optionally with noise."""
    sweep = sweep or _sweep()
    t = np.arange(n) / FS
    tau = 2.0 * range_m / C_LIGHT
    samples = np.exp(1j * sweep_phase_rad(sweep, t - tau))
    if noise:
        generator = rng or np.random.default_rng(0)
        samples = samples + noise * (
            generator.normal(size=n) + 1j * generator.normal(size=n)
        )
    return ReceivedSignal(samples=samples[None, :], fs_hz=FS, sweep=sweep)


# ---------------------------------------------------------------------------
# Closed-form cases
# ---------------------------------------------------------------------------


@pytest.mark.parametrize("k_f", [1.0e6, 2.8e6, 5.6e6, 1.12e7])
def test_beat_frequency_round_trip_is_exact(k_f):
    """Range to beat frequency and back, to floating point.

    ``f_b = 2 K_f R / c`` and its inverse are exact, so the round trip is
    limited only by the arithmetic.  This mapping is what every apparent range
    in the package goes through.
    """
    ranges = np.array([1.0e3, 5.0e4, 3.0e5, 1.0e6, 2.5e6])
    back = range_from_beat_hz(beat_frequency(ranges, k_f), k_f)
    assert np.allclose(back, ranges, rtol=1.0e-13, atol=0.0)


def test_beat_frequency_is_linear_in_range():
    """Doubling the range doubles the beat frequency, exactly."""
    assert beat_frequency(2.0e5, K_F) == pytest.approx(2.0 * beat_frequency(1.0e5, K_F))


@pytest.mark.parametrize("window", ["rect", "hann", "hamming", "blackman"])
def test_range_cell_noise_bandwidth_is_enbw_over_block_duration(window):
    """One range cell admits ``enbw / T`` of noise bandwidth.

    This identity is the sub-watt argument of C1 written as arithmetic: what
    sets the noise per cell is the coherent block duration, not the swept
    bandwidth.
    """
    n, fs = 4096, FS
    duration = n / fs
    _, enbw = window_gains(window_taper(window, n))
    resolution = enbw * C_LIGHT / (2.0 * K_F * duration)
    assert range_cell_noise_bandwidth(resolution, K_F) == pytest.approx(
        enbw / duration, rel=1.0e-12
    )


def test_one_kilometre_cell_admits_about_37_hz():
    """The number quoted throughout the package and in the proposal."""
    assert range_cell_noise_bandwidth(1.0e3, 5.6e6) == pytest.approx(37.36, abs=0.05)


def test_hann_window_gains():
    """Coherent gain 0.5 and equivalent noise bandwidth 1.5, the textbook values."""
    coherent, enbw = window_gains(window_taper("hann", 1024))
    assert coherent == pytest.approx(0.5, rel=1.0e-9)
    assert enbw == pytest.approx(1.5, rel=1.0e-3)


def test_unknown_window_raises_rather_than_falling_back():
    """A misspelled window must not silently become rectangular."""
    with pytest.raises(ValueError, match="unknown window"):
        window_taper("hanning-ish", 64)


@pytest.mark.parametrize("range_m", [1.0e5, 4.0e5, 9.0e5])
def test_transform_puts_a_synthetic_tone_in_the_expected_bin(range_m):
    """A noiseless point target peaks where the mapping says it should.

    Rectangular window, no zero padding, no noise: the peak bin must be the one
    :func:`tss.l3_signal.transform.range_bin_index` predicts, to within the
    half bin the straddling allows.  This couples the de-chirp sign convention,
    the direction-dependent orientation and the bin-to-range mapping in one
    test, which is where those three most often disagree with each other.
    """
    n = 4096
    spectrum = transform(dechirp(_point_target(range_m, n=n)), window="rect", n_fft=n)
    peak = int(np.argmax(spectrum.power_db[0]))
    expected = float(range_bin_index(range_m, FS, n, K_F))
    assert abs(peak - expected) <= 0.5
    # Half a bin is the straddling error and it is an absolute tolerance in
    # metres, not a relative one: at short range half a bin is a large
    # fraction of the range and the mapping is no less exact for it.
    assert spectrum.range_m[peak] == pytest.approx(
        range_m, abs=0.5 * spectrum.bin_spacing_m
    )


def test_transform_level_is_window_independent():
    """A unit-amplitude tone reads 0 dB whatever window is applied.

    The coherent gain is divided out, which is what makes a level comparable
    across window choices.  Without it, switching to a Hann window would move
    every reported level by 6 dB and look like a sensitivity change.
    """
    block = dechirp(_point_target(4.0e5))
    for window in ("rect", "hann", "blackman"):
        spectrum = transform(block, window=window, n_fft=8192)
        assert float(np.max(spectrum.power_db[0])) == pytest.approx(0.0, abs=0.2)


def test_transform_reports_padding_as_spacing_not_resolution():
    """Zero padding reduces the bin spacing and leaves the resolution alone."""
    block = dechirp(_point_target(4.0e5, n=4096))
    plain = transform(block, window="hann", n_fft=4096)
    padded = transform(block, window="hann", n_fft=16384)
    assert padded.bin_spacing_m == pytest.approx(plain.bin_spacing_m / 4.0, rel=1.0e-9)
    assert padded.range_resolution_m == pytest.approx(plain.range_resolution_m, rel=1.0e-9)


def test_down_sweep_lands_at_positive_beat_frequency_too():
    """Orientation removes the direction asymmetry from the transform.

    On a down-sweep the raw mixer product sits on the opposite side of zero.
    :func:`tss.l3_signal.dechirp.dechirp` conjugates on an up-sweep so that an
    echo is at positive beat frequency either way, and getting that wrong is
    how the sign of the Doppler bias, which is C3, goes wrong.
    """
    n = 4096
    down = _sweep(direction=-1)
    spectrum = transform(dechirp(_point_target(4.0e5, n=n, sweep=down)),
                         window="rect", n_fft=n)
    peak = int(np.argmax(spectrum.power_db[0]))
    assert spectrum.range_m[peak] == pytest.approx(4.0e5, rel=2.0e-3)


# ---------------------------------------------------------------------------
# Trace extraction and the baseline tier
# ---------------------------------------------------------------------------


def test_leading_edge_reports_a_non_detection_rather_than_raising():
    """Noise alone gives NaN and a flag, not an exception.

    Contract rule 3 in :mod:`tss.interfaces`: a tier that cannot answer reports
    a masked result with a reason, because a frequency at which the wave was
    evanescent or penetrated is a fact about the ionosphere and not a bug.
    """
    rng = np.random.default_rng(3)
    n = 2048
    noise = (rng.normal(size=n) + 1j * rng.normal(size=n)) * 0.01
    signal = ReceivedSignal(samples=noise[None, :], fs_hz=FS, sweep=_sweep())
    sample = leading_edge(transform(dechirp(signal), window="hann"), threshold_db=15.0)
    assert np.isnan(sample.range_m)
    assert sample.flags.get("no_detection") is True


def test_leading_edge_is_at_or_before_the_peak():
    """The leading edge never sits behind the peak of the return.

    The inversion wants the range at which the echo starts.  The peak sits
    behind it by an amount that depends on the profile and the window, which is
    why the edge is what is extracted; this test is the ordering that makes
    that claim checkable.
    """
    spectrum = transform(dechirp(_point_target(4.0e5, noise=0.01)), window="hann",
                         n_fft=8192)
    sample = leading_edge(spectrum, threshold_db=10.0)
    assert np.isfinite(sample.range_m)
    assert sample.range_m <= sample.peak_range_m + 1.0e-9


def test_baseline_tier_forms_an_ionogram_and_is_registered():
    """The baseline tier registers and produces an ionogram of the right shape."""
    assert "baseline" in available("l3")["l3"]
    tier = BaselineSignalProcessing(n_frequencies=16, pad_factor=2)
    ionogram = tier.form_ionogram(_point_target(4.0e5, n=32768, noise=0.02))
    assert isinstance(ionogram, Ionogram)
    assert ionogram.freq_hz.shape == (16,)
    assert ionogram.virtual_range_m.shape == (16, 1)
    assert ionogram.snr_db.shape == (16, 1)
    detected = np.isfinite(ionogram.virtual_range_m[:, 0])
    assert detected.sum() >= 12
    assert np.allclose(ionogram.virtual_range_m[detected, 0], 4.0e5, rtol=0.03)


def test_baseline_tier_declines_to_claim_mode_or_doppler():
    """The baseline tier says what it did not do rather than saying zero.

    ``doppler_mps`` left ``None`` is the tier reporting that a single sweep
    cannot separate Doppler from range, which is C3.  A zero there would assert
    that it had separated them and found none.  ``mode_discrimination`` False
    is the same statement about INS-08.
    """
    tier = BaselineSignalProcessing(n_frequencies=8)
    ionogram = tier.form_ionogram(_point_target(4.0e5, n=16384))
    assert ionogram.doppler_mps is None
    assert ionogram.residual_range_bias_m is None
    assert ionogram.flags["mode_discrimination"] is False
    assert ionogram.flags["doppler_separated"] is False
    assert tuple(ionogram.modes) == ("unresolved",)


def test_baseline_tier_does_not_mutate_its_input():
    """Contract rule 1: a tier returns new objects and leaves the input alone."""
    signal = _point_target(4.0e5, n=16384)
    before = signal.samples.copy()
    BaselineSignalProcessing(n_frequencies=8).form_ionogram(signal)
    assert np.array_equal(signal.samples, before)


def test_baseline_tier_describe_round_trips_through_json():
    """Contract rule 5: ``describe()`` goes verbatim into the run manifest."""
    import json

    record = BaselineSignalProcessing().describe()
    assert json.loads(json.dumps(dict(record)))["name"] == "baseline"


# ---------------------------------------------------------------------------
# Stubs.  Every one of these asserts the work package is named.
# ---------------------------------------------------------------------------


def test_robust_trace_extract_is_a_stub_naming_wp4():
    with pytest.raises(NotImplementedError, match="WP4"):
        robust_trace_extract([])


@pytest.mark.parametrize(
    "call",
    [
        lambda: sweep_pairing.pair_sweeps(None, None),
        lambda: sweep_pairing.associate_traces(None, None),
        lambda: sweep_pairing.unbiased_range_m(np.zeros(1), np.zeros(1)),
        lambda: sweep_pairing.range_rate_from_pair_mps(
            np.zeros(1), np.zeros(1), np.ones(1), K_F
        ),
        lambda: sweep_pairing.residual_bias_cells(None, 1.0e4),
        lambda: sweep_pairing.SweepPair(),
    ],
)
def test_sweep_pairing_is_a_stub_naming_wp5(call):
    """C3 is unstarted research and every entry point says so.

    The half-sum arithmetic is stated in the module docstring and deliberately
    not offered as a function, because a function that averaged two range
    arrays would be used on unassociated traces and would return a number that
    looked like a measurement.
    """
    with pytest.raises(NotImplementedError, match="WP5"):
        call()


@pytest.mark.parametrize(
    "call",
    [
        lambda: motion_compensation.required_coherent_interval_s(50.0, 1.0e6),
        lambda: motion_compensation.phase_history_rad([], 1.0e6, np.zeros(3)),
        lambda: motion_compensation.compensate_motion(None, [], None),
        lambda: motion_compensation.integrate_coherently([], []),
    ],
)
def test_motion_compensation_is_a_stub_naming_wp5(call):
    with pytest.raises(NotImplementedError, match="WP5"):
        call()


@pytest.mark.parametrize(
    "call",
    [
        lambda: aperture_synthesis.along_track_aperture_m(7600.0, 1.0),
        lambda: aperture_synthesis.synthesize_aperture([], []),
        lambda: aperture_synthesis.along_track_resolution_m(1.0e4, 1.0e7, 4.0e5),
    ],
)
def test_aperture_synthesis_is_a_stub_naming_wp5(call):
    with pytest.raises(NotImplementedError, match="WP5"):
        call()


@pytest.mark.parametrize(
    "call",
    [
        lambda: leakage_cancel.estimate_leakage(None),
        lambda: leakage_cancel.cancel_leakage(None),
        lambda: leakage_cancel.residual_floor_db(None),
        lambda: leakage_cancel.cancellation_margin_db(0.0, 0.0),
    ],
)
def test_leakage_cancel_is_a_stub_naming_wp5(call):
    with pytest.raises(NotImplementedError, match="WP5"):
        call()


@pytest.mark.parametrize(
    "call",
    [
        lambda: vector.characteristic_polarizations(1.0e6, np.zeros(3), np.zeros(3)),
        lambda: vector.separate_modes(None, np.zeros(3), np.zeros(3)),
        lambda: vector.arrival_angle_deg(None, np.zeros(1)),
        lambda: vector.label_ionogram_modes(None, None, np.zeros(3), np.zeros(3)),
    ],
)
def test_vector_receive_is_a_stub_naming_wp11(call):
    with pytest.raises(NotImplementedError, match="WP11"):
        call()


def test_motion_separated_tier_registers_but_refuses_to_run():
    """The campaign shape is fixed before the algorithm exists.

    ``config/campaigns/c3_doppler.yaml`` names this tier, so it must register
    for the campaign to validate, and it must refuse when run.  Both halves
    are the point: ``tss config c3_doppler`` succeeds and ``tss run
    c3_doppler`` says WP5.
    """
    assert "motion_separated" in available("l3")["l3"]
    tier = MotionSeparatedSignalProcessing()
    assert tier.describe()["status"] == "stub"
    with pytest.raises(NotImplementedError, match="WP5"):
        tier.form_ionogram(_point_target(4.0e5))
