"""The delay of an echo follows its trace through the dwell.

An ionospheric echo's group delay depends on frequency, and the sweep moves
through frequency while a step is collected, so the de-chirped beat slides
instead of holding one tone.  These tests hold the L2 baseline to that: the
trace model is exact for what it claims, a flat trace changes nothing, a
sloping trace gives the slide and the smear the slope predicts, and the switch
that turns it off restores one delay per row.
"""

from __future__ import annotations

import numpy as np
import pytest

from tss import sweep_law as law
from tss.constants import C_LIGHT
from tss.interfaces import EchoSet
from tss.l2_instrument.baseline import BaselineInstrument, DelayTrace

A1 = law.stepped_schedule(3.3e5, 2.0e7, 0.0268)
SCALE_HEIGHT_M = 150.0e3
F_N_SAT_HZ = 2.0e6


def _topside_delay_s(f_hz):
    """Two-way group delay to reflection in an exponential topside.

    ``P' = H ln((1 + mu) / (1 - mu))`` with ``mu`` the refractive index at the
    spacecraft, so ``dP'/d ln f = 2 H / mu``.
    """
    mu = np.sqrt(1.0 - (F_N_SAT_HZ / np.asarray(f_hz, dtype=float)) ** 2)
    return 2.0 * SCALE_HEIGHT_M * np.log((1.0 + mu) / (1.0 - mu)) / C_LIGHT


def _echoes(delay_fn, amplitude_db=-60.0, lo=2.05e6, hi=8.0e6):
    centres = law.step_centres_hz(A1)
    f = centres[(centres > lo) & (centres < hi)]
    return EchoSet(
        freq_hz=f,
        group_delay_s=np.asarray(delay_fn(f), dtype=float),
        amplitude_db=np.full(f.size, amplitude_db),
        absorption_db=np.zeros(f.size),
        mode=np.array(["O"] * f.size),
    )


def _step_block(signal, f_hz, skip_s):
    """Samples of the step holding ``f_hz``, after ``skip_s``, with its rate."""
    fs = signal.fs_hz
    i = int(law.step_index(A1, float(law.time_at_frequency_s(A1, f_hz))))
    dwell = A1.step_dwell_s
    lo = int(np.ceil((i * dwell + skip_s) * fs))
    hi = int((i + 1) * dwell * fs)
    edges = law.step_edges_hz(A1)
    return signal.samples[0, lo:hi], (edges[i + 1] - edges[i]) / dwell, i, lo / fs


def test_trace_integral_is_exact_for_the_piecewise_linear_model():
    flat = DelayTrace([1.0e6, 2.0e6, 4.0e6], [3.0e-3, 3.0e-3, 3.0e-3])
    f = np.array([0.5e6, 1.0e6, 1.7e6, 4.0e6, 5.0e6])
    assert np.allclose(flat.delay_s(f), 3.0e-3)
    assert np.allclose(flat.phase_cycles(f), f * 3.0e-3)

    ramp = DelayTrace([1.0e6, 3.0e6], [1.0e-3, 2.0e-3])
    # tau = 1e-3 + 5e-10 (f - 1e6); its integral from 1e6 plus the anchor f0 tau0.
    d = f - 1.0e6
    assert np.allclose(ramp.delay_s(f), 1.0e-3 + 5.0e-10 * d)
    assert np.allclose(ramp.phase_cycles(f), 1.0e3 + 1.0e-3 * d + 2.5e-10 * d**2)


def test_trace_needs_distinct_frequencies():
    with pytest.raises(ValueError, match="distinct"):
        DelayTrace([1.0e6, 1.0e6], [1.0e-3, 2.0e-3])


def test_flat_trace_gives_the_undispersed_samples():
    """With no dispersion to model the two settings must agree sample for sample."""
    echoes = _echoes(lambda f: np.full(np.shape(f), 2.0 * 3.0e5 / C_LIGHT))
    on = BaselineInstrument(quantise=False).receive(echoes, A1, np.random.default_rng(7))
    off = BaselineInstrument(quantise=False, in_dwell_dispersion=False).receive(
        echoes, A1, np.random.default_rng(7)
    )
    assert on.flags["in_dwell_dispersion_modelled"] is True
    assert off.flags["in_dwell_dispersion_modelled"] is False
    scale = np.max(np.abs(off.samples))
    assert np.max(np.abs(on.samples - off.samples)) < 1.0e-6 * scale
    assert on.flags["max_in_dwell_delay_change_s"] == pytest.approx(0.0, abs=1e-12)


def test_beat_slides_with_the_trace():
    """The beat at any instant is the local rate times the delay arriving then."""
    echoes = _echoes(_topside_delay_s)
    signal = BaselineInstrument(quantise=False).receive(echoes, A1, np.random.default_rng(11))
    block, rate, i, t0 = _step_block(signal, 3.0e6, skip_s=1.1 * float(_topside_delay_s(3.2e6)))
    fs = signal.fs_hz
    n = block.size
    quarter = n // 4
    phase = np.unwrap(np.angle(block))
    t = np.arange(n) / fs

    def beat(sl):
        return np.polyfit(t[sl], phase[sl], 1)[0] / (2.0 * np.pi)

    def expected(sl):
        # Emission frequency at the middle of the slice: f(t) less the beat.
        t_mid = t0 + float(np.mean(t[sl]))
        f_now = float(law.frequency_hz(A1, t_mid))
        f_e = f_now
        for _ in range(8):
            f_e = f_now - rate * float(_topside_delay_s(f_e))
        return rate * float(_topside_delay_s(f_e))

    early, late = slice(0, quarter), slice(n - quarter, n)
    assert beat(early) == pytest.approx(expected(early), rel=0.01)
    assert beat(late) == pytest.approx(expected(late), rel=0.01)
    assert beat(late) > 1.03 * beat(early)


def test_dispersion_smears_the_range_spectrum_and_costs_peak():
    echoes = _echoes(_topside_delay_s)
    width_m, peak_db = {}, {}
    for on in (False, True):
        signal = BaselineInstrument(quantise=False, in_dwell_dispersion=on).receive(
            echoes, A1, np.random.default_rng(5)
        )
        block, rate, _, _ = _step_block(signal, 3.0e6, skip_s=1.1 * float(_topside_delay_s(3.2e6)))
        n = block.size
        spectrum = np.abs(np.fft.fft(block * np.hanning(n), 8 * n)) ** 2
        rng_m = C_LIGHT * np.fft.fftfreq(8 * n, 1.0 / signal.fs_hz) / (2.0 * rate)
        keep = (rng_m > 0.0) & (rng_m < 1.5e6)
        spectrum, rng_m = spectrum[keep], rng_m[keep]
        above = rng_m[spectrum > spectrum.max() / 2.0]
        width_m[on] = float(above.max() - above.min())
        peak_db[on] = 10.0 * np.log10(spectrum.max())
        if on:
            # 0.1 H / mu over a 5 per cent step: about 20 km at 3 MHz here.
            change_m = C_LIGHT * signal.flags["max_in_dwell_delay_change_s"] / 2.0
            assert change_m > 15.0e3
    assert width_m[True] > 3.0 * width_m[False]
    assert peak_db[False] - peak_db[True] > 4.0


def test_rows_without_a_single_valued_trace_keep_one_delay():
    """Two rays at each frequency in one mode: no unique trace, so no dispersion."""
    one = _echoes(_topside_delay_s, lo=2.9e6, hi=3.6e6)
    two = EchoSet(
        freq_hz=np.concatenate([one.freq_hz, one.freq_hz]),
        group_delay_s=np.concatenate([one.group_delay_s, 1.2 * one.group_delay_s]),
        amplitude_db=np.concatenate([one.amplitude_db, one.amplitude_db]),
        absorption_db=np.concatenate([one.absorption_db, one.absorption_db]),
        mode=np.array(["O"] * (2 * one.freq_hz.size)),
    )
    signal = BaselineInstrument(quantise=False).receive(two, A1, np.random.default_rng(2))
    assert signal.flags["in_dwell_dispersion_modelled"] is False
    assert signal.flags["rows_dispersed"] == 0
    assert signal.flags["rows_without_trace"] == signal.flags["rows_used"]


def test_modes_are_traced_separately():
    o = _echoes(_topside_delay_s, lo=2.9e6, hi=3.6e6)
    both = EchoSet(
        freq_hz=np.concatenate([o.freq_hz, o.freq_hz]),
        group_delay_s=np.concatenate([o.group_delay_s, 1.2 * o.group_delay_s]),
        amplitude_db=np.concatenate([o.amplitude_db, o.amplitude_db]),
        absorption_db=np.concatenate([o.absorption_db, o.absorption_db]),
        mode=np.array(["O"] * o.freq_hz.size + ["X"] * o.freq_hz.size),
    )
    signal = BaselineInstrument(quantise=False).receive(both, A1, np.random.default_rng(2))
    assert signal.flags["rows_dispersed"] == signal.flags["rows_used"]
    assert signal.flags["rows_without_trace"] == 0


def test_configuration_can_turn_it_off():
    assert BaselineInstrument.from_config({}).in_dwell_dispersion is True
    assert BaselineInstrument.from_config({"in_dwell_dispersion": False}).in_dwell_dispersion is False
    assert BaselineInstrument().describe()["in_dwell_dispersion"] is True
