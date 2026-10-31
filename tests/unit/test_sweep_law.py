"""Unit tests for the stepped A1 sweep law, :mod:`tss.sweep_law`.

The numbers checked here are the ones Requirements Rev G states for the A1
ionogram sweep: 159 steps over 0.33 to 20 MHz, 26.8 ms each, 4.26 s, peak rate
5.6 MHz/s, 211 steps at the 0.025 goal.  The end-to-end test is the one that
matters most: an echo below 3 MHz, where the step rate is up to nine times
slower than the peak, has to land at its true range.  A processor that used the
peak rate for every block would put it at a ninth of that range.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from tss import sweep_law as law
from tss.constants import C_LIGHT
from tss.interfaces import EchoSet, SweepSchedule
from tss.l2_instrument.waveform import (
    schedule_from_config,
    sweep_duration_s,
    sweep_phase_rad,
    validate_schedule,
)
from tss.l3_signal import BaselineSignalProcessing
from tss.l3_signal.base import instantaneous_frequency_hz
from tss.l3_signal.base import sweep_phase_rad as l3_phase

A1 = law.stepped_schedule(3.3e5, 2.0e7, 0.0268)


def test_a1_matches_requirements_rev_g():
    assert law.step_count(A1) == 159
    assert law.active_duration_s(A1) == pytest.approx(4.26, abs=0.005)
    assert A1.k_f_hz_per_s == pytest.approx(5.6e6, rel=0.001)
    assert sweep_duration_s(A1) == pytest.approx(law.active_duration_s(A1))
    goal = law.stepped_schedule(3.3e5, 2.0e7, 0.0268, step_fraction=0.025)
    assert law.step_count(goal) == 211


def test_step_widths_follow_ins02():
    e = law.step_edges_hz(A1)
    w = np.diff(e)
    assert e[0] == 3.3e5 and e[-1] == 2.0e7
    # Every step obeys min(0.05 f, 150 kHz) at its lower edge; the last may be
    # shorter because it is clipped at f_stop.
    assert np.all(w <= np.minimum(0.05 * e[:-1], 150.0e3) * (1 + 1e-12))
    assert np.allclose(w[:-1], np.minimum(0.05 * e[:-2], 150.0e3))


def test_bottom_step_range_cell_is_9_1_km():
    """c / (2 df) at the first step, the Rev G range cell at 0.33 MHz."""
    df = law.step_widths_hz(A1)[0]
    assert C_LIGHT / (2.0 * df) == pytest.approx(9.08e3, rel=0.005)


def test_frequency_is_continuous_and_phase_is_its_integral():
    t = np.linspace(0.0, law.active_duration_s(A1) - 1e-9, 200_001)
    f = law.frequency_hz(A1, t)
    assert np.all(np.diff(f) >= 0.0)
    assert f[0] == pytest.approx(3.3e5)
    assert f[-1] == pytest.approx(2.0e7, rel=1e-6)
    phi = law.phase_rad(A1, t)
    numeric = np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(t))])
    assert np.allclose(phi / (2 * np.pi), numeric, rtol=0, atol=1e-2)


def test_nonlinearity_keeps_step_edges_and_phase_exact():
    s = replace(A1, nonlinearity_ppm=500.0)
    edges_t = np.arange(law.step_count(s)) * s.step_dwell_s
    assert np.allclose(law.frequency_hz(s, edges_t), law.step_edges_hz(s)[:-1])
    t = np.linspace(0.0, 0.2, 400_001)
    f = law.frequency_hz(s, t)
    numeric = np.concatenate([[0.0], np.cumsum(0.5 * (f[1:] + f[:-1]) * np.diff(t))])
    assert np.allclose(law.phase_rad(s, t) / (2 * np.pi), numeric, atol=1e-3)


def test_down_sweep_visits_the_same_steps_in_reverse():
    down = replace(A1, direction=-1)
    total = law.active_duration_s(A1)
    # Offset from the step edges so that no sample sits on a boundary.
    t = (np.arange(5000) + 0.37) * total / 5000.0
    assert np.allclose(law.frequency_hz(down, total - t), law.frequency_hz(A1, t), rtol=1e-9)
    assert law.frequency_hz(down, 0.0) == pytest.approx(2.0e7)


def test_time_at_frequency_inverts_the_law():
    f = np.geomspace(3.4e5, 1.99e7, 57)
    for sweep in (A1, replace(A1, direction=-1)):
        t = law.time_at_frequency_s(sweep, f)
        assert np.allclose(law.frequency_hz(sweep, t), f, rtol=1e-9)


def test_rate_is_local_and_slow_at_the_bottom():
    t_bottom = 0.5 * A1.step_dwell_s
    t_top = law.active_duration_s(A1) - 1.5 * A1.step_dwell_s
    assert law.rate_hz_per_s(A1, t_bottom) == pytest.approx(0.05 * 3.3e5 / 0.0268)
    assert law.rate_hz_per_s(A1, t_top) == pytest.approx(A1.k_f_hz_per_s)


def test_l2_and_l3_helpers_share_the_stepped_law():
    t = np.linspace(0.0, 1.0, 1001)
    assert np.allclose(sweep_phase_rad(A1, t), l3_phase(A1, t))
    assert np.allclose(instantaneous_frequency_hz(A1, t), law.frequency_hz(A1, t))


def test_validate_refuses_a_rate_that_disagrees_with_the_steps():
    bad = replace(A1, k_f_hz_per_s=2.0e6)
    with pytest.raises(ValueError, match="peak step rate"):
        validate_schedule(bad)


def test_config_builds_a1_and_checks_a_stated_rate():
    cfg = {"f_start_hz": 3.3e5, "f_stop_hz": 2.0e7, "step_dwell_s": 0.0268,
           "step_fraction": 0.05, "step_cap_hz": 1.5e5, "sweep_period_s": 0.0}
    s = schedule_from_config(cfg)
    assert law.is_stepped(s) and s.sweep_period_s == pytest.approx(4.2612)
    # Rev G's rounded 5.6 MHz/s is accepted.
    assert schedule_from_config({**cfg, "k_f_hz_per_s": 5.6e6}).k_f_hz_per_s == s.k_f_hz_per_s
    with pytest.raises(ValueError, match="step law"):
        schedule_from_config({**cfg, "k_f_hz_per_s": 1.0e6})


def test_uniform_schedule_is_unchanged():
    u = SweepSchedule(f_start_hz=1.0e6, f_stop_hz=5.0e6, k_f_hz_per_s=5.6e6,
                      sweep_period_s=1.0)
    assert not law.is_stepped(u)
    assert sweep_duration_s(u) == pytest.approx(4.0e6 / 5.6e6)


@pytest.mark.parametrize("range_m", [3.0e5, 6.0e5])
def test_stepped_echo_lands_at_its_range_across_the_band(range_m):
    """L2 to L3 on the A1 law: the range is right below and above 3 MHz."""
    from tss.l2_instrument.baseline import BaselineInstrument

    sweep = law.stepped_schedule(3.3e5, 6.0e6, 0.0268)
    freq = law.step_centres_hz(sweep)
    echoes = EchoSet(
        freq_hz=freq,
        group_delay_s=np.full(freq.size, 2.0 * range_m / C_LIGHT),
        # -60 dB rather than the -90 dB of the uniform test: the sky is some
        # 30 dB hotter at 0.33 MHz than at 3 MHz, and this test is about the
        # range scale, not the sensitivity at the band edge.
        amplitude_db=np.full(freq.size, -60.0),
        absorption_db=np.zeros(freq.size),
        mode=np.array(["O"] * freq.size),
    )
    signal = BaselineInstrument().receive(echoes, sweep, np.random.default_rng(3))
    # The WP3-era fixed rule (WP4 Phase 2): at -60 dB the echo is about 80 dB
    # over the noise at the top of this 6 MHz sweep, and in its short last
    # step, swept at a much lower rate, the energy that arrives during the
    # step transition smears from 400 to 1050 km.  The stated-P_fa detector
    # rightly detects it, and the range jump then fails the chain rule.  This
    # test is about the range scale, so it keeps the rule it was written for;
    # the partial last step is a step-transition item for WP4 Phase 3.
    tier = BaselineSignalProcessing(detector="fixed")

    # The range scale, checked on the spectrum peak of every column.  A block
    # converted at the peak rate would put the bottom-of-band echoes at a ninth
    # of their range, so this is the check that each block uses its own step.
    spectra = tier.spectra(signal)
    assert len(spectra) == law.step_count(sweep)
    for sp in spectra:
        peak = float(sp.range_m[int(np.argmax(sp.power_db[0]))])
        assert abs(peak - range_m) <= sp.range_resolution_m

    ionogram = tier.form_ionogram(signal)
    assert np.allclose(ionogram.freq_hz, freq)
    ranges = ionogram.virtual_range_m[:, 0]
    assert np.all(np.isfinite(ranges))
    # The leading edge is never beyond the echo.  How far before it sits
    # depends on SNR and the window skirt, which is not what this test is about.
    cells = np.asarray(ionogram.flags["range_resolution_m_by_column"])
    assert np.all(ranges - range_m <= cells)
    # The cell is coarse at the bottom and fine at the top, as the step law says.
    assert cells[0] > 5.0 * np.min(cells)
    assert ionogram.flags["range_resolution_m"] == pytest.approx(cells.max())
    assert ionogram.flags["sweep_law"] == "stepped"


def test_campaign_grid_uses_step_centres():
    from tss.studies import _grids

    cfg = {"sweep": {"f_start_hz": 3.3e5, "f_stop_hz": 2.0e7, "step_dwell_s": 0.0268,
                     "step_fraction": 0.05, "step_cap_hz": 1.5e5}}
    _, freq = _grids(cfg)
    assert freq.size == 159
    assert np.allclose(freq, law.step_centres_hz(A1))


class _NoNoise:
    """Stands in for the generator so that L2 adds no noise."""

    def normal(self, loc=0.0, scale=1.0, size=None):
        shape = np.shape(scale) if size is None else size
        return np.zeros(shape)


def _single_step_signal(step: int, range_m: float):
    """L2 output for one echo row filed in one A1 step, noise-free."""
    from tss.l2_instrument.baseline import BaselineInstrument

    sweep = law.stepped_schedule(3.3e5, 6.0e6, 0.0268)
    f = law.step_centres_hz(sweep)[[step]]
    echoes = EchoSet(
        freq_hz=f,
        group_delay_s=np.full(1, 2.0 * range_m / C_LIGHT),
        amplitude_db=np.full(1, -60.0),
        absorption_db=np.zeros(1),
        mode=np.array(["O"]),
    )
    signal = BaselineInstrument(quantise=False).receive(echoes, sweep, _NoNoise())
    return sweep, signal


def test_l2_carries_the_exact_step_transition():
    """In the first tau of a step the beat ramps from K_(i-1) tau to K_i tau."""
    range_m = 1.0e6
    tau = 2.0 * range_m / C_LIGHT
    step = 3
    sweep, signal = _single_step_signal(step, range_m)
    assert signal.flags["step_transition_modelled"] is True
    fs = float(signal.fs_hz)
    x = np.asarray(signal.samples)[0]
    t = np.arange(x.size) / fs
    t_i = step * sweep.step_dwell_s
    k = law.step_rates_hz_per_s(sweep)
    # Instantaneous beat from the phase increment, well away from the edges.
    beat = np.diff(np.unwrap(np.angle(x))) * fs / (2.0 * np.pi)
    mid = 0.5 * (t[1:] + t[:-1])
    steady = (mid > t_i + tau + 1e-3) & (mid < t_i + sweep.step_dwell_s - 1e-3)
    early = (mid > t_i + 0.1 * tau) & (mid < t_i + 0.2 * tau)
    assert np.allclose(beat[steady], k[step] * tau, rtol=1e-4)
    expected = k[step - 1] * tau + (k[step] - k[step - 1]) * (mid[early] - t_i)
    assert np.allclose(beat[early], expected, rtol=1e-4)


def test_step_transition_puts_a_skirt_short_of_the_echo():
    """Whole-dwell spectrum shows the skirt; gating the first tau removes it."""
    range_m = 1.0e6
    tau = 2.0 * range_m / C_LIGHT
    step = 2
    sweep, signal = _single_step_signal(step, range_m)
    fs = float(signal.fs_hz)
    x = np.asarray(signal.samples)[0]
    t = np.arange(x.size) / fs
    t_i = step * sweep.step_dwell_s
    k = float(law.step_rates_hz_per_s(sweep)[step])

    def peak_between(sel, lo_m, hi_m):
        s = x[sel]
        w = np.hanning(s.size)
        spec = np.abs(np.fft.fft(s * w, 16 * s.size)) / w.sum()
        fr = np.fft.fftfreq(spec.size, 1.0 / fs)
        r = C_LIGHT * fr / (2.0 * k)
        band = (r > lo_m) & (r < hi_m)
        return 20 * np.log10(spec[band].max() / spec.max())

    dwell = (t >= t_i) & (t < t_i + sweep.step_dwell_s)
    gated = (t >= t_i + tau) & (t < t_i + sweep.step_dwell_s)
    # 930 to 970 km is outside the main lobe of both windows (cells of 9 and
    # 12 km with Hann) and inside the R/1.05 to R span of the transition.
    whole = peak_between(dwell, 9.3e5, 9.7e5)
    clean = peak_between(gated, 9.3e5, 9.7e5)
    assert whole > -45.0
    assert whole > clean + 6.0


def test_no_echo_before_the_sweep_has_been_transmitted():
    sweep, signal = _single_step_signal(0, 1.0e6)
    fs = float(signal.fs_hz)
    tau = 2.0 * 1.0e6 / C_LIGHT
    x = np.asarray(signal.samples)[0]
    assert np.all(x[: int(0.9 * tau * fs)] == 0)
