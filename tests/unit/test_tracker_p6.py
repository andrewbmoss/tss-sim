"""WP4 Phase 6: the f_N trace floor, stop-at-loss and the onboard starting value."""

from __future__ import annotations

from types import SimpleNamespace

import numpy as np
import pytest

from tss.l3_signal import tracker as tk
from tss.l4_retrieval.base import starting_plasma_frequency_hz


def _image(freq, ranges, n_r=500, grid=2e3, snr=25.0):
    img = np.full((freq.size, n_r), 0.0)
    for i, r in ranges.items():
        img[i, int(r / grid)] = snr
    return img


def test_trace_floor_takes_the_far_trace_below_a_near_start():
    # O starts far out at column 2; X starts at near range at column 6 (the P5 failure)
    f = 1.0e6 * 1.05 ** np.arange(20)
    ranges = {i: 400e3 + 10e3 * i for i in range(2, 20)}
    ranges.update({i: 20e3 + 30e3 * (i - 6) for i in range(6, 12)})
    img = _image(f, {**{i: ranges[i] for i in ranges}})
    for i in range(6, 12):
        img[i, int((20e3 + 30e3 * (i - 6)) / 2e3)] = 25.0
    streak = np.zeros(f.size, dtype=bool)
    old = tk.onboard_priors(f, img, streak, 1.0e6, 2e3, tk.TrackerSettings(fn_prior="near_start"))
    new = tk.onboard_priors(f, img, streak, 1.0e6, 2e3, tk.TrackerSettings(fn_prior="trace_floor"))
    assert old["f_n_method"] == "trace_start" and old["f_n_hz"] > f[5]
    assert new["f_n_method"] == "trace_floor"
    assert f[1] < new["f_n_hz"] < f[2]


def test_trace_floor_keeps_a_good_near_start():
    f = 1.0e6 * 1.05 ** np.arange(20)
    img = _image(f, {i: 10e3 + 20e3 * i for i in range(3, 20)})
    streak = np.zeros(f.size, dtype=bool)
    new = tk.onboard_priors(f, img, streak, 0.8e6, 2e3, tk.TrackerSettings(fn_prior="trace_floor"))
    assert new["f_n_method"] == "trace_start"


def test_unknown_fn_prior_refused():
    with pytest.raises(ValueError):
        tk.TrackerSettings(fn_prior="guess")


def test_stop_at_loss_on_an_x_takeover():
    f = np.arange(10.0)
    r = np.linspace(300e3, 700e3, 10)
    q = np.zeros(10, dtype=int)
    conf = np.ones(10)
    conf[7:] = 0.1
    assert tk.stop_at_loss(f, r, q, conf, tk.TrackerSettings(stop_run=2)) == 7.0
    conf[:] = 1.0
    conf[9] = 0.1        # one low sample at the very end is a run of the rest of the trace
    assert tk.stop_at_loss(f, r, q, conf, tk.TrackerSettings(stop_run=2)) == 9.0
    conf[:] = 1.0
    assert np.isnan(tk.stop_at_loss(f, r, q, conf, tk.TrackerSettings()))


def test_stop_at_loss_drops_a_track_only_tail():
    f = np.arange(10.0)
    r = np.linspace(300e3, 700e3, 10)
    q = np.zeros(10, dtype=int)
    q[8:] = tk.Q_TRACK_ONLY
    assert tk.stop_at_loss(f, r, q, np.ones(10), tk.TrackerSettings()) == 8.0


def test_onboard_starting_value():
    ion = SimpleNamespace(freq_hz=np.array([1.0e6, 1.1e6, 1.2e6]),
                          virtual_range_m=np.array([[100e3], [200e3], [300e3]]),
                          snr_db=np.array([[20.0], [20.0], [20.0]]), modes=("O",), flags={})
    f0, src = starting_plasma_frequency_hz(ion, onboard_hz=0.97e6)
    assert (f0, src) == (0.97e6, "onboard")
    # an estimate at or above the lowest detected frequency is not a plasma frequency at the spacecraft
    f0, src = starting_plasma_frequency_hz(ion, onboard_hz=1.05e6)
    assert (f0, src) == (1.0e6, "lowest_trace")


def test_phase6_defaults():
    s = tk.TrackerSettings()
    assert (s.fn_prior, s.far_start_cost, s.stop_at_loss, s.cusp_up_cap_m) == ("trace_floor", 2.0, False, s.up_cap_m)
    from tss.l4_retrieval.lamination import MagnetoionicLaminationRetrieval
    assert MagnetoionicLaminationRetrieval().start_policy == "onboard"


def test_resonance_below_the_floor_and_a_harmonic_ignored():
    # O trace from column 20 (f about 2.65 MHz); f_N resonance at column 19;
    # a harmonic of f_H rings one column off its nominal frequency at column 5
    f = 1.0e6 * 1.05 ** np.arange(40)
    img = _image(f, {i: 300e3 + 10e3 * i for i in range(20, 40)})
    streak = np.zeros(f.size, dtype=bool)
    streak[[5, 19]] = True
    f_h = f[5] * 1.04            # 1 f_H, more than half a step from column 5
    pri = tk.onboard_priors(f, img, streak, f_h, 2e3, tk.TrackerSettings())
    assert pri["f_n_method"] == "resonance" and pri["f_n_hz"] == f[19]
    old = tk.onboard_priors(f, img, streak, f_h, 2e3, tk.TrackerSettings(fn_prior="near_start"))
    assert old["f_n_hz"] == f[5]      # P4's rule took the lowest streak
