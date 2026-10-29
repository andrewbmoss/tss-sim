"""WP3 Phase 2: transmitter impairments in the samples, against their predictions.

Exit criterion: each impairment switches on and off, and with INS-21 off the
range-correlation credit is lost by the predicted amount.  Densities are
measured from the de-chirped samples with a Welch estimate and compared with
the prediction averaged over the same bins, to 0.5 dB.
"""

from __future__ import annotations

import math

import numpy as np
import pytest
from scipy import signal as sps
from scipy import special

from tss import sweep_law
from tss.config import load_yaml, resolve_path
from tss.constants import C_LIGHT
from tss.interfaces import EchoSet, SweepSchedule
from tss.l2_instrument import impairments as im
from tss.l2_instrument import leakage as lk
from tss.l2_instrument.baseline import BaselineInstrument
from tss.l2_instrument.oscillator import Oscillator
from tss.studies.data_volume import A1

TOL_DB = 0.5


def _uniform(f_start=5.0e6, k=5.6e6, seconds=1.0, direction=1):
    return SweepSchedule(
        f_start_hz=f_start, f_stop_hz=f_start + k * seconds, k_f_hz_per_s=k,
        sweep_period_s=seconds, direction=direction,
    )


def _density_dbc_hz(x, fs, f_lo, f_hi, nperseg=4096):
    """Mean one-sided-offset density (positive offsets) relative to the carrier."""
    f, p = sps.welch(x, fs=fs, window="hann", nperseg=nperseg, return_onesided=False,
                     scaling="density", detrend=False)
    sel = (f >= f_lo) & (f <= f_hi)
    return f[sel], p[sel] / np.mean(np.abs(x) ** 2)


def _mean_db(lin):
    return 10.0 * math.log10(float(np.mean(lin)))


def _leak(imp, path_m, *, sweep=None, fs=100.0e3, n=1 << 16, seed=0, t0=0.01):
    m = lk.LeakageModel.from_coupling_path(path_m)
    return lk.leakage_signal(sweep or _uniform(), fs, m, np.random.default_rng(seed),
                             n_samples=n, t0_s=t0, impairments=imp)


# ----------------------------------------------------------------------------
# The specification
# ----------------------------------------------------------------------------


def test_no_impairment_is_the_phase1_leakage():
    sweep = A1.schedule()
    m = lk.LeakageModel.from_coupling_path(10.0)
    a = lk.leakage_signal(sweep, 78.4e3, m, np.random.default_rng(0), n_samples=2000, t0_s=1.0)
    b = lk.leakage_signal(sweep, 78.4e3, m, np.random.default_rng(0), n_samples=2000, t0_s=1.0,
                          impairments=im.TransmitImpairments.none())
    assert np.allclose(a, b, rtol=0, atol=1e-15)
    assert not im.TransmitImpairments.none().any


def test_the_power_law_reproduces_the_allan_deviation():
    for noise in ("white_fm", "flicker_fm", "random_walk_fm"):
        osc = Oscillator({1.0: 2.0e-9, 2.0: 2.0e-9}, noise_type=noise)
        h, a = im.power_law_h(osc)
        sigma2 = {0: h / 2.0, -1: 2.0 * math.log(2.0) * h, -2: 2.0 * math.pi**2 * h / 3.0}[a]
        assert math.sqrt(sigma2) == pytest.approx(2.0e-9, rel=1e-12)


def test_spur_indices():
    beta = im.pm_spur_beta(-40.0)
    assert special.j1(beta) / special.j0(beta) == pytest.approx(0.01, rel=1e-9)
    assert im.am_spur_index(-40.0) == pytest.approx(0.02)
    with pytest.raises(ValueError):
        im.Spur(offset_hz=1e3, level_dbc=3.0)


def test_the_configuration_block_reads_and_absent_means_none():
    cfg = load_yaml(resolve_path("instrument/baseline_6u.yaml"))
    inst = BaselineInstrument.from_config(cfg)
    # Since WP3 Phase 4 the scored transmitter carries its allocations.
    assert inst.impairments.any and inst.impairments.am_noise_dbc_hz == -140.0   # INS-30
    assert not im.TransmitImpairments.from_config(None).any
    assert not im.TransmitImpairments.from_config({}).any
    imp = im.TransmitImpairments.from_config(
        {"oscillator_noise": True, "am_noise_dbc_hz": -130, "spurs": [{"offset_hz": 1e4, "level_dbc": -70}],
         "doppler": True},
        oscillator=inst.oscillator,
    )
    assert imp.oscillator is inst.oscillator and imp.single_reference and imp.doppler
    assert imp.spurs[0].kind == "pm" and imp.am_noise_dbc_hz == -130.0


# ----------------------------------------------------------------------------
# Phase noise: range correlation and INS-21
# ----------------------------------------------------------------------------


@pytest.mark.parametrize("path_m", [3.0, 30.0])
def test_a_phase_noise_floor_leaves_L_times_4sin2(path_m):
    l_dbc = -90.0
    x = _leak(im.TransmitImpairments(phase_noise_floor_dbc_hz=l_dbc), path_m)
    tau = path_m / C_LIGHT
    for lo, hi in ((3e3, 8e3), (15e3, 25e3), (30e3, 45e3)):
        f, d = _density_dbc_hz(x, 100e3, lo, hi)
        pred = 10 ** (im.dechirped_pm_density_dbc_hz(l_dbc, f, tau) / 10)
        assert _mean_db(d) == pytest.approx(_mean_db(pred), abs=TOL_DB)


def test_an_independent_reference_loses_the_credit_by_the_predicted_amount():
    """INS-21 off: the residual is 2L, not L 4 sin^2(pi f tau)."""
    l_dbc, path = -90.0, 10.0
    tau = path / C_LIGHT
    common = _leak(im.TransmitImpairments(phase_noise_floor_dbc_hz=l_dbc), path)
    indep = _leak(im.TransmitImpairments(phase_noise_floor_dbc_hz=l_dbc, single_reference=False), path)
    f, dc = _density_dbc_hz(common, 100e3, 15e3, 25e3)
    _, di = _density_dbc_hz(indep, 100e3, 15e3, 25e3)
    lost = _mean_db(di) - _mean_db(dc)
    pc = 10 ** (im.dechirped_pm_density_dbc_hz(l_dbc, f, tau) / 10)
    pi_ = 10 ** (im.dechirped_pm_density_dbc_hz(l_dbc, f, tau, single_reference=False) / 10)
    assert lost == pytest.approx(_mean_db(pi_) - _mean_db(pc), abs=TOL_DB)
    assert lost > 40.0
    assert _mean_db(di) == pytest.approx(l_dbc + 10 * math.log10(2.0), abs=TOL_DB)


@pytest.mark.parametrize("single_reference", [True, False])
def test_oscillator_noise_at_an_echo_delay_follows_the_range_correlation(single_reference):
    """Flicker-FM reference noise multiplied to 10 MHz, through a 7 ms path (1050 km)."""
    osc = Oscillator({1.0: 2.0e-9, 2.0: 2.0e-9}, noise_type="flicker_fm")
    imp = im.TransmitImpairments(oscillator=osc, single_reference=single_reference)
    fc, tau, fs = 10.0e6, 7.0e-3, 20.0e3
    sweep = _uniform(f_start=fc, k=1.0, seconds=40.0)       # the carrier barely moves
    m = lk.LeakageModel(coupling_delay_s=tau, isolation_db=0.0)
    x = lk.leakage_signal(sweep, fs, m, np.random.default_rng(3), n_samples=1 << 19, t0_s=0.02,
                          impairments=imp)
    for lo, hi in ((60.0, 130.0), (300.0, 700.0), (1500.0, 3000.0)):
        f, d = _density_dbc_hz(x, fs, lo, hi, nperseg=1 << 14)
        lref = im.oscillator_phase_density_dbc_hz(osc, fc, f)
        pred_db = im.dechirped_pm_density_dbc_hz(0.0, f, tau, single_reference=single_reference) + lref
        assert _mean_db(d) == pytest.approx(_mean_db(10 ** (pred_db / 10)), abs=TOL_DB)


# ----------------------------------------------------------------------------
# AM noise and spurs
# ----------------------------------------------------------------------------


@pytest.mark.parametrize("path_m", [1.0, 30.0])
@pytest.mark.parametrize("single_reference", [True, False])
def test_am_noise_reaches_the_residual_whole(path_m, single_reference):
    """Committee Review 5: AM is not range-correlated, whatever the path or INS-21."""
    l_am = -110.0
    x = _leak(im.TransmitImpairments(am_noise_dbc_hz=l_am, single_reference=single_reference), path_m)
    _, d = _density_dbc_hz(x, 100e3, 3e3, 45e3)
    assert _mean_db(d) == pytest.approx(l_am, abs=TOL_DB)


def test_pm_spurs_are_suppressed_and_am_spurs_are_not():
    fs, n = 100.0e3, 5000                 # 20 Hz bins
    path = 10.0
    tau = path / C_LIGHT
    for kind, f_s in (("pm", 20.0e3), ("am", 12.0e3)):
        imp = im.TransmitImpairments(spurs=(im.Spur(offset_hz=f_s, level_dbc=-40.0, kind=kind),))
        x = _leak(imp, path, fs=fs, n=n)
        w = np.hanning(n)
        spec = np.abs(np.fft.fft(x * w))
        side = 20 * np.log10(spec[int(round(f_s * n / fs))] / (np.sum(w) * np.sqrt(np.mean(np.abs(x) ** 2))))
        if kind == "pm":
            beta_r = 2.0 * im.pm_spur_beta(-40.0) * math.sin(math.pi * f_s * tau)
            expected = 20 * math.log10(special.j1(beta_r) / special.j0(beta_r))
        else:
            expected = -40.0
        assert side == pytest.approx(expected, abs=TOL_DB)


# ----------------------------------------------------------------------------
# Sweep nonlinearity (INS-25)
# ----------------------------------------------------------------------------


def _nl_local_cycles(sweep, t):
    dwell = float(sweep.step_dwell_s)
    i = sweep_law.step_index(sweep, t)
    u = t - i * dwell
    s = u / dwell
    width = sweep_law.rate_hz_per_s(sweep, t) * dwell
    return width * float(sweep.nonlinearity_ppm) * 1e-6 * 4 * dwell * (s**2 / 2 - s**3 / 3)


def test_sweep_nonlinearity_is_range_correlated_only_with_a_common_reference():
    eps_ppm = 2000.0
    base = A1.schedule()
    sweep = sweep_law.stepped_schedule(base.f_start_hz, base.f_stop_hz, base.step_dwell_s,
                                       nonlinearity_ppm=eps_ppm)
    dwell = float(sweep.step_dwell_s)
    fs = 200.0e3
    i = int(sweep_law.step_index(sweep, sweep_law.time_at_frequency_s(sweep, 5.0e6)))
    t0, n = i * dwell + 1e-4, int(0.9 * dwell * fs)
    t = t0 + np.arange(n) / fs
    path = 30.0
    tau = path / C_LIGHT
    m = lk.LeakageModel.from_coupling_path(path)
    ideal = lk.leakage_signal(sweep, fs, m, np.random.default_rng(0), n_samples=n, t0_s=t0,
                              impairments=im.TransmitImpairments())
    common = lk.leakage_signal(sweep, fs, m, np.random.default_rng(0), n_samples=n, t0_s=t0,
                               impairments=im.TransmitImpairments(sweep_nonlinearity=True))
    indep = lk.leakage_signal(sweep, fs, m, np.random.default_rng(0), n_samples=n, t0_s=t0,
                              impairments=im.TransmitImpairments(sweep_nonlinearity=True,
                                                                 single_reference=False))
    c_now, c_then = _nl_local_cycles(sweep, t), _nl_local_cycles(sweep, t - tau)
    got_common = np.angle(common * np.conj(ideal))
    assert np.allclose(got_common, 2 * np.pi * (c_now - c_then), atol=1e-9)
    # With an independent ideal reference the whole in-step error survives (less
    # the constant of the step's accumulated phase, removed by referring to t0).
    got_indep = np.unwrap(np.angle(indep * np.conj(ideal)))
    want = -2 * np.pi * c_then
    assert np.allclose(got_indep - got_indep[0], want - want[0], atol=1e-6)
    # The range correlation credit on this error is the ratio of the two.
    credit_db = 20 * np.log10(np.ptp(want) / np.ptp(2 * np.pi * (c_now - c_then)))
    assert credit_db > 60.0


# ----------------------------------------------------------------------------
# Realizations and the baseline's echoes
# ----------------------------------------------------------------------------


def test_one_realization_serves_the_leakage_and_must_match_its_record():
    imp = im.TransmitImpairments(phase_noise_floor_dbc_hz=-100.0)
    sweep = _uniform()
    real = imp.realize(sweep, 100e3, 4096, np.random.default_rng(1), t0_s=0.01, max_delay_s=1e-6)
    m = lk.LeakageModel.from_coupling_path(10.0)
    a = lk.leakage_signal(sweep, 100e3, m, np.random.default_rng(9), n_samples=4096, t0_s=0.01,
                          impairments=imp, realization=real)
    b = lk.leakage_signal(sweep, 100e3, m, np.random.default_rng(5), n_samples=4096, t0_s=0.01,
                          impairments=imp, realization=real)
    assert np.array_equal(a, b)
    with pytest.raises(ValueError):
        lk.leakage_signal(sweep, 100e3, m, np.random.default_rng(5), n_samples=4000, t0_s=0.01,
                          impairments=imp, realization=real)


def _echo(f_hz, range_m, amp_db=-80.0, doppler_hz=None):
    return EchoSet(
        freq_hz=np.array([f_hz]), mode=np.array(["O"], dtype=object),
        group_delay_s=np.array([2.0 * range_m / C_LIGHT]), amplitude_db=np.array([amp_db]),
        absorption_db=np.array([0.0]), flags=np.array([""], dtype=object), source="test",
        doppler_hz=None if doppler_hz is None else np.array([doppler_hz]),
    )


@pytest.mark.parametrize("direction", [1, -1])
def test_doppler_moves_the_echo_beat_by_minus_fd(direction):
    cfg = load_yaml(resolve_path("instrument/baseline_6u.yaml"))
    sweep = sweep_law.stepped_schedule(0.33e6, 20e6, 0.0268, direction=direction)
    fd = 37.0
    echo = _echo(5.0e6, 600e3, doppler_hz=fd)
    out = {}
    for on in (False, True):
        inst = BaselineInstrument.from_config(
            cfg, impairments=im.TransmitImpairments(doppler=on), quantise=False)
        out[on] = inst.receive(echo, sweep, np.random.default_rng(4))
    s_on, s_off = out[True].samples[0], out[False].samples[0]
    fs = out[True].fs_hz
    i = int(sweep_law.step_index(sweep, sweep_law.time_at_frequency_s(sweep, 5.0e6)))
    lo = int((i * 0.0268 + 0.005) * fs)
    hi = int(((i + 1) * 0.0268 - 0.001) * fs)
    ph = np.unwrap(np.angle(s_on[lo:hi] * np.conj(s_off[lo:hi])))
    slope_hz = np.polyfit(np.arange(lo, hi) / fs, ph, 1)[0] / (2 * np.pi)
    assert slope_hz == pytest.approx(-fd, abs=0.5)
    assert out[True].flags["rows_doppler_applied"] == 1
    assert out[False].flags["rows_doppler_applied"] == 0


def test_the_baseline_carries_every_impairment_on_its_echoes():
    cfg = load_yaml(resolve_path("instrument/baseline_6u.yaml"))
    base = BaselineInstrument.from_config(
        cfg, quantise=False, impairments=im.TransmitImpairments.none(), leakage_enabled=False
    )
    imp = im.TransmitImpairments(
        oscillator=base.oscillator, phase_noise_floor_dbc_hz=-120.0, am_noise_dbc_hz=-120.0,
        spurs=(im.Spur(1.0e4, -60.0), im.Spur(2.0e4, -60.0, kind="am")), sweep_nonlinearity=True,
        doppler=True,
    )
    inst = BaselineInstrument.from_config(cfg, impairments=imp, quantise=False, leakage_enabled=False)
    sweep = A1.schedule()
    echo = _echo(5.0e6, 600e3, doppler_hz=10.0)
    clean = base.receive(echo, sweep, np.random.default_rng(2))
    dirty = inst.receive(echo, sweep, np.random.default_rng(2))
    assert dirty.flags["transmit_impairments"]["oscillator"] is not None
    assert clean.flags["transmit_impairments"]["oscillator"] is None
    assert not np.allclose(clean.samples, dirty.samples)
