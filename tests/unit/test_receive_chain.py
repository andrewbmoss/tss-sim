"""WP3 Phase 3: the receive chain with leakage (front end and leakage high-pass)."""

from __future__ import annotations

import copy
import math

import numpy as np
import pytest

from tss import sweep_law
from tss.config import load_yaml, resolve_path
from tss.constants import C_LIGHT
from tss.interfaces import EchoSet, SweepSchedule
from tss.l2_instrument import leakage as lk
from tss.l2_instrument import receive_chain as rc
from tss.l2_instrument.baseline import BaselineInstrument
from tss.studies.data_volume import A1


def _cfg(leakage=True):
    cfg = copy.deepcopy(load_yaml(resolve_path("instrument/baseline_6u.yaml")))
    cfg["leakage"]["enabled"] = leakage
    return cfg


def _uniform(f_start=5.0e6, k=5.6e6, seconds=2.0):
    return SweepSchedule(f_start_hz=f_start, f_stop_hz=f_start + k * seconds, k_f_hz_per_s=k,
                         sweep_period_s=seconds, direction=1)


# ----------------------------------------------------------------------------
# The high-pass
# ----------------------------------------------------------------------------


@pytest.mark.parametrize("order", [3, 4])
def test_the_discrete_high_pass_is_the_analog_butterworth(order):
    fs, fc = 78.4e3, 1000.0
    hp = rc.LeakageHighpass(order=order, fixed_corner_hz=fc)
    n = int(fs * 2)
    t = np.arange(n) / fs
    for f in (300.0, 700.0, 1000.0, 2000.0, 5000.0):
        x = np.exp(2j * np.pi * f * t)
        y = hp.apply(x, _uniform(), fs)
        got = 20 * np.log10(abs(np.mean(y[n // 2:] * np.conj(x[n // 2:]))))
        assert got == pytest.approx(float(hp.response_db(f, fc)[0]), abs=0.05)
    # The analog response is Butterworth: -3 dB at the corner, 20N dB/decade below.
    assert float(hp.response_db(fc, fc)[0]) == pytest.approx(-3.01, abs=0.01)
    assert float(hp.response_db(fc / 10, fc)[0]) == pytest.approx(-20 * order, abs=0.1)


def test_the_corner_follows_the_step_and_keeps_the_100km_echo():
    hp = rc.LeakageHighpass()
    sweep = A1.schedule()
    rates = sweep_law.step_rates_hz_per_s(sweep)
    corners = hp.corner_hz(rates)
    assert corners.min() == pytest.approx(0.5 * 2 * rates.min() * 1e5 / C_LIGHT)
    assert corners.min() == pytest.approx(205.5, abs=1.0)
    assert corners.max() == pytest.approx(1867, abs=5)
    loss = [-float(hp.response_db(2 * k * 1e5 / C_LIGHT, c)[0]) for k, c in zip(rates, corners, strict=True)]
    assert max(loss) == pytest.approx(0.017, abs=0.001)


def test_switching_the_corner_does_not_shake_the_leakage_loose():
    """Integrator states carry across a change of corner, as capacitor voltages
    would, so a near-DC leakage that the filter has already absorbed stays absorbed."""
    sweep = A1.schedule()
    fs = 78.4e3
    hp = rc.LeakageHighpass()
    t0 = 0.0
    n = int(40 * 0.0268 * fs)                       # 40 steps of 5 % corner changes
    t = t0 + np.arange(n) / fs
    x = np.exp(2j * np.pi * 0.05 * t)               # a carrier at 0.05 Hz, unit power
    y = hp.apply(x, sweep, fs, t0_s=t0)
    p = np.abs(y) ** 2
    boundaries = (np.arange(5, 40) * 0.0268 * fs).astype(int)
    worst = max(p[b:b + 200].max() for b in boundaries)
    assert 10 * np.log10(worst) < -100.0


def test_the_filter_state_round_trips():
    sweep = A1.schedule()
    fs = 78.4e3
    hp = rc.LeakageHighpass()
    x = np.exp(2j * np.pi * 0.3 * np.arange(20000) / fs) + 1e-3 * np.exp(2j * np.pi * 3e3 * np.arange(20000) / fs)
    whole = hp.apply(x, sweep, fs)
    a, st = hp.apply(x[:7000], sweep, fs, return_state=True)
    b = hp.apply(x[7000:], sweep, fs, t0_s=7000 / fs, state=st)
    # The trapezoid at the cut sees the last input repeated instead of the next
    # one, a perturbation of order Q (u[n+1] - u[n]); otherwise the two agree.
    assert np.allclose(np.concatenate([a, b]), whole, atol=1e-4)


# ----------------------------------------------------------------------------
# The front end
# ----------------------------------------------------------------------------


def test_the_rapp_front_end_is_1db_down_at_its_p1db():
    fe = rc.FrontEnd(input_p1db_dbm=1.0, am_pm_deg_at_p1db=3.0)
    r1 = math.sqrt(fe.input_p1db_w)
    g = fe.gain(np.array([r1 * 1e-3, r1]))
    assert 20 * np.log10(abs(g[0])) == pytest.approx(0.0, abs=1e-4)
    assert 20 * np.log10(abs(g[1])) == pytest.approx(-1.0, abs=1e-9)
    assert math.degrees(np.angle(g[1])) == pytest.approx(3.0)
    assert rc.FrontEnd().linear and np.array_equal(rc.FrontEnd().apply(np.array([1 + 1j])), np.array([1 + 1j]))


@pytest.mark.parametrize("backoff_db", [3.0, 10.0])
def test_a_weak_echo_beside_the_leakage_loses_gain_and_gains_an_image_as_predicted(backoff_db):
    fe = rc.FrontEnd(input_p1db_dbm=1.0, am_pm_deg_at_p1db=2.0)
    p_c = fe.input_p1db_w * 10 ** (-backoff_db / 10)
    fs, n, f_b = 10e3, 10000, 1234.0
    t = np.arange(n) / fs
    eps = 1e-4 * math.sqrt(p_c)
    x = math.sqrt(p_c) + eps * np.exp(2j * np.pi * f_b * t)
    y = fe.apply(x)
    spec = np.fft.fft(y) / n
    k = int(round(f_b * n / fs))
    gain = 20 * np.log10(abs(spec[k]) / eps)
    image = 20 * np.log10(abs(spec[-k]) / abs(spec[k]))
    assert gain == pytest.approx(rc.small_signal_gain_db(fe, p_c), abs=0.01)
    assert image == pytest.approx(rc.compression_image_dbc(fe, p_c), abs=0.05)


def test_the_bandpass_equivalent_model_is_the_passband_nonlinearity_after_de_chirp():
    """The plan's check: compression applied to the de-chirped envelope equals a real
    memoryless nonlinearity at radio frequency followed by the de-chirp, because
    the de-chirp is a frequency translation.  tanh at RF against its first-zone
    describing function at baseband."""
    v = 1.0
    th = np.linspace(0.0, 2 * np.pi, 4097)[:-1]

    def describing(r):
        r = np.atleast_1d(r)
        return np.array([2.0 * np.mean(v * np.tanh(ri * np.cos(th) / v) * np.cos(th)) for ri in r])

    rg = np.linspace(0.0, 4.0, 2001)
    dg = describing(rg)
    fe = rc.FrontEnd(am_am=lambda r: np.interp(r, rg, dg))

    fs_rf, dur = 8.0e6, 0.02
    t = np.arange(int(fs_rf * dur)) / fs_rf
    f0, k = 500e3, 5.6e6
    tau = 2.0e-3

    def phi(tt):
        return 2 * np.pi * (f0 * tt + k * tt**2 / 2)

    a_l, a_e = 1.0, 0.03
    s = a_l * np.cos(phi(t)) + a_e * np.cos(phi(t - tau))
    y = v * np.tanh(s / v)
    z = 2.0 * y * np.exp(1j * phi(t))                   # ref - rx convention, as the L2 tier
    zf = np.fft.fft(z)
    ff = np.fft.fftfreq(z.size, 1 / fs_rf)
    zf[np.abs(ff) > 200e3] = 0.0
    bb_rf = np.fft.ifft(zf)

    x_bb = a_l + a_e * np.exp(1j * (phi(t) - phi(t - tau)))
    bb_model = fe.apply(x_bb)

    keep = slice(int(0.003 * fs_rf), int(0.019 * fs_rf), 40)  # skip the edges, decimate
    beat = k * tau
    tt = t[keep]

    w = np.hanning(tt.size)

    def tone(b, f):
        return abs(np.sum(w * b[keep] * np.exp(-2j * np.pi * f * tt)))

    for f in (0.0, beat, -beat):
        got_rf, got_bb = tone(bb_rf, f), tone(bb_model, f)
        assert 20 * np.log10(got_rf / got_bb) == pytest.approx(0.0, abs=0.1)
    # The compression is real: the image is there, about 20 dB below the echo.
    assert -40 < 20 * np.log10(tone(bb_rf, -beat) / tone(bb_rf, beat)) < -5


# ----------------------------------------------------------------------------
# The baseline with leakage on
# ----------------------------------------------------------------------------


def test_leakage_is_off_by_default_and_on_loads_the_mom_bus_profile():
    off = BaselineInstrument.from_config(_cfg(False))
    assert not off.leakage_enabled and not off.highpass_enabled and not off.front_end_enabled
    on = BaselineInstrument.from_config(_cfg(True))
    assert on.leakage_enabled and on.highpass_enabled and on.front_end_enabled
    bus = lk.mom_anchor_table("bus")
    assert np.allclose(lk.isolation_profile_db(bus["f_hz"], on.leakage), bus["iso_crossed_ends_offset_db"])
    assert on.leakage.coupling_delay_s == 3.0e-9
    cfg = _cfg(True)
    cfg["leakage"]["profile"] = "flat"
    flat = BaselineInstrument.from_config(cfg)
    assert flat.leakage.table_f_hz is None and flat.leakage.isolation_db == 40.0


def _echo(f_hz, range_m, amp_db):
    return EchoSet(
        freq_hz=np.array([f_hz]), mode=np.array(["O"], dtype=object),
        group_delay_s=np.array([2.0 * range_m / C_LIGHT]), amplitude_db=np.array([amp_db]),
        absorption_db=np.array([0.0]), flags=np.array([""], dtype=object), source="test",
    )


def test_without_the_high_pass_the_digitiser_drowns_and_with_it_the_echo_survives():
    from tss.registry import build

    sweep = A1.schedule()
    echo = _echo(5.0e6, 600e3, -95.0)
    l3 = build("l3", "baseline")
    off = BaselineInstrument.from_config(_cfg(False))
    on = BaselineInstrument.from_config(_cfg(True))
    no_hp = BaselineInstrument.from_config(_cfg(True), highpass_enabled=False)
    r_off = off.receive(echo, sweep, np.random.default_rng(3))
    r_on = on.receive(echo, sweep, np.random.default_rng(3))
    r_no = no_hp.receive(echo, sweep, np.random.default_rng(3))
    assert r_on.flags["leakage_modelled"] and not r_off.flags["leakage_modelled"]
    assert r_no.flags["n_clipped"] > 0.9 * r_no.flags["n_samples"]       # INS-11 needs the filter
    snr_off = float(np.asarray(l3.snr_at(r_off, [5.0e6], [600e3]))[0])
    snr_on = float(np.asarray(l3.snr_at(r_on, [5.0e6], [600e3]))[0])
    assert snr_on == pytest.approx(snr_off, abs=1.0)


def test_the_study_scores_ins11_on_the_noise_rule():
    from tss.studies import leakage_receive_chain as st

    inst = st._instrument()
    r = st.per_step(inst, A1.schedule())
    assert r["f_hz"].size == sweep_law.step_count(A1.schedule())
    # Rev G's rule asks more than the noise rule at the bottom of the band.
    assert r["required_revg_rule_db"][0] > r["required_noise_rule_db"][0]
    # The filter's rejection of the carrier in the settled part of every step
    # exceeds the noise rule.
    assert np.all(r["achieved_db"] >= r["required_noise_rule_db"])
    assert np.max(r["echo_loss_100km_db"]) < 0.02
