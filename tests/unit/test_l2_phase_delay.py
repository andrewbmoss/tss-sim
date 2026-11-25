"""WP5 Phase 2: the phase path in L2 (EchoSet.phase_delay_s, D8 condition 1).

An echo's carrier phase follows its phase path, ``2 pi f tau_p``, while the
beat and the envelope follow the group delay.  To first order in the
sub-band (stationary phase, ``d(f tau_p)/df = tau_g``) the de-chirped phase is
the group-delay form plus ``2 pi f_k (tau_p - tau_g)``.  An EchoSet without the
field, or with it equal to the group delay, gives the samples it always gave.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from tss import sweep_law
from tss.interfaces import EchoSet
from tss.l2_instrument.baseline import BaselineInstrument, SequenceContext

FC = 5.0e6
TAU_G = 2 * 500e3 / 299_792_458.0


def _sweep():
    return sweep_law.drift_mode(FC)


def _echo(tau_g=TAU_G, tau_p=None, doppler=None, amp_db=-120.0):
    return EchoSet(
        freq_hz=np.array([FC]), mode=np.array(["O"], dtype=object),
        group_delay_s=np.array([tau_g]), amplitude_db=np.array([amp_db]),
        absorption_db=np.array([0.0]),
        doppler_hz=None if doppler is None else np.array([doppler]),
        flags=np.array(["facet:f1"], dtype=object),
        phase_delay_s=None if tau_p is None else np.array([tau_p]),
    )


def _l2(**kw):
    return BaselineInstrument(quantise=False, thermal_noise=False, **kw)


def _phase_at_centre(sig, tau_g, f_d=0.0):
    """Demodulate the facet tone and refer its phase to the sweep centre."""
    sw = sig.sweep
    k = sw.k_f_hz_per_s
    fs = sig.fs_hz
    x = sig.samples[0]
    t = np.arange(x.size) / fs
    live = np.abs(x) > 0
    tc = 0.5 * x.size / fs
    z = np.sum(x[live] * np.exp(-2j * np.pi * (k * tau_g - f_d) * (t[live] - tc)))
    return float(np.angle(z))


def test_drift_mode_sweep_law():
    sw = sweep_law.drift_mode(FC)
    assert sw.f_start_hz == pytest.approx(FC - 75e3) and sw.f_stop_hz == pytest.approx(FC + 75e3)
    assert sw.sweep_period_s == pytest.approx(1 / 45)
    assert sw.k_f_hz_per_s == pytest.approx(150e3 * 45)
    assert sweep_law.drift_mode(FC, repetition_hz=71.0).sweep_period_s == pytest.approx(1 / 71)
    t = sweep_law.drift_mode_times(5)
    assert np.allclose(np.diff(t), 1 / 45) and t[0] == 0.0


def test_absent_or_equal_phase_delay_changes_nothing():
    l2 = BaselineInstrument()
    a = l2.receive(_echo(), _sweep(), np.random.default_rng(3))
    b = l2.receive(_echo(tau_p=TAU_G), _sweep(), np.random.default_rng(3))
    assert np.array_equal(a.samples, b.samples)


def test_carrier_phase_follows_the_phase_path():
    tau_p = TAU_G * 0.93  # phase delay shorter than group delay in a plasma
    sig = _l2().receive(_echo(tau_p=tau_p), _sweep(), np.random.default_rng(0))
    assert sig.flags["n_clipped"] == 0
    sw = _sweep()
    k = sw.k_f_hz_per_s
    tc = 0.5 * sig.samples.shape[1] / sig.fs_hz  # the record's centre
    f0 = sw.f_start_hz
    want = 2 * np.pi * (f0 * TAU_G + k * TAU_G * tc - k * TAU_G**2 / 2 + FC * (tau_p - TAU_G))
    got = _phase_at_centre(sig, TAU_G)
    assert abs(np.angle(np.exp(1j * (got - want)))) < 1e-6
    assert sig.flags["rows_phase_path"] == 1


def test_closing_echo_phase_and_doppler_agree_in_sign():
    # A closing echo: the phase path shortens by v T each sweep, f_d = -f dtau/dt > 0.
    v = 25.0
    T = _sweep().sweep_period_s
    dtau = 2 * v * T / 299_792_458.0
    f_d = FC * 2 * v / 299_792_458.0
    l2 = _l2()
    from tss.l2_instrument.impairments import TransmitImpairments

    l2.impairments = TransmitImpairments(doppler=True)
    s1 = l2.receive(_echo(tau_p=TAU_G, doppler=f_d), _sweep(), np.random.default_rng(0))
    s2 = l2.receive(_echo(tau_p=TAU_G - dtau, doppler=f_d), _sweep(), np.random.default_rng(0))
    p1, p2 = _phase_at_centre(s1, TAU_G, f_d), _phase_at_centre(s2, TAU_G, f_d)
    step = np.angle(np.exp(1j * (p2 - p1)))
    assert step < 0.0
    assert step == pytest.approx(-2 * np.pi * f_d * T, rel=1e-6)
    # Inside one sweep the Doppler turns the phase the same way, by -2 pi f_d T/2
    # from the start to the centre; the phase is continuous across the boundary.
    assert -2 * np.pi * f_d * T / 2 < 0.0


def test_blanking_the_first_delay_of_a_repeated_sweep():
    sig = _l2().receive(_echo(tau_p=TAU_G), _sweep(), np.random.default_rng(0),
                        sequence=SequenceContext(blank_first_delay=True))
    n_tau = int(math.floor(TAU_G * sig.fs_hz))
    assert np.all(sig.samples[0, :n_tau] == 0)
    assert np.all(np.abs(sig.samples[0, n_tau + 2:]) > 0)
    assert sig.flags["first_delay_blanked"] is True
    plain = _l2().receive(_echo(tau_p=TAU_G), _sweep(), np.random.default_rng(0))
    assert np.all(np.abs(plain.samples[0]) > 0)
    assert plain.flags["first_delay_blanked"] is False


def test_thermal_noise_off_is_deterministic():
    a = _l2().receive(_echo(), _sweep(), np.random.default_rng(1))
    b = _l2().receive(_echo(), _sweep(), np.random.default_rng(2))
    assert np.array_equal(a.samples, b.samples)
    assert a.flags["thermal_noise"] is False


def test_facet_rows_never_join_a_dispersion_trace():
    e = EchoSet(
        freq_hz=np.array([FC - 10e3, FC + 10e3]), mode=np.array(["O", "O"], dtype=object),
        group_delay_s=np.array([TAU_G, TAU_G * 1.001]), amplitude_db=np.array([-120.0, -120.0]),
        absorption_db=np.zeros(2), flags=np.array(["facet:a", "facet:b"], dtype=object),
        phase_delay_s=np.array([TAU_G, TAU_G]),
    )
    sig = _l2().receive(e, _sweep(), np.random.default_rng(0))
    assert sig.flags["rows_dispersed"] == 0
