"""WP5 Phase 2: what a sequence adds to L2 (oscillator continuity, specular
fading) and the layer echo beside the facet."""

from __future__ import annotations

import numpy as np
import pytest

from tss import sweep_law
from tss.interfaces import IonosphereState, PlatformState
from tss.l1_propagation import plate_echo as pe
from tss.l2_instrument.baseline import BaselineInstrument
from tss.l2_instrument.fading import FadingModel
from tss.l2_instrument.impairments import TransmitImpairments
from tss.l2_instrument.oscillator import Oscillator
from tss.registry import build

F = 5.0e6
PROF = pe.ExponentialTopside()
DEPTH = PROF.depth_for_mu(F, 0.5)


def _echo_sets(n, t0=-0.2):
    tier = build("l1", "facet", facets=[{"id": "f1", "depth_m": DEPTH, "drift_mps": 25.0}])
    t = sweep_law.drift_mode_times(n, t0_s=t0)
    st = IonosphereState(alt_m=np.linspace(1e5, 6e5, 11), ne=np.full(11, 1e10))
    sets = [tier.echoes(st, PlatformState(alt_m=6e5, lat_deg=0, lon_deg=0, t_s=x), np.array([F]),
                        np.random.default_rng(0)) for x in t]
    return sets, t


def _osc_l2(**kw):
    osc = Oscillator()
    imp = TransmitImpairments(oscillator=osc, doppler=True)
    return BaselineInstrument(quantise=False, thermal_noise=False, impairments=imp, **kw)


def test_one_oscillator_draw_spans_the_sequence():
    sets, t = _echo_sets(6)
    l2 = _osc_l2()
    sweep = sweep_law.drift_mode(F)
    real, offsets = l2.sequence_realization(sets, sweep, t, np.random.default_rng(4))
    recs = l2.receive_sequence(sets, sweep, t, np.random.default_rng(4))
    assert all(r.flags["oscillator_continuous"] for r in recs)
    n = recs[0].samples.shape[1]
    ref = real.ref_phase()
    # Consecutive records sit next to each other on one draw.
    assert offsets[1] - offsets[0] == pytest.approx(round(sweep.sweep_period_s * recs[0].fs_hz), abs=1)
    step_in = np.median(np.abs(np.diff(ref[offsets[0]: offsets[0] + n])))
    step_across = abs(ref[offsets[1]] - ref[offsets[0] + n - 1])
    assert step_across < 20 * step_in + 1e-12
    # The carrier is the sawtooth, not one ramp extrapolated.
    fc = real.carrier_hz_at(np.array([0.0, sweep.sweep_period_s + 1e-6]))
    assert fc[1] == pytest.approx(sweep.f_start_hz, abs=50.0)


def test_without_impairments_no_realization_is_drawn():
    sets, t = _echo_sets(3)
    l2 = BaselineInstrument(quantise=False, thermal_noise=False)
    recs = l2.receive_sequence(sets, sweep_law.drift_mode(F), t, np.random.default_rng(0))
    assert not any(r.flags["oscillator_continuous"] for r in recs)


def _tone(rec, echo, sweep):
    x = rec.samples[0]
    fs = rec.fs_hz
    tt = np.arange(x.size) / fs
    tau, fd = float(echo.group_delay_s[0]), float(echo.doppler_hz[0])
    live = np.abs(x) > 0
    return np.sum(x[live] * np.exp(-2j * np.pi * (sweep.k_f_hz_per_s * tau - fd) * tt[live]))


def test_specular_fading_turns_the_phase_and_keeps_the_amplitude():
    sets, t = _echo_sets(40)
    sweep = sweep_law.drift_mode(F)
    plain = BaselineInstrument(quantise=False, thermal_noise=False)
    faded = BaselineInstrument(quantise=False, thermal_noise=False,
                               fading=FadingModel(tau_i_s=0.5, tau_log_sigma=0.0))
    a = plain.receive_sequence(sets, sweep, t, np.random.default_rng(1))
    b = faded.receive_sequence(sets, sweep, t, np.random.default_rng(1))
    za = np.array([_tone(r, e, sweep) for r, e in zip(a, sets, strict=True)])
    zb = np.array([_tone(r, e, sweep) for r, e in zip(b, sets, strict=True)])
    assert np.allclose(np.abs(zb), np.abs(za), rtol=1e-9)
    dphi = np.angle(zb / za)
    assert np.std(dphi) > 0.01
    f = b[0].flags["fading"]["f1"]
    assert f["tau_i_s"] == 0.5 and f["specular"] == "phase"


def test_fading_is_keyed_to_the_echo_not_the_row():
    sets, t = _echo_sets(10)
    l2 = BaselineInstrument(quantise=False, thermal_noise=False,
                            fading=FadingModel(tau_i_s=0.5, tau_log_sigma=0.0))
    g1 = l2.fading_gains(sets, t, np.random.default_rng(3))
    g2 = l2.fading_gains(sets, t, np.random.default_rng(3))
    assert set(g1) == {"f1"} and np.array_equal(g1["f1"].gain, g2["f1"].gain)


def test_a_diffuse_part_is_refused_here():
    sets, t = _echo_sets(3)
    l2 = BaselineInstrument(fading=FadingModel(tau_i_s=0.5, k_factor=3.0, tau_d_s=0.05))
    with pytest.raises(ValueError, match="diffuse"):
        l2.receive_sequence(sets, sweep_law.drift_mode(F), t, np.random.default_rng(0))


def test_layer_echo_beside_the_facet():
    alt = np.linspace(80e3, 600e3, 521)
    st = IonosphereState(alt_m=alt, ne=1.2e12 * np.exp(-((alt - 300e3) / 80e3) ** 2) + 1e9)
    tier = build("l1", "facet", facets=[{"id": "f1", "depth_m": 200e3}], layer_echo=True)
    e = tier.echoes(st, PlatformState(alt_m=600e3, lat_deg=0, lon_deg=0, t_s=0.0), np.array([F]),
                    np.random.default_rng(0))
    flags = list(e.flags)
    assert "facet:f1" in flags and "layer:O" in flags
    i = flags.index("layer:O")
    assert np.isfinite(e.group_delay_s[i]) and e.phase_delay_s[i] < e.group_delay_s[i]
    sig = BaselineInstrument(quantise=False, thermal_noise=False).receive_sequence(
        [e], sweep_law.drift_mode(F), np.array([0.0]), np.random.default_rng(0))[0]
    assert sig.flags["rows_used"] == 2
