"""WP5 Phase 2: the drift-mode sequence path and its exit criterion.

Exit (WP5 plan, P2): a noise-free sequence through L1 (facet) and L2 reproduces
the plate_echo phase to 0.01 rad and the Doppler to 1e-3.  The phase is taken
from the L2 record of each sweep at its centre, with the beat terms that the
group delay sets removed; it is compared with the field of ``echo_stratified``
(the verified model, called directly) at the same instant.
"""

from __future__ import annotations

from dataclasses import replace

import numpy as np
import pytest

from tss import sweep_law
from tss.interfaces import PlatformState, ReceivedSequence
from tss.l1_propagation import plate_echo as pe
from tss.l1_propagation.doppler import doppler_shift_hz
from tss.l2_instrument.baseline import BaselineInstrument
from tss.pipeline import Pipeline
from tss.registry import build

F = 5.0e6
PROF = pe.ExponentialTopside()
DEPTH = PROF.depth_for_mu(F, 0.5)
DRIFT = 25.0


def _pipeline(**l2kw):
    tiers = {
        "l0": build("l0", "chapman"),
        "l1": build("l1", "facet", facets=[{"id": "f1", "depth_m": DEPTH, "drift_mps": DRIFT}]),
        "l2": BaselineInstrument(quantise=False, thermal_noise=False, **l2kw),
    }
    return Pipeline(tiers, seed=7)


def _run(n=91, t0=-1.0, **l2kw):
    sweep = sweep_law.drift_mode(F)
    t = sweep_law.drift_mode_times(n, t0_s=t0)
    plats = [PlatformState(alt_m=600e3, lat_deg=0.0, lon_deg=0.0) for _ in t]
    res = _pipeline(**l2kw).run_sequence(
        "seq", sweep=sweep, t_s=t, platforms=plats, alt_m=np.linspace(100e3, 600e3, 101),
        epoch="2026-01-01T12:00:00Z",
    )
    return sweep, t, res


def _carrier_phase(rec, echo, sweep):
    """Phase 2 pi f tau_p at the record centre, from the samples."""
    x = rec.samples[0]
    fs = rec.fs_hz
    tt = np.arange(x.size) / fs
    tc = 0.5 * x.size / fs
    tau = float(echo.group_delay_s[0])
    fd = float(echo.doppler_hz[0])
    k = sweep.k_f_hz_per_s
    live = np.abs(x) > 0
    z = np.sum(x[live] * np.exp(-2j * np.pi * (k * tau - fd) * (tt[live] - tc)))
    f0 = sweep.f_start_hz
    beat = 2 * np.pi * (f0 * tau + k * tau * tc - k * tau**2 / 2 - F * tau)
    # What is left is 2 pi f tau_p(t_n) - 2 pi f_d tc ~ 2 pi f tau_p(t_n + tc).
    return float(np.angle(z)) - beat, tc


def test_run_sequence_returns_a_received_sequence():
    sweep, t, res = _run(n=5)
    seq = res.sequence
    assert isinstance(seq, ReceivedSequence)
    assert len(seq.records) == 5 and np.array_equal(seq.t_s, t)
    assert [p.t_s for p in seq.platforms] == pytest.approx(list(t + sweep.sweep_period_s / 2))
    assert all(r.flags["first_delay_blanked"] for r in seq.records)
    assert len(res.echoes) == 5
    assert res.flags["l2_sequence"] is True


def test_exit_noise_free_phase_reproduces_plate_echo_to_0p01_rad():
    sweep, t, res = _run()
    plate = pe.Plate(depth_m=DEPTH, drift_mps=DRIFT)
    errs, mags = [], []
    for rec, echo, tn in zip(res.sequence.records, res.echoes, t, strict=True):
        assert rec.flags["n_clipped"] == 0
        phase, tc = _carrier_phase(rec, echo, sweep)
        ref = pe.echo_stratified(np.array([tn + tc]), F, plate, PROF)
        errs.append(np.angle(np.exp(1j * (phase + np.angle(ref.field[0])))))
        mags.append(abs(ref.field[0]))
    errs, mags = np.abs(errs), np.array(mags)
    # The whole window, to the first nulls (+/-1.05 s here).  L1 is evaluated at
    # each sweep's centre and L2 holds the Doppler fixed about it, so the error
    # is second order about the centre (1e-9 rad measured, 2026-11-25); with
    # L1 at the sweep start it reached 0.02 rad near the nulls.
    lobe = mags >= 0.1 * mags.max()
    assert lobe.sum() >= 60
    assert np.max(errs) < 0.01


def test_exit_doppler_to_1e3():
    sweep, t, res = _run(n=21, t0=-10 / 45 - 1 / 90)
    i0 = 10  # the sweep centred on overflight
    mu = float(np.sqrt(PROF.mu_squared(F, 0.0)))
    want = float(doppler_shift_hz(DRIFT, F, mu_sat=mu))
    # The echo set's Doppler at the overflight sweep.
    assert res.echoes[i0].doppler_hz[0] == pytest.approx(want, rel=1e-3)
    # And the Doppler from the L2 phase steps, symmetric about overflight so
    # the platform's own term cancels.
    ph = np.unwrap([_carrier_phase(r, e, sweep)[0] for r, e in zip(res.sequence.records, res.echoes,
                                                                    strict=True)])
    T = sweep.sweep_period_s
    est = -(ph[i0 + 1] - ph[i0 - 1]) / (2 * np.pi * 2 * T)
    tc = 0.5 * res.sequence.records[0].samples.shape[1] / res.sequence.records[0].fs_hz
    # The pair is centred on the record centre tc, within a sample of the
    # overflight sweep's centre T/2; the platform term there is (tc - T/2) times the slope.
    slope = 2 * F * 7600.0**2 / (299_792_458.0 * float(res.echoes[i0].group_delay_s[0]) * 299_792_458.0 / 2)
    assert est + slope * (tc - T / 2) == pytest.approx(want, rel=1e-3)


def test_exit_holds_on_the_orbit_in_the_l0_profile():
    """The same check on the geometry Andrew chose for Phase 2: an orbit and the
    L0 profile.  The reference is the facet tier itself, evaluated at the record
    centre, since no independent spherical model exists."""
    from tss.interfaces import IonosphereState
    from tss.platform.orbits import CircularOrbit, platform_state_at

    orbit = CircularOrbit(600e3, 60.0, u0_deg=20.0)
    p0 = platform_state_at(orbit, 0.0)
    alt = np.linspace(50e3, 700e3, 1301)
    state = IonosphereState(alt_m=alt, ne=(1.5e6**2 * np.exp((600e3 - alt) / 200e3)) / 80.6163)

    class _L0:
        name = "fixed"
        layer = "l0"

        def state_at(self, lat, lon, epoch, alt_m):
            return state

        def describe(self):
            return {"name": "fixed"}

    facet = build("l1", "facet", geometry="orbit", profile="l0", ref_alt_m=600e3,
                  facets=[{"id": "o", "lat_deg": p0.lat_deg, "lon_deg": p0.lon_deg, "alt_m": 300e3,
                           "drift_mps": DRIFT}])
    p = Pipeline({"l0": _L0(), "l1": facet,
                  "l2": BaselineInstrument(quantise=False, thermal_noise=False)}, seed=3)
    sweep = sweep_law.drift_mode(F)
    t = sweep_law.drift_mode_times(61, t0_s=-30 / 45)
    plats = [platform_state_at(orbit, float(x) + 1 / 90) for x in t]
    res = p.run_sequence("orb", sweep=sweep, t_s=t, platforms=plats, alt_m=alt,
                         epoch="2026-01-01T12:00:00Z", truth_at=(p0.lat_deg, p0.lon_deg))
    errs, mags = [], []
    for rec, echo, tn in zip(res.sequence.records, res.echoes, t, strict=True):
        phase, tc = _carrier_phase(rec, echo, sweep)
        pc = replace(platform_state_at(orbit, float(tn + tc)), t_s=float(tn + tc))
        ref = facet.echoes(state, pc, np.array([F]), np.random.default_rng(0))
        errs.append(abs(np.angle(np.exp(1j * (phase - 2 * np.pi * F * ref.phase_delay_s[0])))))
        mags.append(10 ** (ref.amplitude_db[0] / 20))
    errs, mags = np.array(errs), np.array(mags)
    assert mags.max() > 0 and np.max(errs) < 0.01


def test_truth_is_taken_once_where_asked():
    sweep = sweep_law.drift_mode(F)
    t = sweep_law.drift_mode_times(2)
    plats = [PlatformState(alt_m=600e3, lat_deg=10.0, lon_deg=20.0)] * 2
    res = _pipeline().run_sequence("s", sweep=sweep, t_s=t, platforms=plats,
                                   alt_m=np.linspace(100e3, 600e3, 101), epoch="2026-01-01T12:00:00Z",
                                   truth_at=(11.0, 21.0))
    assert res.flags["truth_at"] == [11.0, 21.0]
    assert res.truth.lat_deg == 11.0 and res.truth.lon_deg == 21.0


def test_per_sweep_fallback_for_an_l2_without_sequences():
    class _Plain:
        name = "plain"
        layer = "l2"

        def __init__(self):
            self.inner = BaselineInstrument(quantise=False, thermal_noise=False)

        def receive(self, echoes, sweep, rng):
            return self.inner.receive(echoes, sweep, rng)

        def describe(self):
            return {"name": "plain"}

    p = _pipeline()
    p.tiers["l2"] = _Plain()
    sweep = sweep_law.drift_mode(F)
    t = sweep_law.drift_mode_times(3)
    res = p.run_sequence("s", sweep=sweep, t_s=t,
                         platforms=[PlatformState(alt_m=600e3, lat_deg=0.0, lon_deg=0.0)] * 3,
                         alt_m=np.linspace(100e3, 600e3, 101), epoch="2026-01-01T12:00:00Z")
    assert res.flags["l2_sequence"] is False
    assert len(res.sequence.records) == 3
