"""The WP3 leakage model (Phase 1) against its analytic forms.

The Phase 1 exit criterion: the model reproduces the single-tone analytics of
Rev G section 5.9 to 0.5 dB.  Each spectral test measures a tone in the
spectrum of the de-chirped leakage, the way the instrument would, rather than
reading the model's own phase back.
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from tss import sweep_law
from tss.config import load_yaml, resolve_path
from tss.constants import C_LIGHT
from tss.interfaces import SweepSchedule
from tss.l2_instrument import leakage as lk
from tss.l2_instrument import reference_echo as rx
from tss.l2_instrument.baseline import BaselineInstrument
from tss.studies.data_volume import A1

TOL_DB = 0.5


@pytest.fixture(scope="module")
def inst():
    return BaselineInstrument.from_config(load_yaml(resolve_path("instrument/baseline_6u.yaml")))


def _uniform(f_start=5.0e6, k=5.6e6, seconds=0.25):
    return SweepSchedule(
        f_start_hz=f_start, f_stop_hz=f_start + k * seconds, k_f_hz_per_s=k,
        sweep_period_s=seconds, direction=1,
    )


def _tone_dbc(x: np.ndarray, fs: float, f_tone: float) -> float:
    """Level of the tone nearest ``f_tone`` against the total power, Hann window.

    The carrier is within a fraction of a hertz of zero and holds all but a
    negligible part of the power, so the total power is the carrier.  The
    record is chosen so the tone sits on a bin, which makes the window's
    scalloping zero.
    """
    n = x.size
    w = np.hanning(n)
    spec = np.fft.fft(x * w)
    k = int(round(f_tone * n / fs)) % n
    return float(20.0 * np.log10(np.abs(spec[k]) / (np.sum(w) * np.sqrt(np.mean(np.abs(x) ** 2)))))


# ----------------------------------------------------------------------------
# The profile
# ----------------------------------------------------------------------------


def test_flat_profile_is_the_isolation_everywhere():
    m = lk.LeakageModel(isolation_db=40.0)
    f = np.linspace(0.33e6, 20e6, 50)
    assert np.allclose(lk.isolation_profile_db(f, m), 40.0)
    assert np.allclose(np.angle(lk.coupling_response(f, m)), 0.0)


def test_ripple_is_peak_to_peak_about_the_base():
    m = lk.LeakageModel(isolation_db=40.0, isolation_ripple_db=6.0, ripple_period_hz=1.0e6)
    iso = lk.isolation_profile_db(np.linspace(0.33e6, 20e6, 20001), m)
    assert iso.max() == pytest.approx(43.0, abs=1e-3)
    assert iso.min() == pytest.approx(37.0, abs=1e-3)


def test_a_resonance_lowers_the_isolation_by_its_depth_and_has_its_q():
    r = lk.Resonance(f0_hz=7.0e6, q=50.0, depth_db=20.0)
    m = lk.LeakageModel(isolation_db=50.0, resonances=(r,))
    assert lk.isolation_profile_db(7.0e6, m) == pytest.approx(30.0, abs=1e-9)
    assert lk.isolation_profile_db(2.0e6, m) == pytest.approx(50.0, abs=0.1)
    # The resonant term alone falls to half power at f0 (1 +/- 1/(2Q)).
    g = 10.0
    edge = 7.0e6 * (1.0 + 1.0 / 100.0)
    excess = np.abs(r.response(edge) - 1.0) / (g - 1.0)
    assert excess == pytest.approx(1.0 / math.sqrt(2.0), rel=1e-9)
    # The coupling phase turns through the resonance.
    ph = np.angle(lk.coupling_response(np.array([6.8e6, 7.2e6]), m))
    assert ph[0] > 0.1 and ph[1] < -0.1


def test_a_table_interpolates_in_db_and_holds_its_ends():
    m = lk.LeakageModel.from_table([1e6, 2e6, 4e6], [60.0, 40.0, 50.0])
    f = np.array([0.5e6, 1.5e6, 3e6, 10e6])
    assert np.allclose(lk.isolation_profile_db(f, m), [60.0, 50.0, 45.0, 50.0])


@pytest.mark.parametrize(
    "kwargs",
    [
        dict(isolation_ripple_db=3.0),
        dict(table_f_hz=(1e6, 2e6)),
        dict(table_f_hz=(2e6, 1e6), table_isolation_db=(40.0, 40.0)),
        dict(coupling_delay_s=-1e-9),
        dict(isolation_ripple_db=-1.0, ripple_period_hz=1e6),
    ],
)
def test_bad_models_are_refused(kwargs):
    with pytest.raises(ValueError):
        lk.LeakageModel(**kwargs)


def test_the_g2_bracket_and_the_baseline_delay():
    assert lk.COUPLING_PATHS_M == (1.0, 3.0, 10.0, 30.0)
    assert lk.LeakageModel.from_coupling_path(10.0).coupling_delay_s == pytest.approx(10.0 / C_LIGHT)
    assert lk.LeakageModel().coupling_delay_s == 3.0e-9
    assert lk.ISOLATION_BRACKETS_DB == (30.0, 40.0, 50.0, 60.0)


def test_the_model_reads_the_instrument_configuration():
    cfg = load_yaml(resolve_path("instrument/baseline_6u.yaml"))
    m = lk.LeakageModel.from_config(cfg["isolation"])
    assert m.isolation_db == 40.0 and m.coupling_delay_s == 3.0e-9


# ----------------------------------------------------------------------------
# Level and phase
# ----------------------------------------------------------------------------


def test_the_leakage_carries_transmit_power_less_the_isolation():
    sweep = A1.schedule()
    fs = 78.4e3
    for iso in lk.ISOLATION_BRACKETS_DB:
        x = lk.leakage_signal(sweep, fs, lk.LeakageModel(isolation_db=iso), np.random.default_rng(0),
                              n_samples=4000, t0_s=1.0)
        assert 10 * np.log10(np.mean(np.abs(x) ** 2)) == pytest.approx(0.0 - iso, abs=1e-6)


def test_leakage_over_echo_reproduces_the_phase0_reference(inst):
    sweep = A1.schedule()
    flat = lk.leakage_over_echo_db(inst, sweep, lk.LeakageModel(isolation_db=40.0))
    ref = rx.reference_levels(inst, sweep)
    for key in ("leakage_over_rmax_echo_db", "leakage_over_weakest_db", "leakage_over_noise_db"):
        assert np.allclose(flat[key], ref[key], atol=1e-9)
    f = flat["f_hz"] / 1e6
    i5 = int(np.argmin(np.abs(f - 5.0)))
    assert flat["leakage_over_rmax_echo_db"][i5] == pytest.approx(97.2, abs=0.1)
    # The binding step of INS-22 (Rev G section 5.9, amended 2026-10-27).
    ib = int(np.argmin(np.abs(f - 1.2628)))
    assert flat["leakage_over_rmax_echo_db"][ib] == pytest.approx(120.55, abs=0.05)


def test_sampled_leakage_over_the_reference_echo_matches_to_half_a_db(inst):
    """The plan's test: the model's leakage, measured from its samples, over the
    Phase 0 echo, against the analytic ratio, at every bracket and three steps."""
    sweep = A1.schedule()
    fs = inst.sample_rate_for(sweep)
    for f_mhz in (1.2628, 5.0, 15.0):
        t0 = float(sweep_law.time_at_frequency_s(sweep, f_mhz * 1e6))
        lev = rx.reference_levels(inst, sweep, f_hz=np.array([f_mhz * 1e6]))
        for iso in lk.ISOLATION_BRACKETS_DB:
            x = lk.leakage_signal(sweep, fs, lk.LeakageModel(isolation_db=iso), np.random.default_rng(1),
                                  n_samples=1000, t0_s=t0)
            measured = 10 * np.log10(np.mean(np.abs(x) ** 2)) - float(lev["echo_rmax_dbw"][0])
            expected = float(lev["leakage_over_rmax_echo_db"][0]) + 40.0 - iso
            assert measured == pytest.approx(expected, abs=TOL_DB)


def test_the_phase_is_the_exact_dechirped_law_across_a_step_transition():
    sweep = A1.schedule()
    fs = 2.0e6
    tau = 30.0 / C_LIGHT
    dwell = float(sweep.step_dwell_s)
    t0 = 20 * dwell - 50.0 / fs                      # straddles the 20th step boundary
    m = lk.LeakageModel.from_coupling_path(30.0, isolation_db=0.0)
    x = lk.leakage_signal(sweep, fs, m, np.random.default_rng(0), n_samples=400, t0_s=t0)
    t = t0 + np.arange(400) / fs
    exact = sweep_law.phase_rad(sweep, t) - sweep_law.phase_rad(sweep, t - tau)
    err = np.angle(x * np.exp(-1j * exact))
    assert np.max(np.abs(err)) < 1e-6
    # The beat is K_i tau inside each step.
    k_before = sweep_law.rate_hz_per_s(sweep, t[0])
    k_after = sweep_law.rate_hz_per_s(sweep, t[-1])
    assert k_after > k_before
    ph = np.unwrap(np.angle(x))
    beat_before = (ph[40] - ph[0]) / (40 / fs) / (2 * np.pi)
    beat_after = (ph[399] - ph[300]) / (99 / fs) / (2 * np.pi)
    assert beat_before == pytest.approx(k_before * tau, rel=1e-3)
    assert beat_after == pytest.approx(k_after * tau, rel=1e-3)


def test_nothing_arrives_before_the_coupling_delay():
    sweep = A1.schedule()
    m = lk.LeakageModel.from_coupling_path(30.0)
    x = lk.leakage_signal(sweep, 1.0e9, m, np.random.default_rng(0), n_samples=200)
    assert np.all(x[:101] == 0) and np.all(x[101:] != 0)


# ----------------------------------------------------------------------------
# Range correlation (Rev G section 5.9)
# ----------------------------------------------------------------------------

_FS = 200.0e3
_T = 0.05                       # 20 Hz bins
_N = int(_FS * _T)


@pytest.mark.parametrize("path_m", lk.COUPLING_PATHS_M)
@pytest.mark.parametrize("f_m", [3740.0, 20000.0, 39220.0])
def test_a_transmit_spur_is_suppressed_by_4_sin2(path_m, f_m):
    dev = 200.0
    fm = lk.SingleToneFM(deviation_hz=dev, rate_hz=f_m)
    m = lk.LeakageModel.from_coupling_path(path_m, isolation_db=40.0)
    sweep = _uniform(seconds=0.1)
    x = lk.leakage_signal(sweep, _FS, m, np.random.default_rng(0), n_samples=_N, t0_s=0.01,
                          tx_phase_error=fm)
    tau = path_m / C_LIGHT
    sideband = _tone_dbc(x, _FS, f_m)
    assert sideband == pytest.approx(lk.dechirped_fm_sideband_dbc(dev, f_m, tau), abs=TOL_DB)
    suppression = sideband - lk.fm_spur_dbc(dev, f_m)
    assert suppression == pytest.approx(float(lk.range_correlation_suppression_db(f_m, tau)), abs=TOL_DB)
    # Rev G's small-angle form, (pi Df tau)^2, holds where f_m tau << 1.
    assert sideband == pytest.approx(lk.small_angle_sideband_dbc(dev, tau), abs=TOL_DB)


def test_an_independent_reference_forfeits_the_credit():
    """INS-21: without a common source the spur survives the de-chirp whole."""
    dev, f_m = 200.0, 20000.0
    fm = lk.SingleToneFM(deviation_hz=dev, rate_hz=f_m)
    m = lk.LeakageModel.from_coupling_path(10.0)
    x = lk.leakage_signal(_uniform(seconds=0.1), _FS, m, np.random.default_rng(0), n_samples=_N,
                          t0_s=0.01, tx_phase_error=fm, reference_phase_error=None)
    assert _tone_dbc(x, _FS, f_m) == pytest.approx(lk.fm_spur_dbc(dev, f_m), abs=TOL_DB)
    common = lk.leakage_signal(_uniform(seconds=0.1), _FS, m, np.random.default_rng(0), n_samples=_N,
                               t0_s=0.01, tx_phase_error=fm)
    lost = _tone_dbc(x, _FS, f_m) - _tone_dbc(common, _FS, f_m)
    assert lost == pytest.approx(-float(lk.range_correlation_suppression_db(f_m, 10.0 / C_LIGHT)), abs=TOL_DB)


def test_range_correlation_holds_inside_one_a1_step():
    """The same check on the stepped law, one 26.8 ms dwell at 5 MHz."""
    sweep = A1.schedule()
    fs = 4 * 78.4e3
    dwell = float(sweep.step_dwell_s)
    i = int(sweep_law.step_index(sweep, sweep_law.time_at_frequency_s(sweep, 5.0e6)))
    n = int(round(dwell * fs))
    f_m = 100.0 / (n / fs)                         # on bin 100
    dev = 300.0
    m = lk.LeakageModel.from_coupling_path(10.0)
    x = lk.leakage_signal(sweep, fs, m, np.random.default_rng(0), n_samples=n, t0_s=i * dwell,
                          tx_phase_error=lk.SingleToneFM(dev, f_m))
    expected = lk.dechirped_fm_sideband_dbc(dev, f_m, 10.0 / C_LIGHT)
    assert _tone_dbc(x, fs, f_m) == pytest.approx(expected, abs=TOL_DB)


# ----------------------------------------------------------------------------
# Departures from flat
# ----------------------------------------------------------------------------


@pytest.mark.parametrize("ripple_db", [0.5, 1.0, 3.0])
def test_an_isolation_ripple_becomes_am_sidebands_at_k_over_period(ripple_db):
    k = 5.6e6
    period = 50.0e3
    fs, seconds = 2000.0, 1.0                      # 1 Hz bins
    f_am = lk.ripple_modulation_hz(k, period)      # 112 Hz
    m = lk.LeakageModel(isolation_db=40.0, isolation_ripple_db=ripple_db, ripple_period_hz=period)
    x = lk.leakage_signal(_uniform(k=k, seconds=1.2), fs, m, np.random.default_rng(0),
                          n_samples=int(fs * seconds), t0_s=0.1)
    # The total power includes the sidebands here; refer to the carrier bin instead.
    n = x.size
    w = np.hanning(n)
    spec = np.abs(np.fft.fft(x * w))
    carrier = spec[0]
    side = spec[int(round(f_am * n / fs))]
    assert 20 * np.log10(side / carrier) == pytest.approx(lk.ripple_sideband_dbc(ripple_db), abs=TOL_DB)


def test_the_leakage_amplitude_follows_a_resonance():
    r = lk.Resonance(f0_hz=5.1e6, q=30.0, depth_db=15.0)
    m = lk.LeakageModel(isolation_db=45.0, resonances=(r,))
    sweep = _uniform(seconds=0.1)
    x = lk.leakage_signal(sweep, 10e3, m, np.random.default_rng(0), n_samples=1000)
    t = np.arange(1000) / 10e3
    f_e = 5.0e6 + 5.6e6 * (t - m.coupling_delay_s)
    expected = -lk.isolation_profile_db(f_e, m)
    got = 20 * np.log10(np.abs(x[1:]))
    assert np.allclose(got, expected[1:], atol=1e-6)


def test_delay_jitter_gives_2_pi_f_sigma_of_phase():
    sigma = 2.0e-12
    f = 20.0e6
    fs, n = 20.0e3, 20000
    m = lk.LeakageModel(delay_jitter_s=sigma, jitter_bandwidth_hz=1000.0)
    sweep = _uniform(f_start=f, k=1.0, seconds=10.0)  # 1 Hz/s: the frequency barely moves
    x = lk.leakage_signal(sweep, fs, m, np.random.default_rng(7), n_samples=n, t0_s=1e-3)
    still = lk.leakage_signal(sweep, fs, lk.LeakageModel(), np.random.default_rng(7), n_samples=n,
                               t0_s=1e-3)
    rms = float(np.sqrt(np.mean(np.angle(x / still) ** 2)))
    assert 20 * np.log10(rms) == pytest.approx(20 * np.log10(lk.jitter_phase_rms_rad(f, sigma)), abs=TOL_DB)


def test_jitter_is_reproducible_and_absent_jitter_draws_nothing():
    sweep = A1.schedule()
    m = lk.LeakageModel(delay_jitter_s=1e-12, jitter_bandwidth_hz=10.0)
    a = lk.leakage_signal(sweep, 78.4e3, m, np.random.default_rng(3), n_samples=500, t0_s=1.0)
    b = lk.leakage_signal(sweep, 78.4e3, m, np.random.default_rng(3), n_samples=500, t0_s=1.0)
    assert np.array_equal(a, b)
    rng = np.random.default_rng(3)
    lk.leakage_signal(sweep, 78.4e3, lk.LeakageModel(), rng, n_samples=500, t0_s=1.0)
    assert rng.normal() == np.random.default_rng(3).normal()


def test_a_retune_transient_decays_from_each_step_boundary():
    sweep = A1.schedule()
    dwell = float(sweep.step_dwell_s)
    tc = 2.0e-3
    m = lk.LeakageModel(isolation_db=40.0, retune_transient_db=6.0, retune_time_constant_s=tc)
    fs = 100.0e3
    t0 = 30 * dwell
    x = lk.leakage_signal(sweep, fs, m, np.random.default_rng(0), n_samples=int(dwell * fs), t0_s=t0)
    level = 20 * np.log10(np.abs(x)) + 40.0
    t = np.arange(x.size) / fs - m.coupling_delay_s
    # The first sample left the antenna just before the boundary, in the previous step.
    assert np.allclose(level[1:], 6.0 * np.exp(-t[1:] / tc), atol=1e-6)
    assert level[0] < 1e-3


# ----------------------------------------------------------------------------
# The moment-method anchor (G1(b))
# ----------------------------------------------------------------------------


def test_the_mom_anchor_loads_as_a_tabulated_profile():
    tab = lk.mom_anchor_table()
    m = lk.mom_anchor_model("crossed_offset")
    assert m.table_f_hz is not None
    assert tab["f_hz"][0] == pytest.approx(0.33e6) and tab["f_hz"][-1] == pytest.approx(20e6)
    assert np.allclose(lk.isolation_profile_db(tab["f_hz"], m), tab["iso_crossed_offset_db"])
    with pytest.raises(ValueError):
        lk.mom_anchor_model(geometry="coaxial")


def test_the_mom_dipole_is_the_balanis_dipole_and_coupling_times_q_sets_the_floor():
    """Below a few MHz the coil loss dominates both antennas, so the isolation of
    two coupled resonators is set by k Q and is flat, not by efficiency."""
    tab = lk.mom_anchor_table()
    low = tab["f_hz"] < 2.0e6
    assert np.ptp(tab["iso_crossed_offset_db"][low]) < 0.1
    assert np.ptp(tab["iso_parallel_db"][low]) < 0.1
    # The transmit efficiency the anchor computes is the tss ShortDipole's to 1 dB.
    from tss.l2_instrument.antenna import ShortDipole
    eta = np.asarray(ShortDipole(length_m=5.0).efficiency_db(tab["f_hz"] / 1e6), dtype=float)
    assert np.max(np.abs(eta - tab["eta_tx_db"])) < 1.0


def test_the_bus_anchor_is_the_default_and_obeys_the_coupled_resonator_law():
    bus = lk.mom_anchor_table("bus")
    m = lk.mom_anchor_model()
    assert np.allclose(lk.isolation_profile_db(bus["f_hz"], m), bus["iso_crossed_ends_offset_db"])
    low = bus["f_hz"] < 3.0e6
    for mount in ("crossed_ends_offset", "parallel_ends"):
        predicted = lk.coupled_resonator_isolation_db(bus[f"k_{mount}"][low], 150.0)
        assert np.max(np.abs(predicted - bus[f"iso_{mount}_db"][low])) < TOL_DB
    # The crossed mount on the bus is below INS-23's 40 dB at the bottom of the band.
    assert np.all(bus["iso_crossed_ends_offset_db"][low] < 40.0)
    assert lk.coupled_resonator_isolation_db(1.3e-4, 150.0) == pytest.approx(40.0, abs=0.3)
