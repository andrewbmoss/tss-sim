"""WP5 Phase 3: the drift tier (``l3/drift``), route one at signal level.

Through L1 (facet) and L2 (baseline, quantisation and thermal noise off).
The pass criterion of D5 option A at noise-free level is drift within 1 m/s.
One test adds white noise to the records at 20 dB per sweep and compares the
spread of t_pk with the parameter-level Monte Carlo (docs/drift-tier.md,
section 6).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from tss import sweep_law
from tss.constants import C_LIGHT
from tss.interfaces import DriftProduct, IonosphereState, PlatformState
from tss.l1_propagation import plate_echo as pe
from tss.l2_instrument.baseline import BaselineInstrument
from tss.l3_signal.drift import DriftProcessing, OnboardPrior, effective_speed_ratio
from tss.pipeline import Pipeline
from tss.platform.orbits import RE, CircularOrbit, platform_state_at
from tss.registry import build

F = 5.0e6
PROF = pe.ExponentialTopside()
DEPTH = PROF.depth_for_mu(F, 0.5)
PRIOR = OnboardPrior.exponential(f_n_sat_hz=1.5e6, scale_height_m=200e3, ref_alt_m=600e3, source="truth")


def _straight(facet, n=91):
    tiers = {"l0": build("l0", "chapman"), "l1": build("l1", "facet", facets=[facet]),
             "l2": BaselineInstrument(quantise=False, thermal_noise=False)}
    sweep = sweep_law.drift_mode(F)
    t = sweep_law.drift_mode_times(n, t0_s=-(n // 2) / 45 - 1 / 90)
    plats = [PlatformState(alt_m=600e3, lat_deg=0.0, lon_deg=0.0) for _ in t]
    return Pipeline(tiers, seed=11).run_sequence("d", sweep=sweep, t_s=t, platforms=plats,
                                                  alt_m=np.linspace(100e3, 600e3, 101),
                                                  epoch="2026-01-01T12:00:00Z")


@pytest.fixture(scope="module")
def straight_run():
    return _straight({"id": "f1", "depth_m": DEPTH, "drift_mps": 25.0})


def test_noise_free_straight_track_drift_within_1_mps(straight_run):
    prod = DriftProcessing(prior=PRIOR).process_sequence(straight_run.sequence)
    assert isinstance(prod, DriftProduct)
    assert list(prod.echo_id) == ["e1"]
    p = {k: v[0] for k, v in prod.products.items()}
    assert p["echo_class"] == "compact_facet"
    assert p["drift_mps"] == pytest.approx(25.0, abs=1.0)
    assert p["geometry"] == "flat" and p["v_eff_mps"] == 7600.0
    # The apparent range at the peak: the echo set's group path at overflight.
    i0 = int(np.argmin(np.abs(straight_run.sequence.t_s + 1 / 90)))
    want = C_LIGHT * float(straight_run.echoes[i0].group_delay_s[0]) / 2
    assert p["range_m"] == pytest.approx(want, abs=150.0)
    assert prod.flags["prior_source"] == "truth"
    assert "TILT_UNFLAGGED" in prod.flags["echo_flags"][0]


def test_the_blanked_interval_is_gated_out_and_counted(straight_run):
    prod = DriftProcessing(prior=PRIOR).process_sequence(straight_run.sequence)
    rec = straight_run.sequence.records[0]
    n = rec.samples.shape[1]
    tau = float(straight_run.echoes[45].group_delay_s[0])
    assert prod.products["blanked_fraction"][0] == pytest.approx((tau * rec.fs_hz + 3) / n, abs=2 / n)
    # The gate's centre fixes the frequency the velocity is converted at (D8 condition 2).
    sw = rec.sweep
    assert prod.products["f_eff_hz"][0] == pytest.approx(
        sw.f_start_hz + sw.k_f_hz_per_s * (prod.products["t_ref_s"][0] - tau), abs=50.0)


def test_contract_names_are_present(straight_run):
    prod = DriftProcessing(prior=PRIOR).process_sequence(straight_run.sequence)
    for k in ("range_m", "los_rate_mps", "drift_mps", "drift_sigma_mps", "t_pk_s", "t_zc_s",
              "echo_class", "platform_doppler_hz", "slope_ratio", "peak_contrast", "along_track_m"):
        assert k in prod.products
    assert "drift_sigma_mps" in prod.flags["not_done"]


def test_extended_plate_is_declined():
    run = _straight({"id": "L", "depth_m": DEPTH, "drift_mps": 25.0, "length_m": 80e3, "taper": "tukey"},
                    n=61)
    p = DriftProcessing(prior=PRIOR).process_sequence(run.sequence)
    assert p.products["echo_class"][0] == "extended"
    assert math.isnan(p.products["drift_mps"][0])


def test_without_a_prior_no_drift_is_formed(straight_run):
    prod = DriftProcessing().process_sequence(straight_run.sequence)
    assert "NO_PRIOR" in prod.flags["echo_flags"][0]
    assert math.isnan(prod.products["drift_mps"][0]) if "drift_mps" in prod.products else True


def test_with_noise_the_envelope_fit_matches_the_parameter_level_monte_carlo():
    """20 dB per sweep (the budget's definition: the per-sweep complex value at the
    envelope's peak over its noise), +/-2.5 envelope widths, the beat given so the
    estimator is isolated.  The spread of t_pk is the route_one Monte Carlo's
    (9.4 ms at this width) within a factor of 1.6.  Before the field was scaled
    to a unit peak it was 66 ms: the fit stopped on absolute tolerances."""
    import dataclasses

    from tss.studies import route_one as r1

    width = 0.926  # r1.envelope_width_s(5 MHz, 529.8 km, 2 km)
    n = int(2 * 2.5 * width * 45) | 1
    run = _straight({"id": "f1", "depth_m": DEPTH, "drift_mps": 25.0}, n=n)
    seq = run.sequence
    tier = DriftProcessing(prior=PRIOR)
    x = np.stack([q.samples[0] for q in seq.records])
    rec = seq.records[0]
    fb = tier._detect(x, rec.fs_hz, rec.sweep.k_f_hz_per_s)[0]
    tt = np.arange(x.shape[1]) / rec.fs_hz
    i0 = int(math.ceil(fb / rec.sweep.k_f_hz_per_s * rec.fs_hz)) + 3
    z = x[:, i0:] @ np.exp(-2j * math.pi * fb * (tt[i0:] - tt[i0:].mean()))
    sig2 = float(np.max(np.abs(z) ** 2)) / ((x.shape[1] - i0) * 100.0)
    tpk, drift = [], []
    for seed in range(20):
        rng = np.random.default_rng(500 + seed)
        recs = [dataclasses.replace(q, samples=q.samples + math.sqrt(sig2 / 2) * (
            rng.standard_normal(q.samples.shape) + 1j * rng.standard_normal(q.samples.shape)))
            for q in seq.records]
        t2 = DriftProcessing(prior=PRIOR)
        t2._detect = lambda *a, fb=fb: [fb]  # noqa: E731
        p = t2.process_sequence(dataclasses.replace(seq, records=recs)).products
        tpk.append(float(p["t_pk_s"][0]))
        drift.append(float(p["drift_mps"][0]))
    tpk = np.array(tpk)
    spread = 1.4826 * np.median(np.abs(tpk - np.median(tpk)))
    want = r1.overflight_time_sigma_s(width, snr_db=20.0, n_trials=200)
    assert spread < 1.6 * want
    assert np.median(drift) == pytest.approx(25.0, abs=1.0)


def _with_noise(seq, sig2, seed, signal=True):
    import dataclasses

    rng = np.random.default_rng(seed)
    recs = [dataclasses.replace(q, samples=(q.samples if signal else 0.0 * q.samples) + math.sqrt(sig2 / 2) * (
        rng.standard_normal(q.samples.shape) + 1j * rng.standard_normal(q.samples.shape))) for q in seq.records]
    return dataclasses.replace(seq, records=recs)


def _sigma2_for(seq, snr_db):
    """Per-sample noise variance giving ``snr_db`` per sweep, the budget's definition."""
    x = np.stack([q.samples[0] for q in seq.records])
    rec = seq.records[0]
    fb = DriftProcessing(detector="summed")._detect(x, rec.fs_hz, rec.sweep.k_f_hz_per_s)[0]
    tt = np.arange(x.shape[1]) / rec.fs_hz
    i0 = int(math.ceil(fb / rec.sweep.k_f_hz_per_s * rec.fs_hz)) + 3
    z = x[:, i0:] @ np.exp(-2j * math.pi * fb * (tt[i0:] - tt[i0:].mean()))
    return float(np.max(np.abs(z) ** 2)) / ((x.shape[1] - i0) * 10 ** (snr_db / 10)), fb


def test_windowed_detector_false_alarm_rate_is_the_stated_one(straight_run):
    """Noise only: the fraction of (window, range bin) cells over threshold is
    the stated P_fa (1e-3 here, for counts), within a factor of 1.5."""
    seq = straight_run.sequence
    over = cells = 0
    for seed in range(10):
        tier = DriftProcessing(prior=PRIOR, pfa=1e-3)
        prod = tier.process_sequence(_with_noise(seq, 1.0, 900 + seed, signal=False))
        d = prod.flags["detection"]
        assert d["method"] == "windowed" and not d["noise_free"]
        over += d["n_cells_over"]
        cells += d["n_cells"]
    assert 1e-3 / 1.5 < over / cells < 1.5e-3


def test_windowed_detector_finds_the_facet_at_10_db_per_sweep(straight_run):
    """10 dB per sweep (budget definition), 21-sweep window, 1e-7 per cell:
    found in 10 of 10 (the summed-spectrum detector found 0 of 5)."""
    seq = straight_run.sequence
    sig2, fb = _sigma2_for(seq, 10.0)
    k = seq.records[0].sweep.k_f_hz_per_s
    for seed in range(10):
        prod = DriftProcessing(prior=PRIOR, pfa=1e-7).process_sequence(_with_noise(seq, sig2, 300 + seed))
        beats = 2.0 * k * np.asarray(prod.products["range_beat_m"], dtype=float) / C_LIGHT
        assert np.min(np.abs(beats - fb)) < 2 * 66.0


def test_detection_report_is_in_echo_order(straight_run):
    """At 1e-4 per cell a false echo usually takes the second slot; the
    per-echo entries of the report follow the echoes (by range), so the
    facet's is the one at its beat, and it is the strongest."""
    seq = straight_run.sequence
    sig2, fb = _sigma2_for(seq, 10.0)
    k = seq.records[0].sweep.k_f_hz_per_s
    seen_two = False
    for seed in range(5):
        prod = DriftProcessing(prior=PRIOR, pfa=1e-4).process_sequence(_with_noise(seq, sig2, 700 + seed))
        d = prod.flags["detection"]
        beats = 2.0 * k * np.asarray(prod.products["range_beat_m"], dtype=float) / C_LIGHT
        assert len(d["excess_db"]) == len(beats) == len(d["window_start_sweep"])
        assert list(beats) == sorted(beats)
        i = int(np.argmin(np.abs(beats - fb)))
        assert d["excess_db"][i] == max(d["excess_db"])
        seen_two |= len(beats) == 2
    assert seen_two


def test_default_false_alarm_rate_is_stated_per_sequence(straight_run):
    """Decisions (m) and (q), 2026-11-27: by default the rate is 1e-2 per
    sequence, divided over the sequence's (window, range bin) cells.  Noise
    only: the report says so, and of 30 sequences at most 2 declare any cell
    (the union bound gives an expectation of at most 0.3; at 1e-4 per cell
    nearly every sequence declares one)."""
    seq = straight_run.sequence
    with_any = 0
    for seed in range(30):
        tier = DriftProcessing(prior=PRIOR)
        assert tier.describe()["pfa_per_sequence"] == 1e-2 and tier.describe()["pfa"] is None
        d = tier.process_sequence(_with_noise(seq, 1.0, 1300 + seed, signal=False)).flags["detection"]
        assert d["pfa_per_sequence"] == 1e-2
        assert d["pfa_per_cell"] * d["n_cells"] == pytest.approx(1e-2, rel=1e-12)
        with_any += d["n_cells_over"] > 0
    assert with_any <= 2


def test_a_range_gate_raises_the_per_cell_rate_and_keeps_the_sequence_rate(straight_run):
    """Fewer range bins searched, the same chance per sequence (decision (q))."""
    seq = _with_noise(straight_run.sequence, 1.0, 1400, signal=False)
    tier = DriftProcessing(prior=PRIOR)
    wide = tier.process_sequence(seq).flags["detection"]
    gated = tier.process_sequence(seq, max_range_m=600e3).flags["detection"]
    assert gated["n_cells"] < wide["n_cells"]
    assert gated["pfa_per_cell"] > wide["pfa_per_cell"]
    assert gated["pfa_per_cell"] * gated["n_cells"] == pytest.approx(1e-2, rel=1e-12)
    assert gated["threshold_factor"] < wide["threshold_factor"]


def test_a_per_cell_rate_overrides_the_sequence_rate(straight_run):
    seq = _with_noise(straight_run.sequence, 1.0, 1401, signal=False)
    d = DriftProcessing(prior=PRIOR, pfa=1e-3).process_sequence(seq).flags["detection"]
    assert d["pfa_per_cell"] == 1e-3 and d["pfa_per_sequence"] is None
    for bad in ({"pfa_per_sequence": 0.0}, {"pfa_per_sequence": 1.0}, {"pfa": 0.0}):
        with pytest.raises(ValueError):
            DriftProcessing(**bad)


def test_default_detector_finds_the_facet_at_10_db_with_no_second_echo(straight_run):
    """10 dB per sweep (budget definition) at the default rate: the facet in 10
    of 10, and the second echo slot empty in at least 9 (at 1e-4 per cell it
    held a false echo in 196 of 200, docs/drift-tier.md section 8)."""
    seq = straight_run.sequence
    sig2, fb = _sigma2_for(seq, 10.0)
    k = seq.records[0].sweep.k_f_hz_per_s
    single = 0
    for seed in range(10):
        prod = DriftProcessing(prior=PRIOR).process_sequence(_with_noise(seq, sig2, 1500 + seed))
        beats = 2.0 * k * np.asarray(prod.products["range_beat_m"], dtype=float) / C_LIGHT
        assert beats.size and np.min(np.abs(beats - fb)) < 2 * 66.0
        single += beats.size == 1
    assert single >= 9


# -- part 3 step 2: the joint estimators in the tier -----------------------------------

JOINT = ("joint", "joint_fading")


def test_the_default_estimator_is_the_pinned_one(straight_run):
    """Decision (t): pinned until the joint fit meets its exit."""
    tier = DriftProcessing(prior=PRIOR)
    assert tier.describe()["estimator"] == "pinned"
    prod = tier.process_sequence(straight_run.sequence)
    assert prod.flags["estimator"] == "pinned"
    assert math.isnan(prod.products["drift_sigma_mps"][0]) and "drift_sigma_mps" in prod.flags["not_done"]
    assert "drift_pinned_mps" not in prod.products
    with pytest.raises(ValueError):
        DriftProcessing(estimator="other")


@pytest.mark.parametrize("estimator", JOINT)
def test_joint_noise_free_straight_track_drift_within_1_mps(straight_run, estimator):
    pinned = DriftProcessing(prior=PRIOR).process_sequence(straight_run.sequence).products
    prod = DriftProcessing(prior=PRIOR, estimator=estimator).process_sequence(straight_run.sequence)
    p = {k: v[0] for k, v in prod.products.items()}
    assert prod.flags["estimator"] == estimator and p["echo_class"] == "compact_facet"
    assert p["drift_mps"] == pytest.approx(25.0, abs=1.0)
    assert p["drift_pinned_mps"] == pytest.approx(float(pinned["drift_mps"][0]), rel=1e-12)
    assert p["joint_ok"] == 1.0 and p["drift_sigma_mps"] > 0.0
    assert p["t_pk_pinned_s"] == pytest.approx(float(pinned["t_pk_s"][0]), rel=1e-12)
    assert p["along_track_m"] == pytest.approx(p["v_eff_mps"] * p["t_pk_s"], rel=1e-12)
    # The identifiability statistics are the pinned test's (decision (s)).
    assert p["slope_ratio"] == pytest.approx(float(pinned["slope_ratio"][0]), rel=1e-12)
    assert "drift_sigma_mps" not in prod.flags["not_done"]
    fading = estimator == "joint_fading"
    # The fitted exponential time has its own name; tau_c_s is docs/drift-mode.md's
    # lag-one product, not formed yet (Andrew, 2026-11-27).
    assert ("tau_c_fit_s" in prod.products) == fading
    assert "tau_c_s" not in prod.products and "tau_c_s" in prod.flags["not_done"]
    assert "TAU_LOWER_BOUND" not in prod.flags["echo_flags"][0]
    if fading:
        # A constant gain over the pass: the fitted time is only a lower bound.
        assert "TAU_FIT_LOWER_BOUND" in prod.flags["echo_flags"][0]


@pytest.mark.parametrize("estimator", JOINT)
def test_joint_orbit_noise_free_drift_within_1_mps(orbit_run, estimator):
    p = DriftProcessing(prior=PRIOR, estimator=estimator).process_sequence(orbit_run.sequence).products
    assert p["joint_ok"][0] == 1.0
    assert p["drift_mps"][0] == pytest.approx(50.0, abs=1.0)


@pytest.mark.parametrize("estimator", JOINT)
def test_joint_extended_plate_is_declined_by_the_pinned_gate(estimator):
    run = _straight({"id": "L", "depth_m": DEPTH, "drift_mps": 25.0, "length_m": 80e3, "taper": "tukey"},
                    n=61)
    p = DriftProcessing(prior=PRIOR, estimator=estimator).process_sequence(run.sequence)
    assert p.products["echo_class"][0] == "extended"
    assert math.isnan(p.products["drift_mps"][0])
    assert "joint_ok" not in p.products or math.isnan(p.products["joint_ok"][0])


def test_the_tiers_noise_variance_is_the_gated_sums(straight_run):
    """The fading fit's sigma^2 from the detector's negative-beat noise: the
    per-sweep value sums the gated samples, then the field is scaled to a unit
    peak (docs/drift-tier.md section 11)."""
    seq = straight_run.sequence
    sig2, _ = _sigma2_for(seq, 20.0)
    n = seq.records[0].samples.shape[1]
    ratio = []
    for seed in range(5):
        p = DriftProcessing(prior=PRIOR, estimator="joint_fading").process_sequence(
            _with_noise(seq, sig2, 1600 + seed)).products
        n_gate = n * (1.0 - float(p["blanked_fraction"][0]))
        want = n_gate * sig2 / 10.0 ** (float(p["power_peak_db"][0]) / 10.0)
        ratio.append(float(p["field_noise_var"][0]) / want)
    assert np.mean(ratio) == pytest.approx(1.0, abs=0.02)


def test_joint_falls_back_to_the_pinned_drift_and_says_so(straight_run, monkeypatch):
    from tss.l3_signal import drift as drift_module
    from tss.l3_signal.novel.route_one_joint import JointEstimate

    monkeypatch.setattr(drift_module, "estimate_route_one_joint", lambda *a, **k: JointEstimate(False, "no"))
    p = DriftProcessing(prior=PRIOR, estimator="joint").process_sequence(straight_run.sequence)
    assert p.products["joint_ok"][0] == 0.0
    assert p.products["drift_mps"][0] == p.products["drift_pinned_mps"][0]
    assert math.isnan(p.products["drift_sigma_mps"][0])
    assert any(f.startswith("JOINT_FAILED_PINNED_USED") for f in p.flags["echo_flags"][0])


def test_fading_fit_without_a_noise_reference_uses_the_pinned_drift(straight_run):
    p = DriftProcessing(prior=PRIOR, estimator="joint_fading", detector="summed").process_sequence(
        straight_run.sequence)
    assert p.products["joint_ok"][0] == 0.0
    assert p.products["drift_mps"][0] == p.products["drift_pinned_mps"][0]
    assert any("no noise" in f for f in p.flags["echo_flags"][0])


# -- part 3 step 2b: the ambiguity check ----------------------------------------------


def test_the_ambiguity_check_is_off_by_default_and_reports_the_margin(straight_run):
    tier = DriftProcessing(prior=PRIOR, estimator="joint_fading")
    assert tier.describe()["ambiguity_margin"] is None
    p = tier.process_sequence(straight_run.sequence).products
    assert p["ambiguity_nll_margin"][0] > 5.0 and "drift_alt_mps" in p
    assert p["drift_mps"][0] == pytest.approx(25.0, abs=1.0)


def test_the_ambiguity_check_passes_an_unambiguous_echo(straight_run):
    prod = DriftProcessing(prior=PRIOR, estimator="joint_fading", ambiguity_margin=5.0).process_sequence(
        straight_run.sequence)
    assert prod.products["drift_mps"][0] == pytest.approx(25.0, abs=1.0)
    assert not any("ambiguous" in f for f in prod.flags["echo_flags"][0])


@pytest.mark.parametrize(("found", "declined"), [(1.0, True), (3.0, False), (math.inf, False), (math.nan, True)])
def test_a_second_maximum_within_the_margin_declines_the_drift(straight_run, monkeypatch, found, declined):
    """Margin 2: a second maximum 1 below the best declines the echo, 3 does
    not; no second maximum does not; a margin that could not be formed does."""
    from tss.l3_signal import drift as drift_module
    from tss.l3_signal.novel.route_one_joint import JointEstimate

    fit = JointEstimate(True, "converged", drift_mps=30.0, drift_sigma_mps=2.0, t0_s=0.01, t_zc_s=0.27,
                        width_s=0.9, tau_c_fit_s=1.0, nll_margin=found, drift_alt_mps=-40.0)
    monkeypatch.setattr(drift_module, "estimate_route_one_joint_fading", lambda *a, **k: fit)
    prod = DriftProcessing(prior=PRIOR, estimator="joint_fading", ambiguity_margin=2.0).process_sequence(
        straight_run.sequence)
    p = {k: v[0] for k, v in prod.products.items()}
    flags = prod.flags["echo_flags"][0]
    assert any(f.startswith("DECLINED: ambiguous") for f in flags) == declined
    for name in ("drift_mps", "drift_sigma_mps", "t_pk_s", "t_zc_s", "along_track_m"):
        assert math.isnan(p[name]) == declined, name
    # What route one saw stays beside it, and the class is the pinned test's.
    assert math.isfinite(p["drift_pinned_mps"]) and p["echo_class"] == "compact_facet"
    assert p["joint_ok"] == 1.0


def test_the_ambiguity_check_belongs_to_the_fading_fit():
    for bad in ({"estimator": "pinned"}, {"estimator": "joint"}):
        with pytest.raises(ValueError):
            DriftProcessing(ambiguity_margin=2.0, **bad)
    with pytest.raises(ValueError):
        DriftProcessing(estimator="joint_fading", ambiguity_margin=-1.0)


# -- the study's recorded measure (Andrew, 2026-11-27: record, decide after the study) ----


def test_the_power_far_from_the_peak_is_recorded_for_every_echo():
    """Raw power more than 1.5 fitted widths from the fitted overflight time over
    the noise: near 1 plus sidelobes for a facet; a declined echo has it too.  It
    needs the detector's noise reference, and samples that far out (a window of
    +/-2.5 widths has them, +/-1 has none)."""
    seq = _straight({"id": "f1", "depth_m": DEPTH, "drift_mps": 25.0}, n=209).sequence
    sig2, _ = _sigma2_for(seq, 20.0)
    noisy = _with_noise(seq, sig2, 1700)
    p = DriftProcessing(prior=PRIOR).process_sequence(noisy).products
    assert 0.3 < p["far_power_over_noise"][0] < 8.0
    declined = DriftProcessing(prior=PRIOR, min_contrast=1e12).process_sequence(noisy).products
    assert math.isnan(declined["drift_mps"][0])
    assert declined["far_power_over_noise"][0] == p["far_power_over_noise"][0]
    # The summed detector (which misses it in this noise) gives no noise reference.
    summed = DriftProcessing(prior=PRIOR, detector="summed").process_sequence(seq).products
    assert math.isnan(summed["far_power_over_noise"][0])
    short = _with_noise(_straight({"id": "f1", "depth_m": DEPTH, "drift_mps": 25.0}, n=61).sequence, sig2, 1)
    assert math.isnan(DriftProcessing(prior=PRIOR).process_sequence(short).products["far_power_over_noise"][0])


def test_the_power_in_the_windows_edges_is_recorded_for_a_patch_too():
    """A patch's power fit fills the window, so no sample is 1.5 fitted widths
    out; the mean power of the window's outer fifths needs no fit and is there
    for every echo (measured 2026-11-27 on the study's trial run)."""
    facet = _straight({"id": "f1", "depth_m": DEPTH, "drift_mps": 25.0}, n=209).sequence
    patch = _straight({"id": "f1", "depth_m": DEPTH, "drift_mps": 25.0,
                       "halo": {"k": 0.0, "n_points": 20, "seed": 3}}, n=209).sequence
    sig2, _ = _sigma2_for(facet, 20.0)
    a = DriftProcessing(prior=PRIOR).process_sequence(_with_noise(facet, sig2, 1800)).products
    b = DriftProcessing(prior=PRIOR).process_sequence(_with_noise(patch, sig2, 1800)).products
    assert 0.3 < a["edge_power_over_noise"][0] < 8.0
    assert b["edge_power_over_noise"][0] > 4.0 * a["edge_power_over_noise"][0]


def test_range_gate_from_the_trace():
    from types import SimpleNamespace

    from tss.l3_signal.drift import range_gate_from_ionogram

    freq = np.array([4.7e6, 4.85e6, 5.0e6, 5.15e6, 5.3e6, 5.45e6])
    vr = np.array([[300e3, np.nan], [310e3, 305e3], [np.nan, 312e3], [330e3, np.nan], [np.nan, np.nan],
                   [np.nan, np.nan]])
    ion = SimpleNamespace(freq_hz=freq, virtual_range_m=vr)
    gate, rep = range_gate_from_ionogram(ion, 5.0e6, margin_m=10e3)
    # Column 0 interpolates to 320 km; column 1 (X, or a swapped label) gives 312: the larger.
    assert gate == pytest.approx(320e3 + 10e3) and rep["reason"] == "gated"
    assert rep["columns_used"] == [0, 1]
    assert range_gate_from_ionogram(ion, 5.0e6, mode_index=1, margin_m=0.0)[0] == pytest.approx(312e3)
    assert range_gate_from_ionogram(ion, 5.4e6)[0] is None  # above the last detected column
    vr[3, 0] = np.nan
    vr[5, 0] = 360e3
    assert range_gate_from_ionogram(ion, 5.2e6)[0] is None  # a gap of three steps


def test_effective_speed_in_a_thin_plasma_is_the_radius_ratio():
    thin = OnboardPrior.exponential(f_n_sat_hz=1e3, scale_height_m=200e3, ref_alt_m=600e3)
    r_p = RE + 600e3
    ratio, depth = effective_speed_ratio(thin, 5e6, 600e3, 300e3, r_p)
    assert depth == pytest.approx(300e3, rel=1e-6)
    assert ratio == pytest.approx(math.sqrt((r_p - 300e3) / r_p), rel=1e-6)


ORBIT = CircularOrbit(600e3, 60.0, u0_deg=20.0)


def _orbit_run(drift):
    p0 = platform_state_at(ORBIT, 0.0)
    alt = np.linspace(50e3, 700e3, 1301)
    state = IonosphereState(alt_m=alt, ne=(1.5e6**2 * np.exp((600e3 - alt) / 200e3)) / 80.6163)

    class _L0:
        name = "fixed"
        layer = "l0"

        def state_at(self, lat, lon, epoch, alt_m):
            return state

        def describe(self):
            return {"name": "fixed"}

    facet = build("l1", "facet", geometry="orbit", ref_alt_m=600e3,
                  facets=[{"id": "o", "lat_deg": p0.lat_deg, "lon_deg": p0.lon_deg, "alt_m": 300e3,
                           "drift_mps": drift}])
    p = Pipeline({"l0": _L0(), "l1": facet,
                  "l2": BaselineInstrument(quantise=False, thermal_noise=False)}, seed=3)
    sweep = sweep_law.drift_mode(F)
    t = sweep_law.drift_mode_times(91, t0_s=-45 / 45 - 1 / 90)
    plats = [platform_state_at(ORBIT, float(x) + 1 / 90) for x in t]
    return p.run_sequence("orb", sweep=sweep, t_s=t, platforms=plats, alt_m=alt,
                          epoch="2026-01-01T12:00:00Z", truth_at=(p0.lat_deg, p0.lon_deg))


@pytest.fixture(scope="module")
def orbit_run():
    return _orbit_run(50.0)


def test_orbit_noise_free_drift_within_1_mps_with_the_effective_speed(orbit_run):
    prod = DriftProcessing(prior=PRIOR).process_sequence(orbit_run.sequence)
    p = {k: v[0] for k, v in prod.products.items()}
    assert p["geometry"] == "orbit" and p["echo_class"] == "compact_facet"
    assert p["v_eff_mps"] < p["v_mps"]
    assert p["facet_depth_m"] == pytest.approx(300e3, rel=0.01)
    assert p["drift_mps"] == pytest.approx(50.0, abs=1.0)


def test_orbit_effective_speed_moves_the_slope_ratio_not_the_drift(orbit_run):
    """Measured 2026-11-25, noise-free: the drift route one forms is the Doppler
    at the power peak plus S (t_bar - t_pk), and the last term, the only one v
    enters, vanishes for a symmetric lobe; the curvature moves the slope ratio
    (the compact-facet test) by (v_eff / v)^2, and the drift only at second
    order (0.06 to 0.28 m/s with an asymmetric envelope, docs/drift-tier.md).
    P2's note that the drift would be 4 % wrong was arithmetic on an assumed
    estimator, and is corrected."""
    a = DriftProcessing(prior=PRIOR).process_sequence(orbit_run.sequence).products
    prod = DriftProcessing(prior=PRIOR, use_effective_speed=False).process_sequence(orbit_run.sequence)
    b = prod.products
    assert float(b["drift_mps"][0]) == pytest.approx(float(a["drift_mps"][0]), abs=0.05)
    k = (float(a["v_eff_mps"][0]) / float(a["v_mps"][0])) ** 2
    assert float(b["slope_ratio"][0]) / float(a["slope_ratio"][0]) == pytest.approx(k, rel=2e-3)
    assert 0.95 < k < 0.96
    assert "FLAT_SPEED" in prod.flags["echo_flags"][0]
