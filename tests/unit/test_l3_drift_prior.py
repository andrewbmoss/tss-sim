"""WP5 Phase 3 part 2: the drift tier's onboard prior from the ionogram chain (Option A).

f_N at the spacecraft comes from the WP4 tracker on an ionogram of the drift
scene's own L0; the profile below is an exponential of an assumed scale height.
"""

from __future__ import annotations

import io
import math
from types import SimpleNamespace

import numpy as np
import pytest

from tss.config import load_campaign
from tss.constants import plasma_frequency
from tss.interfaces import PlatformState
from tss.l3_signal.drift import (
    OnboardPrior,
    _group_integrals,
    effective_speed_ratio,
    onboard_prior_from_ionogram,
)
from tss.pipeline import Pipeline
from tss.platform.orbits import RE
from tss.registry import build

PLAT = PlatformState(alt_m=600e3, lat_deg=0.0, lon_deg=0.0)


def _ion(f_n, lowest=2.0e6, priors=True):
    freq = np.array([1.0e6, 1.5e6, lowest, 3.0e6])
    vr = np.full((4, 2), np.nan)
    vr[2:, 0] = [300e3, 320e3]
    flags = {"tracker": {"priors": {"f_n_hz": f_n, "f_n_method": "trace_start", "trace_floor_hz": 1.9e6}}} \
        if priors else {}
    return SimpleNamespace(freq_hz=freq, virtual_range_m=vr, flags=flags)


def test_a_valid_tracker_estimate_becomes_an_onboard_exponential_prior():
    prior, rep = onboard_prior_from_ionogram(_ion(1.8e6), PLAT, scale_height_m=150e3)
    assert rep["reason"] == "accepted" and prior is not None
    assert prior.source == "onboard"
    assert prior.extra["scale_height_source"] == "assumed"
    assert prior.extra["f_n_method"] == "trace_start"
    assert plasma_frequency(np.interp(600e3, prior.alt_m, prior.ne)) == pytest.approx(1.8e6, rel=1e-6)
    assert plasma_frequency(np.interp(450e3, prior.alt_m, prior.ne)) == pytest.approx(
        1.8e6 * math.exp(0.5), rel=1e-3)


@pytest.mark.parametrize("ion, why", [
    (_ion(2.1e6), "not below the lowest"),
    (_ion(float("nan")), "no plasma frequency"),
    (_ion(None), "no plasma frequency"),
    (_ion(1.8e6, priors=False), "no tracker priors"),
])
def test_an_invalid_estimate_gives_no_prior(ion, why):
    prior, rep = onboard_prior_from_ionogram(ion, PLAT)
    assert prior is None
    assert why in rep["reason"]


def test_the_depth_is_found_when_the_group_path_reaches_past_the_prior_reflection_level():
    """A 100 km scale height at 5 MHz: the group path of a facet 170 km down
    is longer than the depth to reflection, so the first bracket is evanescent.
    The bracket used to be halved until it propagated, which cut it below the
    answer (108 km)."""
    prior = OnboardPrior.exponential(f_n_sat_hz=1.82e6, scale_height_m=100e3, ref_alt_m=600e3)
    r_p = RE + 600e3
    g, _ = _group_integrals(prior, 5e6, 600e3, 170e3, r_p)
    assert not math.isfinite(_group_integrals(prior, 5e6, 600e3, g, r_p)[0])
    _, depth = effective_speed_ratio(prior, 5e6, 600e3, g, r_p)
    assert depth == pytest.approx(170e3, abs=10.0)


def test_the_prior_ionogram_of_a_chapman_scene_gives_f_n_within_a_step():
    """The tracker on the A1 sweep, clutter off: within the step's half-width
    (the trace-start rule takes the midpoint below the first column)."""
    cfg = load_campaign("smoke").config
    from tss.studies import build_instrument, campaign_grids, sweep_from_config

    sweep, _ = sweep_from_config(cfg, io.StringIO())
    alt, freq = campaign_grids(cfg)
    conf = {"name": "leo_600", "orbit": {"altitude_m": 6.0e5, "inclination_deg": 97.4},
            "antenna_count": 1, "antenna_length_m": 5.0}
    l0 = build("l0", "chapman")
    pipe = Pipeline({"l0": l0}, seed=7)
    tiers = {"l1": build("l1", "analytic"), "l2": build_instrument("baseline", cfg, conf), "l3": build("l3", "baseline")}
    res = pipe.run_prior_ionogram("p", tiers=tiers, platform=PLAT, sweep=sweep, freq_hz=freq, alt_m=alt,
                                  epoch="2026-01-01T12:00:00Z")
    prior, rep = onboard_prior_from_ionogram(res.ionogram, PLAT)
    truth = l0.state_at(0.0, 0.0, "2026-01-01T12:00:00Z", alt)
    f_true = float(plasma_frequency(np.interp(600e3, truth.alt_m, truth.ne)))
    assert prior is not None, rep["reason"]
    step = float(np.min(np.abs(freq - f_true))) + 0.5 * float(np.diff(freq)[np.searchsorted(freq, f_true)])
    assert abs(rep["f_n_hz"] - f_true) < step
    assert res.flags["prior_pass"]["clutter"] is False


@pytest.mark.slow
def test_on_an_l0_orbit_scene_the_onboard_prior_moves_the_depth_not_the_drift():
    """Chapman L0 (f_N 1.82 MHz at 600 km, topside scale height 100 km), a 3-D
    plate at 430 km drifting 50 m/s, 5 MHz, noise-free.  The truth prior is the
    L0 profile; the onboard prior is the tracker's f_N with an assumed 200 km
    scale height.  Measured 2026-11-25: drift 49.982 (truth) and 49.992
    (onboard); depth 169.9 and 189.8 km.  With clutter on in the prior pass
    the tracker takes a resonance 2.4 % high and the drift moves 0.4 %."""
    from tss import sweep_law
    from tss.l1_propagation.clutter import ClutterGenerator
    from tss.l2_instrument.baseline import BaselineInstrument
    from tss.l3_signal.drift import DriftProcessing
    from tss.platform.orbits import CircularOrbit, platform_state_at
    from tss.studies import build_instrument, campaign_grids, sweep_from_config

    orbit = CircularOrbit(altitude=600e3, inclination_deg=97.4, name="o")
    p0 = platform_state_at(orbit, 0.0)
    facet = build("l1", "facet", geometry="orbit", profile="l0", ref_alt_m=600e3,
                  facets=[{"id": "o", "lat_deg": p0.lat_deg, "lon_deg": p0.lon_deg, "alt_m": 430e3,
                           "drift_mps": 50.0}])
    pipe = Pipeline({"l0": build("l0", "chapman"), "l1": facet,
                     "l2": BaselineInstrument(quantise=False, thermal_noise=False)}, seed=3)
    t = sweep_law.drift_mode_times(91, t0_s=-45 / 45 - 1 / 90)
    run = pipe.run_sequence("orb", sweep=sweep_law.drift_mode(5e6), t_s=t,
                            platforms=[platform_state_at(orbit, float(x) + 1 / 90) for x in t],
                            alt_m=np.linspace(50e3, 900e3, 1701), epoch="2026-01-01T12:00:00Z",
                            truth_at=(p0.lat_deg, p0.lon_deg))
    cfg = load_campaign("smoke").config
    sweep, _ = sweep_from_config(cfg, io.StringIO())
    alt, freq = campaign_grids(cfg)
    conf = {"name": "leo_600", "orbit": {"altitude_m": 6.0e5, "inclination_deg": 97.4},
            "antenna_count": 1, "antenna_length_m": 5.0}
    tiers = {"l1": build("l1", "analytic"), "l2": build_instrument("baseline", cfg, conf),
             "l3": build("l3", "baseline")}

    def drift(prior):
        p = DriftProcessing(prior=prior).process_sequence(run.sequence)
        return {k: v[0] for k, v in p.products.items()}, p.flags

    truth, _ = drift(OnboardPrior(alt_m=run.truth.alt_m, ne=run.truth.ne, source="truth"))
    ion = pipe.run_prior_ionogram("orb", tiers=tiers, platform=p0, sweep=sweep, freq_hz=freq, alt_m=alt,
                                  epoch="2026-01-01T12:00:00Z").ionogram
    onboard, flags = drift(onboard_prior_from_ionogram(ion, p0)[0])
    assert truth["drift_mps"] == pytest.approx(50.0, abs=0.1)
    assert truth["facet_depth_m"] == pytest.approx(170e3, abs=1e3)
    assert onboard["drift_mps"] == pytest.approx(truth["drift_mps"], abs=0.05)
    assert onboard["facet_depth_m"] - truth["facet_depth_m"] == pytest.approx(20e3, abs=3e3)
    assert flags["prior_source"] == "onboard" and flags["prior"]["scale_height_source"] == "assumed"

    # The range gate from the same ionogram: above the facet's range, and the
    # gated drift is the ungated one (measured: gate 379 km, range 215.65 km).
    from tss.l3_signal.drift import range_gate_from_ionogram

    gate, rep = range_gate_from_ionogram(ion, 5e6)
    assert gate is not None, rep["reason"]
    assert gate > onboard["range_m"] + 50e3
    gated = DriftProcessing(prior=onboard_prior_from_ionogram(ion, p0)[0]).process_sequence(
        run.sequence, max_range_m=gate)
    assert gated.flags["range_gate"] and gated.flags["detection"]["range_bins"] < 400
    assert gated.products["drift_mps"][0] == pytest.approx(onboard["drift_mps"], abs=0.02)

    noisy = {**tiers, "clutter": ClutterGenerator.from_config(load_campaign("baseline").config["clutter"])}
    ion_c = pipe.run_prior_ionogram("orb", tiers=noisy, platform=p0, sweep=sweep, freq_hz=freq, alt_m=alt,
                                    epoch="2026-01-01T12:00:00Z").ionogram
    prior_c, rep_c = onboard_prior_from_ionogram(ion_c, p0)
    assert prior_c is not None, rep_c["reason"]
    with_clutter, _ = drift(prior_c)
    assert with_clutter["drift_mps"] == pytest.approx(truth["drift_mps"], rel=0.01)
