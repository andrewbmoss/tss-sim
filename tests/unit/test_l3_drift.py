"""WP5 Phase 3: the drift tier (``l3/drift``), route one at signal level.

Noise-free, through L1 (facet) and L2 (baseline, quantisation and thermal
noise off).  The pass criterion of D5 option A at noise-free level is drift
within 1 m/s.  Nothing here runs with noise yet (docs/drift-tier.md,
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
