"""Physics and regression checks on the migrated link budget.

Ported from the author's ``LinkBudget/test_link_budget.py``.  Every test
function and every assertion is carried across.  The only change is the import
path: the module under test is now ``tss.l2_instrument.linkbudget`` instead of a
top-level ``tss_link_budget``, and it is bound to the same local name ``lb`` so
that the bodies are unchanged.

These are behavioural checks on the model, not just numerical regressions: each
one asserts something that must be true of the physics, so a future change that
breaks the physics fails here rather than silently producing a plausible number.

They are also the reason the link budget module is kept whole rather than split
along layer lines.  See ``docs/architecture.md``, section "Known seams".
"""

from __future__ import annotations

import math
import sys
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
SRC = REPO / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

from tss.l2_instrument import linkbudget as lb  # noqa: E402


def _cfg(**kw):
    base = dict(sweep_rate=5.6e6, range_resolution=1000.0, h_sat=600e3,
                true_range=300e3, tx_power=1.0, polarization_loss_db=-6.0,
                rx_noise_figure_db=3.0, sky_model="franklin",
                ionosphere=lb.IonosphereProfile(nmf2=1.0e12, hmf2=300e3,
                                                plasma_scale_height=100e3))
    base.update(kw)
    return lb.SounderConfig(**base)


# --- plasma physics ---------------------------------------------------------
def test_plasma_frequency_round_trip():
    for ne in (1e9, 1e10, 1e11, 1e12):
        assert math.isclose(lb.density_from_plasma_frequency(
            lb.plasma_frequency(ne)), ne, rel_tol=1e-12)


def test_plasma_frequency_matches_textbook_constant():
    # fN [Hz] ~ 8.98 sqrt(Ne [m^-3]) is the form used throughout the literature.
    assert abs(lb.PLASMA_CONST - 8.98) < 0.01


# --- sky noise --------------------------------------------------------------
def test_franklin_anchors():
    assert math.isclose(float(lb.sky_temperature(0.6)), 2.0e7, rel_tol=1e-9)
    assert math.isclose(float(lb.sky_temperature(12.0)), 4.0e5, rel_tol=1e-9)


def test_sky_temperature_falls_with_frequency():
    f = np.array([1.0, 3.0, 10.0, 30.0])
    for model in ("franklin", "itu"):
        t = lb.sky_temperature(f, model)
        assert np.all(np.diff(t) < 0)


def test_franklin_and_itu_agree_within_a_few_db_in_midband():
    # Independent references; they should not disagree wildly where both apply.
    for f in (3.0, 5.0, 8.0, 12.0):
        d = abs(10 * math.log10(float(lb.sky_temperature(f, "franklin"))
                                / float(lb.sky_temperature(f, "itu"))))
        assert d < 5.0, f"{f} MHz: {d:.1f} dB apart"
    # They are closest in mid-band (0.4 dB at 5 MHz) and diverge toward 12 MHz
    # (4.2 dB), where Franklin's in-situ values run above the ground-based model.


def test_legacy_references_are_optimistic():
    # The whole point of the correction.
    t_true = float(lb.sky_temperature(5.0, "franklin"))
    assert t_true / float(lb.sky_temperature(5.0, "legacy300")) > 1000
    assert t_true / float(lb.sky_temperature(5.0, "legacy_sni")) > 100


# --- profile and geometry ---------------------------------------------------
def test_reflection_height_decreases_with_frequency():
    cfg = _cfg()
    hs = [lb.reflection_height(f * 1e6, cfg.ionosphere, cfg.h_sat)
          for f in (3.0, 4.0, 5.0, 6.0, 7.0, 8.0)]
    assert all(isinstance(h, float) for h in hs)
    assert all(a > b for a, b in zip(hs, hs[1:], strict=False))


def test_band_edges_are_flagged():
    cfg = _cfg()
    f_lo = cfg.ionosphere.local_plasma_frequency(cfg.h_sat)
    f_hi = cfg.ionosphere.fof2()
    assert lb.link_budget(f_lo * 0.9, cfg).status == "evanescent"
    assert lb.link_budget(f_hi * 1.1, cfg).status == "penetrates"
    assert lb.link_budget(0.5 * (f_lo + f_hi), cfg).status == "ok"
    assert math.isnan(lb.link_budget(f_lo * 0.9, cfg).snr_db)


# --- absorption -------------------------------------------------------------
def test_absorption_is_converged_and_smooth():
    cfg = _cfg()
    for f in (3.0, 5.0, 7.0):
        vals = [lb.absorption_db(f * 1e6, cfg.ionosphere, cfg.h_sat, n_steps=n)
                for n in (500, 3000, 20000)]
        assert max(vals) - min(vals) < 1e-3, f"{f} MHz not converged: {vals}"


def test_absorption_rises_toward_penetration():
    cfg = _cfg()
    f_hi = cfg.ionosphere.fof2() / 1e6
    a = [lb.absorption_db(f * 1e6, cfg.ionosphere, cfg.h_sat)
         for f in np.linspace(3.0, f_hi * 0.99, 12)]
    assert all(b > a_ for a_, b in zip(a, a[1:], strict=False))


def test_legacy_absorption_diverges_where_the_new_model_does_not():
    # This is the term that produced the apparent low-frequency shortfall.
    cfg = _cfg()
    f = 3.0e6                       # inside the usable band for this profile
    new = lb.absorption_db(f, cfg.ionosphere, cfg.h_sat)
    old = lb.absorption_db(f, cfg.ionosphere, cfg.h_sat, legacy_model=True)
    assert new < 1.0 and old > 10.0


# --- the external-noise-limited result --------------------------------------
def test_receive_efficiency_cancels_while_the_sky_dominates():
    """
    The central correction: with an external-noise-limited receiver, dropping
    the radiation efficiency of BOTH antennas by 3 dB costs 3 dB, not 6 dB,
    because the receive antenna couples less sky noise in exactly the same
    proportion as it couples less signal.
    """
    cfg = _cfg()
    ant = lb.Antenna(eff_ref_db=(-14.76, -8.76))     # 3 dB worse everywhere
    cfg_lo = _cfg(tx_antenna=ant, rx_antenna=ant)
    d = lb.link_budget(5e6, cfg).snr_db - lb.link_budget(5e6, cfg_lo).snr_db
    assert 2.9 < d < 3.1, d


def test_receiver_noise_takes_over_at_very_low_efficiency():
    """Below about -30 dB efficiency the receiver, not the sky, sets the floor,
    and the cancellation above stops applying."""
    poor = lb.Antenna(eff_ref_db=(-45.0, -45.0))
    cfg = _cfg(tx_antenna=poor, rx_antenna=poor)
    r = lb.link_budget(20e6, cfg)
    t_rx = (lb.undb(cfg.rx_noise_figure_db) - 1.0) * lb.T0_REF
    assert r.system_temperature_k < 3.0 * t_rx


# --- FMCW processing --------------------------------------------------------
def test_finer_range_cell_improves_snr_and_costs_dwell():
    """
    A result specific to de-chirped FMCW and worth stating: the range cell IS
    the analysis bin, so a FINER range cell is a NARROWER bin and therefore
    LESS noise.  Range resolution and sensitivity are not in tension here, as
    they are for a pulse radar.  What a finer cell costs is dwell time
    (T = 1/B), and therefore along-track resolution.
    """
    fine, coarse = _cfg(range_resolution=1000.0), _cfg(range_resolution=2000.0)
    assert math.isclose(coarse.noise_bandwidth / fine.noise_bandwidth, 2.0,
                        rel_tol=1e-12)
    d = lb.link_budget(5e6, fine).snr_db - lb.link_budget(5e6, coarse).snr_db
    assert math.isclose(d, 10 * math.log10(2.0), abs_tol=1e-9), d
    assert math.isclose(fine.dwell_time / coarse.dwell_time, 2.0, rel_tol=1e-12)


def test_coherent_integration_gain():
    a = _cfg()
    b = _cfg(coherent_sweeps=64)
    d = lb.link_budget(5e6, b).snr_db - lb.link_budget(5e6, a).snr_db
    assert math.isclose(d, 10 * math.log10(64), abs_tol=1e-9)


def test_snr_scales_with_transmit_power():
    a = lb.link_budget(5e6, _cfg(tx_power=0.1)).snr_db
    b = lb.link_budget(5e6, _cfg(tx_power=1.0)).snr_db
    assert math.isclose(b - a, 10.0, abs_tol=1e-9)


def test_required_power_inverts_the_budget():
    cfg = _cfg()
    for target in (5.0, 10.0, 20.0):
        p = lb.required_tx_power(5e6, cfg, target)
        got = lb.link_budget(5e6, _cfg(tx_power=p)).snr_db
        assert math.isclose(got, target, abs_tol=1e-9)


# --- legacy regression ------------------------------------------------------
def test_legacy_port_matches_the_original_matlab():
    """Hand-evaluated from snr_ionosonde.m for Kf=5.6 MHz/s, R=300 km,
    Rmax=500 km, f=8 MHz, Pt=1 W."""
    snr, pwr, n = lb.legacy_snr_ionosonde(5.6e6, 300e3, 500e3, 8e6, 1.0)
    assert abs(snr - 13.58) < 0.02
    assert abs(pwr - (-115.88)) < 0.02
    assert abs(n - (-129.46)) < 0.02


def test_legacy_model_is_pessimistic_overall():
    """The swept-band noise error dominates the too-cold-sky error, so the old
    model understates SNR despite two optimistic assumptions."""
    cfg = _cfg()
    new = lb.link_budget(5e6, cfg).snr_db
    old, _, _ = lb.legacy_snr_ionosonde(cfg.sweep_rate, cfg.true_range, 500e3,
                                        5e6, cfg.tx_power)
    assert old < new - 15.0


# --- migration checks, added here rather than in the original ----------------
def test_constants_come_from_the_package_and_not_from_a_local_copy():
    """The migration replaced local constant definitions with imports.

    Not in the author's original file.  It is here because the one way the
    migration could have gone wrong silently is a constant that was retyped
    rather than imported, and a retyped constant that differs in the last digit
    would move every number in this module by an amount none of the tests above
    would catch.
    """
    from tss import constants, units

    assert lb.C_LIGHT is constants.C_LIGHT
    assert lb.K_BOLTZ is constants.K_BOLTZ
    assert lb.PLASMA_CONST is constants.PLASMA_CONST
    assert lb.T0_REF == constants.T0_K
    assert lb.db is units.db
    assert lb.undb is units.undb


def test_antenna_noise_temperature_agrees_with_the_antenna_module():
    """The two statements of ``T_a = eta T_sky + (1-eta) T_phys`` must agree.

    Not in the author's original file.  ``tss.l2_instrument.antenna`` restates
    the relation because that module is where the argument for it lives, and two
    statements of one relation is exactly the thing that drifts.
    """
    from tss.l2_instrument.antenna import antenna_noise_temperature_k

    ant = lb.Antenna()
    for f_mhz in (0.5, 1.0, 5.0, 20.0):
        t_sky = float(lb.sky_temperature(f_mhz))
        here = float(ant.noise_temperature(f_mhz, t_sky))
        there = float(
            antenna_noise_temperature_k(
                float(ant.efficiency(f_mhz)), t_sky, ant.physical_temperature
            )
        )
        assert math.isclose(here, there, rel_tol=1e-12)


if __name__ == "__main__":
    fns = [(k, v) for k, v in sorted(globals().items()) if k.startswith("test_")]
    bad = 0
    for name, fn in fns:
        try:
            fn()
            print(f"  PASS  {name}")
        except AssertionError as e:
            bad += 1
            print(f"  FAIL  {name}: {e}")
    print("-" * 60)
    print(f"{len(fns) - bad}/{len(fns)} passed" if bad else f"ALL {len(fns)} TESTS PASS")
