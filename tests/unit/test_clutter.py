"""WP4 Phase 1: the clutter generator (tss.l1_propagation.clutter) and its L2 rendering.

The exit criterion of the phase: each term matches its closed form to 0.5 dB,
and the pipeline with clutter off reproduces the Phase 0 baseline.
"""

from __future__ import annotations

import copy
import math
from dataclasses import replace

import numpy as np
import pytest

from tss.interfaces import EchoSet
from tss.l1_propagation import clutter as cl
from tss.l1_propagation.analytic import AnalyticPropagation
from tss.registry import build

TOL_DB = 0.5
C = 299_792_458.0
EPOCH = "2026-01-01T12:00:00Z"


def _db(x):
    return 10.0 * np.log10(x)


def _only(term, **over):
    s = copy.deepcopy(cl.DEFAULTS)
    for t in ("ground_echo", "spread_f", "plasma_resonances", "interference"):
        s[t]["enabled"] = t == term
    s[term].update(over)
    return cl.ClutterGenerator(s)


def _scene(fo_hz=6.0e6, h_sat=500e3, freqs=None):
    from tss.platform.orbits import CircularOrbit, platform_state_at

    l0 = build("l0", "chapman", nmf2=fo_hz**2 / 80.6, hmf2_m=300e3, plasma_scale_height_m=100e3)
    alt = np.arange(90e3, 2.0e6 + 1, 5e3)
    platform = platform_state_at(CircularOrbit(altitude=h_sat, inclination_deg=97.4), 0.0, epoch=EPOCH)
    platform = replace(platform, lat_deg=40.0, lon_deg=0.0)
    state = l0.state_at(platform.lat_deg, platform.lon_deg, EPOCH, alt)
    f = np.asarray(freqs if freqs is not None else np.linspace(0.5e6, 20e6, 80), float)
    l1 = AnalyticPropagation() if hasattr(AnalyticPropagation, "echoes") else build("l1", "analytic")
    echoes = l1.echoes(state, platform, f, np.random.default_rng(0))
    return state, platform, f, echoes


def _clutter_rows(e, kind):
    fl = np.array([str(x) for x in e.flags], dtype=object)
    return np.array([x.startswith("clutter:" + kind) for x in fl], dtype=bool)


# -- configuration ---------------------------------------------------------------


def test_from_config_is_off_when_missing_or_disabled():
    assert cl.ClutterGenerator.from_config(None) is None
    assert cl.ClutterGenerator.from_config({}) is None
    assert cl.ClutterGenerator.from_config({"enabled": False}) is None
    g = cl.ClutterGenerator.from_config({"enabled": True, "spread_f": {"level_db": -10.0}})
    assert g is not None and g.settings["spread_f"]["level_db"] == -10.0
    assert g.settings["spread_f"]["rows_per_echo"] == cl.DEFAULTS["spread_f"]["rows_per_echo"]


def test_default_config_turns_clutter_on_and_smoke_pins_it_off():
    from tss.config import load_campaign
    from tss.studies import _base_config

    assert cl.ClutterGenerator.from_config(_base_config(()).get("clutter")) is not None
    assert cl.ClutterGenerator.from_config(load_campaign("smoke").config.get("clutter")) is None


# -- ground echo -------------------------------------------------------------------


def test_ground_echo_in_a_thin_ionosphere_is_the_free_space_closed_form():
    # foF2 0.3 MHz: every sounding frequency penetrates, and from 5 MHz up (well
    # above f_H, where the X mode's group retardation is still 2 per cent at
    # 2 MHz) the group path is the altitude to a part in a thousand and the
    # F-region loss is nil, so the echo is a mirror at the ground behind the
    # declared losses.
    state, platform, f, e0 = _scene(fo_hz=0.3e6, freqs=np.linspace(5e6, 20e6, 10))
    g = _only("ground_echo")
    e = g.augment(e0, state, platform, f, np.random.default_rng(1))
    sel = _clutter_rows(e, "ground")
    assert np.count_nonzero(sel) == np.count_nonzero(np.array([x == "penetrating" for x in e0.flags]))
    fr, tau, amp = e.freq_hz[sel], e.group_delay_s[sel], e.amplitude_db[sel]
    h = platform.alt_m
    np.testing.assert_allclose(tau, 2 * h / C, rtol=2e-3)
    s = cl.DEFAULTS["ground_echo"]
    expect = (20 * np.log10((C / fr) / (4 * math.pi * 2 * h))
              - s["d_region_a_db_mhz2"] / (fr / 1e6 + s["d_region_f_l_mhz"]) ** 2 + s["surface_reflection_db"])
    assert np.max(np.abs(amp - expect)) < TOL_DB


def test_ground_echo_only_above_penetration_and_later_than_the_vacuum_path():
    state, platform, f, e0 = _scene(fo_hz=6.0e6)
    e = _only("ground_echo").augment(e0, state, platform, f, np.random.default_rng(1))
    sel = _clutter_rows(e, "ground")
    assert np.all(e.freq_hz[sel & (e.mode == "O")] > 6.0e6 * 0.99)
    assert np.all(e.group_delay_s[sel] > 2 * platform.alt_m / C)
    assert set(e.mode[sel]) <= {"O", "X"}


# -- spread F -------------------------------------------------------------------


def test_spread_f_power_offsets_and_frequency_spread():
    state, platform, f, e0 = _scene(fo_hz=6.0e6)
    sp = cl.DEFAULTS["spread_f"]
    e = _only("spread_f").augment(e0, state, platform, f, np.random.default_rng(3))
    sel = _clutter_rows(e, "spread_f")
    n = sp["rows_per_echo"]
    flags0 = np.array([str(x) for x in e0.flags], dtype=object)
    parent = np.isfinite(e0.group_delay_s) & np.array([x not in ("evanescent", "penetrating") for x in flags0])
    # Total power beside each parent echo is level_db below it, to round-off.
    for i in np.nonzero(parent)[0][:20]:
        mine = sel & (e.freq_hz == e0.freq_hz[i]) & (e.mode == e0.mode[i])
        assert np.count_nonzero(mine) == n
        tot = _db(np.sum(10 ** (e.amplitude_db[mine] / 10)))
        assert abs(tot - (e0.amplitude_db[i] + sp["level_db"])) < TOL_DB
        off = (e.group_delay_s[mine] - e0.group_delay_s[i]) * C / 2
        assert np.all((off >= 0) & (off <= sp["extent_km"] * 1e3))
    # The truncated exponential's mean, over every row.
    m, x = sp["mean_km"], sp["extent_km"]
    mean_expect = m - x * math.exp(-x / m) / (1 - math.exp(-x / m))
    offs = []
    for i in np.nonzero(parent)[0]:
        mine = sel & (e.freq_hz == e0.freq_hz[i]) & (e.mode == e0.mode[i])
        offs.extend((e.group_delay_s[mine] - e0.group_delay_s[i]) * C / 2e3)
    assert abs(np.mean(offs) - mean_expect) < 4 * np.std(offs) / math.sqrt(len(offs))
    # Some rows above each mode's last echo (frequency spread).
    top_o = np.max(e0.freq_hz[parent & (e0.mode == "O")])
    above = sel & (e.mode == "O") & (e.freq_hz > top_o)
    assert np.any(above) and np.all(e.freq_hz[above] <= top_o * (1 + sp["beyond_penetration_fraction"]))


# -- plasma resonances -------------------------------------------------------------


def test_resonances_sit_in_the_step_that_contains_them_and_decay_as_stated():
    from tss.l1_propagation import magnetoionic
    from tss.l1_propagation.base import plasma_frequency_at

    state, platform, f, e0 = _scene(fo_hz=6.0e6, freqs=np.linspace(0.33e6, 20e6, 400))
    g = _only("plasma_resonances")
    e = g.augment(e0, state, platform, f, np.random.default_rng(0))
    r = g.settings["plasma_resonances"]
    f_n = float(plasma_frequency_at(state, platform.alt_m))
    f_h = float(np.ravel(magnetoionic.field_at(state, platform.alt_m)[0])[0])
    step = f[1] - f[0]
    for name, f_r, spec in (("fN", f_n, r["fN"]), ("fUH", math.hypot(f_n, f_h), r["fUH"]),
                            ("1fH", f_h, r["nfH"]), ("2fH", 2 * f_h, r["nfH"])):
        sel = _clutter_rows(e, "resonance:" + name)
        assert np.count_nonzero(sel) > 10, name
        assert np.all(np.abs(e.freq_hz[sel] - f_r) <= step / 2 + 1.0), name
        tau = e.group_delay_s[sel]
        expect = spec["level_db"] - 10 / math.log(10) * tau / (spec["decay_ms"] * 1e-3)
        assert np.max(np.abs(e.amplitude_db[sel] - expect)) < TOL_DB
        span = r["floor_db"] * spec["decay_ms"] * 1e-3 / (10 / math.log(10)) * C / 2
        assert np.max(tau) * C / 2 <= min(span, r["extent_km"] * 1e3) + 1.0
        np.testing.assert_allclose(np.diff(np.sort(tau)) * C / 2, r["row_spacing_km"] * 1e3, rtol=1e-9)


# -- interference ----------------------------------------------------------------


def test_carriers_only_in_broadcast_bands_above_fof2():
    state, platform, f, e0 = _scene(fo_hz=8.0e6)
    e = _only("interference").augment(e0, state, platform, f, np.random.default_rng(5))
    c = e.interference
    assert c is not None and c["freq_hz"].size > 0
    assert np.all(c["freq_hz"] > 8.0e6 * 0.999)
    in_band = [any(lo * 1e6 <= x <= hi * 1e6 for lo, hi in cl.BROADCAST_BANDS_MHZ) for x in c["freq_hz"]]
    assert all(in_band)
    assert np.all(np.mod(c["freq_hz"] - np.round(c["freq_hz"] / 5e3) * 5e3, 5e3) < 1e-3)


def test_carrier_baseband_power_and_band_limit():
    from tss.l2_instrument.interference import carrier_baseband

    fs, k = 200e3, 1e6
    t = np.arange(0, 1.0, 1 / fs)
    f_ref = 5e6 + k * t
    phi = 2 * np.pi * (5e6 * t + k * t**2 / 2)
    p_w = 1e-9
    x = carrier_baseband(phi, f_ref, t, 5.5e6, p_w, 0.3, fs)
    fb = f_ref - 5.5e6
    inside = np.abs(fb) < 0.85 * fs / 2
    assert abs(_db(np.mean(np.abs(x[inside]) ** 2)) - _db(p_w)) < 1e-6
    assert np.all(x[np.abs(fb) >= fs / 2] == 0)
    # Its instantaneous frequency is f_ref - f_c (the receiver's sign convention).
    ph = np.unwrap(np.angle(x[inside]))
    inst = np.diff(ph) * fs / (2 * np.pi)
    np.testing.assert_allclose(inst, fb[inside][:-1], atol=k / fs + 1.0)


def test_render_carriers_applies_the_antenna_gain_and_polarization_loss():
    from tss.l2_instrument.interference import render_carriers

    fs, k = 200e3, 1e6
    t = np.arange(0, 1.0, 1 / fs)
    f_ref = 5e6 + k * t
    phi = 2 * np.pi * (5e6 * t + k * t**2 / 2)
    carriers = {"freq_hz": np.array([5.5e6]), "power_dbm_iso": np.array([-70.0]), "phase_rad": np.array([0.0])}
    x, info = render_carriers(carriers, phi_ref_rad=phi, f_ref_hz=f_ref, t_s=t, fs_hz=fs,
                              gain_dbi=lambda f_mhz: -12.0, polarization_loss_db=-3.0)
    inside = np.abs(f_ref - 5.5e6) < 0.85 * fs / 2
    assert info["carriers_rendered"] == 1
    assert abs(_db(np.mean(np.abs(x[inside]) ** 2)) + 30 - (-70.0 - 12.0 - 3.0)) < TOL_DB


# -- L2 and scoring --------------------------------------------------------------


def _instrument():
    from tss.l2_instrument.baseline import BaselineInstrument
    from tss.studies import _base_config, instrument_config

    icfg = instrument_config(_base_config(()))
    icfg["leakage"]["enabled"] = False
    icfg["transmit_impairments"] = {}
    return BaselineInstrument.from_config(icfg)


def _sweep():
    import io

    from tss.studies import _base_config, _sweep_from_config

    return _sweep_from_config(_base_config(()), io.StringIO())[0]


def test_l2_synthesises_a_clutter_row_exactly_as_an_echo():
    inst, sweep = _instrument(), _sweep()
    row = dict(freq_hz=np.array([5.0e6]), mode=np.array(["O"], dtype=object),
               group_delay_s=np.array([3.4e-3]), amplitude_db=np.array([-100.0]),
               absorption_db=np.zeros(1), doppler_hz=np.zeros(1))
    a = inst.receive(EchoSet(**row, flags=np.array([""], dtype=object)), sweep, np.random.default_rng(9))
    b = inst.receive(EchoSet(**row, flags=np.array(["clutter:ground"], dtype=object)), sweep,
                     np.random.default_rng(9))
    np.testing.assert_array_equal(a.samples, b.samples)
    assert b.flags["clutter_rows_used"] == 1 and a.flags["clutter_rows_used"] == 0


def test_l2_adds_listed_carriers_and_nothing_else():
    inst, sweep = _instrument(), _sweep()
    e = EchoSet(freq_hz=np.array([5.0e6]), mode=np.array(["O"], dtype=object),
                group_delay_s=np.array([3.4e-3]), amplitude_db=np.array([-100.0]), absorption_db=np.zeros(1))
    car = {"freq_hz": np.array([12.0e6]), "power_dbm_iso": np.array([-60.0]), "phase_rad": np.array([0.0])}
    a = inst.receive(e, sweep, np.random.default_rng(4))
    b = inst.receive(replace(e, interference=car), sweep, np.random.default_rng(4))
    assert b.flags["interference"]["carriers_rendered"] == 1
    d = np.abs(b.samples[0] - a.samples[0])
    assert np.count_nonzero(d) > 0
    assert np.count_nonzero(d) < 0.05 * d.size  # only while the sweep passes 12 MHz


def test_truth_cells_never_include_clutter():
    from tss.l5_assessment.scoring import truth_echo_cells

    state, platform, f, e0 = _scene(fo_hz=6.0e6)
    gen = cl.ClutterGenerator(copy.deepcopy(cl.DEFAULTS))
    e = gen.augment(e0, state, platform, f, np.random.default_rng(2))
    assert e.freq_hz.size > e0.freq_hz.size
    a, b = truth_echo_cells(e0), truth_echo_cells(e)
    np.testing.assert_array_equal(a["freq_hz"], b["freq_hz"])
    np.testing.assert_array_equal(a["range_m"], b["range_m"])
    np.testing.assert_array_equal(a["cutoff_step"], b["cutoff_step"])


def test_augment_does_not_modify_its_input_and_is_seeded():
    state, platform, f, e0 = _scene(fo_hz=6.0e6)
    before = (e0.freq_hz.copy(), e0.amplitude_db.copy(), list(e0.flags))
    gen = cl.ClutterGenerator(copy.deepcopy(cl.DEFAULTS))
    a = gen.augment(e0, state, platform, f, np.random.default_rng(7))
    b = gen.augment(e0, state, platform, f, np.random.default_rng(7))
    np.testing.assert_array_equal(e0.freq_hz, before[0])
    np.testing.assert_array_equal(e0.amplitude_db, before[1])
    assert list(e0.flags) == before[2]
    np.testing.assert_array_equal(a.group_delay_s, b.group_delay_s)
    np.testing.assert_array_equal(a.interference["freq_hz"], b.interference["freq_hz"])


# -- climatology -----------------------------------------------------------------


def test_climatology_occurrence_and_scaling():
    assert cl._spread_f_probability(5.0, 21.0) == 0.7
    assert cl._spread_f_probability(5.0, 12.0) == 0.02
    assert cl._spread_f_probability(70.0, 12.0) == 0.5
    assert cl._spread_f_probability(40.0, 0.0) == 0.1
    assert cl.cos_solar_zenith(0.0, 180.0, EPOCH) < 0.0  # midnight
    assert cl.cos_solar_zenith(-23.0, 0.0, EPOCH) > 0.99  # noon, subsolar latitude in January
    assert cl.local_time_h(135.0, EPOCH) == 21.0
    assert cl.d_region_loss_db(4.8e6, 100.0, 1.2) == pytest.approx(100.0 / 36.0)


# -- the pipeline -----------------------------------------------------------------


def test_clutter_off_reproduces_the_phase0_baseline():
    import csv
    from pathlib import Path

    from tss.pipeline import build_tiers
    from tss.studies import _base_config
    from tss.studies import wp4_campaign as w

    ref = Path(__file__).resolve().parents[2] / "data" / "reference" / "wp4_baseline" / "wp4_campaign_runs.csv"
    rows = {(r["profile"], r["seed"]): r for r in csv.DictReader(ref.open(encoding="utf-8"))}
    cfg = _base_config(())
    shared = build_tiers({k: v for k, v in cfg["tiers"].items() if k not in ("l0", "l2")})
    got = w.run_one("smoke", 1, config=cfg, configuration="leo_500km_single", clutter=None, shared=shared)
    want = rows[("smoke", "1")]
    for k in ("sci05", "sci07_m", "ins10_db"):  # the CSV carries six decimals
        assert got[k] == pytest.approx(float(want[k]), rel=0, abs=6e-7)
    assert got["per_column"] == want["per_column"]


def test_pipeline_records_the_clutter_it_applied():
    from tss.pipeline import build_tiers
    from tss.studies import _base_config
    from tss.studies import wp4_campaign as w

    cfg = _base_config(("grid.alt_step_m=20000",))
    shared = build_tiers({k: v for k, v in cfg["tiers"].items() if k not in ("l0", "l2")})
    gen = cl.ClutterGenerator.from_config(cfg["clutter"])
    row = w.run_one("smoke", 1, config=cfg, configuration="leo_500km_single", clutter=gen, shared=shared,
                    keep=True)
    res = row["_result"]
    assert row["clutter"] == "on" and row["clutter_rows"] > 0
    assert res.flags["clutter"]["ground_echo_rows"] > 0
    assert res.signal.flags["clutter_rows_used"] == row["clutter_rows"]
