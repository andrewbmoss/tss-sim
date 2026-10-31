"""WP4 Phase 3: the dispersion and step-transition treatments.

The treatments live in tss.l3_signal.ionogram.BaselineSignalProcessing
(``dispersion``, ``transition``) and the narrow-band step law in
tss.sweep_law; the comparison is the study tss.studies.wp4_dispersion.
"""

from __future__ import annotations

import io

import numpy as np
import pytest

from tss import sweep_law
from tss.interfaces import EchoSet
from tss.l3_signal.detect import cfar_detect
from tss.l3_signal.ionogram import BaselineSignalProcessing

C = 299_792_458.0
PFA = 1.0e-4


def _cfg():
    from tss.studies import _base_config

    return _base_config(())


def _sweep(**over):
    from tss.studies import _sweep_from_config

    cfg = _cfg()
    cfg["sweep"] = {**cfg["sweep"], **over}
    return _sweep_from_config(cfg, io.StringIO())[0]


def _instrument(leakage=True, dispersion=True):
    from tss.l2_instrument.baseline import BaselineInstrument
    from tss.studies import instrument_config

    icfg = instrument_config(_cfg())
    icfg["leakage"]["enabled"] = leakage
    if not leakage:
        icfg["transmit_impairments"] = {}
    icfg["in_dwell_dispersion"] = dispersion
    return BaselineInstrument.from_config(icfg)


def _empty():
    return EchoSet(freq_hz=np.array([]), mode=np.array([], dtype=object), group_delay_s=np.array([]),
                   amplitude_db=np.array([]), absorption_db=np.array([]))


# -- option 3: the narrow-band step law -----------------------------------------------


def test_narrow_band_step_law():
    a1 = _sweep()
    nb = _sweep(narrow_lo_hz=2.0e6, narrow_hi_hz=8.0e6, narrow_cap_hz=50.0e3)
    assert sweep_law.step_count(a1) == 159
    e = sweep_law.step_edges_hz(nb)
    w = np.diff(e)
    inside = (e[:-1] >= 2.0e6) & (e[:-1] < 8.0e6)
    assert np.all(w[inside] <= 50.0e3 + 1e-6)
    assert np.all(w[~inside] <= 150.0e3 + 1e-6)
    assert sweep_law.step_count(nb) == 237
    assert sweep_law.active_duration_s(nb) == pytest.approx(237 * 0.0268)
    # The phase stays continuous across the narrow band's edges.
    t = np.linspace(0, sweep_law.active_duration_s(nb), 200001)
    dphi = np.diff(sweep_law.phase_rad(nb, t))
    f = sweep_law.frequency_hz(nb, t[:-1])
    assert np.max(np.abs(dphi - 2 * np.pi * f * (t[1] - t[0]))) < 0.05


# -- the detector under the treatments ---------------------------------------------------


@pytest.mark.parametrize("kw", [{"dispersion": "blocks", "dispersion_blocks": 2},
                                {"dispersion": "blocks", "dispersion_blocks": 4},
                                {"dispersion": "rate_search"}, {"dispersion": "none", "transition": "gate"},
                                {"dispersion": "none", "transition": "matched"}])
def test_measured_pfa_under_each_treatment(kw):
    sweep = _sweep()
    inst = _instrument(leakage=False)
    l3 = BaselineSignalProcessing(**kw)
    hits = cells = 0
    # The rate search's look count is calibrated (about 7 effective looks), so
    # it is checked on six records of its own seeds, not the calibration's.
    seeds = range(200, 206) if kw.get("dispersion") == "rate_search" else range(2)
    for seed in seeds:
        for sp in l3.spectra(inst.receive(_empty(), sweep, np.random.default_rng(seed))):
            r = cfar_detect(sp, **l3.detector_settings())
            ok = sp.range_m <= 1.05e6
            hits += int(r.mask[ok].sum())
            cells += int(ok.sum())
    assert 1 / 1.5 < (hits / cells) / PFA < 1.5, hits / cells


# -- option 2: the rate search -----------------------------------------------------------


def test_rate_search_compresses_a_sliding_beat_and_returns_its_rate():
    from tss.l3_signal.dechirp import Dechirped

    fs, n = 78.4e3, 2100
    dur = n / fs
    t = (np.arange(n) - 0.5 * (n - 1)) / fs
    j0 = 12
    gamma = j0 / dur**2
    f_b = 20e3
    x = np.exp(2j * np.pi * (f_b * t + 0.5 * gamma * t**2))
    sw = _sweep()
    blk = Dechirped(samples=x[None, :], fs_hz=fs, t_s=t + 1.0, sweep=sw, oriented=True,
                    centre_frequency_hz=10e6, k_f_hz_per_s=5.6e6)
    plain = BaselineSignalProcessing(dispersion="none")._transform_block(blk, n)
    rs = BaselineSignalProcessing(dispersion="rate_search")._transform_block(blk, n)
    assert np.max(rs.power_db) > np.max(plain.power_db) + 3.0
    assert np.max(rs.power_db) == pytest.approx(0.0, abs=0.5)  # a unit tone reads 0 dB
    assert rs.rate_hz_per_s2[0] == pytest.approx(gamma, rel=0.01)
    k = int(np.argmax(rs.power_db[0]))
    assert rs.beat_hz[k] == pytest.approx(f_b, abs=2 * fs / n)


def test_rate_search_and_narrow_steps_recover_the_dispersion_loss_on_a_dispersive_topside():
    from dataclasses import replace

    from tss.config import load_yaml, resolve_path
    from tss.studies import _platform_for
    from tss.studies import wp4_campaign as w
    from tss.studies import wp4_dispersion as d

    cfg = _cfg()
    conf = load_yaml(resolve_path("leo_500km_single", subdir="configurations"))
    lat, lt, *_ = w.PROFILES["mid_day_min"]
    pl = replace(_platform_for(conf, w.EPOCH), lat_deg=lat, lon_deg=w._lon_for_local_time(lt))
    tr = {k: d.TREATMENTS[k] for k in ("D0", "D2w", "D3")}
    rows = d.run_scene("mid_day_min", w._l0_for("mid_day_min"), pl, w.EPOCH, config=cfg, seeds=[1],
                       treatments=tr)
    agg = d.aggregate(rows)
    assert agg["D0"]["all"]["loss_db_median"] > 2.0  # the effect is there
    assert agg["D2w"]["all"]["loss_db_median"] < 0.5
    assert agg["D3"]["all"]["loss_db_median"] < 0.5
    assert agg["D2w"]["all"]["width_km_median"] < agg["D0"]["all"]["width_km_median"]
    assert agg["D2w"]["all"]["dpdf_err_median_abs"] < 0.1


# -- option 1: sub-blocks --------------------------------------------------------------


def test_sub_blocks_coarsen_the_cell_by_their_count():
    sweep = _sweep()
    rec = _instrument(leakage=False).receive(_empty(), sweep, np.random.default_rng(0))
    a = BaselineSignalProcessing(dispersion="none").spectra(rec)
    b = BaselineSignalProcessing(dispersion="blocks", dispersion_blocks=4).spectra(rec)
    assert b[50].n_incoherent == 4
    assert b[50].range_resolution_m == pytest.approx(4 * a[50].range_resolution_m, rel=0.01)


# -- the step transition ----------------------------------------------------------------


def test_gate_drops_the_first_2_rmax_over_c():
    sweep = _sweep()
    rec = _instrument(leakage=False).receive(_empty(), sweep, np.random.default_rng(0))
    a = BaselineSignalProcessing(dispersion="none").spectra(rec)
    g = BaselineSignalProcessing(dispersion="none", transition="gate").spectra(rec)
    n_a = a[10].range_resolution_m
    n_g = g[10].range_resolution_m
    gate = 2 * 1.05e6 / C
    assert n_g / n_a == pytest.approx(0.0268 / (0.0268 - gate), rel=0.02)


def test_matched_filter_removes_the_transition_skirt_of_a_long_range_echo():
    # A strong non-dispersive echo at 900 km in a low step, where the rate
    # changes from step to step: in the first tau the beat is the previous
    # step's, which puts a skirt at R K_(i-1) / K_i, about 5 per cent short.
    # Strong enough that the skirt is detected, which is when the hybrid
    # filter looks at it (it re-evaluates detections only).
    sweep = _sweep()
    f = sweep_law.step_centres_hz(sweep)
    j = int(np.argmin(np.abs(f - 1.2e6)))
    e = EchoSet(freq_hz=np.array([f[j]]), mode=np.array(["O"], dtype=object),
                group_delay_s=np.array([2 * 900e3 / C]), amplitude_db=np.array([-50.0]),
                absorption_db=np.zeros(1))
    rec = _instrument(leakage=False, dispersion=False).receive(e, sweep, np.random.default_rng(2))
    acc = BaselineSignalProcessing(dispersion="none").spectra(rec)[j]
    mf = BaselineSignalProcessing(dispersion="none", transition="matched").spectra(rec)[j]
    r = acc.range_m
    skirt = (r > 0.93 * 900e3) & (r < 0.97 * 900e3)
    main = np.abs(r - 900e3) < 2 * acc.range_resolution_m
    assert np.max(mf.power_db[0][skirt]) < np.max(acc.power_db[0][skirt]) - 3.0
    assert np.max(mf.power_db[0][main]) == pytest.approx(np.max(acc.power_db[0][main]), abs=1.0)


def test_the_rate_search_is_the_default_after_gate_g2():
    t = BaselineSignalProcessing()
    assert t.dispersion == "rate_search" and t.rate_search_max_cells == 64 and t.transition == "accept"
    assert t.rate_search_effective_looks == pytest.approx(0.2 * 67)


def test_unknown_treatments_are_refused():
    with pytest.raises(ValueError):
        BaselineSignalProcessing(dispersion="magic")
    with pytest.raises(ValueError):
        BaselineSignalProcessing(transition="skip")


# -- the ISIS-1 scene -----------------------------------------------------------------


def test_measured_topside_keeps_the_data_and_adds_a_peak():
    import datetime as dt
    from types import SimpleNamespace

    from tss.studies.wp4_dispersion import MeasuredTopside

    prof = SimpleNamespace(alt_km=np.array([600.0, 550.0, 500.0, 450.0, 400.0, 360.0]),
                           ne_m3=np.array([6e10, 8e10, 1.2e11, 1.9e11, 3.0e11, 4.2e11]),
                           extra={"hmf2_iri_km": 300.0, "fof2_iri_mhz": 7.0}, epoch=dt.datetime(1969, 2, 2),
                           lat_deg=10.0, lon_deg=-73.0, sat_alt_km=600.0)
    l0 = MeasuredTopside(prof)
    alt = np.arange(90e3, 2.0e6 + 1, 5e3)
    st = l0.state_at(10.0, -73.0, "1969-02-02T03:00:00Z", alt)
    for h, n in zip(prof.alt_km, prof.ne_m3, strict=True):
        assert np.interp(h * 1e3, st.alt_m, st.ne) == pytest.approx(n, rel=0.02)
    assert st.alt_m[np.argmax(st.ne)] == pytest.approx(300e3, abs=5e3)
    assert np.max(st.ne) >= 1.15 * 4.2e11 * 0.99
