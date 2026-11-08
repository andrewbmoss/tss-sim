"""WP4 Phase 4: the trace tracker (tss.l3_signal.tracker) and its use in the L3 tier."""

from __future__ import annotations

import copy
from dataclasses import replace

import numpy as np
import pytest

from tss.l3_signal import tracker as tr


def _run(profile, clutter, seed=1, l3=None):
    from tss.config import config_dir, load_yaml
    from tss.pipeline import build_tiers
    from tss.studies import wp4_campaign as w

    cfg = load_yaml(config_dir() / "default.yaml")
    if l3 is not None:
        cfg["tiers"]["l3"] = l3
    c = copy.deepcopy(cfg)
    shared = build_tiers({k: v for k, v in c["tiers"].items() if k not in ("l0", "l2")})
    return w.run_one(profile, seed, config=c, configuration="leo_500km_single", clutter=w._clutter_for(c, clutter),
                     shared=shared, keep=True)


def test_settings_refuse_unknown_policies():
    with pytest.raises(ValueError):
        tr.TrackerSettings(blend_policy="average")
    with pytest.raises(ValueError):
        tr.TrackerSettings(cusp_policy="ignore")


def test_viterbi_follows_a_weak_rising_trace_and_not_noise():
    rng = np.random.default_rng(3)
    s = tr.TrackerSettings()
    n_c, n_r = 60, 300
    z = rng.exponential(size=(n_c, n_r))
    truth = (10 + 2 * np.arange(n_c)).astype(int)
    z[np.arange(n_c), truth] += 10 ** 0.8 * rng.exponential(size=n_c)  # 8 dB fading: under the 9.7 dB cell threshold
    img = 10 * np.log10(z)
    start = np.zeros((n_c, n_r))
    path, total = tr.viterbi_path(tr._llr(img, s), np.ones_like(img, dtype=bool), start, s.grid_m, s)
    hit = np.abs(path - truth) <= 1
    # a fading 8 dB echo: columns where it fades under the bias are skipped as gaps
    assert np.mean(hit[path >= 0]) > 0.9 and np.count_nonzero(path >= 0) >= 30
    noise = 10 * np.log10(rng.exponential(size=(n_c, n_r)))
    _, t_noise = tr.viterbi_path(tr._llr(noise, s), np.ones_like(noise, dtype=bool), start, s.grid_m, s)
    assert t_noise < s.min_track_score


def test_cusp_fit_recovers_the_critical_frequency():
    fc = 10.0e6
    f = np.linspace(7.0e6, 9.8e6, 12)
    r = 120e3 + 60e3 * tr._cusp_shape(f / fc)
    fit = tr.cusp_fit(f, r)
    assert fit["fo_hz"] == pytest.approx(fc, rel=0.01)
    # a trailing point that flattens (an X sample, a ground cell) is left out
    fit2 = tr.cusp_fit(np.append(f, 9.95e6), np.append(r, r[-1] + 5e3))
    assert fit2["n_dropped"] == 1 and fit2["fo_hz"] == pytest.approx(fc, rel=0.01)


def test_x_trace_predicted_from_the_true_o_trace_matches_the_true_x_trace():
    from tss.l4_retrieval.lamination import geomagnetic_field_model
    from tss.l5_assessment.scoring import truth_echo_cells
    from tss.studies import wp4_campaign as w

    row = _run("mid_day_min", "off")
    res = row["_result"]
    lat, lt, *_ = w.PROFILES["mid_day_min"]
    from tss.config import load_yaml, resolve_path
    from tss.studies import _platform_for
    conf = load_yaml(resolve_path("leo_500km_single", subdir="configurations"))
    pl = replace(_platform_for(conf, w.EPOCH), lat_deg=float(lat), lon_deg=w._lon_for_local_time(lt))
    alt = pl.alt_m
    fld, _ = geomagnetic_field_model(pl.lat_deg, pl.lon_deg, pl.epoch or w.EPOCH, alt + 2e4)
    cells = truth_echo_cells(res.echoes)
    f = np.asarray(cells["freq_hz"])
    r = np.asarray(cells["range_m"])
    m = np.array([str(v) for v in cells["mode"]])
    cut = np.asarray(cells["cutoff_step"], dtype=bool)
    o, x = (m == "O") & ~cut, (m == "X") & ~cut
    st = res.truth
    f_n = float(np.sqrt(80.6 * np.interp(alt, st.alt_m, st.ne)))
    pred = tr.predict_x_range_m(f[o], r[o], f_n, alt, fld, f[x])
    ok = np.isfinite(pred) & (f[x] < 0.9 * np.max(f[o]))
    assert np.count_nonzero(ok) >= 10
    assert np.median(np.abs(pred[ok] - r[x][ok])) < 3.0e3


@pytest.mark.parametrize("profile, clutter, min_complete, max_label_err",
                         [("mid_day_min", "off", 0.95, 0.02), ("mid_day_max", "on", 0.8, 0.1)])
def test_tracker_labels_o_and_x_on_mid_latitude_noon(profile, clutter, min_complete, max_label_err):
    from tss.studies import wp4_campaign as w

    row = _run(profile, clutter)
    ion = row["_result"].ionogram
    assert tuple(ion.modes) == ("O", "X")
    assert ion.flags["extractor"] == "tracker" and ion.flags["mode_labelling"] == "model_self_consistent"
    assert ion.flags["mode_discrimination"] is False  # INS-08 needs the vector receive (WP11)
    # The reason says the labels are a model, not that there are none (code review F7).
    reason = ion.flags["mode_discrimination_reason"]
    assert "labelled by model self-consistency" in reason and "not measured" in reason
    m = w.trace_measures(ion, row["_result"].echoes)
    assert m["completeness_o"] >= min_complete
    assert m["labelling_error"] <= max_label_err
    assert m["n_spurious"] <= 2


def test_the_ground_echo_no_longer_takes_the_trace():
    # P1: with clutter on, the WP3-era extractor followed the ground echo
    # (O completeness 0.00).  The tracker masks it and reports foF2 from it.
    row = _run("mid_day_max", "on")
    info = row["_result"].ionogram.flags["tracker"]
    fo = info["ground"]["fo_from_ground_hz"]
    assert fo == pytest.approx(row["fof2_mhz"] * 1e6, rel=0.05)
    assert row["completeness_o"] >= 0.7


def test_onboard_plasma_frequency_from_the_resonances():
    row = _run("mid_night_min", "on")
    pri = row["_result"].ionogram.flags["tracker"]["priors"]
    st = row["_result"].truth
    from tss.studies import wp4_campaign as w  # noqa: F401
    alt = 5.0e5
    f_n = float(np.sqrt(80.6 * np.interp(alt, st.alt_m, st.ne)))
    assert pri["f_n_method"].startswith("resonance")
    assert pri["f_n_hz"] == pytest.approx(f_n, rel=0.05)


def test_legacy_extractor_is_kept():
    from tss.l3_signal.ionogram import BaselineSignalProcessing

    with pytest.raises(ValueError):
        BaselineSignalProcessing(extractor="peak")
    assert BaselineSignalProcessing(extractor="legacy").extractor == "legacy"
