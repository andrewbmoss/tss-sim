"""WP4 Phase 2: the detector (tss.l3_signal.detect).

Exit criteria of the phase: the measured P_fa on noise-only records is within a
factor of 1.5 of nominal, and the detection probability against SNR matches
theory to 0.5 dB.
"""

from __future__ import annotations

import io
import math
from dataclasses import replace

import numpy as np
import pytest

from tss.l3_signal import detect as d
from tss.l3_signal.transform import transform

PFA = 1.0e-4


class _Block:
    """The fields of a de-chirped block that :func:`transform` reads."""

    def __init__(self, x, fs=78.4e3, k=5.6e6):
        self.samples = np.atleast_2d(x)
        self.fs_hz = fs
        self.k_f_hz_per_s = k
        self.centre_frequency_hz = 10e6
        self.sweep = None


def _noise(rng, n, sigma2=1.0):
    return np.sqrt(sigma2 / 2) * (rng.standard_normal(n) + 1j * rng.standard_normal(n))


def _spectrum(x, pad=2):
    return transform(_Block(x), window="hann", n_fft=pad * np.size(x))


def _tone(n, bin_index, amp):
    t = np.arange(n)
    return amp * np.exp(2j * np.pi * bin_index * t / n)


# -- the threshold multiplier -------------------------------------------------------


def test_os_alpha_tends_to_the_known_noise_threshold():
    a = d.os_alpha(PFA, 100_000, 50_000)
    assert a == pytest.approx(-math.log(PFA) / math.log(2), rel=0.01)


def test_os_alpha_matches_monte_carlo_for_a_short_reference():
    rng = np.random.default_rng(1)
    n, k, pfa = 24, 12, 1e-2
    a = d.os_alpha(pfa, n, k)
    ref = np.sort(rng.exponential(size=(400_000, n)), axis=1)[:, k - 1]
    cut = rng.exponential(size=400_000)
    assert np.mean(cut > a * ref) == pytest.approx(pfa, rel=0.1)


# -- false alarms ----------------------------------------------------------------


def test_measured_pfa_on_white_noise_is_within_1p5_of_nominal():
    rng = np.random.default_rng(2)
    hits = cells = 0
    for _ in range(150):
        sp = _spectrum(_noise(rng, 2048))
        r = d.cfar_detect(sp, pfa=PFA)
        hits += int(r.mask.sum())
        cells += r.mask.size
    assert cells * PFA > 30
    assert 1 / 1.5 < (hits / cells) / PFA < 1.5


def test_measured_pfa_on_full_l2_noise_records_with_leakage():
    from tss.interfaces import EchoSet
    from tss.l2_instrument.baseline import BaselineInstrument
    from tss.l3_signal.ionogram import BaselineSignalProcessing
    from tss.studies import _base_config, _sweep_from_config, instrument_config

    cfg = _base_config(())
    sweep = _sweep_from_config(cfg, io.StringIO())[0]
    inst = BaselineInstrument.from_config(instrument_config(cfg))  # leakage on, the default
    assert inst.leakage_enabled
    l3 = BaselineSignalProcessing()
    empty = EchoSet(freq_hz=np.array([]), mode=np.array([], dtype=object), group_delay_s=np.array([]),
                    amplitude_db=np.array([]), absorption_db=np.array([]))
    hits = cells = near = 0
    for seed in range(2):
        for sp in l3.spectra(inst.receive(empty, sweep, np.random.default_rng(seed))):
            r = d.cfar_detect(sp, **l3.detector_settings())
            ok = sp.range_m <= 1.05e6
            hits += int(r.mask[ok].sum())
            cells += int(ok.sum())
            near += int(r.mask[sp.range_m < 50e3].sum())
    assert 1 / 1.5 < (hits / cells) / PFA < 1.5
    # The near range is searched (no exclusion) and the leakage residual there
    # stays below the threshold.
    assert near <= 2


# -- detection probability ------------------------------------------------------------


def _pd_theory(snr, thr_over_mean):
    """Steady target in complex Gaussian noise, known noise: Marcum Q_1."""
    from scipy.stats import ncx2

    return float(ncx2.sf(2 * thr_over_mean, 2, 2 * snr))


def _snr_for(pd, thr_over_mean):
    lo, hi = 1e-3, 1e4
    for _ in range(100):
        mid = math.sqrt(lo * hi)
        if _pd_theory(mid, thr_over_mean) < pd:
            lo = mid
        else:
            hi = mid
    return 10 * math.log10(mid)


def test_detection_probability_matches_theory_to_half_a_db():
    rng = np.random.default_rng(3)
    n, b = 2048, 300
    # Cell SNR: the peak-cell power of the tone over the mean noise power per cell.
    noise_cell = np.mean([np.mean(10 ** (_spectrum(_noise(rng, n)).power_db[0][400:1600] / 10)) for _ in range(20)])
    tone_cell = 10 ** (_spectrum(_tone(n, b, 1.0)).power_db[0][2 * b] / 10)
    thr = -math.log(PFA)  # known-noise threshold over the mean
    snrs_db = np.arange(6.0, 18.01, 1.0)
    pd = []
    for s in snrs_db:
        amp = math.sqrt(10 ** (s / 10) * noise_cell / tone_cell)
        hit = 0
        for _ in range(300):
            r = d.cfar_detect(_spectrum(_tone(n, b, amp) + _noise(rng, n)), pfa=PFA)
            hit += bool(r.mask[2 * b])
        pd.append(hit / 300)
    pd = np.asarray(pd)
    for target in (0.5, 0.9):
        measured = float(np.interp(target, pd, snrs_db))
        assert abs(measured - _snr_for(target, thr)) < 0.5, (target, measured, _snr_for(target, thr))


# -- the noise reference and the veto ----------------------------------------------------


def test_the_noise_reference_ignores_returns_and_clutter():
    rng = np.random.default_rng(4)
    n = 2048
    x = _noise(rng, n)
    clean = d.noise_reference(_spectrum(x))
    # A 40 dB echo and a continuum of returns over a third of the range.
    cont = sum(_tone(n, b, 3.0 * math.exp(-(b - 300) / 200)) * np.exp(1j * rng.uniform(0, 6.3))
               for b in range(300, 700))
    dirty = _spectrum(x + _tone(n, 200, 100.0) + cont)
    ref = d.noise_reference(dirty)
    assert ref["source"] == "negative_half"
    assert abs(ref["mean_db"] - clean["mean_db"]) < 0.1
    col = 10 * np.log10(np.median(10 ** (dirty.power_db[0] / 10)) / math.log(2))
    assert col - ref["mean_db"] > 3.0  # the column median is raised by the continuum


def test_the_local_veto_removes_sidelobes():
    rng = np.random.default_rng(5)
    n = 2048
    sp = _spectrum(_tone(n, 400, 100.0) + _noise(rng, n))  # about 71 dB over the noise cell
    with_veto = d.cfar_detect(sp, pfa=PFA)
    without = d.cfar_detect(sp, pfa=PFA, local_pfa=None, sidelobe_margin_db=None)
    dist = np.abs(np.arange(sp.power_db.shape[1]) - 800)
    skirt = (dist > 4) & (dist < 80)  # beyond the Hann main lobe (+-2 bins, +-4 padded)
    assert without.mask[skirt].sum() > 5  # the sidelobes pass the global test
    assert with_veto.mask[skirt].sum() == 0
    assert with_veto.mask[dist <= 2].all()


def test_the_veto_keeps_o_beside_a_stronger_x():
    # O and X on a scalar receiver: a few decibels apart and a few tens of
    # kilometres (here 12 cells) apart in range.
    rng = np.random.default_rng(9)
    n = 2048
    noise_cell = 1.5 / n
    x_amp = math.sqrt(10 ** 3.0 * noise_cell)  # 30 dB
    o_amp = math.sqrt(10 ** 2.4 * noise_cell)  # 24 dB, 6 dB weaker
    hit_o = hit_x = 0
    for _ in range(50):
        sp = _spectrum(_tone(n, 400, x_amp) + _tone(n, 412, o_amp) * np.exp(1j * rng.uniform(0, 6.3))
                       + _noise(rng, n))
        r = d.cfar_detect(sp, pfa=PFA)
        hit_x += bool(r.mask[798:803].any())
        hit_o += bool(r.mask[822:827].any())
    assert hit_x == 50 and hit_o == 50


def test_sidelobe_envelope_of_hann():
    env, null = d.sidelobe_envelope("hann", 1024)
    over = env.size / 512
    assert null == pytest.approx(2.0, abs=0.05)
    assert 20 * np.log10(env[int(2.2 * over)]) == pytest.approx(-31.5, abs=0.3)
    assert np.all(np.diff(env) <= 0)


def test_search_limit():
    rng = np.random.default_rng(6)
    sp = _spectrum(_tone(2048, 900, 30.0) + _noise(rng, 2048))
    r_all = d.cfar_detect(sp, pfa=PFA)
    lim = float(sp.range_m[1700])
    r_lim = d.cfar_detect(sp, pfa=PFA, max_range_m=lim)
    assert r_all.mask[1795:1805].any() and not r_lim.mask[sp.range_m > lim].any()


def test_without_a_negative_half_it_falls_back_and_says_so():
    rng = np.random.default_rng(7)
    sp = replace(_spectrum(_noise(rng, 2048)), negative_power_db=None)
    assert d.cfar_detect(sp, pfa=PFA).source == "positive_half_median"


# -- the tier ----------------------------------------------------------------------


def test_the_tier_defaults_to_the_stated_pfa_detector():
    from tss.l3_signal.ionogram import BaselineSignalProcessing

    t = BaselineSignalProcessing()
    s = t.describe()
    assert t.detector == "cfar" and t.pfa == PFA and t.snr_floor == "noise"
    assert t.max_range_m == 1.05e6
    assert BaselineSignalProcessing(detector="fixed").snr_floor == "column"
    with pytest.raises(ValueError):
        BaselineSignalProcessing(detector="peak")
    assert s is not None


def test_snr_at_against_noise_and_against_the_column():
    from tss.interfaces import EchoSet
    from tss.l2_instrument.baseline import BaselineInstrument
    from tss.l3_signal.ionogram import BaselineSignalProcessing
    from tss.studies import _base_config, _sweep_from_config, instrument_config

    cfg = _base_config(())
    sweep = _sweep_from_config(cfg, io.StringIO())[0]
    icfg = instrument_config(cfg)
    icfg["leakage"]["enabled"] = False
    icfg["transmit_impairments"] = {}
    inst = BaselineInstrument.from_config(icfg)
    f = 5.0e6
    # One echo at 300 km and a continuum from 400 to 900 km in the same step.
    taus = [2 * 300e3 / 3e8] + list(2 * np.arange(400e3, 900e3, 3e3) / 3e8)
    amps = [-95.0] + [-100.0] * (len(taus) - 1)
    e = EchoSet(freq_hz=np.full(len(taus), f), mode=np.array(["O"] * len(taus), dtype=object),
                group_delay_s=np.asarray(taus), amplitude_db=np.asarray(amps),
                absorption_db=np.zeros(len(taus)),
                flags=np.array([""] + ["clutter:resonance:x"] * (len(taus) - 1), dtype=object))
    rec = inst.receive(e, sweep, np.random.default_rng(8))
    l3 = BaselineSignalProcessing()
    s_noise = float(l3.snr_at(rec, [f], [300e3])[0])
    s_col = float(l3.snr_at(rec, [f], [300e3], floor="column")[0])
    assert s_noise > s_col + 3.0  # the continuum raises the column floor, not the noise reference
