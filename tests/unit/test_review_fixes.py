"""The fixes of the 2026-10-17 code review: sample rate, instrument file, INS-10.

The velocity frame has its own tests in ``test_platform.py`` and the contract
suite.  These cover the other three:

* the L2 digitiser rate is an instrument setting (INS-06), not 2.5 times the
  largest beat among the simulated echoes;
* campaigns build L2 from the instrument file, with the configuration's
  antenna and receive overrides, and the matching-network loss is counted once,
  in the antenna efficiency;
* INS-10 is scored on the SNR measured at every truth echo, which can fail,
  rather than on the minimum over detected cells, which could not.
"""

from __future__ import annotations

import dataclasses
import io
import json

import numpy as np
import pytest

import tss  # noqa: F401
from tss.config import load_campaign, load_yaml, resolve_path
from tss.constants import C_LIGHT
from tss.interfaces import EchoSet
from tss.l2_instrument import BaselineInstrument
from tss.l2_instrument.baseline import DEFAULT_DESIGN_RANGE_M
from tss.l5_assessment import check_requirements, compute_metrics
from tss.l5_assessment.scoring import truth_echo_cells, truth_echo_metrics
from tss.pipeline import Pipeline, build_tiers
from tss.registry import build
from tss.studies import (
    _grids,
    _platform_for,
    _sweep_from_config,
    build_instrument,
    instrument_config,
    run_campaign,
)

EPOCH = "2026-01-01T12:00:00Z"


def _smoke():
    config = load_campaign("smoke").config
    sweep, _ = _sweep_from_config(config, io.StringIO())
    alt_m, freq_hz = _grids(config)
    return config, sweep, alt_m, freq_hz


def _configuration(name: str):
    return load_yaml(resolve_path(name, subdir="configurations"))


def _echo(freq_hz, range_m, mode="O"):
    f = np.atleast_1d(np.asarray(freq_hz, dtype=float))
    r = np.atleast_1d(np.asarray(range_m, dtype=float))
    return EchoSet(
        freq_hz=f,
        mode=np.array([mode] * f.size, dtype=object),
        group_delay_s=2.0 * r / C_LIGHT,
        amplitude_db=np.full(f.size, -110.0),
        absorption_db=np.zeros(f.size),
        flags=np.array([""] * f.size, dtype=object),
    )


# ---------------------------------------------------------------------------
# The sample rate
# ---------------------------------------------------------------------------


def test_the_sample_rate_is_the_ins06_rate_whatever_the_echoes():
    _, sweep, _, _ = _smoke()
    l2 = BaselineInstrument()
    near = l2.receive(_echo(5.0e6, 3.0e5), sweep, np.random.default_rng(0))
    far = l2.receive(_echo(5.0e6, 9.0e5), sweep, np.random.default_rng(0))
    expected = 4.0 * sweep.k_f_hz_per_s * DEFAULT_DESIGN_RANGE_M / C_LIGHT
    assert near.fs_hz == far.fs_hz
    assert near.fs_hz == pytest.approx(expected, rel=1e-8)
    assert near.fs_hz == pytest.approx(78.4e3, rel=1e-3)          # INS-06 at 5.6 MHz/s
    assert near.flags["unambiguous_range_m"] >= DEFAULT_DESIGN_RANGE_M


def test_the_unambiguous_range_reaches_ins04_through_the_whole_chain():
    config, sweep, alt_m, freq_hz = _smoke()
    platform = _platform_for(_configuration("leo_500km_single"), EPOCH)
    tiers = build_tiers({k: v for k, v in config["tiers"].items() if k != "l2"})
    tiers["l2"] = build_instrument("baseline", config, _configuration("leo_500km_single"))
    result = Pipeline(tiers, seed=1).run_configuration(
        "leo_500km_single", platform=platform, sweep=sweep, freq_hz=freq_hz,
        alt_m=alt_m, epoch=EPOCH,
    )
    ins04 = result.assessment.requirements["INS-04"]
    assert ins04["achieved"] >= 1.05e6
    assert ins04["verdict"] in ("threshold", "goal")


def test_a_configured_sample_rate_wins():
    _, sweep, _, _ = _smoke()
    signal = BaselineInstrument(sample_rate_hz=1.2e5).receive(
        _echo(5.0e6, 3.0e5), sweep, np.random.default_rng(0)
    )
    assert signal.fs_hz == 1.2e5
    assert signal.flags["sample_rate_rule"] == "configured sample_rate_hz"


def test_an_echo_beyond_the_unambiguous_range_is_counted():
    _, sweep, _, _ = _smoke()
    signal = BaselineInstrument().receive(
        _echo([10.0e6, 12.0e6], [6.0e5, 1.4e6]), sweep, np.random.default_rng(0)
    )
    assert signal.flags["rows_beyond_unambiguous_range"] == 1


@pytest.mark.parametrize("bad", [{"sample_rate_hz": 0.0}, {"design_range_m": -1.0}])
def test_a_nonsensical_rate_is_refused(bad):
    with pytest.raises(ValueError):
        BaselineInstrument(**bad)


# ---------------------------------------------------------------------------
# The instrument file
# ---------------------------------------------------------------------------


def test_the_matching_loss_is_counted_once_in_the_antenna_efficiency():
    config, _, _, _ = _smoke()
    l2 = build_instrument("baseline", config, _configuration("leo_500km_single"))
    record = l2.describe()
    assert record["tx_chain_efficiency_db"] == 0.0
    assert record["stated_efficiency_db"] is None


def test_the_antenna_efficiency_follows_the_sounding_frequency():
    """Decision of 2026-10-17: the baseline flies the physical model, matching
    coil included, at every frequency, not the flat -20 dB it once stated."""
    from tss.l2_instrument import ShortDipole

    config, _, _, _ = _smoke()
    l2 = build_instrument("baseline", config, _configuration("leo_500km_single"))
    assert l2.describe()["efficiency_source"] == "ShortDipole physical model"
    f_mhz = np.array([0.33, 1.0, 3.0, 10.0, 20.0])
    eta = l2.efficiency_db(f_mhz)
    assert np.allclose(eta, ShortDipole(length_m=5.0).efficiency_db(f_mhz))
    assert np.all(np.diff(eta) > 0.0)
    assert eta[0] < -40.0 and eta[-1] > -1.0


def test_a_configuration_that_changes_the_antenna_uses_the_physical_model():
    config, _, _, _ = _smoke()
    cfg = instrument_config(config, _configuration("leo_500km_small_antenna"))
    assert cfg["antenna"]["length_m"] == 1.0
    assert "radiation_efficiency_db" not in cfg["antenna"]
    assert any("1 m" in note for note in cfg["configuration_overrides"])
    small = build_instrument("baseline", config, _configuration("leo_500km_small_antenna"))
    record = small.describe()
    assert record["efficiency_source"] == "ShortDipole physical model"
    assert record["antenna"]["length_m"] == 1.0
    assert record["configuration_notes"]
    full = build_instrument("baseline", config, _configuration("leo_500km_single"))
    assert small.efficiency_db(1.0) < full.efficiency_db(1.0) - 10.0


def test_a_stated_flat_efficiency_is_dropped_when_the_length_changes(tmp_path):
    """A flat efficiency stated for one antenna does not carry to another."""
    import yaml

    config, _, _, _ = _smoke()
    cfg = instrument_config(config)
    cfg["antenna"]["radiation_efficiency_db"] = -20.0
    path = tmp_path / "flat.yaml"
    path.write_text(yaml.safe_dump(cfg), encoding="utf-8")
    out = instrument_config(
        {**config, "instrument": str(path)}, _configuration("leo_500km_small_antenna")
    )
    assert "radiation_efficiency_db" not in out["antenna"]
    assert any("physical model" in note for note in out["configuration_overrides"])


def test_the_vector_configuration_has_three_receive_channels():
    config, _, _, _ = _smoke()
    l2 = build_instrument("baseline", config, _configuration("leo_500km_vector"))
    assert l2.describe()["receiver"]["channels"] == 3


def test_a_campaign_records_the_instrument_it_flew(tmp_path):
    assert run_campaign("smoke", out=tmp_path, stream=io.StringIO()) in (0, 1)
    manifest = json.loads((tmp_path / "manifest.json").read_text(encoding="utf-8"))
    l2 = manifest["tiers"]["l2_by_configuration"]["leo_500km_single"]
    assert l2["configuration"] == "baseline_6u/leo_500km_single"
    assert l2["tx_chain_efficiency_db"] == 0.0


# ---------------------------------------------------------------------------
# INS-10 at every truth echo
# ---------------------------------------------------------------------------


def test_truth_echo_cells_keep_only_echoes_and_map_delay_to_range():
    echoes = dataclasses.replace(
        _echo([1.0e6, 2.0e6, 3.0e6], [2.0e5, 3.0e5, 4.0e5]),
        flags=np.array(["evanescent", "", "penetrating"], dtype=object),
    )
    cells = truth_echo_cells(echoes)
    assert cells["freq_hz"].tolist() == [2.0e6]
    assert cells["range_m"] == pytest.approx([3.0e5])
    assert cells["mode"] == ["O"]


def test_snr_is_measured_at_an_echo_and_not_where_there_is_none():
    _, sweep, _, _ = _smoke()
    l3 = build("l3", "baseline")
    echoes = _echo([6.0e6], [4.0e5])
    signal = BaselineInstrument(stated_efficiency_db=0.0).receive(
        echoes, sweep, np.random.default_rng(3)
    )
    snr = l3.snr_at(signal, [6.0e6, 6.0e6, 6.0e6], [4.0e5, 7.0e5, 2.0e6])
    assert snr[0] > 20.0                       # the echo
    assert snr[1] < 10.0                       # an empty cell: noise only
    assert np.isnan(snr[2])                    # beyond the unambiguous range


def test_ins10_fails_when_a_truth_echo_is_below_10_db():
    record = {"freq_hz": np.array([2.0e6, 4.0e6]), "range_m": np.array([3e5, 4e5]),
              "mode": ["O", "X"], "snr_db": np.array([25.0, 6.0])}
    metrics = truth_echo_metrics(record)
    assert metrics["snr_db_min_at_truth_echoes"] == 6.0
    assert metrics["snr_db_min_at_truth_echoes_x"] == 6.0
    assert metrics["snr_db_min_at_truth_echoes_freq_hz"] == 4.0e6
    verdict = check_requirements({"snr_db_min_at_truth_echoes": 6.0})["INS-10"]
    assert verdict["verdict"] == "fail"


def test_an_unmeasured_echo_is_counted_and_not_scored():
    metrics = truth_echo_metrics({"snr_db": np.array([np.nan, 15.0]), "mode": ["O", "O"],
                                  "freq_hz": np.array([1e6, 2e6])})
    assert metrics["n_truth_echoes"] == 2
    assert metrics["n_truth_echoes_unmeasured"] == 1
    assert metrics["snr_db_min_at_truth_echoes"] == 15.0


def test_ins10_is_not_evaluated_without_the_measurement():
    from tss.interfaces import DensityProfile, Ionogram, IonosphereState

    alt = np.linspace(2.0e5, 5.0e5, 40)
    truth = IonosphereState(alt_m=alt, ne=1.0e11 * np.exp(-(alt - 2.0e5) / 1.0e5))
    ion = Ionogram(freq_hz=np.array([1.0e6, 2.0e6]), virtual_range_m=np.array([[1e5], [2e5]]),
                   snr_db=np.array([[30.0], [30.0]]), modes=("unresolved",))
    metrics = compute_metrics(DensityProfile(alt_m=alt, ne=truth.ne), truth, ion)
    assert metrics["snr_db_min"] == 30.0      # the old metric, still a diagnostic
    record = check_requirements(metrics)["INS-10"]
    assert record["verdict"] == "not_evaluated"
    assert "snr_db_min_at_truth_echoes" in record["reason"]


def test_the_pipeline_scores_ins10_at_the_truth_echoes():
    config, sweep, alt_m, freq_hz = _smoke()
    platform = _platform_for(_configuration("leo_500km_single"), EPOCH)
    tiers = build_tiers({k: v for k, v in config["tiers"].items() if k != "l2"})
    tiers["l2"] = build_instrument("baseline", config, _configuration("leo_500km_single"))
    result = Pipeline(tiers, seed=1).run_configuration(
        "leo_500km_single", platform=platform, sweep=sweep, freq_hz=freq_hz,
        alt_m=alt_m, epoch=EPOCH,
    )
    errors = result.assessment.errors
    assert errors["n_truth_echoes"] > errors["n_detected"]   # both modes, detected or not
    ins10 = result.assessment.requirements["INS-10"]
    assert ins10["achieved"] == errors["snr_db_min_at_truth_echoes"]
    # And SCI-09 is now scored, from the platform the pipeline passes to L5.
    sci09 = result.assessment.requirements["SCI-09"]
    assert sci09["achieved"] == pytest.approx(4.26 * 7.69e3, rel=0.02)


def test_ins10_reports_each_modes_last_step_below_cutoff_separately():
    """G5c (2026-10-29): the cusp echo is reported, not scored."""
    from tss.interfaces import EchoSet

    f = np.array([3e6, 4e6, 5e6, 3e6, 4e6, 5e6, 6e6, 6e6])
    mode = np.array(["O", "O", "O", "X", "X", "X", "O", "X"], dtype=object)
    flags = np.array(["", "", "", "", "", "", "penetrating", ""], dtype=object)
    delay = np.full(f.size, 4e-3)
    echoes = EchoSet(freq_hz=f, mode=mode, group_delay_s=delay, amplitude_db=np.zeros(f.size),
                     absorption_db=np.zeros(f.size), flags=flags, source="test")
    cells = truth_echo_cells(echoes)
    # O penetrates at 6 MHz, so its 5 MHz echo is the cutoff step; X does not
    # penetrate inside the set, so none of its echoes is.
    assert list(cells["cutoff_step"]) == [False, False, True, False, False, False, False]
    snr = np.array([20.0, 18.0, 4.0, 22.0, 21.0, 19.0, 15.0])
    m = truth_echo_metrics({**cells, "snr_db": snr})
    assert m["snr_db_min_at_truth_echoes"] == 15.0
    assert m["snr_db_min_at_cutoff_steps"] == 4.0
    assert m["snr_db_min_at_truth_echoes_incl_cutoff"] == 4.0
    assert m["n_cutoff_step_echoes"] == 1
    assert m["snr_db_min_at_truth_echoes_o"] == 18.0
