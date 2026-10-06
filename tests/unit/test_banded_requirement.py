"""Requirements stated per frequency band (INS-05, Rev G amended 2026-09-30)."""

from __future__ import annotations

import pytest

from tss.requirements import Requirement, load_requirements

BANDS = (
    {"f_max_hz": 5.0e5, "threshold": 1.5e4, "goal": 1.0e4},
    {"f_min_hz": 5.0e5, "f_max_hz": 7.0e5, "threshold": 1.0e4, "goal": 1.0e4},
    {"f_min_hz": 7.0e5, "threshold": 1.0e4, "goal": 5.0e3},
)


def _req() -> Requirement:
    return Requirement(id="INS-05", title="Virtual range resolution", sense="max",
                       unit="m", threshold=1.0e4, goal=5.0e3, source="derived",
                       derivation="test", bands=BANDS)


def test_each_band_is_scored_against_its_own_limits():
    rec = _req().evaluate({"freq_hz": [3.4e5, 6.0e5, 5.0e6], "value": [1.36e4, 9.0e3, 1.5e3]})
    verdicts = {b["band"]: b["verdict"] for b in rec["bands"]}
    assert verdicts == {"0-0.5 MHz": "threshold", "0.5-0.7 MHz": "goal", "0.7-inf MHz": "goal"}
    assert rec["verdict"] == "threshold"
    assert rec["achieved"] == pytest.approx(1.36e4)
    assert rec["threshold"] == 1.5e4 and rec["band"] == "0-0.5 MHz"


def test_a_value_that_passes_low_fails_high():
    """13.6 km passes below 0.5 MHz but is a failure above it."""
    rec = _req().evaluate({"freq_hz": [3.4e5, 8.0e5], "value": [1.2e4, 1.2e4]})
    assert rec["verdict"] == "fail"
    assert rec["band"] == "0.7-inf MHz"


def test_band_edges_are_lower_inclusive():
    rec = _req().evaluate({"freq_hz": [5.0e5], "value": [1.2e4]})
    assert rec["verdict"] == "fail"


def test_scalar_is_scored_against_the_strictest_limits():
    rec = _req().evaluate(1.2e4)
    assert rec["verdict"] == "fail"
    assert _req().evaluate(4.0e3)["verdict"] == "goal"


def test_empty_bands_are_not_evaluated():
    rec = Requirement(id="X", title="x", sense="max", threshold=1.0, source="s",
                      bands=({"f_min_hz": 1e7, "threshold": 1.0},)).evaluate(
        {"freq_hz": [1e6], "value": [0.5]})
    assert rec["verdict"] == "not_evaluated"


def test_shipped_ins05_is_banded_and_valid():
    reqs = load_requirements()
    assert reqs.validate() == []
    ins05 = reqs["INS-05"]
    assert len(ins05.bands) == 3
    assert ins05.evaluate({"freq_hz": [3.3e5], "value": [1.364e4]})["verdict"] == "threshold"
