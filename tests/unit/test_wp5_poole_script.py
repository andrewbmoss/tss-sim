"""WP5 Phase 3 part 3: Poole's comparator on the drift sequences (``scripts/wp5_poole.py``)."""

from __future__ import annotations

import importlib.util
import math
from pathlib import Path

import pytest


@pytest.fixture(scope="module")
def script():
    path = Path(__file__).resolve().parents[2] / "scripts" / "wp5_poole.py"
    spec = importlib.util.spec_from_file_location("wp5_poole", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


@pytest.fixture(scope="module")
def clean(script):
    return script.trial(("nadir", 0), estimators=("pinned",), levels=(None,))


def test_noise_free_pooles_pairs_follow_the_line_of_sight_rate(clean):
    """Between successive sweeps Poole's phase velocity is the combined
    line-of-sight rate L1 put into the echo, the platform's part included."""
    (row,) = clean
    assert row["detected"]
    assert row["poole_pair_rms_mps"] < 1.5
    # The same observable as the tier's phase steps.
    assert row["poole_minus_tier_rms_mps"] < 1.0


def test_noise_free_the_lobe_mean_of_pooles_velocity_is_near_the_drift(clean):
    """Over a whole symmetric lobe the platform's part averages out."""
    (row,) = clean
    assert row["err_poole_mean"] == pytest.approx(0.0, abs=1.0)
    assert abs(row["err_poole_peak"]) < 8.0
    assert row["err_pinned"] == pytest.approx(0.0, abs=1.0)


def test_an_extended_plate_gets_a_poole_velocity_and_no_route_one_drift(script):
    (row,) = script.trial(("extended", 0), estimators=("pinned",), levels=(None,))
    assert math.isfinite(row["err_poole_mean"]) and math.isnan(row["err_pinned"])


def test_the_summary_keeps_each_quantity_apart(script, clean):
    out = script.summarise(clean)[("nadir", None)]
    assert out["n"] == 1 and out["detected"] == 1
    assert out["err_poole_mean"]["n"] == 1 and out["err_pinned"]["n"] == 1
