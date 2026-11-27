"""WP5 Phase 3 part 3: statistics that might spot a rough patch (``scripts/wp5_rough_test.py``).

Report only: nothing in route one or the drift tier uses them.
"""

from __future__ import annotations

import importlib.util
import math
from pathlib import Path

import numpy as np
import pytest


@pytest.fixture(scope="module")
def script():
    path = Path(__file__).resolve().parents[2] / "scripts" / "wp5_rough_test.py"
    spec = importlib.util.spec_from_file_location("wp5_rough_test", path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_the_statistics_of_a_model_echo_and_of_noise(script):
    """A model facet in unit noise: the residual is the noise (ratio near 1),
    small beside the echo, and nothing is left outside the lobe."""
    from tss.l3_signal.novel.route_one_joint import estimate_route_one_joint, joint_model

    kw = {"f_hz": 4.988e6, "group_path_m": 529.8e3, "mu_sat": 0.9537, "v_s_mps": 7600.0}
    t = (np.arange(209) - 104) / 45.0
    rng = np.random.default_rng(1)
    s = joint_model(t, t0_s=0.0, drift_mps=25.0, width_s=0.926, **kw)
    y = 10 * s + (rng.standard_normal(t.size) + 1j * rng.standard_normal(t.size)) / math.sqrt(2)
    fit = estimate_route_one_joint(t, y, t0_start_s=0.0, drift_start_mps=25.0, **kw)
    got = script.statistics(t, y, fit, noise_var=1.0, **kw)
    assert got["residual_over_noise"] == pytest.approx(1.0, abs=0.3)
    assert got["residual_over_echo"] < 0.1
    assert got["outside_over_noise"] == pytest.approx(1.0, abs=0.3)
    # Far from the lobe only the sidelobes and noise remain; a steady background shows there.
    assert 0.7 < got["far_over_noise"] < 4.0
    raised = script.statistics(t, y + 5.0 * np.exp(2j * np.pi * rng.uniform(size=t.size)), fit, noise_var=1.0,
                               **kw)
    assert raised["far_over_noise"] > got["far_over_noise"] + 15.0


def test_a_threshold_rejects_the_stated_share_of_the_calibration_set(script):
    calib = list(range(100))
    assert script.threshold(calib, 0.05, high_is_suspect=True) == pytest.approx(94.05)
    assert np.mean(script.rejected(calib, script.threshold(calib, 0.05, True), True)) == pytest.approx(0.05)
    low = script.threshold(calib, 0.10, high_is_suspect=False)
    assert np.mean(script.rejected(calib, low, False)) == pytest.approx(0.10)


def test_one_realisation_gives_the_echo_and_its_statistics(script):
    (row,) = script.trial(("no halo", 0), levels=(20.0,))
    assert row["formed"] and abs(row["err_joint_mps"]) < 10.0
    for name in script.STATISTICS:
        assert math.isfinite(row[name])
    (patch,) = script.trial(("r2 K0", 0), levels=(20.0,))
    if patch["formed"]:
        assert patch["residual_over_echo"] > row["residual_over_echo"]


def test_the_cases_cover_true_reflectors_and_patches(script):
    assert list(script.CASES) == ["no halo", "fading 0.5 s", "fading 2 s", "r2 K10", "r2 K3", "r2 K0"]
    assert script.CALIBRATION_SEED0 >= 1000
