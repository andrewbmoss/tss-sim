"""The pipeline's optional hooks are declared and recorded (code review F5, 2026-11-08).

Four things the pipeline uses are not part of the six tier protocols: L3's
``snr_at`` (INS-10 at the truth echoes, with a ``floor`` keyword for signal to
noise plus clutter), L5's ``l2_metrics`` keyword (INS-11, INS-23, INS-24,
INS-29 to INS-31), L1's ``group_index_used`` (the gate on SCI-05 and SCI-07)
and the clutter generator.  Until this fix they were found by ``getattr`` and
signature probes, so a renamed method turned requirements into "not evaluated"
without an error.  They are now declared in :mod:`tss.interfaces`, and every run
records which of them it found.
"""

from __future__ import annotations

import io
import json

import numpy as np

import tss  # noqa: F401
from tss.config import load_yaml, resolve_path
from tss.interfaces import (
    ClutterAugmenter,
    MeasuresTruthEchoSNR,
    ReportsGroupIndex,
    SweepSchedule,
    accepts_keyword,
    hook_support,
)
from tss.l1_propagation.clutter import ClutterGenerator
from tss.pipeline import Pipeline, build_tiers
from tss.studies import _platform_for, run_campaign

_UNIFORM = SweepSchedule(
    f_start_hz=1.0e6, f_stop_hz=2.0e6, k_f_hz_per_s=5.6e6, sweep_period_s=1.0 / 5.6
)
_EPOCH = "2026-01-01T12:00:00Z"
_DEFAULT = {"l0": "chapman", "l1": "analytic", "l2": "baseline", "l3": "baseline",
            "l4": "lamination_magnetoionic", "l5": "standard"}


def _run(tiers):
    platform = _platform_for(
        load_yaml(resolve_path("leo_500km_single", subdir="configurations")), _EPOCH
    )
    return Pipeline(tiers, seed=3).run_configuration(
        "leo_500km_single", platform=platform, sweep=_UNIFORM,
        freq_hz=np.linspace(1.05e6, 1.95e6, 16), alt_m=np.arange(9.0e4, 6.0e5, 2.0e4),
        epoch=_EPOCH,
    )


def test_the_default_tiers_provide_every_hook():
    tiers = build_tiers(_DEFAULT)
    assert isinstance(tiers["l1"], ReportsGroupIndex)
    assert isinstance(tiers["l3"], MeasuresTruthEchoSNR)
    assert accepts_keyword(tiers["l3"].snr_at, "floor")
    assert accepts_keyword(tiers["l5"].score, "l2_metrics")
    clutter = ClutterGenerator({"enabled": True})
    assert isinstance(clutter, ClutterAugmenter)
    assert hook_support(tiers, clutter) == {
        "l1_group_index": True, "l3_snr_at": True, "l3_snr_floor": True,
        "l5_l2_metrics": True, "clutter": True,
    }
    assert hook_support(tiers)["clutter"] is None


def test_a_run_records_the_hooks_it_used():
    result = _run(build_tiers(_DEFAULT))
    assert result.flags["hooks"] == {
        "l1_group_index": True, "l3_snr_at": True, "l3_snr_floor": True,
        "l5_l2_metrics": True, "clutter": None,
    }


class _L3WithoutSnr:
    """An L3 tier that forms ionograms but has no ``snr_at``."""

    def __init__(self, inner) -> None:
        self._inner = inner

    def form_ionogram(self, signal, platform=None):
        return self._inner.form_ionogram(signal, platform)

    def describe(self):
        return {"name": "no_snr"}


def test_an_l3_without_snr_at_is_recorded_and_leaves_ins10_not_evaluated():
    tiers = build_tiers(_DEFAULT)
    tiers["l3"] = _L3WithoutSnr(tiers["l3"])
    assert not isinstance(tiers["l3"], MeasuresTruthEchoSNR)
    result = _run(tiers)
    assert result.flags["hooks"]["l3_snr_at"] is False
    assert result.flags["hooks"]["l3_snr_floor"] is False
    assert result.truth_echo_snr is None
    assert result.assessment.requirements["INS-10"]["verdict"] == "not_evaluated"


def test_the_campaign_record_carries_the_hooks(tmp_path):
    assert run_campaign("smoke", out=tmp_path, stream=io.StringIO()) in (0, 1)
    record = json.loads((tmp_path / "campaign.json").read_text(encoding="utf-8"))
    hooks = record["per_configuration"]["leo_500km_single"]["hooks"]
    assert hooks["l3_snr_at"] is True and hooks["l5_l2_metrics"] is True
    assert hooks["clutter"] is None  # smoke is pinned clutter-off
