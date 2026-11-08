"""A uniform sweep with the leakage on (code review F2, 2026-11-08).

``config/instrument/baseline_6u.yaml`` turns the transmit-to-receive leakage
on, and the leakage metrics then check the INS-11 high-pass rejection segment
by segment.  For the worst segment they converted its mid-time to a frequency
with :func:`tss.sweep_law.frequency_hz`, which knows only the stepped A1 law,
so any uniform sweep (``step_dwell_s = 0``) on the baseline instrument raised
``ValueError: a stepped sweep needs step_fraction or step_cap_hz``.  The rest
of the method already branched on the law; this one call did not.

These tests hold the baseline instrument to rule 3 on a uniform sweep: it
returns a record and its metrics, and a campaign on a uniform sweep runs to an
exit status.  A uniform sweep has one high-pass segment, so the frequency
reported for the INS-11 margin is the mid-sweep frequency.
"""

from __future__ import annotations

import io

import numpy as np
import pytest

import tss  # noqa: F401
from tss.config import load_campaign, load_yaml, resolve_path
from tss.interfaces import EchoSet, SweepSchedule
from tss.studies import build_instrument, run_campaign

#: A short uniform sweep, 1 to 2 MHz at 5.6 MHz/s, quick to synthesise.
_UNIFORM = SweepSchedule(
    f_start_hz=1.0e6, f_stop_hz=2.0e6, k_f_hz_per_s=5.6e6, sweep_period_s=1.0 / 5.6
)


def _no_echoes() -> EchoSet:
    e = np.empty(0, dtype=float)
    return EchoSet(
        freq_hz=e, mode=np.empty(0, dtype=object), group_delay_s=e,
        amplitude_db=e, absorption_db=e,
    )


@pytest.fixture(scope="module")
def record():
    config = load_campaign("smoke").config
    configuration = load_yaml(resolve_path("leo_500km_single", subdir="configurations"))
    instrument = build_instrument("baseline", config, configuration)
    return instrument, instrument.receive(_no_echoes(), _UNIFORM, np.random.default_rng(0))


def test_the_baseline_instrument_receives_a_uniform_sweep_with_leakage_on(record):
    instrument, signal = record
    assert instrument.leakage_enabled, "the instrument file is expected to turn the leakage on"
    assert signal.flags["leakage_modelled"]
    assert signal.flags["sweep_law"] == "uniform"
    assert np.all(np.isfinite(signal.samples))


def test_the_leakage_metrics_of_a_uniform_sweep_are_reported(record):
    _, signal = record
    m = signal.flags["l2_metrics"]
    assert np.isfinite(m["dynamic_range_db_at_digitiser"])
    assert np.isfinite(m["ins11_rejection_margin_db"])
    # One high-pass segment: the reported frequency is the mid-sweep one.
    f = m["ins11_rejection_margin_freq_hz"]
    assert _UNIFORM.f_start_hz <= f <= _UNIFORM.f_stop_hz
    assert f == pytest.approx(0.5 * (_UNIFORM.f_start_hz + _UNIFORM.f_stop_hz), rel=1e-3)
    # INS-24 is scored per step and says why it is not scored here.
    assert "stepped sweep" in m["ins24_not_scored_reason"]


def test_a_campaign_on_a_uniform_sweep_runs_to_an_exit_status(tmp_path):
    status = run_campaign(
        "smoke", out=tmp_path, overrides=["sweep.step_dwell_s=0"], stream=io.StringIO()
    )
    assert status in (0, 1)
    assert (tmp_path / "manifest.json").is_file()
