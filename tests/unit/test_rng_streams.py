"""Separate random streams for the L2 sub-models and for L1 (code review F3, 2026-11-08).

Until this fix one generator served every draw of a record in turn: the
transmit impairments, the leakage, the RF canceller and then the thermal
noise, all in L2, after L1 had been handed the same generator.  Switching a
sub-model on or off changed how many draws came before the thermal noise, so
the noise realisation changed with it, and a paired on/off comparison was not
a comparison of the sub-model alone.

Now the impairments, the leakage and the canceller each draw from a child
stream spawned from the generator L2 is given (which does not advance it), the
thermal noise draws from that generator itself, and the pipeline gives L1 a
stream of its own.  These tests hold the instrument and the pipeline to that.
"""

from __future__ import annotations

import numpy as np

import tss  # noqa: F401
from tss.config import load_campaign, load_yaml, resolve_path
from tss.interfaces import EchoSet, SweepSchedule
from tss.l2_instrument import BaselineInstrument
from tss.l2_instrument.impairments import TransmitImpairments
from tss.pipeline import Pipeline, build_tiers
from tss.studies import _platform_for, build_instrument, instrument_config

#: A short uniform sweep, 1 to 2 MHz at 5.6 MHz/s, quick to synthesise.
_UNIFORM = SweepSchedule(
    f_start_hz=1.0e6, f_stop_hz=2.0e6, k_f_hz_per_s=5.6e6, sweep_period_s=1.0 / 5.6
)
_EPOCH = "2026-01-01T12:00:00Z"


def _no_echoes() -> EchoSet:
    e = np.empty(0, dtype=float)
    return EchoSet(
        freq_hz=e, mode=np.empty(0, dtype=object), group_delay_s=e,
        amplitude_db=e, absorption_db=e,
    )


def _smoke_config():
    return load_campaign("smoke").config


def test_switching_the_transmit_impairments_leaves_the_noise_alone():
    """With no echo and no leakage the impairments have nothing to act on, so
    the record must not change when they are switched on."""
    imp = TransmitImpairments.from_config(
        instrument_config(_smoke_config()).get("transmit_impairments")
    )
    assert imp.needs_realization, "the baseline impairments are expected to draw"
    on = BaselineInstrument(impairments=imp).receive(_no_echoes(), _UNIFORM, np.random.default_rng(4))
    off = BaselineInstrument().receive(_no_echoes(), _UNIFORM, np.random.default_rng(4))
    np.testing.assert_array_equal(on.samples, off.samples)


def test_the_same_seed_gives_the_same_record_with_every_sub_model_on():
    inst = build_instrument(
        "baseline", _smoke_config(),
        load_yaml(resolve_path("leo_500km_single", subdir="configurations")),
    )
    assert inst.leakage_enabled and inst.impairments.needs_realization
    a = inst.receive(_no_echoes(), _UNIFORM, np.random.default_rng(9))
    b = inst.receive(_no_echoes(), _UNIFORM, np.random.default_rng(9))
    np.testing.assert_array_equal(a.samples, b.samples)
    assert a.flags["l2_substreams"] is True


class _DrawingL1:
    """An L1 tier that draws from its generator before delegating."""

    group_index_used = None

    def __init__(self, inner) -> None:
        self.inner = inner

    def echoes(self, state, platform, freq_hz, rng):
        rng.normal(size=1000)
        return self.inner.echoes(state, platform, freq_hz, rng)

    def describe(self):
        return {"name": "drawing"}


def _l2_record(l1_wraps: bool) -> np.ndarray:
    tiers = build_tiers({"l0": "chapman", "l1": "analytic", "l2": "baseline",
                         "l3": "baseline", "l4": "lamination_magnetoionic", "l5": "standard"})
    if l1_wraps:
        tiers["l1"] = _DrawingL1(tiers["l1"])
    platform = _platform_for(
        load_yaml(resolve_path("leo_500km_single", subdir="configurations")), _EPOCH
    )
    result = Pipeline(tiers, seed=3).run_configuration(
        "leo_500km_single", platform=platform, sweep=_UNIFORM,
        freq_hz=np.linspace(1.05e6, 1.95e6, 16), alt_m=np.arange(9.0e4, 6.0e5, 2.0e4),
        epoch=_EPOCH,
    )
    return result.signal.samples


def test_an_l1_that_draws_does_not_change_the_l2_record():
    np.testing.assert_array_equal(_l2_record(True), _l2_record(False))


def test_l1_and_l2_draw_from_different_streams():
    p = Pipeline({}, seed=1)
    a = p._l1_rng("leo_500km_single").normal(size=8)
    b = p._rng("leo_500km_single").normal(size=8)
    assert not np.array_equal(a, b)
    # Each is stable for its seed and configuration.
    np.testing.assert_array_equal(a, Pipeline({}, seed=1)._l1_rng("leo_500km_single").normal(size=8))

