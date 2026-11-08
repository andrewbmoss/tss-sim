"""Independent thermal noise on every receive channel (code review F1, 2026-11-08).

Until this fix the extra channels of a vector receive configuration were built
from the first channel's samples after its thermal noise had been added, plus
a second draw.  Channels 2..N therefore carried twice the noise power of
channel 1 (+3 dB) and noise correlated with it (about 0.7), while the comment
beside the loop said "the same signal with an independent noise realisation".

These tests hold the receiver to that comment:

* every channel has the same noise power;
* the noise of any two channels is uncorrelated;
* channel 1 of a three-channel record is bit-identical to the record of the
  same instrument with one channel, so the fix leaves every single-channel
  result as it was.

The records are echo-free or carry one weak echo, with leakage off, so the
only thing in them that differs between channels is the receiver noise.
"""

from __future__ import annotations

import itertools

import numpy as np
import pytest

import tss  # noqa: F401
from tss.interfaces import EchoSet, SweepSchedule
from tss.l2_instrument import BaselineInstrument
from tss.l2_instrument.receiver import Receiver

#: A short uniform sweep, 1 to 2 MHz at 5.6 MHz/s: about 14 000 samples at the
#: INS-06 rate, enough for the statistics below and quick to run.
_SWEEP = SweepSchedule(
    f_start_hz=1.0e6, f_stop_hz=2.0e6, k_f_hz_per_s=5.6e6, sweep_period_s=1.0 / 5.6
)


def _instrument(channels: int) -> BaselineInstrument:
    return BaselineInstrument(receiver=Receiver(channels=channels))


def _no_echoes() -> EchoSet:
    e = np.empty(0, dtype=float)
    return EchoSet(
        freq_hz=e, mode=np.empty(0, dtype=object), group_delay_s=e,
        amplitude_db=e, absorption_db=e,
    )


def _one_echo() -> EchoSet:
    return EchoSet(
        freq_hz=np.array([1.5e6]), mode=np.array(["O"], dtype=object),
        group_delay_s=np.array([2.0e-3]), amplitude_db=np.array([-40.0]),
        absorption_db=np.array([0.0]),
    )


def test_every_receive_channel_has_the_same_noise_power():
    x = _instrument(3).receive(_no_echoes(), _SWEEP, np.random.default_rng(7)).samples
    assert x.shape[0] == 3
    p = np.mean(np.abs(x) ** 2, axis=1)
    rel_db = 10.0 * np.log10(p / p[0])
    # One-sigma scatter of a power estimate over ~14 000 complex samples is
    # about 0.04 dB; the defect this guards against was +3 dB.
    assert np.all(np.abs(rel_db) < 0.2), rel_db


def test_the_noise_of_any_two_receive_channels_is_uncorrelated():
    x = _instrument(3).receive(_no_echoes(), _SWEEP, np.random.default_rng(11)).samples
    for i, j in itertools.combinations(range(x.shape[0]), 2):
        rho = np.abs(np.vdot(x[i], x[j])) / np.sqrt(
            np.vdot(x[i], x[i]).real * np.vdot(x[j], x[j]).real
        )
        # The scatter of |rho| for independent noise is about 1/sqrt(N), 0.008
        # here; the defect gave 0.7 between channel 1 and the others.
        assert rho < 0.05, (i, j, rho)


@pytest.mark.parametrize("echoes", [_no_echoes, _one_echo], ids=["no_echo", "one_echo"])
def test_channel_one_is_the_single_channel_record(echoes):
    """Adding channels changes nothing in channel 1: same draws, same arithmetic."""
    single = _instrument(1).receive(echoes(), _SWEEP, np.random.default_rng(3)).samples
    vector = _instrument(3).receive(echoes(), _SWEEP, np.random.default_rng(3)).samples
    assert single.shape[0] == 1
    np.testing.assert_array_equal(vector[0], single[0])


def test_the_echo_is_in_every_receive_channel():
    """The extra channels still carry the signal, not noise alone."""
    sig = _instrument(3).receive(_one_echo(), _SWEEP, np.random.default_rng(5)).samples
    quiet = _instrument(3).receive(_no_echoes(), _SWEEP, np.random.default_rng(5)).samples
    gained = np.mean(np.abs(sig) ** 2, axis=1) / np.mean(np.abs(quiet) ** 2, axis=1)
    assert np.all(gained > 1.5), gained
