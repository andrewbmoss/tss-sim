"""Layer: L3
WP: WP5
Claims: C4
Status: stub
Requirements: INS-11, INS-23, INS-24

Adaptive cancellation of transmit-to-receive leakage.  Not implemented.  WP5.

Why this is the claim the concept now turns on
------------------------------------------------
The waveform purity argument of :mod:`tss.l3_signal.dechirp` shows that range
correlation removes most of what an oscillator specification would otherwise
have to deliver.  A shared reference between the transmitter and the de-chirp
mixer suppresses oscillator phase noise at offset ``f_m`` by
``4 sin^2(pi f_m tau)``, which for a 300 km echo is 18 dB at 10 Hz offset and
38 dB at 1 Hz.  The close-in noise, which is the expensive part of an
oscillator, is exactly the part that is removed.

The leakage path gets none of that relief.  Its two-way delay is essentially
zero, so its suppression factor is essentially zero as well, and it carries the
oscillator's full noise into the receiver.  It is also tens of decibels
stronger than any echo, because the transmitter is running continuously while
the receiver listens, which is what continuous-wave operation means.  So the
hard problem did not go away when range correlation was credited; it moved from
the waveform to the isolation, and this module is where it now sits.

What the requirements ask
--------------------------
INS-23 asks for 40 dB of transmit-to-receive isolation at threshold and 60 dB
at goal, which is what the antennas, the physical separation and the front end
deliver before any processing.  INS-24 asks that the de-chirped leakage
residual sit at least 10 dB below the weakest echo INS-10 admits, 20 dB at
goal.  INS-11 asks for 60 dB of instantaneous receiver dynamic range at
threshold and 90 dB at goal, which is what the receiver must have in order that
the leakage and the weakest echo can be present at the same time without the
front end compressing.

The gap between what the hardware gives and what the processing must supply is
the subject of this module.  Note the order in which they apply: the isolation
of INS-23 acts before the receiver, so it is what keeps the front end out of
compression, and no amount of processing can substitute for it.  Cancellation
acts after digitisation and can only remove what was linearly recorded.

Intended approach
-----------------
The leakage is a near-zero-delay copy of the transmitted sweep, so after
de-chirping it appears as a strong component at near-zero beat frequency with a
skirt set by the window.  Three levels of treatment, in increasing order of
what they demand:

1.  **Exclusion.**  Ignore the near range.  This is what
    :func:`tss.l3_signal.trace_extract.leading_edge` does with its
    ``exclude_below_m`` argument.  It is not cancellation and it loses the near
    range, which at the top of the band is where the topside echo is.
2.  **Subtraction of an estimated reference.**  Estimate the leakage transfer
    function, which is the complex coupling between the transmit and receive
    antennas as a function of frequency, and subtract the predicted leakage.
    The coupling varies across a 0.33 to 20 MHz band by a large factor because
    the antennas are electrically small at the bottom of it, so the estimate is
    a frequency-dependent function and not a constant.
3.  **Adaptive cancellation.**  Update the estimate continuously against the
    received data, which is what is needed if the coupling changes in flight.
    It does change: the antenna impedance depends on the surrounding plasma,
    and the plasma changes along the orbit.  That dependence is the reason the
    cancellation has to be adaptive rather than calibrated once on the bench.

Open questions
--------------
1.  **How deep cancellation can go before the residual is set by something
    other than the estimate.**  Cancellation cannot remove what the receiver
    did not record linearly.  Front-end compression, quantisation and the
    receiver's own intermodulation all put a floor under the residual that no
    adaptive filter reaches, and where that floor is decides whether INS-24 is
    achievable at all.
2.  **Whether the residual is stationary enough to adapt to.**  The plasma
    impedance changes on the orbit timescale, which is slow, but a depletion
    crossing is fast.  An adaptive canceller that is still converging across a
    depletion crossing loses the measurement exactly where the science is.
3.  **The interaction with C3.**  Cancellation modifies the signal before the
    Doppler separation sees it.  If it introduces a frequency-dependent phase,
    it introduces an apparent range bias, and that bias would be
    indistinguishable from the Doppler bias the pairing of
    :mod:`tss.l3_signal.novel.sweep_pairing` is trying to remove.

This is a hardware claim with a modelling component, not a modelling claim.
``docs/verification.md`` says so under V4 and the dissertation should present it
that way: the residual numbers that C4 turns on are ultimately measured on the
bench in WP9, and what this module can contribute is how much of the remaining
gap processing can close.

Raises
------
NotImplementedError
    On every call.  WP5.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

import numpy as np

from tss.interfaces import ReceivedSignal

__all__ = [
    "estimate_leakage",
    "cancel_leakage",
    "residual_floor_db",
    "cancellation_margin_db",
]

_MESSAGE = (
    "WP5: adaptive leakage cancellation is not implemented.  This is claim C4, the "
    "claim the concept now turns on: the leakage path has essentially zero delay, "
    "so it gets no range-correlation credit and carries the oscillator's full noise "
    "into the receiver while being tens of decibels stronger than any echo.  "
    "Requirements INS-23, INS-24 and INS-11.  The measured residual is a WP9 bench "
    "result; see docs/verification.md section V4 and docs/work-packages.md rows WP5 "
    "and WP9."
)


def estimate_leakage(
    signal: ReceivedSignal,
    *,
    isolation_db: float | np.ndarray | None = None,
) -> np.ndarray:
    """Estimate the leakage component of a received record.  Not implemented.

    Intended to return the complex leakage waveform, same shape as
    ``signal.samples``, predicted from a frequency-dependent coupling estimate.
    ``isolation_db`` is the INS-23 isolation profile when one is known from L2;
    passing it does not make the estimate a calibration, because the coupling
    phase matters as much as its magnitude and the isolation profile carries
    only the magnitude.

    Raises
    ------
    NotImplementedError
        Always.  WP5.
    """
    del signal, isolation_db
    raise NotImplementedError(_MESSAGE)


def cancel_leakage(
    signal: ReceivedSignal,
    *,
    method: str = "adaptive",
    settings: Mapping[str, Any] | None = None,
) -> ReceivedSignal:
    """Return a record with the leakage removed.  Not implemented.

    Intended to return a new :class:`~tss.interfaces.ReceivedSignal` with the
    cancellation applied and the achieved cancellation depth recorded in
    ``flags``, so that INS-24 is scored against a measured residual rather than
    against an assumed one.  The input is not modified.

    Raises
    ------
    NotImplementedError
        Always.  WP5.
    """
    del signal, method, settings
    raise NotImplementedError(_MESSAGE)


def residual_floor_db(
    signal: ReceivedSignal,
    *,
    exclude_above_m: float | None = None,
) -> float:
    """De-chirped residual leakage floor in dB.  Not implemented.

    This is the quantity INS-24 is scored against, as the margin between the
    weakest echo INS-10 admits and the residual.  It cannot be computed from a
    simulation alone: the floor is set by front-end compression, quantisation
    and receiver intermodulation as much as by the canceller, and those are
    measured on the bench in WP9.

    Raises
    ------
    NotImplementedError
        Always.  WP5.
    """
    del signal, exclude_above_m
    raise NotImplementedError(_MESSAGE)


def cancellation_margin_db(
    residual_db: float | np.ndarray,
    weakest_echo_db: float | np.ndarray,
) -> np.ndarray:
    """Margin of the weakest admitted echo over the residual.  Not implemented.

    The arithmetic is a subtraction.  It is not offered as a function because
    it would be used with an assumed residual, and an INS-24 verdict computed
    from an assumed residual is the failure mode the requirement was written to
    prevent.

    Raises
    ------
    NotImplementedError
        Always.  WP5.
    """
    del residual_db, weakest_echo_db
    raise NotImplementedError(_MESSAGE)
