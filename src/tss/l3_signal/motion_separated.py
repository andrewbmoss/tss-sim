"""Layer: L3
WP: WP5
Claims: C3
Status: stub
Requirements: INS-09, INS-15, SCI-23

The motion-separated L3 tier.  Registered, not implemented.

This is the tier ``config/campaigns/c3_doppler.yaml`` names.  It is registered
so that the campaign validates, and it raises when it is run.  That is
deliberate and it is worth saying why rather than leaving it to be discovered.

The campaign exists before the algorithm on purpose.  The shape of the answer
C3 has to produce is fixed by the requirement set: a residual apparent-range
bias under one INS-05 range cell, a velocity resolution of 50 m/s at threshold
and 25 m/s at goal from SCI-23, and a coherent interval of 1 to 2 s from
INS-15, scored per configuration.  Writing the campaign first fixes what the
result must look like before there is a result, which is the first governing
principle of Section 3.1 applied to the author's own work: the number a result
is scored against has to exist before the result does, and that discipline is
easier to keep when the campaign file was written first.

``tss config c3_doppler`` therefore succeeds, because the tier is registered,
and ``tss run c3_doppler`` fails with a message naming WP5, because the
algorithm is not written.  Those are the correct two behaviours and the
registry's contract, stated in :mod:`tss.registry` and in ``CONTRIBUTING.md``,
is what makes them possible: a stub raises when it is called and not when it is
built.

What the tier would compose
----------------------------
1.  :func:`tss.l3_signal.novel.leakage_cancel.cancel_leakage` on the record,
    because everything downstream works on what survives the leakage.
2.  The baseline chain of :class:`tss.l3_signal.ionogram.BaselineSignalProcessing`
    on the up-sweep and on the down-sweep separately, giving two ionograms whose
    traces are displaced in opposite directions by the same amount.
3.  :func:`tss.l3_signal.novel.sweep_pairing.pair_sweeps` to associate the two
    traces and take the half-sum and half-difference.
4.  :func:`tss.l3_signal.novel.motion_compensation.integrate_coherently` over
    the interval INS-15 asks for, when a velocity resolution rather than a
    single velocity estimate is wanted.

Step 3 is the principal contribution and step 1 is the claim the concept turns
on.  Neither exists, so neither does this.

Why it is not a thin wrapper that raises immediately
------------------------------------------------------
It could have been written to raise in its constructor, which would be simpler.
It raises in ``form_ionogram`` instead so that ``tss tiers`` lists it, so that
the campaign validator accepts it, and so that a reader can see from
:meth:`describe` what settings the tier will take when it exists.  Recording
those settings now means the first run that does work has them in its manifest,
rather than having them added later and absent from the earlier runs.

Raises
------
NotImplementedError
    On every call to ``form_ionogram``.  WP5.
"""

from __future__ import annotations

from typing import Any

from tss.interfaces import Ionogram, PlatformState, ReceivedSignal
from tss.l3_signal.base import SignalProcessingBase
from tss.registry import register

__all__ = ["MotionSeparatedSignalProcessing"]

_MESSAGE = (
    "WP5: the motion-separated L3 tier is not implemented.  It composes "
    "tss.l3_signal.novel.leakage_cancel (C4), the baseline chain on the up and down "
    "sweeps separately, tss.l3_signal.novel.sweep_pairing (C3, the principal "
    "contribution) and tss.l3_signal.novel.motion_compensation.  The tier is "
    "registered so that config/campaigns/c3_doppler.yaml validates before the "
    "algorithms exist, which is deliberate; see the module docstring and "
    "docs/work-packages.md row WP5."
)


@register("l3", "motion_separated")
class MotionSeparatedSignalProcessing(SignalProcessingBase):
    """Range and Doppler separated using the platform motion.  Stub.

    Parameters
    ----------
    window
        Taper for the underlying transforms, passed through to the baseline
        chain.
    coherent_interval_s
        Coherent processing interval.  INS-15 asks for 1 s at threshold and 2 s
        at goal, and the default is the threshold.
    sweep_pairs
        How up-sweeps and down-sweeps are paired.  ``"up_down"`` is the
        alternating schedule ``config/campaigns/c3_doppler.yaml`` names.
    leakage_cancellation
        Whether the leakage canceller runs before the chain.  Defaults on,
        because a C3 result taken with the leakage still present would be a
        result about a different instrument.
    """

    def __init__(
        self,
        *,
        window: str = "hann",
        coherent_interval_s: float = 1.0,
        sweep_pairs: str = "up_down",
        leakage_cancellation: bool = True,
    ) -> None:
        super().__init__(window=window)
        self.coherent_interval_s = float(coherent_interval_s)
        self.sweep_pairs = str(sweep_pairs)
        self.leakage_cancellation = bool(leakage_cancellation)

    def _settings(self) -> dict[str, Any]:
        return {
            "window": self.window,
            "coherent_interval_s": self.coherent_interval_s,
            "sweep_pairs": self.sweep_pairs,
            "leakage_cancellation": self.leakage_cancellation,
            "status": "stub",
        }

    def form_ionogram(
        self,
        signal: ReceivedSignal,
        platform: PlatformState | None = None,
    ) -> Ionogram:
        """Form a Doppler-separated ionogram.  Not implemented.

        Raises
        ------
        NotImplementedError
            Always.  WP5.
        """
        del signal, platform
        raise NotImplementedError(_MESSAGE)
