"""Layer: L3
WP: WP4
Claims: C1, C2
Status: partial
Requirements: INS-01, INS-04, INS-05, INS-07, INS-10

The baseline L3 tier: de-chirp, transform, extract, assemble an ionogram.

This is the signal chain of WP4 with nothing novel in it.  It does what a
conventional chirp sounder processor does, so that the novel tiers of
:mod:`tss.l3_signal.novel` have something to be compared against and so that
the rest of the pipeline has an ionogram to invert.

How a sweep becomes an ionogram
--------------------------------
The transmitted sweep visits every frequency in the band once.  The processor
cuts the received record into contiguous blocks, de-chirps each block against
the sweep that was actually sent, transforms it, and files the result under the
instantaneous transmitted frequency at the centre of that block.  Stacking the
blocks against that frequency is the ionogram.

The block length is the one free choice and it sets everything else:

* the coherent block duration ``T = n_block / fs``;
* the range resolution, ``c enbw / (2 K_f T)``, so a longer block gives finer
  range cells;
* the noise bandwidth per cell, ``1 / T`` times the window's equivalent noise
  bandwidth, so a longer block also gives a lower noise floor;
* the frequency span the block covers, ``K_f T``, and therefore how far the
  sounding frequency moves while the block is being collected.

The first three want a long block and the last wants a short one.  The block
cannot be so long that the sounding frequency moves across a feature of the
trace, and it cannot be so short that the range cell exceeds INS-05.  At the
baseline ``K_f`` of 5.6 MHz/s, a 1 km range cell needs ``T`` of about 27 ms,
during which the sounding frequency moves 150 kHz.  That is the INS-02 step
cap, which is not a coincidence: the frequency step law and the range
resolution are the same trade seen from two sides.

What this tier does not do
---------------------------
It does not separate O from X.  A single receive channel sees the two traces
superposed and there is no information in one channel with which to label them.
The ionogram therefore carries one mode column, named ``"unresolved"``, and
``flags["mode_discrimination"]`` is False.  L5 scores INS-08 against that flag
and reports a failure, which is the correct report: the baseline configuration
does not meet INS-08 and the vector receive of WP11 is what would.  Filling two
columns with the same trace would score INS-08 as met and would be a fabricated
result.

It does not populate ``doppler_mps``.  A single sweep cannot separate the
Doppler shift from the range, which is the whole of C3.  Leaving the field
``None`` is the tier saying so, and L5 records INS-09 as ``not_evaluated``
rather than as a pass or a failure.  :mod:`tss.l3_signal.motion_separated` is
the tier that will populate it.

It does not populate ``residual_range_bias_m`` either, for the same reason.  A
tier that has not estimated the Doppler does not know how much range bias it
has left in, and writing a zero there would assert that it left none.

Where the hooks are
--------------------
Two operations are marked in the code as hooks, because they are the two places
a fidelity change enters and a subclass or a future tier will want to replace
them without touching the block loop:

``_detect_and_extract``
    Turns the stack of range spectra into one trace sample per frequency.  The
    baseline delegates to :func:`tss.l3_signal.trace_extract.extract_trace`,
    which is the clean-data extractor.  The robust extractor of WP4 replaces
    exactly this call.
``_estimate_doppler``
    Returns ``(doppler_mps, residual_range_bias_m)`` for the assembled trace.
    The baseline returns ``(None, None)``.  The up and down sweep pairing of
    WP5 replaces exactly this call, and it is the only place C3 touches the
    baseline chain.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from typing import Any

import numpy as np

from tss.interfaces import Ionogram, PlatformState, ReceivedSignal
from dataclasses import replace as _replace

from tss import sweep_law
from tss.l3_signal.base import SignalProcessingBase, sweep_duration_s
from tss.l3_signal.dechirp import dechirp
from tss.l3_signal.trace_extract import TraceSample, extract_trace
from tss.l3_signal.transform import RangeSpectrum, transform
from tss.registry import register

__all__ = ["BaselineSignalProcessing", "UNRESOLVED_MODE"]

#: The mode label used when the receiver cannot tell O from X.  A scalar
#: receive chain always uses it.  ``"O"`` would be a guess and the guess is
#: right about half the time near the cusp.
UNRESOLVED_MODE = "unresolved"


@register("l3", "baseline")
class BaselineSignalProcessing(SignalProcessingBase):
    """Conventional de-chirp, transform and leading-edge ionogram formation.

    Parameters
    ----------
    window
        Taper applied before the transform.  See
        :data:`tss.l3_signal.transform.WINDOWS`.  Hann by default, because a
        topside ionogram has a strong near echo whose rectangular-window
        sidelobes would sit on top of the weak far structure.
    n_frequencies
        How many blocks the sweep is cut into, and therefore how many columns
        the ionogram has.  The block length follows from the record length.
        Ignored for a stepped sweep, which is cut one block per step so that
        every block has a single sweep rate and so a single range scale.
    pad_factor
        Zero-padding factor for the transform.  Reduces the bin spacing without
        changing the range resolution, which reduces the straddling error on
        the leading edge at the cost of a longer transform.
    threshold_db
        Detection threshold above the estimated noise floor.
    min_run
        Cells that must be consecutively over threshold before a detection is
        declared.
    exclude_below_m
        Near range excluded from both the floor estimate and the detection, to
        keep the transmit leakage out.  A blunt substitute for cancelling it;
        see :mod:`tss.l3_signal.novel.leakage_cancel`.
    max_jump_m
        Continuity limit for :func:`tss.l3_signal.trace_extract.connect_trace`.
        ``None`` disables the rule, which is the default.
    chain_max_gap, chain_max_jump_m
        Settings for :func:`tss.l3_signal.trace_extract.longest_connected_trace`,
        which keeps the one continuous trace and rejects isolated detections
        below the start of the trace and above penetration.  On by default,
        with up to two empty columns bridged and 150 km of range change allowed
        per column.  ``chain_max_gap=None`` disables it.
    """

    def __init__(
        self,
        *,
        window: str = "hann",
        n_frequencies: int = 64,
        pad_factor: int = 2,
        threshold_db: float = 10.0,
        min_run: int = 2,
        exclude_below_m: float = 0.0,
        max_jump_m: float | None = None,
        chain_max_gap: int | None = 2,
        chain_max_jump_m: float = 1.5e5,
    ) -> None:
        super().__init__(window=window)
        if int(n_frequencies) < 1:
            raise ValueError("n_frequencies must be at least 1")
        if int(pad_factor) < 1:
            raise ValueError("pad_factor must be at least 1")
        self.n_frequencies = int(n_frequencies)
        self.pad_factor = int(pad_factor)
        self.threshold_db = float(threshold_db)
        self.min_run = int(min_run)
        self.exclude_below_m = float(exclude_below_m)
        self.max_jump_m = None if max_jump_m is None else float(max_jump_m)
        self.chain_max_gap = None if chain_max_gap is None else int(chain_max_gap)
        self.chain_max_jump_m = float(chain_max_jump_m)

    # -- manifest -----------------------------------------------------------

    def _settings(self) -> dict[str, Any]:
        return {
            "window": self.window,
            "n_frequencies": self.n_frequencies,
            "pad_factor": self.pad_factor,
            "threshold_db": self.threshold_db,
            "min_run": self.min_run,
            "exclude_below_m": self.exclude_below_m,
            "max_jump_m": self.max_jump_m,
            "chain_max_gap": self.chain_max_gap,
            "chain_max_jump_m": self.chain_max_jump_m,
            "mode_discrimination": False,
            "doppler_separation": False,
            "status": "partial",
        }

    # -- the block loop -----------------------------------------------------

    def _block_plan(self, signal: ReceivedSignal) -> tuple[int, int, bool]:
        """Usable sample count, block length and whether the record was cut.

        The record is cut into ``n_frequencies`` contiguous blocks of equal
        length.  Any samples left over at the end, at most ``n_frequencies - 1``
        of them, are dropped rather than being folded into the last block,
        because a last block of a different length would have a different range
        resolution and a different noise floor from every other column.

        A record longer than one sweep is truncated to one sweep.  Processing
        several sweeps into one ionogram is a coherent integration question and
        it belongs to :mod:`tss.l3_signal.novel.motion_compensation`, not here.

        The plan is computed in one place so that the block duration reported
        in the ionogram flags is the block duration the transforms actually
        used.  Deriving it twice is how the reported range resolution and the
        achieved range resolution come apart.
        """
        samples = np.atleast_2d(np.asarray(signal.samples))
        n_total = int(samples.shape[1])
        fs = float(signal.fs_hz)
        if fs <= 0.0:
            raise ValueError("ReceivedSignal.fs_hz must be positive")

        n_sweep = int(round(sweep_duration_s(signal.sweep) * fs))
        n_use = min(n_total, n_sweep) if n_sweep > 0 else n_total
        n_block = n_use // self.n_frequencies
        if n_block < 2:
            raise ValueError(
                f"{n_use} usable samples cut into {self.n_frequencies} blocks gives "
                f"{n_block} samples per block; a block needs at least two"
            )
        return n_use, n_block, n_total > n_use

    def _block_starts(self, signal: ReceivedSignal) -> tuple[list[int], int]:
        """Start index of every block and the common block length.

        A uniform sweep is cut into ``n_frequencies`` equal contiguous blocks.
        A stepped sweep is cut at its step boundaries, one block per step.  Each
        block starts on the first whole sample of its step and is one sample
        shorter than the step, so that no block reaches into the next step and
        picks up a second rate.  Steps that lie beyond the record, which happens
        when L2 capped the record length, are left out.
        """
        if not sweep_law.is_stepped(signal.sweep):
            _, n_block, _ = self._block_plan(signal)
            return [b * n_block for b in range(self.n_frequencies)], n_block
        samples = np.atleast_2d(np.asarray(signal.samples))
        n_total = int(samples.shape[1])
        fs = float(signal.fs_hz)
        per_step = float(signal.sweep.step_dwell_s) * fs
        n_block = int(np.floor(per_step)) - 1
        if n_block < 2:
            raise ValueError(
                f"a {signal.sweep.step_dwell_s:g} s step at {fs:g} Hz gives "
                f"{n_block} samples per block; a block needs at least two"
            )
        starts = []
        for b in range(sweep_law.step_count(signal.sweep)):
            start = int(np.ceil(b * per_step))
            if start + n_block > n_total:
                break
            starts.append(start)
        if not starts:
            raise ValueError("the record is shorter than one step of the stepped sweep")
        return starts, n_block

    def spectra(self, signal: ReceivedSignal) -> list[RangeSpectrum]:
        """De-chirp and transform the record, one block per ionogram column."""
        starts, n_block = self._block_starts(signal)
        stepped = sweep_law.is_stepped(signal.sweep)
        if stepped:
            # File each column at its step's nominal centre.  The sample-time
            # centre of a block wanders by a fraction of a sample from step to
            # step, and the A1 law sits exactly on INS-02's limit, so that
            # wander alone would score as a step law violation.
            centres = sweep_law.step_centres_hz(signal.sweep)
            if int(signal.sweep.direction) < 0:
                centres = centres[::-1]

        out: list[RangeSpectrum] = []
        for b, start in enumerate(starts):
            block = dechirp(signal, start=start, stop=start + n_block, t0_s=0.0, orient=True)
            if stepped:
                block = _replace(block, centre_frequency_hz=float(centres[b]))
            out.append(
                transform(
                    block,
                    window=self.window,
                    n_fft=n_block * self.pad_factor,
                )
            )
        return out

    # -- hook: detection and trace extraction -------------------------------
    #
    # This is the first of the two replaceable hooks named in the module
    # docstring.  The robust extractor of WP4 replaces this method and nothing
    # else in the tier.

    def _detect_and_extract(self, spectra: Sequence[RangeSpectrum]) -> list[TraceSample]:
        """Stack of range spectra to one trace sample per sounding frequency."""
        return extract_trace(
            list(spectra),
            channel=0,
            threshold_db=self.threshold_db,
            min_run=self.min_run,
            exclude_below_m=self.exclude_below_m,
            max_jump_m=self.max_jump_m,
            chain_max_gap=self.chain_max_gap,
            chain_max_jump_m=self.chain_max_jump_m,
        )

    # -- hook: Doppler separation -------------------------------------------
    #
    # The second replaceable hook.  The up and down sweep pairing of WP5
    # replaces this method.  The baseline has nothing to say, and says nothing
    # rather than saying zero.

    def _estimate_doppler(
        self,
        samples: Sequence[TraceSample],
        signal: ReceivedSignal,
        platform: PlatformState | None,
    ) -> tuple[np.ndarray | None, np.ndarray | None]:
        """Doppler in m/s and residual range bias in m, or ``(None, None)``.

        A single sweep carries the Doppler shift and the range in one beat
        frequency and cannot separate them.  Returning ``None`` is the tier
        reporting that, and L5 then records INS-09 as ``not_evaluated``.
        """
        del samples, signal, platform
        return None, None

    # -- the L3 contract ----------------------------------------------------

    def form_ionogram(
        self,
        signal: ReceivedSignal,
        platform: PlatformState | None = None,
    ) -> Ionogram:
        """Form an ionogram from a received record.  Satisfies L3SignalProcessing.

        Returns a new :class:`~tss.interfaces.Ionogram`.  ``signal`` is not
        modified.  Frequencies with no detection carry NaN in
        ``virtual_range_m`` and ``snr_db``, which is what contract rule 3 in
        :mod:`tss.interfaces` asks for: a frequency at which the wave was
        evanescent or penetrated is a result, not a failure.
        """
        spectra = self.spectra(signal)
        samples = self._detect_and_extract(spectra)
        doppler_mps, residual_bias_m = self._estimate_doppler(samples, signal, platform)

        freq_hz = np.array([s.freq_hz for s in samples], dtype=float)
        virtual_range_m = np.array([s.range_m for s in samples], dtype=float)[:, None]
        snr_db = np.array([s.snr_db for s in samples], dtype=float)[:, None]

        detected = np.isfinite(virtual_range_m[:, 0])
        first = spectra[0]
        stepped = sweep_law.is_stepped(signal.sweep)
        if stepped:
            _, n_block_samples = self._block_starts(signal)
            n_sweep = int(round(sweep_duration_s(signal.sweep) * float(signal.fs_hz)))
            truncated = int(np.atleast_2d(np.asarray(signal.samples)).shape[1]) > n_sweep
        else:
            _, n_block_samples, truncated = self._block_plan(signal)
        block_duration_s = n_block_samples / float(signal.fs_hz)
        resolution_by_column = np.array([sp.range_resolution_m for sp in spectra], dtype=float)
        unambiguous_by_column = np.array([sp.unambiguous_range_m for sp in spectra], dtype=float)

        flags: dict[str, Any] = {
            # Instrument-side numbers L5 scores requirements against.  For a
            # uniform sweep every column has the same values.  For a stepped
            # sweep they vary with the step rate, and the scored values are the
            # worst column: the coarsest resolution (the bottom of the band)
            # and the shortest unambiguous range (the fastest step).
            "range_resolution_m": float(np.max(resolution_by_column)),
            "bin_spacing_m": float(max(sp.bin_spacing_m for sp in spectra)),
            "noise_bandwidth_hz": float(first.noise_bandwidth_hz),
            "unambiguous_range_m": float(np.min(unambiguous_by_column)),
            "sweep_law": "stepped" if stepped else "uniform",
            "range_resolution_m_by_column": resolution_by_column.tolist(),
            "ionogram_period_s": float(signal.sweep.sweep_period_s),
            # The band the transmitted sweep covered.  INS-01a and INS-01b are
            # scored against these, not against the detected edges, which the
            # ionosphere sets (the local plasma frequency and fxF2).
            "sweep_min_frequency_hz": float(
                min(signal.sweep.f_start_hz, signal.sweep.f_stop_hz)
            ),
            "sweep_max_frequency_hz": float(
                max(signal.sweep.f_start_hz, signal.sweep.f_stop_hz)
            ),
            "block_duration_s": float(block_duration_s),
            "window": str(first.window),
            "enbw_bins": float(first.enbw_bins),
            # What the tier did not do, stated so that L5 does not have to
            # infer it from a missing field.
            "mode_discrimination": False,
            "mode_discrimination_reason": (
                "scalar receive; O and X labelling needs the vector receive of WP11 "
                "(INS-08, C5)"
            ),
            "doppler_separated": False,
            "doppler_separation_reason": (
                "a single sweep cannot separate the Doppler shift from the range; "
                "that is claim C3 and it lives in tss.l3_signal.novel.sweep_pairing (WP5)"
            ),
            "n_frequencies": int(freq_hz.size),
            "n_detected": int(np.count_nonzero(detected)),
            "n_rejected_disconnected": int(
                sum(1 for s in samples if s.flags.get("rejected_disconnected"))
            ),
            "input_dechirped": bool(signal.flags.get("dechirped", False)),
            "record_truncated_to_one_sweep": bool(truncated),
        }
        if np.any(detected):
            flags["min_detected_frequency_hz"] = float(np.min(freq_hz[detected]))
            flags["max_detected_frequency_hz"] = float(np.max(freq_hz[detected]))

        return Ionogram(
            freq_hz=freq_hz,
            virtual_range_m=virtual_range_m,
            snr_db=snr_db,
            modes=(UNRESOLVED_MODE,),
            doppler_mps=doppler_mps,
            aoa_deg=None,
            residual_range_bias_m=residual_bias_m,
            flags=flags,
            source=f"tss.l3_signal.ionogram.{type(self).__name__}",
        )

    def describe(self) -> Mapping[str, Any]:
        return super().describe()
