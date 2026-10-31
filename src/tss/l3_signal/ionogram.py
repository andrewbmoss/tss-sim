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
    exactly this call.  Since WP4 Phase 2 its cell decisions come from the
    stated-P_fa detector of :mod:`tss.l3_signal.detect` (``detector="cfar"``,
    the default; ``"fixed"`` is the WP3-era 10 dB rule).
``_estimate_doppler``
    Returns ``(doppler_mps, residual_range_bias_m)`` for the assembled trace.
    The baseline returns ``(None, None)``.  The up and down sweep pairing of
    WP5 replaces exactly this call, and it is the only place C3 touches the
    baseline chain.
"""

from __future__ import annotations

from collections.abc import Mapping, Sequence
from dataclasses import replace as _replace
from typing import Any

import numpy as np

from tss import sweep_law
from tss.interfaces import Ionogram, PlatformState, ReceivedSignal
from tss.l3_signal.base import SignalProcessingBase, sweep_duration_s
from tss.l3_signal.dechirp import dechirp
from tss.l3_signal.trace_extract import TraceSample, extract_trace, noise_floor_db
from tss.l3_signal.transform import RangeSpectrum, transform
from tss.registry import register

_C_LIGHT = 299_792_458.0

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
    detector
        ``"cfar"`` (default since WP4 Phase 2, gate G1): the cell test of
        :mod:`tss.l3_signal.detect`, with ``pfa`` per range cell against the
        negative-beat noise reference and sidelobe blanking under stronger
        peaks.  ``"fixed"``: ``threshold_db`` over the median of the column,
        the WP3-era rule.
    pfa, local_pfa, local_ref_cells, guard_cells, local_order, sidelobe_margin_db
        The detector's settings (decision G1: 1e-4 per cell).  ``local_pfa``
        turns on the optional ordered-statistic veto (16 cells a side beyond 2
        guard cells, at the 0.75 statistic); it is off by default because it
        masks echoes spread in range by dispersion.  Sidelobe blanking at a
        6 dB margin is on.
    max_range_m
        Search limit, INS-04's 1050 km by default.  ``None`` searches to the
        unambiguous range.
    dispersion
        WP4 Phase 3, the treatment of dispersion within a dwell
        (Technical Volume section 5.3): ``"none"`` (one transform of the whole
        dwell, the baseline), ``"blocks"`` (``dispersion_blocks`` sub-blocks
        transformed and added in power, option 1) or ``"rate_search"`` (the
        block de-chirped again against beat chirp rates of ``j / T^2`` for
        ``j`` from -2 to ``rate_search_max_cells``, and the transform with the
        highest peak kept for the column, option 2).  Option 3, narrower steps, is a sweep setting, not an L3 one.
    transition
        WP4 Phase 3, the treatment of the first ``tau`` of each step:
        ``"accept"`` (the baseline, the decision of 2026-09-30), ``"gate"``
        (drop the first ``2 R_max / c`` of every dwell) or ``"matched"`` (the
        hybrid matched filter: in the steps whose rate differs from the step
        before, re-evaluate the cells around each detection by a matched filter
        on the exact de-chirped history of an echo at that delay).
    snr_floor
        What :meth:`snr_at` measures against: ``"noise"`` (the detector's
        negative-beat reference, decision G1) or ``"column"`` (the median of the
        column, which includes clutter: the WP3-era measure).  Defaults to
        ``"noise"`` with the CFAR detector and ``"column"`` with the fixed one.
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
        detector: str = "cfar",
        pfa: float = 1.0e-4,
        local_pfa: float | None = None,
        local_ref_cells: int = 16,
        guard_cells: int = 2,
        local_order: float = 0.75,
        sidelobe_margin_db: float | None = 6.0,
        max_range_m: float | None = 1.05e6,
        snr_floor: str | None = None,
        dispersion: str = "none",
        dispersion_blocks: int = 4,
        rate_search_max_cells: int = 32,
        transition: str = "accept",
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
        if detector not in ("cfar", "fixed"):
            raise ValueError(f"unknown detector {detector!r}; expected 'cfar' or 'fixed'")
        self.detector = str(detector)
        self.pfa = float(pfa)
        self.local_pfa = None if local_pfa is None else float(local_pfa)
        self.local_ref_cells = int(local_ref_cells)
        self.guard_cells = int(guard_cells)
        self.local_order = float(local_order)
        self.sidelobe_margin_db = None if sidelobe_margin_db is None else float(sidelobe_margin_db)
        self.max_range_m = None if max_range_m is None else float(max_range_m)
        floor = snr_floor or ("noise" if self.detector == "cfar" else "column")
        if floor not in ("noise", "column"):
            raise ValueError(f"unknown snr_floor {snr_floor!r}; expected 'noise' or 'column'")
        self.snr_floor = floor
        if dispersion not in ("none", "blocks", "rate_search"):
            raise ValueError(f"unknown dispersion {dispersion!r}; expected 'none', 'blocks' or 'rate_search'")
        if transition not in ("accept", "gate", "matched"):
            raise ValueError(f"unknown transition {transition!r}; expected 'accept', 'gate' or 'matched'")
        self.dispersion = str(dispersion)
        self.dispersion_blocks = max(int(dispersion_blocks), 1)
        self.rate_search_max_cells = max(int(rate_search_max_cells), 0)
        self.transition = str(transition)

    def detector_settings(self) -> dict[str, Any]:
        """Keyword arguments for :func:`tss.l3_signal.detect.cfar_detect`."""
        return {"pfa": self.pfa, "local_pfa": self.local_pfa, "local_ref_cells": self.local_ref_cells,
                "guard_cells": self.guard_cells, "local_order": self.local_order,
                "sidelobe_margin_db": self.sidelobe_margin_db,
                "max_range_m": self.max_range_m, "exclude_below_m": self.exclude_below_m}

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
            "detector": self.detector,
            "pfa": self.pfa if self.detector == "cfar" else None,
            "local_pfa": self.local_pfa if self.detector == "cfar" else None,
            "local_ref_cells": self.local_ref_cells,
            "guard_cells": self.guard_cells,
            "local_order": self.local_order,
            "sidelobe_margin_db": self.sidelobe_margin_db if self.detector == "cfar" else None,
            "max_range_m": self.max_range_m,
            "snr_floor": self.snr_floor,
            "dispersion": self.dispersion,
            "dispersion_blocks": self.dispersion_blocks if self.dispersion == "blocks" else None,
            "rate_search_max_cells": self.rate_search_max_cells if self.dispersion == "rate_search" else None,
            "transition": self.transition,
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

        gate = 0
        if self.transition == "gate":
            reach = self.max_range_m if self.max_range_m is not None else 1.05e6
            gate = int(np.ceil(2.0 * reach / _C_LIGHT * float(signal.fs_hz)))
            if n_block - gate < 16:
                raise ValueError("gating the first 2 R_max / c leaves no dwell to transform")
        transition_steps = self._transition_steps(signal) if self.transition == "matched" else None
        out: list[RangeSpectrum] = []
        for b, start in enumerate(starts):
            s0, n = start + gate, n_block - gate
            block = dechirp(signal, start=s0, stop=s0 + n, t0_s=0.0, orient=True)
            if stepped:
                block = _replace(block, centre_frequency_hz=float(centres[b]))
            sp = self._transform_block(block, n)
            if transition_steps is not None and transition_steps[b]:
                sp = self._matched_refine(sp, signal, s0, n)
            out.append(sp)
        return out

    # -- WP4 Phase 3: dispersion and step-transition treatments ----------------

    def _transform_block(self, block: Any, n: int) -> RangeSpectrum:
        """One block to a range spectrum under the configured dispersion treatment."""
        n_fft = n * self.pad_factor
        if self.dispersion == "none":
            return transform(block, window=self.window, n_fft=n_fft)
        if self.dispersion == "blocks":
            m = min(self.dispersion_blocks, max(n // 16, 1))
            sub = n // m
            pos = neg = None
            first = None
            for j in range(m):
                piece = _replace(block, samples=block.samples[:, j * sub:(j + 1) * sub],
                                 t_s=block.t_s[j * sub:(j + 1) * sub])
                s = transform(piece, window=self.window, n_fft=n_fft)
                first = s if first is None else first
                p = 10.0 ** (np.asarray(s.power_db) / 10.0)
                q = 10.0 ** (np.asarray(s.negative_power_db) / 10.0)
                pos = p if pos is None else pos + p
                neg = q if neg is None else neg + q
            with np.errstate(divide="ignore"):
                return _replace(first, power_db=10.0 * np.log10(pos / m),
                                negative_power_db=10.0 * np.log10(neg / m), n_incoherent=m)
        # Rate search: de-chirp again against beat chirp rates j / T^2 about the
        # block centre, so that an echo whose beat slides by j cells in the dwell
        # is compressed back into one cell at the step-centre range.
        fs = float(block.fs_hz)
        t_rel = (np.arange(n, dtype=float) - 0.5 * (n - 1)) / fs
        dur = n / fs
        # One rate per column, the one whose transform has the highest peak in
        # the search range.  A per-cell maximum over rates was tried first: the
        # mismatched rates spread an echo ahead of itself and the maximum keeps
        # that spread, so the leading edge read early (2026-10-30).
        base = transform(block, window=self.window, n_fft=n_fft)
        rng = np.asarray(base.range_m, dtype=float)
        search = rng <= (self.max_range_m if self.max_range_m is not None else rng[-1])
        best, best_peak, best_gamma = base, float(np.max(np.asarray(base.power_db)[0][search])), 0.0
        js = [j for j in range(-2, self.rate_search_max_cells + 1) if j != 0]
        for j in js:
            gamma = j / dur**2
            chirp = np.exp(-1j * np.pi * gamma * t_rel**2)
            s = transform(_replace(block, samples=block.samples * chirp[None, :]),
                          window=self.window, n_fft=n_fft)
            peak = float(np.max(np.asarray(s.power_db)[0][search]))
            if peak > best_peak:
                best, best_peak, best_gamma = s, peak, gamma
        # The noise reference stays the unchirped transform's (white noise is
        # unchanged by the chirp); the selection over rates is a choice among
        # looks, which the detector's per-look P_fa accounts for.
        return _replace(best, negative_power_db=base.negative_power_db, n_looks=len(js) + 1,
                        rate_hz_per_s2=np.full(rng.size, best_gamma))

    def _transition_steps(self, signal: ReceivedSignal) -> np.ndarray:
        """Per block, whether its step's rate differs from the step before it."""
        if not sweep_law.is_stepped(signal.sweep):
            return np.zeros(self.n_frequencies, dtype=bool)
        starts, _ = self._block_starts(signal)
        dwell = float(signal.sweep.step_dwell_s)
        k = np.asarray(sweep_law.rate_hz_per_s(signal.sweep, (np.arange(len(starts)) + 0.5) * dwell), dtype=float)
        prev = np.concatenate([[k[0]], k[:-1]])
        return np.abs(k - prev) > 1e-6 * np.maximum(k, 1.0)

    def _matched_refine(self, sp: RangeSpectrum, signal: ReceivedSignal, start: int, n: int) -> RangeSpectrum:
        """The hybrid matched filter: re-evaluate the cells around each detection.

        Each candidate delay ``tau`` gets ``|sum w x conj(s_tau)|^2 / (sum w)^2``
        with ``s_tau(t) = exp(j [phi(t) - phi(t - tau)])`` from the transmitted
        law, which carries the step transition exactly (zero before ``tau`` from
        the start of the sweep).  The scaling is the transform's, so a unit echo
        reads 0 dB and the noise per cell is unchanged.  Non-dispersive by
        construction: near f_N and the cusp the delay moves in the dwell and the
        filter is mismatched (recorded drawback, 2026-09-30).
        """
        from tss.l3_signal.detect import cfar_detect
        from tss.l3_signal.trace_extract import _runs
        from tss.l3_signal.transform import window_taper

        res = cfar_detect(sp, **self.detector_settings())
        runs = _runs(res.mask)
        if not runs:
            return sp
        fs = float(signal.fs_hz)
        x = np.atleast_2d(np.asarray(signal.samples))[0, start:start + n].astype(complex)
        t = (start + np.arange(n, dtype=float)) / fs
        w = window_taper(self.window, n)
        norm = float(np.sum(w)) ** 2
        phi_t = sweep_law.phase_rad(signal.sweep, t)
        if int(signal.sweep.direction) < 0:
            x = np.conj(x)
            phi_t = -phi_t
        margin = int(np.ceil(4 * self.pad_factor))
        power = np.array(sp.power_db, dtype=float, copy=True)
        rng = np.asarray(sp.range_m, dtype=float)
        done = np.zeros(rng.size, dtype=bool)
        for a, b in runs:
            for i in range(max(a - margin, 0), min(b + margin, rng.size)):
                if done[i]:
                    continue
                done[i] = True
                tau = 2.0 * rng[i] / _C_LIGHT
                live = t >= tau
                ph = sweep_law.phase_rad(signal.sweep, np.clip(t - tau, 0.0, None))
                if int(signal.sweep.direction) < 0:
                    ph = -ph
                s = np.where(live, np.exp(1j * (phi_t - ph)), 0.0)
                val = abs(np.sum(w * x * np.conj(s))) ** 2 / norm
                power[0, i] = 10.0 * np.log10(val) if val > 0.0 else -np.inf
        return _replace(sp, power_db=power)

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
            exclude_above_m=None,
            max_jump_m=self.max_jump_m,
            chain_max_gap=self.chain_max_gap,
            chain_max_jump_m=self.chain_max_jump_m,
            detector=self.detector,
            detector_settings=self.detector_settings(),
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

    def snr_at(
        self,
        signal: ReceivedSignal,
        freq_hz: Sequence[float] | np.ndarray,
        range_m: Sequence[float] | np.ndarray,
        *,
        half_width_cells: float = 2.0,
        floor: str | None = None,
    ) -> np.ndarray:
        """Signal to noise ratio per range cell at given (frequency, range) cells.

        A measurement, not a detection: for each cell the column whose step
        contains ``freq_hz`` is taken, its noise floor is estimated exactly as
        the detector estimates it, and the result is the peak power within
        ``half_width_cells`` range cells of ``range_m`` above that floor.  It is
        what INS-10 is scored on at every echo the truth returned, which a
        detected-cell minimum cannot be, because a cell is detected only if it
        is already above the 10 dB threshold.

        NaN where the frequency is outside every column or the range is beyond
        the column's unambiguous range, which is where the echo has folded and
        its range cell cannot be measured.  ``signal`` is not modified.

        ``floor`` (default: the tier's ``snr_floor``) chooses the floor:
        ``"noise"``, the detector's negative-beat reference, which is noise and
        carriers and no clutter (decision G1, 2026-10-30), or ``"column"``, the
        median of the column's positive half plus the median-to-mean factor,
        which a clutter continuum in the column raises (the WP3-era measure,
        reported beside INS-10 as signal to noise plus clutter).
        """
        which = floor or self.snr_floor
        if which not in ("noise", "column"):
            raise ValueError(f"unknown floor {floor!r}; expected 'noise' or 'column'")
        f = np.atleast_1d(np.asarray(freq_hz, dtype=float))
        r = np.atleast_1d(np.asarray(range_m, dtype=float))
        if f.shape != r.shape:
            raise ValueError("freq_hz and range_m must be parallel")
        out = np.full(f.shape, np.nan)
        if f.size == 0:
            return out
        spectra = self.spectra(signal)
        centres = np.array([sp.centre_frequency_hz for sp in spectra], dtype=float)
        order = np.argsort(centres)
        edges = np.concatenate((
            [centres[order[0]] - 0.5 * (centres[order[1]] - centres[order[0]])
             if centres.size > 1 else -np.inf],
            0.5 * (centres[order[1:]] + centres[order[:-1]]),
            [centres[order[-1]] + 0.5 * (centres[order[-1]] - centres[order[-2]])
             if centres.size > 1 else np.inf],
        ))
        floors: dict[int, float] = {}
        for i, (fi, ri) in enumerate(zip(f, r, strict=True)):
            if not (np.isfinite(fi) and np.isfinite(ri)):
                continue
            slot = int(np.searchsorted(edges, fi, side="right")) - 1
            if not 0 <= slot < centres.size:
                continue
            j = int(order[slot])
            sp = spectra[j]
            if ri > sp.unambiguous_range_m:
                continue
            row = np.asarray(sp.power_db, dtype=float)[0]
            rng = np.asarray(sp.range_m, dtype=float)
            if j not in floors:
                if which == "noise":
                    from tss.l3_signal.detect import noise_reference

                    floors[j] = float(noise_reference(sp)["mean_db"])
                else:
                    floors[j] = noise_floor_db(row, rng, exclude_below_m=self.exclude_below_m)
            near = np.abs(rng - ri) <= float(half_width_cells) * sp.range_resolution_m
            if not np.any(near):
                continue
            out[i] = float(np.max(row[near]) - floors[j])
        return out

    def describe(self) -> Mapping[str, Any]:
        return super().describe()
