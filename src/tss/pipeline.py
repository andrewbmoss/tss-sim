"""The L0 to L5 orchestrator.

The pipeline owns the sequencing and nothing else.  It builds the tier
implementations a campaign names, runs one configuration end to end, and hands
the result to the assessment tier.  It contains no physics, so a change of
fidelity in any tier requires no change here.

Two properties are deliberate.  First, every intermediate product is retained
when ``keep_intermediates`` is set, because verification level V1 compares a
layer against the one beneath it and cannot do that if only the final profile
survives.  Second, the RNG is created once from the campaign seed and split per
configuration, so two configurations in the same campaign do not share a noise
realisation and a single configuration reruns identically.

Between L1 and L2 the pipeline applies the clutter generator of WP4 Phase 1
(:mod:`tss.l1_propagation.clutter`) when it is given one.  It is not a tier: it
adds returns that are not echoes of the profile to whatever echo set L1 made,
so it works with any L1 tier.  It draws from a generator of its own, derived
from the same seed and configuration with a further spawn key, so clutter on or
off does not by itself redraw L2's noise.  The default is no clutter; the
campaign runner passes one when the configuration's ``clutter`` block is on.
"""

from __future__ import annotations

import zlib
from collections.abc import Mapping
from dataclasses import dataclass, field, replace
from typing import Any

import numpy as np

from tss.interfaces import (
    Assessment,
    DensityProfile,
    EchoSet,
    Ionogram,
    IonosphereState,
    PlatformState,
    ReceivedSignal,
    SweepSchedule,
)
from tss.registry import build

__all__ = ["Pipeline", "RunResult", "build_tiers"]


def _accepts(fn: Any, name: str) -> bool:
    """Whether ``fn`` takes a keyword ``name``."""
    import inspect

    try:
        params = inspect.signature(fn).parameters
    except (TypeError, ValueError):
        return False
    return name in params or any(p.kind is inspect.Parameter.VAR_KEYWORD for p in params.values())


@dataclass
class RunResult:
    """Everything one configuration produced."""

    configuration: str
    assessment: Assessment | None = None
    truth: IonosphereState | None = None
    echoes: EchoSet | None = None
    signal: ReceivedSignal | None = None
    ionogram: Ionogram | None = None
    retrieved: DensityProfile | None = None
    timing_s: dict[str, float] = field(default_factory=dict)
    flags: dict[str, Any] = field(default_factory=dict)
    truth_echo_snr: dict[str, Any] | None = None


def build_tiers(tiers: Mapping[str, Any]) -> dict[str, Any]:
    """Instantiate every tier a campaign names.

    A tier entry is either a bare name (``l0: chapman``) or a mapping with a
    ``name`` and keyword arguments (``l0: {name: iri2020, topside: NeQuick}``).
    """
    built: dict[str, Any] = {}
    for layer, spec in tiers.items():
        if isinstance(spec, str):
            built[layer] = build(layer, spec)
        elif isinstance(spec, Mapping):
            kwargs = {k: v for k, v in spec.items() if k != "name"}
            built[layer] = build(layer, spec["name"], **kwargs)
        else:
            raise TypeError(f"tier spec for {layer} must be a name or a mapping")
    return built


class Pipeline:
    """Runs configurations through L0 to L5."""

    def __init__(self, tiers: Mapping[str, Any], *, seed: int = 0,
                 keep_intermediates: bool = True, clutter: Any = None) -> None:
        self.tiers = dict(tiers)
        self.seed = int(seed)
        self.keep_intermediates = keep_intermediates
        self.clutter = clutter

    def _clutter_rng(self, configuration: str) -> np.random.Generator:
        """The clutter generator's own stream (spawn key 0xC1 after the configuration's)."""
        mixed = np.random.SeedSequence(
            entropy=self.seed,
            spawn_key=(zlib.crc32(configuration.encode("utf-8")), 0xC1),
        )
        return np.random.default_rng(mixed)

    def _rng(self, configuration: str) -> np.random.Generator:
        """A generator unique to (seed, configuration) and stable across runs.

        The configuration name is folded in with CRC-32 rather than the built-in
        ``hash``, which Python salts per process (``PYTHONHASHSEED``).  With
        ``hash`` the same seed drew different noise in every process, so no two
        runs of a campaign were the same run.
        """
        mixed = np.random.SeedSequence(
            entropy=self.seed,
            spawn_key=(zlib.crc32(configuration.encode("utf-8")),),
        )
        return np.random.default_rng(mixed)

    def run_configuration(
        self,
        configuration: str,
        *,
        platform: PlatformState,
        sweep: SweepSchedule,
        freq_hz: np.ndarray,
        alt_m: np.ndarray,
        epoch: str,
    ) -> RunResult:
        """Run one configuration end to end."""
        import time

        rng = self._rng(configuration)
        result = RunResult(configuration=configuration)
        t: dict[str, float] = {}

        t0 = time.perf_counter()
        truth = self.tiers["l0"].state_at(
            platform.lat_deg, platform.lon_deg, epoch, alt_m
        )
        t["l0"] = time.perf_counter() - t0

        t0 = time.perf_counter()
        echoes = self.tiers["l1"].echoes(truth, platform, freq_hz, rng)
        t["l1"] = time.perf_counter() - t0

        if self.clutter is not None:
            t0 = time.perf_counter()
            echoes = self.clutter.augment(echoes, truth, platform, freq_hz, self._clutter_rng(configuration))
            t["clutter"] = time.perf_counter() - t0
            result.flags["clutter"] = dict(getattr(self.clutter, "last_report", {}) or {})

        t0 = time.perf_counter()
        signal = self.tiers["l2"].receive(echoes, sweep, rng)
        t["l2"] = time.perf_counter() - t0

        t0 = time.perf_counter()
        ionogram = self.tiers["l3"].form_ionogram(signal, platform)
        t["l3"] = time.perf_counter() - t0

        # INS-10 is scored at every echo the truth returned (2026-10-17).  The
        # cells come from the echo set and the measurement from L3, which is
        # the layer that knows how a range cell's SNR is measured; the
        # pipeline only hands one to the other.  A tier without the
        # measurement leaves INS-10 not evaluated.
        truth_echo_snr: dict[str, Any] | None = None
        measure = getattr(self.tiers["l3"], "snr_at", None)
        if callable(measure):
            from tss.l5_assessment.scoring import truth_echo_cells

            t0 = time.perf_counter()
            cells = truth_echo_cells(echoes)
            truth_echo_snr = {
                **cells,
                "snr_db": np.asarray(measure(signal, cells["freq_hz"], cells["range_m"])),
            }
            # WP4 Phase 2 (G1): INS-10 is scored against the noise floor; the
            # same cells measured against the column, which a clutter continuum
            # raises, are reported beside it as signal to noise plus clutter.
            if _accepts(measure, "floor"):
                truth_echo_snr["snr_db_column"] = np.asarray(
                    measure(signal, cells["freq_hz"], cells["range_m"], floor="column"))
            t["l3_truth_echo_snr"] = time.perf_counter() - t0

        t0 = time.perf_counter()
        retrieved = self.tiers["l4"].invert(ionogram, platform)
        t["l4"] = time.perf_counter() - t0

        # L5 scores a profile against SCI-05 and SCI-07 only when the echoes
        # and the inversion treat the magnetic field the same way.  The
        # retrieval cannot know how its echoes were made, so that is recorded
        # here, where both tiers are in hand.  This is bookkeeping, not physics.
        forward = getattr(self.tiers["l1"], "group_index_used", None)
        if forward is not None:
            retrieved = replace(
                retrieved,
                flags={**dict(retrieved.flags or {}), "forward_group_index": str(forward)},
            )

        # Requirement evidence that only L2 can give (WP3 Phase 4): the
        # leakage scores and the declared values, carried on the record.  Only
        # an assessment tier that takes them is given them.
        l5_kwargs: dict[str, Any] = {"platform": platform, "truth_echo_snr": truth_echo_snr}
        l2_metrics = (signal.flags or {}).get("l2_metrics") if signal is not None else None
        if l2_metrics is not None and _accepts(self.tiers["l5"].score, "l2_metrics"):
            l5_kwargs["l2_metrics"] = l2_metrics

        t0 = time.perf_counter()
        assessment = self.tiers["l5"].score(
            retrieved, truth, ionogram, configuration, **l5_kwargs
        )
        t["l5"] = time.perf_counter() - t0

        result.assessment = assessment
        result.timing_s = t
        result.truth_echo_snr = truth_echo_snr
        if self.keep_intermediates:
            result.truth = truth
            result.echoes = echoes
            result.signal = signal
            result.ionogram = ionogram
            result.retrieved = retrieved
        return result

    def describe(self) -> dict[str, Any]:
        """Tier descriptions for the run manifest."""
        out: dict[str, Any] = {"seed": self.seed}
        for layer, tier in self.tiers.items():
            out[layer] = dict(tier.describe()) if hasattr(tier, "describe") else {
                "name": getattr(tier, "name", type(tier).__name__)
            }
        out["clutter"] = dict(self.clutter.describe()) if self.clutter is not None else {"enabled": False}
        return out
