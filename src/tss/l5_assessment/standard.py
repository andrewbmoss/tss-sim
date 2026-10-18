"""Layer: L5
WP: WP0, WP7
Claims: C1, C2, C3, C4, C5
Status: implemented
Requirements: INS-02, INS-04, INS-05, INS-07, INS-08, INS-09, INS-10, INS-15, SCI-09, SCI-23

The standard L5 tier: compute the metrics, score them, return an Assessment.

It composes :mod:`tss.l5_assessment.scoring` and
:mod:`tss.l5_assessment.requirements_check` and adds nothing of its own.  That
is the whole design: the arithmetic and the judging live in modules that can be
tested without a pipeline, and the tier is the thin thing that the pipeline
calls.

The returned :class:`~tss.interfaces.Assessment` carries the configuration
name, because Section 1.3 of the proposal reports per configuration and a
result without its configuration is not a result.  The campaign validator
refuses a campaign that names no configurations for the same reason, so a run
that reaches this tier always has one.

``Assessment.errors`` carries the truth-versus-retrieval comparison and the
achieved instrument numbers, including the ones no requirement is scored
against.  Nothing computed is dropped merely because the requirement schema
cannot express it; see the INS-01 note in
:mod:`tss.l5_assessment.requirements_check`.

``Assessment.flags`` carries the summary verdict counts, the requirement set
revision, and any flags the ionogram or the retrieved profile brought with them,
prefixed so that their origin stays visible.  Two of those are load bearing:
``l4.starting_value`` says which of the three starting-value routes the
inversion used, and ``l3.mode_discrimination`` says whether the O and X traces
were labelled.  Both change how a number should be read.
"""

from __future__ import annotations

from collections.abc import Mapping
from typing import Any

from tss.interfaces import (
    Assessment,
    DensityProfile,
    Ionogram,
    IonosphereState,
    PlatformState,
)
from tss.l5_assessment.base import AssessmentBase, summarise_verdicts
from tss.l5_assessment.requirements_check import CARRY_IN_ERRORS, check_requirements
from tss.l5_assessment.scoring import compute_metrics
from tss.registry import register

__all__ = ["StandardAssessment"]


@register("l5", "standard")
class StandardAssessment(AssessmentBase):
    """Score a retrieval against truth and against the requirement set.

    Parameters
    ----------
    requirements
        Name of the requirement set, resolved by
        :func:`tss.config.resolve_path`.
    platform
        Optional platform state used for the along-track sampling metric.  The
        :class:`tss.interfaces.L5Assessment` protocol does not pass one to
        ``score``, so a campaign that wants SCI-09 scored supplies it here when
        the tier is built, and it is recorded in :meth:`describe`.
    extra_metrics
        Metrics computed elsewhere, merged in before scoring.  This is how an
        L2-owned number such as the achieved leakage isolation reaches the
        requirements check without this layer pretending to have measured it.
        A campaign that uses it should say so in its own record.
    """

    def __init__(
        self,
        *,
        requirements: str = "requirements/requirements",
        platform: PlatformState | None = None,
        extra_metrics: Mapping[str, Any] | None = None,
    ) -> None:
        super().__init__(requirements=requirements)
        self.platform = platform
        self.extra_metrics = dict(extra_metrics) if extra_metrics else {}

    def _settings(self) -> dict[str, Any]:
        return {
            "requirements": self.requirements_name,
            "has_platform": self.platform is not None,
            "extra_metrics": sorted(self.extra_metrics),
        }

    def score(
        self,
        retrieved: DensityProfile,
        truth: IonosphereState,
        ionogram: Ionogram,
        configuration: str,
        *,
        platform: PlatformState | None = None,
        truth_echo_snr: Mapping[str, Any] | None = None,
    ) -> Assessment:
        """Return a scored :class:`~tss.interfaces.Assessment`.  Nothing is mutated.

        ``platform``, when given, is used for the along-track metric in place of
        the one the tier was built with; the pipeline passes it.
        ``truth_echo_snr`` is the SNR L3 measured at every truth echo, which
        INS-10 is scored on; without it INS-10 is not evaluated.
        """
        if not configuration:
            raise ValueError(
                "configuration must be named.  Section 1.3 of the proposal reports "
                "per configuration and a result without its configuration is not a "
                "result."
            )

        metrics = compute_metrics(
            retrieved,
            truth,
            ionogram,
            platform=platform if platform is not None else self.platform,
            extra=self.extra_metrics,
            truth_echo_snr=truth_echo_snr,
        )
        records = check_requirements(metrics, self.requirements)

        errors = {
            key: metrics[key]
            for key in CARRY_IN_ERRORS
            if metrics.get(key) is not None
        }

        flags: dict[str, Any] = {
            "summary": summarise_verdicts(records),
            "requirements_revision": self.requirements.revision,
            "requirements_source": str(self.requirements.source),
            "n_requirements": len(self.requirements),
            # How the scored set sits inside the full catalogue, and every
            # catalogue requirement this run does not score, with the reason.
            "coverage": self.requirements.coverage(),
            "not_scored": self.requirements.not_scored_records(),
        }
        for prefix, source in (("l3", ionogram.flags), ("l4", retrieved.flags)):
            for key, value in (source or {}).items():
                flags[f"{prefix}.{key}"] = value

        return Assessment(
            configuration=str(configuration),
            requirements=records,
            errors=errors,
            flags=flags,
            source=f"tss.l5_assessment.standard.{type(self).__name__}",
        )
