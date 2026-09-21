"""Layer: L5
WP: WP0, WP7
Claims: C1, C2, C3, C4, C5
Status: implemented
Requirements: -

Shared helpers for the L5 assessment tiers.

There is no scoring arithmetic here.  The error metrics are in
:mod:`tss.l5_assessment.scoring` and the requirement mapping is in
:mod:`tss.l5_assessment.requirements_check`.  What lives here is the loading
and caching of the requirement set, the description machinery, and the summary
counting that both the assessment tier and the report renderer need.

One convention, stated once.  A tier holds the name of the requirement set it
scores against rather than a loaded object, and loads it on first use.  The
name goes into :meth:`AssessmentBase.describe` and therefore into the run
manifest, so a stored assessment can always be traced back to the exact
requirement file it was scored against, which is what
``docs/requirements-traceability.md`` asks for.  Holding a loaded object in the
tier instead would put an unserialisable thing in the manifest and would make
the revision invisible.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Mapping
from typing import Any

from tss.requirements import RequirementSet, load_requirements

__all__ = ["AssessmentBase", "VERDICTS", "summarise_verdicts"]

#: The four verdicts, in the order a report shows them.  ``goal`` and
#: ``threshold`` are both met; the distinction is which level.
VERDICTS: tuple[str, ...] = ("goal", "threshold", "fail", "not_evaluated")


def summarise_verdicts(records: Mapping[str, Mapping[str, Any]]) -> dict[str, int]:
    """Count verdicts, with every verdict present even when its count is zero.

    Zero counts are kept so that a summary line has the same shape every time
    and a reader can see that nothing failed rather than having to notice that
    ``fail`` is missing.
    """
    counts = Counter(str(r.get("verdict", "not_evaluated")) for r in records.values())
    return {verdict: int(counts.get(verdict, 0)) for verdict in VERDICTS}


class AssessmentBase:
    """Common base for L5 tiers.

    A subclass registers itself with :func:`tss.registry.register`, which sets
    ``name`` and ``layer`` on the class, and implements ``score`` to satisfy
    :class:`tss.interfaces.L5Assessment`.

    Parameters
    ----------
    requirements
        Name of the requirement set to score against, resolved by
        :func:`tss.config.resolve_path`.  The default is the set the campaign
        configuration also defaults to.
    """

    name: str = "base"
    layer: str = "l5"

    def __init__(self, *, requirements: str = "requirements/requirements") -> None:
        self.requirements_name = str(requirements)
        self._requirements: RequirementSet | None = None

    @property
    def requirements(self) -> RequirementSet:
        """The requirement set, loaded on first use and cached thereafter."""
        if self._requirements is None:
            self._requirements = load_requirements(self.requirements_name)
        return self._requirements

    def _settings(self) -> dict[str, Any]:
        return {"requirements": self.requirements_name}

    def describe(self) -> Mapping[str, Any]:
        """JSON-serialisable identity and settings, for the run manifest."""
        record: dict[str, Any] = {
            "layer": "l5",
            "name": getattr(self, "name", type(self).__name__),
            "class": type(self).__name__,
            "module": type(self).__module__,
        }
        record.update(self._settings())
        try:
            record["requirements_revision"] = self.requirements.revision
            record["n_requirements"] = len(self.requirements)
        except Exception as exc:  # pragma: no cover - only when the file is missing
            record["requirements_revision"] = f"unavailable: {exc}"
        return record
