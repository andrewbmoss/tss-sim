"""The requirement set, made machine readable.

The first governing principle of Section 3.1 is that no simulated performance
figure is scored pass or fail until the number it is scored against exists and
is traced to a source.  This module is that principle in code.  The requirement
set lives in ``config/requirements/requirements.yaml``, one entry per ID, each
carrying a threshold, a goal, a sense, a unit, the claim it serves, the work
package that validates it, and the source it is traced to.  L5 scores against
this file and nothing else.

A consequence worth stating: a requirement with ``source: derived`` and no
``derivation`` reference will fail :func:`validate`.  An untraceable number
cannot be used to score a result.

Keeping the YAML in step with ``TSS_Science_Requirements_RevE.docx`` and
``TSS_SRM_RevE.xlsx`` is a manual step for now.  ``scripts/import_srm.py``
exists to reduce that to a diff review rather than a retype.
"""

from __future__ import annotations

from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal

from tss.config import load_yaml, resolve_path

__all__ = ["Requirement", "RequirementSet", "load_requirements", "Verdict"]

Verdict = Literal["goal", "threshold", "fail", "not_evaluated"]

#: ``sense`` says which direction is better.  ``max`` means the achieved value
#: must be at or below the stated number (a limit on error, mass, power);
#: ``min`` means at or above (SNR, isolation, coverage); ``range`` means the
#: achieved value must lie inside ``[lower, upper]``; ``bool`` is a capability
#: that is either demonstrated or not.
SENSES = ("min", "max", "range", "bool")


@dataclass(frozen=True)
class Requirement:
    """One requirement with its threshold, goal and traceability."""

    id: str
    title: str
    sense: str
    unit: str = ""
    threshold: Any = None
    goal: Any = None
    claims: tuple[str, ...] = ()
    work_packages: tuple[str, ...] = ()
    source: str = ""
    derivation: str = ""
    notes: str = ""

    def evaluate(self, achieved: Any) -> dict[str, Any]:
        """Score an achieved value and return a record for the Assessment."""
        record: dict[str, Any] = {
            "id": self.id,
            "title": self.title,
            "achieved": achieved,
            "threshold": self.threshold,
            "goal": self.goal,
            "unit": self.unit,
            "sense": self.sense,
            "claims": list(self.claims),
            "source": self.source,
        }
        if achieved is None:
            record["verdict"] = "not_evaluated"
            record["margin"] = None
            return record

        verdict, margin = self._verdict(achieved)
        record["verdict"] = verdict
        record["margin"] = margin
        return record

    def _verdict(self, achieved: Any) -> tuple[Verdict, float | None]:
        if self.sense == "bool":
            met = bool(achieved)
            return ("goal" if met else "fail"), None

        if self.sense == "range":
            lo, hi = self.threshold
            inside = lo <= achieved <= hi
            if not inside:
                return "fail", float(min(achieved - lo, hi - achieved))
            if self.goal is not None:
                glo, ghi = self.goal
                if glo <= achieved <= ghi:
                    return "goal", float(min(achieved - glo, ghi - achieved))
            return "threshold", float(min(achieved - lo, hi - achieved))

        if self.sense == "min":
            if self.goal is not None and achieved >= self.goal:
                return "goal", float(achieved - self.goal)
            if self.threshold is not None and achieved >= self.threshold:
                return "threshold", float(achieved - self.threshold)
            base = self.threshold if self.threshold is not None else 0.0
            return "fail", float(achieved - base)

        if self.sense == "max":
            if self.goal is not None and achieved <= self.goal:
                return "goal", float(self.goal - achieved)
            if self.threshold is not None and achieved <= self.threshold:
                return "threshold", float(self.threshold - achieved)
            base = self.threshold if self.threshold is not None else 0.0
            return "fail", float(base - achieved)

        raise ValueError(f"{self.id}: unknown sense {self.sense!r}")


class RequirementSet(Mapping[str, Requirement]):
    """An immutable, indexable collection of requirements."""

    def __init__(self, requirements: Iterable[Requirement], *, revision: str = "",
                 source: Path | None = None) -> None:
        self._by_id = {r.id: r for r in requirements}
        self.revision = revision
        self.source = source

    def __getitem__(self, key: str) -> Requirement:
        return self._by_id[key]

    def __iter__(self):
        return iter(self._by_id)

    def __len__(self) -> int:
        return len(self._by_id)

    def for_claim(self, claim: str) -> list[Requirement]:
        """Every requirement that a given claim (C1..C5) is scored against."""
        return [r for r in self._by_id.values() if claim in r.claims]

    def for_work_package(self, wp: str) -> list[Requirement]:
        return [r for r in self._by_id.values() if wp in r.work_packages]

    def score(self, achieved: Mapping[str, Any]) -> dict[str, dict[str, Any]]:
        """Score a mapping of requirement ID to achieved value.

        Requirements absent from ``achieved`` are reported as
        ``not_evaluated`` rather than dropped, so the assessment always shows
        the full set and it is obvious what a campaign did not reach.
        """
        return {rid: req.evaluate(achieved.get(rid)) for rid, req in self._by_id.items()}

    def validate(self) -> list[str]:
        """Return a list of problems, empty if the set is well formed."""
        problems: list[str] = []
        for rid, req in self._by_id.items():
            if req.sense not in SENSES:
                problems.append(f"{rid}: sense {req.sense!r} not in {SENSES}")
            if req.sense != "bool" and req.threshold is None:
                problems.append(f"{rid}: no threshold")
            if not req.source:
                problems.append(f"{rid}: no source, so it cannot be used to score")
            if req.source == "derived" and not req.derivation:
                problems.append(f"{rid}: source is 'derived' but no derivation is named")
            if req.sense == "range":
                for field_name in ("threshold", "goal"):
                    val = getattr(req, field_name)
                    if val is not None and (not isinstance(val, (list, tuple)) or len(val) != 2):
                        problems.append(f"{rid}: {field_name} must be [lower, upper] for a range")
        return problems


def load_requirements(name: str | Path = "requirements/requirements") -> RequirementSet:
    """Load the requirement set from YAML."""
    path = resolve_path(name)
    data = load_yaml(path)
    entries = data.get("requirements", {})
    reqs = [
        Requirement(
            id=rid,
            title=body.get("title", ""),
            sense=body.get("sense", "max"),
            unit=body.get("unit", ""),
            threshold=body.get("threshold"),
            goal=body.get("goal"),
            claims=tuple(body.get("claims", ()) or ()),
            work_packages=tuple(body.get("work_packages", ()) or ()),
            source=body.get("source", ""),
            derivation=body.get("derivation", ""),
            notes=body.get("notes", ""),
        )
        for rid, body in entries.items()
    ]
    return RequirementSet(reqs, revision=data.get("revision", ""), source=path)
