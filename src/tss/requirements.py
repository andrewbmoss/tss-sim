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

_RANK = {"fail": 0, "threshold": 1, "goal": 2}


def _band_label(band: Mapping[str, Any]) -> str:
    lo = band.get("f_min_hz")
    hi = band.get("f_max_hz")
    lo_s = "0" if lo is None else f"{float(lo) / 1e6:g}"
    hi_s = "inf" if hi is None else f"{float(hi) / 1e6:g}"
    return f"{lo_s}-{hi_s} MHz"


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
    #: Frequency bands with their own threshold and goal, for a requirement
    #: stated per band (INS-05 in Rev G, amended 2026-09-30).  Each band is a
    #: mapping with ``f_min_hz`` (inclusive, absent = no lower edge),
    #: ``f_max_hz`` (exclusive, absent = no upper edge), ``threshold`` and
    #: ``goal``.  The top-level ``threshold`` and ``goal`` are then the
    #: strictest band's, used when an achieved value carries no frequencies.
    bands: tuple[Mapping[str, Any], ...] = ()

    def _band_values(self, achieved: Mapping[str, Any]):
        freq = [float(x) for x in achieved.get("freq_hz", [])]
        vals = [float(x) for x in achieved.get("value", [])]
        if len(freq) != len(vals):
            raise ValueError(f"{self.id}: freq_hz and value must be parallel")
        for band in self.bands:
            lo = band.get("f_min_hz")
            hi = band.get("f_max_hz")
            sel = [
                v for f, v in zip(freq, vals)
                if (lo is None or f >= float(lo)) and (hi is None or f < float(hi))
            ]
            yield band, sel

    def _evaluate_banded(self, achieved: Mapping[str, Any]) -> dict[str, Any]:
        """Score a per-frequency achieved value band by band.

        Each band is scored on its worst value (the largest for sense max, the
        smallest for sense min).  The requirement's verdict is the worst band
        verdict, and the record's ``achieved``, ``threshold`` and ``goal`` are
        those of the governing band, the one with that verdict and the least
        relative margin, so the one-line report stays readable.
        """
        per_band: list[dict[str, Any]] = []
        for band, values in self._band_values(achieved):
            if not values:
                continue
            worst = max(values) if self.sense == "max" else min(values)
            sub = Requirement(
                id=self.id, title=self.title, sense=self.sense, unit=self.unit,
                threshold=band.get("threshold"), goal=band.get("goal"),
            )
            verdict, margin = sub._verdict(worst)
            ref = band.get("threshold") or 1.0
            per_band.append({
                "band": _band_label(band),
                "achieved": worst,
                "threshold": band.get("threshold"),
                "goal": band.get("goal"),
                "verdict": verdict,
                "margin": margin,
                "relative_margin": None if margin is None else margin / abs(float(ref)),
                "n_values": len(values),
            })
        record: dict[str, Any] = {
            "id": self.id, "title": self.title, "unit": self.unit, "sense": self.sense,
            "claims": list(self.claims), "source": self.source, "bands": per_band,
        }
        if not per_band:
            record.update(achieved=None, threshold=self.threshold, goal=self.goal,
                          verdict="not_evaluated", margin=None)
            return record
        gov = min(per_band, key=lambda b: (_RANK[b["verdict"]], b["relative_margin"]))
        record.update(achieved=gov["achieved"], threshold=gov["threshold"], goal=gov["goal"],
                      verdict=gov["verdict"], margin=gov["margin"], band=gov["band"])
        return record

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

        if self.bands and isinstance(achieved, Mapping):
            return self._evaluate_banded(achieved)

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
            for band in req.bands:
                if not isinstance(band, Mapping) or band.get("threshold") is None:
                    problems.append(f"{rid}: every band needs a threshold")
                    continue
                lo, hi = band.get("f_min_hz"), band.get("f_max_hz")
                if lo is not None and hi is not None and float(hi) <= float(lo):
                    problems.append(f"{rid}: band {_band_label(band)} is empty")
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
            bands=tuple(body.get("bands", ()) or ()),
        )
        for rid, body in entries.items()
    ]
    return RequirementSet(reqs, revision=data.get("revision", ""), source=path)
