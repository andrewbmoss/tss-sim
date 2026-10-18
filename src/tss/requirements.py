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

The catalogue and the scored set
--------------------------------
The SRM workbook (``TSS_SRM_RevG.xlsx``) is the requirement set of record.
``scripts/import_srm.py`` transcribes every requirement in it into
``config/requirements/srm_catalog.yaml``, a generated file that carries the
workbook's wording and is never edited by hand.  ``requirements.yaml`` names
that catalogue under ``catalog:`` and then does two things with it:

* it carries, as scored entries, the requirements a simulation run can score,
  with the numbers, senses and units the catalogue states in prose; and
* it lists every other catalogue ID under ``not_scored``, each with a category
  from :data:`NOT_SCORED_CATEGORIES` and a reason.

:meth:`RequirementSet.validate` fails if a catalogue ID is in neither place, in
both, or if a scored ID is missing from the catalogue.  So the scored subset can
be smaller than the full set, and it is, but it cannot shrink silently: a
requirement that is not scored says why, and the L5 report prints it.

Taking in a new SRM revision is then a diff review rather than a retype: rerun
``import_srm.py``, read the diff of the catalogue, and amend the scored entries
whose text changed.
"""

from __future__ import annotations

from collections import Counter
from collections.abc import Iterable, Mapping
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Literal

from tss.config import load_yaml, resolve_path

__all__ = [
    "Requirement",
    "RequirementSet",
    "NotScored",
    "NOT_SCORED_CATEGORIES",
    "load_requirements",
    "Verdict",
]

Verdict = Literal["goal", "threshold", "fail", "not_evaluated"]

#: ``sense`` says which direction is better.  ``max`` means the achieved value
#: must be at or below the stated number (a limit on error, mass, power);
#: ``min`` means at or above (SNR, isolation, coverage); ``range`` means the
#: achieved value must lie inside ``[lower, upper]``; ``bool`` is a capability
#: that is either demonstrated or not.
SENSES = ("min", "max", "range", "bool")

#: Why a catalogue requirement is not scored by a simulation run.
#:
#: ``flowed_down``  scored through the derived requirements it names in
#:                  ``scored_through``, so scoring it again would count one
#:                  measurement twice.
#: ``pending``      scorable in simulation, but no tier produces its metric
#:                  yet; ``work_packages`` names who will.
#: ``analysis``     verified by analysis, inspection or a stand-alone study,
#:                  not by a per-run score.
#: ``bench``        verified on hardware or in orbit.
#: ``mission``      a mission-level requirement, verified by orbit design and
#:                  propagation rather than by an instrument run.
#: ``software``     a requirement on tss-sim itself, verified by its tests.
#: ``tbd``          the value is not yet set.
NOT_SCORED_CATEGORIES: tuple[str, ...] = (
    "flowed_down",
    "pending",
    "analysis",
    "bench",
    "mission",
    "software",
    "tbd",
)

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
                v for f, v in zip(freq, vals, strict=True)
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


@dataclass(frozen=True)
class NotScored:
    """A catalogue requirement that a simulation run does not score, and why."""

    id: str
    category: str
    reason: str
    scored_through: tuple[str, ...] = ()
    work_packages: tuple[str, ...] = ()

    def as_record(self, catalog_entry: Mapping[str, Any] | None = None) -> dict[str, Any]:
        """A JSON-serialisable record for the assessment and the report."""
        entry = catalog_entry or {}
        return {
            "id": self.id,
            "category": self.category,
            "reason": self.reason,
            "scored_through": list(self.scored_through),
            "work_packages": list(self.work_packages),
            "requirement": str(entry.get("requirement", "")),
            "verification_methods": list(entry.get("verification_methods", []) or []),
        }


class RequirementSet(Mapping[str, Requirement]):
    """An immutable, indexable collection of scored requirements.

    Iterating, indexing and :meth:`score` cover the scored entries only.  The
    catalogue the set belongs to, and the entries it deliberately does not
    score, are carried beside them in :attr:`catalog` and :attr:`not_scored`.
    """

    def __init__(
        self,
        requirements: Iterable[Requirement],
        *,
        revision: str = "",
        source: Path | None = None,
        catalog: Mapping[str, Mapping[str, Any]] | None = None,
        catalog_info: Mapping[str, Any] | None = None,
        not_scored: Iterable[NotScored] = (),
    ) -> None:
        self._by_id = {r.id: r for r in requirements}
        self.revision = revision
        self.source = source
        self.catalog: dict[str, Mapping[str, Any]] = dict(catalog or {})
        self.catalog_info: dict[str, Any] = dict(catalog_info or {})
        self.not_scored: dict[str, NotScored] = {n.id: n for n in not_scored}

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

    def not_scored_records(self) -> dict[str, dict[str, Any]]:
        """Every deliberately unscored catalogue requirement, in catalogue order."""
        order = list(self.catalog) or sorted(self.not_scored)
        return {
            rid: self.not_scored[rid].as_record(self.catalog.get(rid))
            for rid in order
            if rid in self.not_scored
        }

    def coverage(self) -> dict[str, Any]:
        """How much of the catalogue the scored set covers, and how the rest is held."""
        by_category = Counter(n.category for n in self.not_scored.values())
        unaccounted = sorted(set(self.catalog) - set(self._by_id) - set(self.not_scored))
        return {
            "catalog_revision": str(self.catalog_info.get("revision", "")),
            "catalog_workbook": str(self.catalog_info.get("source_workbook", "")),
            "catalog_sha256": str(self.catalog_info.get("source_sha256", "")),
            "n_catalog": len(self.catalog),
            "n_scored": len(self._by_id),
            "n_not_scored": len(self.not_scored),
            "not_scored_by_category": {
                cat: int(by_category.get(cat, 0)) for cat in NOT_SCORED_CATEGORIES
            },
            "unaccounted": unaccounted,
        }

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

        for rid, entry in self.not_scored.items():
            if entry.category not in NOT_SCORED_CATEGORIES:
                problems.append(
                    f"{rid}: not_scored category {entry.category!r} not in {NOT_SCORED_CATEGORIES}"
                )
            if not entry.reason.strip():
                problems.append(f"{rid}: not_scored without a reason")
            if rid in self._by_id:
                problems.append(f"{rid}: both scored and listed under not_scored")
            if entry.category == "flowed_down":
                if not entry.scored_through:
                    problems.append(f"{rid}: flowed_down must name what it is scored through")
                for target in entry.scored_through:
                    if target not in self._by_id:
                        problems.append(
                            f"{rid}: scored_through {target}, which is not a scored requirement"
                        )
            if entry.category == "pending" and not entry.work_packages:
                problems.append(f"{rid}: pending must name the work package that will score it")

        if self.catalog:
            for rid in self._by_id:
                if rid not in self.catalog:
                    problems.append(f"{rid}: scored but not in the catalogue")
            for rid in self.not_scored:
                if rid not in self.catalog:
                    problems.append(f"{rid}: listed under not_scored but not in the catalogue")
            for rid in self.coverage()["unaccounted"]:
                problems.append(
                    f"{rid}: in the catalogue but neither scored nor listed under not_scored"
                )
        elif self.not_scored:
            problems.append("not_scored is given but no catalogue is named")
        return problems


def _not_scored_from(block: Mapping[str, Any]) -> list[NotScored]:
    out: list[NotScored] = []
    for rid, body in (block or {}).items():
        body = body or {}
        out.append(
            NotScored(
                id=str(rid),
                category=str(body.get("category", "")),
                reason=str(body.get("reason", "")),
                scored_through=tuple(body.get("scored_through", ()) or ()),
                work_packages=tuple(body.get("work_packages", ()) or ()),
            )
        )
    return out


def load_requirements(name: str | Path = "requirements/requirements") -> RequirementSet:
    """Load the requirement set from YAML, with its catalogue when one is named."""
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

    catalog: dict[str, Mapping[str, Any]] = {}
    catalog_info: dict[str, Any] = {}
    if data.get("catalog"):
        doc = load_yaml(resolve_path(str(data["catalog"])))
        catalog = dict(doc.get("requirements", {}) or {})
        catalog_info = {k: v for k, v in doc.items() if k != "requirements"}

    return RequirementSet(
        reqs,
        revision=data.get("revision", ""),
        source=path,
        catalog=catalog,
        catalog_info=catalog_info,
        not_scored=_not_scored_from(data.get("not_scored", {})),
    )
