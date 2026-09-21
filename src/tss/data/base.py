"""Layer: data
WP: WP1, WP2, WP8, WP12
Claims: -
Status: implemented
Requirements: -

The :class:`Source` descriptor and the registry that holds them.

A source is an external data set the simulator can consume.  The descriptor
carries what a reader of the dissertation needs in order to check the work:
where the data came from, what licence it is under, how it should be cited, how
often it changes, and which work packages depend on it.  None of that is
optional.  A source whose licence or citation is unknown is a source that
cannot appear in the dissertation, so the constructor refuses to build one.

Registration mirrors :mod:`tss.registry`.  A fetcher module decorates its
fetch function and the decorator records the descriptor::

    @register_source(
        name="spdf_isis_ne",
        description="ISIS-1 hand-scaled topside Ne profiles",
        url="https://spdf.gsfc.nasa.gov/pub/data/isis/",
        licence="Public domain (US Government work)",
        citation="...",
        cadence="static archive",
        work_packages=("WP1", "WP8"),
    )
    def fetch_isis_ne(*, force: bool = False, limit: int | None = None):
        ...

Every fetcher takes ``force`` and ``limit`` as keyword-only arguments and
returns the list of local paths it produced, newest first is not required but
the order must be stable.  ``force`` means ignore the cache.  ``limit`` means
stop after that many files or records, and a fetcher that cannot honour a limit
should say so in its docstring rather than silently ignore it.
"""

from __future__ import annotations

from collections.abc import Callable, Mapping, Sequence
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

__all__ = [
    "Source",
    "SOURCES",
    "register_source",
    "get_source",
    "UnknownSource",
    "FetchError",
]


FetchFn = Callable[..., "list[Path]"]


class UnknownSource(KeyError):
    """Raised when a name does not match any registered source."""


class FetchError(RuntimeError):
    """Raised when a fetch fails for a reason the caller should see.

    This is deliberately distinct from the transport library's own exception
    types so that callers do not have to import requests to catch a failure.
    """


@dataclass(frozen=True)
class Source:
    """One external data set, with the provenance needed to cite it.

    Attributes
    ----------
    name
        Registry key.  Lowercase with underscores.  This is what ``tss fetch``
        takes on the command line.
    description
        One line, shown by ``tss fetch list``.
    url
        The landing page or archive root a human should visit.  Not necessarily
        the file URL, which several of these sources only expose through an
        API.
    licence
        The licence or rights statement, in the words the provider uses.  "Public
        domain (US Government work)" is a licence for this purpose.  Empty is
        not allowed.
    citation
        How the source asks to be cited, or the best available reference.  Empty
        is not allowed.  This string is copied into the run manifest.
    cadence
        How often the upstream data changes, which is what decides a sensible
        cache lifetime.  Free text such as "static archive", "daily", "live".
    work_packages
        The work packages that consume this source.  Used to answer "what breaks
        if this archive goes away".
    fetch
        The callable that downloads or locates the data.
    acknowledgement
        Any acknowledgement text the provider asks users to carry, separate from
        the citation.  Several US Government archives are public domain but ask
        for acknowledgement anyway, and honouring that is free.
    notes
        Candid remarks about limitations, rate limits, or known gaps.
    """

    name: str
    description: str
    url: str
    licence: str
    citation: str
    cadence: str
    work_packages: Sequence[str]
    fetch: FetchFn
    acknowledgement: str = ""
    notes: str = ""
    extra: Mapping[str, Any] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if not self.name or not self.name.strip():
            raise ValueError("a source needs a name")
        if not self.licence.strip():
            raise ValueError(
                f"source {self.name!r} has no licence; a data set whose rights are "
                "unknown cannot be used in the dissertation"
            )
        if not self.citation.strip():
            raise ValueError(
                f"source {self.name!r} has no citation; every external input is "
                "cited in the manifest, so an uncitable source is not usable"
            )
        if not callable(self.fetch):
            raise TypeError(f"source {self.name!r} has a non-callable fetch")

    def describe(self) -> dict[str, Any]:
        """A JSON-serialisable record of this source, for the run manifest."""
        return {
            "name": self.name,
            "description": self.description,
            "url": self.url,
            "licence": self.licence,
            "citation": self.citation,
            "cadence": self.cadence,
            "work_packages": list(self.work_packages),
            "acknowledgement": self.acknowledgement,
            "notes": self.notes,
            "extra": dict(self.extra),
        }


SOURCES: dict[str, Source] = {}


def register_source(
    *,
    name: str,
    description: str,
    url: str,
    licence: str,
    citation: str,
    cadence: str,
    work_packages: Sequence[str],
    acknowledgement: str = "",
    notes: str = "",
    extra: Mapping[str, Any] | None = None,
) -> Callable[[FetchFn], FetchFn]:
    """Decorator registering a fetch function as a named source.

    The decorated function is returned unchanged so that it stays directly
    importable and directly testable.  Re-registering the same name with the
    same function is a no-op, which keeps a double import from being an error;
    re-registering it with a different function raises.
    """

    def _decorate(fn: FetchFn) -> FetchFn:
        source = Source(
            name=name,
            description=description,
            url=url,
            licence=licence,
            citation=citation,
            cadence=cadence,
            work_packages=tuple(work_packages),
            fetch=fn,
            acknowledgement=acknowledgement,
            notes=notes,
            extra=dict(extra or {}),
        )
        existing = SOURCES.get(name)
        if existing is not None and existing.fetch is not fn:
            raise ValueError(
                f"source {name!r} is already registered to "
                f"{getattr(existing.fetch, '__qualname__', existing.fetch)!r}"
            )
        SOURCES[name] = source
        return fn

    return _decorate


def get_source(name: str) -> Source:
    """Look up a source, with an error that lists what is registered."""
    try:
        return SOURCES[name]
    except KeyError as exc:
        known = ", ".join(sorted(SOURCES)) or "nothing"
        raise UnknownSource(
            f"no data source named {name!r}; registered: {known}"
        ) from exc
