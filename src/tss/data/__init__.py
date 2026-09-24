"""Layer: data
WP: WP1, WP2, WP8, WP12
Claims: -
Status: implemented
Requirements: INS-01

External data acquisition for the simulator.

Four work packages need something from outside the package and each of them
needs it for a different reason.

WP1 needs ionosphere truth to check L0 against.  The Alouette and ISIS topside
sounders are the only substantial body of measured topside electron density
profiles that exists, and the hand-scaled CRC profiles from ISIS-1 give the
local plasma frequency at the spacecraft directly.  That is the quantity the
lower band edge of INS-01 is built from, so the same file serves as truth for
WP1 and as the measurement INS-01 is argued against.  See :mod:`tss.data.spdf_isis`.

WP2 needs an absorption reference.  The data.gov collection "Measurements of
the absorptive properties of the ionosphere" is the reference the integrated
absorption model is scored against, and the WSPR spot record is an independent
cross-check on its diurnal and frequency behaviour.  See
:mod:`tss.data.datagov_absorption` and :mod:`tss.data.wspr`.

WP8 and V4 need archived real ionograms, which live in the same SPDF archive as
the profiles but as a much larger and less uniform product.  This layer records
where they are and does not bulk-download them.

WP1 and WP8 also need present-day truth that the ISIS archive cannot give.
Madrigal supplies incoherent scatter radar profiles through the F2 peak and into
the topside, and DMSP in-situ density at about 850 km; GIRO DIDBase supplies
the F2 peak from ground Digisondes; IRTAM supplies assimilated global peak maps
for event-driven L0 truth.  See :mod:`tss.data.madrigal`, :mod:`tss.data.giro`
and :mod:`tss.data.irtam`.

WP12 needs a sub-megahertz sky noise reference.  The WSPR record does not reach
below 137 kHz and is a poor proxy there, so what this layer provides for WP12
is the band behaviour above 1.8 MHz as context; the sub-megahertz reference
itself comes from the Franklin and MacLean tabulation handled elsewhere.  That
gap is real and is stated here rather than papered over.

Interface.  :data:`SOURCES` maps a name to a :class:`~tss.data.base.Source`
descriptor carrying licence, citation, cadence and the work packages that
depend on it.  :func:`fetch` runs one source's fetcher and returns the local
paths.  ``tss fetch list`` prints the registry and ``tss fetch <name>`` calls
:func:`fetch`.

Import cost.  Importing this package must stay cheap, because
:mod:`tss.cli` imports it to list sources.  Nothing here imports requests,
pandas or h5py at module scope; every one of those is imported inside the
function that needs it.  A machine with none of them installed can still run
``tss fetch list`` and can still import every module in this package.

Offline behaviour.  Nothing in this package performs a network request at
import time, and nothing does so unless a fetcher or a query function is
called.  The test suite never calls either.
"""

from __future__ import annotations

from pathlib import Path
from typing import Any

# Importing the fetcher modules is what populates SOURCES, exactly as importing
# the layer packages is what populates the tier registry.  These modules are
# import-cheap by construction; see the note on import cost above.
from tss.data import (  # noqa: E402,F401
    cache,
    datagov_absorption,
    giro,
    irtam,
    local_models,
    madrigal,
    spdf_isis,
    wspr,
)
from tss.data.base import (
    SOURCES,
    FetchError,
    Source,
    UnknownSource,
    get_source,
    register_source,
)

__all__ = [
    "SOURCES",
    "Source",
    "FetchError",
    "UnknownSource",
    "fetch",
    "get_source",
    "register_source",
    "describe_sources",
    "cache",
    "giro",
    "irtam",
    "local_models",
    "madrigal",
]


def fetch(name: str, *, force: bool = False, limit: int | None = None, **kwargs: Any) -> list[Path]:
    """Fetch one named source and return the local paths it produced.

    ``force`` ignores the cache and re-downloads.  ``limit`` caps how much is
    fetched, where a source can honour it; a source that cannot says so in its
    own docstring and ignores it rather than pretending.  Extra keyword
    arguments are passed through to the source's fetcher, which is how
    :func:`tss.data.spdf_isis.fetch_isis_ne` takes a ``spacecraft`` argument
    from a caller that wants one.

    Raises :class:`~tss.data.base.UnknownSource` for a name that is not
    registered, and :class:`~tss.data.base.FetchError` for a fetch that failed
    for a reason the caller should see.

    This function makes network requests.  It is never called from the test
    suite.
    """
    source = get_source(name)
    paths = source.fetch(force=force, limit=limit, **kwargs)
    return [Path(p) for p in paths]


def describe_sources() -> dict[str, dict[str, Any]]:
    """Every source descriptor as JSON, for the run manifest.

    A run that consumed external data records what it consumed and under what
    terms.  Copying the whole registry rather than only the sources used is
    deliberate and cheap: it means a manifest also records what was available
    and not taken.
    """
    return {name: src.describe() for name, src in sorted(SOURCES.items())}
