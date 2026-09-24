"""Layer: data
WP: WP1
Claims: -
Status: implemented (fetch and parse; evaluation in tss.l0_ionosphere.irtam_maps)
Requirements: -

IRTAM, the IRI-based Real-Time Assimilative Model, as coefficient files.

What it provides.  IRTAM, run at the Lowell GIRO Data Center, assimilates the
GIRO Digisonde network into IRI and publishes the result as a set of
coefficients every 15 minutes for foF2, hmF2, B0 and B1.  The coefficients
replace IRI's climatological CCIR/URSI maps for that instant, so evaluating them
gives a global, data-driven map of the F2 peak and bottomside shape for a real
day.

What it is good for here.  Event-driven L0 truth.  IRI gives a monthly-median
ionosphere; IRTAM gives the ionosphere that actually existed at 12:00 UT on a
given day, with the peak where the Digisondes saw it.  That is the input for
simulating a pass over a real storm or a real equatorial anomaly, and it is the
right reference for any conjunction study against GIRO, because the two share
their underlying measurements.  It is not independent of GIRO, and it must not
be used to validate against GIRO.

Two routes to the same files.

* The static archive at ``ulcar.uml.edu/GAMBIT/GambitCoefficients/`` holds
  files named ``IRTAM_<char>_COEFFS_<YYYYMMDD>_<HHMMSS>.ASC`` in
  ``COEFFS_<YYYY>/<YYYY_MM_DD>/`` directories, at 15-minute cadence, for 2000
  through 2025 as of 2026-09-23.  B0 and B1 are named ``B0in`` and ``B1in`` there.
* The servlet at ``lgdc.uml.edu/rix/gambit-coeffs?charName=<char>&time=<ISO>``
  computes the same file on request and is the route for times the archive
  does not yet hold.  It serves foF2, hmF2, B0 and B1 to a guest, as checked
  on 2026-09-23, although the web form in front of it offers only three.  It
  covered every 2026 epoch tried through 2026-09-21, answered HTTP 500 for
  2026-09-22 and for "now", and answered 429 once probed a dozen times in a few
  minutes.  Treat the last day or two as not yet available, and keep requests
  sparse.

:func:`fetch_coefficients` tries the archive first, since static files cost the
provider least, and falls back to the servlet on a 404.

File format.  A ``#`` header between ``START_HEADER`` and ``END_HEADER`` with the
characteristic and its unit, the time of validity, the assimilative engine, the
expansion basis and its lengths, and the list of assimilated stations; then
the coefficients as whitespace-separated floats, four per line.  The basis is
stated as ``JonesGallet_LinTrend`` with lengths ``14(temporal) x 76(spatial)``,
so 1064 numbers.

Layout of the 1064 numbers.  Not documented in the file.  The layout used
here is the one in PyIRTAM (Forsythe et al. 2024, Space Weather,
doi:10.1029/2024SW003965), the published Python evaluator from the Naval
Research Laboratory, which reads the file as follows.  The first 988 numbers
are a CCIR-style block of 13 temporal by 76 spatial coefficients stored in
Fortran order, so the temporal index varies fastest.  The last 76 numbers are
one additional temporal term, the linear trend the header's ``_LinTrend``
refers to, which PyIRTAM inserts as temporal row 1.  The resulting
``(14, 76)`` matrix therefore has row order: mean, trend, then the twelve
diurnal harmonic terms.  :attr:`IrtamCoefficients.matrix` returns exactly
that.  It was checked by reading PyIRTAM 0.0.7's ``IRTAM_read_files``, not by
evaluating a map against IRI.

Evaluation.  This module fetches and parses.  Turning coefficients into maps
of foF2, hmF2, B0 and B1 is :mod:`tss.l0_ionosphere.irtam_maps`, which
reproduces PyIRTAM to rounding error; see ``docs/irtam-validation.md``.

Privacy.  The servlet writes the requester's IP address into the header.  The
cached file keeps it; fixtures redact it.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from datetime import date, datetime, timedelta, timezone
from pathlib import Path
from typing import Any

from tss.data import cache
from tss.data.base import FetchError, register_source

__all__ = [
    "SOURCE",
    "ARCHIVE_ROOT",
    "SERVLET",
    "CHARACTERISTICS",
    "CADENCE_MIN",
    "IrtamCoefficients",
    "archive_url",
    "servlet_url",
    "round_to_cadence",
    "parse_coefficients",
    "fetch_coefficients",
    "fetch_irtam_sample",
]

SOURCE = "irtam"
ARCHIVE_ROOT = "https://ulcar.uml.edu/GAMBIT/GambitCoefficients/"
SERVLET = "https://lgdc.uml.edu/rix/gambit-coeffs"
CADENCE_MIN = 15
MIN_INTERVAL_S = 5.0
"""Minimum gap between requests, seconds.  The servlet answered 429 at about three
seconds spacing on 2026-09-23."""

RECENT_LAG_DAYS = 3
"""How far back the default sample reaches.  The servlet could not compute the
most recent day on 2026-09-23 (HTTP 500), while three days back was reliable."""

CHARACTERISTICS: dict[str, dict[str, str | None]] = {
    "foF2": {"archive": "foF2", "servlet": "foF2", "unit": "MHz"},
    "hmF2": {"archive": "hmF2", "servlet": "hmF2", "unit": "km"},
    "B0": {"archive": "B0in", "servlet": "B0", "unit": "km"},
    "B1": {"archive": "B1in", "servlet": "B1", "unit": ""},
}
"""Characteristic to its archive file name, its servlet name, and its unit.

``None`` for the servlet means the servlet does not offer it to a guest.
"""


@dataclass(frozen=True)
class IrtamCoefficients:
    """One IRTAM coefficient file, parsed.

    ``coefficients`` is the flat vector exactly as the file lists it.
    ``n_temporal`` and ``n_spatial`` are the basis lengths the header declares,
    and :attr:`matrix` reshapes the vector under the inferred layout described
    in the module docstring.
    """

    characteristic: str
    unit: str
    valid_time: datetime
    engine: str
    basis: str
    n_temporal: int
    n_spatial: int
    grid: str
    stations: tuple[str, ...]
    generated_by: str
    coefficients: Any
    header: tuple[str, ...] = ()

    @property
    def matrix(self) -> Any:
        """``(n_temporal, n_spatial)`` coefficient matrix in PyIRTAM's layout.

        Row 0 is the mean term, row 1 the linear-trend term (the last
        ``n_spatial`` numbers in the file), rows 2 onward the diurnal harmonic
        terms.  Columns are the Jones and Gallet spatial functions.  See the
        module docstring for where this layout comes from.
        """
        import numpy as np

        nt, ns = self.n_temporal, self.n_spatial
        main = self.coefficients[: (nt - 1) * ns].reshape((nt - 1, ns), order="F")
        trend = self.coefficients[(nt - 1) * ns :]
        out = np.empty((nt, ns), dtype=float)
        out[0] = main[0]
        out[1] = trend
        out[2:] = main[1:]
        return out


# ---------------------------------------------------------------------------
# URLs
# ---------------------------------------------------------------------------


def _as_utc(value: str | date | datetime) -> datetime:
    if isinstance(value, datetime):
        return value.astimezone(timezone.utc) if value.tzinfo else value.replace(tzinfo=timezone.utc)
    if isinstance(value, date):
        return datetime(value.year, value.month, value.day, tzinfo=timezone.utc)
    if isinstance(value, str):
        dt = datetime.fromisoformat(value.strip().replace("Z", "+00:00"))
        return dt.astimezone(timezone.utc) if dt.tzinfo else dt.replace(tzinfo=timezone.utc)
    raise TypeError(f"cannot read a {type(value).__name__} as a time")


def round_to_cadence(t: str | date | datetime) -> datetime:
    """Round a time to the nearest 15-minute IRTAM epoch, in UTC.

    Halfway cases round up.  The caller is told which epoch was used through
    :attr:`IrtamCoefficients.valid_time`, so rounding is never silent.
    """
    dt = _as_utc(t)
    minutes = dt.hour * 60 + dt.minute + (dt.second + dt.microsecond / 1e6) / 60.0
    epoch = int((minutes + CADENCE_MIN / 2.0) // CADENCE_MIN) * CADENCE_MIN
    base = datetime(dt.year, dt.month, dt.day, tzinfo=timezone.utc)
    return base + timedelta(minutes=epoch)


def _check_char(characteristic: str) -> dict[str, str | None]:
    try:
        return CHARACTERISTICS[characteristic]
    except KeyError as exc:
        raise ValueError(
            f"unknown IRTAM characteristic {characteristic!r}; known: "
            f"{', '.join(CHARACTERISTICS)}"
        ) from exc


def archive_url(characteristic: str, t: str | date | datetime) -> str:
    """The static archive URL for one characteristic at one epoch."""
    name = _check_char(characteristic)["archive"]
    e = round_to_cadence(t)
    return (
        f"{ARCHIVE_ROOT}COEFFS_{e:%Y}/{e:%Y_%m_%d}/"
        f"IRTAM_{name}_COEFFS_{e:%Y%m%d}_{e:%H%M%S}.ASC"
    )


def servlet_url(characteristic: str, t: str | date | datetime) -> str:
    """The servlet URL for one characteristic at one epoch.

    Raises ``ValueError`` for a characteristic the servlet does not serve.
    """
    name = _check_char(characteristic)["servlet"]
    if name is None:
        raise ValueError(f"the IRTAM servlet does not serve {characteristic}; use the archive")
    e = round_to_cadence(t)
    return f"{SERVLET}?charName={name}&time={e:%Y-%m-%dT%H:%M}"


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

_CHAR = re.compile(r"Ionospheric Characteristic:\s*(\S+)\s*(?:\[([^\]]*)\])?")
_LENGTHS = re.compile(r"Basis Lengths:\s*(\d+)\s*\(temporal\)\s*x\s*(\d+)\s*\(spatial\)")


def parse_coefficients(text: str | bytes) -> IrtamCoefficients:
    """Parse an IRTAM coefficient file from either route.

    Checks that the number of coefficients equals the product of the declared
    basis lengths, because a truncated download is otherwise indistinguishable
    from a valid file.
    """
    import numpy as np

    if isinstance(text, bytes):
        text = text.decode("utf-8", "replace")
    lines = text.splitlines()
    if not any("START_HEADER" in ln for ln in lines[:3]):
        snippet = " ".join(text.split())[:160]
        raise FetchError(f"not an IRTAM coefficient file: {snippet!r}")

    header: list[str] = []
    numbers: list[str] = []
    in_header = True
    for ln in lines:
        s = ln.strip()
        if not s:
            continue
        if in_header:
            header.append(s)
            if "END_HEADER" in s:
                in_header = False
            continue
        if s.startswith("#"):
            continue
        numbers.extend(s.split())

    joined = "\n".join(header)
    m = _CHAR.search(joined)
    if not m:
        raise FetchError("IRTAM header has no 'Ionospheric Characteristic' line")
    char, unit = m.group(1), (m.group(2) or "").strip()

    def field(label: str) -> str:
        for h in header:
            body = h.lstrip("#").strip()
            if body.startswith(label):
                return body[len(label):].strip(" :")
        return ""

    valid = field("Time of Validity")
    try:
        valid_time = _as_utc(valid)
    except ValueError as exc:
        raise FetchError(f"cannot read IRTAM time of validity {valid!r}") from exc

    lm = _LENGTHS.search(joined)
    if not lm:
        raise FetchError("IRTAM header has no 'Basis Lengths' line")
    nt, ns = int(lm.group(1)), int(lm.group(2))

    stations: list[str] = []
    grab = False
    for h in header:
        body = h.lstrip("#").strip()
        if body.startswith("Assimilated stations"):
            grab = True
            continue
        if grab:
            if body.startswith("Requested by") or "END_HEADER" in body or not body:
                break
            stations.extend(body.split())

    try:
        coeffs = np.array([float(x) for x in numbers], dtype=float)
    except ValueError as exc:
        raise FetchError(
            f"IRTAM file has an unreadable coefficient ({exc}); the download is "
            "probably truncated"
        ) from exc
    if coeffs.size != nt * ns:
        raise FetchError(
            f"IRTAM file has {coeffs.size} coefficients but its header declares "
            f"{nt} x {ns} = {nt * ns}; the download is probably truncated"
        )

    return IrtamCoefficients(
        characteristic=char,
        unit=unit,
        valid_time=valid_time,
        engine=field("Assimilative Engine"),
        basis=field("Expansion Basis"),
        n_temporal=nt,
        n_spatial=ns,
        grid=field("Earth Grid"),
        stations=tuple(stations),
        generated_by=field("Generated by"),
        coefficients=coeffs,
        header=tuple(header),
    )


# ---------------------------------------------------------------------------
# Fetching
# ---------------------------------------------------------------------------


def fetch_coefficients(
    characteristic: str,
    t: str | date | datetime,
    *,
    route: str = "auto",
    force: bool = False,
) -> tuple[IrtamCoefficients, Path]:
    """Fetch, cache and parse one coefficient file.

    ``route`` is ``"archive"``, ``"servlet"`` or ``"auto"``.  Auto tries the
    archive and falls back to the servlet only on a 404, which is how the
    archive says it does not hold that epoch yet.  Any other failure is
    raised, because a 5xx from the archive is not evidence the servlet will
    do better.

    Coefficients for a past epoch do not change, so a cached file is never
    stale.  The one exception is a servlet file for an epoch within the last
    day, which may be recomputed as late data arrive; it is re-fetched after
    an hour.
    """
    if route not in ("auto", "archive", "servlet"):
        raise ValueError(f"route must be 'auto', 'archive' or 'servlet', not {route!r}")
    _check_char(characteristic)
    epoch = round_to_cadence(t)
    recent = datetime.now(timezone.utc) - epoch < timedelta(days=1)

    order = ["archive", "servlet"] if route == "auto" else [route]
    if CHARACTERISTICS[characteristic]["servlet"] is None and "servlet" in order:
        order = [r for r in order if r != "servlet"] or order

    def _where(r: str) -> tuple[str, Path, float]:
        url = archive_url(characteristic, epoch) if r == "archive" else servlet_url(characteristic, epoch)
        path = cache.cached_path(SOURCE, url, ".ASC" if r == "archive" else ".txt")
        max_age = (1.0 / 24.0) if (r == "servlet" and recent) else float("inf")
        return url, path, max_age

    # Every route's cache is checked before any route touches the network, so an
    # epoch the archive lacks does not cost an archive 404 on every later call.
    if not force:
        for r in order:
            _url, path, max_age = _where(r)
            if cache.is_fresh(path, max_age) and cache.verify_cached(path):
                return parse_coefficients(path.read_bytes()), path

    last_error: FetchError | None = None
    for r in order:
        url, path, _max_age = _where(r)
        try:
            response = cache.polite_get(url, timeout=120.0, min_interval_s=MIN_INTERVAL_S)
        except FetchError as exc:
            last_error = exc
            if route == "auto" and r == "archive" and " 404 " in f" {exc} ":
                continue
            raise
        body = response.content
        parsed = parse_coefficients(body)  # refuse to cache a malformed file
        path.write_bytes(body)
        headers = getattr(response, "headers", {}) or {}
        cache.record_meta(
            path,
            url=url,
            source=SOURCE,
            etag=headers.get("ETag"),
            last_modified=headers.get("Last-Modified"),
            extra={"route": r, "characteristic": characteristic, "epoch": epoch.isoformat()},
        )
        return parsed, path

    raise last_error or FetchError(f"no IRTAM route could serve {characteristic} at {epoch}")


@register_source(
    name=SOURCE,
    description="IRTAM assimilative IRI coefficients (foF2, hmF2, B0, B1) every 15 minutes",
    url="https://giro.uml.edu/IRTAM/",
    licence=(
        "Distributed by the Lowell GIRO Data Center; derived from GIRO measurements, "
        "which are released under CC-BY-NC-SA 4.0 with the GIRO Rules of the Road"
    ),
    citation=(
        "Galkin, I. A., B. W. Reinisch, X. Huang and D. Bilitza (2012), Assimilation of "
        "GIRO data into a real-time IRI, Radio Sci., 47, RS0L07, doi:10.1029/2011RS004952; "
        "and Galkin, I. A., et al. (2020), Assimilation of sparse continuous near-Earth "
        "weather measurements by NECTAR model morphing, Space Weather, "
        "doi:10.1029/2020SW002463."
    ),
    cadence="every 15 minutes; static archive 2000-2025, servlet for recent epochs",
    work_packages=("WP1",),
    acknowledgement="IRTAM coefficients courtesy of the Lowell GIRO Data Center.",
    notes=(
        "Evaluated into maps by tss.l0_ionosphere.irtam_maps (verified against PyIRTAM; "
        "see docs/irtam-validation.md). IRTAM assimilates GIRO data, so it is not an "
        "independent check on GIRO stations it assimilated."
    ),
    extra={"archive_root": ARCHIVE_ROOT, "servlet": SERVLET, "cadence_min": CADENCE_MIN},
)
def fetch_irtam_sample(
    *,
    force: bool = False,
    limit: int | None = None,
    time: str | None = None,
    characteristics: tuple[str, ...] = ("foF2", "hmF2", "B0", "B1"),
    route: str = "auto",
) -> list[Path]:
    """Fetch one epoch of IRTAM coefficients for each characteristic.

    Without arguments this takes 12:00 UT three days ago.  Not yesterday: on
    2026-09-23 the servlet answered HTTP 500 for every characteristic at every
    epoch tried on 2026-09-22 and for "now", while every epoch through
    2026-09-21 worked, so the most recent day or two are not reliably
    available.  ``limit`` caps the number of characteristics fetched.
    """
    if time is None:
        today = datetime.now(timezone.utc).date()
        time = datetime(today.year, today.month, today.day, 12, tzinfo=timezone.utc) - timedelta(
            days=RECENT_LAG_DAYS
        )
    chars = list(characteristics)[: limit if limit is not None else None]
    return [fetch_coefficients(c, time, route=route, force=force)[1] for c in chars]
