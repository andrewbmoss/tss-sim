"""Layer: data
WP: WP1, WP8
Claims: -
Status: implemented
Requirements: -

GIRO DIDBase, scaled ionogram characteristics from the global Digisonde network.

What it provides.  The Lowell GIRO Data Center holds the ionograms and the
scaled characteristics of roughly 120 Digisondes.  The characteristics that
matter here are foF2, hmF2, B0 and B1, which are the F2 peak and the IRI
bottomside shape parameters, plus fmin.  Each record carries the autoscaling
confidence score (CS) and a two-character qualifying/descriptive code per value.

What it is good for here.  Ground truth for the peak only.  A topside sounder
retrieves the profile from the spacecraft down to hmF2, and the bottom of that
retrieval is exactly the foF2 and hmF2 a ground ionosonde measures from below.
A satellite pass near a Digisonde is therefore a direct, independent check on
the retrieved peak.  It says nothing about the topside shape above hmF2, which
an ionosonde cannot see; that comes from :mod:`tss.data.madrigal`.

The Western Pacific stations (Guam, Kwajalein, Okinawa, Wake, Darwin,
Townsville, Learmonth, Niue, Norfolk) are listed in :data:`WESTERN_PACIFIC`
because the bistatic OTH side study needs them as ionospheric state near the
JORN coverage.

The endpoint.  ``https://lgdc.uml.edu/fastchar/getbest`` is a GET service that
takes the station, a comma-separated list of characteristic names and a UTC
window, and returns the plain-text table.  It is what the DIDBase web form at
``giro.uml.edu/didbase/scaled.php`` redirects to: the form answers a POST with a
302 to this URL, translating its positional option indices into names on the
way.  Calling the service directly saves the redirect and removes the index
table the form route needed, which broke silently if the form's option list
was ever reordered.  Checked on 2026-09-23: for the same query the two routes
return identical bodies apart from the generation timestamp.

The older GET service ``lgdc.uml.edu/common/DIDBGetValues``, which most
published scripts use, answers 404 "resource not available"; it should not be
reintroduced.

Licence.  The response header states that all GIRO measurements are released
under CC-BY-NC-SA 4.0, that the Lowell GIRO Data Center Rules of the Road
apply, and that the data provider of each station must be acknowledged.  The
non-commercial clause is fine for the dissertation.  Whether it covers
Navy-sponsored use in the bistatic side study is a question for the Lowell GIRO
Data Center, not something to assume.

Privacy.  The server writes the requester's public IP address into the
response header.  The cached file keeps the bytes as served, because raw bytes
are never edited, and ``data/external/`` is not committed.  Checked-in fixtures
have that line redacted.
"""

from __future__ import annotations

import re
from collections.abc import Sequence
from dataclasses import dataclass, field
from datetime import date, datetime, timedelta, timezone
from pathlib import Path
from typing import Any

from tss.data import cache
from tss.data.base import FetchError, register_source

__all__ = [
    "ENDPOINT",
    "SOURCE",
    "CHARACTERISTICS",
    "SERVER_NAMES",
    "DEFAULT_CHARACTERISTICS",
    "WESTERN_PACIFIC",
    "MAX_WINDOW_DAYS",
    "MIN_INTERVAL_S",
    "GiroQueryError",
    "ScaledCharacteristics",
    "build_query",
    "request_url",
    "parse_scaled",
    "query_scaled",
    "fetch_giro_sample",
]

ENDPOINT = "https://lgdc.uml.edu/fastchar/getbest"
SOURCE = "giro_didbase"

MIN_INTERVAL_S = 3.0
"""Minimum gap between two queries, seconds.  Three times the package default."""

MAX_WINDOW_DAYS = 31
"""Longest window sent in one request.  Longer windows are split into chunks
and each chunk is cached separately, so an interrupted collection resumes."""

SERVER_NAMES: tuple[str, ...] = (
    "foF2", "foF1", "foE", "foEs", "fbEs", "foEa", "foP", "fxI", "MUF(D)", "M(D)",
    "h`F2", "h`F", "h`E", "h`Es", "h`Ea", "h`P", "hmF2", "hmF1", "hmE", "zhalfNm",
    "yF2", "yF1", "yE", "scaleF2", "B0", "B1", "D1", "TEC", "FF", "FE", "QF", "QE",
    "fmin", "fminF", "fminE", "fminEs", "foF2p",
)
"""Characteristic names as the service spells them, read from the DIDBase form's
redirect on 2026-09-23.  The virtual heights use a backtick, not an apostrophe.

``Type Es`` is offered by the form but left out here: it is a classification,
not a number, the service drops it from the table when asked, and its space
would break the whitespace-separated column header anyway.
"""

CHARACTERISTICS: dict[str, str] = {n: n for n in SERVER_NAMES}
CHARACTERISTICS.update({
    "MUFD": "MUF(D)", "MD": "M(D)",
    "h'F2": "h`F2", "h'F": "h`F", "h'E": "h`E", "h'Es": "h`Es", "h'Ea": "h`Ea", "h'P": "h`P",
})
"""Accepted characteristic name to the name the service expects.

Every server spelling is accepted as itself, and a few spellings that are easier
to type map onto it: ``MUFD`` and ``MD`` for ``MUF(D)`` and ``M(D)``, and an
apostrophe for the backtick in the virtual heights.  Column names in a result
are always the server's spelling, because they are read from the response.
"""

DEFAULT_CHARACTERISTICS: tuple[str, ...] = ("foF2", "hmF2", "B0", "B1", "fmin")

WESTERN_PACIFIC: dict[str, str] = {
    "GU513": "Guam",
    "KJ609": "Kwajalein",
    "OK426": "Okinawa",
    "WA619": "Wake Island",
    "DW41K": "Darwin",
    "TV51R": "Townsville",
    "LM42B": "Learmonth",
    "ND61R": "Niue",
    "NI63_": "Norfolk Island",
    "CS31K": "Cocos Island",
}
"""Digisondes relevant to the Western Pacific bistatic side study, by URSI code."""

_URSI = re.compile(r"^[A-Z0-9_]{5}$")
_MISSING = {"---", "--", "-", ""}


class GiroQueryError(FetchError):
    """DIDBase answered, but with an error status other than "no data"."""


@dataclass(frozen=True)
class ScaledCharacteristics:
    """One DIDBase response, parsed.

    ``values`` maps a characteristic name, as the server spells it in the
    column header (``foF2``, ``hmF2``, ``B0``...), to a float array with NaN
    where the server wrote ``---``.  ``qd`` holds the matching two-character
    qualifying/descriptive codes.  ``cs`` is the autoscaling confidence score:
    0 to 100, 999 for a manually scaled record, -1 if unknown.

    ``units`` maps each characteristic to the unit the header states, or an
    empty string for a dimensionless one such as B1.  ``status`` is ``"ok"``
    or the server's no-data message; an empty result is an answer, not an
    error.
    """

    ursi_code: str
    station: str
    lat_deg: float
    lon_deg: float
    instrument: str
    times: Any
    cs: Any
    values: dict[str, Any]
    qd: dict[str, Any]
    units: dict[str, str] = field(default_factory=dict)
    status: str = "ok"
    header: tuple[str, ...] = ()

    def __len__(self) -> int:
        return len(self.times)

    def good(self, min_cs: int = 70) -> Any:
        """Boolean mask of records whose confidence score is at least ``min_cs``.

        A manually scaled record (CS 999) always passes.  An unknown score (-1)
        never does.  70 is a common autoscaling threshold in the Digisonde
        literature, not a GIRO recommendation; pick it deliberately for a study.
        """
        import numpy as np

        cs = np.asarray(self.cs)
        return (cs >= min_cs) | (cs == 999)


# ---------------------------------------------------------------------------
# Request building
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


def _check_ursi(code: str) -> str:
    code = str(code).strip().upper()
    if not _URSI.match(code):
        raise ValueError(
            f"{code!r} is not an URSI station code; these are five characters, "
            "e.g. 'MHJ45' for Millstone Hill or 'GU513' for Guam"
        )
    return code


def build_query(
    ursi_code: str,
    start: str | date | datetime,
    end: str | date | datetime,
    characteristics: Sequence[str] = DEFAULT_CHARACTERISTICS,
    *,
    dmuf_km: int = 3000,
) -> list[tuple[str, str]]:
    """The query parameters for one request, as ``(name, value)`` pairs.

    Split out so it can be tested without a network.  Times are UTC, sent as
    ``YYYY/MM/DD hh:mm:ss``.  Unknown characteristic names raise here rather
    than being sent: the service would drop them with a warning, and if none
    were known it would fail the whole query.  Duplicates are sent once.
    """
    code = _check_ursi(ursi_code)
    t0, t1 = _as_utc(start), _as_utc(end)
    if t1 <= t0:
        raise ValueError(f"end {t1.isoformat()} is not after start {t0.isoformat()}")
    if not characteristics:
        raise ValueError("at least one characteristic is needed")
    unknown = [c for c in characteristics if c not in CHARACTERISTICS]
    if unknown:
        raise ValueError(
            f"unknown characteristic(s) {', '.join(unknown)}; known: "
            f"{', '.join(CHARACTERISTICS)}"
        )

    names = list(dict.fromkeys(CHARACTERISTICS[c] for c in characteristics))
    return [
        ("ursiCode", code),
        ("charName", ",".join(names)),
        ("DMUF", str(int(dmuf_km))),
        ("fromDate", t0.strftime("%Y/%m/%d %H:%M:%S")),
        ("toDate", t1.strftime("%Y/%m/%d %H:%M:%S")),
    ]


def _windows(t0: datetime, t1: datetime, days: int) -> list[tuple[datetime, datetime]]:
    out = []
    step = timedelta(days=days)
    cur = t0
    while cur < t1:
        nxt = min(cur + step, t1)
        out.append((cur, nxt))
        cur = nxt
    return out


# ---------------------------------------------------------------------------
# Parsing
# ---------------------------------------------------------------------------

_LOCATION = re.compile(
    r"GEO\s*\(\s*([-\d.]+)\s*([NS])\s+([-\d.]+)\s*([EW])\s*\)\s*,\s*URSI-Code\s+(\S+)\s*,\s*(.+)$"
)
_UNIT = re.compile(r"^#\s*([^\s\[=]+)\s*(?:\[([^\]]*)\])?\s*=")


def parse_scaled(text: str | bytes) -> ScaledCharacteristics:
    """Parse a DIDBase response body.

    Raises :class:`GiroQueryError` when the server reports an error other than
    "no measurement data", for example an unknown characteristic.  Returns an
    empty result, with ``status`` set to the server's message, when the station
    simply has no data in the window.
    """
    import numpy as np

    if isinstance(text, bytes):
        text = text.decode("utf-8", "replace")
    lines = text.splitlines()
    header = [ln for ln in lines if ln.startswith("#")]
    body = [ln for ln in lines if ln.strip() and not ln.startswith("#")]

    if not header or "DIDBase" not in text[:400]:
        snippet = " ".join(text.split())[:160]
        raise GiroQueryError(f"response does not look like a DIDBase table: {snippet!r}")

    ursi, station, lat, lon, instrument = "", "", float("nan"), float("nan"), ""
    units: dict[str, str] = {}
    status = "ok"
    col_line = ""
    for ln in header:
        content = ln[1:].strip()
        if content.startswith("Location:"):
            m = _LOCATION.search(content)
            if m:
                lat = float(m.group(1)) * (1 if m.group(2) == "N" else -1)
                lon = float(m.group(3)) * (1 if m.group(4) == "E" else -1)
                ursi, station = m.group(5), m.group(6).strip()
        elif content.startswith("Instrument:"):
            instrument = content.split(":", 1)[1].strip()
        elif content.startswith("STATUS:"):
            status = content.split(":", 1)[1].strip()
        elif content.startswith("Time"):
            col_line = content
        else:
            m = _UNIT.match(ln)
            if m and m.group(1) not in ("CS",):
                units[m.group(1)] = (m.group(2) or "").strip()

    if status != "ok" and status.upper().startswith("ERROR"):
        if "no measurement data" in status.lower():
            return ScaledCharacteristics(
                ursi_code=ursi, station=station, lat_deg=lat, lon_deg=lon,
                instrument=instrument, times=np.array([], dtype="datetime64[ms]"),
                cs=np.array([], dtype=int), values={}, qd={}, units=units,
                status=status, header=tuple(header),
            )
        raise GiroQueryError(f"DIDBase refused the query for {ursi or 'station'}: {status}")

    if not col_line:
        raise GiroQueryError("DIDBase response has no column header line")

    cols = col_line.split()
    # Layout: Time CS name QD name QD ...
    if len(cols) < 2 or cols[0] != "Time" or cols[1] != "CS":
        raise GiroQueryError(f"unexpected column header {col_line!r}")
    names = cols[2::2]
    if any(q != "QD" for q in cols[3::2]):
        raise GiroQueryError(f"unexpected column header {col_line!r}")

    times: list[Any] = []
    cs: list[int] = []
    vals: dict[str, list[float]] = {n: [] for n in names}
    qds: dict[str, list[str]] = {n: [] for n in names}
    for ln in body:
        parts = ln.split()
        if len(parts) != 2 + 2 * len(names):
            raise GiroQueryError(f"row has {len(parts)} fields, expected {2 + 2 * len(names)}: {ln!r}")
        times.append(np.datetime64(parts[0].rstrip("Z"), "ms"))
        cs.append(int(parts[1]))
        for i, n in enumerate(names):
            raw = parts[2 + 2 * i]
            vals[n].append(float("nan") if raw in _MISSING else float(raw))
            qds[n].append(parts[3 + 2 * i])

    return ScaledCharacteristics(
        ursi_code=ursi,
        station=station,
        lat_deg=lat,
        lon_deg=lon,
        instrument=instrument,
        times=np.array(times, dtype="datetime64[ms]"),
        cs=np.array(cs, dtype=int),
        values={n: np.array(v, dtype=float) for n, v in vals.items()},
        qd={n: np.array(v, dtype=object) for n, v in qds.items()},
        units=units,
        status=status,
        header=tuple(header),
    )


# ---------------------------------------------------------------------------
# Fetching
# ---------------------------------------------------------------------------


def request_url(params: Sequence[tuple[str, str]]) -> str:
    """The full request URL, encoded as the service's own redirect encodes it."""
    from urllib.parse import urlencode

    return ENDPOINT + "?" + urlencode(list(params), safe=",")


def query_scaled(
    ursi_code: str,
    start: str | date | datetime,
    end: str | date | datetime,
    characteristics: Sequence[str] = DEFAULT_CHARACTERISTICS,
    *,
    force: bool = False,
    max_age_days: float = 30.0,
    max_chunks: int | None = None,
) -> tuple[ScaledCharacteristics, list[Path]]:
    """Query one station over a window and return the merged result and files.

    The window is split into chunks of at most :data:`MAX_WINDOW_DAYS`.  Each
    chunk is one GET, cached under ``data/external/giro_didbase/`` with the
    usual sidecar, so a long collection that is interrupted picks up where it
    stopped.  ``max_chunks`` caps how many chunks are fetched; the result then
    covers only the start of the window and says so in ``status``.

    Chunk boundaries are inclusive at both ends on the server, so a record that
    falls exactly on a boundary would appear twice.  Duplicate timestamps are
    dropped when the chunks are merged.

    This function makes network requests unless every chunk is cached.
    """
    import numpy as np

    t0, t1 = _as_utc(start), _as_utc(end)
    windows = _windows(t0, t1, MAX_WINDOW_DAYS)
    truncated = max_chunks is not None and len(windows) > max_chunks
    if max_chunks is not None:
        windows = windows[: max(0, int(max_chunks))]

    parts: list[ScaledCharacteristics] = []
    paths: list[Path] = []
    for w0, w1 in windows:
        params = build_query(ursi_code, w0, w1, characteristics)
        url = request_url(params)
        path = cache.cached_path(SOURCE, url, ".txt")
        if force or not (cache.is_fresh(path, max_age_days) and cache.verify_cached(path)):
            # The URL is sent pre-encoded rather than through requests' params,
            # so the bytes on the wire are exactly the URL in the sidecar.
            response = cache.polite_get(url, timeout=180.0, min_interval_s=MIN_INTERVAL_S)
            body = response.content
            if not body.strip():
                raise FetchError(f"DIDBase returned an empty body for {ursi_code} {w0} to {w1}")
            path.write_bytes(body)
            cache.record_meta(path, url=url, source=SOURCE, extra={"params": dict(params)})
        parts.append(parse_scaled(path.read_bytes()))
        paths.append(path)

    merged = _merge(parts, np)
    if len(merged) == 0:
        # An empty answer still carries the columns that were asked for, so a
        # caller can index r.values["foF2"] and get an empty array rather than
        # a KeyError on a station that simply had no data that day.
        from dataclasses import replace

        merged = replace(
            merged,
            values={CHARACTERISTICS[c]: np.array([], dtype=float) for c in characteristics},
            qd={CHARACTERISTICS[c]: np.array([], dtype=object) for c in characteristics},
        )
    if truncated:
        merged = _replace_status(merged, f"truncated to the first {max_chunks} chunk(s)")
    return merged, paths


def _replace_status(r: ScaledCharacteristics, status: str) -> ScaledCharacteristics:
    from dataclasses import replace

    return replace(r, status=status)


def _merge(parts: list[ScaledCharacteristics], np: Any) -> ScaledCharacteristics:
    data = [p for p in parts if len(p)]
    if not data:
        if parts:
            return parts[0]
        raise FetchError("no chunks were queried")
    first = data[0]
    names = list(first.values)
    times = np.concatenate([p.times for p in data])
    order = np.argsort(times, kind="stable")
    times = times[order]
    keep = np.ones(len(times), dtype=bool)
    keep[1:] = times[1:] != times[:-1]
    idx = order[keep]

    def cat(get: Any) -> Any:
        return np.concatenate([get(p) for p in data])[idx]

    return ScaledCharacteristics(
        ursi_code=first.ursi_code,
        station=first.station,
        lat_deg=first.lat_deg,
        lon_deg=first.lon_deg,
        instrument=first.instrument,
        times=times[keep],
        cs=cat(lambda p: p.cs),
        values={n: cat(lambda p, n=n: p.values[n]) for n in names},
        qd={n: cat(lambda p, n=n: p.qd[n]) for n in names},
        units=first.units,
        status="ok",
        header=first.header,
    )


@register_source(
    name=SOURCE,
    description="GIRO DIDBase scaled Digisonde characteristics (foF2, hmF2, B0, B1)",
    url="https://giro.uml.edu/didbase/scaled.php",  # the human-facing page
    licence=(
        "CC-BY-NC-SA 4.0, as stated in every DIDBase response; Lowell GIRO Data Center "
        "Rules of the Road apply (https://ulcar.uml.edu/DIDBase/RulesOfTheRoadForDIDBase.htm)"
    ),
    citation=(
        "Reinisch, B. W. and I. A. Galkin (2011), Global Ionospheric Radio Observatory "
        "(GIRO), Earth Planets Space, 63, 377-381, doi:10.5047/eps.2011.03.001. "
        "Data from the Lowell GIRO Data Center DIDBase."
    ),
    cadence="live; stations report every 5 to 15 minutes, archive back to the 1990s",
    work_packages=("WP1", "WP8"),
    acknowledgement=(
        "Ionosonde data from the Lowell GIRO Data Center. The data provider of each "
        "station used must be acknowledged by name, as the Rules of the Road require."
    ),
    notes=(
        "Truth for the F2 peak only; an ionosonde cannot see the topside. "
        "Non-commercial licence: confirm with the Lowell GIRO Data Center before using "
        "these data in Navy-sponsored work. The response header carries the requester's "
        "IP address; cached files keep it, fixtures redact it."
    ),
    extra={"endpoint": ENDPOINT, "max_window_days": MAX_WINDOW_DAYS},
)
def fetch_giro_sample(
    *,
    force: bool = False,
    limit: int | None = None,
    station: str = "MHJ45",
    start: str | None = None,
    end: str | None = None,
    characteristics: Sequence[str] = DEFAULT_CHARACTERISTICS,
) -> list[Path]:
    """Fetch scaled characteristics for one station and window, and cache them.

    Without arguments this fetches one day at Millstone Hill (MHJ45), the most
    recent complete UTC day, so ``tss fetch giro_didbase`` does something
    sensible.  ``limit`` caps the number of 31-day chunks requested.
    """
    if start is None or end is None:
        today = datetime.now(timezone.utc).date()
        t1 = datetime(today.year, today.month, today.day, tzinfo=timezone.utc)
        start = start or (t1 - timedelta(days=1)).isoformat()
        end = end or t1.isoformat()
    _result, paths = query_scaled(
        station, start, end, characteristics, force=force, max_chunks=limit
    )
    return paths
