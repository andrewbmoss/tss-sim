"""Layer: data
WP: WP2, WP12
Claims: -
Status: implemented
Requirements: -

wspr.live, the queryable archive of WSPR reception reports.

What WSPR spots are.  WSPR is a weak-signal beacon mode.  A transmitter sends a
two-minute frame carrying its callsign, grid square and power at a known slot
time; every receiver that decodes it reports the decode with a signal-to-noise
ratio referred to a 2.5 kHz bandwidth.  Several thousand stations do this
continuously across the HF bands, so the aggregate is a dense, free,
continuously updated record of which HF paths were open, at what frequency, at
what hour, at what received SNR, for a transmit power that is usually stated to
within a decibel.

What it is good for here.  WP2 models integrated D-region absorption.  An
absorption model that is wrong by several decibels in the diurnal sense will
disagree with the WSPR record in a way that is visible without any calibration
at all: the ratio of night to day spot counts on a given band, the lowest band
that carries spots at local noon, and the drop in median SNR through sunrise are
all consequences of absorption and all measurable from these data.  Used that
way, WSPR is an independent check on the diurnal and frequency dependence of the
absorption model, and on band behaviour generally, without needing an
instrument.

What it is not.  WSPR is ground-to-ground oblique propagation.  The rays leave
a ground antenna at a low elevation angle, pass twice through the D region on
each hop, and reflect from the bottomside.  A topside sounder looks downward
from above the F2 peak through a completely different geometry and a single
one-way D-region traverse only in the limit where the wave penetrates, which it
mostly does not.  Nothing measured here validates a topside result.  What it
provides is a bound and a sanity check: an absorption model that cannot
reproduce the ground-to-ground record at the same frequencies and the same
solar zenith angles is wrong in a way that will also be wrong topside, while
one that does reproduce it has passed a necessary and not a sufficient test.
That asymmetry is the whole of the claim being made from these data and the
dissertation should say so plainly wherever it uses them.

The endpoint.  wspr.live exposes a ClickHouse HTTP interface at
``https://db1.wspr.live/`` which takes the SQL in a ``query`` parameter.  The
convention is to append ``FORMAT CSVWithNames`` or ``FORMAT JSONEachRow``.  The
table of reception reports is ``wspr.rx``.

Respecting the site.  This is one person's server carrying a large database as a
public service.  Every query here carries an enforced ``LIMIT``, the default
inter-request interval is well above the package default, and no function in
this module loops over a date range issuing a query per day.  Anything that
wants more than a few hundred thousand rows should take the published archive
files that wspr.live provides for bulk download instead of the live endpoint;
see :data:`BULK_NOTE`.  A query that is too expensive is more likely to be
killed than answered, so a narrow query is also the fast one.
"""

from __future__ import annotations

from collections.abc import Iterable, Sequence
from datetime import date, datetime, timedelta, timezone
from pathlib import Path
from typing import Any

from tss.data import cache
from tss.data.base import FetchError, register_source

__all__ = [
    "ENDPOINT",
    "TABLE",
    "DEFAULT_LIMIT",
    "MAX_LIMIT",
    "MIN_INTERVAL_S",
    "SPOT_COLUMNS",
    "BANDS_MHZ",
    "BULK_NOTE",
    "query",
    "spots",
    "build_spots_sql",
    "quote",
    "fetch_wspr_sample",
]


ENDPOINT = "https://db1.wspr.live/"
TABLE = "wspr.rx"

DEFAULT_LIMIT = 100_000
"""Rows returned when the caller does not say.  Large enough for a day of one
band, small enough that a mistake does not cost the server much."""

MAX_LIMIT = 1_000_000
"""Hard ceiling.  A request above this is clamped, not refused, and the clamp is
reported in the frame's ``attrs``.  Anything needing more belongs in the bulk
archive."""

MIN_INTERVAL_S = 5.0
"""Minimum gap between live queries, seconds.  Five times the package default."""

BULK_NOTE = (
    "wspr.live publishes monthly archive files for bulk use.  Take those rather "
    "than paging the live endpoint.  The live endpoint is for narrow, targeted "
    "queries, and a job that would issue more than a handful of them should be "
    "rewritten against the archive files."
)

SPOT_COLUMNS: tuple[str, ...] = (
    "time",
    "band",
    "rx_sign",
    "rx_lat",
    "rx_lon",
    "tx_sign",
    "tx_lat",
    "tx_lon",
    "distance",
    "azimuth",
    "frequency",
    "power",
    "snr",
    "drift",
)
"""The columns of ``wspr.rx`` this module selects.

``frequency`` is in Hz as reported, ``power`` is the transmitter's stated power
in dBm, ``snr`` is in dB referred to 2.5 kHz, ``distance`` is in km and
``band`` is the integer band code described at :data:`BANDS_MHZ`.  The table has
more columns than these; the selection is deliberate because ``SELECT *`` on a
wide columnar table costs the server far more than it costs us.
"""

BANDS_MHZ: dict[int, float] = {
    -1: 0.1375,
    0: 0.4746,
    1: 1.8366,
    3: 3.5686,
    5: 5.2887,
    7: 7.0386,
    10: 10.1387,
    14: 14.0956,
    18: 18.1046,
    21: 21.0946,
    24: 24.9246,
    28: 28.1246,
}
"""WSPR band code to dial frequency in MHz.

The code is the band in megahertz rounded down, with 0 for 630 m and -1 for
2200 m.  Only 1.8 MHz through 28 MHz overlap the 1 to 20 MHz sounding band that
INS-01 covers, and only the lower half of that is useful as an absorption check,
because above about 15 MHz daytime absorption is small and the spot record is
dominated by whether the path supported the frequency at all.
"""


def quote(value: Any) -> str:
    """Render a Python value as a ClickHouse SQL literal.

    Strings are escaped for backslash and single quote and wrapped in single
    quotes.  Datetimes are rendered as ``'YYYY-MM-DD hh:mm:ss'`` in UTC.
    Numbers and booleans pass through as literals.  Nothing else is accepted,
    because a value whose type this function does not recognise is exactly the
    case where string formatting turns into an injection.

    This is the only place in the module that builds a literal.  There is no
    parameter binding over the ClickHouse HTTP interface as it is used here, so
    quoting correctly is the whole of the defence, and every caller goes
    through this function.
    """
    if isinstance(value, bool):
        return "1" if value else "0"
    if isinstance(value, (int, float)):
        if value != value or value in (float("inf"), float("-inf")):
            raise ValueError(f"cannot quote non-finite number {value!r}")
        return repr(value)
    if isinstance(value, datetime):
        dt = value.astimezone(timezone.utc) if value.tzinfo else value.replace(tzinfo=timezone.utc)
        return "'" + dt.strftime("%Y-%m-%d %H:%M:%S") + "'"
    if isinstance(value, date):
        return "'" + value.strftime("%Y-%m-%d 00:00:00") + "'"
    if isinstance(value, str):
        escaped = value.replace("\\", "\\\\").replace("'", "\\'")
        return "'" + escaped + "'"
    raise TypeError(
        f"refusing to quote a {type(value).__name__} into SQL; pass a str, a number, "
        "a date or a datetime"
    )


def _as_datetime(value: str | date | datetime) -> datetime:
    """Accept an ISO string, a date or a datetime and return an aware datetime."""
    if isinstance(value, datetime):
        return value if value.tzinfo else value.replace(tzinfo=timezone.utc)
    if isinstance(value, date):
        return datetime(value.year, value.month, value.day, tzinfo=timezone.utc)
    if isinstance(value, str):
        text = value.strip().replace("Z", "+00:00")
        try:
            dt = datetime.fromisoformat(text)
        except ValueError as exc:
            raise ValueError(
                f"cannot read {value!r} as a time; use ISO 8601, e.g. "
                "'2024-03-01' or '2024-03-01T12:00:00Z'"
            ) from exc
        return dt if dt.tzinfo else dt.replace(tzinfo=timezone.utc)
    raise TypeError(f"cannot read a {type(value).__name__} as a time")


def _clamp_limit(limit: int | None) -> tuple[int, bool]:
    """Return the limit actually used and whether it was clamped."""
    if limit is None:
        return DEFAULT_LIMIT, False
    n = int(limit)
    if n <= 0:
        raise ValueError(f"limit must be positive, got {limit!r}")
    if n > MAX_LIMIT:
        return MAX_LIMIT, True
    return n, False


def build_spots_sql(
    start: str | date | datetime,
    end: str | date | datetime,
    *,
    band: int | Sequence[int] | None = None,
    tx_sign: str | None = None,
    rx_sign: str | None = None,
    limit: int | None = None,
    columns: Sequence[str] = SPOT_COLUMNS,
    fmt: str | None = None,
) -> str:
    """Build the SQL for a spot query.

    Split out from :func:`spots` so that it can be tested without any network.
    Every literal goes through :func:`quote`, the column list is checked against
    :data:`SPOT_COLUMNS` rather than interpolated freely, and a ``LIMIT`` is
    always emitted.  Callsigns are upper-cased because the table stores them
    that way, and a lower-case argument silently matching nothing is a worse
    failure than a surprising case change.

    The time window is half-open, ``start <= time < end``, so that consecutive
    windows tile without double-counting a spot on the boundary.
    """
    t0 = _as_datetime(start)
    t1 = _as_datetime(end)
    if t1 <= t0:
        raise ValueError(f"end {t1.isoformat()} is not after start {t0.isoformat()}")

    unknown = [c for c in columns if c not in SPOT_COLUMNS]
    if unknown:
        raise ValueError(
            f"unknown column(s) {', '.join(unknown)}; this module selects from "
            f"{', '.join(SPOT_COLUMNS)}"
        )
    if not columns:
        raise ValueError("at least one column must be selected")

    where: list[str] = [
        f"time >= {quote(t0)}",
        f"time < {quote(t1)}",
    ]

    if band is not None:
        if isinstance(band, (int,)) and not isinstance(band, bool):
            where.append(f"band = {quote(int(band))}")
        elif isinstance(band, Iterable) and not isinstance(band, (str, bytes)):
            codes = [int(b) for b in band]
            if not codes:
                raise ValueError("band was an empty sequence")
            where.append("band IN (" + ", ".join(quote(b) for b in codes) + ")")
        else:
            raise TypeError(f"band must be an int or a sequence of ints, got {band!r}")

    if tx_sign is not None:
        where.append(f"tx_sign = {quote(str(tx_sign).strip().upper())}")
    if rx_sign is not None:
        where.append(f"rx_sign = {quote(str(rx_sign).strip().upper())}")

    n, _clamped = _clamp_limit(limit)

    sql = (
        "SELECT " + ", ".join(columns)
        + f" FROM {TABLE}"
        + " WHERE " + " AND ".join(where)
        + " ORDER BY time"
        + f" LIMIT {n}"
    )
    if fmt:
        sql += f" FORMAT {fmt}"
    return sql


def query(
    sql: str,
    *,
    fmt: str = "CSVWithNames",
    timeout: float = 120.0,
    min_interval_s: float = MIN_INTERVAL_S,
    cache_result: bool = True,
    force: bool = False,
    max_age_days: float = 1.0,
) -> Any:
    """Run one SQL query against wspr.live and return a pandas DataFrame.

    ``fmt`` is appended as ``FORMAT <fmt>`` unless the SQL already ends in a
    FORMAT clause.  ``CSVWithNames`` and ``JSONEachRow`` are both handled;
    anything else is sent as given and returned as a single-column frame of raw
    text, because guessing at a parser for an unknown format would be worse than
    handing the caller the bytes.

    Results are cached by the hash of the query text under
    ``data/external/wspr/``, with the usual sidecar.  A repeated query inside
    ``max_age_days`` is answered from disk and never touches the network, which
    is both polite and what makes a notebook re-runnable on a train.

    This function makes a network request.  It is never called by the test
    suite and never called at import.
    """
    import pandas as pd

    text = sql.strip().rstrip(";")
    lowered = text.lower()
    if " format " not in lowered[-40:]:
        text = f"{text} FORMAT {fmt}"
        used_fmt = fmt
    else:
        used_fmt = text.rsplit(None, 1)[-1]

    suffix = ".csv" if used_fmt.lower().startswith("csv") else ".txt"
    if used_fmt.lower().startswith("json"):
        suffix = ".jsonl"
    path = cache.cached_path("wspr", text, suffix)

    if cache_result and not force and cache.is_fresh(path, max_age_days):
        raw = path.read_bytes()
    else:
        response = cache.polite_get(
            ENDPOINT,
            params={"query": text},
            timeout=timeout,
            min_interval_s=min_interval_s,
        )
        raw = response.content
        if cache_result:
            path.write_bytes(raw)
            cache.record_meta(
                path,
                url=ENDPOINT,
                source="wspr",
                extra={"query": text, "format": used_fmt},
            )

    return _parse(raw, used_fmt, pd)


def _parse(raw: bytes, fmt: str, pd: Any) -> Any:
    """Turn a ClickHouse response body into a DataFrame."""
    import io

    low = fmt.lower()
    if low.startswith("csv"):
        frame = pd.read_csv(io.BytesIO(raw))
    elif low.startswith("json"):
        frame = pd.read_json(io.BytesIO(raw), lines=True)
    else:
        return pd.DataFrame({"raw": [raw.decode("utf-8", "replace")]})

    if "time" in frame.columns:
        frame["time"] = pd.to_datetime(frame["time"], utc=True, errors="coerce")
    if "band" in frame.columns:
        frame["freq_mhz_nominal"] = frame["band"].map(BANDS_MHZ)
    frame.attrs["endpoint"] = ENDPOINT
    frame.attrs["bulk_note"] = BULK_NOTE
    return frame


def spots(
    start: str | date | datetime,
    end: str | date | datetime,
    *,
    band: int | Sequence[int] | None = None,
    tx_sign: str | None = None,
    rx_sign: str | None = None,
    limit: int | None = DEFAULT_LIMIT,
    columns: Sequence[str] = SPOT_COLUMNS,
    fmt: str = "CSVWithNames",
    **query_kwargs: Any,
) -> Any:
    """Fetch WSPR reception reports over a time window.

    The window is half-open in UTC.  ``band`` takes a WSPR band code or a
    sequence of them, see :data:`BANDS_MHZ`.  ``tx_sign`` and ``rx_sign`` are
    callsigns, quoted and upper-cased.  ``limit`` is always applied and is
    clamped to :data:`MAX_LIMIT`.

    The returned frame carries the query in ``attrs['sql']``, so a figure made
    from it can record exactly what was asked.  A result whose row count equals
    the limit is a truncated result and ``attrs['truncated']`` says so; treat
    that as a reason to narrow the window rather than to raise the limit.
    """
    sql = build_spots_sql(
        start,
        end,
        band=band,
        tx_sign=tx_sign,
        rx_sign=rx_sign,
        limit=limit,
        columns=columns,
    )
    n, clamped = _clamp_limit(limit)
    frame = query(sql, fmt=fmt, **query_kwargs)
    frame.attrs["sql"] = sql
    frame.attrs["limit"] = n
    frame.attrs["limit_clamped"] = clamped
    frame.attrs["truncated"] = len(frame) >= n
    return frame


@register_source(
    name="wspr_live",
    description="WSPR reception reports from wspr.live, for HF path availability and SNR",
    url=ENDPOINT,
    licence=(
        "Reception reports are contributed by amateur operators and redistributed by "
        "wspr.live for public use; attribution to wspr.live and WSPRnet expected"
    ),
    citation=(
        "WSPR spot database, wspr.live (ClickHouse mirror of WSPRnet reception "
        "reports), https://wspr.live/"
    ),
    cadence="live, updated continuously",
    work_packages=("WP2", "WP12"),
    acknowledgement="WSPR spots courtesy of wspr.live and the WSPRnet contributors.",
    notes=BULK_NOTE + (
        "  WSPR is ground-to-ground oblique propagation and is a bound and a sanity "
        "check on the WP2 absorption model, not a topside validation."
    ),
    extra={"table": TABLE, "default_limit": DEFAULT_LIMIT, "max_limit": MAX_LIMIT},
)
def fetch_wspr_sample(
    *,
    force: bool = False,
    limit: int | None = None,
    start: str | None = None,
    end: str | None = None,
    band: int | Sequence[int] | None = (3, 7, 10, 14),
) -> list[Path]:
    """Fetch one sample day of spots on the bands WP2 uses, and cache it.

    This exists so that ``tss fetch wspr_live`` does something sensible without
    arguments.  It is a sample, not a data set: the real WP2 cross-check picks
    its own windows.  ``start`` and ``end`` default to the most recent complete
    UTC day.

    Returns the cached CSV path.  Nothing here downloads more than one day.
    """
    if start is None or end is None:
        today = datetime.now(timezone.utc).date()
        t1 = datetime(today.year, today.month, today.day, tzinfo=timezone.utc)
        t0 = t1 - timedelta(days=1)
        start = start or t0.isoformat()
        end = end or t1.isoformat()

    sql = build_spots_sql(start, end, band=band, limit=limit)
    text = f"{sql} FORMAT CSVWithNames"
    path = cache.cached_path("wspr", text, ".csv")
    if not force and cache.is_fresh(path, 7.0):
        return [path]

    response = cache.polite_get(
        ENDPOINT, params={"query": text}, timeout=180.0, min_interval_s=MIN_INTERVAL_S
    )
    body = response.content
    if not body.strip():
        raise FetchError(
            f"wspr.live returned an empty body for the window {start} to {end}; "
            "check the dates and the band codes"
        )
    path.write_bytes(body)
    cache.record_meta(
        path, url=ENDPOINT, source="wspr_live", extra={"query": text}
    )
    return [path]
