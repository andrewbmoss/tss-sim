"""Layer: data
WP: WP1, WP8
Claims: -
Status: implemented
Requirements: INS-01

Madrigal (CEDAR), incoherent scatter radar profiles and DMSP in-situ density.

What it provides.  The CEDAR Madrigal database at MIT Haystack holds, among
much else, the incoherent scatter radar record (Millstone Hill, Jicamarca,
Arecibo, Poker Flat, Resolute Bay, EISCAT) and the DMSP SSIES in-situ plasma
record.  Two products matter here.

* ISR electron density profiles.  A zenith single-pulse Millstone Hill profile
  runs from about 150 km to about 800 km, through the F2 peak and well into the
  topside, with a per-point error estimate.  That is topside profile shape,
  measured, which no ionosonde and no GNSS product can give.  It is the truth a
  retrieved topside profile is scored against, at a handful of sites and on
  the days the radars ran.
* DMSP in-situ electron density at about 840 to 870 km, at one-second
  cadence.  That is the local plasma frequency at a spacecraft, the quantity
  that sets the INS-01 lower band edge, measured in the current solar cycle
  rather than in the ISIS era.

The endpoint.  Madrigal 3 exposes plain-text web services under the site root:
``getInstrumentsService.py``, ``getExperimentsService.py``,
``getExperimentFilesService.py``, ``isprintService.py`` and ``getMadfile.cgi``.
They are what the ``madrigalWeb`` client calls; this module calls them directly
through :func:`tss.data.cache.polite_get` so that caching, sidecars, throttling
and the User-Agent are the same as for every other source, and so that
``madrigalWeb`` is not a dependency.

Two behaviours of the services that the parsers here exist to handle:

* ``isprintService.py`` takes one ``parms=`` per parameter.  A single
  comma-joined ``parms`` value is rejected.
* It answers HTTP 200 even when it fails, with a Python traceback in the body.
  :func:`parse_isprint` looks for that and raises rather than returning an
  empty table.
* It answers a window that matches nothing with the single line "No records
  were selected with the filters above", also HTTP 200.  That is an answer,
  and :func:`parse_isprint` returns an empty table for it.

Identity.  Madrigal requires the user's name, email and affiliation on every
data request; the listings of instruments, experiments and files do not need
them.  They are read from ``TSS_MADRIGAL_NAME``, ``TSS_MADRIGAL_EMAIL`` and
``TSS_MADRIGAL_AFFILIATION`` and never written anywhere the repository keeps:
the URL recorded in a cache sidecar has them removed, because sidecars are
committed.  A data request without them raises :class:`MadrigalIdentityError`
with instructions rather than sending an anonymous request.
"""

from __future__ import annotations

import html
import os
import re
from collections.abc import Sequence
from dataclasses import dataclass
from datetime import date, datetime, timedelta, timezone
from pathlib import Path
from typing import Any
from urllib.parse import quote, urlencode

from tss.data import cache
from tss.data.base import FetchError, register_source

__all__ = [
    "SITE",
    "SOURCE",
    "INSTRUMENTS",
    "KINDAT",
    "MadrigalIdentityError",
    "MadrigalServiceError",
    "Instrument",
    "Experiment",
    "ExperimentFile",
    "IsprintTable",
    "identity",
    "parse_instruments",
    "parse_experiments",
    "parse_files",
    "parse_isprint",
    "build_isprint_params",
    "redact_url",
    "list_instruments",
    "list_experiments",
    "list_files",
    "isprint",
    "profiles",
    "fetch_madrigal_sample",
]

SITE = "https://cedar.openmadrigal.org/"
SOURCE = "madrigal"
MIN_INTERVAL_S = 2.0

INSTRUMENTS: dict[str, int] = {
    "jicamarca": 10,
    "arecibo_linefeed": 20,
    "arecibo_gregorian": 21,
    "millstone": 30,
    "millstone_steerable": 31,
    "millstone_zenith": 32,
    "poker_flat": 61,
    "eiscat_tromso_uhf": 72,
    "eiscat_tromso_vhf": 74,
    "resolute_north": 91,
    "eiscat_svalbard": 95,
    "gnss_tec": 8000,
    "dmsp": 8100,
}
"""Madrigal instrument codes for the instruments this project uses.

Taken from ``getInstrumentsService.py`` on 2026-09-23.  The full list has about
300 entries; :func:`list_instruments` returns it.
"""

KINDAT: dict[str, int] = {
    "mlh_zenith_single_pulse": 3420,
    "mlh_zenith_alternating_code": 3425,
    "mlh_combined_basic": 3410,
    "mlh_gridded": 13300,
    "dmsp_f16_1s": 10116,
    "dmsp_f17_1s": 10117,
    "dmsp_f18_1s": 10118,
}
"""Kind-of-data codes seen in the listings, for selecting a file within an
experiment.  Millstone's zenith single-pulse product is the one with the
longest topside altitude coverage."""

_NO_RECORDS = "No records were selected with the filters above"
"""What isprintService.py returns, with HTTP 200, when the window matches nothing."""

_IDENTITY_ENV = ("TSS_MADRIGAL_NAME", "TSS_MADRIGAL_EMAIL", "TSS_MADRIGAL_AFFILIATION")
_IDENTITY_KEYS = ("user_fullname", "user_email", "user_affiliation")


class MadrigalIdentityError(FetchError):
    """A data request was attempted without the identity Madrigal requires."""


class MadrigalServiceError(FetchError):
    """A Madrigal service answered with an error in the body, whatever the status."""


# ---------------------------------------------------------------------------
# Records
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class Instrument:
    name: str
    code: int
    mnemonic: str
    lat_deg: float
    lon_deg: float
    alt_km: float
    category: str


@dataclass(frozen=True)
class Experiment:
    """One row of ``getExperimentsService.py``.

    ``url`` is the experiment's page on the site, and ``path`` the part of it
    after ``madtoc/``, which is what identifies the experiment directory.
    ``is_local`` is False for an experiment held at another Madrigal site;
    such an experiment's files cannot be listed from here.
    """

    id: int
    url: str
    name: str
    site_id: int
    site_name: str
    inst_code: int
    inst_name: str
    start: datetime
    end: datetime
    is_local: bool
    pi: str
    pi_email: str

    @property
    def path(self) -> str:
        return self.url.split("madtoc/", 1)[-1]


@dataclass(frozen=True)
class ExperimentFile:
    """One row of ``getExperimentFilesService.py``.

    ``category`` is 1 default, 2 variant, 3 history, 4 real-time.  ``name`` is
    the full server path, which is what ``isprintService.py`` takes.
    """

    name: str
    kindat: int
    description: str
    category: int
    status: str
    permission: int
    doi: str

    @property
    def basename(self) -> str:
        return self.name.rsplit("/", 1)[-1]


@dataclass(frozen=True)
class IsprintTable:
    """An isprint result: one float array per requested parameter.

    ``columns`` preserves the order the parameters were requested in, which is
    the only thing that names the columns, since the service writes no header
    line.  Missing values are NaN.  ``times`` is built from the YEAR, MONTH,
    DAY, HOUR, MIN and SEC columns when they were requested.
    """

    columns: tuple[str, ...]
    data: dict[str, Any]
    times: Any = None
    file: str = ""

    def __len__(self) -> int:
        return len(next(iter(self.data.values()))) if self.data else 0

    def __getitem__(self, key: str) -> Any:
        return self.data[key.upper()]


# ---------------------------------------------------------------------------
# Identity and URLs
# ---------------------------------------------------------------------------


def identity() -> dict[str, str]:
    """The Madrigal identity from the environment, or an explanatory error."""
    values = [os.environ.get(k, "").strip() for k in _IDENTITY_ENV]
    missing = [k for k, v in zip(_IDENTITY_ENV, values, strict=True) if not v]
    if missing:
        raise MadrigalIdentityError(
            "Madrigal requires a name, email and affiliation on every data request. "
            f"Set {', '.join(missing)} in the environment (not in any file in the "
            "repository) and try again."
        )
    return dict(zip(_IDENTITY_KEYS, values, strict=True))


def redact_url(url: str) -> str:
    """``url`` with the identity parameters' values replaced by ``[redacted]``."""
    return re.sub(
        r"(user_fullname|user_email|user_affiliation)=[^&]*", r"\1=[redacted]", url
    )


def _encode(pairs: Sequence[tuple[str, str]]) -> str:
    """URL-encode pairs the way Madrigal's services expect.

    Madrigal treats a literal ``+`` as a space in names, and ``madrigalWeb``
    encodes both spaces and plus signs as ``%2B`` in the identity fields.  The
    filter string keeps ``:`` unescaped, as ``madrigalWeb`` does.
    """
    out = []
    for k, v in pairs:
        if k in ("user_fullname", "user_affiliation"):
            out.append(f"{k}={quote(v.replace(' ', '+'), safe='')}")
        elif k == "user_email":
            out.append(f"{k}={quote(v, safe='@')}")
        elif k == "filters":
            out.append(f"{k}={quote(v, safe=':')}")
        else:
            out.append(urlencode([(k, v)], safe="/,"))
    return "&".join(out)


def _as_utc(value: str | date | datetime) -> datetime:
    if isinstance(value, datetime):
        return value.astimezone(timezone.utc) if value.tzinfo else value.replace(tzinfo=timezone.utc)
    if isinstance(value, date):
        return datetime(value.year, value.month, value.day, tzinfo=timezone.utc)
    if isinstance(value, str):
        dt = datetime.fromisoformat(value.strip().replace("Z", "+00:00"))
        return dt.astimezone(timezone.utc) if dt.tzinfo else dt.replace(tzinfo=timezone.utc)
    raise TypeError(f"cannot read a {type(value).__name__} as a time")


# ---------------------------------------------------------------------------
# Parsers
# ---------------------------------------------------------------------------


def _check_service_error(text: str, service: str) -> None:
    low = text.lower()
    if "traceback" in low or "error occurred" in low or re.search(r"^\w*error:", text, re.M) \
            or "raise " in text or "<html" in low[:200]:
        msg = " ".join(html.unescape(text).split())[:300]
        raise MadrigalServiceError(f"{service} returned an error: {msg}")


def parse_instruments(text: str | bytes) -> list[Instrument]:
    """Parse ``getInstrumentsService.py``: name, code, mnemonic, lat, lon, alt, category, ..."""
    text = text.decode("utf-8", "replace") if isinstance(text, bytes) else text
    _check_service_error(text, "getInstrumentsService.py")
    out = []
    for line in text.splitlines():
        items = line.split(",")
        if len(items) < 7:
            continue
        try:
            out.append(Instrument(
                name=items[0], code=int(items[1]), mnemonic=items[2],
                lat_deg=float(items[3]), lon_deg=float(items[4]), alt_km=float(items[5]),
                category=items[6],
            ))
        except ValueError:
            continue
    return out


def parse_experiments(text: str | bytes) -> list[Experiment]:
    """Parse ``getExperimentsService.py``.

    Field order, as ``madrigalWeb`` reads it and as checked against a live
    response on 2026-09-23: id, url, name, site id, site name, instrument code,
    instrument name, start y m d h m s, end y m d h m s, local flag, PI, PI
    email, then a timestamp and an access code.  Lines with fewer than 20
    fields are skipped, as ``madrigalWeb`` skips them.
    """
    text = text.decode("utf-8", "replace") if isinstance(text, bytes) else text
    _check_service_error(text, "getExperimentsService.py")
    out = []
    for line in text.splitlines():
        it = line.split(",")
        if len(it) < 20:
            continue
        s = [int(x) for x in it[7:13]]
        e = [int(x) for x in it[13:19]]
        out.append(Experiment(
            id=int(it[0]), url=it[1], name=it[2], site_id=int(it[3]), site_name=it[4],
            inst_code=int(it[5]), inst_name=it[6],
            start=datetime(*s, tzinfo=timezone.utc), end=datetime(*e, tzinfo=timezone.utc),
            is_local=it[19].strip() == "1",
            pi=it[20] if len(it) > 20 else "", pi_email=it[21] if len(it) > 21 else "",
        ))
    return out


def parse_files(text: str | bytes, *, include_non_default: bool = False) -> list[ExperimentFile]:
    """Parse ``getExperimentFilesService.py``.

    By default this keeps what ``madrigalWeb`` keeps: the default files
    (category 1) if the experiment has any, otherwise the real-time ones
    (category 4).  Preliminary history versions are dropped unless
    ``include_non_default`` is set, because scoring against a superseded
    calibration is an easy mistake to make silently.
    """
    text = text.decode("utf-8", "replace") if isinstance(text, bytes) else text
    _check_service_error(text, "getExperimentFilesService.py")
    rows = []
    for line in text.splitlines():
        it = line.split(",")
        if len(it) < 6:
            continue
        rows.append(ExperimentFile(
            name=it[0], kindat=int(it[1]), description=it[2], category=int(it[3]),
            status=it[4], permission=int(it[5]), doi=it[6] if len(it) > 6 else "",
        ))
    if include_non_default:
        return rows
    has_default = any(r.category == 1 for r in rows)
    keep = 1 if has_default else 4
    return [r for r in rows if r.category == keep]


def parse_isprint(text: str | bytes, parms: Sequence[str], *, file: str = "") -> IsprintTable:
    """Parse an ``isprintService.py`` body for the requested ``parms``.

    The body has no header line, so the columns are named by ``parms`` in
    order, and a row whose field count differs from ``len(parms)`` is an error
    rather than something to guess about.  ``nan`` and ``missing`` are read as
    NaN.
    """
    import numpy as np

    text = text.decode("utf-8", "replace") if isinstance(text, bytes) else text
    _check_service_error(text, "isprintService.py")
    names = tuple(p.upper() for p in parms)
    if _NO_RECORDS in text:
        # An empty selection is an answer, not an error: the file exists and
        # simply has nothing in the requested window.
        empty = np.array([], dtype=float)
        return IsprintTable(columns=names, data={n: empty.copy() for n in names},
                            times=np.array([], dtype="datetime64[s]"), file=file)
    rows = []
    for line in text.splitlines():
        parts = line.split()
        if not parts:
            continue
        if len(parts) != len(names):
            raise MadrigalServiceError(
                f"isprint row has {len(parts)} fields but {len(names)} parameters were "
                f"requested: {line.strip()[:120]!r}"
            )
        rows.append([float("nan") if p.lower() in ("nan", "missing") else float(p) for p in parts])
    arr = np.array(rows, dtype=float).reshape(-1, len(names))
    data = {n: arr[:, i] for i, n in enumerate(names)}

    times = None
    need = ("YEAR", "MONTH", "DAY", "HOUR", "MIN", "SEC")
    if all(k in data for k in need) and len(arr):
        times = np.array(
            [
                np.datetime64(
                    f"{int(y):04d}-{int(mo):02d}-{int(d):02d}T{int(h):02d}:{int(mi):02d}:{int(s):02d}",
                    "s",
                )
                for y, mo, d, h, mi, s in zip(*(data[k] for k in need), strict=True)
            ]
        )
    return IsprintTable(columns=names, data=data, times=times, file=file)


def build_isprint_params(
    file: str,
    parms: Sequence[str],
    start: str | date | datetime | None = None,
    end: str | date | datetime | None = None,
    *,
    extra_filters: str = "",
    ident: dict[str, str] | None = None,
) -> list[tuple[str, str]]:
    """The query pairs for one isprint call.

    ``parms`` become one ``parms=`` pair each.  ``start`` and ``end`` become the
    isprint ``date1 time1 date2 time2`` filter; ``extra_filters`` is appended
    verbatim in isprint's own syntax, for example ``"filter=gdalt,300,900"``.
    """
    if not parms:
        raise ValueError("isprint needs at least one parameter")
    filt = []
    if start is not None:
        t0 = _as_utc(start)
        filt.append(f"date1={t0:%m/%d/%Y} time1={t0:%H:%M:%S}")
    if end is not None:
        t1 = _as_utc(end)
        filt.append(f"date2={t1:%m/%d/%Y} time2={t1:%H:%M:%S}")
    if extra_filters:
        filt.append(extra_filters.strip())
    pairs = [("file", file)]
    pairs += [("parms", p.strip().upper()) for p in parms]
    pairs.append(("filters", " ".join(filt)))
    if ident:
        pairs += [(k, ident[k]) for k in _IDENTITY_KEYS]
    return pairs


# ---------------------------------------------------------------------------
# Service calls
# ---------------------------------------------------------------------------


def _get_cached(
    service: str,
    pairs: Sequence[tuple[str, str]],
    *,
    suffix: str,
    max_age_days: float,
    force: bool,
    timeout: float = 120.0,
) -> Path:
    """GET one service call through the cache, recording a redacted URL."""
    url = f"{SITE}{service}?{_encode(pairs)}"
    public = redact_url(url)
    # The cache key is the redacted URL, so the identity never reaches a filename
    # and two people fetching the same data share a cache entry.
    path = cache.cached_path(SOURCE, public, suffix)
    if not force and cache.is_fresh(path, max_age_days) and cache.verify_cached(path):
        return path
    response = cache.polite_get(url, timeout=timeout, min_interval_s=MIN_INTERVAL_S)
    body = response.content
    _check_service_error(body.decode("utf-8", "replace"), service)
    path.write_bytes(body)
    cache.record_meta(path, url=public, source=SOURCE, extra={"service": service})
    return path


def list_instruments(*, force: bool = False) -> list[Instrument]:
    """Every instrument the site knows.  Needs no identity."""
    path = _get_cached("getInstrumentsService.py", [], suffix=".txt", max_age_days=90.0, force=force)
    return parse_instruments(path.read_bytes())


def list_experiments(
    instrument: int | str,
    start: str | date | datetime,
    end: str | date | datetime,
    *,
    force: bool = False,
) -> list[Experiment]:
    """Experiments for one instrument overlapping a UTC window.  Needs no identity.

    ``instrument`` is a code or a key of :data:`INSTRUMENTS`.
    """
    code = INSTRUMENTS[instrument] if isinstance(instrument, str) else int(instrument)
    t0, t1 = _as_utc(start), _as_utc(end)
    if t1 <= t0:
        raise ValueError(f"end {t1.isoformat()} is not after start {t0.isoformat()}")
    pairs = [("code", str(code))]
    for tag, t in (("start", t0), ("end", t1)):
        pairs += [
            (f"{tag}year", str(t.year)), (f"{tag}month", str(t.month)),
            (f"{tag}day", str(t.day)), (f"{tag}hour", str(t.hour)),
            (f"{tag}min", str(t.minute)), (f"{tag}sec", str(t.second)),
        ]
    pairs.append(("local", "1"))
    path = _get_cached("getExperimentsService.py", pairs, suffix=".txt", max_age_days=30.0, force=force)
    return parse_experiments(path.read_bytes())


def list_files(
    experiment: int | Experiment, *, include_non_default: bool = False, force: bool = False
) -> list[ExperimentFile]:
    """The files of one experiment.  Needs no identity."""
    exp_id = experiment.id if isinstance(experiment, Experiment) else int(experiment)
    path = _get_cached(
        "getExperimentFilesService.py", [("id", str(exp_id))],
        suffix=".txt", max_age_days=30.0, force=force,
    )
    return parse_files(path.read_bytes(), include_non_default=include_non_default)


def isprint(
    file: str | ExperimentFile,
    parms: Sequence[str],
    start: str | date | datetime | None = None,
    end: str | date | datetime | None = None,
    *,
    extra_filters: str = "",
    force: bool = False,
) -> tuple[IsprintTable, Path]:
    """Extract selected parameters from one file over a window.  Needs identity.

    This is the preferred way to take data from Madrigal: the server does the
    selection and returns a few kilobytes of text instead of a file of tens of
    megabytes.  Results are cached as served, keyed on the redacted request.
    """
    name = file.name if isinstance(file, ExperimentFile) else str(file)
    pairs = build_isprint_params(name, parms, start, end, extra_filters=extra_filters, ident=identity())
    path = _get_cached("isprintService.py", pairs, suffix=".txt", max_age_days=float("inf"),
                       force=force, timeout=300.0)
    return parse_isprint(path.read_bytes(), parms, file=name), path


def profiles(
    table: IsprintTable,
    *,
    alt: str = "GDALT",
    ne: str = "NE",
    dne: str | None = "DNE",
    max_rel_err: float | None = None,
) -> list[dict[str, Any]]:
    """Group an ISR isprint table into one profile per timestamp.

    Each profile is a dict with ``time``, ``alt_km``, ``ne_m3`` and, if
    ``dne`` was requested, ``dne_m3``.  Points with NaN density are dropped.
    ``max_rel_err`` also drops points whose relative error ``dne/ne`` exceeds
    it.  That is necessary and not sufficient.  In the Millstone fixture the
    22:01:03 profile has a point at 645 km reading 4.7e11 m^-3, about three
    times its neighbours, with a stated error of only 2 %; an error-bar filter
    keeps it.  Topside ISR points need a physical plausibility check as well
    (monotone decrease above hmF2 within the noise, for instance) before they
    are used as truth, and that check belongs to the study, not to this reader.
    Altitudes are sorted ascending.
    """
    import numpy as np

    if table.times is None:
        raise ValueError("profiles() needs YEAR, MONTH, DAY, HOUR, MIN and SEC in the table")
    t = table.times
    out = []
    for stamp in np.unique(t):
        m = t == stamp
        a, n = table[alt][m], table[ne][m]
        good = np.isfinite(n)
        d = table[dne][m] if dne and dne.upper() in table.data else None
        if d is not None and max_rel_err is not None:
            with np.errstate(divide="ignore", invalid="ignore"):
                good &= np.isfinite(d) & (d / n <= max_rel_err)
        order = np.argsort(a[good])
        prof = {"time": stamp, "alt_km": a[good][order], "ne_m3": n[good][order]}
        if d is not None:
            prof["dne_m3"] = d[good][order]
        out.append(prof)
    return out


@register_source(
    name=SOURCE,
    description="CEDAR Madrigal: ISR topside Ne profiles and DMSP in-situ Ne",
    url=SITE,
    licence=(
        "CEDAR Madrigal data are openly available; the CEDAR Rules of the Road apply, "
        "requiring acknowledgement of the instrument and its funding and, where "
        "appropriate, contact with the instrument PI"
    ),
    citation=(
        "Data obtained from the CEDAR Madrigal database, MIT Haystack Observatory "
        "(https://cedar.openmadrigal.org/). Cite each file by the DOI in its "
        "getExperimentFilesService listing, and the instrument reference its PI supplies."
    ),
    cadence="archive, growing; experiment files move from preliminary to final",
    work_packages=("WP1", "WP8"),
    acknowledgement=(
        "Acknowledge the Madrigal database and each instrument used, e.g. Millstone Hill "
        "(MIT Haystack Observatory, NSF) and DMSP SSIES (Boston College, AFRL)."
    ),
    notes=(
        "Data requests need TSS_MADRIGAL_NAME, TSS_MADRIGAL_EMAIL and "
        "TSS_MADRIGAL_AFFILIATION; these are never written to the repository. ISR sites "
        "and campaign days are sparse. isprintService answers 200 on failure; errors are "
        "detected in the body."
    ),
    extra={"site": SITE, "instruments": dict(INSTRUMENTS)},
)
def fetch_madrigal_sample(
    *,
    force: bool = False,
    limit: int | None = None,
    file: str | None = None,
    parms: Sequence[str] | None = None,
    start: str | None = None,
    end: str | None = None,
    extra_filters: str = "",
) -> list[Path]:
    """Fetch one isprint extract and cache it.

    Without arguments this takes ten minutes of the Millstone Hill zenith
    single-pulse product on 2023-09-20 from 22:00 UT, which is the fixture's
    source: a topside profile to 800 km.  ``limit`` is not meaningful for a
    single extract and is ignored.
    """
    del limit
    if file is None:
        file = "/opt/openmadrigal/madroot/experiments5/2023/mlh/20sep23/mlh230920i.004.hdf5"
        start = start or "2023-09-20T22:00:00"
        end = end or "2023-09-20T22:10:00"
    if parms is None:
        parms = ("YEAR", "MONTH", "DAY", "HOUR", "MIN", "SEC", "GDLAT", "GLON", "GDALT", "NE", "DNE")
    if start is None or end is None:
        t1 = datetime.now(timezone.utc).replace(minute=0, second=0, microsecond=0)
        start = start or (t1 - timedelta(hours=1)).isoformat()
        end = end or t1.isoformat()
    _table, path = isprint(file, parms, start, end, extra_filters=extra_filters, force=force)
    return [path]
