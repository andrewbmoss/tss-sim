"""Layer: data
WP: WP1, WP8
Claims: -
Status: implemented
Requirements: INS-01

NASA SPDF Alouette/ISIS topside sounder archive.

The Alouette and ISIS spacecraft flew topside sounders from 1962 to 1984 and
are the only large body of real topside ionograms and real topside electron
density profiles in existence.  Two products matter here.

The hand-scaled electron density profiles, held under
``topside_sounder/crc_ne_profile_ascii/``, were scaled at the Communications
Research Centre in Ottawa by Jackson's method.  ``isis1_ne.asc`` is 8.4 MB and
holds ISIS-1 profiles from 1969 to 1980 in a fixed-width ASCII format.  The profile count is not asserted here because this repository has never had the file to count; read_ne_ascii reports it at read time.  Each profile runs downward
from the spacecraft, so its first sample is the local density at the satellite.
That single number is what sets the lower band edge of INS-01, which is why
this file is a direct test of the requirement and not an indirect one.  WP1
uses the same profiles as truth to check the L0 models against.

The digitised ionograms under the same archive are what WP8 and V4 need in
order to run the retrieval against a real measurement rather than a simulated
one.  Those are large binary products and this module does not fetch them
automatically; :data:`IONOGRAM_DIRS` records where they live and the status
page explains which segments have been digitised.

Rights.  The data are a US Government work and are in the public domain.  The
archive nonetheless asks that it be acknowledged, and the CRC scaling effort
should be credited, so the source descriptor carries an acknowledgement string
which the run manifest copies.

Offline use.  The author already holds ``isis1_ne.asc`` locally.  The fetcher
therefore looks for a local copy before it considers downloading, and
``TSS_ISIS_NE`` names an explicit path.  Nothing here downloads anything unless
the local copies are all missing.
"""

from __future__ import annotations

import os
from collections.abc import Iterator, Sequence
from dataclasses import dataclass, field
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any

import numpy as np

from tss.data import cache
from tss.data.base import FetchError, register_source

__all__ = [
    "BASE_URL",
    "NE_DIR_URL",
    "STATUS_URL",
    "IONOGRAM_DIRS",
    "NE_FILES",
    "HDR",
    "PLASMA",
    "Profile",
    "fetch_isis_ne",
    "read_ne_ascii",
    "iter_profiles",
    "local_ne_path",
]


BASE_URL = "https://spdf.gsfc.nasa.gov/pub/data/isis/"
NE_DIR_URL = BASE_URL + "topside_sounder/crc_ne_profile_ascii/"
STATUS_URL = "https://spdf.gsfc.nasa.gov/isis/isis-status.html"

IONOGRAM_DIRS: tuple[str, ...] = (
    BASE_URL + "topside_sounder/",
    BASE_URL + "topside_sounder/ionogram_images/",
)
"""Where the digitised ionogram products sit.

Listed rather than fetched.  These are large and the archive's coverage is
uneven by segment, so WP8 picks specific passes by hand after reading the status
page.  Guessing at filenames here would be inventing a layout.
"""

NE_FILES: dict[str, str] = {
    "isis1": "isis1_ne.asc",
    "isis2": "isis2_ne.asc",
    "alouette1": "alouette1_ne.asc",
    "alouette2": "alouette2_ne.asc",
}
"""Per-spacecraft hand-scaled profile files in ``crc_ne_profile_ascii/``.

Only ``isis1`` is confirmed present and used; the other three names follow the
archive's own convention and the fetcher reports a clear error rather than a
silent empty result if one of them is not there.
"""

ISIS_MAX_AGE_DAYS = float("inf")
"""The archive is closed.  A cached copy never goes stale."""


# ---------------------------------------------------------------------------
# Fixed-width record layout, from the archive's 00readme.txt
# ---------------------------------------------------------------------------
#
# This block and the two helpers below are carried over unchanged in behaviour
# from the author's crc_read.py.  The format is a Fortran fixed-width record
# and the field widths are not negotiable, so it is reproduced rather than
# re-derived.
#
#     line 1                M            I5     number of records
#     per record header     itype        I1     spacecraft id
#                           iquality     I2     0 best .. 10 worst
#                           date         I6     yyddd
#                           time         I6     hhmmss
#                           long         F6.1   deg
#                           lat          F6.2   deg
#                           dip          F6.2   deg (IGRF)
#                           L            F5.2
#                           xhi          F6.1   solar zenith at 100 km
#                           Rz12         I3
#                           IG12         I3
#                           rtec         F6.2
#                           lnNmF2_iri   F7.4   IRI F2 peak, ln(cm^-3)
#                           hmF2_iri     F5.1   IRI F2 peak height, km
#                           TEC_iri      F6.2
#                           n            I3     number of profile points
#     heights               14I5   h_km * 10
#     densities             10I7   ln(Ne / cm^-3) * 1e5

HDR: list[tuple[str, int, type]] = [
    ("itype", 1, int), ("iquality", 2, int), ("date", 6, int),
    ("time", 6, int), ("lon", 6, float), ("lat", 6, float),
    ("dip", 6, float), ("L", 5, float), ("xhi", 6, float),
    ("rz12", 3, int), ("ig12", 3, int), ("rtec", 6, float),
    ("ln_nmf2_iri", 7, float), ("hmf2_iri", 5, float),
    ("tec_iri", 6, float), ("n", 3, int),
]

PLASMA = 8.98
"""fN[Hz] = 8.98 * sqrt(Ne[m^-3]).

Kept at the value crc_read.py used so that results reproduce bit for bit.  The
package constant in :mod:`tss.constants` carries more digits; the difference is
below the scaling precision of the archive and is not worth a discrepancy.
"""


def _fixed(line: str, spec: Sequence[tuple[str, int, type]]) -> dict[str, Any]:
    """Split one fixed-width header line.  Unparseable fields become None."""
    out: dict[str, Any] = {}
    i = 0
    for name, w, cast in spec:
        s = line[i:i + w].strip()
        i += w
        try:
            out[name] = cast(s) if s else None
        except ValueError:
            out[name] = None
    return out


def _packed(
    lines: Sequence[str], count: int, width: int, start: int
) -> tuple[list[int], int]:
    """Read ``count`` fixed-width integers running across wrapped lines."""
    vals: list[int] = []
    li, ci = start, 0
    while len(vals) < count:
        line = lines[li]
        if ci >= len(line.rstrip("\n")):
            li += 1
            ci = 0
            continue
        chunk = line[ci:ci + width]
        ci += width
        s = chunk.strip()
        if s == "":
            continue
        vals.append(int(s))
    if ci >= len(lines[li].rstrip("\n")):
        li += 1
    else:
        li += 1          # a value block always ends its line
    return vals, li


def _read_records(
    path: str | Path,
    *,
    want_profiles: bool = True,
    limit: int | None = None,
) -> tuple[int, list[dict[str, Any]]]:
    """The crc_read.py record loop, unchanged in what it computes.

    Returns the declared record count from line 1 and the parsed records.  The
    declared count and the parsed count are returned separately on purpose: they
    differ when ``limit`` is set, and a mismatch without a limit is a signal the
    file is truncated.
    """
    with open(path) as f:
        lines = f.readlines()
    m = int(lines[0].strip())
    recs: list[dict[str, Any]] = []
    li = 1
    while li < len(lines) and (limit is None or len(recs) < limit):
        if not lines[li].strip():
            li += 1
            continue
        h = _fixed(lines[li], HDR)
        li += 1
        n = h["n"]
        if not n or n < 2:
            continue
        hts, li = _packed(lines, n, 5, li)
        den, li = _packed(lines, n, 7, li)
        h["h_km"] = np.asarray(hts, float) / 10.0
        h["ne_cm3"] = np.exp(np.asarray(den, float) / 1e5)
        h["h_sat_km"] = h["h_km"][0]
        h["ne_sat_cm3"] = h["ne_cm3"][0]
        h["fn_sat_mhz"] = PLASMA * np.sqrt(h["ne_sat_cm3"] * 1e6) / 1e6
        h["nmf2_iri_cm3"] = np.exp(h["ln_nmf2_iri"]) if h["ln_nmf2_iri"] else np.nan
        h["fof2_iri_mhz"] = PLASMA * np.sqrt(h["nmf2_iri_cm3"] * 1e6) / 1e6
        # observed peak: deepest point of the profile (profiles stop at/near hmF2)
        h["ne_bot_cm3"] = h["ne_cm3"][-1]
        h["h_bot_km"] = h["h_km"][-1]
        h["fo_bot_mhz"] = PLASMA * np.sqrt(h["ne_bot_cm3"] * 1e6) / 1e6
        yy = h["date"] // 1000
        h["year"] = 1900 + yy if yy >= 50 else 2000 + yy
        h["doy"] = h["date"] % 1000
        hh = h["time"] // 10000
        h["ut_hour"] = hh + (h["time"] // 100 % 100) / 60.0
        h["lt_hour"] = (h["ut_hour"] + h["lon"] / 15.0) % 24.0
        if not want_profiles:
            h.pop("h_km")
            h.pop("ne_cm3")
        recs.append(h)
    return m, recs


def _epoch(rec: dict[str, Any]) -> datetime | None:
    """UTC timestamp from the record's year, day of year and hhmmss.

    Derived, not read: the file carries ``yyddd`` and ``hhmmss`` separately and
    this assembles them.  Returns None rather than raising if either is absent,
    because a record with a damaged date is still a usable profile for the
    altitude statistics.
    """
    try:
        year = int(rec["year"])
        doy = int(rec["doy"])
        t = int(rec["time"])
    except (KeyError, TypeError, ValueError):
        return None
    hh, mm, ss = t // 10000, (t // 100) % 100, t % 100
    try:
        base = datetime(year, 1, 1, tzinfo=timezone.utc)
        return base + timedelta(days=doy - 1, hours=hh, minutes=mm, seconds=ss)
    except (ValueError, OverflowError):
        return None


# ---------------------------------------------------------------------------
# Public reader
# ---------------------------------------------------------------------------


def read_ne_ascii(
    path: str | Path,
    *,
    want_profiles: bool = True,
    limit: int | None = None,
    strict: bool = False,
) -> Any:
    """Read a CRC hand-scaled Ne profile file into a pandas DataFrame.

    One row per observation.  Column semantics are exactly those of the author's
    crc_read.py, which in turn follow the archive's own 00readme.txt:

    ==================  =========================================================
    ``itype``           spacecraft id
    ``iquality``        scaling quality, 0 best to 10 worst
    ``date``            yyddd as stored
    ``time``            hhmmss as stored
    ``lon``, ``lat``    degrees
    ``dip``             IGRF dip angle, degrees
    ``L``               McIlwain L
    ``xhi``             solar zenith angle at 100 km, degrees
    ``rz12``, ``ig12``  the twelve-month indices
    ``rtec``            scaled TEC above the profile's lowest point
    ``ln_nmf2_iri``     IRI F2 peak density, ln(cm^-3)
    ``hmf2_iri``        IRI F2 peak height, km
    ``tec_iri``         IRI TEC
    ``n``               number of profile samples
    ``h_km``            profile heights, km, downward from the satellite
    ``ne_cm3``          profile densities, cm^-3, parallel to ``h_km``
    ``h_sat_km``        satellite height, which is ``h_km[0]``
    ``ne_sat_cm3``      local density at the satellite, ``ne_cm3[0]``
    ``fn_sat_mhz``      local plasma frequency at the satellite, MHz
    ``nmf2_iri_cm3``    IRI F2 peak density, cm^-3
    ``fof2_iri_mhz``    IRI foF2, MHz
    ``ne_bot_cm3``      density at the profile's deepest sample
    ``h_bot_km``        height of that sample
    ``fo_bot_mhz``      plasma frequency there, MHz
    ``year``            four-digit year, with the 50 pivot for the two-digit field
    ``doy``             day of year
    ``ut_hour``         UT hour, fractional to the minute
    ``lt_hour``         local solar time, UT plus lon/15, wrapped to 24
    ``epoch``           derived UTC timestamp, added by this module
    ==================  =========================================================

    ``h_km`` and ``ne_cm3`` are object columns holding numpy arrays, because the
    profiles have different lengths.  Pass ``want_profiles=False`` to drop them
    and get a purely scalar table, which is much cheaper for the whole
    record file.

    ``strict=True`` raises if the number of records parsed does not match the
    count the file declares on its first line and no ``limit`` was given.  That
    is off by default because a truncated file is still readable and saying so
    is more useful than refusing.

    pandas is imported inside the function so that importing :mod:`tss.data`
    does not pull it in.
    """
    import pandas as pd

    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(
            f"no ISIS Ne profile file at {p}.  Run `tss fetch spdf_isis_ne`, or "
            f"set TSS_ISIS_NE to a local copy of a file from {NE_DIR_URL}"
        )

    declared, recs = _read_records(p, want_profiles=want_profiles, limit=limit)
    if strict and limit is None and declared != len(recs):
        raise FetchError(
            f"{p} declares {declared} records but {len(recs)} parsed; the file "
            "looks truncated"
        )

    for rec in recs:
        rec["epoch"] = _epoch(rec)

    frame = pd.DataFrame.from_records(recs)
    frame.attrs["source_file"] = str(p)
    frame.attrs["declared_records"] = declared
    frame.attrs["archive"] = NE_DIR_URL
    meta = cache.read_meta(p)
    if meta:
        frame.attrs["sha256"] = meta.get("sha256")
        frame.attrs["retrieved_utc"] = meta.get("retrieved_utc")
    return frame


@dataclass(frozen=True)
class Profile:
    """One scaled topside observation.

    ``alt_km`` and ``ne_cm3`` run downward from the spacecraft, so ``alt_km[0]``
    is ``sat_alt_km`` and ``ne_cm3[0]`` is the local density there.  This is the
    shape :class:`tss.interfaces.IonosphereState` wants reversed, and
    :meth:`ascending` does that without mutating anything.
    """

    epoch: datetime | None
    lat_deg: float
    lon_deg: float
    sat_alt_km: float
    alt_km: np.ndarray
    ne_cm3: np.ndarray
    quality: int | None = None
    spacecraft: int | None = None
    extra: dict[str, Any] = field(default_factory=dict)

    @property
    def ne_m3(self) -> np.ndarray:
        """Density in SI, which is what every tier above L0 expects."""
        return self.ne_cm3 * 1.0e6

    @property
    def alt_m(self) -> np.ndarray:
        return self.alt_km * 1.0e3

    @property
    def plasma_freq_sat_mhz(self) -> float:
        """Local plasma frequency at the spacecraft, the INS-01 quantity."""
        return float(PLASMA * np.sqrt(self.ne_cm3[0] * 1.0e6) / 1.0e6)

    def ascending(self) -> tuple[np.ndarray, np.ndarray]:
        """``(alt_m, ne_m3)`` sorted by increasing altitude."""
        order = np.argsort(self.alt_m)
        return self.alt_m[order], self.ne_m3[order]


def iter_profiles(
    path: str | Path | None = None,
    *,
    limit: int | None = None,
    max_quality: int | None = None,
) -> Iterator[Profile]:
    """Yield one :class:`Profile` per observation.

    ``path`` defaults to whatever :func:`local_ne_path` finds.  ``max_quality``
    keeps only records whose ``iquality`` is at or below the given value, where
    0 is the best scaling and 10 the worst.  Filtering happens after ``limit``
    is applied to the parse, so a limit is a limit on records read and not on
    records yielded; that keeps the read cost predictable on the full file.

    This is a generator over an already-parsed list rather than a streaming
    parser.  The whole file is 8.4 MB and reading it takes under a second, so
    streaming would buy complexity and nothing else.
    """
    p = Path(path) if path is not None else local_ne_path()
    if p is None:
        raise FileNotFoundError(
            "no local ISIS Ne profile file found; run `tss fetch spdf_isis_ne` "
            "or set TSS_ISIS_NE"
        )
    _declared, recs = _read_records(p, want_profiles=True, limit=limit)
    for rec in recs:
        if max_quality is not None:
            q = rec.get("iquality")
            if q is None or q > max_quality:
                continue
        yield Profile(
            epoch=_epoch(rec),
            lat_deg=float(rec["lat"]) if rec.get("lat") is not None else float("nan"),
            lon_deg=float(rec["lon"]) if rec.get("lon") is not None else float("nan"),
            sat_alt_km=float(rec["h_sat_km"]),
            alt_km=rec["h_km"],
            ne_cm3=rec["ne_cm3"],
            quality=rec.get("iquality"),
            spacecraft=rec.get("itype"),
            extra={
                "dip_deg": rec.get("dip"),
                "L": rec.get("L"),
                "solar_zenith_deg": rec.get("xhi"),
                "rz12": rec.get("rz12"),
                "ig12": rec.get("ig12"),
                "ut_hour": rec.get("ut_hour"),
                "lt_hour": rec.get("lt_hour"),
                "hmf2_iri_km": rec.get("hmf2_iri"),
                "fof2_iri_mhz": rec.get("fof2_iri_mhz"),
                "source_file": str(p),
            },
        )


# ---------------------------------------------------------------------------
# Locating and fetching
# ---------------------------------------------------------------------------


def local_ne_path(spacecraft: str = "isis1") -> Path | None:
    """Find an already-present copy of the Ne profile file, or None.

    Search order: ``TSS_ISIS_NE``, then the cache directory, then the places the
    author's own working tree keeps it.  The last of those is a convenience for
    one machine and is checked last so that it never shadows a deliberate copy.
    """
    filename = NE_FILES.get(spacecraft, NE_FILES["isis1"])

    env = os.environ.get("TSS_ISIS_NE")
    if env:
        p = Path(env).expanduser()
        if p.is_dir():
            p = p / filename
        if p.is_file():
            return p

    d = cache.source_dir("spdf_isis", create=False)
    if d.is_dir():
        direct = d / filename
        if direct.is_file():
            return direct
        for candidate in sorted(d.glob(f"{Path(filename).stem}-*{Path(filename).suffix}")):
            return candidate

    for candidate in _AUTHOR_PATHS:
        p = Path(candidate) / filename
        if p.is_file():
            return p
    return None


_AUTHOR_PATHS: tuple[str, ...] = (
    r"C:\Users\Andrew\PhD\Dissertation\ISISValidation",
    "/mnt/user-data/uploads/PhD/Dissertation/ISISValidation",
)
"""Working copies on the author's machines.  Checked last, never written to."""


@register_source(
    name="spdf_isis_ne",
    description="ISIS-1 hand-scaled topside electron density profiles (CRC, via NASA SPDF)",
    url=NE_DIR_URL,
    licence="Public domain (US Government work); acknowledgement requested",
    citation=(
        "Communications Research Centre topside sounder Ne profiles, ISIS-1, scaled "
        "by Jackson's method; distributed by NASA Space Physics Data Facility, "
        "https://spdf.gsfc.nasa.gov/pub/data/isis/topside_sounder/crc_ne_profile_ascii/"
    ),
    cadence="static archive, closed mission",
    work_packages=("WP1", "WP8"),
    acknowledgement=(
        "Data provided by the NASA Space Physics Data Facility (SPDF).  Topside "
        "profiles were hand-scaled at the Communications Research Centre, Ottawa."
    ),
    notes=(
        "ISIS-1 hand-scaled profiles, 1969 to 1980, about 8.4 MB.  The scaled selection is not "
        "a uniform orbit-track sample: median |lat| is about 19 degrees against an "
        "80 degree inclination track, so distribution statistics confound model "
        "error with sampling and must be paired observation by observation.  "
        "Digitised ionograms for WP8 are a separate, much larger product and are "
        "not fetched here; see the status page."
    ),
    extra={"status_page": STATUS_URL, "ionogram_dirs": list(IONOGRAM_DIRS)},
)
def fetch_isis_ne(
    *,
    force: bool = False,
    limit: int | None = None,
    spacecraft: str = "isis1",
) -> list[Path]:
    """Locate or download the hand-scaled Ne profile file.

    A local copy wins unless ``force`` is set, because the archive is closed and
    a byte-identical file is a byte-identical file.  When a local copy is used
    and has no sidecar, one is written so that the run manifest can still record
    the digest.

    ``limit`` is accepted for interface consistency and ignored: the archive
    ships one file per spacecraft and there is nothing to truncate at fetch
    time.  Pass ``limit`` to :func:`read_ne_ascii` instead.
    """
    if spacecraft not in NE_FILES:
        raise FetchError(
            f"unknown spacecraft {spacecraft!r}; known: {', '.join(sorted(NE_FILES))}"
        )
    filename = NE_FILES[spacecraft]

    if not force:
        existing = local_ne_path(spacecraft)
        if existing is not None:
            if cache.read_meta(existing) is None and _is_writable(existing.parent):
                cache.record_meta(
                    existing,
                    url=NE_DIR_URL + filename,
                    source="spdf_isis_ne",
                    extra={"note": "located locally, not downloaded"},
                )
            return [existing]

    url = NE_DIR_URL + filename
    dest = cache.source_dir("spdf_isis") / filename
    path = cache.download(
        url,
        dest,
        source="spdf_isis_ne",
        force=force,
        max_age_days=ISIS_MAX_AGE_DAYS,
        min_interval_s=2.0,
    )
    return [path]


def _is_writable(directory: Path) -> bool:
    """True if a sidecar can be written next to a file we only read."""
    return os.access(directory, os.W_OK)
