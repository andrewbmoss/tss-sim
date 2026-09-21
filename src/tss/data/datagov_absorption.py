"""Layer: data
WP: WP2
Claims: -
Status: partial
Requirements: -

The data.gov data set "Measurements of the absorptive properties of the
ionosphere".

WP2 computes integrated D-region absorption along the sounding path and has to
be scored against something measured.  This data set is the reference chosen
for that: a NOAA-lineage collection of ionospheric absorption measurements
catalogued at data.gov.  The catalogue entry is at

    https://catalog.data.gov/dataset/measurements-of-the-absorptive-properties-of-the-ionosphere

Discovery rather than hardcoding.  data.gov is a CKAN catalogue, and a CKAN
catalogue moves its files.  Hardcoding a resource URL produces a fetcher that
works today and 404s in a year with no explanation.  This module therefore asks
``/api/3/action/package_show`` for the package, reads the resource list, and
downloads what it finds.  The package metadata is cached with its own sidecar,
so a run manifest records which revision of the catalogue entry was used as well
as which files came out of it.

Status is "partial" and honestly so.  The exact column layout of the delivered
resources has not been verified against a download, because this module was
written with no network access.  Rather than invent a schema, the reader
introspects whatever columns arrive, maps them through
:data:`COLUMN_ALIASES`, and raises :class:`AbsorptionSchemaError` naming the
columns it actually found when it cannot make the mapping.  The first person to
run the fetcher will get either a tidy frame or an error that tells them exactly
what to add to the alias table.  That is a deliberate trade: a loud, specific
failure is worth more than a quiet, wrong parse.

:func:`reference_absorption` is the accessor WP2 scores against.  Read its
docstring before using it; the interpolation it does across frequency carries an
assumption that has to be checked once the real table is in hand.
"""

from __future__ import annotations

import json
from collections.abc import Mapping, Sequence
from dataclasses import dataclass
from pathlib import Path
from typing import Any

import numpy as np

from tss.data import cache
from tss.data.base import FetchError, register_source

__all__ = [
    "CATALOG_URL",
    "CKAN_API",
    "DATASET_ID",
    "COLUMN_ALIASES",
    "REQUIRED_FIELDS",
    "AbsorptionSchemaError",
    "AbsorptionTable",
    "package_metadata",
    "resource_urls",
    "fetch_absorption",
    "read_absorption",
    "reference_absorption",
]


DATASET_ID = "measurements-of-the-absorptive-properties-of-the-ionosphere"
CATALOG_URL = f"https://catalog.data.gov/dataset/{DATASET_ID}"
CKAN_API = "https://catalog.data.gov/api/3/action/package_show"

METADATA_MAX_AGE_DAYS = 30.0
DATA_MAX_AGE_DAYS = 365.0

_DOWNLOADABLE_SUFFIXES = (".csv", ".txt", ".dat", ".tsv", ".asc", ".zip", ".json")


class AbsorptionSchemaError(ValueError):
    """Raised when the delivered file's columns cannot be mapped.

    The message always names the columns that were actually present, because the
    fix is to add them to :data:`COLUMN_ALIASES` and that cannot be done without
    knowing what they are.
    """


# ---------------------------------------------------------------------------
# Schema normalisation
# ---------------------------------------------------------------------------

COLUMN_ALIASES: dict[str, tuple[str, ...]] = {
    "frequency_hz": (
        "frequency_hz", "freq_hz", "frequency", "freq", "f_hz",
        "frequency_mhz", "freq_mhz", "f_mhz", "frequency_khz", "freq_khz",
        "sounding_frequency", "operating_frequency",
    ),
    "absorption_db": (
        "absorption_db", "absorption", "abs_db", "l_db", "loss_db",
        "absorption_decibels", "non_deviative_absorption_db", "a_db",
    ),
    "time_utc": (
        "time_utc", "time", "datetime", "date_time", "timestamp", "epoch",
        "observation_time", "ut",
    ),
    "date": ("date", "obs_date", "yyyymmdd"),
    "lat_deg": ("lat_deg", "lat", "latitude", "geographic_latitude"),
    "lon_deg": ("lon_deg", "lon", "long", "longitude", "geographic_longitude"),
    "station": ("station", "station_name", "site", "observatory", "station_id"),
    "solar_zenith_deg": (
        "solar_zenith_deg", "solar_zenith_angle", "chi", "sza", "zenith_angle",
    ),
    "quality": ("quality", "qc", "flag", "quality_flag"),
}
"""Expected tidy name to the upstream spellings that map onto it.

Matching is case-insensitive and ignores spaces, hyphens and underscores.  A
frequency column whose name carries a unit (``freq_mhz``, ``freq_khz``) is
converted to hertz on read; a bare ``frequency`` is assumed to be in megahertz
and that assumption is recorded in the frame's ``attrs`` so it can be checked.
"""

REQUIRED_FIELDS: tuple[str, ...] = ("frequency_hz", "absorption_db")
"""Without a frequency and an absorption there is nothing for WP2 to score."""

_FREQ_SCALE: dict[str, float] = {
    "hz": 1.0,
    "khz": 1.0e3,
    "mhz": 1.0e6,
}


def _canon(name: str) -> str:
    """Normalise a column name for matching."""
    return "".join(ch for ch in str(name).lower() if ch.isalnum())


_ALIAS_LOOKUP: dict[str, str] = {
    _canon(alias): tidy
    for tidy, aliases in COLUMN_ALIASES.items()
    for alias in aliases
}


def _map_columns(columns: Sequence[str]) -> dict[str, str]:
    """Map delivered column names onto tidy names.  Unmapped columns are dropped.

    The first delivered column that maps onto a tidy name wins.  A file with two
    plausible absorption columns is a file whose schema needs reading by a human,
    and picking the second silently would hide that.
    """
    mapping: dict[str, str] = {}
    taken: set[str] = set()
    for col in columns:
        tidy = _ALIAS_LOOKUP.get(_canon(col))
        if tidy is None or tidy in taken:
            continue
        mapping[col] = tidy
        taken.add(tidy)
    return mapping


def _frequency_unit(original_name: str) -> str:
    """Guess the unit of a frequency column from its name.

    Returns one of ``hz``, ``khz``, ``mhz``.  A name with no unit in it returns
    ``mhz``, which is the convention ionospheric absorption tables overwhelmingly
    use, and the caller records that it was a guess.
    """
    c = _canon(original_name)
    if c.endswith("khz") or "khz" in c:
        return "khz"
    if c.endswith("mhz") or "mhz" in c:
        return "mhz"
    if c.endswith("hz") or "hz" in c:
        return "hz"
    return "mhz"


# ---------------------------------------------------------------------------
# CKAN discovery
# ---------------------------------------------------------------------------


def package_metadata(*, force: bool = False, dataset_id: str = DATASET_ID) -> dict[str, Any]:
    """Fetch and cache the CKAN ``package_show`` record for the data set.

    Returns the ``result`` object.  Raises :class:`~tss.data.base.FetchError` if
    CKAN reports failure, which it does with ``success: false`` and an HTTP 200,
    so the status code alone is not enough to tell.
    """
    path = cache.cached_path("datagov_absorption", f"{CKAN_API}?id={dataset_id}", ".json")
    if not force and cache.is_fresh(path, METADATA_MAX_AGE_DAYS):
        payload = json.loads(path.read_text(encoding="utf-8"))
    else:
        response = cache.polite_get(
            CKAN_API, params={"id": dataset_id}, timeout=60.0, min_interval_s=2.0
        )
        payload = response.json()
        path.write_text(json.dumps(payload, indent=2, sort_keys=True), encoding="utf-8")
        cache.record_meta(
            path,
            url=CKAN_API,
            source="datagov_absorption",
            extra={"dataset_id": dataset_id},
        )

    if not isinstance(payload, dict) or not payload.get("success"):
        error = payload.get("error") if isinstance(payload, dict) else None
        raise FetchError(
            f"CKAN package_show for {dataset_id!r} did not succeed: {error!r}.  "
            f"Check the catalogue entry at {CATALOG_URL}"
        )
    result = payload.get("result")
    if not isinstance(result, dict):
        raise FetchError(f"CKAN returned no result object for {dataset_id!r}")
    return result


def resource_urls(
    metadata: Mapping[str, Any] | None = None,
    *,
    suffixes: Sequence[str] = _DOWNLOADABLE_SUFFIXES,
    force: bool = False,
) -> list[dict[str, str]]:
    """The downloadable resources in the package, as ``{name, url, format}``.

    HTML landing pages and API endpoints are filtered out by extension and by
    declared format, because a CKAN resource list routinely mixes them in with
    the data files and downloading a landing page produces a confusing file
    rather than an error.
    """
    meta = metadata if metadata is not None else package_metadata(force=force)
    out: list[dict[str, str]] = []
    for res in meta.get("resources", []) or []:
        url = str(res.get("url") or "").strip()
        if not url:
            continue
        fmt = str(res.get("format") or "").strip().lower()
        tail = url.split("?", 1)[0].lower()
        if fmt in {"html", "web page", "api"}:
            continue
        if not tail.endswith(tuple(suffixes)) and fmt not in {
            "csv", "text", "txt", "zip", "json", "tsv", "ascii",
        }:
            continue
        out.append(
            {
                "name": str(res.get("name") or Path(tail).name or "resource"),
                "url": url,
                "format": fmt or Path(tail).suffix.lstrip("."),
                "description": str(res.get("description") or ""),
            }
        )
    return out


@register_source(
    name="datagov_absorption",
    description="Measurements of the absorptive properties of the ionosphere (data.gov)",
    url=CATALOG_URL,
    licence="US Government work, public domain (see the catalogue entry for the exact rights statement)",
    citation=(
        "Measurements of the absorptive properties of the ionosphere, data.gov "
        f"catalogue entry {DATASET_ID}, {CATALOG_URL}"
    ),
    cadence="static archive, catalogue metadata may be revised",
    work_packages=("WP2",),
    acknowledgement=(
        "Data obtained through data.gov; original collection by NOAA National "
        "Centers for Environmental Information and its predecessors."
    ),
    notes=(
        "Resource URLs are discovered through the CKAN API rather than hardcoded.  "
        "The column schema of the delivered files has not been verified against a "
        "download, so the reader introspects the delivered columns and raises a "
        "named error rather than assuming a layout."
    ),
    extra={"ckan_api": CKAN_API, "dataset_id": DATASET_ID},
)
def fetch_absorption(*, force: bool = False, limit: int | None = None) -> list[Path]:
    """Discover and download the data set's resources.

    ``limit`` caps the number of resources downloaded, which is useful when the
    package turns out to hold many files and only the first is wanted for a look.
    Returns the local paths in catalogue order.
    """
    meta = package_metadata(force=force)
    resources = resource_urls(meta)
    if not resources:
        raise FetchError(
            f"the CKAN record for {DATASET_ID!r} lists no downloadable resources.  "
            f"Open {CATALOG_URL} and check whether the files moved."
        )
    if limit is not None:
        resources = resources[: max(0, int(limit))]

    paths: list[Path] = []
    for res in resources:
        dest = cache.cached_path("datagov_absorption", res["url"], _suffix_for(res))
        paths.append(
            cache.download(
                res["url"],
                dest,
                source="datagov_absorption",
                force=force,
                max_age_days=DATA_MAX_AGE_DAYS,
                min_interval_s=2.0,
                extra_meta={
                    "resource_name": res["name"],
                    "resource_format": res["format"],
                    "dataset_id": DATASET_ID,
                },
            )
        )
    return paths


def _suffix_for(resource: Mapping[str, str]) -> str:
    """A file suffix for a resource, from its URL then its declared format."""
    tail = resource["url"].split("?", 1)[0]
    suffix = Path(tail).suffix
    if suffix and len(suffix) <= 6:
        return suffix.lower()
    fmt = (resource.get("format") or "").strip().lower()
    return f".{fmt}" if fmt and fmt.isalnum() else ".dat"


# ---------------------------------------------------------------------------
# Reading
# ---------------------------------------------------------------------------


def read_absorption(path: str | Path, **read_kwargs: Any) -> Any:
    """Read one downloaded resource into a tidy DataFrame.

    Tidy means one row per measurement with at least ``frequency_hz`` and
    ``absorption_db``, plus whatever of ``time_utc``, ``lat_deg``, ``lon_deg``,
    ``station``, ``solar_zenith_deg`` and ``quality`` the file carries.  Columns
    that do not map are dropped, and the dropped names are recorded in
    ``attrs['unmapped_columns']`` so nothing disappears silently.

    The delimiter is sniffed for ``.txt`` and ``.dat`` files and taken as a comma
    for ``.csv``.  If the file cannot be parsed as a table at all, or if the
    required fields are missing after mapping, this raises
    :class:`AbsorptionSchemaError` naming every column it found.

    ``attrs['frequency_unit_assumed']`` is True when the frequency column's name
    carried no unit and megahertz was assumed.  Check that before trusting a
    scored result.
    """
    import pandas as pd

    p = Path(path)
    if not p.is_file():
        raise FileNotFoundError(f"no absorption resource at {p}; run `tss fetch datagov_absorption`")

    suffix = p.suffix.lower()
    try:
        if suffix in {".csv", ".tsv"}:
            sep = "\t" if suffix == ".tsv" else ","
            raw = pd.read_csv(p, sep=sep, **read_kwargs)
        elif suffix == ".json":
            raw = pd.read_json(p)
        else:
            raw = pd.read_csv(p, sep=None, engine="python", **read_kwargs)
    except Exception as exc:
        raise AbsorptionSchemaError(
            f"could not read {p} as a table: {exc}.  The delivered format may not be "
            "delimited text; open the file and add a branch to read_absorption()."
        ) from exc

    found = [str(c) for c in raw.columns]
    mapping = _map_columns(found)
    missing = [f for f in REQUIRED_FIELDS if f not in mapping.values()]
    if missing:
        raise AbsorptionSchemaError(
            f"{p.name}: cannot find {', '.join(missing)} among the delivered columns.  "
            f"Columns actually present: {', '.join(found) or '(none)'}.  "
            f"Add the right spelling to COLUMN_ALIASES in {__name__}."
        )

    tidy = raw.rename(columns=mapping)[list(mapping.values())].copy()

    freq_source = next(orig for orig, t in mapping.items() if t == "frequency_hz")
    unit = _frequency_unit(freq_source)
    tidy["frequency_hz"] = pd.to_numeric(tidy["frequency_hz"], errors="coerce") * _FREQ_SCALE[unit]
    tidy["absorption_db"] = pd.to_numeric(tidy["absorption_db"], errors="coerce")

    for col in ("lat_deg", "lon_deg", "solar_zenith_deg"):
        if col in tidy.columns:
            tidy[col] = pd.to_numeric(tidy[col], errors="coerce")
    if "time_utc" in tidy.columns:
        tidy["time_utc"] = pd.to_datetime(tidy["time_utc"], utc=True, errors="coerce")

    before = len(tidy)
    tidy = tidy.dropna(subset=["frequency_hz", "absorption_db"])

    tidy.attrs["source_file"] = str(p)
    tidy.attrs["dataset_id"] = DATASET_ID
    tidy.attrs["catalog_url"] = CATALOG_URL
    tidy.attrs["column_mapping"] = mapping
    tidy.attrs["unmapped_columns"] = [c for c in found if c not in mapping]
    tidy.attrs["frequency_unit"] = unit
    tidy.attrs["frequency_unit_assumed"] = _frequency_unit_was_guessed(freq_source)
    tidy.attrs["rows_dropped_nonnumeric"] = before - len(tidy)
    meta = cache.read_meta(p)
    if meta:
        tidy.attrs["sha256"] = meta.get("sha256")
        tidy.attrs["url"] = meta.get("url")
    return tidy


def _frequency_unit_was_guessed(name: str) -> bool:
    """True if the frequency column's name carried no unit at all."""
    c = _canon(name)
    return not (c.endswith("hz") or "hz" in c)


# ---------------------------------------------------------------------------
# The WP2 accessor
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class AbsorptionTable:
    """Frequency-sorted reference absorption, ready to interpolate.

    ``freq_hz`` is strictly increasing and ``absorption_db`` is the aggregate
    (median by default) of the measurements in each frequency bin.  ``count``
    says how many measurements each point rests on, because a reference point
    built from three observations should not be scored against as if it were
    built from three thousand.
    """

    freq_hz: np.ndarray
    absorption_db: np.ndarray
    count: np.ndarray
    statistic: str = "median"
    source_file: str = ""
    selection: Mapping[str, Any] | None = None


def _build_table(
    frame: Any,
    *,
    statistic: str = "median",
    bins: int | Sequence[float] | None = None,
    selection: Mapping[str, Any] | None = None,
) -> AbsorptionTable:
    """Aggregate a tidy frame into an :class:`AbsorptionTable`."""
    import pandas as pd

    f = np.asarray(frame["frequency_hz"], dtype=float)
    a = np.asarray(frame["absorption_db"], dtype=float)

    if bins is None:
        key = f
    else:
        edges = np.asarray(bins, dtype=float) if not isinstance(bins, int) else np.linspace(
            np.nanmin(f), np.nanmax(f), int(bins) + 1
        )
        idx = np.clip(np.digitize(f, edges) - 1, 0, len(edges) - 2)
        key = 0.5 * (edges[idx] + edges[idx + 1])

    grouped = pd.DataFrame({"f": key, "a": a}).groupby("f")["a"]
    agg = getattr(grouped, statistic)()
    counts = grouped.count()
    order = np.argsort(np.asarray(agg.index, dtype=float))
    return AbsorptionTable(
        freq_hz=np.asarray(agg.index, dtype=float)[order],
        absorption_db=np.asarray(agg.values, dtype=float)[order],
        count=np.asarray(counts.values, dtype=float)[order],
        statistic=statistic,
        source_file=str(frame.attrs.get("source_file", "")),
        selection=dict(selection or {}),
    )


def reference_absorption(
    freq_hz: float | Sequence[float] | np.ndarray,
    *,
    table: AbsorptionTable | None = None,
    frame: Any | None = None,
    path: str | Path | None = None,
    statistic: str = "median",
    bins: int | Sequence[float] | None = None,
    extrapolate: str = "f_squared",
) -> np.ndarray:
    """Reference absorption in dB at the given frequencies, for WP2 to score against.

    Supply exactly one of ``table``, ``frame`` or ``path``; with none of them the
    function reads the first cached resource it finds.

    Inside the range of the reference table the value is interpolated linearly in
    ``log(frequency)`` against ``absorption_db``, because non-deviative
    absorption falls roughly as an inverse power of frequency and a linear
    interpolation in frequency across a decade would be visibly wrong.

    Outside the range the behaviour depends on ``extrapolate``:

    ``"f_squared"`` (the default)
        Scale the nearest endpoint by ``(f_end / f)**2``, the classical
        non-deviative dependence well above the gyrofrequency.  This is a model
        assumption and not a measurement, and it is the weakest part of this
        function.  WP2 should prefer to report only over the measured range and
        treat anything outside it as unscored.
    ``"nan"``
        Return NaN outside the measured range, which is the honest answer.
    ``"clip"``
        Hold the endpoint value.

    Returns an array the shape of ``freq_hz``, or a scalar-shaped array for a
    scalar input.  The caller should check
    ``AbsorptionTable.count`` before drawing a conclusion from a sparsely
    sampled frequency.
    """
    if table is None:
        if frame is None:
            if path is None:
                path = _first_cached_resource()
            frame = read_absorption(path)
        table = _build_table(frame, statistic=statistic, bins=bins)

    f = np.atleast_1d(np.asarray(freq_hz, dtype=float))
    ft = table.freq_hz
    at = table.absorption_db
    if ft.size == 0:
        raise AbsorptionSchemaError("the reference absorption table is empty")
    if ft.size == 1:
        out = np.full(f.shape, at[0], dtype=float)
        return out if np.ndim(freq_hz) else out

    out = np.interp(np.log(f), np.log(ft), at, left=np.nan, right=np.nan)

    below = f < ft[0]
    above = f > ft[-1]
    if extrapolate == "nan":
        pass
    elif extrapolate == "clip":
        out[below] = at[0]
        out[above] = at[-1]
    elif extrapolate == "f_squared":
        with np.errstate(divide="ignore", invalid="ignore"):
            out[below] = at[0] * (ft[0] / f[below]) ** 2
            out[above] = at[-1] * (ft[-1] / f[above]) ** 2
    else:
        raise ValueError(
            f"unknown extrapolate mode {extrapolate!r}; use 'f_squared', 'nan' or 'clip'"
        )
    return out


def _first_cached_resource() -> Path:
    """The first already-downloaded resource, or a clear error saying to fetch."""
    d = cache.source_dir("datagov_absorption", create=False)
    if d.is_dir():
        for candidate in sorted(d.iterdir()):
            if candidate.is_file() and candidate.suffix.lower() in {".csv", ".txt", ".dat", ".tsv"}:
                return candidate
    raise FileNotFoundError(
        "no cached absorption resource; run `tss fetch datagov_absorption` first.  "
        f"Catalogue entry: {CATALOG_URL}"
    )
