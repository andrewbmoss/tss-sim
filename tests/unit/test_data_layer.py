"""Unit tests for the external data acquisition layer.

Every test here runs offline.  Nothing in this file opens a socket, and the two
tests that exercise HTTP behaviour do it against a fake session object rather
than against a server.  That is not only politeness towards the archives: a test
that depends on a remote host is a test that fails for reasons which have
nothing to do with the code, and a suite that fails for unrelated reasons stops
being read.

The fixture ``tests/fixtures/isis1_ne_sample.asc`` was written by hand to match
the CRC ``crc_ne_profile_ascii`` fixed-width Fortran layout, because the real
``isis1_ne.asc`` was not present on the machine where this module was written.
It is seven synthetic profiles with plausible topside numbers and exercises the
awkward parts of the format: header fields packed with no separators, negative
longitudes and latitudes butted up against their neighbours, and height and
density blocks that wrap across lines at 14 and 10 values respectively with a
short final line.  It is not real data and no scientific claim may be made from
it.  Replacing it with the first forty lines of the real file would be strictly
better and the tests would not need changing.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[2]
SRC = REPO / "src"
if str(SRC) not in sys.path:
    sys.path.insert(0, str(SRC))

FIXTURE = REPO / "tests" / "fixtures" / "isis1_ne_sample.asc"


# ---------------------------------------------------------------------------
# The registry
# ---------------------------------------------------------------------------


def test_sources_is_populated() -> None:
    from tss.data import SOURCES

    assert SOURCES, "importing tss.data must register at least one source"
    for expected in ("spdf_isis_ne", "wspr_live", "datagov_absorption"):
        assert expected in SOURCES, f"{expected} did not register"


def test_every_source_carries_licence_and_citation() -> None:
    """A source without rights and a reference cannot appear in the dissertation."""
    from tss.data import SOURCES

    for name, src in SOURCES.items():
        assert src.licence.strip(), f"{name} has no licence"
        assert src.citation.strip(), f"{name} has no citation"
        assert src.url.strip(), f"{name} has no url"
        assert src.description.strip(), f"{name} has no description"
        assert src.work_packages, f"{name} names no work package"
        assert callable(src.fetch), f"{name} has a non-callable fetch"


def test_source_refuses_to_be_built_without_a_licence() -> None:
    from tss.data.base import Source

    with pytest.raises(ValueError, match="licence"):
        Source(
            name="nameless",
            description="d",
            url="https://example.invalid/",
            licence="   ",
            citation="c",
            cadence="static",
            work_packages=("WP1",),
            fetch=lambda **_: [],
        )

    with pytest.raises(ValueError, match="citation"):
        Source(
            name="nameless",
            description="d",
            url="https://example.invalid/",
            licence="public domain",
            citation="",
            cadence="static",
            work_packages=("WP1",),
            fetch=lambda **_: [],
        )


def test_describe_sources_is_json_serialisable() -> None:
    from tss.data import describe_sources

    json.dumps(describe_sources())


def test_unknown_source_names_what_is_registered() -> None:
    from tss.data.base import UnknownSource, get_source

    with pytest.raises(UnknownSource, match="registered"):
        get_source("no_such_source")


def test_data_package_does_not_import_heavy_dependencies() -> None:
    """Importing tss.data must stay cheap enough for ``tss fetch list``.

    requests, pandas and h5py are imported inside the functions that need them.
    This checks the modules rather than sys.modules, because another test or
    another import may legitimately have pulled pandas in already.
    """
    import ast

    package = SRC / "tss" / "data"
    offenders: list[str] = []
    for path in sorted(package.glob("*.py")):
        tree = ast.parse(path.read_text(encoding="utf-8"))
        for node in tree.body:  # module scope only
            names: list[str] = []
            if isinstance(node, ast.Import):
                names = [a.name.split(".")[0] for a in node.names]
            elif isinstance(node, ast.ImportFrom) and node.level == 0:
                names = [(node.module or "").split(".")[0]]
            for n in names:
                if n in {"requests", "pandas", "h5py", "matplotlib", "scipy"}:
                    offenders.append(f"{path.name} imports {n} at module scope")
    assert not offenders, "; ".join(offenders)


# ---------------------------------------------------------------------------
# The cache
# ---------------------------------------------------------------------------


@pytest.fixture()
def tmp_cache(tmp_path, monkeypatch):
    """Point the cache at a temporary directory for the duration of one test."""
    monkeypatch.setenv("TSS_DATA_CACHE", str(tmp_path / "external"))
    from tss.data import cache

    return cache


def test_cache_root_follows_the_environment(tmp_cache, tmp_path) -> None:
    assert tmp_cache.cache_root() == tmp_path / "external"
    d = tmp_cache.source_dir("demo")
    assert d.is_dir() and d.name == "demo"


def test_cached_path_is_stable_and_distinguishes_queries(tmp_cache) -> None:
    a = tmp_cache.cached_path("demo", "https://example.invalid/data.csv?x=1", ".csv")
    b = tmp_cache.cached_path("demo", "https://example.invalid/data.csv?x=2", ".csv")
    again = tmp_cache.cached_path("demo", "https://example.invalid/data.csv?x=1", ".csv")

    assert a == again, "the same key must map to the same path"
    assert a != b, "different keys must not collide"
    assert a.suffix == ".csv"
    assert "data.csv" in a.name or "data" in a.name


def test_record_meta_and_read_meta_round_trip(tmp_cache) -> None:
    from tss.provenance import file_digest

    payload = b"absorption,frequency\n1.0,3.5e6\n"
    path = tmp_cache.source_dir("demo") / "sample.csv"
    path.write_bytes(payload)

    meta_file = tmp_cache.record_meta(
        path,
        url="https://example.invalid/sample.csv",
        etag='W/"abc123"',
        last_modified="Wed, 21 Oct 2015 07:28:00 GMT",
        source="demo",
        extra={"note": "written by a test"},
    )

    assert meta_file.name == "sample.csv.meta.json"
    assert meta_file.is_file()

    meta = tmp_cache.read_meta(path)
    assert meta is not None
    assert meta["url"] == "https://example.invalid/sample.csv"
    assert meta["etag"] == 'W/"abc123"'
    assert meta["last_modified"] == "Wed, 21 Oct 2015 07:28:00 GMT"
    assert meta["bytes"] == len(payload)
    assert meta["sha256"] == file_digest(path)
    assert meta["retrieved_utc"].endswith("+00:00")
    assert meta["extra"]["note"] == "written by a test"
    assert "Utah State University" in meta["user_agent"]

    assert tmp_cache.verify_cached(path) is True
    path.write_bytes(payload + b"tampered\n")
    assert tmp_cache.verify_cached(path) is False


def test_read_meta_returns_none_when_there_is_no_sidecar(tmp_cache) -> None:
    path = tmp_cache.source_dir("demo") / "orphan.csv"
    path.write_bytes(b"x\n")
    assert tmp_cache.read_meta(path) is None


def test_read_meta_survives_a_corrupt_sidecar(tmp_cache) -> None:
    path = tmp_cache.source_dir("demo") / "bad.csv"
    path.write_bytes(b"x\n")
    tmp_cache.meta_path(path).write_text("{not json", encoding="utf-8")
    assert tmp_cache.read_meta(path) is None


def test_is_fresh(tmp_cache) -> None:
    missing = tmp_cache.source_dir("demo") / "absent.csv"
    assert tmp_cache.is_fresh(missing, 1.0) is False

    empty = tmp_cache.source_dir("demo") / "empty.csv"
    empty.write_bytes(b"")
    assert tmp_cache.is_fresh(empty, 1.0) is False

    fresh = tmp_cache.source_dir("demo") / "fresh.csv"
    fresh.write_bytes(b"x\n")
    tmp_cache.record_meta(fresh, url="https://example.invalid/x", source="demo")
    assert tmp_cache.is_fresh(fresh, 1.0) is True
    assert tmp_cache.is_fresh(fresh, float("inf")) is True
    assert tmp_cache.is_fresh(fresh, 0.0) is False

    stale_meta = json.loads(tmp_cache.meta_path(fresh).read_text(encoding="utf-8"))
    stale_meta["retrieved_utc"] = "1999-01-01T00:00:00+00:00"
    tmp_cache.meta_path(fresh).write_text(json.dumps(stale_meta), encoding="utf-8")
    assert tmp_cache.is_fresh(fresh, 1.0) is False


def test_user_agent_names_the_project_and_the_institution(monkeypatch) -> None:
    from tss.data import cache

    monkeypatch.delenv("TSS_USER_AGENT", raising=False)
    ua = cache.user_agent()
    assert "tss-sim" in ua
    assert "Utah State University" in ua
    assert "@" not in ua, "no personal contact detail belongs in a request header"

    monkeypatch.setenv("TSS_USER_AGENT", "something-else/1.0")
    assert cache.user_agent() == "something-else/1.0"


# ---------------------------------------------------------------------------
# polite_get
# ---------------------------------------------------------------------------


class _Response:
    """The smallest thing polite_get treats as a response."""

    def __init__(self, status_code: int, body: bytes = b"", reason: str = "") -> None:
        self.status_code = status_code
        self.content = body
        self.reason = reason
        self.headers: dict[str, str] = {}

    @property
    def text(self) -> str:
        return self.content.decode("utf-8", "replace")


class _Session:
    """A fake requests session that records calls and never touches a socket."""

    def __init__(self, *responses: _Response) -> None:
        self._responses = list(responses)
        self.calls: list[dict] = []

    def get(self, url, params=None, headers=None, timeout=None, stream=False):
        self.calls.append(
            {"url": url, "params": params, "headers": dict(headers or {}),
             "timeout": timeout, "stream": stream}
        )
        if not self._responses:
            raise AssertionError("polite_get made more requests than the test allowed")
        return self._responses.pop(0)


@pytest.fixture()
def no_sleep(monkeypatch):
    """Make backoff and throttling instantaneous."""
    from tss.data import cache

    monkeypatch.setattr(cache.time, "sleep", lambda _s: None)
    return cache


@pytest.mark.parametrize("status", [400, 401, 403, 404, 410, 422, 429])
def test_polite_get_does_not_retry_a_4xx(no_sleep, status) -> None:
    """A client error is an answer.  Asking again will not change it."""
    from tss.data.base import FetchError

    session = _Session(_Response(status, b"no", reason="Client Error"))
    with pytest.raises(FetchError) as excinfo:
        no_sleep.polite_get(
            "https://example.invalid/thing",
            session=session,
            max_retries=5,
            min_interval_s=0.0,
        )

    assert len(session.calls) == 1, f"a {status} was retried {len(session.calls)} times"
    assert str(status) in str(excinfo.value)
    assert "not retrying" in str(excinfo.value)


def test_polite_get_retries_a_5xx_then_succeeds(no_sleep) -> None:
    session = _Session(
        _Response(503, b"busy"),
        _Response(502, b"busy"),
        _Response(200, b"payload"),
    )
    response = no_sleep.polite_get(
        "https://example.invalid/thing",
        session=session,
        max_retries=4,
        min_interval_s=0.0,
    )
    assert response.status_code == 200
    assert response.content == b"payload"
    assert len(session.calls) == 3


def test_polite_get_gives_up_after_max_retries(no_sleep) -> None:
    from tss.data.base import FetchError

    session = _Session(*[_Response(503, b"busy") for _ in range(4)])
    with pytest.raises(FetchError, match="giving up"):
        no_sleep.polite_get(
            "https://example.invalid/thing",
            session=session,
            max_retries=3,
            min_interval_s=0.0,
        )
    assert len(session.calls) == 4


def test_polite_get_sends_the_user_agent_and_never_disables_tls(no_sleep) -> None:
    session = _Session(_Response(200, b"ok"))
    no_sleep.polite_get(
        "https://example.invalid/thing", session=session, min_interval_s=0.0
    )
    call = session.calls[0]
    assert "Utah State University" in call["headers"]["User-Agent"]
    assert "verify" not in call, "polite_get must never pass a verify argument"


# ---------------------------------------------------------------------------
# WSPR
# ---------------------------------------------------------------------------


def test_quote_escapes_a_callsign() -> None:
    from tss.data.wspr import quote

    assert quote("W1AW") == "'W1AW'"
    assert quote("O'BRIEN") == "'O\\'BRIEN'"
    assert quote("back\\slash") == "'back\\\\slash'"
    assert quote(7) == "7"

    with pytest.raises(TypeError):
        quote(["W1AW"])


def test_spots_sql_quotes_callsigns_and_enforces_a_limit() -> None:
    from tss.data.wspr import DEFAULT_LIMIT, MAX_LIMIT, TABLE, build_spots_sql

    sql = build_spots_sql(
        "2024-03-01",
        "2024-03-02",
        band=7,
        tx_sign="w1aw",
        rx_sign="ve3abc",
        limit=250,
    )

    assert f"FROM {TABLE}" in sql
    assert "tx_sign = 'W1AW'" in sql
    assert "rx_sign = 'VE3ABC'" in sql
    assert "band = 7" in sql
    assert "time >= '2024-03-01 00:00:00'" in sql
    assert "time < '2024-03-02 00:00:00'" in sql
    assert sql.rstrip().endswith("LIMIT 250")

    # A limit is always emitted, even when the caller asks for none.
    default = build_spots_sql("2024-03-01", "2024-03-02")
    assert f"LIMIT {DEFAULT_LIMIT}" in default

    # And it is clamped rather than honoured without bound.
    clamped = build_spots_sql("2024-03-01", "2024-03-02", limit=MAX_LIMIT * 10)
    assert f"LIMIT {MAX_LIMIT}" in clamped


def test_spots_sql_refuses_an_injected_callsign() -> None:
    """A quote in a callsign is escaped, not passed through."""
    from tss.data.wspr import build_spots_sql, quote

    nasty = "W1AW' OR 1=1; DROP TABLE wspr.rx --"
    sql = build_spots_sql("2024-03-01", "2024-03-02", tx_sign=nasty)

    # The text survives, but only as the contents of one escaped string literal,
    # so the semicolon and the comment marker are data and not syntax.
    assert quote(nasty.upper()) in sql
    assert "\\'" in sql, "the injected quote must have been escaped"
    assert "OR 1=1" not in sql.replace(quote(nasty.upper()), ""), (
        "nothing from the callsign may appear outside the literal"
    )
    assert sql.rstrip().split()[-2] == "LIMIT"


def test_spots_sql_rejects_a_backwards_window_and_unknown_columns() -> None:
    from tss.data.wspr import build_spots_sql

    with pytest.raises(ValueError, match="not after"):
        build_spots_sql("2024-03-02", "2024-03-01")

    with pytest.raises(ValueError, match="unknown column"):
        build_spots_sql("2024-03-01", "2024-03-02", columns=("time", "secret"))

    with pytest.raises(ValueError, match="positive"):
        build_spots_sql("2024-03-01", "2024-03-02", limit=0)


def test_spots_sql_takes_a_band_list() -> None:
    from tss.data.wspr import build_spots_sql

    sql = build_spots_sql("2024-03-01", "2024-03-02", band=[3, 7, 14])
    assert "band IN (3, 7, 14)" in sql


def test_wspr_docstring_states_the_limitation() -> None:
    """The ground-to-ground caveat is load bearing and must not be edited out."""
    import tss.data.wspr as wspr

    doc = (wspr.__doc__ or "").lower()
    assert "ground-to-ground" in doc
    assert "not a topside validation" in doc or "nothing measured here validates" in doc
    assert "archive" in doc


# ---------------------------------------------------------------------------
# ISIS reader
# ---------------------------------------------------------------------------


def test_fixture_exists() -> None:
    assert FIXTURE.is_file(), f"missing fixture {FIXTURE}"


def test_read_ne_ascii_parses_the_fixture() -> None:
    from tss.data.spdf_isis import read_ne_ascii

    frame = read_ne_ascii(FIXTURE)

    assert len(frame) == 7
    assert frame.attrs["declared_records"] == 7

    for col in (
        "itype", "iquality", "date", "time", "lon", "lat", "dip", "L", "xhi",
        "rz12", "ig12", "rtec", "ln_nmf2_iri", "hmf2_iri", "tec_iri", "n",
        "h_km", "ne_cm3", "h_sat_km", "ne_sat_cm3", "fn_sat_mhz",
        "nmf2_iri_cm3", "fof2_iri_mhz", "ne_bot_cm3", "h_bot_km", "fo_bot_mhz",
        "year", "doy", "ut_hour", "lt_hour", "epoch",
    ):
        assert col in frame.columns, f"missing column {col}"

    first = frame.iloc[0]
    assert first["itype"] == 1
    assert first["date"] == 69120
    assert first["time"] == 123045
    assert first["n"] == 18
    assert first["lon"] == pytest.approx(250.5)
    assert first["lat"] == pytest.approx(45.20)
    assert first["year"] == 1969
    assert first["doy"] == 120
    assert first["ut_hour"] == pytest.approx(12.5)

    # The two-digit year pivots at 50, so 69 is 1969 and not 2069.
    assert set(frame["year"]) <= {1969, 1970, 1971, 1972}

    # A negative longitude butted against a negative latitude still splits.
    second = frame.iloc[1]
    assert second["lon"] == pytest.approx(-112.3)
    assert second["lat"] == pytest.approx(-22.10)
    assert second["dip"] == pytest.approx(-41.80)


def test_read_ne_ascii_profile_semantics() -> None:
    """The first sample of a profile is the spacecraft, which is what INS-01 needs."""
    from tss.data.spdf_isis import PLASMA, read_ne_ascii

    frame = read_ne_ascii(FIXTURE)
    for _, row in frame.iterrows():
        h = np.asarray(row["h_km"], dtype=float)
        ne = np.asarray(row["ne_cm3"], dtype=float)

        assert len(h) == row["n"]
        assert len(ne) == row["n"]

        # Profiles run downward from the satellite.
        assert h[0] == pytest.approx(row["h_sat_km"])
        assert ne[0] == pytest.approx(row["ne_sat_cm3"])
        assert h[-1] == pytest.approx(row["h_bot_km"])
        assert h[0] > h[-1]

        # Topside density increases downward towards the F2 peak.
        assert ne[-1] > ne[0]

        expected = PLASMA * np.sqrt(row["ne_sat_cm3"] * 1e6) / 1e6
        assert row["fn_sat_mhz"] == pytest.approx(expected)


def test_read_ne_ascii_can_drop_the_profiles_and_can_limit() -> None:
    from tss.data.spdf_isis import read_ne_ascii

    scalar = read_ne_ascii(FIXTURE, want_profiles=False)
    assert "h_km" not in scalar.columns
    assert "ne_cm3" not in scalar.columns
    assert len(scalar) == 7

    short = read_ne_ascii(FIXTURE, limit=3)
    assert len(short) == 3


def test_read_ne_ascii_reports_a_missing_file() -> None:
    from tss.data.spdf_isis import read_ne_ascii

    with pytest.raises(FileNotFoundError, match="TSS_ISIS_NE"):
        read_ne_ascii(FIXTURE.with_name("definitely-not-here.asc"))


def test_iter_profiles_yields_usable_records() -> None:
    from tss.data.spdf_isis import Profile, iter_profiles

    profiles = list(iter_profiles(FIXTURE))
    assert len(profiles) == 7
    assert all(isinstance(p, Profile) for p in profiles)

    p = profiles[0]
    assert p.epoch is not None
    assert p.epoch.year == 1969
    assert p.epoch.hour == 12 and p.epoch.minute == 30 and p.epoch.second == 45
    assert p.lat_deg == pytest.approx(45.20)
    assert p.lon_deg == pytest.approx(250.5)
    assert p.sat_alt_km == pytest.approx(880.0)

    alt_m, ne_m3 = p.ascending()
    assert np.all(np.diff(alt_m) > 0), "ascending() must sort by altitude"
    assert alt_m[-1] == pytest.approx(p.sat_alt_km * 1e3)
    assert ne_m3[-1] == pytest.approx(p.ne_cm3[0] * 1e6)
    assert p.plasma_freq_sat_mhz > 0.0


def test_iter_profiles_filters_on_quality() -> None:
    from tss.data.spdf_isis import iter_profiles

    good = list(iter_profiles(FIXTURE, max_quality=1))
    assert good, "the fixture should contain at least one well-scaled profile"
    assert all(p.quality is not None and p.quality <= 1 for p in good)
    assert len(good) < 7


# ---------------------------------------------------------------------------
# data.gov absorption
# ---------------------------------------------------------------------------


def test_absorption_reader_names_the_columns_it_found(tmp_path) -> None:
    """An unmappable schema must fail loudly and say what it saw."""
    from tss.data.datagov_absorption import AbsorptionSchemaError, read_absorption

    path = tmp_path / "mystery.csv"
    path.write_text("alpha,beta,gamma\n1,2,3\n4,5,6\n", encoding="utf-8")

    with pytest.raises(AbsorptionSchemaError) as excinfo:
        read_absorption(path)

    message = str(excinfo.value)
    assert "alpha" in message and "beta" in message and "gamma" in message
    assert "COLUMN_ALIASES" in message


def test_absorption_reader_normalises_expected_names(tmp_path) -> None:
    from tss.data.datagov_absorption import read_absorption

    path = tmp_path / "abs.csv"
    path.write_text(
        "Freq_MHz,Absorption (dB),Station\n"
        "2.5,14.2,BOULDER\n"
        "5.0,4.1,BOULDER\n"
        "10.0,1.0,BOULDER\n",
        encoding="utf-8",
    )
    tidy = read_absorption(path)

    assert list(tidy.columns) == ["frequency_hz", "absorption_db", "station"]
    assert tidy["frequency_hz"].iloc[0] == pytest.approx(2.5e6)
    assert tidy["absorption_db"].iloc[0] == pytest.approx(14.2)
    assert tidy.attrs["frequency_unit"] == "mhz"
    assert tidy.attrs["frequency_unit_assumed"] is False


def test_reference_absorption_interpolates_and_flags_extrapolation(tmp_path) -> None:
    from tss.data.datagov_absorption import read_absorption, reference_absorption

    path = tmp_path / "abs.csv"
    path.write_text(
        "freq_mhz,absorption_db\n2.0,16.0\n4.0,4.0\n8.0,1.0\n", encoding="utf-8"
    )
    frame = read_absorption(path)

    inside = reference_absorption([2.0e6, 8.0e6], frame=frame)
    assert inside[0] == pytest.approx(16.0)
    assert inside[1] == pytest.approx(1.0)

    outside = reference_absorption([1.0e6], frame=frame, extrapolate="nan")
    assert np.isnan(outside[0])

    clipped = reference_absorption([1.0e6], frame=frame, extrapolate="clip")
    assert clipped[0] == pytest.approx(16.0)

    scaled = reference_absorption([1.0e6], frame=frame, extrapolate="f_squared")
    assert scaled[0] == pytest.approx(16.0 * 4.0)

    with pytest.raises(ValueError, match="unknown extrapolate"):
        reference_absorption([1.0e6], frame=frame, extrapolate="magic")


# ---------------------------------------------------------------------------
# Local model discovery
# ---------------------------------------------------------------------------


def test_probe_reports_every_model_without_raising() -> None:
    """probe() must answer on a machine with none of the models installed.

    ``trial=False`` is used here so the suite never spends a minute starting
    MATLAB.  The trial path is the default in production and is argued for in
    the module docstring; what is tested here is that the inventory is total and
    JSON-serialisable.
    """
    from tss.data.local_models import probe

    status = probe(trial=False)
    for name in ("pharlap", "iri2016", "nequick2", "matlab"):
        assert name in status
        assert isinstance(status[name]["available"], bool)
        if not status[name]["available"]:
            assert status[name]["reason"], f"{name} is unavailable without a reason"
    json.dumps(status)


def test_require_raises_something_actionable(monkeypatch) -> None:
    from tss.data import local_models
    from tss.data.local_models import MissingModel, ModelStatus

    monkeypatch.setattr(local_models, "_probe_cache", {})
    monkeypatch.setitem(
        local_models._PROBES,
        "pharlap",
        lambda trial: ModelStatus(
            name="pharlap",
            available=False,
            reason="set TSS_PHARLAP_PATH to the unpacked PHaRLAP root",
        ),
    )
    with pytest.raises(MissingModel, match="TSS_PHARLAP_PATH"):
        local_models.require("pharlap", trial=False)


def test_require_rejects_an_unknown_name() -> None:
    from tss.data.local_models import require

    with pytest.raises(KeyError, match="known"):
        require("not_a_model", trial=False)


def test_local_models_docstring_records_the_exist_lesson() -> None:
    """The reasoning for trial calls is the point of the module, not decoration."""
    import tss.data.local_models as lm

    doc = lm.__doc__ or ""
    assert "exist(" in doc
    assert "trial" in doc.lower()
