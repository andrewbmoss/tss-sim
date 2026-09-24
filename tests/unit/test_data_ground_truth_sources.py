"""Tests for the Madrigal, GIRO DIDBase and IRTAM sources.

Every fixture under ``tests/fixtures/{madrigal,giro,irtam}/`` is a real response
captured from the live service on 2026-09-23, trimmed where noted and with the
requester's IP address redacted from GIRO headers.  No test here touches the
network: fetch paths are driven through a fake transport patched over
``tss.data.cache.polite_get``.
"""

from __future__ import annotations

import json
import re
from pathlib import Path

import numpy as np
import pytest

REPO = Path(__file__).resolve().parents[2]
FIX = REPO / "tests" / "fixtures"


class _Resp:
    def __init__(self, body: bytes, status: int = 200) -> None:
        self.content = body
        self.status_code = status
        self.headers: dict[str, str] = {}
        self.reason = ""

    @property
    def text(self) -> str:
        return self.content.decode("utf-8", "replace")


@pytest.fixture()
def tmp_cache(tmp_path, monkeypatch):
    monkeypatch.setenv("TSS_DATA_CACHE", str(tmp_path / "external"))
    return tmp_path / "external"


# ---------------------------------------------------------------------------
# Fixture hygiene
# ---------------------------------------------------------------------------

_IPV4 = re.compile(r"\b(?:\d{1,3}\.){3}\d{1,3}\b")


@pytest.mark.parametrize("sub", ["giro", "irtam", "madrigal"])
def test_fixtures_carry_no_ip_address_or_personal_identity(sub) -> None:
    files = sorted((FIX / sub).iterdir())
    assert files, f"no fixtures under {sub}"
    for f in files:
        text = f.read_text(encoding="utf-8", errors="replace")
        assert not _IPV4.search(text), f"{f.name} contains an IP address"
        assert "sdl.usu.edu" not in text and "user_email" not in text, f.name


def test_every_new_source_is_registered_with_licence_and_citation() -> None:
    from tss.data import SOURCES

    for name in ("madrigal", "giro_didbase", "irtam"):
        src = SOURCES[name]
        assert src.licence.strip() and src.citation.strip()
        assert "WP1" in src.work_packages


# ---------------------------------------------------------------------------
# GIRO DIDBase
# ---------------------------------------------------------------------------


def test_didbase_parses_the_millstone_fixture() -> None:
    from tss.data.giro import parse_scaled

    r = parse_scaled((FIX / "giro" / "didbase_MHJ45_20240301.txt").read_bytes())
    assert r.ursi_code == "MHJ45" and r.station == "MILLSTONE HILL"
    assert r.lat_deg == pytest.approx(42.6) and r.lon_deg == pytest.approx(288.5)
    assert "DPS-4D" in r.instrument
    assert len(r) == 17
    assert list(r.values) == ["foF2", "hmF2", "B0", "B1", "fmin"]
    assert r.units == {"foF2": "MHz", "hmF2": "km", "B0": "km", "B1": "", "fmin": "MHz"}
    assert r.times[0] == np.datetime64("2024-03-01T00:00:00", "ms")
    assert r.times[-1] == np.datetime64("2024-03-01T02:00:00", "ms")
    assert r.values["foF2"][0] == pytest.approx(8.400)
    assert r.values["hmF2"][0] == pytest.approx(270.9)
    # 00:45 has foF2 but no profile: '---' is NaN and the QD code is '__'.
    assert r.values["foF2"][6] == pytest.approx(8.350)
    assert np.isnan(r.values["hmF2"][6]) and r.qd["hmF2"][6] == "__"
    assert r.cs[6] == 0


def test_didbase_confidence_mask_passes_manual_scaling() -> None:
    from tss.data.giro import parse_scaled

    r = parse_scaled((FIX / "giro" / "didbase_MHJ45_20240301.txt").read_bytes())
    assert r.good(70).sum() == int(((r.cs >= 70)).sum())
    from dataclasses import replace

    manual = replace(r, cs=np.full(len(r), 999))
    assert manual.good(100).all()


def test_didbase_no_data_is_an_empty_answer_not_an_error() -> None:
    from tss.data.giro import parse_scaled

    r = parse_scaled((FIX / "giro" / "didbase_nodata.txt").read_text())
    assert len(r) == 0
    assert "No measurement data" in r.status


def test_didbase_unknown_characteristic_is_an_error() -> None:
    from tss.data.giro import GiroQueryError, parse_scaled

    with pytest.raises(GiroQueryError, match="No known requested characteristic"):
        parse_scaled((FIX / "giro" / "didbase_badchar.txt").read_text())


def test_didbase_rejects_something_that_is_not_a_table() -> None:
    from tss.data.giro import GiroQueryError, parse_scaled

    with pytest.raises(GiroQueryError):
        parse_scaled("<html><body>Service unavailable</body></html>")


def test_didbase_query_uses_names_and_the_service_date_format() -> None:
    from tss.data.giro import build_query

    q = dict(build_query("mhj45", "2024-03-01T00:00", "2024-03-01T02:00", ["foF2", "hmF2", "B0", "B1"]))
    assert q == {
        "ursiCode": "MHJ45",
        "charName": "foF2,hmF2,B0,B1",
        "DMUF": "3000",
        "fromDate": "2024/03/01 00:00:00",
        "toDate": "2024/03/01 02:00:00",
    }


def test_didbase_request_url_matches_the_verified_request() -> None:
    """The exact URL that returned the fixture's body on 2026-09-23."""
    from tss.data.giro import build_query, request_url

    url = request_url(build_query("MHJ45", "2024-03-01T00:00", "2024-03-01T02:00",
                                  ["foF2", "hmF2", "B0", "B1", "fmin"]))
    assert url == (
        "https://lgdc.uml.edu/fastchar/getbest?ursiCode=MHJ45&charName=foF2,hmF2,B0,B1,fmin"
        "&DMUF=3000&fromDate=2024%2F03%2F01+00%3A00%3A00&toDate=2024%2F03%2F01+02%3A00%3A00"
    )


def test_didbase_friendly_names_map_to_server_spellings() -> None:
    from tss.data.giro import CHARACTERISTICS, SERVER_NAMES, build_query, request_url

    q = dict(build_query("MHJ45", "2024-03-01", "2024-03-02", ["MUFD", "M(D)", "h'F2", "h`F2", "foF2p"]))
    assert q["charName"] == "MUF(D),M(D),h`F2,foF2p"  # duplicates sent once
    assert "charName=MUF%28D%29,M%28D%29,h%60F2,foF2p" in request_url(build_query(
        "MHJ45", "2024-03-01", "2024-03-02", ["MUFD", "MD", "h'F2", "foF2p"]))
    assert set(CHARACTERISTICS.values()) == set(SERVER_NAMES)
    assert "Type Es" not in SERVER_NAMES and "TypeEs" not in CHARACTERISTICS


@pytest.mark.parametrize(
    "args,match",
    [
        (("MILLSTONE", "2024-03-01", "2024-03-02"), "URSI"),
        (("MHJ45", "2024-03-02", "2024-03-01"), "not after"),
        (("MHJ45", "2024-03-01", "2024-03-02", ["foF3"]), "unknown characteristic"),
    ],
)
def test_didbase_form_refuses_bad_requests(args, match) -> None:
    from tss.data.giro import build_query

    with pytest.raises(ValueError, match=match):
        build_query(*args)


def test_didbase_western_pacific_codes_are_well_formed() -> None:
    from tss.data.giro import WESTERN_PACIFIC, build_query

    for code in WESTERN_PACIFIC:
        build_query(code, "2024-01-01", "2024-01-02")


def test_didbase_query_chunks_caches_and_does_not_refetch(tmp_cache, monkeypatch) -> None:
    from tss.data import cache, giro

    body = (FIX / "giro" / "didbase_MHJ45_20240301.txt").read_bytes()
    calls: list[str] = []

    def fake_get(url, **kw):
        assert "params" not in kw, "the URL is sent pre-encoded, not rebuilt by requests"
        calls.append(url)
        return _Resp(body)

    monkeypatch.setattr(cache, "polite_get", fake_get)

    result, paths = giro.query_scaled("MHJ45", "2024-01-01", "2024-03-11")
    # 70 days at 31 days per chunk is three requests, one GET each, no redirect.
    assert len(calls) == 3 and len(paths) == 3
    assert all(u.startswith("https://lgdc.uml.edu/fastchar/getbest?ursiCode=MHJ45&") for u in calls)
    assert "fromDate=2024%2F01%2F01+00%3A00%3A00&toDate=2024%2F02%2F01+00%3A00%3A00" in calls[0]
    # The same fixture served three times collapses to one set of timestamps.
    assert len(result) == 17
    for p, u in zip(paths, calls, strict=True):
        meta = json.loads(cache.meta_path(p).read_text())
        assert meta["source"] == "giro_didbase"
        assert meta["url"] == u, "the sidecar records exactly the URL that was sent"
        assert meta["extra"]["params"]["charName"] == "foF2,hmF2,B0,B1,fmin"

    giro.query_scaled("MHJ45", "2024-01-01", "2024-03-11")
    assert len(calls) == 3, "a cached window must not be fetched again"

    _r, capped = giro.query_scaled("MHJ45", "2024-01-01", "2024-03-11", max_chunks=1, force=True)
    assert len(capped) == 1 and len(calls) == 4
    assert "truncated" in _r.status


def test_didbase_source_notes_flag_the_non_commercial_licence() -> None:
    from tss.data import SOURCES

    src = SOURCES["giro_didbase"]
    assert "CC-BY-NC-SA" in src.licence
    assert "Navy" in src.notes


# ---------------------------------------------------------------------------
# IRTAM
# ---------------------------------------------------------------------------

ARCHIVE = FIX / "irtam" / "IRTAM_foF2_COEFFS_20121017_233000.ASC"
SERVLET_FILE = FIX / "irtam" / "servlet_foF2_20240301T1200.txt"


@pytest.mark.parametrize("path", [ARCHIVE, SERVLET_FILE])
def test_irtam_parses_both_routes(path) -> None:
    from tss.data.irtam import parse_coefficients

    c = parse_coefficients(path.read_bytes())
    assert c.characteristic == "foF2" and c.unit == "MHz"
    assert (c.n_temporal, c.n_spatial) == (14, 76)
    assert c.coefficients.shape == (1064,)
    assert c.basis == "JonesGallet_LinTrend"
    assert c.grid == "46 lats x 45 lons"
    assert "MHJ45" in c.stations
    assert all(len(s) == 5 for s in c.stations)


def test_irtam_archive_header_values() -> None:
    from datetime import datetime, timezone

    from tss.data.irtam import parse_coefficients

    c = parse_coefficients(ARCHIVE.read_bytes())
    assert c.valid_time == datetime(2012, 10, 17, 23, 30, tzinfo=timezone.utc)
    assert c.engine == "NECTAR v0.2A"
    assert len(c.stations) == 37
    assert c.coefficients[0] == pytest.approx(8.90280741)


def test_irtam_matrix_follows_the_pyirtam_layout() -> None:
    """Layout checked against PyIRTAM 0.0.7 IRTAM_read_files on both fixtures.

    First 988 numbers: 13 x 76 in Fortran order (temporal index fastest).
    Last 76: the trend term, inserted as row 1.
    """
    from tss.data.irtam import parse_coefficients

    c = parse_coefficients(ARCHIVE.read_bytes())
    v, m = c.coefficients, c.matrix
    assert m.shape == (14, 76)
    np.testing.assert_array_equal(m[1], v[988:])
    np.testing.assert_array_equal(m[0], v[0:988:13])
    np.testing.assert_array_equal(m[2], v[1:988:13])
    np.testing.assert_array_equal(m[13], v[12:988:13])
    # The physical signature that makes the layout credible: a large mean term,
    # a trend term three orders of magnitude smaller.
    assert np.median(np.abs(m[0])) > 100 * np.median(np.abs(m[1]))


def test_irtam_refuses_a_truncated_file() -> None:
    from tss.data.base import FetchError
    from tss.data.irtam import parse_coefficients

    text = ARCHIVE.read_text()
    with pytest.raises(FetchError, match="truncated"):
        parse_coefficients(text[: len(text) // 2])


def test_irtam_refuses_an_html_error_page() -> None:
    from tss.data.base import FetchError
    from tss.data.irtam import parse_coefficients

    with pytest.raises(FetchError, match="not an IRTAM"):
        parse_coefficients("<!doctype html><html>404</html>")


def test_irtam_urls_and_rounding() -> None:
    from tss.data.irtam import archive_url, round_to_cadence, servlet_url

    assert archive_url("B0", "2012-10-17T23:37") == (
        "https://ulcar.uml.edu/GAMBIT/GambitCoefficients/COEFFS_2012/2012_10_17/"
        "IRTAM_B0in_COEFFS_20121017_233000.ASC"
    )
    assert archive_url("foF2", "2012-10-17T23:30:00Z").endswith(
        "IRTAM_foF2_COEFFS_20121017_233000.ASC"
    )
    assert servlet_url("hmF2", "2024-03-01T11:53") == (
        "https://lgdc.uml.edu/rix/gambit-coeffs?charName=hmF2&time=2024-03-01T12:00"
    )
    assert servlet_url("B1", "2024-03-01T12:00").endswith("charName=B1&time=2024-03-01T12:00")
    assert round_to_cadence("2024-03-01T23:53").isoformat() == "2024-03-02T00:00:00+00:00"
    assert round_to_cadence("2024-03-01T12:07:29").minute == 0
    assert round_to_cadence("2024-03-01T12:07:30").minute == 15
    with pytest.raises(ValueError, match="unknown IRTAM"):
        archive_url("NmF2", "2024-03-01")


def test_irtam_auto_route_falls_back_to_the_servlet_only_on_404(tmp_cache, monkeypatch) -> None:
    from tss.data import cache, irtam
    from tss.data.base import FetchError

    seen: list[str] = []

    def fake_get(url, **kw):
        seen.append(url)
        if "ulcar.uml.edu" in url:
            raise FetchError(f"GET {url} returned 404 Not Found; not retrying a client error.")
        return _Resp(SERVLET_FILE.read_bytes())

    monkeypatch.setattr(cache, "polite_get", fake_get)
    coeffs, path = irtam.fetch_coefficients("foF2", "2024-03-01T12:00")
    assert len(seen) == 2 and "gambit-coeffs" in seen[1]
    assert coeffs.valid_time.year == 2024
    meta = json.loads(cache.meta_path(path).read_text())
    assert meta["extra"]["route"] == "servlet"

    # Cached: a second call for a past epoch makes no request.
    irtam.fetch_coefficients("foF2", "2024-03-01T12:00")
    assert len(seen) == 2


def test_irtam_does_not_fall_back_on_a_server_error(tmp_cache, monkeypatch) -> None:
    from tss.data import cache, irtam
    from tss.data.base import FetchError

    def fake_get(url, **kw):
        raise FetchError(f"GET {url} returned 503; giving up after 5 attempts")

    monkeypatch.setattr(cache, "polite_get", fake_get)
    with pytest.raises(FetchError, match="503"):
        irtam.fetch_coefficients("foF2", "2012-10-17T23:30")


def test_irtam_does_not_cache_a_malformed_body(tmp_cache, monkeypatch) -> None:
    from tss.data import cache, irtam
    from tss.data.base import FetchError

    monkeypatch.setattr(cache, "polite_get", lambda url, **kw: _Resp(b"<html>oops</html>"))
    with pytest.raises(FetchError):
        irtam.fetch_coefficients("foF2", "2012-10-17T23:30", route="archive")
    assert not any(tmp_cache.rglob("*.ASC"))


# ---------------------------------------------------------------------------
# Madrigal
# ---------------------------------------------------------------------------

MLH_PARMS = ("YEAR", "MONTH", "DAY", "HOUR", "MIN", "SEC", "GDLAT", "GLON", "GDALT", "NE", "DNE")
DMSP_PARMS = ("YEAR", "MONTH", "DAY", "HOUR", "MIN", "SEC", "GDLAT", "GLON", "GDALT", "NE")


@pytest.fixture()
def madrigal_identity(monkeypatch):
    monkeypatch.setenv("TSS_MADRIGAL_NAME", "Test Person")
    monkeypatch.setenv("TSS_MADRIGAL_EMAIL", "test.person@example.org")
    monkeypatch.setenv("TSS_MADRIGAL_AFFILIATION", "Example University")


def test_madrigal_experiments_parse_the_live_field_order() -> None:
    from datetime import datetime, timezone

    from tss.data.madrigal import parse_experiments

    exps = parse_experiments((FIX / "madrigal" / "experiments_dmsp_20240229_0302.txt").read_bytes())
    assert len(exps) == 6
    e = exps[2]
    assert e.id == 100288979 and e.inst_code == 8100 and e.site_name == "CEDAR"
    assert e.start == datetime(2024, 3, 1, tzinfo=timezone.utc)
    assert e.end == datetime(2024, 3, 2, tzinfo=timezone.utc)
    assert e.is_local and e.pi == "Keith Groves"
    assert e.path == "experiments3/2024/dms/01mar24"


def test_madrigal_files_keep_only_final_defaults() -> None:
    from tss.data.madrigal import KINDAT, parse_files

    text = (FIX / "madrigal" / "files_mlh_100004470.txt").read_bytes()
    files = parse_files(text)
    assert {f.category for f in files} == {1}
    assert {f.status for f in files} == {"Final"}
    zen = [f for f in files if f.kindat == KINDAT["mlh_zenith_single_pulse"]]
    assert [f.basename for f in zen] == ["mlh230920i.004.hdf5"]
    assert zen[0].doi.startswith("https://w3id.org/cedar")
    assert len(parse_files(text, include_non_default=True)) == 11


def test_madrigal_instruments_parse() -> None:
    from tss.data.madrigal import INSTRUMENTS, parse_instruments

    inst = {i.code: i for i in parse_instruments((FIX / "madrigal" / "instruments_sample.txt").read_bytes())}
    assert inst[INSTRUMENTS["millstone"]].lat_deg == pytest.approx(42.619)
    assert inst[INSTRUMENTS["dmsp"]].category == "Satellite Instruments"


def test_madrigal_isprint_parses_dmsp_with_nan() -> None:
    from tss.data.madrigal import parse_isprint

    t = parse_isprint((FIX / "madrigal" / "isprint_dmsp_f16_20240301.txt").read_bytes(), DMSP_PARMS)
    assert len(t) == 20
    assert np.isnan(t["NE"][:7]).all()
    assert t["NE"][7] == pytest.approx(5.41610e10)
    assert t["GDALT"][0] == pytest.approx(872.20)
    assert t.times[0] == np.datetime64("2024-03-01T00:00:00")
    assert t["ne"] is t["NE"]


def test_madrigal_isprint_error_body_raises_even_on_http_200() -> None:
    from tss.data.madrigal import MadrigalServiceError, parse_isprint

    with pytest.raises(MadrigalServiceError, match="Illegal parameter"):
        parse_isprint((FIX / "madrigal" / "isprint_error.txt").read_bytes(), DMSP_PARMS)


def test_madrigal_isprint_column_count_mismatch_raises() -> None:
    from tss.data.madrigal import MadrigalServiceError, parse_isprint

    with pytest.raises(MadrigalServiceError, match="fields"):
        parse_isprint((FIX / "madrigal" / "isprint_dmsp_f16_20240301.txt").read_bytes(), MLH_PARMS)


def test_madrigal_profiles_group_the_millstone_topside() -> None:
    from tss.data.madrigal import parse_isprint, profiles

    t = parse_isprint((FIX / "madrigal" / "isprint_mlh_zenith_20230920.txt").read_bytes(), MLH_PARMS)
    profs = profiles(t)
    assert len(profs) == 2
    clean = profs[1]
    assert clean["alt_km"][0] < 200 and clean["alt_km"][-1] > 800
    peak = clean["alt_km"][np.argmax(clean["ne_m3"])]
    assert 270 < peak < 310
    assert np.all(np.diff(clean["alt_km"]) > 0)


def test_madrigal_error_bar_filter_is_not_sufficient_on_its_own() -> None:
    """The documented limitation: the 645 km spike survives a 50 % error cut."""
    from tss.data.madrigal import parse_isprint, profiles

    t = parse_isprint((FIX / "madrigal" / "isprint_mlh_zenith_20230920.txt").read_bytes(), MLH_PARMS)
    first = profiles(t, max_rel_err=0.5)[0]
    i = int(np.argmin(np.abs(first["alt_km"] - 645.0)))
    assert first["ne_m3"][i] > 4e11, "the spike should still be there"


def test_madrigal_isprint_params_repeat_parms_and_build_the_date_filter() -> None:
    from tss.data.madrigal import build_isprint_params

    pairs = build_isprint_params("/x/f.hdf5", ["year", "ne"], "2024-03-01T00:00:00", "2024-03-01T00:01:00")
    assert [v for k, v in pairs if k == "parms"] == ["YEAR", "NE"]
    assert dict(pairs)["filters"] == "date1=03/01/2024 time1=00:00:00 date2=03/01/2024 time2=00:01:00"


def test_madrigal_encoding_matches_the_request_verified_live() -> None:
    """The exact form of the isprint URL that returned data on 2026-09-23."""
    from tss.data import madrigal

    pairs = madrigal.build_isprint_params(
        "/opt/openmadrigal/madroot/experiments3/2024/dms/01mar24/dms_20240301_16s1.001.hdf5",
        ["YEAR", "NE"], "2024-03-01T00:00:00", "2024-03-01T00:01:00",
        ident={"user_fullname": "A B. C", "user_email": "a@b.org", "user_affiliation": "U S"},
    )
    q = madrigal._encode(pairs)
    assert "file=/opt/openmadrigal/madroot/experiments3/2024/dms/01mar24/dms_20240301_16s1.001.hdf5" in q
    assert "&parms=YEAR&parms=NE&" in q
    assert "filters=date1%3D03%2F01%2F2024%20time1%3D00:00:00%20date2%3D03%2F01%2F2024%20time2%3D00:01:00" in q
    assert "user_fullname=A%2BB.%2BC" in q
    assert "user_email=a@b.org" in q
    assert "user_affiliation=U%2BS" in q


def test_madrigal_identity_missing_is_an_explained_error(monkeypatch) -> None:
    from tss.data.madrigal import MadrigalIdentityError, identity

    for k in ("TSS_MADRIGAL_NAME", "TSS_MADRIGAL_EMAIL", "TSS_MADRIGAL_AFFILIATION"):
        monkeypatch.delenv(k, raising=False)
    with pytest.raises(MadrigalIdentityError, match="TSS_MADRIGAL_EMAIL"):
        identity()


def test_madrigal_isprint_never_writes_the_identity_to_the_cache(
    tmp_cache, monkeypatch, madrigal_identity
) -> None:
    from tss.data import cache, madrigal

    sent: list[str] = []

    def fake_get(url, **kw):
        sent.append(url)
        return _Resp((FIX / "madrigal" / "isprint_mlh_zenith_20230920.txt").read_bytes())

    monkeypatch.setattr(cache, "polite_get", fake_get)
    table, path = madrigal.isprint("/x/mlh230920i.004.hdf5", MLH_PARMS, "2023-09-20T22:00", "2023-09-20T22:10")
    assert len(table) == 76

    # The request carried the identity...
    assert "test.person@example.org" in sent[0]
    # ...and nothing on disk does: not the file name, not the sidecar.
    for f in tmp_cache.rglob("*"):
        if f.is_file():
            assert "example.org" not in f.name and "Test" not in f.name
    meta_text = cache.meta_path(path).read_text()
    assert "example.org" not in meta_text and "Test" not in meta_text
    assert "user_email=[redacted]" in json.loads(meta_text)["url"]

    madrigal.isprint("/x/mlh230920i.004.hdf5", MLH_PARMS, "2023-09-20T22:00", "2023-09-20T22:10")
    assert len(sent) == 1, "an extract is cached permanently"


def test_madrigal_listing_needs_no_identity(tmp_cache, monkeypatch) -> None:
    from tss.data import cache, madrigal

    for k in ("TSS_MADRIGAL_NAME", "TSS_MADRIGAL_EMAIL", "TSS_MADRIGAL_AFFILIATION"):
        monkeypatch.delenv(k, raising=False)
    sent: list[str] = []

    def fake_get(url, **kw):
        sent.append(url)
        return _Resp((FIX / "madrigal" / "experiments_dmsp_20240229_0302.txt").read_bytes())

    monkeypatch.setattr(cache, "polite_get", fake_get)
    exps = madrigal.list_experiments("dmsp", "2024-02-29", "2024-03-03")
    assert len(exps) == 6
    assert "code=8100" in sent[0] and "startyear=2024" in sent[0] and "user_" not in sent[0]


def test_madrigal_service_error_is_not_cached(tmp_cache, monkeypatch, madrigal_identity) -> None:
    from tss.data import cache, madrigal

    monkeypatch.setattr(
        cache, "polite_get",
        lambda url, **kw: _Resp((FIX / "madrigal" / "isprint_error.txt").read_bytes()),
    )
    with pytest.raises(madrigal.MadrigalServiceError):
        madrigal.isprint("/x/f.hdf5", ["YEAR"])
    assert not [p for p in tmp_cache.rglob("*.txt")]


def test_didbase_empty_answer_still_has_the_requested_columns(tmp_cache, monkeypatch) -> None:
    """Found live: Guam (GU513) had no data on 2024-03-01, and the result had no columns."""
    from tss.data import cache, giro

    body = (FIX / "giro" / "didbase_nodata.txt").read_bytes()
    monkeypatch.setattr(cache, "polite_get", lambda url, **kw: _Resp(body))
    r, _paths = giro.query_scaled("GU513", "2024-03-01", "2024-03-02", ["foF2", "hmF2"])
    assert len(r) == 0 and "No measurement data" in r.status
    assert r.values["foF2"].shape == (0,) and r.values["hmF2"].shape == (0,)
    assert r.good().shape == (0,)
    # Columns are keyed by the server's spelling, whatever alias was requested.
    r2, _ = giro.query_scaled("GU513", "2024-03-01", "2024-03-02", ["MUFD", "h'F2"])
    assert set(r2.values) == {"MUF(D)", "h`F2"}


def test_madrigal_isprint_no_records_is_an_empty_table() -> None:
    """Captured live: a window outside the file's span returns this line with HTTP 200."""
    from tss.data.madrigal import parse_isprint

    t = parse_isprint((FIX / "madrigal" / "isprint_no_records.txt").read_bytes(), DMSP_PARMS)
    assert len(t) == 0
    assert t.columns == DMSP_PARMS
    assert t["NE"].shape == (0,) and t.times.shape == (0,)
