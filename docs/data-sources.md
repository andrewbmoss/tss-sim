# External data sources

Three external data sources and a set of locally installed models feed this
repository.  This document says what each one provides, where it comes from,
what its licence and citation requirements are, how it is cached, how to be a
good citizen when fetching it, what is checked in versus downloaded, and which
work package consumes it.

The implementation is `src/tss/data/`, one module per source, with the shared
registry and cache in `base.py` and `cache.py`.  Registered source names are
what `tss fetch` takes on the command line:

```bash
tss fetch list                       # every source with its one-line description
tss fetch spdf_isis_ne
tss fetch wspr_live --limit 20000
tss fetch datagov_absorption --force
```

Two rules apply throughout.

**Tests never hit the network.**  Not once, not with a skip guard, not "only in
CI".  Every test that needs external data reads a small checked-in extract from
`tests/fixtures/`.  A fetcher is tested against a recorded response, not against
the live endpoint.  The reason is partly reproducibility and partly courtesy: a
test suite that hammers a public archive on every commit is a test suite that
gets the author blocked.

**Downloaded data is not checked in.**  `data/external/` is in `.gitignore`.
What is checked in is the fetcher, the fixture extract, and the manifest entry
recording the SHA-256 of every file a run consumed.  A result is reproducible
because the manifest says exactly which bytes went into it, not because the
bytes are in the repository.

## Cache layout

Every external file lands under `data/external/<source>/` with a
`<file>.meta.json` sidecar beside it.  The cache root can be moved with the
`TSS_DATA_CACHE` environment variable, which is what lets a large archive live
on a different drive without editing any configuration.

```
data/external/
  spdf_isis_ne/
    isis1_ne.asc                 raw bytes, exactly as served, never edited
    isis1_ne.asc.meta.json
    ...
  wspr_live/
    <query-hash>.csv
    <query-hash>.csv.meta.json
  datagov_absorption/
    ...
```

The sidecar records the URL the file came from, the UTC time of retrieval, the
server's ETag or Last-Modified if either was offered, the SHA-256 of the bytes
and the byte count.  That is the minimum needed for `tss.provenance.Manifest` to
state which version of an upstream data set produced a figure, and for a reader
to tell whether their copy is the same one.  `verify_cached` re-checks the hash.

Filenames are content-addressed in the weak sense that they carry a hash of the
request, not of the response.  Hashing the response would mean a file has no
name until it has been downloaded, which makes resume and cache lookup awkward.
The response hash goes in the sidecar instead.

Raw bytes are never edited.  A parsing fix re-parses from the cached file, which
means a parser bug can be corrected without re-downloading and without any doubt
about what was downloaded.

Each source carries its provenance as data on the `Source` record in
`src/tss/data/base.py`: `licence`, `citation`, `acknowledgement`, `cadence`,
`work_packages` and candid `notes`.  Neither `licence` nor `citation` may be
empty.  The citation string is copied into the run manifest, so a figure carries
its attribution without anyone having to remember it at writing-up time.

## NASA SPDF Alouette and ISIS topside sounder archive

Registered as **`spdf_isis_ne`**.  Module `src/tss/data/spdf_isis.py`.

**What it provides.**  The Alouette and ISIS spacecraft flew topside sounders
from 1962 to 1984 and are the only large body of real topside ionograms and real
topside electron density profiles in existence.  Two products matter here.

The hand-scaled electron density profiles under
`topside_sounder/crc_ne_profile_ascii/` were scaled at the Communications
Research Centre in Ottawa by Jackson's method.  `isis1_ne.asc` is 8.4 MB and
holds ISIS-1 profiles from 1969 to 1980 in a fixed-width ASCII format.  The profile count is not asserted here because this repository has never had the file to count; read_ne_ascii reports it at read time.  Each profile runs downward from
the spacecraft, so its first sample is the local electron density at the
satellite.  That single number is what sets the lower band edge of INS-01, which
makes this file a direct test of the requirement rather than an indirect one.
The delivered WP8 validation pairs 5,118 of these profiles against the model.

The digitised ionograms in the same archive are what a full V4 retrieval
comparison needs, in order to run the inversion against a real measurement
rather than a simulated one.  Those are large binary products and the module
does not fetch them automatically; `IONOGRAM_DIRS` records where they live.

**Where.**  https://spdf.gsfc.nasa.gov/isis/isis-status.html and the archive
directories it links to.

**Licence and citation.**  Public domain as a United States Government work,
with acknowledgement requested.  Cite the CRC topside sounder Ne profiles as
distributed by the NASA Space Physics Data Facility, and acknowledge both SPDF
and the Communications Research Centre hand-scaling.  Record the retrieval date,
because the holdings have grown over time as further films were digitised.  The
exact strings are on the `Source` record and reach the manifest automatically.

**Cache.**  `data/external/spdf_isis_ne/`.  A static archive from a closed
mission, so a complete cache never needs refreshing and `tss fetch` does nothing
when the checksums match.

**Rate limits and etiquette.**  No published rate limit.  `polite_get` sets a
User-Agent naming the project and the institution, enforces a minimum interval
between requests to the same host, retries 5xx with exponential backoff and does
not retry 4xx at all, because a 404 is an answer.  Requests are serial.  Do not
parallelise them; the archive is a public service run on a modest budget.

**Checked in versus downloaded.**  Downloaded: everything under
`data/external/spdf_isis_ne/`.  Checked in: `tests/fixtures/spdf_isis/`, a few
hundred lines of `isis1_ne.asc` in the real format, enough to test the parser
and the profile comparison offline.

**Consumed by.**  WP1, where the scaled profiles feed the `tabulated` L0 tier
and serve as truth to check the L0 models against.  WP8 and verification level
V4, for the archive comparison.  Scores INS-01.

**Known limitation.**  Alouette and ISIS were pulse sounders.  The archive
contains no range-Doppler coupling and no transmit-to-receive leakage, because
neither exists in a pulse system with a switched receiver.  V4 therefore
validates the ionospheric and retrieval half of the chain and cannot touch the
FMCW half.  See `docs/verification.md`.

## wspr.live

Registered as **`wspr_live`**.  Module `src/tss/data/wspr.py`.

**What it provides.**  WSPR is a weak-signal beacon mode.  A transmitter sends a
two-minute frame carrying its callsign, grid square and power at a known slot
time, and every receiver that decodes it reports the decode with a
signal-to-noise ratio referred to a 2.5 kHz bandwidth.  Several thousand
stations do this continuously across the HF bands, so the aggregate is a dense,
free, continuously updated record of which paths were open, at what frequency,
at what hour, at what received SNR, for a transmit power usually stated to
within a decibel.

Used here in two ways.  For WP2, it is an independent check on the diurnal and
frequency dependence of the absorption model: the ratio of night to day spot
counts on a band, the lowest band carrying spots at local noon, and the drop in
median SNR through sunrise are all consequences of absorption and all measurable
from these data without any calibration.  For WP12, it is a large body of
distributed receive measurements from which a noise reference can be
constrained.

**Where.**  wspr.live exposes a ClickHouse HTTP endpoint.  Queries are sent as
SQL over HTTPS and results come back in the requested format.

**Licence and citation.**  Reception reports are contributed by amateur
operators and redistributed by wspr.live for public use, with attribution to
wspr.live and WSPRnet expected.  Cite the wspr.live spot database and note the
query date range.  Spot records carry callsigns, which are public identifiers;
the repository has no reason to retain them past aggregation and drops them once
the aggregates a study consumes have been built.

**Cache.**  `data/external/wspr_live/`, keyed by a hash of the query text plus
its date range.  The query SQL is stored in the sidecar alongside its result, so
a result can always be traced back to the exact SQL that produced it.  Cadence
is live, so a cached query is a snapshot and its retrieval time matters.

**Rate limits and etiquette.**  This is the source where etiquette matters most,
because it is a single volunteer-run endpoint and a careless query can ask it
for hundreds of millions of rows.  wspr.live explicitly asks that bulk users
take the published archive files rather than hammer the live endpoint, and the
module carries that note on its `Source` record.  The conventions:

* Always bound the query by time.  Never issue a bare `SELECT ... FROM wspr.rx`
  without a date predicate.
* For anything beyond a few weeks, take the published archive files instead of
  querying.  That is what the provider asks for and it is faster.
* Aggregate server-side.  If a study needs a daily median, ask ClickHouse for
  the daily median rather than pulling raw spots and computing it locally.
* Fetch in day-sized or week-sized chunks and cache each chunk, so an
  interrupted collection resumes rather than restarts.
* One request at a time, with a minimum interval between them.  No parallel
  fetching.
* Prefer a long collection at a polite rate over a short one at an impolite
  rate.  WP12 has a long lead time for exactly this reason, which is why
  `docs/work-packages.md` recommends starting the collection early.

**Checked in versus downloaded.**  Downloaded: everything under
`data/external/wspr_live/`.  Checked in: `tests/fixtures/wspr/`, one recorded
query response of a few hundred rows, used to test the parser, the aggregation
and the chunking logic offline.

**Consumed by.**  WP2 as an absorption sanity check, and WP12 for the
sub-megahertz noise reference, and so indirectly C1.

**Known limitation, stated early rather than discovered late.**  WSPR is
ground-to-ground oblique propagation, and it does not operate below the HF
amateur bands.  It cannot by itself supply a measurement at 0.33 MHz.  What it
can supply is a bound and a constraint at the frequencies where it does operate,
which narrows the extrapolation below 0.6 MHz rather than removing it.  If WP12
concludes with a constrained extrapolation rather than a measurement, the band
edge result stays conditional and has to be written up that way.

## data.gov ionospheric absorption measurements

Registered as **`datagov_absorption`**.  Module
`src/tss/data/datagov_absorption.py`.

**What it provides.**  Measurements of the absorptive properties of the
ionosphere, used as an independent check on the absorption model delivered in
WP2.

**Where.**
https://catalog.data.gov/dataset/measurements-of-the-absorptive-properties-of-the-ionosphere
Resource URLs are discovered through the CKAN API rather than hard-coded,
because data.gov resource links move and the catalogue entry is the stable
identifier.

**Licence and citation.**  United States Government work, public domain; the
catalogue entry carries the exact rights statement.  Cite the data.gov catalogue
entry by its dataset identifier and acknowledge the originating collection,
which is NOAA National Centers for Environmental Information and its
predecessors.  Record the catalogue metadata alongside the data, because entries
are revised and the catalogue page is not a stable copy of the data.

**Cache.**  `data/external/datagov_absorption/`, with the CKAN catalogue record
stored beside the data so the version is recoverable.  A static archive whose
catalogue metadata may be revised.

**Rate limits and etiquette.**  A small, one-time download through the CKAN API.
Fetch once and cache.  There is no reason to poll.

**Checked in versus downloaded.**  Downloaded: the dataset and its catalogue
record.  Checked in: a small extract under `tests/fixtures/absorption/`,
sufficient to test the parser and one comparison case.

**Consumed by.**  WP2, and WP8 at verification level V3, as one of the
independent absorption references alongside Barry 1971.

## Local model installations

Module `src/tss/data/local_models.py`.  Nothing here downloads anything.
PHaRLAP is not redistributable, NeQuick is distributed under its own terms,
MATLAB is licensed, and the IRI bindings need a Fortran compiler.  What the
module does is answer one question honestly for each of them: can L0 actually
call it right now on this machine.

Availability is decided by a trial call, not by the presence of a file.  That is
deliberate.  Testing for PHaRLAP by checking that a file exists has already
produced a false negative once on the author's machine, and a false negative
here costs an afternoon.  Having one place that answers the availability
question means a tier does not carry its own detection code, and a campaign can
fail at configuration time with a message a person can act on, rather than at
hour three of a run with an import error.

**IRI-2020.**  The International Reference Ionosphere, the primary L0 truth
model.  Reached from Python through the installed bindings and from MATLAB
through the batch generator.  Consumed by WP1, and by the band-edge study over
1.2 million profiles in WP0 and WP1.  Cite the IRI-2020 reference publication
and state the options used, the topside option in particular, since the topside
formulation is exactly what a topside sounder study is sensitive to.  The
`IonosphereState.source` field carries that option string into every manifest.

**NeQuick 2.**  A second, independent topside formulation.  Consumed by WP1, and
used at verification level V2 to show whether a result depends on the topside
model rather than on the instrument.  Not yet wired in; this is the outstanding
half of WP1.

**PHaRLAP 4.7.4.**  Numerical HF ray tracing, installed at
`C:\Users\Andrew\PhD\Chartier Verification\pharlap_4.7.4`.  MATLAB only, which
is one of the two reasons MATLAB is a peer implementation rather than a
convenience; see `docs/adr/0002-python-primary-matlab-peer.md`.  Consumed by WP2
for the ray-traced L1 tier and by V2 as the cross-fidelity reference against the
analytic tier.  Licensed by Defence Science and Technology Group and not
redistributed.  Its terms govern; check them before publishing derived results.

Paths to these installations live in the `paths:` block of
`config/default.yaml`.  No module contains a literal path.  Any test or study
needing one carries `@pytest.mark.pharlap`, `@pytest.mark.matlab` or
`@pytest.mark.data` and is skipped when the dependency is absent.  A skipped
test reports as skipped rather than passing quietly, so a run on a machine
without PHaRLAP does not look like a clean full run.

## Adding a source

1. Write a module under `src/tss/data/` and register it with
   `@register_source`, filling in every field.  `licence` and `citation` may not
   be empty, and `notes` is the place for candid remarks about limitations and
   gaps.  Importing the module is what populates `SOURCES`, exactly as importing
   a layer package populates the tier registry.
2. Fetch through `polite_get` and write through the cache helpers, so the
   sidecar, the checksum and the rate limiting come for free.
3. Check in a small representative extract under `tests/fixtures/<name>/` and
   write the parser tests against it.  No network in tests.
4. Add a section to this document under the same headings used above: what it
   provides, where, licence and citation, cache, rate limits and etiquette,
   checked in versus downloaded, which work package consumes it, and any known
   limitation.
5. Make sure the study that consumes it calls `Manifest.add_input` for every
   file it reads, so the checksums and the citation reach the manifest.
