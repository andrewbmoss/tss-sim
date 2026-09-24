# External data sources

Six external data sources and a set of locally installed models feed this
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
tss fetch giro_didbase               # one day at Millstone Hill
tss fetch irtam                      # one epoch, all four characteristics
tss fetch madrigal                   # needs TSS_MADRIGAL_* set; see below
```

The first three sources were chosen for the proposal.  Madrigal, GIRO DIDBase
and IRTAM were added on 2026-09-23 after an evaluation of candidate ground-truth
sources (Dissertation project note `data-sources-evaluation.md`).  They answer a
gap the first three leave: present-day truth for the topside profile, the F2
peak and the density at the spacecraft, where ISIS supplies only 1969 to 1980.

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

## CEDAR Madrigal: incoherent scatter radar and DMSP

Registered as **`madrigal`**.  Module `src/tss/data/madrigal.py`.

**What it provides.**  Two products from the CEDAR Madrigal database.
Incoherent scatter radar electron density profiles, of which the Millstone Hill
zenith single-pulse product (kindat 3420) is the most useful: it runs from about
150 km through the F2 peak to about 800 km with a per-point error estimate.
That is measured topside profile shape, which no ionosonde and no GNSS product
provides.  And DMSP SSIES in-situ electron density at about 850 km at
one-second cadence, which is the local plasma frequency at a spacecraft in the
current solar cycle, the modern counterpart to the first sample of an ISIS
profile.

**Where.**  https://cedar.openmadrigal.org/ through its plain-text services
`getInstrumentsService.py`, `getExperimentsService.py`,
`getExperimentFilesService.py` and `isprintService.py`.  These are what the
`madrigalWeb` client calls.  The module calls them directly through
`polite_get` so caching, sidecars and throttling match every other source and
`madrigalWeb` is not a dependency.  `isprint` asks the server to extract only
the parameters and time window wanted, which returns kilobytes of text instead
of an HDF5 file of tens of megabytes.

**Identity.**  Madrigal requires a name, email and affiliation on every data
request.  Listings do not.  Set them in the environment, never in a file in the
repository:

```bash
export TSS_MADRIGAL_NAME="Andrew B. Moss"
export TSS_MADRIGAL_EMAIL="..."
export TSS_MADRIGAL_AFFILIATION="Utah State University"
```

On Windows use `setx` or the conda environment's `env_vars`.  A data request
without them raises `MadrigalIdentityError` rather than going out anonymously.
The URL recorded in a sidecar has all three values replaced by `[redacted]`,
and the cache key is built from the redacted URL, so the identity reaches
neither a filename nor a committed sidecar.

**Licence and citation.**  CEDAR Madrigal data are openly available under the
CEDAR Rules of the Road: acknowledge the database, each instrument and its
funding, and consult the instrument PI where a result depends heavily on their
data.  Each file listing carries a DOI-style `w3id.org/cedar` URL, which is the
citation for that file.  The PI of each experiment is in the listing.

**Cache.**  `data/external/madrigal/`.  Listings expire after 30 days
(instruments after 90), because files move from preliminary to final.  An
isprint extract is cached permanently, since it names a specific file version.

**Rate limits and etiquette.**  No published limit.  Two seconds between
requests.  Use `isprint` with a narrow window rather than downloading whole
files.  `list_files` keeps only the final default files unless asked otherwise,
because scoring against a superseded calibration is an easy silent mistake.

**Behaviours the parser exists to handle.**  `isprintService.py` takes one
`parms=` per parameter and rejects a comma-joined list.  It answers HTTP 200 on
failure with a Python traceback in the body; `parse_isprint` detects that and
raises, and an error body is never cached.  A window that matches nothing
returns "No records were selected with the filters above", also with HTTP 200,
and is read as an empty table.  The output has no header line, so
columns are named by the requested parameters in order and a row with the
wrong field count is an error.

**Checked in versus downloaded.**  Checked in under `tests/fixtures/madrigal/`,
all captured live on 2026-09-23: a DMSP experiment listing, two file listings,
an instrument sample, 20 s of DMSP F16, two Millstone Hill zenith profiles from
2023-09-20 22:01 UT, the error body returned for a malformed request, and the
body returned for a window with no records.

**Consumed by.**  WP1 as topside truth for the L0 models, and WP8 at V4 as the
modern complement to the ISIS comparison.  DMSP bears on INS-01.

**Known limitations.**  ISR sites are few and radars run on campaign days, so
coverage is sparse and mid- or high-latitude except for Jicamarca.  DMSP flies
at about 850 km in dawn-dusk and similar orbits, higher than a 500 km cubesat,
so its density is an upper-altitude, fixed-local-time sample.  Error bars alone
do not clean ISR topside data: the Millstone fixture's first profile has a point
at 645 km three times its neighbours with a stated error of 2 %.  A plausibility
check belongs in the study that uses the profiles.

## GIRO DIDBase: scaled Digisonde characteristics

Registered as **`giro_didbase`**.  Module `src/tss/data/giro.py`.

**What it provides.**  Scaled characteristics from about 120 Digisondes:
foF2, hmF2, B0, B1 and fmin by default, with the autoscaling confidence score
(0 to 100, 999 for manual scaling) and a qualifying code per value.  Millstone
Hill reports every 7.5 minutes.

**What it is truth for.**  The F2 peak only.  A topside retrieval ends at hmF2
and an ionosonde measures foF2 and hmF2 from below, so a satellite pass near a
Digisonde checks the bottom of the retrieved profile independently.  Nothing
above hmF2 is visible to an ionosonde.

**Where.**  The GET service `https://lgdc.uml.edu/fastchar/getbest`, with
`ursiCode`, a comma-separated `charName` list, `DMUF` and a UTC window as
`fromDate`/`toDate` in `YYYY/MM/DD hh:mm:ss`.  It returns the plain-text table
directly.  This is where the DIDBase web form
(https://giro.uml.edu/didbase/scaled.php) redirects; calling it directly saves
the redirect and uses characteristic names rather than the form's positional
option indices, which would have broken silently had the form been reordered.
Both routes return the same body (checked 2026-09-23).  The service spells the
virtual heights with a backtick (`h`F2`) and uses `MUF(D)` and `M(D)`;
`CHARACTERISTICS` also accepts `h'F2`, `MUFD` and `MD`.  `Type Es` is not
supported: it is a classification, and the service drops it from the table.
The older GET service `lgdc.uml.edu/common/DIDBGetValues`, used by most
published scripts, answered 404 "resource not available" on 2026-09-23 and
should not be reintroduced.

**Licence and citation.**  Every response states that GIRO measurements are
released under CC-BY-NC-SA 4.0, that the Lowell GIRO Data Center Rules of the
Road apply, and that the data provider of each station used must be
acknowledged.  Cite Reinisch and Galkin (2011), Earth Planets Space 63, 377,
and name each station's provider.  The non-commercial clause is fine for the
dissertation; confirm with the Lowell GIRO Data Center before using these data
in Navy-sponsored work.

**Cache.**  `data/external/giro_didbase/`, one file per 31-day chunk, keyed on
the request URL.  The sidecar records that URL, exactly as sent, and the parameters.

**Rate limits and etiquette.**  Three seconds between requests.  Windows longer
than 31 days are split into chunks and each is cached, so a long collection
resumes rather than restarts.

**Privacy.**  The server writes the requester's public IP address into the
response header.  Cached files keep it because raw bytes are never edited, and
`data/external/` is not committed except for sidecars, which do not contain the
body.  Fixtures have the line redacted.

**Checked in versus downloaded.**  `tests/fixtures/giro/`: two hours at
Millstone Hill on 2024-03-01, the no-data response, and the unknown-characteristic
error response.

**Consumed by.**  WP1 and WP8, for peak conjunctions.  `WESTERN_PACIFIC` lists
the stations relevant to the bistatic side study (Guam, Kwajalein, Okinawa,
Wake, Darwin, Townsville, Learmonth, Niue, Norfolk, Cocos).

**Known limitations.**  Autoscaled values vary in quality and the confidence
score is the only guard; `ScaledCharacteristics.good(min_cs)` applies it, and
the threshold should be chosen and stated per study.  Station coverage is
uneven and thin over oceans.

## IRTAM: assimilative IRI coefficients

Registered as **`irtam`**.  Module `src/tss/data/irtam.py`.

**What it provides.**  The IRI-based Real-Time Assimilative Model assimilates
the GIRO network into IRI every 15 minutes and publishes coefficients for foF2,
hmF2, B0 and B1 that replace IRI's climatological maps for that instant.  That
is event-driven L0 truth: the ionosphere of a real day rather than a monthly
median.

**Where.**  Two routes to the same files.  The static archive
`https://ulcar.uml.edu/GAMBIT/GambitCoefficients/COEFFS_<YYYY>/<YYYY_MM_DD>/`
holds 2000 through 2025, with B0 and B1 named `B0in` and `B1in`.  The servlet
`https://lgdc.uml.edu/rix/gambit-coeffs?charName=<char>&time=<ISO>` computes the
file on request and serves all four characteristics.  `fetch_coefficients`
tries the archive first and falls back to the servlet on a 404 only.  On
2026-09-23 the archive ended at 2025, the servlet served every 2026 epoch tried
up to 2026-09-21 and answered HTTP 500 for 2026-09-22 and for "now".  The
default `tss fetch irtam` therefore takes 12:00 UT three days ago, not
yesterday.  The servlet also answered 429 after about a dozen requests in a few
minutes, so a campaign should take the archive wherever it can and pace servlet
requests at no less than the five-second minimum.

**Coefficient layout.**  1064 numbers, declared as 14 temporal by 76 spatial
Jones and Gallet terms.  The layout follows PyIRTAM (Forsythe et al. 2024,
Space Weather, doi:10.1029/2024SW003965): the first 988 are 13 by 76 in Fortran
order, and the last 76 are the linear-trend term inserted as temporal row 1.
`IrtamCoefficients.matrix` reproduces PyIRTAM 0.0.7's `IRTAM_read_files`
exactly on both fixtures.

**Licence and citation.**  Distributed by the Lowell GIRO Data Center and
derived from GIRO measurements, so the GIRO terms carry over.  Cite Galkin et
al. (2012), Radio Science 47, RS0L07, doi:10.1029/2011RS004952, and Galkin et al.
(2020), Space Weather, doi:10.1029/2020SW002463.

**Cache.**  `data/external/irtam/`.  A past epoch never changes, so its file is
never stale; a servlet file for an epoch within the last day is re-fetched
after an hour, since late data may change it.  Every route's cache is checked
before any network request.

**Rate limits and etiquette.**  Five seconds between requests.  Prefer the
archive, which is static files.  A malformed or truncated body is refused
before it is cached.

**Checked in versus downloaded.**  `tests/fixtures/irtam/`: one archive file
(foF2, 2012-10-17 23:30 UT) and one servlet file (foF2, 2024-03-01 12:00 UT,
IP redacted).

**Consumed by.**  WP1, for event-driven L0.

**Evaluation.**  `tss.l0_ionosphere.irtam_maps` evaluates the coefficients into
foF2, hmF2, B0 and B1 at any points, with modip from IGRF.  It reproduces
PyIRTAM to rounding error; verification and a first validation against GIRO
are in `docs/irtam-validation.md`.

**Known limitations.**  IRTAM gives the peak and the bottomside shape only; an
L0 tier also needs a topside, which is not built yet.  IRTAM assimilates GIRO,
so it is not an independent check on GIRO and must not be scored against GIRO
stations it assimilated; each file's header lists them.

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
