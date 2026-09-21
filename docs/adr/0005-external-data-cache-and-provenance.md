# ADR 0005  External data cache and provenance

Status: accepted
Date: 2025
Supersedes: nothing

## Context

The repository consumes four external data sources and three local model
installations.  The NASA SPDF Alouette and ISIS archive supplies the 5,118
hand-scaled profiles that verification level V4 rests on.  wspr.live supplies
the spot data that WP12 will use to constrain a sub-megahertz noise reference.
The data.gov absorption dataset is an independent check on the WP2 absorption
model.  IRI-2020, NeQuick 2 and PHaRLAP are installed locally and licensed
separately.

Four problems follow from that.

External sources change.  The SPDF holdings have grown as further films were
digitised, data.gov catalogue entries are revised, and wspr.live accumulates
rows continuously.  A result that says "validated against the ISIS archive"
without saying which bytes is not reproducible, and re-running it a year later
against a larger archive gives a different answer for reasons that have nothing
to do with the model.

Some of these data sets are large.  Checking them into git would make the
repository unusable and would not solve the provenance problem anyway, because
the reader would still have to trust that the checked-in copy was what the
server sent.

Public archives are a shared resource.  A test suite that fetches from SPDF or
wspr.live on every commit is impolite and will eventually be blocked, and the
wspr.live ClickHouse endpoint in particular can be asked for hundreds of
millions of rows by accident.

And a fetcher is code that needs testing, but it cannot be tested against a live
endpoint without reintroducing the first three problems.

## Decision

**A cache under `data/external/`**, one subdirectory per source, with raw
downloads kept exactly as received under `raw/` and never edited, and derived or
parsed forms in sibling directories regenerated from `raw/`.

**An `index.json` per source**, written by the fetcher, recording the URL or the
exact query text, the fetch timestamp in UTC, and the byte count and SHA-256 of
every file.  For wspr.live the query SQL is stored alongside its result, so a
result can be traced to the query that produced it.

**Checksums reach the run manifest.**  A study calls
`Manifest.add_input(path, role)` for every external file it reads, and
`src/tss/provenance.py` records the path, role, size and SHA-256.  A result is
reproducible because the manifest names the bytes, not because the bytes are in
the repository.

**Nothing under `data/external/` is checked in.**  It is in `.gitignore`.  What
is checked in is the fetcher module under `src/tss/data/`, a small
representative extract under `tests/fixtures/<source>/`, and the manifest
entries recording what a run consumed.

**Tests never hit the network.**  Not once, not with a skip guard, not only in
continuous integration.  Fetchers are tested against recorded responses and
parsers against the checked-in extracts.  A test needing a cold cache carries
`@pytest.mark.data` and is skipped, visibly, when the cache is absent.

**Fetch etiquette is part of the fetcher, not of the operator.**  Serial
requests with a delay, a descriptive `User-Agent` naming the project and a
contact address, conditional requests against the cached index where the server
supports them, time-bounded and server-side-aggregated queries for wspr.live,
and chunked collection so an interrupted fetch resumes rather than restarts.
`tss fetch <source>` does nothing when the cache is complete and checksums
match, unless `--force` is given.

**Local model installations are configured, not hard-coded.**  Paths to
IRI-2020, NeQuick 2 and PHaRLAP live in the `paths:` block of
`config/default.yaml`, and no module contains a literal path.  These models are
not redistributed; each user obtains and licenses them.

## Consequences

A result can be traced to the exact bytes it consumed, and a rerun that produces
a different answer can be diagnosed by comparing manifest checksums rather than
by guessing whether the archive changed.

The repository stays small and the sources stay unstressed.  A cold-cache
machine fetches once, politely, and thereafter reads from disk.

Because `raw/` is never edited, a parsing bug is fixed by correcting the parser
and regenerating the derived directory.  No re-download is needed and there is
no doubt about what was downloaded.  That has already paid for itself once with
the ISIS scaled profiles.

The no-network rule means the fetchers themselves are only lightly tested.  The
parsing, chunking and caching logic runs against fixtures in continuous
integration, but the part that actually talks to the server is exercised only
when the author runs a real fetch.  A change in the server's response format
would not be caught by the test suite.  That is a real gap and the mitigation is
thin: the fetchers validate what they receive against the expected schema and
fail loudly rather than writing a malformed cache.

Ignoring `data/external/` means a fresh clone cannot reproduce a V4 result
without first fetching several gigabytes, which takes time and requires the
archive still to be up.  If SPDF were to go offline, the checked-in fixture
extracts would not be enough to rerun the full comparison.  The honest statement
is that reproducibility of the V4 results depends on an external service
remaining available, and that a separate archival copy of the exact files used
should be deposited with the dissertation rather than relied on from the live
archive.

Configuring model paths rather than hard-coding them means a machine without
PHaRLAP still imports the package and runs everything else, with the affected
tests skipped rather than failing.  The cost is that a missing path shows up as
a skip, and a run on a machine missing three of them looks superficially like a
clean run.  The verification report counts and names its skips at the top for
this reason.

A side benefit to storing the wspr.live query text alongside its result is that
the queries themselves become reviewable.  A noise reference built from an
ill-posed query is a common way to get a wrong answer confidently, and having
the SQL in the cache next to the rows makes that checkable by someone other than
the author.
