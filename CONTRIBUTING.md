# Contributing to tss-sim

This is a single-author research codebase supporting a dissertation, so the
conventions here are aimed at keeping results traceable and reviewable rather
than at coordinating a large team.  They still apply to every change, including
changes made by the author, because the reason for them is reproducibility and
not politeness.

## Units and naming

SI everywhere, with a suffix on the name whenever a quantity is carried in
anything else.  A bare `alt` is ambiguous, so the codebase uses `alt_m` and
`alt_km` and never `alt`.  The suffixes in use are listed in `src/tss/units.py`:

```
_m   metres        _km   kilometres     _s    seconds
_hz  hertz         _mhz  megahertz      _mps  metres per second
_k   kelvin        _w    watts          _dbm  dBm
_db  decibels      _deg  degrees        _rad  radians
```

Write conversions using the helpers in `tss.units` rather than inline factors.
`db` and `undb` are power ratios; `db_amplitude` is the 20 log10 form.  The link
budget uses power ratios throughout, so if you find yourself writing a bare
`20 * np.log10` somewhere in L2, check that you meant it.

Array attribute names are plural or unsuffixed only when the interface in
`src/tss/interfaces.py` says so.  Do not invent a parallel spelling for a field
that already exists on a dataclass there.

## One tier per module

Each fidelity implementation lives in its own module under its layer package.
`src/tss/l0_ionosphere/chapman.py` holds the alpha-Chapman tier and nothing
else.  A module that holds two tiers is acceptable only when they share a
substantial private core, and in that case say so in the module docstring.

The layer packages themselves are thin.  A layer `__init__.py` imports its tier
modules so that registration happens, and re-exports the names another layer
legitimately needs.  It should not contain physics.

There is one standing exception, the migrated link budget, documented as a known
seam in `docs/architecture.md`.  Do not add new exceptions without recording the
reason in that document.

## Register via decorator

A tier makes itself available by decorating its class:

```python
from tss.registry import register

@register("l2", "baseline")
class BaselineInstrument:
    ...
```

The decorator sets `name` and `layer` on the class.  Do not assign those
attributes yourself, and do not add entries to the registry dictionaries
directly.  A duplicate name for the same layer raises at import time, which is
what you want: a silently shadowed tier would produce results attributed to the
wrong implementation.

Every tier must implement `describe()` returning a JSON-serialisable mapping of
its identity and settings.  The manifest records that mapping verbatim, so if a
setting changes a result, it belongs in `describe()`.

A registered tier whose implementation is still a stub raises
`NotImplementedError` when it is called, not when it is built.  That is
deliberate, because it lets `tss config <campaign>` validate a campaign that
names a tier before that tier is finished.

## No heavy imports at module scope

`import tss` imports every layer package in order to populate the registry.  It
must stay cheap and must work on a machine that has none of the optional
dependencies installed.  PHaRLAP, the MATLAB engine, the IRI and NeQuick
bindings, `h5py`, `requests` and `matplotlib` are all imported inside the
function or method that needs them.

```python
def echoes(self, state, platform, freq_hz, rng):
    from tss.l1_propagation._pharlap import trace_rays   # lazy, optional
    ...
```

`numpy` and `pyyaml` are the exceptions.  They are hard dependencies and may be
imported at module scope.

A test that guards this lives in `tests/test_import_cost.py`.  If you add a
module-scope import of an optional package, that test fails.

## Docstring metadata block

Every module under `src/tss` carries a metadata block at the end of its module
docstring.  It is how a reader gets from a file to the work package and the
claim it serves without going through the proposal.  The fields are required
even when the answer is "none".

```python
"""Alpha-Chapman ionosphere with a parametrised topside scale height.

A closed-form profile used for V1 limiting-case tests and for the fast tier of
the trade study, where an IRI call per grid point would be too slow.

Layer:        L0
WP:           WP1
Claims:       C1, C2
Status:       delivered
Requirements: INS-01, INS-02
"""
```

* `Layer` is one of `L0` to `L5`, or `support` for modules such as
  `units.py`, `config.py`, `provenance.py` and `registry.py`.
* `WP` is one or more of `WP0` to `WP12`.
* `Claims` lists the claims the module contributes evidence to, or `none`.
* `Status` is one of `stub`, `in progress`, `delivered`, `migrated`.  Use
  `migrated` for code brought in from the author's earlier work whose tests
  came with it.
* `Requirements` lists the requirement IDs the module produces an achieved
  value for, or `none`.  These must exist in
  `config/requirements/requirements.yaml`.

`tests/test_metadata.py` parses these blocks and fails on a missing field, an
unknown work package, or a requirement ID that is not in the requirement set.
The check is cheap and it catches the common drift, which is a requirement being
renumbered in the document and not in the code.

## Tests colocated by layer

`tests/` mirrors the layer packages:

```
tests/
  conftest.py              shared fixtures, including the seeded RNG
  fixtures/                small checked-in data, used instead of the network
  test_interfaces.py       contract tests every tier is run through
  test_metadata.py         docstring metadata block check
  l0_ionosphere/
  l1_propagation/
  l2_instrument/
  l3_signal/
  l4_retrieval/
  l5_assessment/
  verify/                  V2 to V4 drivers that are cheap enough to run in CI
```

A new tier gets tests in the directory for its layer.  At minimum it is run
through the shared contract tests in `test_interfaces.py`, which check that the
tier does not mutate its inputs, that it returns the declared object type, that
`describe()` round-trips through JSON, and that two calls with the same seeded
generator give identical output.

Markers:

* `@pytest.mark.slow` for anything over a couple of seconds.  Convergence
  studies and large IRI sweeps carry it.
* `@pytest.mark.matlab` for a test that needs MATLAB on the path.  These are
  skipped when MATLAB is absent, which is the normal case in CI.
* `@pytest.mark.pharlap` likewise for PHaRLAP.
* `@pytest.mark.data` for a test that needs a cached external data set under
  `data/external/`.  Skipped when the cache is cold.

No test may reach the network.  Use the fixtures under `tests/fixtures/`.  If a
new external source is added, check in a small representative extract alongside
the fetcher and document both in `docs/data-sources.md`.

Running them:

```bash
pytest -q                       # everything available on this machine
pytest -q -m "not slow"         # the fast suite, which is what CI runs
pytest -q tests/l2_instrument   # one layer
pytest -q --cov=tss --cov-report=term-missing
tss verify --level V1           # the verification drivers rather than pytest
```

The fast suite is expected to stay under about a minute.  If a change pushes it
past that, mark the offending test slow rather than letting the suite rot.

## Commit messages

One logical change per commit.  Subject line in the imperative, under 72
characters, prefixed with the scope:

```
l2: add transmit-to-receive leakage model with configurable isolation

Implements the coupling path between the transmit and receive antennas as a
frequency-dependent isolation profile, so that INS-23 can be scored rather than
assumed.  The de-chirped residual is not yet modelled; that is the remaining
piece of WP3 and is tracked in docs/work-packages.md.

Refs: WP3, C4, INS-23, INS-24
```

Scopes are the layer package names (`l0` to `l5`), plus `config`, `req`,
`data`, `matlab`, `docs`, `tests`, `ci` and `scripts`.

End the body with a `Refs:` line naming the work packages, claims and
requirement IDs the change touches.  It is the cheapest form of traceability
available and it makes the dissertation write-up much easier, because the
history can be queried by claim.

Commits that change a number appearing in the dissertation must say so in the
body, and must say which figure or table.

## Keeping MATLAB and Python in step

The rules are in `docs/matlab-python-bridge.md`.  The short form:

1. Python is the primary implementation.  MATLAB is a peer for tiers where its
   toolboxes or PHaRLAP win, and for migrated code that already works.
2. MATLAB functions are named `tss_something.m` in a flat `matlab/` tree.  There
   is no `+tss` package namespace.  The reason is recorded in
   `docs/adr/0002-python-primary-matlab-peer.md`.
3. When a function exists on both sides, one of them is the reference and the
   other is a transcription.  Say which in both docstrings.
4. A transcription is checked against its reference by the transcription
   verifier, which runs both over a shared set of inputs written to
   `tests/fixtures/bridge/` and compares to machine precision.  Agreement to
   floating point round-off is the standard, not agreement to plotting
   accuracy.
5. Changing a reference without rerunning the verifier against its transcription
   is a defect.  The verifier is part of V2 and `tss verify --level V2` runs it
   when MATLAB is available.  It is skipped, loudly, when MATLAB is absent.
6. Do not fix a discrepancy by loosening the tolerance.  Find out which side is
   wrong and record the finding in the commit body.

## Documentation that must move with the code

* A new tier: add it to the tier lists in `docs/architecture.md`.
* A new external source: add a section to `docs/data-sources.md` and a fetcher
  under `src/tss/data/`.
* A change in requirement thresholds: change
  `config/requirements/requirements.yaml`, and follow the procedure in
  `docs/requirements-traceability.md` so the YAML, the Rev E document and the
  SRM spreadsheet do not drift apart.
* A decision with consequences that outlive the change: write an ADR under
  `docs/adr/`, numbered sequentially.  ADRs are short.  Context, Decision,
  Consequences, one page.
* A change in work package status: update `docs/work-packages.md` and the
  summary table in `README.md` together.  Those two tables disagreeing is the
  most common documentation defect in this repository.
