# The MATLAB and Python bridge

Python is the primary implementation.  MATLAB is a peer for the tiers where its
toolboxes or PHaRLAP win, and for migrated code that already works and already
has tests.  The two sides exchange data through files on disk.  This document
says how.

The decision itself, and the reasoning behind it, is recorded in
`docs/adr/0002-python-primary-matlab-peer.md`.  This document is the operating
manual.

## Why files rather than the MATLAB Engine API

The MATLAB Engine API for Python is the obvious alternative and it was not
chosen.  Three reasons.

**Hermetic continuous integration.**  The test runner has no MATLAB and no
MATLAB licence.  With the engine API, any code path that crosses into MATLAB is
untestable on the runner, which in practice means it is untested.  With files,
the Python side is tested against checked-in `.mat` and `.npz` fixtures produced
once by MATLAB, so the Python half of every crossing runs in CI.

**Reproducible intermediate artefacts.**  A file handed across the boundary is a
thing that can be checksummed, recorded in a run manifest, archived and rerun
against six months later.  An in-process engine call leaves nothing behind.
Since the reproducibility promise of the dissertation is that every figure is
regenerable from stored configuration, an intermediate that evaporates is a hole
in that promise.  The PHaRLAP ray sets in particular take long enough to
generate that keeping them is worth the disk.

**Decoupled execution.**  A MATLAB-side generation step can run overnight on the
machine that has the licence while the Python analysis runs anywhere.  The
1.2 million profile IRI batch used by the band-edge study is generated this way.
It would be awkward to drive that from a Python process holding an engine
handle for hours.

**The drawback** is real and should be stated.  A file-based handoff is slower
per call, so it cannot be used inside a tight loop; the granularity has to be
one batch, not one profile.  It also makes a MATLAB-side tier a two-step
operation rather than a function call, which means a campaign that selects a
MATLAB tier is not a single command.  A side benefit to this choice is that the
same files serve as the V2 transcription fixtures without any extra work, since
the inputs and outputs are already on disk in a documented format.

## Handoff formats

Three formats, chosen by what is being carried.

**`.mat` v7.3** for anything MATLAB produces that Python reads: ray sets, IRI
batches, link budget reference outputs.  Version 7.3 is required, not
preferred.  It is HDF5 underneath, so `h5py` reads it without `scipy.io`
fallbacks and without the 2 GB limit of the older format.  Save with
`save(path, 'varname', '-v7.3')`.  The Python reader is
`tss.bridge.read_mat` and it imports `h5py` lazily.

A v7.3 file transposes array dimensions relative to MATLAB's own view, because
HDF5 is row-major and MATLAB is column-major.  `read_mat` transposes back, so
Python sees the array with the shape MATLAB reported.  Do not do that transpose
a second time in the caller.

**`.npz`** for anything Python produces that MATLAB reads, and for
Python-to-Python intermediates.  MATLAB reads `.npz` through the small reader in
`matlab/tss_read_npz.m`, which handles the subset of the format the repository
uses: numeric arrays, real and complex, no object arrays and no pickled
objects.  Always save with `numpy.savez` and explicit keyword names, never
positional, because the MATLAB reader looks names up.

**`.csv`** for small tabular results that a human will look at: requirement
scoring tables, comparison summaries, the per-profile outputs of the band-edge
study.  CSV is used where the number of rows is modest and where being able to
open the file in a spreadsheet is worth more than being able to store an array.
Write a header row.  Use the SI column names with their unit suffixes, the same
names the code uses, so that a column in a CSV and a field in a dataclass are
obviously the same quantity.

Nothing is handed across in `.mat` v7 or earlier, in `.txt`, or in a pickle.
Pickle in particular is excluded because it is neither readable by MATLAB nor
safe to load from an archive.

## Naming convention

```
data/bridge/<producer>/<study>_<what>_<tag>.<ext>
```

* `<producer>` is `matlab` or `python`, naming the side that wrote the file.
* `<study>` is the study or work package the file belongs to, for example
  `bandedge`, `linkbudget`, `raytrace_v2`.
* `<what>` names the contents in the same vocabulary the code uses: `profiles`,
  `rays`, `echoes`, `snr`, `reference`.
* `<tag>` is a short discriminator.  Use the configuration name when there is
  one, otherwise the eight-character config hash from the manifest.  Do not use
  a date; the manifest carries the date and a date in the name invites files
  that differ only by when they were made.
* The extension follows the format rules above.

Examples:

```
data/bridge/matlab/bandedge_profiles_iri2020_r12-105.mat
data/bridge/matlab/raytrace_v2_rays_vertical_3f0a91c7.mat
data/bridge/python/linkbudget_reference_baseline_6u.npz
data/bridge/python/bandedge_summary_iri2020_r12-105.csv
```

Every bridge file gets a sidecar `<name>.json` written by the producing side,
carrying the producer, the version of the code that made it, the inputs it
consumed, and the same config hash that appears in the run manifest.  A bridge
file without a sidecar is treated as untrusted and a study that reads one warns.

`data/bridge/` is not checked in.  The sidecars alone are small enough that they
could be, and they are not either, for consistency: the manifest in `results/`
is the checked-in record.

## Path setup: `startup_tss.m`

`matlab/startup_tss.m` is the single entry point on the MATLAB side.  Run it
from the repository root, or from anywhere, since it locates the root from its
own file location rather than from `pwd`:

```matlab
>> run('C:\path\to\tss-sim\matlab\startup_tss.m')
```

It does five things:

1. Finds the repository root from `mfilename('fullpath')` and stores it in a
   `tss_root` function so that no MATLAB function hard-codes a path.
2. Adds `matlab/` and its subdirectories to the path, excluding anything under
   a directory named `private` or starting with a dot.
3. Reads `config/default.yaml` for the `paths:` block and adds the PHaRLAP
   installation to the path from there, rather than from a literal path in the
   script.  On the author's machine that resolves to
   `C:\Users\Andrew\PhD\Chartier Verification\pharlap_4.7.4`.  If the path is
   absent it warns and continues, so the non-PHaRLAP functions still work.
4. Sets the PHaRLAP environment variables the toolkit expects for its reference
   data directory.
5. Creates `data/bridge/matlab/` if it does not exist, and prints a one-line
   summary saying which components were found.

It does not change `pwd`, does not call `clear all`, and does not save the path.
A MATLAB session that has not run it will fail on the first `tss_` call with a
function-not-found error, which is the intended failure: a half-configured path
is worse than an obvious one.

## Which tiers are MATLAB-side, and why

Two things are MATLAB-only, and the reason in both cases is that the tool is.

**PHaRLAP ray tracing.**  PHaRLAP 4.7.4 is a MATLAB toolkit.  There is no Python
binding, and reimplementing three-dimensional numerical ray tracing in Python in
order to avoid a language boundary would be a large piece of work whose main
output would be a second thing to verify.  So the ray-traced L1 tier generates
its ray sets in MATLAB, writes them as `.mat` v7.3, and the Python tier
`tss.l1_propagation.raytrace` reads them back and presents them through the
ordinary `L1Propagation` protocol.  From the pipeline's point of view it is a
tier like any other; the difference is that it requires a generation step to
have been run first, and it says so with a clear error if the ray set it needs
is not in `data/bridge/`.

**IRI batch generation at scale.**  The band-edge study consumes 1.2 million
IRI-2020 profiles.  The MATLAB IRI interface on the author's machine is the one
that has been exercised at that volume, and regenerating the set through a
different binding would invalidate a completed result for no scientific gain.
`matlab/tss_iri_batch.m` writes the profiles as `.mat` v7.3 and the
`matlab_iri` L0 tier reads them.

**Everything else is Python-primary.**  The link budget is the interesting case,
because it exists on both sides.  Python is the reference; `matlab/tss_link_budget.m`
is the transcription, kept because it predates the Python version and because
having an independent implementation is itself a verification asset.  Both
docstrings say which is which.

## The transcription verifier

When a function exists on both sides, one is the reference and the other is a
transcription.  The verifier checks that the transcription agrees with the
reference to machine precision over a shared set of inputs.  It is part of
verification level V2.

The procedure:

1. The Python side generates a set of input cases, chosen to cover the
   parameter space including the edges: band edges, zero and maximum
   absorption, the smallest and largest antenna sizes, and the degenerate cases
   where a term goes to zero.  It writes them to
   `tests/fixtures/bridge/<function>_inputs.npz`.
2. The Python reference is evaluated over those inputs and the outputs are
   written to `tests/fixtures/bridge/<function>_python.npz`.  Both files are
   checked in.  They are small, and checking them in is what lets the Python
   half of the comparison run in CI.
3. `matlab/tss_verify_transcription.m` reads the inputs, evaluates the MATLAB
   transcription, and writes
   `data/bridge/matlab/<function>_matlab.mat`.
4. `tss verify --level V2` reads both output sets and compares them.

The comparison standard is floating point round-off, not plotting accuracy.  The
default tolerances are a relative tolerance of 1e-12 and an absolute tolerance
scaled to the magnitude of the quantity, and a discrepancy larger than that is a
finding rather than noise.  Quantities that legitimately cannot meet that, an
iterative solve with its own convergence tolerance being the usual case, declare
their own tolerance in the verifier's case table together with a one-line
justification.  There are currently two such declarations and both are in the
absorption integration.

Three rules go with it:

* Changing a reference without rerunning the verifier against its transcription
  is a defect.
* A discrepancy is not fixed by loosening the tolerance.  Find out which side is
  wrong, fix that side, and record the finding in the commit body.  A
  transcription that disagrees with its reference has usually found a real bug
  in one of them, and that is most of the value of having two implementations.
* When MATLAB is absent, `tss verify --level V2` skips the transcription cases
  and says so loudly in its output and in the written report.  A verification
  report that silently omits the cases it could not run is worse than no report.

## Day-to-day

```matlab
% MATLAB side: generate and hand across
run('matlab/startup_tss.m')
tss_iri_batch('config/campaigns/bandedge.yaml')          % writes .mat v7.3
tss_raytrace_vertical('config/campaigns/v2_raytrace.yaml')
tss_verify_transcription('tss_link_budget')
```

```bash
# Python side: consume, run, verify
tss config bandedge                  # confirm the campaign resolves
tss study band_edge --strict         # reads data/bridge/matlab/, writes results/
tss verify --level V2 --write
```

If a Python study reports that a bridge file is missing, the fix is to run the
corresponding MATLAB generator, not to loosen the study.  The error message
names the expected path and the generator that produces it.
