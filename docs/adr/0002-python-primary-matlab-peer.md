# ADR 0002  Python primary, MATLAB peer, bridged by files

Status: accepted
Date: 2025
Supersedes: nothing

## Context

Both languages are already in use.  The link budget exists in Python and in
MATLAB, with 21 unit tests, and the MATLAB version came first.  The IRI batch
that produced the 1.2 million profile band-edge study was generated in MATLAB.
PHaRLAP 4.7.4, which is needed for the ray-traced propagation tier, is a MATLAB
toolkit with no Python binding.

Three options were considered.

Making everything MATLAB would keep PHaRLAP and the existing code in one place.
It would also mean no free continuous integration runner, since the test machine
has no MATLAB licence, a weaker package and test ecosystem, and a codebase
harder for anyone without a licence to run or review.

Making everything Python would mean reimplementing three-dimensional numerical
ray tracing, and regenerating a completed 1.2 million profile study through a
different IRI binding.  The main output of that work would be a second thing to
verify, for no scientific gain.

Keeping both means a language boundary, which has to be managed.

Given that both are staying, the coupling mechanism is the second decision.  The
MATLAB Engine API for Python is the obvious choice and has the property that any
code path crossing into MATLAB is untestable on a runner without MATLAB, which
in practice means untested.  It also leaves no artefact behind, so an
intermediate cannot be checksummed, recorded in a manifest or rerun against.

## Decision

**Python is the primary implementation.**  New work goes in Python unless there
is a specific reason it cannot.  The package is `tss`, and the layer contracts,
registry, configuration, requirements, provenance and pipeline are Python only.

**MATLAB is a peer**, not a convenience wrapper, for tiers where its toolboxes
or PHaRLAP win, and for migrated code that already works and already has tests.
Two things are MATLAB-only and will stay that way: PHaRLAP ray tracing, because
the toolkit is, and IRI batch generation at the scale the band-edge study needs,
because that is the binding the completed result used.

**The two sides exchange data through files**, in `.mat` v7.3, `.npz` and
`.csv`, under `data/bridge/` with a naming convention and a JSON sidecar per
file.  Not through the MATLAB Engine API.  The reasons are hermetic continuous
integration, reproducible intermediate artefacts that can be checksummed and
archived, and the ability to run a long MATLAB generation step independently of
the Python analysis.

**MATLAB functions use a `tss_` prefix in a flat `matlab/` tree**, for example
`tss_link_budget.m`, `tss_iri_batch.m`, `tss_raytrace_vertical.m`.  There is no
`+tss` package namespace.  A package namespace would be the more modern MATLAB
style and would give real scoping.  It is rejected because the author's existing
code is already written with the prefix convention, `tss_link_budget.m` among
it, and adopting `+tss` would mean rewriting every call site in migrated code to
`tss.link_budget(...)`.  That rewrite would touch code covered by passing tests
in order to change nothing about behaviour, which is the kind of change that
introduces defects for no benefit.  Consistency with what exists wins here.

**Where a function exists on both sides, one is the reference and the other is a
transcription**, stated in both docstrings, and the transcription is checked
against the reference to machine precision by the transcription verifier, which
is part of verification level V2.  For the link budget, Python is the reference
and MATLAB is the transcription, even though the MATLAB version came first.

## Consequences

The Python half of every language crossing runs in continuous integration,
because the fixtures are checked in and the tests read them from disk.  MATLAB
is needed only to regenerate a fixture or to run the V2 transcription cases,
and when it is absent those cases report as skipped rather than silently
passing.

Bridge files are durable artefacts.  A ray set that took hours to generate is on
disk, checksummed, recorded in the manifest, and reusable.  The PHaRLAP output
in particular is expensive enough that keeping it is worth the space.

Having two independent implementations of the link budget is a verification
asset rather than duplication, and the transcription check is the level that
most often finds real bugs, because two implementations written at different
times make different mistakes.

The drawback is that a file handoff is slow per call, so it cannot be used
inside a loop.  The granularity has to be one batch, not one profile.  A
MATLAB-side tier is therefore a two-step operation, and a campaign selecting one
is not a single command: the generation step runs first and the Python study
fails with a clear message naming the missing file and its generator if it has
not.

The prefix convention means MATLAB has no real namespace, so `tss_` is the only
thing preventing a collision with a toolbox function.  That is weaker than a
package and it is accepted.  It also means `startup_tss.m` has to add
directories to the path explicitly, where a package would be found by its parent
directory alone.

Keeping two implementations of anything costs maintenance.  Changing a reference
without rerunning the verifier against its transcription is a defect, and the
verifier only runs where MATLAB is available, which is the author's machine.
The mitigation is that the verifier is part of V2 and that a verification report
states loudly which cases it could not run.  The mitigation is not complete, and
the honest position is that the MATLAB side is checked less often than the
Python side.

A side benefit to the file bridge is that the inputs and outputs it writes serve
directly as the V2 transcription fixtures, with no extra work, because they are
already on disk in a documented format.
