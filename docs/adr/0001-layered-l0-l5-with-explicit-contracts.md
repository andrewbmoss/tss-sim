# ADR 0001  Layered L0 to L5 architecture with explicit contracts

Status: accepted
Date: 2025
Supersedes: nothing

## Context

The simulator has to support an error budget per layer, not just an end-to-end
profile error.  The dissertation argues that a sub-watt FMCW topside sounder is
feasible, and an examiner is entitled to ask which part of the chain contributes
which part of the error.  A single end-to-end model produces a final number and
cannot answer that.

There is a second pressure.  Each part of the chain needs more than one fidelity
of implementation.  Propagation needs a fast analytic tier for large campaign
sweeps and a slow ray-traced tier for geometry the analytic tier cannot reach.
The ionosphere needs a closed-form Chapman layer for limiting-case tests,
an empirical model for realism, and stored profiles for comparison against
archived measurements.  The instrument needs a closed-form link budget and a
sample-level receiver model.  Without a structure, those alternatives turn into
flags and branches inside one large function, and the combination that produced
a given figure becomes hard to state.

A third pressure is the schedule.  Several parts of the chain do not exist yet
and will be written over the next two years, some by the author and possibly
some by collaborators.  Whatever structure is chosen has to let an unwritten
layer be represented by a stub without blocking work on the layers around it.

## Decision

Six layers, L0 to L5, with explicit contracts.

* **L0** ionosphere state, **L1** propagation, **L2** instrument, **L3** signal
  processing, **L4** retrieval, **L5** assessment.
* The objects that pass between layers are frozen dataclasses defined in
  `src/tss/interfaces.py`.  That module contains no physics.
* Each layer defines a `Protocol` that its implementations satisfy.  An
  implementation is called a tier.
* A tier registers itself under a layer and a name with `tss.registry`, and a
  campaign selects tiers by name.  The pipeline in `src/tss/pipeline.py` owns
  sequencing and nothing else.
* Four contract rules apply to every tier: inputs are never mutated; units are
  SI unless the attribute name carries a suffix; a tier that cannot answer
  returns a masked result with a reason in `flags` rather than raising, because
  physical impossibility is a result and not an error; randomness is drawn only
  from an injected `numpy.random.Generator`.
* Every tier implements `describe()`, returning a JSON-serialisable record that
  goes into the run manifest.

## Consequences

The error at L5 can be attributed to a layer, because truth at L0 is exact and
every intermediate product is retained when `keep_intermediates` is set.  That
is the property the whole design exists for.

A new fidelity tier is a new module and a decorator.  It does not touch the
pipeline, and a campaign can select it by name.  A tier that is registered but
unimplemented raises when it is called rather than when it is built, so a
campaign naming future work still validates.

Reproducibility follows from the RNG rule.  One generator is created from the
campaign seed and split per configuration, so two configurations do not share a
noise realisation and a single configuration reruns identically from its
recorded seed.

The masked-result rule means a frequency sweep can run straight through the band
edges.  Below the local plasma frequency at the spacecraft the wave is
evanescent and there is no measurement to report; the run continues and the row
is flagged.  The drawback is that a caller who ignores `flags` will treat a
masked result as a measured one, so every consumer of an `EchoSet` or an
`Ionogram` has to check.  That is a real hazard and it is accepted because the
alternative, raising, would make a band-edge study impossible to write.

The costs are interface discipline and boilerplate.  Six contracts is more code
than three loose functions, and the contracts constrain what a tier may do.
Adding a field to a data object touches every tier that produces it.  A tier
that wants to pass something the contract does not carry has either to extend
the contract or to abuse `flags` or `extra`, and the second is tempting.

There is one known violation of the layering, recorded as a seam in
`docs/architecture.md`: the migrated link budget module is kept whole under L2
and the L0 and L1 packages re-export pieces of it.  The layering is therefore
a description of most of the code rather than all of it, and that has to be
stated wherever the architecture is presented.

A side benefit to this choice is that the verification plan falls out of it.
V1 compares a layer against the one beneath it, and V2 compares two tiers at the
same layer.  Neither comparison would be expressible without the layer boundary
being a real thing in the code.
