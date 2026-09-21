# Verification and validation

Verification runs at four levels.  Each answers a different question, each has
reach that the others do not, and each has limits that have to be stated because
a verification claim without its limits is not a verification claim.

| Level | Question it answers | Where it lives |
| --- | --- | --- |
| V1 | Is the code internally consistent and does it behave correctly in cases where the answer is known in closed form? | `tests/`, driven by `pytest` and by `tss verify --level V1` |
| V2 | Do two independent implementations of the same physics agree? | `src/tss/verify/v2_crossfidelity.py`, `matlab/tss_verify_transcription.m` |
| V3 | Does the model reproduce published results? | `src/tss/verify/v3_literature.py`, `tests/verify/` |
| V4 | Does the model reproduce archived real measurements? | `src/tss/verify/v4_archive.py`, `src/tss/studies/isis_validation.py` |

`tss verify --level all --write` runs everything available on the machine and
writes a report under `results/verification/`, with a manifest.  Cases that
could not run, because MATLAB or PHaRLAP or a cached data set is absent, appear
in the report as skipped with the reason.  They are never silently omitted.

## V1  Internal consistency

**What it covers.**

*Unit tests.*  The ordinary kind, per module, colocated by layer under
`tests/`.  These check the arithmetic, the handling of edge inputs, the shapes
and dtypes of returned arrays, and the behaviour at the boundaries of the valid
domain.

*Contract tests.*  Every registered tier is run through the shared tests in
`tests/test_interfaces.py`, which check the four contract rules from
`src/tss/interfaces.py`: inputs are not mutated, the declared object type comes
back, `describe()` round-trips through JSON, and two calls with the same seeded
generator give bitwise identical output.  A new tier gets these for free by
being registered.

*Closed-form limiting cases.*  Where the physics has an analytic answer, the
code is checked against it.  Vertical incidence in a parabolic layer gives a
known virtual height.  Absorption with collision frequency set to zero gives
zero absorption.  A Chapman layer at its peak gives a known plasma frequency.
Free space with the ionosphere removed gives the geometric range.  The
de-chirped beat frequency for a single stationary target gives the known
apparent range exactly, with no window and no noise.  These cases are where a
sign error or a factor of two is actually caught; a test that only checks that
the output has the right shape catches almost nothing.

*Convergence.*  Numerical results are checked to converge as the grid refines,
the sample rate rises or the integration step shrinks, and the operating value
is chosen with margin over the point where the answer stops moving.
Convergence tests carry `@pytest.mark.slow` and are excluded from the fast
suite.

*Conservation and consistency between layers.*  The SNR that L3 measures for a
given configuration must match what the L2 link budget predicts for the same
configuration, within the loss the processing is known to introduce.  A
mismatch means one of the two is wrong, and this comparison is one of the more
productive V1 tests because it couples a closed-form calculation to a
sample-level simulation.

**What it cannot reach.**  V1 checks that the code does what the author
intended.  It cannot check that the intention was right.  A model built on a
wrong physical assumption passes V1 completely, because every test was written
by the same person who made the assumption.  That is the whole reason V2 to V4
exist.

**Where the tests live.**  `tests/`, mirroring the layer packages, plus
`tests/test_interfaces.py` and `tests/test_metadata.py`.  Run with `pytest -q`,
or `pytest -q -m "not slow"` for the fast suite.

## V2  Cross-fidelity

**What it covers.**

*Analytic against ray-traced.*  The analytic L1 tier and the PHaRLAP ray-traced
tier are run over the same L0 truth at vertical incidence, where the two should
agree, and the group delays and amplitudes are compared.  Agreement establishes
that the analytic tier is usable for the large sweeps where the ray tracer is
too slow.  Disagreement tells you which approximation in the analytic tier is
doing the damage.

*Independent re-implementation of ported code.*  Where a function exists in both
Python and MATLAB, the transcription is checked against its reference to machine
precision over a shared input set.  The procedure and the tolerance standard are
in `docs/matlab-python-bridge.md`.  This is the level that most often finds real
bugs, because two implementations written at different times by the same author
make different mistakes.

*Model against model at L0.*  IRI-2020 against NeQuick 2 over the same epochs
and locations.  This is not a correctness test, since neither is truth.  It is a
sensitivity test: a result that changes materially between the two topside
formulations is a result that depends on the topside model rather than on the
instrument, and it has to be reported that way.  This case is blocked on the
NeQuick 2 tier, which is the outstanding half of WP1.

*Fidelity ladder at L2.*  The closed-form link budget against the sample-level
receiver model, which is the same comparison described under V1 seen from the
other side.

**What it cannot reach.**  Two implementations that share an assumption agree
with each other and are both wrong together.  The analytic and ray-traced tiers
both use the same refractive index formulation; if that formulation is wrong for
the topside, V2 will not notice.  Neither will a transcription check, which by
construction compares an implementation against a copy of itself.

**Where it lives.**  `src/tss/verify/v2_crossfidelity.py` and
`matlab/tss_verify_transcription.m`, driven by `tss verify --level V2`.  The
cheap cases also run under `pytest` in `tests/verify/`.  Cases needing MATLAB or
PHaRLAP are skipped, loudly, when those are absent.

## V3  Literature replication

**What it covers.**  The model is set up to reproduce a published result, and
the output is compared against what was published.  This is the level that tests
assumptions rather than arithmetic, because the published result was produced by
someone who did not share the author's assumptions.

The replication set:

* **Barry 1971**, as a reference for the HF chirp sounder technique and the
  de-chirp signal processing, and for the range resolution relationships that
  INS-05 rests on.
* **Poole 1985**, for chirp sounder ionogram formation and the treatment of the
  sweep in relation to the ionospheric trace.
* **Zawdie 2017**, for topside sounding and modelled topside ionograms, which is
  the closest published analogue to what L1 through L3 produce.
* **Franklin and MacLean 1969**, for the external noise reference the link
  budget uses.  This is the reference that runs out below about 0.6 MHz, which
  is the origin of the extrapolation that WP12 exists to remove.
* **Chartier 2022**, for the PHaRLAP-based propagation comparison.  The PHaRLAP
  installation used here comes from the Chartier verification work on the
  author's machine, which means the ray-tracing configuration can be checked
  against a known-good setup rather than assembled from the manual.

Each replication is a study under `src/tss/studies/` with a stored
configuration, so it reruns on demand rather than being a one-time check.  Each
carries a note recording what agreed, what did not, and by how much.  A
replication that does not fully agree is still useful and is reported as a
partial agreement with the discrepancy quantified.  Removing a replication
because it disagrees would be the wrong response.

**What it cannot reach.**  Published results carry their own errors and their
own undocumented assumptions.  Digitising a figure from a 1971 paper introduces
a reading error of its own, and several of these replications compare against
digitised figures rather than against tabulated numbers, so agreement to a few
per cent is the practical ceiling.  V3 also cannot cover anything nobody has
published, which is most of what is novel here: there is no published
sub-watt CubeSat topside FMCW sounder to replicate.

## V4  Archived real data

**What it covers.**  Digitised Alouette and ISIS ionograms from the NASA SPDF
archive, described in `docs/data-sources.md`.  The delivered case pairs 5,118
hand-scaled ISIS-1 profiles against the model.

The comparison runs the L0 truth model for the time and place of each archived
sounding, propagates it through L1, forms the trace, inverts it at L4, and
compares the retrieved profile against the hand-scaled profile from the archive.
It exercises the ionospheric modelling, the propagation physics, the virtual
height to true height relationship and the inversion, against measurements made
by a real topside sounder in a real ionosphere.

This is the strongest evidence available for the ionospheric half of the chain,
and the ISIS-1 comparison is a substantial part of what makes C1 defensible
above 0.6 MHz.

**What it cannot reach, and this is the important limitation.**

Alouette and ISIS were pulse sounders.  They transmitted a short pulse and
switched the receiver on afterwards.  Consequently:

* **The archive contains no range-Doppler coupling.**  In a pulse system range
  and Doppler are separated by the measurement geometry itself, and a Doppler
  shift does not appear as a range error.  In an FMCW system they are coupled
  through the sweep rate, and separating them with one transmit and receive pair
  is exactly claim C3.  There is nothing in the archive that bears on it.
* **The archive contains no transmit-to-receive leakage.**  A pulse sounder with
  a switched receiver does not have the transmitter running while it listens.
  The continuous-wave leakage path, the isolation requirement INS-23 and the
  de-chirped residual floor INS-24 have no analogue in the archived data.  Claim
  C4 cannot be tested against it at all.
* Sweep nonlinearity, the coherent processing interval of INS-15 and the
  oscillator stability of INS-26 are likewise FMCW-specific and absent from the
  archive.

So V4 validates the ionospheric and retrieval half of the chain, L0, L1, L4 and
the parts of L5 that score profile accuracy.  It does not and cannot validate
the FMCW half.

**How the FMCW half is validated instead.**  By V1 through V3, and by the bench
prototype.  V1 covers the de-chirp arithmetic and the bin-to-range mapping
against closed-form cases.  V2 covers the waveform and receiver implementations
against each other across the language boundary.  V3 covers the chirp technique
against Barry 1971 and Poole 1985.  None of those puts a real transmitter next
to a real receiver, so the leakage and residual numbers that C4 turns on are
ultimately measured on the bench, which is WP9.  The honest statement is that
C4 is a hardware claim with a modelling component, not a modelling claim, and
the dissertation should present it that way.

**Where it lives.**  `src/tss/verify/v4_archive.py` and
`src/tss/studies/isis_validation.py`, with the archive parser in
`src/tss/data/spdf_isis.py`.  The test-suite portion runs against
`tests/fixtures/spdf_isis/` and never touches the network.  The full 5,118-profile
comparison needs the cached archive and carries `@pytest.mark.data`.

## Reading a verification report

The report written by `tss verify --write` has one section per level and, within
each, one row per case with the verdict, the quantity compared, the observed
discrepancy and the tolerance.  Four verdicts:

* **pass**, the case ran and agreed within tolerance.
* **fail**, the case ran and did not agree.  A fail is never suppressed and
  never resolved by widening a tolerance without a recorded justification.
* **skipped**, the case could not run and the reason is given, usually a missing
  MATLAB, PHaRLAP or data cache.  A report with many skips is a weak report and
  says so at the top.
* **partial**, used at V3 only, where a replication agrees in part.  The
  discrepancy is quantified and the note explains it.

The report carries the requirement set revision and the git commit, and it is
written with a manifest like any other result.  A verification report from a
dirty working tree records that it was, which is why publication runs use
`--strict`.

## Current coverage

| Layer | V1 | V2 | V3 | V4 |
| --- | --- | --- | --- | --- |
| L0 ionosphere | Delivered | Blocked on the NeQuick 2 tier, WP1 | Partial | Delivered through the ISIS-1 comparison |
| L1 propagation | Delivered for absorption | Blocked on the PHaRLAP tier, WP2 | Barry 1971 and Chartier 2022 outstanding | Delivered through the ISIS-1 comparison |
| L2 instrument | Delivered for the link budget, 21 tests | Delivered for the link budget transcription | Franklin and MacLean 1969 delivered above 0.6 MHz | Not reachable, see above |
| L3 signal processing | Not started, WP4 | Not started | Barry 1971 and Poole 1985 outstanding | Not reachable, see above |
| L4 retrieval | Not started, WP6 | Not started | Zawdie 2017 outstanding | Delivered through the ISIS-1 comparison |
| L5 assessment | Partial | Not applicable | Not applicable | Partial |

The gaps line up with the work package sequence in `docs/work-packages.md`.  The
largest single gap is that L3 has no verification at any level, because L3 does
not exist yet, and that is WP4.
