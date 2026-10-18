# Architecture

This is the authoritative interface document for `tss-sim`.  It describes what
each layer is responsible for, what object goes in, what object comes out, which
tier implementations exist or are planned, the fidelity tiers available at each
layer, and the seams that are known to be imperfect.

The code that this document describes is `src/tss/interfaces.py`.  Where the two
disagree, the code is right and this document is stale; say so in the commit
that fixes it.

## Why layers

The simulator generates a known ionosphere, puts it through the instrument, and
then tries to recover it.  Because the truth at L0 is exact, every discrepancy
at L5 can be attributed to the layer that introduced it.  That attribution is
the point.  A monolithic end-to-end model would produce the same final numbers
and would not tell you whether an error came from the absorption model, the
receiver noise figure, the trace extractor or the inversion.

The cost is interface discipline.  Six layers with stable contracts is more code
than three layers without them, and the contracts constrain what a tier is
allowed to do.  That cost is accepted because the dissertation has to defend
error budgets per layer, not just a final profile error.

## Contract rules

These apply to every tier implementation at every layer.  They are stated in the
module docstring of `src/tss/interfaces.py` and are enforced by the shared
contract tests in `tests/test_interfaces.py`.

**1.  Inputs are never mutated.**  A tier returns new objects.  The data objects
are frozen dataclasses, which stops accidental field assignment but does not
stop in-place modification of the numpy arrays they hold.  The contract test
checks array contents before and after a call.  If a tier needs to modify an
input array, it copies first.

**2.  Units are SI unless the attribute name carries a suffix.**  `_km`, `_mhz`,
`_db`, `_dbm`, `_deg`, `_mps` and the rest are listed in `src/tss/units.py`.
Fields without a suffix are SI: `ne` is m^-3, `nu` is s^-1, `freq_hz` is hertz,
`group_delay_s` is seconds, `b_field_t` is tesla.

**3.  Physical impossibility is a result, not an error.**  A tier that cannot
answer for a given input returns a result with the relevant entries masked and a
reason recorded in `flags`, rather than raising.  The band edges are the
motivating case.  Below the local plasma frequency at the spacecraft the wave is
evanescent and there is no measurement to report; that is a fact about the
ionosphere, not a bug, and a campaign that sweeps frequency across the band must
be able to run straight through it.  `EchoSet.flags` carries `"evanescent"` and
`"penetrating"` per row for exactly this.  Exceptions are reserved for
programming errors and for configurations that cannot be built.

**4.  Randomness comes from an injected generator.**  A tier draws only from the
`numpy.random.Generator` passed to it.  No tier calls `numpy.random` at module
level, seeds anything itself, or reads the clock.  The pipeline creates one
generator from the campaign seed and splits it per configuration, so two
configurations in the same campaign do not share a noise realisation and a
single configuration reruns identically from its recorded seed.

**5.  Every tier describes itself.**  `describe()` returns a JSON-serialisable
mapping that goes verbatim into the run manifest.  Any setting that can change a
result belongs in it.

## Data objects

All are frozen dataclasses in `src/tss/interfaces.py`.

| Object | Produced by | Carries |
| --- | --- | --- |
| `IonosphereState` | L0 | `alt_m`, `ne`, `nu`, `te`, `b_field_t`, location, epoch, `source`, `extra` |
| `PlatformState` | campaign / orbit | `alt_m`, `lat_deg`, `lon_deg`, `velocity_ned_mps` (local north-east-down), `attitude_quat`, `epoch` |
| `EchoSet` | L1 | parallel arrays of `freq_hz`, `mode`, `group_delay_s`, `amplitude_db`, `absorption_db`, optional `doppler_hz`, `aoa_deg`, `reflection_alt_m`, `flags` |
| `SweepSchedule` | campaign / L2 | `f_start_hz`, `f_stop_hz`, `k_f_hz_per_s`, `sweep_period_s`, `direction`, `nonlinearity_ppm` |
| `ReceivedSignal` | L2 | `samples` as `(n_channels, n_samples)` complex, `fs_hz`, the `sweep` actually transmitted, `noise_floor_dbm_hz`, `flags` |
| `Ionogram` | L3 | `freq_hz`, `virtual_range_m` and `snr_db` as `(n_freq, n_modes)`, optional `doppler_mps`, `aoa_deg`, `residual_range_bias_m` |
| `DensityProfile` | L4 | `alt_m`, `ne`, optional `ne_sigma`, `alt_sigma_m` |
| `Assessment` | L5 | `configuration`, `requirements` keyed by requirement ID, `errors`, `flags` |

Two details are worth calling out because they are easy to get wrong.

`PlatformState` is a first-class input to L1 and to the L3 novel algorithms, not
a detail buried inside an orbit propagator.  C3 rests on the platform's motion
being large and precisely known, so the motion has to be visible at the
interface where the algorithms that exploit it live.

`ReceivedSignal` carries the sweep that was actually transmitted, not the sweep
that was requested.  L3 de-chirps against what was sent.  This matters because
INS-21's single-reference architecture means the same object is the de-chirp
reference, which is what makes the range-correlation credit claimed in Chapter 5
of the Technical Volume legitimate rather than assumed.

## L0  Ionosphere state

**Responsibility.**  Produce the truth ionosphere at a location and epoch on a
given altitude grid.  Nothing in L0 knows that a radar exists.

**In.**  `lat_deg`, `lon_deg`, `epoch`, `alt_m`.
**Out.**  `IonosphereState`.

A backend that cannot supply electron-neutral collision frequency directly
should derive one from an atmosphere model rather than leave `nu` as `None`,
because L1 absorption needs it.  `b_field_t` is optional in the dataclass but is
needed by the magnetoionic splitting in L1 and by the X-mode cutoff used in the
band-edge study, so in practice any tier used for a scored run supplies it.

`extra` carries scalar diagnostics such as `nmf2`, `hmf2_km` and
`solar_zenith_deg`.  Studies may read them.  No tier may depend on a key being
present.

| Tier | Module | Fidelity | Status |
| --- | --- | --- | --- |
| `chapman` | `src/tss/l0_ionosphere/chapman.py` | Closed-form alpha-Chapman with a parametrised topside scale height.  Fast, smooth, used for V1 limiting cases and the fast tier of the trade study. | Delivered |
| `iri2020` | `src/tss/l0_ionosphere/iri2020.py` | Empirical climatology, the reference truth for most campaigns. | Delivered |
| `tabulated` | `src/tss/l0_ionosphere/tabulated.py` | Reads a stored profile, including the hand-scaled ISIS-1 profiles used by V4. | Delivered |
| `matlab_iri` | `src/tss/l0_ionosphere/matlab_iri.py` | Reads the `.mat` output of the MATLAB IRI batch generator.  Used for the 1.2 million profile band-edge set, which is generated MATLAB-side. | Delivered |
| `nequick2` | `src/tss/l0_ionosphere/nequick.py` | NeQuick 2, giving a second independent topside formulation.  Needed for V2 cross-fidelity at L0. | Planned, WP1 |
| `perturbed` | `src/tss/l0_ionosphere/perturbations.py` | Wraps another L0 tier and superposes travelling disturbances, gradients and depletions. | Planned, WP1 |
| `bubble` | `src/tss/l0_ionosphere/bubbles.py` | Equatorial plasma depletion structures. | Planned, WP10, goal level |

**Fidelity ladder at L0.**  `chapman` for analytic tractability, `iri2020` or
`nequick2` for realism, `perturbed` on top of either for structure, `tabulated`
for real measured profiles.  The ladder is used by V2: an answer that changes
materially between `iri2020` and `nequick2` is an answer that depends on the
topside model rather than on the instrument, and has to be reported that way.

## L1  Propagation

**Responsibility.**  Turn a truth state plus a platform state into the echoes
that return to the platform.  Group delay, amplitude, absorption, magnetoionic
splitting and Doppler.

**In.**  `IonosphereState`, `PlatformState`, `freq_hz`, `rng`.
**Out.**  `EchoSet`.

Rows are one per (frequency, mode, ray).  The arrays are parallel and
one-dimensional so that an analytic tier and a ray-traced tier can return the
same object with different row counts.  A ray-traced tier that finds three rays
at one frequency returns three rows; the analytic tier returns one.

| Tier | Module | Fidelity | Status |
| --- | --- | --- | --- |
| `analytic` | `src/tss/l1_propagation/analytic.py` | Vertical incidence, Appleton-Hartree phase refractive index integrated to the reflection height.  No off-vertical geometry. | Delivered |
| `absorption` | `src/tss/l1_propagation/absorption.py` | Non-deviative and deviative absorption from `ne` and `nu`.  Composed into the analytic tier. | Delivered and tested, WP2 |
| `magnetoionic` | `src/tss/l1_propagation/magnetoionic.py` | O and X mode separation and the X-mode cutoff shift. | Delivered |
| `raytrace` | `src/tss/l1_propagation/raytrace.py` | Three-dimensional numerical ray tracing through PHaRLAP, including off-vertical returns and horizontal gradients.  MATLAB-side, read back through the file bridge. | Planned, WP2 |
| `doppler` | `src/tss/l1_propagation/doppler.py` | Platform-motion Doppler per ray, which is the input C3 has to work from. | Planned, WP2 and WP5 |

**Fidelity ladder at L1.**  `analytic` is fast enough to sweep large campaigns
and is exactly the thing V2 checks against `raytrace`.  The analytic tier cannot
produce off-vertical echoes or ground and oblique returns at all, so any result
that depends on clutter geometry is out of its reach and must wait for WP2.

## L2  Instrument

**Responsibility.**  Everything between the echoes arriving at the spacecraft
and complex baseband samples leaving the receiver.  Sweep generation, antenna
response, front end, noise, oscillator phase noise, transmit-to-receive leakage.

**In.**  `EchoSet`, `SweepSchedule`, `rng`.
**Out.**  `ReceivedSignal`.

This is where the instrument's real impairments live, and it is the layer that
decides whether the concept closes.  The electrically small antennas give poor
radiation efficiency at the bottom of the band, the sub-watt transmitter leaves
little margin, and continuous-wave operation means the transmitter is running
while the receiver is listening.

| Tier | Module | Fidelity | Status |
| --- | --- | --- | --- |
| `linkbudget` | `src/tss/l2_instrument/linkbudget.py` | Power budget from transmitter through antennas, path, absorption and external noise to SNR per range cell.  Migrated with its 21 unit tests. | Delivered, WP3 |
| `baseline` | `src/tss/l2_instrument/receiver.py` | Sweep synthesis, mixing, filtering, quantisation, thermal and external noise, phase noise against INS-26. | In progress, WP3 |
| `waveform` | `src/tss/l2_instrument/waveform.py` | Sweep law generation including the INS-02 frequency step law and sweep nonlinearity. | In progress, WP3 |
| `antenna` | `src/tss/l2_instrument/antenna.py` | Electrically small antenna efficiency, pattern and impedance against frequency. | Delivered as part of the link budget, WP3 |
| `leakage` | `src/tss/l2_instrument/leakage.py` | Transmit-to-receive coupling as a frequency-dependent isolation profile, scoring INS-23, plus the de-chirped residual that scores INS-24. | Not started.  **Critical path**, WP3, C4 |
| `vector` | `src/tss/l2_instrument/vector_array.py` | Three-channel vector receive for INS-18. | Planned, WP11, goal level |

**Fidelity ladder at L2.**  The link budget answers "is there enough signal" in
closed form without generating samples.  The baseline receiver generates samples
and therefore can show what the processing actually has to cope with.  Both are
kept, because the closed-form answer is what makes the sample-level answer
checkable: if the measured SNR out of L3 does not match the link budget
prediction for the same configuration, one of the two is wrong, and that
comparison is a V1 test.

## L3  Signal processing

**Responsibility.**  De-chirp against the transmitted sweep, transform, map bin
index to apparent range, detect traces and form the ionogram.

**In.**  `ReceivedSignal`, optionally `PlatformState`.
**Out.**  `Ionogram`.

`virtual_range_m` is `(n_freq, n_modes)` with NaN where no trace was detected.
`snr_db` is parallel to it.  `doppler_mps` is populated only when the
implementation actually resolved Doppler, which is the whole point of C3; a tier
that leaves it `None` is saying it did not, and L5 will record INS-09 as
`not_evaluated` rather than as a failure.

| Tier | Module | Fidelity | Status |
| --- | --- | --- | --- |
| `baseline` | `src/tss/l3_signal/dechirp.py`, `transform.py`, `ionogram.py`, `detect.py` | De-chirp, windowed transform, bin-to-range mapping, CFAR-style trace extraction. | Not started.  Next package, WP4 |
| `updown` | `src/tss/l3_signal/doppler_decouple.py` | Motion-based separation of Doppler from range using up and down sweep pairing on one transmit and receive pair. | Not started, WP5, C3, principal contribution |
| `leakage_suppress` | `src/tss/l3_signal/leakage_suppress.py` | Adaptive cancellation of the de-chirped leakage residual. | Not started, WP5, C4 |
| `aoa` | `src/tss/l3_signal/aoa.py` | Arrival angle and polarization mode separation from the vector channels. | Not started, WP11, C5, goal level |

**Bin index to apparent range.**  This mapping is where the FMCW-specific
errors enter and it deserves its own note.  The de-chirped beat frequency maps
to apparent range through the sweep rate, so a sweep rate error, a sweep
nonlinearity and a Doppler shift all appear as range errors that are
indistinguishable in a single sweep.  `Ionogram.residual_range_bias_m` exists to
carry the part of that error the processing knows it has not removed, so L4 and
L5 can propagate it rather than silently absorbing it into the profile.

## L4  Retrieval

**Responsibility.**  Invert virtual height to true height by lamination, giving
electron density against true altitude with an uncertainty.

**In.**  `Ionogram`, `PlatformState`.
**Out.**  `DensityProfile`.

Topside inversion runs downward from the spacecraft.  The starting density is
the local plasma density at the spacecraft, which in a real mission comes from
an in-situ probe or from the lowest usable sounding frequency, and in the
simulator can be taken from truth to isolate the inversion error from the
starting-value error.  Which of those was done is recorded in
`DensityProfile.flags`, because the two give materially different error budgets
and a result that does not say which is not usable.

| Tier | Module | Fidelity | Status |
| --- | --- | --- | --- |
| `lamination` | `src/tss/l4_retrieval/lamination.py` | Standard laminar inversion with a chosen density law per lamina. | Not started, WP6 |
| `error_budget` | `src/tss/l4_retrieval/error_budget.py` | Propagates range, frequency and starting-value uncertainty into `ne_sigma` and `alt_sigma_m`. | Not started, WP6 |

L4 is gated by L3.  Until WP4 delivers an ionogram, there is nothing to invert
and an L4 tier would be tested only against synthetic ionograms, which tests the
inversion arithmetic and nothing about the instrument.

## L5  Assessment

**Responsibility.**  Score the retrieval against the L0 truth it started from,
and score the run against the requirement set, per configuration.

**In.**  `DensityProfile`, `IonosphereState`, `Ionogram`, `configuration`.
**Out.**  `Assessment`.

`Assessment.requirements` maps a requirement ID to a record carrying the
achieved value, the threshold and goal it was scored against, the verdict and
the margin.  The requirement set comes from
`config/requirements/requirements.yaml` and from nowhere else; see
`docs/requirements-traceability.md`.

`configuration` is mandatory.  Section 1.3 of the proposal reports results per
configuration, and a result without its configuration is not a result.  The
campaign validator refuses a campaign that names no configurations for the same
reason.

| Tier | Module | Fidelity | Status |
| --- | --- | --- | --- |
| `score` | `src/tss/l5_assessment/score.py` | Truth-versus-retrieval error metrics and requirement scoring. | In progress, WP0 and WP7 |
| `feasibility` | `src/tss/l5_assessment/feasibility.py` | Aggregates many runs into a feasibility map over the configuration space. | Not started, WP7 |

Requirements absent from the achieved mapping are reported as `not_evaluated`
rather than dropped, so an assessment always shows the full set and it is
obvious what a campaign did not reach.  That choice makes assessments more
verbose and it is worth the verbosity, because a silently shortened requirement
table reads as a clean sweep.

## Known seams

A seam is a place where the implementation does not match the architecture this
document describes.  They are listed rather than hidden, because an
undocumented seam is how a reviewer loses confidence in everything else.

### The migrated link budget is not split along layer lines

`src/tss/l2_instrument/linkbudget.py` is kept whole.  It arrived from the
author's earlier work as a single Python module with a MATLAB peer and 37
passing unit tests, and it contains, mixed together:

* transmitter power, antenna efficiency, receive chain noise figure and external
  noise, which properly belong to **L2**;
* the free-space and spreading terms and the absorption call, which properly
  belong to **L1**;
* the plasma frequency and reflection condition it uses to decide whether there
  is an echo at all, which properly belongs to **L0**.

The L0 and L1 packages therefore re-export the pieces that belong to them rather
than owning them:

```python
# src/tss/l1_propagation/__init__.py
from tss.l2_instrument.linkbudget import absorption_db, path_loss_db   # seam
```

**Why it is left this way.**  Splitting the module means moving code that is
covered by 21 passing tests, and the tests are written against the module's
current call signatures.  Doing that before WP4 lands would trade real, existing
coverage for tidiness at exactly the moment the codebase most needs a trustworthy
L2, because WP3 leakage work is the critical path and it builds on this module.
The risk of the split is not that it is hard, it is that a subtle transcription
error in the move would not be caught, since the tests would be moving at the
same time.

**What it costs.**  The import direction is wrong: L1 depends on L2, which
violates the rule that a layer may reference `interfaces` but never a sibling
layer.  A reader following the layer structure will not find the absorption
implementation where the structure says it should be.  And the link budget
cannot currently be run with an L0 tier other than the one it assumes
internally.

**When it is fixed.**  The physical split is scheduled for the start of WP4.
The order is: freeze the link budget tests, extract the L1 terms with the tests
duplicated on both sides until they agree, extract the L0 terms, then delete the
re-exports and the duplicated tests.  Until then the re-export lines carry a
`# seam` comment and this section is the explanation.

### The analytic L1 tier cannot produce off-vertical echoes

The analytic tier is vertical-incidence only.  Ground returns, oblique
propagation and the horizontal gradients that produce spread-F style traces are
outside it.  Any result about clutter, about arrival angle, or about the effect
of horizontal structure is therefore not addressable until the PHaRLAP tier
lands in WP2.  Results produced with the analytic tier state that limitation.

### Doppler is generated at L1 but exploited at L3

`EchoSet.doppler_hz` is filled by L1 from the platform velocity and the
reflection geometry.  The algorithms that separate Doppler from range live at
L3 and need `PlatformState`, which is why `form_ionogram` takes it as an
optional argument.  The seam is that L3 therefore has to reconstruct geometry
that L1 already knew.  This is deliberate: giving L3 the L1 geometry directly
would let the decoupling algorithm cheat, and the whole of C3 is the claim that
the algorithm can recover it from the received signal alone.  The cost is some
duplicated geometry code.

### Sub-megahertz noise is extrapolated

The external noise reference used by the link budget is measured down to about
0.6 MHz.  Below that it is extrapolated.  The band edge of 0.33 MHz therefore
rests on an extrapolation, and C1 is stated as resolved only where a measured
reference exists.  WP12 exists to replace the extrapolation with a measurement.
Until it does, any L5 verdict on INS-01 coverage below 0.6 MHz carries a flag
saying the noise reference was extrapolated.
