# Work packages

Thirteen work packages, WP0 to WP12.  The table gives the title, the current
status, the claims each package serves, the modules that implement it and the
next action.  Module paths that do not yet exist are marked planned; they are
listed anyway because the path is part of the plan and a reader should be able
to find where a piece of work will land.

Status values are the same four used in module metadata blocks: complete,
partial, not started, and goal level for the packages that are desirable rather
than required.

## The table

| WP | Title | Status | Claims | Modules | Next action |
| --- | --- | --- | --- | --- | --- |
| WP0 | Close the science requirements | Complete.  Requirement set at Rev E. | C1 C2 C3 C4 C5 | `src/tss/requirements.py`, `config/requirements/requirements.yaml`, `scripts/import_srm.py` | Keep Rev E in step with the SRM spreadsheet.  Re-open only if a claim forces a threshold change, and if so bump to Rev F rather than editing Rev E in place. |
| WP1 | Ionosphere truth model and perturbation library | Partial.  IRI-2020 in use; NeQuick 2 and the depletion generator outstanding. | C1 C2 | `src/tss/l0_ionosphere/iri2020.py`, `chapman.py`, `tabulated.py`, `matlab_iri.py`; planned `nequick.py`, `perturbations.py` | Implement the NeQuick 2 tier so that V2 has a second independent topside formulation.  Then the perturbation library: travelling disturbances, horizontal gradients, depletions. |
| WP2 | Propagation and echo model | Partial.  Absorption delivered and tested; ray-traced tier and off-vertical echoes outstanding. | C1 C2 C3 | `src/tss/l1_propagation/analytic.py`, `absorption.py`, `magnetoionic.py`; planned `raytrace.py`, `doppler.py`; MATLAB `matlab/tss_raytrace_vertical.m` | Stand up the PHaRLAP tier behind the L1 protocol and check it against the analytic tier at vertical incidence.  That comparison is the V2 case for L1. |
| WP3 | Instrument model and link budget | Link budget complete.  Leakage and de-chirped residual outstanding. **Critical path.** | C1 C4 | `src/tss/l2_instrument/linkbudget.py` (migrated, 21 tests), `receiver.py`, `waveform.py`, `antenna.py`; planned `leakage.py` | Model transmit-to-receive coupling as a frequency-dependent isolation profile so INS-23 can be scored, then the de-chirped leakage residual for INS-24.  Nothing downstream is trustworthy until the residual floor is known. |
| WP4 | Baseline signal processing chain | Not started.  Next package to start. | C1 C2 C4 | Planned `src/tss/l3_signal/dechirp.py`, `transform.py`, `ionogram.py`, `detect.py` | De-chirp against the transmitted sweep, transform, map bin index to apparent range, extract traces.  Do the link budget layer split described in `docs/architecture.md` at the same time. |
| WP5 | Novel algorithms | Not started. | C3 C4 | Planned `src/tss/l3_signal/doppler_decouple.py`, `leakage_suppress.py` | Motion-based Doppler and range decoupling from up and down sweep pairing on one transmit and receive pair.  This is the principal contribution and it cannot start until WP4 gives it an ionogram to work on. |
| WP6 | Retrieval and error budget | Not started. | C2 | Planned `src/tss/l4_retrieval/lamination.py`, `error_budget.py` | Laminar inversion from the spacecraft downward, then propagate range, frequency and starting-density uncertainty into the retrieved profile. |
| WP7 | Trade study and design closure | Not started. | C1 C2 C4 | `src/tss/l5_assessment/score.py`; planned `feasibility.py`, `src/tss/studies/trade.py` | Sweep the configuration space and produce feasibility maps per configuration.  Needs WP3 and WP6 to give it real numbers rather than assumed ones. |
| WP8 | Verification and validation | Partial.  Link budget suite and archive comparison delivered. | C1 C2 | `src/tss/verify/`, `tests/`, `src/tss/data/spdf_isis.py`, `src/tss/studies/isis_validation.py` | Extend beyond the link budget suite and the ISIS-1 archive comparison as each new tier lands.  Add the V2 transcription verifier runs for every MATLAB peer. |
| WP9 | Prototype specification | Not started. | C1 C4 | Planned `src/tss/studies/prototype_spec.py`, `docs/prototype/` | Turn the closed design into a bench specification: software-defined radio choice, isolation budget, antenna and test plan.  Depends on WP3 leakage results. |
| WP10 | Plasma bubble modelling | Not started, goal level. | C2 | Planned `src/tss/l0_ionosphere/bubbles.py` | Equatorial depletion structures as an L0 perturbation, to test retrieval under strong horizontal structure. |
| WP11 | Vector receive and arrival angle | Not started, goal level. | C5 | Planned `src/tss/l2_instrument/vector_array.py`, `src/tss/l3_signal/aoa.py` | Three-channel receive model for INS-18, then arrival angle and O/X separation from the channel set. |
| WP12 | Sub-megahertz sky noise reference | Not started, newly created. | C1 | `src/tss/data/wspr.py` (delivered), planned `src/tss/studies/sky_noise.py` | Establish a measured external noise reference below 0.6 MHz, or establish that none exists and state the band edge result as conditional.  This is the package that decides whether C1 can be closed at the bottom of the band. |

## Claim coverage

Read the other way round, so that it is obvious which packages a claim depends
on and where a claim has no work behind it yet.

| Claim | Work packages | Status of the evidence |
| --- | --- | --- |
| C1 sub-watt from small antennas gives usable ionograms | WP1, WP2, WP3, WP4, WP8, WP12 | Affirmative where a measured noise reference exists.  The link budget and the ISIS-1 archive comparison carry it above 0.6 MHz.  Below that, WP12 is the whole of the evidence and it has not started. |
| C2 range, density precision and along-track resolution together | WP1, WP2, WP4, WP6, WP7, WP10 | Open and quantified.  The trade space is understood; the numbers need WP6 to produce a retrieval error and WP7 to sweep it. |
| C3 motion separates Doppler from range | WP2, WP5 | Open, unstarted.  Principal contribution.  Nothing is implemented, and WP5 is blocked on WP4. |
| C4 leakage is suppressible with commercial SDR hardware | WP3, WP5, WP9 | Open.  This is now the claim the concept turns on.  WP3 has to produce a residual leakage floor before WP5 can attempt to suppress it. |
| C5 vector receive gives arrival angle and mode separation | WP11 | Open, goal level.  Not started.  If the schedule tightens, this is the first thing to drop. |

## Sequencing and the critical path

The critical path runs WP3 leakage, then WP4, and WP4 gates WP5, WP6 and
therefore everything at L4 and beyond.

```
  WP0 ──────────────────────────────────────────────┐ (requirements, done)
                                                    │
  WP1 ──┐                                           │
        ├──► WP2 ──┐                                │
  WP12 ─┘          │                                │
                   ├──► WP3 leakage ──► WP4 ──┬──► WP5 ──► (C3, C4)
  WP3 link budget ─┘      CRITICAL            │
        (done)                                ├──► WP6 ──► WP7 ──► WP9
                                              │
                                              └──► WP8 extends alongside
  WP10, WP11  goal level, off the critical path
```

**Why WP3 leakage is critical rather than WP4.**  The instrument transmits
continuously while it receives.  Until the transmit-to-receive isolation and the
de-chirped leakage residual are modelled, the noise floor that L3 has to detect
traces against is unknown, and it is very likely dominated by leakage rather than
by external noise across much of the band.  Building a signal processing chain
against a thermal-noise-only floor would produce detection performance that does
not survive contact with the real floor, and every downstream number would have
to be redone.  The order is therefore: establish the floor, then build the
processing that has to live with it.

This is also why C4 is now the claim the concept turns on.  C1 was the obvious
risk at the start, and the link budget together with the ISIS-1 archive
comparison has largely answered it above 0.6 MHz.  What replaced it is the
question of whether commercial software-defined radio hardware can hold enough
isolation, and whether processing can suppress what remains.  If C4 fails, the
sub-watt result of C1 does not matter, because the leakage rather than the
signal budget sets the detection floor.

**Why WP4 is the next package to start.**  It is the widest gate.  WP5 needs an
ionogram to apply its algorithms to.  WP6 needs an ionogram to invert.  WP7
needs a retrieval to score.  L4 does not exist in any meaningful form without
L3.  WP4 is also the moment to do the link budget layer split described in
`docs/architecture.md`, because that is when the L1 and L2 boundaries are being
exercised anyway and when a transcription error in the move would be caught by
the new L3 tests rather than sitting undetected.

**WP12 sits in parallel and is time-sensitive.**  It does not block anything
downstream, but it decides whether the band-edge result, which is already
published within the proposal as 0.33 MHz, can be stated unconditionally.  The
wspr.live data acquisition has a long lead time if a substantial history is
wanted, so starting the collection early costs little and buys options.  If WP12
finds no usable measurement below 0.6 MHz, the honest outcome is to state the
band edge as conditional on extrapolated noise, and that statement has to be
written into the results chapter rather than discovered at the defence.

**Goal-level packages.**  WP10 and WP11 are off the critical path and serve C5
and the harder parts of C2.  They are worth doing if the schedule allows.  The
decision point is when WP7 closes: if the trade study shows the baseline
configuration closes comfortably, there is room for them, and if it does not,
the effort goes into WP3 and WP5 instead.
