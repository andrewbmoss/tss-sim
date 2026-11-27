# tss-sim: where things stand

Updated 2026-11-27, during WP5 Phase 3 (curved-Earth facets, the sheet wrap,
the drift tier, the onboard prior, the windowed detector with its rate per
sequence, the Cramer-Rao harness, the joint estimators, the power-cut
experiment, the diffuse halo and Poole's comparator done; the study next; the
work on the fading outliers stopped for now).
WP0, WP3 and WP4 are complete; WP5 is in progress (Phases 0 to 2 done, Phase 3
under way, see the last section) and WP5 and WP6 are on the critical path; WP9 is optional and depends
on funding.  The history of what was built and why is in the git log,
the phase log below and the project notes; this file says what the code does
today, what it does not, and what comes next.

## Baseline

* Tag `proposal-baseline` marks `f103fc8`, the code the proposal and the
  Technical Volume cite: 474 tests, `tss run smoke` 0 goal, 8 threshold,
  2 fail (INS-04, INS-08), 10 not evaluated.  The INS-04 failure there was a
  simulation artefact; see the review fixes below.
* Run everything in the `tss` conda environment (`conda activate tss`).  Plain
  `python` in PowerShell is a stand-alone 3.13 without the dependencies.
* `pytest` runs the unit and contract suites; `ruff check src tests scripts` is
  clean; `tss verify` lists what verification has and has not run.
* The remote is `origin` = https://github.com/andrewbmoss/tss-sim.git.

## What each layer is today

| Layer | Implemented tiers | Stubs | What it does not yet model |
|---|---|---|---|
| L0 | chapman, analytic, iri (IRI-2020, `model: pharlap` by default or `pyiri`), iri2020 (PyIRI 0.1.7); matlab_iri (PHaRLAP IRI-2020); iri_wrapper (iri2016, not built) | perturbations (WP10) | NeQuick 2, IRTAM as a tier, depletions, storm days |
| L1 | analytic: vertical incidence, magnetoionic group index for O and X, absorption, Doppler; clutter (WP4: ground echo, plasma resonances, spread F; on by default) | raytrace (WP2, PHaRLAP) | Z branch, lateral deviation, off-vertical echoes, ray-traced clutter |
| L2 | baseline, built from `config/instrument/baseline_6u.yaml`: de-chirped baseband from the A1 stepped sweep at the INS-06 digitiser rate, in-dwell dispersion, sky noise, receiver gain, clipping, quantisation; transmit-to-receive leakage, transmit impairments (`transmit_impairments`), the receive chain (front end, per-step high-pass) and the passive canceller (WP3; defaults AM -140 dBc/Hz and an ideal 26 dB canceller); HF broadcast carriers (WP4 clutter) | - | real vector channels (WP11), measured hardware behaviour (WP9) |
| L3 | baseline (WP4): per-step transform with the chirp-rate search (64 cells, a rate per detection); the stated-P_fa detector (1e-4 per cell); the trace tracker (track-before-detect, O and X labelled by self-consistency, ground-echo and streak masks, cusp model, per-sample sigma, confidence and quality; legacy threshold-and-connect kept); SNR measurement at given cells; Poole's stationary-phase range estimator | - | O/X by measurement (WP11; labelling is model-based), the drift mode over a sequence of sweeps (WP5), measured traces with steep starts and ledges |
| L4 | lamination (field-free), lamination_magnetoionic (Budden 12.13, the default); starting value from the tracker's onboard f_N (`start_policy: onboard`, WP4) | - | uncertainty propagation, depletion crossings, reference-field error |
| L5 | standard: scores the Rev G set, reports what it does not score | feasibility (WP7) | NmF2/hmF2 errors, 1 sigma coverage, Kp classes |

The C3 work (route one, route two, the field-aligned sheet, the tilt map, the
arrival-angle options) runs as studies at parameter level, outside the L0-L5
pipeline.  The pipeline itself runs one sweep at one epoch.

## Step 1 (2026-10-17): requirements and contract tests

* `scripts/import_srm.py` transcribes the SRM workbook into
  `config/requirements/srm_catalog.yaml`: all 70 Rev G requirements
  (SCI 28, INS 28, MIS 10, SW 4).  Generated; never edit it by hand.
* `requirements.yaml` scores 22 of them and lists the other 48 under
  `not_scored` with a category and a reason.  `tss requirements` fails if any
  catalogue ID is in neither place.  SCI-15 (proposal Table 3.1, missing until
  then) and INS-03 were added.  The L5 report ends with every not-scored
  requirement and why.
* To take in a new SRM revision: rerun `import_srm.py`, review the diff of the
  catalogue, amend the scored entries.  A unit test fails if the catalogue and
  the workbook in `..\Dissertation` disagree.
* `tests/contract/test_tier_contracts.py` holds every implemented tier to the
  four rules of `tss/interfaces.py`, and fails if a tier is registered that it
  does not cover.

## Step 1 addendum (2026-10-17): the review fixes

Decisions by Andrew: the matching-network loss belongs in the antenna
efficiency, and the antenna efficiency follows the sounding frequency; INS-10
is scored at every truth echo; fix before closing WP0.

1. **Velocity frame.**  `PlatformState.velocity_ned_mps` (was
   `velocity_ecef_mps`) is local north-east-down; `tss.platform.orbits`
   converts with `ned_from_ecef`.  Before, the Earth-fixed components were read
   as NED and a circular orbit had a 7.5 km/s vertical closing speed in L1.
2. **Sample rate.**  L2 runs at `receiver.sample_rate_hz`, or by default the
   INS-06 rate `4 K_f R / c` at the peak sweep rate with `R` the INS-04
   threshold (78.4 kHz for A1).  Before, it was 2.5 times the deepest
   simulated echo's beat, so INS-04 always scored 1.25 times that echo and an
   echo-free record was sampled at 2.5 Hz.  Echoes beyond the unambiguous
   range fold, as in the instrument, and are counted.
3. **Instrument file.**  Campaigns build L2 with `build_instrument` from the
   instrument file plus each configuration's overrides (antenna length,
   vector receive).  `transmit.efficiency_db` is 0 dB, and the antenna
   efficiency, matching-coil loss included, is the ShortDipole physical
   model's at every sounding frequency: -46 dB at 0.33 MHz rising to -0.5 dB
   at 20 MHz for the 5 m antenna, where the file stated a flat -20 dB.  A
   configuration's antenna length changes it (`leo_500km_small_antenna`
   flies 1 m).  The pipeline now passes the platform to L5, so SCI-09 is
   scored.
4. **INS-10.**  Scored as the worst SNR per range cell that L3 measures at
   every echo the truth returned, both modes, detected or not
   (`snr_db_min_at_truth_echoes`).  The old minimum over detected cells could
   not fail and is kept as a diagnostic.

Smoke after the fixes: 0 goal, 11 threshold, 1 fail (INS-08), 10 not
evaluated.  INS-04 meets threshold at 1050 km; SCI-09 is 32.8 km; INS-10
meets threshold at 19.2 dB, the X echo at 6.64 MHz just below penetration
(O minimum 19.7 dB), over 44 truth echoes of which 25 were detected.  With
the flat -20 dB efficiency the same run failed INS-10 at 9.4 dB.
`field_inversion` through the whole chain gives 4.5 to 5.0 km for the
magnetoionic inversion and 1.3 to 17.9 km field-free, within 0.1 km of the
numbers it gave before, because the old L2 defaults already used the physical
antenna model.  Technical Volume Table 7.1 (`tab:fieldinversion`) was updated
to these numbers on 2026-10-17.

## IRI-2020 through PyIRI (2026-10-17)

PyIRI 0.1.7 is installed in the `tss` environment (and listed in
`environment.yml` and the `iri` extra) and the `l0: iri2020` tier is bound to
its spherical-harmonic library with IRI-2020's default options: URSI foF2,
SHU-2015 hmF2, ABT-2009 bottomside, NeQuick-2 topside.  It needs a stated daily
F10.7 (`f107`); the `baseline` campaign states 120.  Against PHaRLAP's Fortran
IRI-2020 at matched solar level (`docs/pyiri-validation.md`) NmF2 agrees to
about 3 per cent and hmF2 to about 2 km, but the topside density at 500 to
900 km is only 0.44 to 0.51 of IRI-2020's, so the plasma frequency at the
spacecraft is about 35 per cent lower at 900 km.  Traced on 2026-10-18 (`docs/topside-model-campaign.md`): PyIRI 0.1.7 applies the NeQuick G thickness step (H0 = k B2bot / v) on top of the NeQuick 2 shape factor, which IRI-2020 does not, so its H0 is 0.63-0.66 of IRI-2020's; ISIS-1 sides with IRI-2020. Do not use `iri2020` for band-edge or topside-shape work until corrected.  The
`iri2016` wrapper is installed but not built (it needs ninja on first call), so
`matlab_iri` is the only Fortran IRI route that works today.

Re-run against ISIS-1 with each epoch's own indices (R12 = -1, so IRI reads R12, IG12, daily and
81-day F10.7 and ap from its index files, foF2 storm model on; PyIRI at the daily F10.7):
`scripts/topside_campaign.py isis-epoch`, `docs/topside-model-campaign.md`.  At solar maximum nothing
changes (IRI-2020 within 0.07 dB, PyIRI within 0.31 dB, fitted H0 still 0.93 of the IRI law).  At
solar minimum the user-R12 route's derived IG12 (9 for R12 15) is well below the file's (17), which
put the binned densities 16% low: worth checking in the band-edge study.

## IRI implementation as a configuration option (2026-10-18)

Campaigns name `l0: {name: iri, f107: ...}`.  The `model` key picks the code:
`pharlap` (the default; Fortran IRI-2020 in the PHaRLAP mex via MATLAB, i.e.
`matlab_iri`) or `pyiri` (`iri2020`).  The solar index is stated once, `f107`
or `r12`, and converted with IRI's own relation (63.75 + 0.728 R12 +
0.00089 R12^2), so `f107: 120` runs PHaRLAP at R12 = 71.1 and IRI reports
F10.7 = 120 back.  Other keys go to the chosen implementation, and an option
it does not take is an error.  `label`, `provenance` and `describe()` record
the model.  `baseline` now runs `model: pharlap`, so it needs MATLAB and
PHaRLAP; `smoke` and the contract suite still need neither.  The `iri2020`
and `matlab_iri` tier names are unchanged.  Tests: `tests/unit/test_l0_iri.py`
(the PHaRLAP run is `-m pharlap`).

The switch changes the baseline result.  `tss run baseline` (strict, clean tree)
on PHaRLAP against the same campaign with `--set tiers.l0.model=pyiri`:

| Configuration | INS-10 worst SNR at truth echoes, PHaRLAP | PyIRI | Fails, PHaRLAP | PyIRI |
|---|---|---|---|---|
| leo_500km_single | 12.3 dB | 25.9 dB | INS-08 | INS-08 |
| leo_800km_single | 7.0 dB | 18.7 dB | INS-08, INS-10 | INS-08 |
| leo_500km_small_antenna | 4.1 dB | 9.5 dB | INS-08, INS-10 | INS-08, INS-10 |

So at 800 km INS-10 fails once the topside is IRI-2020's.  The PyIRI column
reproduces the pre-change baseline exactly.  Traced in WP3 Phase 0
(2026-10-20, `docs/reference-echo.md` §6): the minimum is the X echo in the
last step below fxF2 on both codes, and the topside does not move it.  The
median SNRs agree to 0.1 dB.  PHaRLAP's foF2 is 0.35 % lower, which puts that
step 10 kHz below fxF2 instead of 48 kHz, at 873 km of group range.  Moving
foF2 by one step moves the minimum through 5.6-13.2 dB on PHaRLAP.

## WP0 closure, Phase 0 and Phase 1 (2026-10-19)

Plan: project note `claude/wp0-closure-plan.md`.  Branch `wp0-closure`.

* Phase 0: the Rev G band-edge run (MinFreqStudy, 2026-09-07) is frozen in
  `data/reference/band_edge` with source hashes.  PHaRLAP's index files allow
  epoch-driven runs from 1958-01-01 to 2024-12-09 (`apf107.dat` ends there).
* Phase 1: `tss study min_frequency` is now the MinFreqStudy harness, ported
  (`src/tss/studies/min_frequency.py`, plans in `config/studies`, doc
  `docs/band-edge.md`).  The placeholder it replaces ran 400-800 km on a coarse
  grid and silently skipped IRI.  It writes the harness's summary columns plus
  SCI-02 columns, the ensemble (minimum across topsides), completeness curves
  and per-sample scalars.
* The `iri` tier and `MatlabIRIBackend` gained an epoch mode (`epoch_indices:
  true`, R12 = -1), with the epochs checked against the index files before
  MATLAB starts, and the extras now carry IG12, 81-day F10.7, daily Ap and the
  foF2 storm ratio.
* Fixed on the way: `iri_options` never reached MATLAB intact (the JSON string
  lost its quotes on the `-batch` command line); they now go as a JSON file.
  Nothing before today passed non-empty options, so no earlier result is
  affected.
* Tests: `tests/unit/test_band_edge.py` (the harness's 18 checks, the ensemble
  direction, orbit sampling against the harness's own draws, the plan loader,
  an end-to-end stand-in run, and two `-m pharlap` cases that match the
  reference to 1e-4).  Suite 632 passed, 1 skipped; `-m pharlap` 4 passed.
* Reference regression: the full 60-case rerun reproduces the Rev G summary to
  2.4e-7 relative in every percentile; governing value O1 / 900 km / solar min /
  NeQuick 0.3342 MHz (95 %), 0.2700 MHz (99 %).  Details in docs/band-edge.md.
* Phase 2 (2026-10-20): plan `band_edge_solarmin_epoch`, epoch-driven solar
  minima 1986, 1996, 2008-09, 2019-20 (foF2 storm on) and 2008-09 storm off;
  `scripts/band_edge_compare.py`.  Governing O1 / 900 km / NeQuick: 0.355,
  0.330, 0.282, 0.279 MHz at 95 % (Rev G 0.334); 0.212 MHz at 99 % in 2019-20
  (Rev G 0.270).  The Rev G level is a 1996-like minimum; the deep minima lower
  the edge 11-16 % (95 %) and 16-22 % (99 %) everywhere.  Storm model effect
  0.0-0.3 %.  Within INS-01a's -20 % tolerance, which Rev G reserved for the
  models' 19 % tail overstatement against ISIS-1.  G1 decided (Andrew,
  2026-10-20): keep the Rev G 0.33 MHz; record the deep-minimum value in
  INS-01a's notes at close-out.  Results
  kept in `data/reference/band_edge/phase2_solarmin_epoch`.
* Phase 3 (2026-10-20): SCI-02 threshold from the kept hmF2 samples
  (`tss.studies.sci02_margin`, `scripts/sci02_threshold.py`, outputs in
  `data/reference/sci02_threshold`), plus two storm plans (`sci02_storm_epoch`,
  `sci02_storm_m3000`).  Met over 600-900 km at every solar level; hmF2
  0.1/99.9 percentiles inside 200-500 km everywhere.  At 500 km, solar max,
  0.03-0.07 % of samples have the peak above the spacecraft (MIS-01 "at all
  times" not met); IRI's hmF2 has no storm uplift, so storm time is not
  verified.  G3 decision (Andrew): MIS-01 restated as a fraction of samples,
  99.9 % (`MIS01_ABOVE_PEAK_FRACTION`); worst cell 99.93 %, met everywhere.
* Phase 4 (2026-10-20): plan `band_edge_goal_1400`, 120 cases to 2000 km, with
  the O+/H+ transition extra and completeness at a 0.33 MHz start
  (`goal_metrics`); `scripts/goal_knockons.py`; outputs in
  `data/reference/band_edge_goal_1400`.  Governing 95 % edge 0.246 MHz at
  1000 km to 0.164 MHz at 1400 km (O1, deep minimum); a 0.33 MHz start covers
  85.6 % to 65.3 % (mean completeness 0.97 to 0.89).  O+/H+ median 820 km (O1)
  and 680 km (O2) at solar min, 1135-1290 km at solar max.  At 1400 km:
  INS-04 1800 km, INS-06 134 kHz (A1), INS-13 580 MB/day, +4.7 dB spreading,
  INS-10 -5 dB (provisional).  G2 (Andrew): requirements stay as they are
  unless the follow-on studies show something significant; the exceedances go
  into the TV and dissertation.
* Phase 5 (2026-10-20): `scripts/isis_high_altitude.py`, 5,233 ISIS-1 CRC
  profiles at 900-1400 km against epoch IRI-2020 (three topsides); outputs in
  `data/reference/isis_high_altitude`.  NeQuick median +0.3 to +0.7 dB, H0 on
  the IRI law; but the modelled 5th percentile is 23-46 % high at solar max
  (realised coverage of the INS-01a-style edge 86-90 %), mostly poleward of
  70 deg (ratio up to 2.1; confounded with CRC quality 6-7).  Bound on Phase 4:
  coverage at 0.33 MHz, 1400 km, deep min, 65 -> 45 % (O1), 75 -> 61 % (O2).
  G2 confirmed (Andrew): the results show the datasets are inadequate rather
  than that the requirements must change; requirements stay, MIS-01 stays at
  99.9 % of samples; the discussion goes into the TV and dissertation.
* Phase 7.1-7.3 and 7.7 (2026-10-20): `tss study data_volume`
  (`src/tss/studies/data_volume.py`, `docs/data-volume.md`): products per
  observing mode, threshold and goal product sets, reference duty cycles for O1
  and O2, A1 and A2, 900 and 1400 km.  A1 / O1 / 900 km: 109 MB/day threshold,
  435 goal; Rev G's 16.7 kB, 215 / 339 MB/day and trace figures reproduced;
  INS-19 vector factor x2.0 (traces), x2.3 (threshold set).  Drafts of INS-13,
  a new INS-28 and a proposed MIS-11 (G5) are in the doc.  Contacts, link
  budgets and processing budgets are out of scope (Andrew).
* `tss_iri_profiles.m` now writes the extras with nine significant figures.
* **WP0 closed (2026-10-20, Phase 8).**  Andrew accepted MIS-11 as drafted
  (G5).  `Dissertation/TSS_Science_Requirements_RevG.docx` and
  `TSS_SRM_RevG.xlsx` carry the 2026-10-20 amendment: MIS-01 as 99.9 % of
  samples, INS-13 and INS-19 restated, INS-28 and MIS-11 added, INS-01a and
  SCI-02 annotated, new sections 5.12 and 5.13, Open Items 17 and 18, SRM
  Derivations section 24 (live formulas, 108.6 / 434.6 MB/day as tss-sim).
  Catalogue regenerated (72 requirements); requirements.yaml lists INS-28
  (bench) and MIS-11 (mission).  Technical Volume updated (band-edge chapter,
  INS-13 section, Table 7.1, verification matrix).  Originals backed up in
  `%TEMP%/tss_wp0/revG_backup_2026-10-20`.
* G4 (Andrew, 2026-10-19): ground stations USU/Logan, KSAT SvalSat
  (78.230 N, 15.408 E) and INPE Cuiaba (15.553 S, 56.073 W).  Reference radios:
  UHF AAC Clyde Space PULSAR-UTRX (9.6 kbps GFSK, 27-33 dBm, 3-5.1 W Tx);
  S-band Innoflight SCR-104 (up to 4.5 Mbps, BPSK/QPSK/OQPSK, up to 2 W RF);
  X-band Innoflight SCR-106 (7.9-8.5 GHz Tx, up to 105 Mbps per SatCatalog,
  132 Mbps per SatNow, 20 mW-2.5 W RF, 5-30 W DC, 250 g).  Public datasheet
  values only; the Innoflight datasheets are on request.  Capacity is to be
  set from a link budget per station, not from the maximum rate.

## Open findings of the 2026-10-17 review

* L3 searches every bin out to each column's unambiguous range, about
  9300 km in the slowest steps, with a fixed 10 dB threshold, so noise alone
  can raise a detection.  WP4's robust extractor (CFAR and a search limit at
  the design range) is the fix.
* L2-owned numbers (INS-11, INS-14, INS-23, INS-24, INS-26) have no route into
  L5 yet.
* The oscillator is built but never applied to the samples, and
  `EchoSet.doppler_hz` is ignored by L2: the samples carry no phase noise and
  no Doppler.
* `config/campaigns/c3_doppler.yaml` still asks for coherent intervals of 1 and
  2 s and up/down sweep pairing, both withdrawn; `motion_separated` describes
  the single-pair model the route studies replaced.  (Both retired 2026-11-24,
  WP5 decision D9, ff9ecf1.)
* No configuration sits at the 900 km governing case of INS-01a and MIS-01.

## WP3 Phase 0 (2026-10-20): the reference echo

Plan: project note `claude/wp3-plan.md`.  Doc: `docs/reference-echo.md`.

* `tss.l2_instrument.reference_echo.reference_levels(instrument, sweep)` gives
  the following per step:
  * the noise in one range cell;
  * the weakest echo INS-10 admits (noise + 10 dB; INS-24's reference);
  * the free-space echo from 1050 km;
  * the leakage carrier at INS-23's isolation, and the echo band
    (100-1050 km beats).

  It uses the same terms the L2 baseline puts into the samples.  A synthetic
  echo through L2 and L3 measures the predicted SNR to within 1 dB
  (`tests/unit/test_reference_echo.py`, study `reference_echo`, outputs in
  `data/reference/wp3_reference_echo`).
* The SRM's "leakage over echo = 122 dB - isolation" is the inherited
  `SounderConfig` budget: 300 km, legacy -6.8 dBi antenna.  On the instrument
  at 1050 km it is 97.2 dB at 5 MHz at 40 dB isolation.  Against the weakest
  admitted echo it is 107.5-128.8 dB over the band (111.0 dB at 5 MHz).
  Carried into Rev G's rules:
  * INS-11's high-pass rejection is 47 dB at 5 MHz, not 32 dB (93 dB at the
    band edge).  Referred to the digitiser-band noise, it is 40 dB (57 dB).
  * INS-22's spur requirement against the weakest echo, at 10 m, is -82.6 dBc
    at the worst step, against Rev G's -70 dBc.

  The requirement documents were corrected on 2026-10-27 (below).
* The `link_budget` study's default `SounderConfig` is the inherited design,
  and the SRM's INS-10 status (38.5/33.1 dB) comes from it.  The study now
  also reports the instrument: 27.8 dB median at 1050 km in free space, below
  10 dB under 1.26 MHz on the 5 m antenna.
* `baseline_6u.yaml`'s -20 dB was already resolved and tested on 2026-10-17.
  The L2 baseline carries the Balanis reactance.
* The INS-10 minimum on PHaRLAP is the penetration cusp (above).
* `scripts/ins10_trace.py` (`echoes`, `cusp`) reproduces the trace.  The
  PHaRLAP outputs are in `data/reference/wp3_reference_echo/ins10_trace`.

## Next, in order

1. WP3, L2: the leakage model (Phase 1) and the de-chirped residual against
   INS-24, measured against `reference_echo`; phase noise and Doppler into the
   samples; a route for L2-owned numbers into L5.
2. WP4, L3: robust trace extraction, O/X labelling (INS-08), the Barry and
   Poole reference former.
3. Gate M2: SCI-05 and SCI-07 scored against truth with both modes present,
   over seeds.
4. Beside 1-3: the MATLAB bridge (`tss_yaml_lite.m`, `tss_call_python.m`), the
   PHaRLAP L1 tier, NeQuick 2 and IRTAM L0 tiers, the depletion generator,
   WP12's noise reference.
5. WP5: the drift-mode sequence path and its L3 tier (route one at signal
   level; `motion_separated` was retired 2026-11-24), then WP6,
   WP7, the V2-V4 campaign, and the goal-level packages.

## Rev G amended 2026-10-27: the leakage-to-echo ratio and INS-22

Rev G section 5.9 and SRM Derivations section 16 took leakage over echo as
122 dB less the isolation.  That is the inherited budget's 300 km echo on the
legacy -6.8 dBi antenna.  On the instrument the free-space echo from 1050 km is
137.2 dB below the transmitter at 5 MHz and 160.5 dB below it at 1.26 MHz, the
lowest step at which it still meets INS-10.

* INS-22 is now -82 dBc threshold and -92 dBc goal, for A1 and A2, over each
  step's echo band.  It binds at 1.26 MHz, where the rule asks -81.1 dBc at
  40 dB isolation and a 10 m path.  At 5 MHz the rule asks -65.5 dBc.
* The A2 values of INS-22 (-86.5 dBc) and INS-23 (56.5 dB) are withdrawn.  A2's
  fast steps are where the ratio is smallest, and its binding step is A1's.
* INS-11's high-pass rejection is 47 dB at 5 MHz (was 32 dB).
* `config/requirements/srm_catalog.yaml` is regenerated from the workbook.
  The Technical Volume (sections 5.6, the WP9 bench numbers and the BOM) and
  its leakage figure are updated.  The Proposal is not.

## WP3 Phase 1 (2026-10-27): the leakage model

Doc: `docs/leakage.md`.  Decisions G1 (parametric isolation brackets, one
moment-method anchor) and G2 (coupling paths 1, 3, 10, 30 m; baseline 3 ns).

* `tss.l2_instrument.leakage.leakage_signal` makes the de-chirped leakage in the
  baseline's units and phase convention, with the exact step transitions, a
  frequency-dependent isolation profile (flat, ripple, resonances, a table,
  retuning transients), delay jitter, and a transmit phase error with the
  INS-21 common-reference switch.
* `tests/unit/test_leakage.py` (40 tests): the spur suppression
  4 sin^2(pi f_m tau), the sideband (pi Df tau)^2, the ripple's AM sidebands,
  2 pi f sigma_tau, and leakage over the Phase 0 echo all hold to 0.011 dB
  against the 0.5 dB exit criterion.
* `data/reference/wp3_isolation_mom`: two 5 m dipoles, coil Q 150, free space.
  Below a few MHz the isolation is set by coupling coefficient times Q, not by
  efficiency, and is flat: 29.0 dB crossed with a 3-8 cm feed offset
  (k = 4.7e-4), 0.3 dB parallel.  On this model INS-23's 40 dB needs
  k < 1.3e-4.  No bus in the model.
* Not yet: the leakage is not added to the baseline's samples, and there is no
  receiver compression or leakage high-pass (Phase 3); no residual or score
  (Phase 4).

## WP3, 2026-10-27 continued: the bus anchor and Phase 2

* `scripts/matlab/wp3_isolation_mom_bus.m`: two 5 m dipoles on a closed
  conducting 6U box, full MoM.  Isolation is the coupled-resonator law
  `-20 log10(kQ/2)` (`leakage.coupled_resonator_isolation_db`).  Crossed with a
  3/8 cm feed offset: 36.6 dB (k = 2.0e-4); parallel 0.4 dB.  Angle error alone
  on a symmetric crossed pair: 0.1 deg gives 36-40 dB, 1 deg 14-21 dB.  INS-23's
  40 dB needs k < 1.3e-4.  `mom_anchor_model()` defaults to the bus mount.
* `tss.l2_instrument.impairments` (Phase 2): `TransmitImpairments` (oscillator
  noise from the Allan deviation, a phase-noise floor, AM noise, PM and AM
  spurs, sweep nonlinearity, Doppler; INS-21 as `single_reference`) and one
  `ImpairmentRealization` per record shared by the echoes and the leakage.
  The L2 baseline applies them to the echoes; `leakage_signal` takes them.
  `tests/unit/test_impairments.py` (19): every prediction within 0.2 dB; the
  INS-21 credit lost is 50.43 dB against 50.47 dB predicted.
* Configuration: `transmit_impairments` in `baseline_6u.yaml`, all off, so no
  earlier result moves.

## WP3 Phase 3 (2026-10-28): leakage in the samples, front end, high-pass

* `leakage.enabled` in `baseline_6u.yaml` (false by default until Phase 4)
  adds the leakage at the receiver input on the MoM bus isolation, sharing the
  record's impairment realization.  `front_end` (Rapp, +1 dBm P1dB) and
  `leakage_highpass` (Butterworth 4, corner 0.5 of each step's 100 km beat)
  follow it.
* `tss.l2_instrument.receive_chain`: `FrontEnd`, `LeakageHighpass`,
  `small_signal_gain_db`, `compression_image_dbc`.  Bandpass equivalence
  checked against a passband tanh to 0.1 dB.
* Study `leakage_receive_chain` (`data/reference/wp3_receive_chain`): INS-11 on
  the noise rule met at every step (least margin 38.6 dB cold, 92 dB primed).
  Cold start loses step 1 for 25 ms and clips the converter by 50 dB there.
  The front end is 7.6 dB under P1dB below 3 MHz (WP9-2 asks 10), with an echo
  image at -41 dBc.

## WP3 Phases 4 to 7 and close-out (2026-10-28 to 2026-10-29)

* Phase 4 (f75e811): `leakage.dechirped_residual` and `score_ins23_ins24`;
  `l2_metrics` on every record; L5 maps INS-11, -14, -23, -24, -26.  Leakage on
  by default (MoM bus isolation, primed).  Finding: the transmitter's AM noise,
  carried on the leakage and not removed by range correlation, sets the INS-24
  pedestal; at -133 dBc/Hz INS-24 failed by 26.7 dB at 0.338 MHz.
* Phase 5 (a37ad40, db35b7f): `leakage.RFCanceller` (amplitude, phase and delay
  errors; retune track / step / sweep; added noise) and study `leakage_budget`.
  Allocation: passive canceller >= 26 dB with a continuously moved weight
  (INS-31) and AM <= -140 dBc/Hz (INS-30); INS-29 (leakage-borne noise <=
  k T_sys) added.  Default: AM -140, ideal 26 dB canceller.
* INS-23 relaxed to 35 dB (0931733); INS-10 restated with each mode's cutoff
  step reported separately (G5c, e881594).
* Phase 6 (8df56b7, 79c059f): study `wp9_bench`, the predictions and pass
  limits of the WP9 conducted bench (paths E0/E1/E2 = 39/61/69 dB emulating
  INS-23 + INS-31 by attenuation).  WP9 Bench Test Specification Rev A in the
  Dissertation folder.  WP9's conducted bench is required from 2026-10-29.
* Phase 7 (4667679, 5fed541): `leakage_budget` scores INS-10 without the
  cutoff-step echoes, and its reference run is regenerated (ideal instrument
  24.9 dB; no echo below 10 dB at any power in the grid).  The Technical
  Volume's leakage section is rewritten around the result, with two new
  figures from `wp3_isolation_mom` and `wp3_leakage_budget`.  Rev G and SRM
  status fields point to the bench tests; `srm_catalog.yaml` re-imported.
* **WP3 closed (2026-10-29).**  C4 is answered in simulation; its two hardware
  allocations (INS-30, INS-31) are measured by WP9.

## WP4 Phases 0 and 1 (2026-10-30)

* Phase 0 (af28d58, 7ba8814): study `wp4_campaign`, six Chapman ionospheres
  (mid-latitude noon at solar max and min, mid-latitude night, equatorial
  evening, high-latitude night, the smoke ionosphere) x five seeds, scored on
  SCI-05, SCI-07, INS-05, INS-10 and the O-trace completeness (with the X,
  false, missed and spurious counts and the range bias).  Baseline reference
  run `data/reference/wp4_baseline`: SCI-05 0.176 and SCI-07 18.8 km median,
  O completeness 0.14 (the chain follows X above f_x on half the O columns and
  reads 5 to 15 km short).
* Phase 1 (9fdeb53): the clutter generator, `tss.l1_propagation.clutter`
  (`docs/clutter.md`): ground echo above penetration, plasma resonances at the
  spacecraft (f_N, f_UH, n f_H), spread F, broadcast carriers above foF2
  (rendered by `tss.l2_instrument.interference`); the compression image listed.
  Flagged EchoSet rows plus `EchoSet.interference`; truth scoring ignores them.
  On by default (`clutter` block in `config/default.yaml`), pinned off in the
  smoke campaign, optional climatology.  The second hop is left out.  Clutter
  off reproduces the baseline in all 30 runs.  Clutter on
  (`data/reference/wp4_p1_clutter`): O completeness 0.00, INS-10 median 5.9 dB.
  The ground echo takes the extractor (P4); the resonances and a carrier in an
  X column lower INS-10 (gate G1: score INS-10 against noise or against
  noise plus clutter).
* Phase 2 (e841a41, d2f6d3e, dfe8d2c): the detector, `tss.l3_signal.detect`
  (`docs/detector.md`), now the default L3 cell test (`detector: cfar`; the
  smoke campaign is pinned to `fixed`).  G1: P_fa 1e-4 per cell against a
  noise reference from the echo-free negative-beat half (ordered-statistic
  median, Rohling's threshold); sidelobe blanking under a stronger peak's
  window envelope plus 6 dB; the ordered-statistic veto optional and off (it
  masked range-extended echoes); search to 1050 km, no near-range exclusion.
  Measured P_fa 0.8e-4 to 1.25e-4; Pd within 0.2 dB of Marcum Q; about 12
  false alarms per ionogram reach the extractor.  INS-10 against the noise
  floor, signal to noise plus clutter beside it: with clutter on the INS-10
  median is 18.1 dB (P1 measure: 5.9).  Neutral on the trace with clutter off;
  the ground echo still takes the extractor (P4).
* Phase 3 (0feb6ac, 196e895, fe20462, 5bd021a, 2d6854b): in-dwell dispersion
  and the step transition (`docs/dispersion.md`, study `wp4_dispersion`,
  `data/reference/wp4_p3_dispersion`).  Treatments in the baseline L3 tier:
  `dispersion: none | blocks | rate_search`, `transition: accept | gate |
  matched`, and the narrow step law (`narrow_lo_hz`, `narrow_hi_hz`,
  `narrow_cap_hz`).  Echo level over 18 scenes (6 Chapman, 12 ISIS-1) x O, X x
  3 seeds, and through the chain, clutter off.  G2: the rate search, span
  widened to 64 cells, now the default (`rate_search_max_cells: 64`; smoke
  pinned to `dispersion: none`); loss 0.10 dB median and 1.24 dB p90 for echoes
  crossing up to 64 cells (94 %), the rest to P4's cusp model.  G3: accept the
  skirt (gating fails INS-05; the matched filter equals accept in the chain).
  The detector keeps 1e-4 per cell: gamma threshold over the median for
  sub-blocks; an effective look count for the rate search (about 13 at 64).
  G6: draft INS-32, a performance condition on the processing, for Rev G at
  P6.  Clutter off (`wp4_p3_clutter_off`): SCI-05 0.145, SCI-07 17.4 km,
  INS-10 median 20.2 dB, 3 fails.  P2's mid-latitude-night INS-10 fail was
  dispersion, not sky noise.  Rate capture (one rate per column, set by the
  strongest return) causes the 3 remaining fails and, with clutter on, lets
  the resonances and ground echo take the rate: rate per detection is a P4
  item.  Clutter on (`wp4_p3_clutter`): the ground echo still takes the
  extractor.  L3 run time about 6 times D0's.
* Phase 4 (0c8ec62, 1cbbcc3, fb9a671, 62076e6, and the docs commit): the trace
  tracker, `tss.l3_signal.tracker` (`docs/tracker.md`), the default extractor
  of the baseline L3 tier (`extractor: legacy` kept; smoke pinned to it).  The
  rate search assigns a rate per detection (`rate_assign`, 15 effective looks).
  Onboard priors (IGRF f_H; f_N from the resonances or the trace start),
  streak and ground-echo masks (foF2 from the ground echo), track-before-detect
  (Viterbi, evidence min(z,1000) - 6, background flattening, 12 dB shadow
  test), O and X ordered and labelled by self-consistency (X predicted from the
  inverted O trace, `predict_x_range_m`), leading edges refined, per-sample
  sigma, confidence and quality bits, a cusp model.  G4: reject + flag.
  `wp4_campaign` reports the labelling error and `coverage_gap_m`.  Reference
  runs `wp4_p4_clutter` (SCI-05 0.063, SCI-07 7.9 km, all goal or threshold,
  labelling error 5.7 %) and `wp4_p4_clutter_off` (0.031, 5.1 km, 0.4 %).
* Phase 5 (58f9013, cd82baa, 40b8171, and the docs commit): acceptance
  (`docs/acceptance.md`).  `tss.l3_signal.poole_range`, Poole's three-cell
  estimator (phase velocity at the window centroid; R' reads 2.5 V* T long, a
  term Poole drops).  Study `fmcw_acceptance` (`data/reference/wp4_p5_acceptance`):
  items A1-A17 from Poole (1985), Poole and Evans (1985) and Barry (1971) at
  their parameters, scored on G5: 12 reproduced, 2 differences explained (A5a,
  A6), A8 reported, 5 not simulated, none open.  Study `isis1_retrieval`
  (`data/reference/wp4_p5_isis1`): 7,790 ISIS-1 records (quality 4-9) x clutter
  on/off through the P4 chain, resumable from a JSON-lines checkpoint (662 CPU
  h).  SCI-05 0.078 / 0.092 median, 28 % fail, concentrated at low f_N (49 %
  under 1.5 MHz, 2 % over 4 MHz); error over the measured span 6.4 % median;
  labelling error 0.16 median.  Open for P6.  `wp4_dispersion.isis_scenes` kept
  the worse-scaled records (index >= 6; 0 is best): default now 4-6, reruns
  `wp4_p5_dispersion_q4to6` and `_q6to10`; G2 and G3 stand.
* Phase 6 (b6dfbcc, 6e4fc9b, 847f4c5, and the docs commit): the M2 campaign
  (`docs/m2.md`).  `isis1_retrieval` scores SCI-05/07 over the measured span
  beside the full profile, with the M2 rule (median at threshold, 80 % of runs
  at threshold).  Baseline (P4 chain, `wp4_m2_base_*`): ISIS-1 fails, 73.1 /
  65.3 % clutter on.  Tracker: `fn_prior: trace_floor` (f_N bounded by the
  lowest persistent trace start; far O starts allowed; the f_N resonance
  taken at or below the floor), `far_start_cost`, options `stop_at_loss` and
  `cusp_up_cap_m` (off).  L4: `start_policy: onboard` (the tracker's f_N).
  Final (`wp4_m2_isis1`, `wp4_m2_chapman_*`): ISIS-1 84.6 / 89.3 % on,
  81.8 / 87.6 % off; Chapman all goal.  M2 met.  `scripts/isis_model_check.py`
  (`wp4_m2_models`): IRI-2020 against ISIS-1 and the coverage of the data.
  The `l4_retrieval.base` starting-value docstring had the bias sign wrong.

## WP5 Phase 0 (2026-11-24): stubs retired, the contract, the harness

* Decisions and plan: project notes `wp5-plan-draft`,
  `wp5-decisions-and-p0-start-2026-11-24`.  C3 by three routes; route one is
  the best solution with one receive channel; M3 is route one alone at signal
  level.
* D9 (ff9ecf1): `motion_compensation`, `aperture_synthesis`, `leakage_cancel`,
  the `motion_separated` tier, `c3_doppler.yaml` and the `doppler_separation`
  study retired; `novel/sweep_pairing` kept as the Phase 5 INS-09 stub.
* `docs/drift-mode.md`: the single-node output contract (every echo, compact
  facets with the several-frequency lever, field-aligned sheets, the extended
  layer, the ionogram-mode INS-09 correction), SCI-15 scored against the
  plasma's own rate (`requirements.yaml` notes), the D8 and C4 screens.
* Study `wp5_campaign` (fe5dd5d), the WP5 harness: 115 scenes by echo class
  (c3_verify's grids A 56, B 32, C 12, D 7; eight field-aligned sheets, two in
  the plasma) through the parameter-level code, and the route-one budget with
  the `route_one` study's seeds.  Each scene has its own seed (campaign seed and
  a checksum of its id) in `manifest.seeds`.  Expected values are pinned in
  `EXPECTED` from the TV; mismatches are listed, not tolerated.  Rows use the
  contract's names; products with no parameter-level logic are listed in
  `not_evaluated`.  Sheets are at 5 MHz and 500 km only (the forward model
  fixes them); plasma facets are placed by `mu_f`, with `P'` recorded.
* Reference run `data/reference/wp5_baseline` (438f124; from fe5dd5d, clean,
  `--strict`): reproduced.  A 56/56; B 0/32 flagged, error 0.9994-1.001 of the
  predicted tilt bias; C 12/12 declined; D compact to 2 km, extended from 4 km;
  sheets to TV Table `c3sheet`; budget within 0.4 % of TV Table `c3r1budget`.
  About 7 minutes; the plasma extended plates (about 20 s each) and the two
  plasma sheets (about 70 s each) dominate.
* Next: Phase 1, the decorrelation-time source or bound and the L2 fading
  model (gate G1), then Phase 2, the drift-mode sequence path.

## WP5 Phase 1 (2026-11-24 evening): decorrelation time and fading

* `docs/decorrelation.md`; project note `wp5-p1-decorrelation`.  No orbit
  measurement exists; ground-based values (ITU-R F.1487, Digisonde practice,
  Ganguly 2001) bound the intrinsic term only.
* G1 (Andrew): ground literature at low confidence and high variance; diffuse
  decorrelation tied to the facet (`tau_d / W = 0.79 D / L`); intrinsic times
  0.1, 0.5, 1.0 s; the route-one 0.1 s case moved to the extended reflector.
* `l2_instrument/fading.py`: Rician two-time-scale gain (specular part a phase wander by default, so K holds for any intrinsic time; `specular="rayleigh"` for the route-one gain), exp or gauss ACF exact at any time,
  per-echo log-normal intrinsic time, halo-derived diffuse time and coherence
  bandwidth, per-echo spawned streams.  Standalone until Phase 2.
* `route_one`: `FADING_CASES` is none and 0.5 s (seed codes kept, numbers
  unchanged); `EXTENDED_FADING_CASES` holds 0.1 s.
* `wp5_campaign`: group `E_budget`, the two budget points a facet can reach in
  the reference profile (5 and 10 MHz at 500 km); the other four are listed
  as unreachable.  `wp5_baseline` re-run (236bb86, from 6a5f9f2, `--strict`):
  117 scenes reproduced, the 115 Phase 0 rows and the budget unchanged,
  E_budget 2/2 at 0.016 m/s.
* Next: Phase 2, the drift-mode sequence path, which wires `fading` in.

## WP5 Phase 2 (2026-11-25): the drift-mode sequence path

* `docs/drift-sequence.md`; project note `wp5-p2-sequence-path`.  Decisions
  (Andrew, 2026-11-25): `EchoSet.phase_delay_s`; the orbit and the L0 profile
  now; the specular fading gain, oscillator continuity and the layer echo in
  this phase; a 150 kHz chirp at 45 Hz (`sweep_law.drift_mode`).
* `Pipeline.run_sequence`: L0 once, L1 per sweep at the sweep's centre, L2
  `receive_sequence` (Doppler inside each sweep, the first group delay
  blanked, one impairment draw over the sequence, specular fading per facet).
  L2 carrier phase follows the phase path (D8 condition 1).
* L1 tier `facet` (`l1_propagation/facet.py`): the plate model per sweep on a
  straight track or on the orbit by Earth flattening, in an exponential or a
  tabulated L0 topside; the layer echo on request.
  `plate_echo.echo_track` is the general entry; `echo_free_space` and
  `echo_stratified` are unedited and reproduced exactly.
* Exit met: noise-free L2 phase against `echo_stratified` to 1.3e-9 rad over
  the window (0.01 rad asked) and the Doppler to 1e-3; smoke and
  `wp4_campaign` unchanged.  The orbit + L0 chain agrees with the tier itself
  to 6.6e-9 rad (a self-consistency check; no independent spherical reference).
* For Phase 3: on the orbit the range curvature near overflight is
  `v^2 r_f / r_p / D`, 4.3 % under the flat value at 300 km below a
  600 km orbit (about 1.1 m/s on 25 m/s), so route one needs
  `v_eff^2 = v^2 r_f / r_p`; the blanked first
  `tau_g` costs 15 to 30 % of each record.  Not done: `sheet_echo`, the
  diffuse halo, any study driving `run_sequence`.
* Next: Phase 3, the drift tier and the Cramer-Rao harness.

## WP5 Phase 3 (2026-11-25, part 1): curved Earth, the sheet wrap, the drift tier

* `docs/drift-tier.md`; project notes `wp5-p3-decisions-2026-11-25`,
  `wp5-p3-part1`.  Decisions (a) to (h) of 2026-11-25 09:01.
* `l1_propagation/sphere.py`: each ray flattened in its own plane through
  the Earth's centre; held to Bouguer's law in radius to 1e-9 m.
* Facet tier: 3-D plates (default on the orbit; `plate_model`), sheets
  (`kind: sheet`; straight track = `sheet_echo.echo_track`, orbit = the 3-D
  engine with IGRF dip and declination).  Sheet exit met: 2.0e-5 rad, Doppler
  3.4e-7 (at 4208661); `field_aligned_verify` unchanged.
* L2: an echo set's geometry arrives half a delay later
  (`geometry_delay_fraction` 0.5 in a sequence); Phase 2's convention biased
  route one by `v^2 / (mu_sat c)`, 0.2 m/s.
* L3 tier `drift` (`l3/drift`): route one at signal level, first delay gated
  out, `f_eff` conversion, platform Doppler removed, explicit N(h) prior.
  Noise-free: 25 -> 24.997 (straight), 50 -> 49.996 (orbit, 3-D plate),
  -60 -> -60.35 (facet 2 km along track: mostly the envelope fit's
  absolute tolerances, fixed in part 2: -59.90 / -59.83 at HEAD), extended
  declined.
  **Noise-free only**: with noise the detector (spectrum summed over the
  whole sequence, 13 dB) finds nothing at 10 dB per sweep; open for the study.
* Correction to Phase 2: noise-free, to first order, `v_eff` does not move the
  drift (0.0004 m/s; 0.06 to 0.28 m/s with an asymmetric envelope); the
  curvature moves the slope ratio by `(v_eff/v)^2` (0.956).
* Next: the diffuse halo (scatterer field), the Cramer-Rao harness, the study
  `wp5_route_one` with the SNR grid, Poole's comparator.

## WP5 Phase 3 (2026-11-25 night, part 2): prior, detector, bound

* Committed c300fdd..60365b8; on origin/master (confirmed 2026-11-27).  `docs/drift-tier.md`
  sections 6 to 9; project note `wp5-p3-part2`.  Decisions (i) to (l), Option A.
* Fix: the drift tier scales the field to a unit peak before route one's
  envelope fit (scipy's absolute tolerances stopped it early on
  physical-scale powers, about -66 dB: `t_pk` spread 66 ms at 20 dB instead
  of about 10).
* Noise (beat given, 200 seeds, straight track): 20 dB 200/200 formed, drift
  sd 1.31 m/s; 10 dB 92/200 formed, sd 6.67 m/s.  The 10 dB declines are
  route one's free-slope test, on an ideal field too (87/200); no earlier
  study ran it with noise, so the budget's 10 dB point overstates
  availability by about two (measured at one point: 5 MHz, 530 km, 2 km
  facet, +/-2.5 W).  The estimator is pinned and unchanged.  The same fix
  explains most of part 1's -60.35 (now -59.90 / -59.83); `c3_verify`'s
  fields are near unit scale and unaffected.
* Onboard prior (Option A): `onboard_prior_from_ionogram` (the tracker's
  f_N, an assumed scale height; NO_PRIOR when invalid) and
  `Pipeline.run_prior_ionogram` (the scene's own L0; additive).  Chapman orbit
  scene: f_N 1.8203 against 1.8224 MHz, drift 49.992 against 49.982 (truth
  prior); the assumed scale height moves the depth (20 km), not the drift.
  Over 20 seeds: NO_PRIOR 0/20; f_N -0.11 % (clutter off) and +2.38 %
  (clutter on, a resonance) in every seed, a systematic error, not noise.
  Fixed `effective_speed_ratio`'s bracket (108 km returned for 170).
* Detector: windowed (21 sweeps), threshold for 1e-4 per cell from the
  negative-beat half; finds the facet at 10 dB in every trial, but at 1e-4
  per cell nearly every sequence also declares a false echo (1 in 13 forms a
  false drift); 1e-7 per cell removes them at no cost at 10 dB.  Default
  kept at 1e-4 as decided; the choice is Andrew's.  Range gate from the
  trace (`range_gate_from_ionogram`, larger of the O and X columns).
* Bound: `novel/drift_crb.py` (deterministic and Rayleigh-fading Fisher
  matrices; anchors 0.42 m/s and singular, Stein's FDOA bound per Moss 2014).
  Straight scene: 0.735 m/s at 20 dB, 2.33 at 10 dB; the estimator's RMS
  error is 1.6 and 2.8 times the bound (plan exit: within 20 % of it).
* Next: the halo, the study `wp5_route_one`, Poole's comparator; Andrew to
  decide the per-cell P_fa and whether route one's slope test changes.

## WP5 Phase 3, part 3 (2026-11-26): the false-alarm rate per sequence

* Branch `claude/wp5-p3-detector-rate`.  Decisions (m), (q) and (u)
  (Dissertation `Plans\WP5_plan.md`, P3, "Plan for part 3").
* Detector: `pfa_per_sequence` 1e-2 by default, divided over the (window,
  range bin) cells of each sequence (7.5e-8 per cell on the 209-sweep scene);
  `pfa` per cell still overrides it; the report carries both.
* 200 noise-only sequences: 2 with any cell over (L = 21), no false drift
  (15 at 1e-4 per cell).  10 dB per sweep: the facet found in 200/200, 2
  false echoes, 80 drifts on the facet and 1 on a false echo; the facet's
  drifts are identical to those at 1e-4 per cell on the same seeds.  7 dB:
  47, 50, 50 of 50 found at L = 11, 21, 41.  `docs/drift-tier.md` section 10.
* Fast suite 1067 passed (about 12 minutes on this machine), ruff clean,
  smoke unchanged before and after.
* No number cited in the TV or dissertation changes; the TV sentence on why
  the rate differs from WP4's is due at P3's TV update (section 8.20).
* Next: the joint estimator (`estimator="joint"` and `"joint_fading"`,
  timebox 8 January 2027), then the power-cut experiment, the halo, Poole's
  comparator and the study `wp5_route_one`.

## WP5 Phase 3, part 3 step 2 (2026-11-26): the joint estimators

* Branch `claude/wp5-p3-joint-estimator`, from step 1.  Decisions (o), (r),
  (s), (t) and the fading-aware fit (Andrew, 2026-11-27).
  `docs/drift-tier.md` section 11.
* New `l3_signal/novel/route_one_joint.py`: a deterministic fit (`joint`) and
  a Rayleigh-fading maximum-likelihood fit with `tau_c` fitted
  (`joint_fading`), both behind the pinned free-slope test.  Tier option
  `estimator`, **default still `pinned`**; `drift_pinned_mps` kept;
  `drift_sigma_mps` and `tau_c_s` formed.
* Two things needed beyond the prototype: the fading fit's noise is floored at
  3e-3 of the peak, the `sinc^2` model's own error (told the noise is 1e-6,
  45 of 56 noise-free facets miss by over 1 m/s); and the search starts from
  six cells at each width (from the six best cells overall it missed the
  lowest likelihood in 13 to 46 % of echoes).
* Noise-free within 0.43 m/s (straight and orbit); 56/56 and the tilt bias to
  1 % on the `c3_verify` fields.
* RMS over the bound, 200 seeds, without fading (20, 15, 10 dB): pinned 1.60,
  1.79, 2.46; joint 0.96, 0.97, 0.91; fading 1.12, 1.06, 0.98.  Under fading
  (20 and 10 dB at 0.5 s; 20 and 10 dB at 2 s): pinned 5.2, 3.8, 8.8, 4.3;
  joint 6.5, 6.5, 9.5, 4.2; fading 1.08, 1.78, 1.37, 2.90, with 0.8, 5.6, 4.0,
  8.2 % outliers and 1.03 to 1.25 without them.  Signal level without fading:
  joint 1.10 and 1.11 at 20 and 10 dB.
* `joint`'s sigma is 0.05 to 0.11 of its spread under fading: not to be quoted
  there.  About 1.5 s per echo for `joint_fading`, 10 ms for `joint`.
* Fast suite 1092 passed (about 14 minutes), the two full fading grids in the
  slow suite passed, ruff clean, smoke unchanged.
* The pinned estimator is unchanged, so `wp5_baseline` stands; no number cited
  in the TV or dissertation changes.
* Default (Andrew, 2026-11-27, amending decision (t)): leave `pinned` for
  now; switch to `joint_fading`, not `joint`, when appropriate.
* For Andrew, open: whether the explained fading shortfall is accepted; and
  `tau_c_s` here is the fitted exponential time, not `docs/drift-mode.md`'s
  lag-one definition.
* Next: the power-cut experiment (n), the halo, Poole's comparator, the study
  `wp5_route_one`.

## WP5 Phase 3, part 3 step 2b (2026-11-26): the fading outliers

* Branch `claude/wp5-p3-fading-outliers`, from master c152931 (steps 1 and 2
  merged).  Andrew, 2026-11-27: keep working on the fading shortfall inside the
  8 January timebox (2B); the 3e-3 noise floor accepted; the fitted time
  renamed `tau_c_fit_s` (flag `TAU_FIT_LOWER_BOUND`), `tau_c_s` kept for the
  `docs/drift-mode.md` definition, not formed yet.
* Ambiguity check: the fading fit reports `nll_margin` and `drift_alt_mps`;
  tier option `ambiguity_margin`, **off by default**, declines the echo.
  `docs/drift-tier.md` section 12.
* Pooled over 523, 305, 719, 305 echoes (20 and 10 dB at 0.5 s; 20 and 10 dB at
  2 s), margin off: outliers 1.1, 5.6, 3.3, 6.6 %, RMS 1.16, 1.58, 1.36, 2.11
  of the bound.  Margin 2: kept 100, 91, 99, 94 %; outliers 1.1, 2.9, 3.4,
  5.6 %; RMS 1.16, 1.33, 1.35, 1.61.  No cost without fading up to margin 8.
  **No margin brings all four within 20 %.**
* A fade hiding the envelope does go with the misses at tau_c 2 s (AUC 0.78 to
  0.83, from the true gain), weakly at 0.5 s; measures formed from the data do
  not find it.  The stated sigma ranks outliers well at 20 dB (AUC 0.86, 0.87).
  Averaging the likelihood's maxima gains little.
* Decided (Andrew, 2026-11-27): the check stays off (so availability stays
  equal across the estimators, decision (s)); the outlier work stops for now;
  M3 reports the outlier fraction beside the RMS; the combined margin and sigma
  flag is left untested.
* Fast suite 1102 passed, the two slow fading grids passed, ruff clean, smoke
  unchanged.
* The default estimator stays `pinned`; no number cited in the TV or
  dissertation changes.
* Next: the power-cut experiment (n), the halo, Poole's comparator, the study
  `wp5_route_one`.

## WP5 Phase 3, part 3 steps 3 and 4 (2026-11-27): power cut, halo, Poole

* Branch `claude/wp5-p3-halo-comparator`, from master ef29d24 (PR #4 merged).
  Decisions (n) and (e); the halo as option A and Poole as option A (velocity
  only), Andrew, 2026-11-27.
* Power cut (`docs/drift-tier.md` section 13, `scripts/wp5_power_cut.py`):
  `estimate_route_one(main_lobe=...)`, default 0.1 unchanged.  0.25 forms more
  drifts only at 10 dB without fading (40 to 53 %), fewer everywhere else; 0.5
  far fewer; both pass a 4 km plate as compact.  Nothing adopted.
* Halo (section 14): plate `halo: {k, ratio, n_points, seed}`, scatterers as
  `halo:<id>` rows; K the facet over the halo's mean power at overflight.
  48 draws per case through the chain: at K 10 the noise-free drift error goes
  from 0.1 to 15-17 m/s (pinned, joint), at K 3 to 33-63; `joint_fading` is the
  worst with the facet present, except at 4 lengths.
* **Finding for Andrew: with no facet (K = 0), which should be declined, route
  one forms a drift on half the draws, wrong by about 100 m/s.**  A random
  patch sweeps at the ephemeris slope, so the free-slope test cannot reject it
  (shown without L1 or L2; pinned by a test).  This contradicts the expectation
  in `docs/decorrelation.md`; the TV's text on a random halo as an extended
  reflector should be read with it (TV not edited).
* Poole (section 15, `poole_range.velocity_series`, `scripts/wp5_poole.py`):
  pair by pair the same observable as the tier's phase steps (0.3 m/s apart
  without noise); its lobe mean is 1.3-1.7 m/s at 20 dB against 0.7-0.9 for
  `joint`.  The range half cannot be formed from a one-frequency sequence.
* Fast suite 1139 passed (about 15 minutes), ruff clean, smoke unchanged.
* No number cited in the TV or dissertation changes; `wp5_baseline` untouched
  (route one's default is unchanged).
* For Andrew: how route one should treat a diffuse or rough reflector (before
  or in the study), and which halo cases the study carries.
* Next: the study `wp5_route_one` (its run time to be put to Andrew first).

## WP5 Phase 3, rough patches: a short look (2026-11-27)

* Andrew chose option (b): a bounded, report-only look before the study.
  `scripts/wp5_rough_test.py`, `docs/drift-tier.md` section 16.  Nothing in
  route one or the tier changed.
* Thresholds set to reject 5 % of unfaded reflectors.  The joint fit's residual
  catches the patches but also 50 to 100 % of faded reflectors.  At 20 dB a
  contrast below 56 catches all patches (K = 0) and loses 5 to 22 % of faded
  reflectors; at 10 dB the power far from the peak catches 85 % and loses 16 %.
* Thresholds depend on SNR, and both rely on the halo staying on through the
  pass (option A).
* The study's halo cases: none, K = 10 and 3 at 2 lengths, and K = 0 (put to
  Andrew 2026-11-27, no objection).
* Decided (Andrew, 2026-11-27): record only.  `wp5_route_one` records contrast
  and far-from-peak power for every echo and reports availability and error
  against halo strength; the rejection rule is chosen afterwards.
* The study is not started: its run time goes to Andrew first.

## WP5 Phase 3, the study `wp5_route_one` (2026-11-27)

* Branch `claude/wp5-p3-study`, from master 53859b9 (PR #5 merged).  Andrew:
  "full study now", the approved layout with 100 seeds, no instrument speed-up.
* `src/tss/studies/wp5_route_one.py`: five blocks (`docs/drift-tier.md`
  section 17), checkpointed per job and resumable; one row per echo and
  estimator.  The tier records `far_power_over_noise`.
* The reference run goes to `data/reference/wp5_route_one`, committed
  separately.
* Halo lengths (plan decision (dd), "two halo lengths"): 2 and 4 facet lengths
  at K = 10 and 3, the patch alone at 2; 600 halo jobs.  Andrew chose 4 as the
  second length ("go 4").  Estimated run about 12 h on 14 processes.

