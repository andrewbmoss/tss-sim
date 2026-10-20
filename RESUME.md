# tss-sim: where things stand

Updated 2026-10-17, at the start of the simulation build that follows the
proposal (Rev D).  The history of what was built and why is in the git log and
in the project notes; this file says what the code does today, what it does
not, and what comes next.

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
| L1 | analytic: vertical incidence, magnetoionic group index for O and X, absorption, Doppler | raytrace (WP2, PHaRLAP) | Z branch, lateral deviation, off-vertical echoes |
| L2 | baseline, built from `config/instrument/baseline_6u.yaml`: de-chirped baseband from the A1 stepped sweep at the INS-06 digitiser rate, in-dwell dispersion, sky noise, receiver gain, clipping, quantisation | leakage (WP3, critical path) | leakage, oscillator phase noise in the samples, Doppler in the samples, sweep nonlinearity, real vector channels |
| L3 | baseline: per-step transform, leading-edge trace, connected-trace filter, SNR measurement at given cells | motion_separated (WP5) | O/X labelling, robust extraction (CFAR, a search limit at the design range), Doppler separation |
| L4 | lamination (field-free), lamination_magnetoionic (Budden 12.13, the default) | - | uncertainty propagation, depletion crossings, reference-field error |
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

So at 800 km INS-10 fails once the topside is IRI-2020's.  Not yet traced
(which echo sets the minimum, and whether it is the denser topside raising the
plasma frequency at the spacecraft).  The PyIRI column reproduces the
pre-change baseline exactly.

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
  INS-10 -5 dB (provisional).  G2 next, after Phase 5.
* `tss_iri_profiles.m` now writes the extras with nine significant figures.
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
  the single-pair model the route studies replaced.
* No configuration sits at the 900 km governing case of INS-01a and MIS-01.
* `docs/work-packages.md` still says Rev E in places.

## Next, in order

1. WP0 closure: the band-edge harness at 1000-1400 km (SCI-02 goal); INS-13
   downlink and on-board processing requirements.
2. WP3, L2: the leakage model and the de-chirped residual against INS-24;
   phase noise and Doppler into the samples; a route for L2-owned numbers into
   L5.
3. WP4, L3: robust trace extraction, O/X labelling (INS-08), the Barry and
   Poole reference former.
4. Gate M2: SCI-05 and SCI-07 scored against truth with both modes present,
   over seeds.
5. Beside 2-4: the MATLAB bridge (`tss_yaml_lite.m`, `tss_call_python.m`), the
   PHaRLAP L1 tier, NeQuick 2 and IRTAM L0 tiers, the depletion generator,
   WP12's noise reference.
6. WP5: `motion_separated` as a real tier over a sequence of sweeps, then WP6,
   WP7, the V2-V4 campaign, and the goal-level packages.
