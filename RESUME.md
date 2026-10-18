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
| L0 | chapman, analytic; iri2020 / iri_wrapper / matlab_iri need external software | perturbations (WP10) | NeQuick 2, IRTAM as a tier, depletions, storm days |
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
