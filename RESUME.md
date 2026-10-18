# tss-sim: where things stand

Updated 2026-10-17, at the start of the simulation build that follows the
proposal (Rev D).  This file replaces the delivery note of 2026-09-21.  The
history of what was built and why is in the git log and in the project notes;
this file says what the code does today, what it does not, and what comes next.

## Baseline

* Tag `proposal-baseline` marks `f103fc8`, the code the proposal and the
  Technical Volume cite: 474 tests, `tss run smoke` 0 goal, 8 threshold,
  2 fail (INS-04, INS-08), 10 not evaluated.
* Run everything in the `tss` conda environment (`conda activate tss`).  Plain
  `python` in PowerShell is a stand-alone 3.13 without the dependencies.
* `pytest` runs the unit and contract suites; `ruff check src tests scripts` is
  clean; `tss verify` lists what verification has and has not run.

## What each layer is today

| Layer | Implemented tiers | Stubs | What it does not yet model |
|---|---|---|---|
| L0 | chapman, analytic; iri2020 / iri_wrapper / matlab_iri need external software | perturbations (WP10) | NeQuick 2, IRTAM as a tier, depletions, storm days |
| L1 | analytic: vertical incidence, magnetoionic group index for O and X, absorption, Doppler | raytrace (WP2, PHaRLAP) | Z branch, lateral deviation, off-vertical echoes |
| L2 | baseline: de-chirped baseband from the A1 stepped sweep, in-dwell dispersion, sky noise, receiver gain, clipping, quantisation | leakage (WP3, critical path) | leakage, oscillator phase noise in the samples, Doppler in the samples, sweep nonlinearity, real vector channels, a configured digitiser rate |
| L3 | baseline: per-step transform, leading-edge trace, connected-trace filter | motion_separated (WP5) | O/X labelling, robust extraction, Doppler separation |
| L4 | lamination (field-free), lamination_magnetoionic (Budden 12.13, the default) | - | uncertainty propagation, depletion crossings, reference-field error |
| L5 | standard: scores the Rev G set, reports what it does not score | feasibility (WP7) | NmF2/hmF2 errors, 1 sigma coverage, Kp classes |

The C3 work (route one, route two, the field-aligned sheet, the tilt map, the
arrival-angle options) runs as studies at parameter level, outside the L0-L5
pipeline.  The pipeline itself runs one sweep at one epoch.

## Requirements (step 1, 2026-10-17)

* `scripts/import_srm.py` transcribes the SRM workbook into
  `config/requirements/srm_catalog.yaml`: all 70 Rev G requirements
  (SCI 28, INS 28, MIS 10, SW 4).  Generated; never edit it by hand.
* `requirements.yaml` scores 22 of them and lists the other 48 under
  `not_scored` with a category and a reason (8 flowed_down, 9 pending,
  7 analysis, 4 bench, 15 mission, 4 software, 1 tbd).  `tss requirements`
  fails if any catalogue ID is in neither place.
* New scored entries: SCI-15 (in proposal Table 3.1, missing until now) and
  INS-03.  The L5 report now ends with every not-scored requirement and why.
* To take in a new SRM revision: rerun `import_srm.py`, review the diff of the
  catalogue, amend the scored entries.  A unit test fails if the catalogue and
  the workbook in `..\Dissertation` disagree.

## Contract tests (step 1, 2026-10-17)

`tests/contract/test_tier_contracts.py` holds every implemented tier to the
four rules of `tss/interfaces.py`, and fails if a tier is registered that it
does not cover.  It found three departures:

1. **Fixed:** an echo set with no rows made L2 raise a TypeError.
2. **Recorded (strict xfail):** `tss.platform.orbits` gives the velocity in
   Earth-fixed axes, but L1 and L5 read it as north-east-down.  In every
   campaign L1 computes a 7.5 km/s vertical closing speed for a circular orbit.
   Harmless today only because L2 ignores Doppler.
3. **Recorded (strict xfail):** L2 sets its sample rate from the deepest
   simulated echo, so an echo-free record is sampled at 2.5 Hz and L3 raises.
   The same rule makes INS-04's scored value 1.25 times the deepest echo: the
   INS-04 failure in every smoke run is an artefact, not an instrument result.

## Other findings of the 2026-10-17 code review

* `BaselineInstrument.from_config` is never called.  Campaigns build L2 with
  its defaults, so `config/instrument/baseline_6u.yaml` and the antenna fields
  of the configuration files (length, count, vector receive) have no effect.
  `leo_500km_small_antenna` runs with the 5 m antenna; the trade study's
  antenna axis would change nothing.  The "20 dB double count" question has
  never affected a pipeline run, because neither -20 dB is applied.
* The pipeline does not pass the platform to L5, so SCI-09 is never scored in
  a campaign; and L2-owned numbers (INS-11, INS-14, INS-23, INS-24, INS-26)
  have no route into L5.
* The oscillator is built but never applied to the samples, and
  `EchoSet.doppler_hz` is ignored by L2: the samples carry no phase noise and
  no Doppler.
* `config/campaigns/c3_doppler.yaml` still asks for coherent intervals of 1 and
  2 s and up/down sweep pairing, both withdrawn; `motion_separated` describes
  the single-pair model the route studies replaced.
* No configuration sits at the 900 km governing case of INS-01a and MIS-01.
* `docs/work-packages.md`, `docs/architecture.md` and
  `config/instrument/baseline_6u.yaml` still say Rev E in places.

## Next, in order

1. Step 1 remainder: push to GitHub once the remote exists.
2. WP0 closure: the band-edge harness at 1000-1400 km (SCI-02 goal); INS-13
   downlink and on-board processing requirements.
3. WP3, L2: fix the velocity frame and give L2 a configured digitiser rate
   (removes both xfails); call `from_config` so the instrument file and the
   configurations take effect; then the leakage model and the de-chirped
   residual against INS-24; phase noise and Doppler into the samples.
4. WP4, L3: robust trace extraction, O/X labelling (INS-08), the Barry and
   Poole reference former.
5. Gate M2: SCI-05 and SCI-07 scored against truth with both modes present,
   over seeds.
6. Beside 3-5: the MATLAB bridge (`tss_yaml_lite.m`, `tss_call_python.m`), the
   PHaRLAP L1 tier, NeQuick 2 and IRTAM L0 tiers, the depletion generator,
   WP12's noise reference.
7. WP5: `motion_separated` as a real tier over a sequence of sweeps, then WP6,
   WP7, the V2-V4 campaign, and the goal-level packages.
