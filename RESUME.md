# tss-sim build state — 2026-09-21

Delivered to `C:\Users\Andrew\PhD\tss-sim`.  170 files: 92 Python, 24 MATLAB,
15 Markdown, 11 YAML.

Status: `ruff check` clean, 239 tests pass, `tss run smoke` completes the full
L0 to L5 pipeline and produces a scored assessment, and all eight studies run or
refuse for a stated reason.

## First thing to do

This is not a git repository yet.  There is no shell on your machine from this
session, so the three commands are yours to run:

    cd C:\Users\Andrew\PhD\tss-sim
    git init && git add . && git commit -m "Initial import of the M&S simulator"

`.gitignore`, `.gitattributes`, `LICENSE` (BSD-3), `CITATION.cff` and a GitHub
Actions workflow are already in place, so `git remote add origin ...` and push is
all that stands between this and GitHub.

## What works today

- **Core**: `interfaces.py` (the L0–L5 contracts), `registry.py`, `config.py`,
  `requirements.py`, `provenance.py`, `pipeline.py`, `cli.py`, `constants.py`,
  `units.py`
- **L0** tiers: `chapman`, `analytic`, `iri2020`, `iri_wrapper`, `matlab_iri`
- **L1** tiers: `analytic` (real: absorption integral, magnetoionic splitting,
  Doppler forward model), `raytrace` (stub, PHaRLAP via the MATLAB bridge)
- **L2**: the migrated link budget whole, plus antenna, noise, waveform,
  receiver, oscillator; `leakage.py` is the critical-path stub
- **L3**: `baseline` tier real (de-chirp, transform, ionogram, trace extract);
  `motion_separated` registered and stubbed so the C3 campaign validates
- **L4**: `lamination` real, 162 m height error at 80 laminae against the
  closed-form exponential topside
- **L5**: `standard` real, scoring against `config/requirements/requirements.yaml`
- **Data**: SPDF ISIS, wspr.live, data.gov absorption, local model discovery
- **Studies**: link_budget, antenna_scale, sky_noise, min_frequency,
  isis_validation, transcription all run; doppler_separation and trade_space are
  documented stubs
- **Docs**: README, CONTRIBUTING, architecture, work-packages,
  requirements-traceability, data-sources, matlab-python-bridge, verification,
  and five ADRs

## Remaining, in priority order

1. **MATLAB gaps** — `matlab/bridge/` needs `tss_yaml_lite.m` and
   `tss_call_python.m`; `matlab/studies/` and `matlab/tests/` are empty;
   `matlab/README.md` is missing.  `matlab/tss/` (19 files) and
   `matlab/bridge/` (4 files) are done.
2. **`scripts/`** — only `fetch_data.py` exists.  Needs `run_study.py`,
   `make_figures.py`, `check_requirements.py` and `import_srm.py` (the last
   imports the remaining Rev E requirements from the .docx and .xlsx).
3. **Requirement set is a partial import** — 17 of the claim-bearing entries from
   Table 3.1.  The full Rev E set is 27 measurement, 26 instrument and 9 mission
   requirements.
4. **Tests** — `tests/contract/`, `tests/replication/` and `tests/archive/` are
   empty directories.  All 239 current tests are unit tests.
5. **WP4 then WP5** — the actual research.  `doppler_separation` is the study
   that resolves C3 and it is the one thing with no evidence of tractability
   behind it.

## Four things needing your judgement, not more code

1. **The 5 m antenna does not hold sky dominance at the band edge.**  (Updated
   2026-10-05.)  The short-dipole reactance in `antenna.py` was corrected to the
   Balanis form; the model inherited from `small_antenna.py` used twice the
   terminal capacitance and was 3 dB optimistic.  At radius 5 mm and coil Q 150
   the 5 m margin at 0.33 MHz is now +0.5 dB against 2.6e7 K and -2.5 dB against
   1.3e7 K, and about +7 dB at 0.6 MHz.  The old *0.5 dB against +10 dB*
   disagreement was one model quoted at different frequencies and skies.  What
   needs judgement is whether 5 m stays a candidate or 10 m becomes the floor.
2. **Possible 20 dB double count in `config/instrument/baseline_6u.yaml`.**
   `transmit.efficiency_db: -20.0` is commented as matching-network dominated,
   but a matching-coil loss is already inside
   `antenna.radiation_efficiency_db: -20.0`.  The baseline tier applies both
   because that is what the config says, and flags it in `describe()`.
3. **INS-01 cannot be scored as written.**  Its `sense: range` compares one
   achieved value against `[lower, upper]`, and its goal band `[0.26, 20] MHz`
   contains its threshold band `[0.33, 20] MHz`, so a sounder whose lowest
   usable frequency was 5 MHz would score `goal`.  L5 deliberately leaves it
   `not_evaluated` and records why.  The fix is in the requirement set: split it
   into a lowest-usable-frequency with `sense: max` and a highest with
   `sense: min`.
4. **Open item 12 is load-bearing for INS-26.**  Whether the range-correlation
   credit applies is unsettled, so `oscillator.py` scores the uncredited number
   and reports the credited one as a note only.  For the baseline the credit
   computes to 38 dB, well outside the 24–48x range recorded in the Technical
   Volume, which is itself a reason to settle the offset convention first.

## Two corrections already applied

- The "37 unit tests" on the link budget was wrong.  The author's
  `test_link_budget.py` has 19 test functions; the port has 21 (two
  migration-specific checks were added and are labelled as such).  Every
  occurrence in code and docs now says 21.
- A subagent wrote "38,953 ISIS-1 profiles" into `spdf_isis.py` and
  `docs/data-sources.md` without access to the file.  The count is now not
  asserted; `read_ne_ascii` reports it at read time.
- `tests/fixtures/isis1_ne_sample.asc` is **synthetic**, hand-built from the
  field widths in `crc_read.py` because the real 8.4 MB file was never staged.
  It round-trips through the unmodified `crc_read.read_crc`.  Replacing it with
  the first ~40 lines of the real file needs no test changes.
