# Clutter (WP4 Phase 1)

Status: implemented 2026-10-30.  Module `src/tss/l1_propagation/clutter.py`,
with L2 rendering of carriers in `src/tss/l2_instrument/interference.py`.  Plan
and decisions: project note `claude/wp4-plan.md` (decision D3, and the four P1
decisions of 2026-10-30).

## Why

Until Phase 1 every tss-sim ionogram held one clean O trace and one clean X
trace.  A trace extractor tested on that is tested against an ionosphere that
does not exist.  Alouette and ISIS topside ionograms carry other returns, and
the generator adds the ones that dominate them.

## The terms

| Term | What it is | How it is made | Flag |
|---|---|---|---|
| Ground echo | Above each mode's penetration frequency, the wave goes through the F layer to the Earth and back | Group path integrated through the profile (magnetoionic where the field is known) plus the vacuum below the grid; spreading over twice the group path; F-region loss integrated through the profile; D-region loss A/(f + f_L)^2; surface reflection | `clutter:ground` |
| Plasma resonances | Ringing of the plasma at the spacecraft at f_N, f_UH and n f_H, excited by the transmitter | Rows every 3 km of range in the step containing the resonance, power falling as exp(-tau/tau_d), to 30 dB down or 1000 km | `clutter:resonance:<name>` |
| Spread F | Oblique returns from field-aligned irregularities beyond the vertical echo | 12 rows per echo at offsets drawn from an exponential (mean 30 km, truncated at 100 km), total 6 dB below the echo; also up to 5 % above each mode's penetration frequency | `clutter:spread_f` |
| Broadcast interference | HF broadcast carriers, which reach orbit only above foF2 | Three carriers per ITU broadcast band on a 5 kHz raster, -70 dBm isotropic, 6 dB log-normal spread; rendered by L2 at the receiver input | (carrier list on `EchoSet.interference`) |
| Compression image | The image WP3's front end produces when the leakage compresses it | Not generated: L2 produces it (`receive_chain.py`), for example with the canceller off. Listed so the extractor is tested against it | — |

The second-hop echo of the original plan is left out.  A path of twice the group
range needs a reflector at the spacecraft, which is negligible, so it is not
what clutters a topside ionogram (decision of 2026-10-30).

Strengths are declared, representative values, not measurements.  The module
docstring says where each comes from, and `config/default.yaml` carries them.

## How it enters the pipeline

**Echo-like terms (ground echo, resonances, spread F).**
- They are extra `EchoSet` rows, on the sounding-frequency grid the echo set
  already uses, flagged `clutter:<kind>`.
- L2 synthesises them exactly as it synthesises an echo: stepped-sweep phase,
  step transition and Doppler.
- A clutter row never joins a mode's in-dwell dispersion trace.
- Clutter rows are synthesised without in-dwell dispersion.  This matters for
  the ground echo just above penetration, and it is a stated simplification.

**Carriers.**
- Carriers are not echoes.  They ride on `EchoSet.interference` (`freq_hz`,
  `power_dbm_iso`, `phase_rad`).
- L2 adds them at the receiver input, after the echoes and before the leakage,
  noise and receive chain.
- A carrier is de-chirped to `sqrt(P) exp(j[phi_ref(t) - 2 pi f_c t - theta])`.
- It is band-limited to `|f_ref - f_c| < fs/2`, with a raised-cosine edge over
  the outer 10 %.
- Its power is the isotropic power plus the antenna gain at f_c (efficiency
  included) plus the polarization loss.

**Truth scoring.**
- `truth_echo_cells` excludes every clutter row.
- So INS-10's cells, and the trace measures of `wp4_campaign`, are the profile's
  echoes only.
- Clutter can still lower the SNR measured in those cells.  That is the point.

**The pipeline.**
- `Pipeline(..., clutter=...)` applies the generator between L1 and L2.
- The generator draws from its own stream: the same seed and configuration,
  spawn key `0xC1`.
- The default is no clutter.

**Campaigns.**
- `run_campaign` builds the generator from the configuration's `clutter`
  block.
- That block is **on** in `config/default.yaml`.
- Turn it off with `--set clutter.enabled=false`, or turn off any one term with
  its own `enabled`.

**Pinned off.**
- The smoke campaign is pinned off, because WP3's cited numbers come from it.
- Studies written before WP4 construct their own pipelines without a generator,
  so they reproduce unchanged.

**Climatology.**
- `climatology: true` (off by default) draws spread-F occurrence from latitude
  and local time.
- It scales the D-region loss with cos(chi)^0.75 of the solar zenith angle.
- It weakens the carriers by 10 dB between 06 and 18 LT.

## Verification

`tests/unit/test_clutter.py` has 16 tests.  The checks are:

- **Ground echo, thin ionosphere.** The echo meets the free-space mirror closed
  form to 0.5 dB, with its delay at 2h/c to 0.2 %.  This holds from 5 MHz up;
  at 2 MHz the X mode's group retardation is still 2 %.
- **Spread F.**
  - The power beside each echo is level_db below it.
  - The offsets lie within the extent, and their mean is the truncated
    exponential's.
  - The frequency spread above penetration is bounded.
- **Resonances.**
  - Each sits in the step that contains it.
  - Each decays as stated, to 0.5 dB.
  - The spacing and extent are as declared.
- **Carriers.**
  - Only in broadcast bands, above foF2, on the raster.
  - The de-chirped power and instantaneous frequency are exact.
  - The antenna gain and polarization loss are applied, to 0.5 dB.
- **L2 rendering.**
  - A clutter row produces samples identical to an unflagged echo row.
  - Carriers change the record only while the sweep passes them.
- **Scoring and reproducibility.**
  - Truth cells are unchanged by clutter.
  - The generator is seeded and does not modify its input.
- **Clutter off reproduces the Phase 0 baseline.** One run is checked in the
  test.  The whole 30-run campaign also matches `data/reference/wp4_baseline`
  in every column but run time.

## What it does to the WP3-era chain

`data/reference/wp4_p1_clutter`: `wp4_campaign`, clutter on, 6 profiles x 5
seeds, against the clutter-off baseline `data/reference/wp4_baseline`.

| | Baseline (off) | Clutter on |
|---|---|---|
| SCI-05 median (worst) | 0.176 (1.17) | 1.48 (36.9) |
| SCI-05 verdicts | 15 fail, 10 threshold, 5 goal | 28 fail, 2 threshold |
| SCI-07 | 14 fail, 16 threshold | 8 fail, 22 not evaluated |
| INS-10 median (worst) | 14.4 dB (-13.9) | 5.9 dB (-16.4) |
| INS-10 fails | 5 of 30 | 20 of 30 |
| O-trace completeness, median | 0.14 | 0.00 |
| Spurious extractions, median | 0 | 40 |

Attribution, with each term on alone (one seed; `eq_evening_max`,
`mid_day_min`, `high_night`):

**The ground echo** destroys the trace.
- The extractor keeps the longest connected trace, takes the ground echo above
  foF2, and rejects the ionospheric trace as disconnected.
- Completeness falls to 0, with 34 to 93 spurious extractions.
- INS-10 is unchanged.

**The resonances** lower INS-10.
- The f_N and f_UH resonances fall in the first echo columns, because echoes
  begin just above f_N at the spacecraft.
- The n f_H resonances fall in low-frequency echo columns.
- The floor L3 measures there rises: INS-10 falls 20.6 dB on
  `eq_evening_max`, 5.4 dB on `mid_day_min` and 11.5 dB on `high_night`.

**Interference** costs 14.4 dB of INS-10 on `eq_evening_max`.
- A carrier at 12.04 MHz is above foF2 (12.0 MHz) but below fxF2, so it sits in
  an X-echo column.

**Spread F** has little effect at the declared level.

These are findings for the phases that follow, not defects of the generator:
- **P4:** the extractor must reject the ground echo, for example by its
  consistency with the penetration frequency and its near-constant range.
- **P2 and gate G1:** INS-10 is a sensitivity requirement, against noise.
  Whether it is scored with clutter off, or against a floor that excludes
  declared clutter columns, is a decision for G1.  As run, it measures
  signal-to-(noise + clutter).

## Cost

Clutter adds about 2,300 rows to an A1 echo set: about 1,500 resonance rows,
600 spread-F rows and 150 ground-echo rows.  L2 synthesises each row in its own
step window, so a run takes about 7 s instead of 1.3 s.  If the P3 campaigns
need it faster, the resonance rows of one step can be synthesised as a single
filtered-noise block.
