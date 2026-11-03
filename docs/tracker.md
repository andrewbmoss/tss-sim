# The trace tracker (WP4 Phase 4)

Status: implemented 2026-10-31.  The module is `src/tss/l3_signal/tracker.py`;
it is the default extractor of the baseline L3 tier (`extractor: tracker`), and
the threshold-and-connect extractor stays as `extractor: legacy` (the smoke
campaign is pinned to it).  Plan and decisions: project notes
`claude/wp4-plan.md` and `claude/wp4-p4-tracker.md`.

## Decisions (Andrew, 2026-10-31)

| Question | Decision |
|---|---|
| Architecture | **Track-before-detect, O and X labelled by self-consistency** |
| What the tracker may know | **Onboard information only**: the spacecraft position (so the IGRF field) and the ionogram itself |
| Gate G4, blended O and X | **Build the options, decide on data** |
| P3's rate capture | **Fix it first**: a rate per detection |

## Step 0: a rate per detection (P3 hand-off)

P3's rate search kept one chirp rate per column, chosen by the column's
highest peak, so a weaker echo with a different slope was transformed at the
wrong rate (rate capture).  `ionogram.assign_rates` now keeps every rate's
transform and starts the column from the best rate; the local maxima of the
per-cell maximum over rates are then accepted, strongest first, when they lie
outside every stronger echo's smear wedge (`|dr| <= 3|dj|/4 + 2` cells) and
are not more than 30 dB under it within its skirt; each accepted echo's
neighbourhood takes its own rate.  `rate_hz_per_s2` is now per cell.

- A weak undispersed echo beside a strong one sliding 40 cells reads
  -20.0 dB (its true level) against -30.0 dB per column.
- The detector keeps 1e-4 per cell with 15 effective looks (0.22 per rate),
  calibrated on operational records (leakage and canceller on): 0.96 to 1.05
  of the stated P_fa on two sets of six records; 0.63 to 0.72 with leakage off.
- `rate_assign: column` keeps the P3 rule.

## The tracker

1. **Onboard priors.**  f_H at the spacecraft from the field model at its
   position.  The local plasma frequency f_N from the ionogram: the plasma
   resonances ring through a whole column from zero range, so the lowest such
   column that is not a gyroharmonic is f_N (checked against f_UH when both are
   seen); with no resonance, the lowest strong near-range echo that a trace
   follows.  The X cut-off at the spacecraft is
   `f_x = f_H/2 + sqrt(f_H^2/4 + f_N^2)`.
2. **The image.**  Each column's power over its noise floor, max-pooled onto a
   2 km grid to 1050 km, with the range dependence of the background divided
   out (median over columns), and a 12 dB shadow test: a cell more than 12 dB
   under the strongest cell in the 100 km nearer than it is noise to the
   tracker (spread F, sidelobes).
3. **Masks.**  Streak columns (resonances, broadcast carriers).  The ground
   echo: tracked from the top of the band down (its range rises towards
   foF2), accepted when it reaches the top 15 % of the band and is mostly over
   the threshold, started at the first point after which it is; every cell at
   or beyond it (30 km margin) is masked.  Its branches give foF2 and fxF2,
   which cap the O and X searches.
4. **Track-before-detect.**  A Viterbi search for the path that maximises the
   summed evidence `min(z, 1000) - 6` (`z` the SNR as a power ratio), range
   non-decreasing with frequency (4 km tolerance), rising at most
   `max(60 km, 0.35 r)` per column (150 km cap), up to 3 columns skipped at a
   cost of 1 each; the path's ends are trimmed to its maximum-sum run.  The
   bias of 6 makes a path through noise lose evidence: the search picks, in
   effect, the best of a few tens of cells per column.
5. **O and X by self-consistency.**  At one frequency the X echo is nearer
   than the O echo, so the two paths are found with that ordering, each
   started near its own cut-off.  The O path is inverted (the magnetoionic
   lamination, with the field model) and the X trace it implies is predicted
   (`predict_x_range_m`: the X reflection height solves
   `f_p^2 = f (f - f_H)`, the virtual range integrates the X group index with
   the square-root singularity taken out); the X path is gated on the
   prediction (`max(5 km, 0.03 r)` plus half its change per column), the O
   path kept out of it where X was not found, and the two iterated.
6. **Samples.**  Each path point is refined to its detector run: the peak
   nearest the path, the run split at a 6 dB valley in front of it, and the
   edge at the first cell within 15 dB of the peak (for an echo under about
   25 dB SNR, the detector's first cell as before).  A track-only point (no
   cell over the threshold) is kept and flagged.  Range sigma:
   `sqrt((cell / 2 sqrt(SNR))^2 + cell^2 / 12)`.
7. **Cusp model.**  The top of the O trace is fitted with a parabolic layer's
   virtual range near its critical frequency, `a + b x ln((1+x)/(1-x))/2`,
   `x = f / foF2`, dropping up to three trailing points that flatten instead of
   rising; O samples at or above the fitted foF2 are removed; outliers are
   flagged (`cusp_policy: flag`) or replaced by the fit (`replace`).
8. **Gate G4** (`blend_policy`): where a column's O and X samples are within
   two range cells, `reject` drops the O sample, `inflate` raises its sigma to
   the overlap, and `separate` fits two echoes to the column's power (after
   Paul 1983).

**Quality bits** per sample (`flags["sample_quality"]`): 1 track-only,
2 blended, 4 spread, 8 cusp, 16 compressed at the edge of the rate-search span
(N > 64), 32 labelling confidence under 0.5, 64 rejected.  The ionogram also
carries `range_sigma_m`, `label_confidence` and the tracker's report
(`flags["tracker"]`: priors, ground echo, X consistency, cusp fit, iterations).

**Labelling is a model, not a measurement.**  `mode_discrimination` stays
False and INS-08 still fails on a scalar receiver; WP11's vector receive is
what measures polarisation.

## Results (`wp4_campaign`, 30 runs; `data/reference/wp4_p4_clutter`, `wp4_p4_clutter_off`)

At 62076e6, with G4 set to reject + flag.

| | P0 baseline, clutter off | P3, clutter on | **P4, clutter on** | **P4, clutter off** |
|---|---|---|---|---|
| SCI-05 median (worst) | 0.176 (1.17) | 2.34 (22.6) | **0.063 (0.114)** | **0.031 (0.065)** |
| SCI-05 verdicts | 15 fail | 26 fail | **25 goal, 5 threshold** | **30 goal** |
| SCI-07 km median (worst) | 18.8 | 51.6 (106.6) | **7.9 (12.1)** | **5.1 (12.2)** |
| SCI-07 verdicts | 14 fail | 5 fail, 21 not evaluated | **24 goal, 6 threshold** | **20 goal, 10 threshold** |
| O completeness, median | 0.14 | 0.00 | **0.67** | **0.95** |
| Labelling error, pooled | (O on X in 51 % of O columns) | — | **5.7 %** | **0.4 %** |
| Lowest retrieved height above hmF2, median (worst) | — | — | 27 km (116) | 27 km (59) |
| INS-10 fails | 5 | 5 | 5 | 2 |

**The exit criterion is met.** The criterion is SCI-05 and SCI-07 on the O trace, with clutter on, at or better than the P0 baseline. Every run is at goal or threshold on both. The labelling error is reported above.

**Where the gain comes from.** It comes from no longer mixing X samples, the ground echo, and spread F into the O inversion. It does not come from scoring a shorter profile. SCI-05 and SCI-07 are scored only from the true peak up to the lowest height the retrieval reaches (`coverage_gap_m`). On the same runs, P3 and P4 reach comparably close to hmF2:
- mid-latitude noon: 22–24 km under both;
- high-latitude night: 37 km under P4, against 58 km under P3;
- mid-latitude night: 59 km under P4, against 27 km under P3. P4 drops the cusp samples here (see the open items).

**INS-10.**
- **Clutter on:** all 5 fails are mid-latitude night, as in P2 and P3. The echo just above f_N shares its column and its range with the f_N resonance.
- **Clutter off:** the 2 fails left are rate capture inside the X echo's true smear wedge, at 8.1 and 9.8 dB. A rate per detection clears the third (11.3 dB). The other two cannot be separated by rate assignment.

## Gate G4 (Andrew, 2026-10-31): reject + flag

Clutter on, 30 runs per variant:

| blend_policy, cusp_policy | SCI-05 median (worst) | SCI-07 km median (worst) | Labelling error |
|---|---|---|---|
| **reject, flag** | 0.064 (0.111) | 7.9 (12.1) | 6.2 % |
| inflate, flag | 0.064 (0.110) | 7.9 (13.5) | 6.8 % |
| separate, flag | 0.064 (0.110) | 7.9 (13.5) | 6.8 % |
| reject, replace | 0.063 (0.111) | 7.9 (12.5) | 6.2 % |

All four have the same medians. Topside O and X traces stay more than two cells apart except near their starts, so blended samples are rare. Reject has the best worst cases and is the simplest to defend.

This comparison was run at 1cbbcc3, before the wedge revert of 62076e6. The revert changes the rate search, not the G4 rule.

## Open items

- **Close traces under heavy spread F.** Mid-latitude noon at solar minimum and the smoke profile have labelling errors of 10 % and 12 %. Their O and X traces are 13–20 km apart, and the clutter model puts spread F behind every echo. The O samples there read between the two truths. The "X bias" column of `wp4_campaign` is those samples, not the X trace.
- **Coverage.**
  - **Mid-latitude night:** the retrieval stops 59 km above hmF2, because the cusp fit removes the top O samples. Look at the trailing-point rule.
  - **High-latitude night, one run with clutter on:** the retrieval stops 116 km above hmF2.
- **f_N from the trace start** can fail (high-latitude night with clutter off). The tracker then falls back to the bottom of the band. The trace is still found.
- **Gate G6** (draft INS-32): re-verify the echo-level numbers at P6 with the rate-per-detection default and O and X together.
- **Cost:** the tracker adds about 0.5–1 s per record (reported, not gated).
