# The detector (WP4 Phase 2)

Status: implemented 2026-10-30.  The module is `src/tss/l3_signal/detect.py`,
and it is the default cell test of the baseline L3 tier.  Plan and decisions:
project notes `claude/wp4-plan.md` and `claude/wp4-p2-detector.md`.

## Decisions (gate G1, 2026-10-30)

| Question | Decision |
|---|---|
| False-alarm probability | **1e-4 per range cell** |
| Noise estimate | **Echo-free negative-beat half**, with a local check against a strong return's sidelobes |
| Integration | **Single-cell CFAR here; track-before-detect in P4's tracker** |
| INS-10 | **SNR against the noise floor** (noise and carriers, not clutter); signal to noise plus clutter reported beside it |

## The method

**The noise reference.**
- The receiver samples complex baseband. After the de-chirp's orientation,
  every echo beats at a positive frequency, so the negative half of each
  block's transform carries noise and carriers and no echo or clutter.
- The transform now keeps that half (`RangeSpectrum.negative_power_db`).
- The reference is its median, an ordered statistic.
  - It excludes 4 cells around zero beat, where the leakage image sits.
  - It is robust to a few contaminated cells: returns folded from beyond the
    unambiguous range, and in hardware the I/Q-imbalance image of strong
    returns.
- The threshold is `alpha * median`. `alpha` comes from the ordered-statistic
  false-alarm law (Rohling 1983),
  `P_fa = prod_{i<k} (N - i)/(N - i + alpha)`.
  - `N` is the number of independent cells: the bins over the padding factor
    and the window's noise bandwidth, about 700 per column.
  - `k = N/2`.
- With `N` that large the test is within a few hundredths of a decibel of the
  known-noise (Neyman-Pearson) test. The threshold is 9.7 dB over the mean
  noise.

**The local check: sidelobe blanking.**
- A Hann window's first sidelobe is 31.5 dB down. A return more than about
  41 dB over the noise therefore has sidelobes that pass the global test: the
  ground echo, a resonance, a strong X echo.
- Any cell under a stronger detected peak's sidelobe envelope plus 6 dB is
  blanked.
- The envelope is the upper bound of the window's transform over every
  fractional bin offset, so it holds wherever the echo falls in its bin.
- A real return more than 6 dB above the skirt survives. An O echo beside a
  stronger X echo is kept.

**The optional ordered-statistic veto (off by default).**
- Built first as the local check: 16 cells a side, 2 guard cells, the 0.75
  statistic, at 1e-2.
- Echoes smeared over 50 to 100 km by dispersion within the dwell (near f_N and
  the cusp) fill their own reference window and were vetoed. On the
  equatorial-evening and mid-latitude-night profiles of `wp4_campaign` that
  lost whole traces.
- This is the self-masking of an extended target that an ordered-statistic
  CFAR is known for. The option is `local_pfa`.

**Search limits.**
- The search runs from 0 to INS-04's 1050 km.
- The WP3-era near-range exclusion is gone. L2's leakage high-pass and the
  canceller keep the residual below the threshold; the test below shows it.

**Leading edge.** The first run of at least 2 bins (`min_run`) over the
detector, as before. The tracker of Phase 4 replaces it.

**The fixed rule.** The WP3-era test, 10 dB over the median of the column, is
still available as `detector: fixed`. The smoke campaign is pinned to it so
that WP3's cited numbers reproduce.

## What is stated, and what was measured

All rates are at `P_fa` 1e-4 per range cell. "Per cell" is the marginal
probability that a noise-only bin is declared. Zero padding correlates
adjacent bins without changing it.

| Check | Result |
|---|---|
| Measured `P_fa`, white noise | 0.8e-4 (exit criterion: within a factor of 1.5 of nominal) |
| Measured `P_fa`, full L2 noise-only records with leakage, canceller and AM at the WP3 allocations | 1.25e-4 |
| False alarms below 50 km on the same records | none |
| Positive half against negative half, noise floor | agree to 0.0 dB (no leakage) and 0.2 dB (leakage) |
| SNR for Pd = 0.5, steady target | 9.46 dB measured, 9.40 dB theory (Marcum Q, known noise); 0.06 dB apart |
| SNR for Pd = 0.9 | 11.94 dB measured, 11.75 dB theory; 0.19 dB apart (exit criterion: 0.5 dB) |
| Pd at INS-10's 10 dB | 0.60 measured, 0.62 theory |
| False alarms per ionogram | about 92,000 independent cells within 1050 km; about 26 distinct over-threshold runs, and about 12 runs of the 2 bins the extractor needs, in about 7 % of columns |

The tests are in `tests/unit/test_detect.py` (14 tests). They also check:
- the noise reference is unmoved by a 40 dB echo and a clutter continuum, which
  raise the column median by more than 3 dB;
- the sidelobes of a 71 dB return are removed, and an O echo 6 dB under an X
  echo 12 cells away is kept;
- a range-extended echo is detected, and the optional veto would mask it;
- the search limit holds;
- the fallback, when there is no negative half, says so.

## INS-10 (G1)

`snr_at` measures against the detector's noise reference by default
(`floor="noise"`).

The pipeline also measures every truth echo against the column
(`floor="column"`, which clutter raises). L5 reports the minimum of that as
`snr_db_min_at_truth_echoes_vs_column`, signal to noise plus clutter. It is
reported, not scored.

The Rev G verification note for INS-10 will state the detector and its
`P_fa`. That goes in with the other Rev G amendments in P6.

## Effect on the chain

`wp4_campaign`, 6 profiles x 5 seeds:
- `data/reference/wp4_p2_clutter`: clutter on;
- `data/reference/wp4_p2_clutter_off`: clutter off.

Compared with:
- `data/reference/wp4_baseline`: P0, fixed rule, clutter off;
- `data/reference/wp4_p1_clutter`: P1, fixed rule, clutter on.

| | P0 (fixed, off) | P2 (CFAR, off) | P1 (fixed, on) | P2 (CFAR, on) |
|---|---|---|---|---|
| INS-10 median (worst) | 14.4 (-13.9) | 14.7 (-14.1) | 5.9 (-16.4) | **18.1 (-13.7)** |
| INS-10 fails | 5 | 5 | 20 | **5** |
| Signal to noise plus clutter, median | — | 14.4 | — | 5.9 |
| SCI-05 median | 0.176 | 0.177 | 1.48 | 0.874 |
| SCI-05 fails | 15 | 15 | 28 | 22 |
| O-trace completeness, median | 0.14 | 0.12 | 0.00 | 0.00 |
| Spurious extractions, median | 0 | 0 | 40 | 22 |

**Clutter off.** The detector is neutral on the trace, as it should be. Its
threshold is 9.7 dB over a floor that matches the column median to 0.2 dB, so
it is the old 10 dB rule with a stated `P_fa`.

**Clutter on.**
- INS-10 is noise-limited again, which is what G1 asked for. The 5 failures are
  the mid-latitude-night profile, which fails with or without clutter: about
  -14 dB at the bottom of the band.  (Corrected in P3, 2026-10-31: this was
  attributed to the antenna and sky noise; it is in-dwell dispersion of the
  echo just above f_N at the spacecraft, which the rate search recovers.  See
  `docs/dispersion.md`.)
- The trace is still lost to the ground echo. That is P4's work.

**An observation for P4.** With the optional veto on, the high-latitude-night
profile's O completeness rose from 0.00 to 0.83. Its X echoes are more spread in
range than its O echoes, the veto suppressed them, and the leading edge fell on
O. That is accidental mode discrimination and not a detector setting to adopt.
But echo extent is a candidate feature for the P4 tracker's O/X hypotheses.

## Comparison point (Poole and Evans 1985)

Poole and Evans declared an echo when at most 16 of the highest spectral values
of a 1024-point transform exceeded an amplitude threshold, a sparsity rule
against broadband interference. P5 (item A8) compares it with this detector.
