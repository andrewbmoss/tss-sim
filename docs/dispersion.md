# In-dwell dispersion and the step transition (WP4 Phase 3)

Status: implemented 2026-10-31.  The treatments live in the baseline L3 tier,
`src/tss/l3_signal/ionogram.py` (`BaselineSignalProcessing`), with detector
support in `src/tss/l3_signal/detect.py` and the narrow step law in
`src/tss/l2_instrument/sweep_law.py`.  The comparison is the study
`wp4_dispersion` (`src/tss/studies/wp4_dispersion.py`).  Plan and decisions:
project notes `claude/wp4-plan.md` and `claude/wp4-p3-dispersion.md`.

## The problem

Within one dwell the group path of an echo changes with frequency, so its beat
frequency slides: an echo with slope dP'/df crosses

    N = (2 Δf² / c) · dP'/df

range cells in a step of width Δf.  Near the local plasma frequency at the
spacecraft and near each penetration cusp, dP'/df is large and one transform of
the whole dwell smears the echo over tens of cells.  The smear lowers the peak
(an SNR loss) and widens the echo (an early leading edge).  Separately, at a
step boundary the de-chirp reference changes rate before the echo does, and the
first round-trip time of every dwell carries a transition skirt.

The TV (section 5.3 and Ch. 2) asked WP4 to compare three dispersion
treatments and three step-transition treatments, choose one of each, and
propose a dispersion condition for Rev G.

## Decisions (Andrew, 2026-10-30 and 2026-10-31)

| Question | Decision |
|---|---|
| How to score | **Echo level, independent of the extractor, plus the chain**, clutter off |
| Narrower steps | **Only where the smear needs it** (2 to 8 MHz) |
| One treatment or a hybrid | **One treatment band-wide** for each of dispersion and the transition |
| The matched transition option | **The hybrid matched filter** (re-evaluate detections in steps whose rate changes) |
| **G2**, dispersion | **D2, the rate search** |
| **G3**, step transition | **Accept the skirt** (T0) |
| Rate-search span | **Widened to 64 cells**; echoes beyond it go to P4's cusp model |
| **G6**, the Rev G condition | **A performance condition on the processing** (draft below) |

## The treatments

All are options of the baseline L3 tier, set in the `l3` block of a campaign
(`dispersion:` and `transition:`) or by keyword.

**Dispersion.**

| Name | Setting | What it does |
|---|---|---|
| D0 | `dispersion: none` | One transform of the whole dwell (the WP3-era chain). |
| D1x2, D1x4 | `dispersion: blocks`, `dispersion_blocks: 2` or `4` | The dwell cut into sub-blocks, each transformed and added in power.  Each sub-block crosses fewer cells, at the cost of coherent gain and a coarser cell. |
| D2 | `dispersion: rate_search`, `rate_search_max_cells: 32` | De-chirp the block again against beat chirp rates γ = j/T², j = −2..32, about the block centre; keep, for the column, the rate whose transform has the highest peak.  An echo whose beat slides j cells is compressed back into one cell at the step-centre range.  The chosen rate gives dP'/df = γc/(2K²), which goes to WP6. |
| **D2w** (default) | `dispersion: rate_search` (span 64) | The same over j = −2..64. |
| D3 | sweep keys `narrow_lo_hz: 2e6`, `narrow_hi_hz: 8e6`, `narrow_cap_hz: 50e3` | A 50 kHz step cap between 2 and 8 MHz (150 kHz elsewhere), sized so that the B2 reference profile crosses about two cells at 3 MHz.  237 steps and a 6.35 s sweep instead of 159 and 4.26 s. |

A per-cell maximum over rates was tried before the per-column rate: the
mismatched rates spread an echo ahead of itself, the maximum kept that spread,
and leading edges read early.

**Step transition.**

| Name | Setting | What it does |
|---|---|---|
| **T0** (default) | `transition: accept` | Keep the whole dwell and accept the skirt. |
| T1 | `transition: gate` | Drop the first 2R_max/c (7 ms) of every dwell. |
| T2 | `transition: matched` | The hybrid matched filter: in the 45 of 159 steps whose rate changes, re-evaluate the cells around each detection against the exact transmitted beat history, `s_tau(t) = exp(j[phi(t) − phi(t − tau)])`, which carries the transition. |

## Keeping the stated P_fa

Each treatment changes the noise statistic the detector sees, and P2's
1e-4 per cell was kept for all of them.

- **Sub-blocks (D1).**  The column is a sum of m independent exponentials, so
  the threshold comes from the gamma law, integrated over the sampling
  distribution of the median reference (`detect.gamma_median_alpha`).  The
  first form, which treated the median as exact, gave 1.8e-4 for D1x4.
- **The rate search (D2).**  Choosing the best of several rates is a choice
  among looks: the detector takes a per-look P_fa of 1 − (1 − P_fa)^(1/L).  The
  rates are strongly correlated, so L is an effective count, calibrated on
  noise-only A1 records: about 7 for the 32-cell span and 12 to 14 for the
  64-cell span (default `0.2 × (span + 3)`, 13).  Counting every rate gave a
  measured 2.9e-5.
- **D1 weakness, recorded.**  With leakage on, D1x4 measures 1.6e-4.  The
  positive and negative halves differ by about 0.2 dB with the leakage
  residual present, and the gamma tail amplifies that.  D1 was not chosen.

Tests: `tests/unit/test_dispersion_treatments.py` (P_fa per treatment, the
compression of a sliding beat, D2w and D3 recovering the loss, sub-block
coarsening, gating, the matched filter removing the skirt, the default).

## The comparison (`data/reference/wp4_p3_dispersion`)

**Scenes.**  The six Chapman ionospheres of `wp4_campaign` and 12 ISIS-1
hand-scaled topside profiles (stratified by latitude, local time and
altitude), each completed below its lowest point by an alpha-Chapman peak.
18 scenes, O and X run separately, three seeds.

**Measures, per truth echo.**  The reference is the same record with
dispersion switched off in L2 (same noise draw) under D0, so a treatment's
loss includes what it costs an undispersed echo as well as what it recovers.
Loss, width within 3 dB, leading-edge error, detection, the column's cell
(INS-05), and for the rate search the relative error of its dP'/df.

**Echo level.**

| Treatment | Loss dB median (p90) | Width km median | Edge error km median | Detected | Loss near f_N / cusp dB |
|---|---|---|---|---|---|
| D0 | 4.1 (9.2) | 4.0 | −5.9 | 0.981 | 8.2 / 8.2 |
| D1x2 | 6.3 (9.3) | 4.0 | −8.2 | 0.987 | 8.1 / 8.5 |
| D1x4 | 8.4 (12.0) | 9.5 | −11.4 | 0.982 | 11.3 / 11.5 |
| D2 | 0.1 (5.4) | 1.5 | −2.8 | 0.989 | 1.1 / 3.3 |
| **D2w** | **0.1 (2.0)** | **1.5** | **−2.8** | **0.992** | **0.3 / 0.6** |
| D3 | 0.3 (4.4) | 4.5 | −5.3 | 0.975 | 5.0 / 5.6 |
| T1 | 3.2 (8.2) | 3.4 | −3.8 | 0.982 | 6.5 / 7.6 |
| T2 | 4.1 (9.2) | 4.0 | −6.0 | 0.981 | 8.2 / 8.2 |

**D2w by cells crossed, N.**

| N | Share of echoes | Loss median / p90 dB |
|---|---|---|
| ≤ 32 | | 0.08 / 1.04 |
| ≤ 64 | 94 % | 0.10 / 1.24 (p95 2.65) |
| 32 to 64 | | 0.41 / 2.13 |
| 64 to 128 | | 3.54 / 8.51 |
| > 128 | | 8.68 (median) |

For N ≤ 64: leading-edge |error| p90 6.6 km, dP'/df |error| median 6.5 %,
detected 0.996.

**The chain** (`wp4_campaign`, clutter off).

| Treatment | SCI-05 median | SCI-07 km median | INS-10 dB median (worst) | INS-05 worst km | Run s per record |
|---|---|---|---|---|---|
| D0 | 0.177 | 18.9 | 14.7 (−14.1) | 13.6 | 0.9 |
| D1x4 | 0.208 | 22.5 | 15.9 (9.8) | 54.5 | 1.0 |
| D2 | 0.147 | 17.4 | 20.2 (−4.6) | 13.6 | 3.3 |
| **D2w** | **0.146** | **17.4** | **20.2 (9.2)** | **13.6** | **5.4** |
| D3 | 0.172 | 18.2 | 17.9 (11.7) | 13.6 | 1.4 |
| T1 | 0.153 | 18.9 | 16.3 (−9.9) | 18.5 | 0.9 |
| T2 | 0.177 | 18.9 | 14.7 (−14.1) | 13.6 | 1.2 |

## Why these choices

- **D2w.**  The only treatment that recovers the loss band-wide without
  changing the sweep: 0.1 dB median, 2.0 dB p90, and 0.3 and 0.6 dB near f_N
  and the cusp.  It narrows the echo (1.5 km against 4.0) and moves the
  leading edge 3 km closer to truth.  Its dP'/df is a by-product for WP6.
- **D1** loses on every measure: the coherent gain given up costs more than
  the smear saved, and the coarser cell fails INS-05 (D1x4: 54.5 km).
- **D3** fixes 2 to 8 MHz only, costs 78 more steps and 2.1 s of sweep
  (INS-07; the sweep spans about 45 km along track instead of 30, against
  SCI-09), and leaves the low band and the cusp untreated.
- **T0 (G3).**  T1 coarsens the cell to 18.5 km and fails INS-05.  T2 is
  identical to T0 on every echo-level and chain measure: the skirt reaches
  the detector only for very strong echoes (its unit test needs an echo at
  −50 dB to see it), which P2's `test_sweep_law` finding was.  T2 removes the
  skirt in that case and stays available as an option.

## The cusp and the 64-cell span

The 32-cell span left a residual near the cusp (3.3 dB); that came from echoes
crossing more than 32 cells, not from curvature of the beat within the dwell.
Widening to 64 brings the cusp loss to 0.6 dB.  The 6 % of echoes with N > 64
(near f_N at the spacecraft and at the cusps) still lose 3.5 dB median or
more.  Decision: **they go to P4's cusp model**, which fits the trace near
penetration rather than thresholding it, and they are reported separately.

## The mid-latitude-night INS-10 fail was dispersion

P2 attributed the five mid-latitude-night INS-10 fails (about −14 dB) to the
antenna and sky noise at the bottom of the band.  That was wrong.  The worst
echo is the one just above f_N at the spacecraft (2.27 MHz, 3 km range), which
crosses so many cells under D0 that it reads −13 dB; D2w compresses it to
+18 dB.  With clutter off the profile's INS-10 median goes from −13.7 to
+9.7 dB, and the campaign's fails from 5 to 3.  The three left are the next
item.

## The rate is chosen per column: capture

D2 keeps one rate for each column, the rate of the column's highest peak.  A
weaker echo in the same column with a different slope is then transformed at
the wrong rate.

- **O and X, clutter off.**  Over the 30-run campaign (1050 scored echoes),
  D2w gains 4.2 dB median over D0 per echo; 3 % of echoes lose more than 1 dB
  against D0 and 0.9 % more than 3 dB (at most 4.1 dB).  All are O echoes in
  columns whose rate is set by a stronger echo (mid-latitude noon at solar
  minimum, 2.5 to 2.6 MHz; mid-latitude night, 3.34 MHz).  The echo-level
  study, which runs O and X separately, does not see this.
- **It sets the three INS-10 fails left with clutter off.**  All three are the
  mid-latitude-night O echo at 3.34 MHz and 459 km, which D2w reads at 8.1 to
  9.7 dB against 11.2 to 12.9 dB under D0, which passes.  Choosing the rate
  per detection should clear them.
- **Clutter on.**  Over the same 1050 echoes, clutter costs D2w 0.07 dB
  median; 1 % of echoes lose more than 3 dB against D2w with clutter off.  The
  serious case is the echo just above f_N at the spacecraft on mid-latitude
  night (2.27 MHz, 3 km): it shares its column and its range with the f_N
  resonance, which sets the rate, and it reads −9 to −13 dB against +17 to
  +18 dB with clutter off.  That is the profile's INS-10 fail with clutter on
  (−10.8 dB median, as in P2).  The other large losses (equatorial-evening X
  echoes at 12.04 MHz, 10 to 13 dB) come from the broadcast carrier in that
  column and hit D0 equally.
- **For P4.**  Choose the rate per detection, not per column: peaks in the
  range-rate plane, with the cells inside a stronger peak's smear wedge
  (|Δr| ≤ |Δj|/2 plus a guard) assigned to that peak, then recalibrate the
  P_fa.  That clears the clutter-off cases.  The f_N echo under its resonance
  overlaps it in range, so it separates only if the two differ enough in rate;
  otherwise it is a resonance-rejection problem.  The rate is also an echo
  feature for the tracker: the ground echo and the resonances have rates
  unlike the trace's.

## The chain on the new default (`data/reference/wp4_p3_clutter_off`, `wp4_p3_clutter`)

`wp4_campaign`, 30 runs, D2w and T0 with the P2 detector.

| | P2, clutter off | P3, clutter off | P2, clutter on | P3, clutter on |
|---|---|---|---|---|
| SCI-05 median (worst) | 0.177 | 0.145 (1.13) | 0.874 | 2.34 (22.6) |
| SCI-05 fails | 15 | 15 | 22 | 26 |
| SCI-07 km median | 18.9 | 17.4 | | 51.6 |
| INS-10 dB median (worst) | 14.7 (−14.1) | 20.2 (8.1) | 18.1 (−13.7) | 20.1 (−12.7) |
| INS-10 fails | 5 | 3 | 5 | 5 |
| O completeness median | 0.12 | 0.16 | 0.00 | 0.00 |
| Run s per record | 0.9 | 5.3 | 6.3 | 11.0 |

- **Clutter off:** every measure improves or holds.  The three INS-10 fails
  left (8.1 to 9.7 dB) are rate capture, above.
- **Clutter on:** the ground echo still takes the extractor (O completeness
  0.00 in both), so SCI-05 and SCI-07 are not meaningful until P4.  SCI-05 is
  worse than P2 on mid-latitude noon at solar minimum and on the smoke profile
  and better on high-latitude night; the extractor that produces it is
  replaced in P4, so the difference is not analysed further.

## Draft Rev G condition (G6)

> **INS-32 Dispersion within a dwell.**  The ionogram processing shall hold the
> dispersion loss to no more than 1 dB median and 3 dB at the 90th percentile
> for echoes that cross up to 64 range cells in one dwell, N = (2Δf²/c)·dP'/df,
> and shall report dP'/df for each detected echo.  Echoes crossing more than 64
> cells (near the local plasma frequency and the penetration cusps) are
> handled by the cusp model and reported separately.
> *Verification: simulation (`tss study wp4_dispersion`).*

Measured: 0.10 dB median and 1.24 dB p90 at N ≤ 64.  The condition goes into
Rev G at P6; the proposal is not changed.

## Cost

Processing budgets are outside the scope of this work (WP0 decision); the cost
is reported, not gated.  D2w transforms each block 68 times: the L3 run time
per record is about 6 times D0's (5.4 s against 0.9 s on the development
machine).  An FFT-based chirp-rate transform, or a coarse-to-fine rate search,
would cut it if it ever matters.

## Items passed on

- **P4:** the cusp model for N > 64; rate capture (rate per detection);
  the rate as a tracker feature; ground-echo and resonance rejection.
- **P6:** INS-32 into Rev G; TV section 5.3 (plan becomes result) and the Ch. 2
  step transition.
- **WP6:** dP'/df from the rate search, |error| median 6.5 % for N ≤ 64.
