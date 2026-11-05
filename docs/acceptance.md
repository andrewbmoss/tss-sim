# FM-CW acceptance and the measured-ionosphere check (WP4 Phase 5)

Status: implemented 2026-11-03.  Studies `fmcw_acceptance` and `isis1_retrieval`;
module `src/tss/l3_signal/poole_range.py`.  Plan and decisions: project notes
`claude/wp4-plan.md` and `claude/wp4-p5-acceptance.md`.

## Decisions (Andrew, 2026-11-03; G5 2026-10-31)

| Question | Decision |
|---|---|
| Test level | **Waveform + L2/L3**: the chain's de-chirp, transform and detector at the papers' own sweep parameters, on synthetic echoes (stated delays, Doppler, dispersion, noise).  No ground vertical-incidence geometry is built. |
| Scope beyond Poole's items and Barry's table | The ISIS-1 measured check, A8 (the echo criterion) and A9 (the group-velocity method).  Barry's gating loss is not simulated. |
| ISIS-1 | Every usable record, clutter on and off; quality index 4 to 9, reported by quality. |
| Gate G5 tolerances | Arithmetic within 1 %; statistical items reproduced when the published value is inside the Monte Carlo's 95 % CI (A6 half-width at 20 dB at most 0.5 deg; measured 0.13); Barry items within 10 %.  Outcomes: reproduced, difference explained, not reproduced (open); hardware items listed as not simulated, with the reason. |

The papers: Poole (1985), *Advanced sounding 1. The FM-CW alternative*,
Radio Sci. 20(6) ("paper 1"); Poole and Evans (1985), *Advanced sounding 2.
First results from an advanced chirp ionosonde*, Radio Sci. 20(6) ("paper 2");
Barry (1971), *A low-power vertical-incidence ionosonde*, IEEE Trans. Geosci.
Electron. GE-9(2).  Their text is not in this repository.

## Poole's estimator (`poole_range`)

Three de-chirped cells, the first two at the same start frequency and the third
`df` higher, T apart.  On the chain's window and scaling:

* phase velocity from the cell 1 to 2 phase (eqs. 72-73), referred to the
  cells' centre frequency (a tapered transform's phase belongs to the window
  centroid; the paper's reference is the cell start);
* group range by stationary phase (eqs. 42, 70-71):
  `R' = c (N + phi_R / 2 pi) / (2 df)`, `N` from the beat frequency.

The algebra drops terms of order `V* T`: evaluated at the cells' centroids
`phi_R` carries `5 df V* T / c` beyond eq. 42, so `R'` reads `2.5 V* T` long
(+0.31 km at 250 m/s), and eq. 59's velocity bound moves by about 2 %.  This is
a reference estimator for WP4 and WP5 (the C3 comparator), not the chain's
range measure.

## `fmcw_acceptance`: the items

`data/reference/wp4_p5_acceptance`, 2000 Monte Carlo trials per point.

| Item | Statement (source) | Published | Measured | Outcome |
|---|---|---|---|---|
| A1 | The chain driven at paper 2 Table 1: k 100 kHz/s, T 0.5 s, 3 cells, df 5 kHz, 75 kHz steps, 50 kHz/s overall | runs | echo in 3/3 cells; 50 kHz/s overall | reproduced |
| A2 | Range cell c/(2kT) (paper 1 eq. 27) | 3.00 km | 3.00 km (Hann cell 4.50 km) | reproduced |
| A3 | f_D = 2kR'/c to 1 % at V* = 30 m/s (k 50 kHz/s, R' 300 km, 5 MHz; eqs. 26-27) | 1.00 % | 1.002 % | reproduced |
| A4a | Ambiguity interval c/(2 df) (eqs. 48, 70) | 29.98 km | 29.98 km | reproduced |
| A4b | f_R1 = k/df (eq. 51) | 20 Hz | 20.0 Hz | reproduced |
| A4c | Poole's R' as a reference estimator, against the chain's leading edge (eqs. 47-52, 70-71) | method | R' error +0.005 km mean, 0.17 km worst over 150-600 km; leading edge -4.5 km | reproduced |
| A5a | Unambiguous Doppler about 300 m/s at 5 MHz, df 5 kHz (eq. 59) | 300 m/s | 305.9 m/s | **difference explained**: eq. 59 drops the `V* T` terms; with them the bound is 305.9 m/s.  With the beat frequency read at the bin centre, as the hardware did, it is 328 m/s (the half bin, f' in eq. 59). |
| A5b | df down to 1 kHz for 500 m/s (eq. 59) | at least 500 m/s | 1527 m/s | reproduced |
| A6 | RMS phase error under 4 deg above 20 dB (eqs. 76-82, Fig. 3) | 4.0 deg | 4.14 deg [4.00, 4.27] | **difference explained**: Poole's own closed form gives 4.06 deg at 20 dB, inside the interval; "4 deg" is that value rounded. |
| A7 | V* rms 0.56 m/s for a 4 deg phase error at T 0.5 s, 5.9 MHz (paper 2) | 0.56 m/s | 0.556 [0.538, 0.573] | reproduced (with 4 deg on the phase difference; 20 dB per cell gives 5.7 deg on the difference and 0.80 m/s) |
| A8 | Echo criterion: at most the 16 highest values over an amplitude criterion, 1024-point FFT (paper 2) | comparison | same threshold: P_fa 7.5e-5 (known mean) vs 8.9e-5 (P2 detector); Pd equal within 0.01 at 6-14 dB; a 40-bin spread echo keeps 16 vs 39 bins | reported |
| A9 | V' = V* + f dV*/df from a swept ionogram, within 20 % of the direct dR'/dt (paper 2 eqs. 1-6, Figs. 4-5) | -44.5 vs -55.1 m/s (20 %) | -51.4 vs -50.2 m/s (+2.4 %); truth -50.0 | reproduced |
| A10 | Time-delay resolution 45 us from 30 kHz coherence at 14 kHz/s (Barry Table I) | 45 us | 47.9 us | reproduced |
| A11 | Virtual height range 0-1000 km with 100 Hz receiver bandwidth (Barry Table I) | 1000 km | 1071 km | reproduced |
| A12 | The chain driven at Barry's parameters | runs | echo detected | reproduced |
| A13 | 8 W peak, 3 W average | | | not simulated: transmitter hardware |
| A14 | Logarithmic sweep 49 s/octave, 0.5-30 MHz | | 289 s per ionogram | not simulated: a schedule property with nothing stated to check |
| A15 | Gating costs 3 dB and half the echoes | | | not simulated: excluded by decision; tss-sim is CW with a canceller (WP3) |
| A16 | PRF randomly varied 350-800 pps | | | not simulated: the T/R gating pattern (see A15) |
| A17 | No echoes below 1.6 MHz; ionograms comparable with pulse sondes | | | not simulated: antenna and bottomside geometry |

Outcomes: 12 reproduced, 2 difference explained, 1 reported, 5 not simulated,
none open.

**A9's layer.**  A parabolic bottomside layer (base 250 km, semi-thickness
100 km) whose critical frequency rises so that the group height at 6.4 MHz
falls at 50 m/s from 376 km; the layer changes shape, so V* differs from V' and
varies with frequency, as in paper 2 (V* -9.2 against paper 2's -8.7 m/s,
dV*/df -6.6e-6 against -5.6e-6 m/s/Hz).  A layer whose critical frequency
changes at a constant rate does not have a constant V', and then the swept and
fixed experiments, minutes apart, measure different things: that was the first
version, and it disagreed by 35 % while each estimate matched the truth at its
own time.

**Phase error against SNR (A6).**

| SNR dB | MC rms deg [95 % CI] | Poole eq. 82 deg |
|---|---|---|
| 10 | 13.35 [12.90, 13.78] | 12.98 |
| 15 | 7.04 [6.81, 7.26] | 7.23 |
| 20 | 4.14 [4.00, 4.27] | 4.06 |
| 25 | 2.21 [2.14, 2.28] | 2.28 |
| 30 | 1.29 [1.25, 1.33] | 1.28 |
| 40 | 0.41 [0.40, 0.43] | 0.41 |

## `isis1_retrieval`: the measured ionosphere

`data/reference/wp4_p5_isis1` (commit 40b8171): 7,790 ISIS-1 records
(iquality 4 to 9, spacecraft at 500-1100 km, at least 6 points), each with
clutter on and off, so 15,580 runs.  No run failed.  662 CPU hours on 17
workers.  The per-run rows are in `isis1_retrieval_records.csv`.

| | Clutter off | Clutter on |
|---|---|---|
| SCI-05 median (p90) | 0.078 (0.57) | 0.092 (0.51) |
| SCI-05 verdicts: goal / threshold / fail | 4680 / 954 / 2151 | 4278 / 1323 / 2118 |
| SCI-07 km median (p90) | 13.0 (56) | 15.5 (51) |
| SCI-07 verdicts: goal / threshold / fail | 2879 / 2591 / 2315 | 1977 / 3001 / 2716 |
| N_e error over the measured span, median (p90) | 0.064 (0.30) | 0.082 (0.28) |
| Runs within 10 % / 15 % over the measured span | 71 % / 82 % | 62 % / 79 % |
| Coverage gap km median (p90) | 27 (88) | 18 (126) |
| O completeness median | 0.43 | 0.39 |
| Labelling error, median per run | 0.16 | 0.20 |

**This does not reproduce the Chapman campaign.**  On the six Chapman profiles
every P4 run is at goal or threshold on SCI-05 and SCI-07; on measured
ionospheres 28 % fail SCI-05 and 30 to 35 % fail SCI-07.  The median is at goal,
and over the heights ISIS-1 measured the median error is 6.4 %.

What the failures follow (clutter off; clutter on is the same shape):

| Plasma frequency at the spacecraft | Runs | SCI-05 fail | Within 15 % over the span | Labelling error, median |
|---|---|---|---|---|
| under 1.5 MHz | 2295 | 49 % | 71 % | 0.21 |
| 1.5-2.5 MHz | 2570 | 28 % | 80 % | 0.22 |
| 2.5-4 MHz | 1795 | 16 % | 88 % | 0.14 |
| 4 MHz and above | 1130 | 2 % | 99 % | 0.07 |

and, through it, night over day, low foF2 and high latitude (high-latitude
night 45 % fail, low-latitude day 14 %).  Scaling quality matters too: 18 %
fail at index 4, 58 % at index 9.

**Diagnosis on a sample** (16 records rerun with the tracker's report: 8
failing and 4 passing under 1.5 MHz, 4 failing above 3 MHz):

* The f_N prior is not the discriminator.  Under 1.5 MHz the tracker found no
  f_N (no resonance column, no trace start) in 7 of the 8 failing records and
  in all 4 passing ones; above 3 MHz, where it was found from the trace start,
  it was within 0.4 MHz.  The hypothesis carried from P4 is not confirmed.
* The failing retrievals read N_e high, 10 to 60 % over most of the span and up
  to a factor of 2 at its bottom.  An X sample inverted as O reads high by
  `f / (f - f_H)`, 1.2 to 1.3 at a few MHz with f_H near 1 MHz, which matches
  the bulk of the bias.  Failing runs carry about three times the samples on
  neither trace (median 3 against 1).
* The labelling error is 0.16 median per run on ISIS-1 against 0.004 pooled
  on the Chapman campaign with clutter off, and falls as f_N rises (0.21 to
  0.07).  At low f_N the two traces start close together (f_x - f_N is under
  f_H), which is where self-consistency has least to work with.
* 38 % of the SCI-05 failures (clutter off) are within 15 % over the measured
  span: their error is below the lowest scaled height, where the truth is the
  Chapman cap fitted to IRI's hmF2 (`MeasuredTopside`).  Part of the failure
  rate is the truth construction, not the chain.

**Open, for P6 and the dissertation**: the self-consistency labelling under
measured topside gradients at low f_N (f_N near or below f_H, night and high
latitude), and a truth below the scaled profile that does not depend on IRI's
hmF2.  The chain was not changed to fit these numbers.

## The ISIS-1 quality filter

The archive's scaling index runs from 0 (best) to 10 (worst), 4 to 10 in these
files with 6 the most likely (`00readme.txt`).  `wp4_dispersion.isis_scenes`
kept index 6 and above, the worse-scaled records, so the 12 ISIS-1 scenes of
the P3 echo-level study (`data/reference/wp4_p3_dispersion`) were drawn from
quality 6 to 9.  The default is now 4 to 6; `isis_qualities: [6, 7, 8, 9, 10]`
reproduces the old set.

**Effect on P3 (rerun at cd82baa, both sets; `data/reference/wp4_p5_dispersion_q4to6`
and `wp4_p5_dispersion_q6to10`).**  At the same commit the old set reproduces the
P3 echo-level numbers exactly, so the difference below is the filter alone.  The
chain part of the study uses only the Chapman scenes and is identical in both.

| Treatment | Loss dB median (p90), q 6-10 | q 4-6 | Near f_N / cusp dB, q 6-10 | q 4-6 |
|---|---|---|---|---|
| D0 | 4.1 (9.2) | 3.8 (9.0) | 8.2 / 8.2 | 8.0 / 7.7 |
| D2 (span 32) | 0.1 (5.4) | 0.1 (5.0) | 1.1 / 3.2 | 1.5 / 1.6 |
| **D2w (span 64)** | **0.1 (2.0)** | **0.1 (2.0)** | **0.3 / 0.6** | **0.4 / 0.7** |
| D3 | 0.3 (4.4) | 0.3 (4.4) | 5.0 / 5.6 | 4.5 / 4.5 |
| T1 | 3.2 (8.2) | 3.0 (7.9) | 6.5 / 7.6 | 6.2 / 7.3 |

G2 (the rate search over 64 cells) and G3 (accept the skirt) stand: D2w is
unchanged at 0.1 dB median and 2.0 dB at the 90th percentile, and still meets
the G6 draft (1 dB / 3 dB).  `data/reference/wp4_p3_dispersion` is kept as the
record the decisions were made on.
