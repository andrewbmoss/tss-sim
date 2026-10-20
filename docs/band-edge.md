# The band-edge study (`tss study min_frequency`)

The lower edge of the sounding band (INS-01a) and the profile completeness behind
SCI-02, from IRI-2020 sampled along the orbit.  Ported on 2026-10-19 from
`Dissertation/MinFreqStudy`, the harness that set INS-01a in Rev G, so that the
WP0 closure work runs, is tested and leaves a manifest inside this repository.

## Running it

    tss study min_frequency                                   # plan band_edge_reference
    tss study min_frequency --set min_frequency.plan=<plan>   # any plan in config/studies
    tss study min_frequency --set min_frequency.samples=2000  # override one key

A plan (`config/studies/*.yaml`) states altitudes, orbit classes (O1 80 deg, O2
30 deg), solar cases, topsides, samples per case, seed, epoch window and the
altitude grid.  A solar case is either

* `{name: min, mode: r12, r12: 10}`: a user R12, the harness's route.  The mex
  derives F10.7 and IG12 from R12 and the foF2 storm model is off; or
* `{name: min_2008, mode: epoch, epoch_start: "2008-07-01T00:00:00",
  duration_days: 365.25, iri_options: {foF2_storm: "on"}}`: R12 = -1, so IRI
  takes R12, IG12, daily and 81-day F10.7 and ap for every sample from its own
  index files.  The span is checked before MATLAB starts
  (`tss.l0_ionosphere.backends.iri_index_span`; 1958-01-01 to 2024-12-09 for
  PHaRLAP 4.7.4).

The IRI route is the `iri` L0 tier with `model: pharlap`.  `backend: analytic`
runs the stand-in, for tests only.

## Outputs

| File | Content |
|---|---|
| `minfreq_summary.csv` | One row per case.  The harness's columns in its order, then `solar_mode`, `epoch_start`, `duration_days`, `fn_sat_p01_MHz`, `hmf2_p001_km`, `hmf2_p999_km`, `frac_hsat_minus_hmf2_lt_100km`, `samples_file`. |
| `completeness_curves.csv` | Coverage and altitude/density completeness against f_min, per case. |
| `ensemble.csv` | Per orbit class, altitude and solar case: the minimum across topsides at 95 and 99 per cent, the unsafe maximum, and the governing model. |
| `samples/*.npz` | Per-sample scalars: time, position, local time, fN at the spacecraft, Ne there, NmF2, hmF2, and the IRI extras (`EXTRA_COLUMNS`, including the indices of an epoch run).  Profiles are not kept. |
| `plan.json`, `manifest.json`, `summary.md` | The resolved plan and the run record. |

## Holding it to the Rev G run

`data/reference/band_edge` keeps the 2026-09-07 MinFreqStudy output that set
INS-01a.  `tests/unit/test_band_edge.py` carries the harness's 18 checks, the
ensemble-direction test, a fixture test that `tss.platform.orbits` draws the
same samples as the harness's `orbits.py`, and a `-m pharlap` test that reruns
two cases and matches the reference to 1e-4.  The full 60-case comparison is a
run of the reference plan (about an hour) and is recorded below.

## Reference regression (2026-10-19)

`tss study min_frequency` (plan `band_edge_reference`, 60 cases, 1.2 million
samples, 77 min on the development machine) against
`data/reference/band_edge/minfreq_summary.csv`:

* The same 60 cases in the same order, with identical sample counts and no
  dropped samples.
* Every percentile column (f_min at 90, 95 and 99 per cent; fN at the
  spacecraft at the 5th, 50th and 95th percentile and its median) agrees to
  within 2.4e-7 relative.  foF2 and hmF2 medians and the band ratio agree to
  within 1.8e-6, the precision of the extras file (six significant figures).
* The ensemble governing value is unchanged: O1, 900 km, solar minimum,
  NeQuick, **0.3342 MHz** at 95 per cent (0.2700 MHz at 99 per cent); the
  unsafe maximum combination would give 0.4766 MHz.
* The completeness curves agree to the four decimals they are written with in
  9,336 of 9,600 values.  The rest differ by at most 0.0075, all in cases with
  a few samples whose hmF2 is at or above the spacecraft (8 of 20,000 in O2,
  500 km, solar maximum, IRI-2001, where hmF2 reaches 536 km).  For those the
  completeness span is clamped to 1 m and the value turns on the sixth
  significant figure of hmF2, which the extras file rounds and the harness's
  file did not.  They are a property of the samples where SCI-02's threshold is
  at issue, not a change in arithmetic, and Phase 3 treats them as such.

## Phase 2: solar minimum on the epoch's own indices (2026-10-20)

Plan `band_edge_solarmin_epoch`, 150 cases (3 h 6 min): the twelve months
centred on each solar minimum inside PHaRLAP's index files, with IRI taking
R12, IG12, F10.7 and ap for every sample from `ig_rz.dat` and `apf107.dat`
(foF2 storm model on), and 2008-09 again with every storm model off.  Compared
with the reference run by `scripts/band_edge_compare.py`; the summary,
ensemble and comparison are kept in `data/reference/band_edge/phase2_solarmin_epoch`.

Indices IRI ran on (median over the samples) and the governing cell, O1 at
900 km, NeQuick in every window:

| Window | R12 | IG12 | F10.7 | 95 % edge | vs Rev G | 99 % edge | vs Rev G |
|---|---|---|---|---|---|---|---|
| Rev G (user R12 10) | 10 | 2.1 | 71.1 | 0.334 | | 0.270 | |
| 1986-03 to 1987-02 | 13.8 | 8.5 | 71.1 | 0.355 | +6.4 % | 0.281 | +4.0 % |
| 1995-11 to 1996-10 | 8.7 | 0.9 | 70.5 | 0.330 | -1.2 % | 0.254 | -6.0 % |
| 2008-05 to 2009-04 | 2.1 | -9.1 | 68.0 | 0.282 | -15.7 % | 0.215 | -20.3 % |
| 2019-06 to 2020-05 | 2.1 | -10.3 | 69.2 | 0.279 | -16.4 % | 0.212 | -21.3 % |

* The Rev G solar level is a 1996-like minimum.  The two most recent minima
  were deeper (IG12 about -10 against 2.1) and lower the band edge by 11 to
  16 per cent at 95 per cent coverage and 16 to 22 per cent at 99 per cent, in
  every cell, both orbit classes, all altitudes.
* The storm model does not matter at solar minimum: with it on, the 2008-09
  edge is 0.0 to 0.3 per cent higher than with every storm model off.  The
  whole shift is the solar index.
* NeQuick is the minimum model in every cell, as in Rev G.
* Decision rule of the WP0 plan: the lowest epoch-driven value, 0.279 MHz, is
  within INS-01a's -20 per cent tolerance (floor 0.264 MHz).  The 99 per cent
  goal point, 0.212 MHz, is within -20 per cent of the 0.26 MHz goal (floor
  0.208).
* ISIS-1 at solar minimum (61 CRC records, 1976-77, paired with IRI in epoch
  mode on 2026-10-18) cannot test the deep minima, whose IG12 is far below
  1976-77's.  At its own epochs the modelled 5th percentile is above the
  measured one: 0.775 against 0.652 MHz at 500-600 km (23 records) and 0.639
  against 0.284 MHz at 600-700 km (31 records); 4 of the 61 records lie below
  0.33 MHz and 3 below 0.28 MHz.
* So the tolerance is now carrying two things.  Rev G set it for the models'
  overstatement of the low tail, about 19 per cent against ISIS-1; the deep
  minimum alone now takes 16 per cent of it.
* Gate G1 decided (Andrew, 2026-10-20): keep the Rev G 0.33 MHz.  INS-01a's
  threshold (0.33 MHz) and goal (0.26 MHz) stand; its notes are to record the
  0.28 MHz deep-minimum edge, the roughly 90 per cent coverage a 0.33 MHz start
  gives in 2008-09 and 2019-20 conditions, and that the tolerance is nearly used.

## Phase 3: SCI-02 threshold, the spacecraft's height above the F2 peak (2026-10-20)

SCI-02's threshold is coverage from the spacecraft down to and including the F2
peak, for h_sat 500-900 km and hmF2 anywhere in 200-500 km; MIS-01 adds that
the spacecraft must be above the peak at all times sampled.  Phase 3 reads that
off the per-sample hmF2 the band-edge runs keep, with
`tss.studies.sci02_margin` and `scripts/sci02_threshold.py`:

    python scripts/sci02_threshold.py OUT revG=<band_edge_reference> \
        solarmin_epoch=<band_edge_solarmin_epoch> storm_epoch=<sci02_storm_epoch> \
        storm_m3000=<sci02_storm_m3000>

Outputs are kept in `data/reference/sci02_threshold` (94 cells, 1.88 million
samples).  hmF2 is a bottomside quantity and is identical, to the metre, across
the three topside models in every cell of every run, so one topside per cell is
read.

Two runs were added for the storm check (3.3), because a user R12 turns IRI's
storm models off and the Phase 2 runs are all at solar minimum:

* `sci02_storm_epoch`: epoch indices, foF2 storm model on, 500 and 600 km, O1
  and O2, over the cycle 23 maximum (2000-07 to 2001-06, median R12 111) and
  densely over the Halloween storms (2003-10-15, 41 days, 1,128 samples per
  cell at Kp 7-9), each with hmF2_storm off and on.  16 cases, 21 minutes.
* `sci02_storm_m3000`: the same 500 km cells with IRI's older hmF2 route
  (M3000F2), the only one through which the foF2 storm model can reach hmF2.
  8 cases, 10 minutes.

Results:

| h_sat | solar | samples below the peak | within 50 km | within 100 km | smallest margin |
|---|---|---|---|---|---|
| 500 km | maximum (R12 120; 2000-01; 2003) | 0-0.07 % | 0-1.2 % | 0.4-8.4 % | -36 km |
| 500 km | minimum (R12 10; 1986-2019 epochs) | 0 | 0 | 0 | 108 km |
| 600 km | maximum | 0 | 0 | 0-0.08 % | 68 km |
| 700-900 km | any | 0 | 0 | 0 | 182 km |

* **Envelope.** Across all 94 cells the 0.1 and 99.9 percentiles of hmF2 lie
  inside 200-500 km (200.7 km in O1 / 800 km / 2008-09; 496 km in O2 / 500 km /
  2000-01).  The extremes are 194 km (solar minimum, 2008-09) and 536 km (solar
  maximum, O2); at most 0.065 per cent of a cell's samples fall below 200 km
  and 0.075 per cent above 500 km.  At 194 km the 900 km true range is 706 km,
  so SCI-08's 1.5 times is 1059 km against INS-04's 1050 km: 9 km, in under 0.1
  per cent of solar-minimum samples.
* **Above the peak.** Met everywhere from 600 to 900 km, at every solar level.
  At 500 km it is met at solar minimum (margin at least 108 km) but not "at all
  times" at solar maximum: 0.03-0.07 per cent of samples have the peak above
  the spacecraft, by up to 36 km, and 4-8 per cent are within 100 km.  These
  are afternoon and evening samples: in O2 / 500 km / R12 120, 18-21 LT has
  0.28 per cent below the peak and 21 per cent within 100 km.
* **Storms.** IRI cannot verify SCI-02 under storm conditions.  With the
  default hmF2 model (Shubin-COSMIC) hmF2 is identical with hmF2_storm off and
  on, and flat across Kp bins (median 303-323 km from Kp 0-3 to Kp 7-9 at solar
  maximum), because that model does not take foF2.  Through M3000F2 the storm
  model lowers hmF2 slightly (0.1 percentile 232 to 215 km, minimum 211 km, at
  Kp 7-9 in the Halloween window) and raises it nowhere.  Neither route has the
  storm-time uplift by penetration electric fields, which is the case that
  would put the peak above a 500 km spacecraft; that needs measured hmF2
  (COSMIC, Digisonde) rather than IRI.
* **Verdict.** SCI-02's threshold is met over 600-900 km, and its 200-500 km
  hmF2 envelope holds at the 0.1/99.9 per cent level, in IRI climatology at
  every solar level.  At 500 km and solar maximum, MIS-01's "above the F2 peak
  at all times sampled" fails in 0.03-0.07 per cent of samples, and the
  storm-time case is not verified.
* **G3 decision (Andrew, 2026-10-20): state MIS-01 as a fraction of samples.**
  The fraction is set at 99.9 per cent, the same level as the 0.1/99.9
  percentile reading of the hmF2 envelope (`MIS01_ABOVE_PEAK_FRACTION` in
  `tss.studies.sci02_margin`; the verdict reports `above_peak_fraction_met`).
  Every cell meets it; the worst, O2 / 500 km / 2000-01, is at 99.93 per cent,
  so the margin is 0.03 per cent of samples, and it is the margin a storm-time
  uplift IRI does not model would have to fit inside.  The wording enters the
  requirements document with the WP0 close-out.

## Phase 4: SCI-02's goal, h_sat 1000-1400 km (2026-10-20)

Plan `band_edge_goal_1400`, 120 cases (2 h 49 min): 1000 to 1400 km in 100 km
steps, O1 and O2, the three topsides, solar minimum as Rev G ran it (R12 10)
and on the 2008-09 and 2019-20 epochs (foF2 storm on), and solar maximum (R12
120), on a 150-2000 km grid.  Two things are new in the harness for it: the
MATLAB generator writes IRI's O+/H+ transition height for every sample (the
`h_oplus_hplus_km` extra: the lowest height above hmF2 where H+ reaches O+,
interpolated in log(H+/O+)), and every summary row carries the coverage and
altitude completeness of a sweep that starts at INS-01a's 0.33 MHz.
`scripts/goal_knockons.py` turns the run into the tables below; outputs are
kept in `data/reference/band_edge_goal_1400`.

**Band edge and what a 0.33 MHz start retrieves** (ensemble at the minimum of
the per-topside percentiles; coverage and completeness for the worst topside;
NeQuick governs every cell):

| h_sat | governing 95 % edge, MHz | 99 % | coverage at 0.33 MHz | mean altitude completeness | 5th percentile |
|---|---|---|---|---|---|
| 1000 km | 0.246 (O1, 2008-09) | 0.186 | 85.6 % | 0.970 | 0.740 |
| 1100 km | 0.217 (O1, 2019-20) | 0.168 | 80.0 % | 0.951 | 0.648 |
| 1200 km | 0.194 (O1, 2019-20) | 0.153 | 74.7 % | 0.932 | 0.577 |
| 1300 km | 0.180 (O1, 2019-20) | 0.138 | 70.5 % | 0.915 | 0.530 |
| 1400 km | 0.164 (O1, 2019-20) | 0.128 | 65.3 % | 0.892 | 0.477 |

* The edge keeps falling with altitude, by 7-12 per cent per 100 km.  The
  deep minima sit 15-17 per cent below Rev G's R12 = 10 at every altitude (at
  1400 km, O1: 0.196 MHz at R12 10, 0.167 and 0.164 MHz in 2008-09 and
  2019-20), the same shift Phase 2 found at 500-900 km.
* At solar maximum a 0.33 MHz start retrieves essentially the whole topside
  at every goal altitude (coverage at least 99.4 per cent; 95 per cent edges
  0.50-0.92 MHz).  The goal's shortfall is a solar-minimum one, and O2 fares
  better than O1 (1400 km, 2019-20: coverage 74.9 per cent, edge 0.206 MHz).
* So with the sweep starting at the threshold edge (G2 option A), a 1400 km
  spacecraft at a deep minimum gets a complete topside in about two profiles
  in three, keeps 89 per cent of the h_sat - hmF2 span on average, and in one
  profile in twenty loses more than half of it.

**O+/H+ transition** (identical across the three topsides to 2e-5 km, so read
from NeQuick; percentage = samples with the transition below the spacecraft;
a transition above the 2000 km grid counts as not below):

| class, solar | transition p05 / p50 / p95, km | below at 1000 km | 1200 km | 1400 km | above 2000 km |
|---|---|---|---|---|---|
| O1, minimum (all three) | about 540 / 820 / 1450 | 68-69 % | 81 % | 91-92 % | 1 % |
| O2, minimum (all three) | about 527 / 682 / 990 | 95-96 % | over 99.5 % | 100 % | 0 |
| O1, maximum | 805-811 / 1280-1290 / 1880-1890 | 19 % | 31 % | 45 % | 27 % |
| O2, maximum | about 790 / 1135 / 1600-1630 | 39 % | 56 % | 80 % | 0.5 % |

* At solar minimum IRI puts the transition well below the 1000-1400 km
  "transition level" SCI-02's basis cites: the median is 820 km (O1) and 680
  km (O2), so a 900 km spacecraft already sits above it in most samples.  At
  solar maximum it rises to a 1135-1290 km median, and 1400 km is above it in
  45 per cent (O1) to 80 per cent (O2) of samples.  The altitude buys the
  transition at solar maximum, which is also where a 0.33 MHz start costs
  nothing.
* IRI's ion composition (RBV-2010 and TTS-2005) above 1000 km is the least
  constrained part of the model; Phase 5's ISIS-1 comparison is about the
  electron density, not the composition.

**Knock-ons at 1400 km** (Rev G derivations with h_sat = 1400 km in place of
900 km; `goal_knockons.json`):

| item | Rev G (900 km) | 1400 km |
|---|---|---|
| INS-04 / SCI-08, 1.5 x (h_sat - 200 km) | 1050 km | 1800 km |
| INS-06, 4 K_f R_max / c, A1 (5.6 MHz/s) | 78.4 kHz | 134 kHz |
| INS-06, A2 (37.3 MHz/s) | 523 kHz | 896 kHz |
| INS-13 ionogram, 10 km bins, 8 bit | 16.7 kB (105 bins) | 28.6 kB (180 bins) |
| INS-13 per day at the A1 sweep (4.26 s) | 339 MB | 580 MB |
| INS-13 per day at INS-07's 6.7 s | 215 MB | 369 MB (401 MB from a 0.164 MHz start) |
| INS-03, A1 steps to 20 MHz | 159 from 0.33 MHz | 159; 173 from 0.164 MHz (sweep 4.64 s) |
| Spreading, two-way, against 900 km | - | +4.7 dB (peak at 200 km), +5.1 dB (median hmF2) |
| INS-20 half-wave length | 454 m at 0.33 MHz | 913 m at 0.164 MHz |
| INS-10 SNR, provisional (analytic link budget, peak at 300 km) | median 33.6, min 27.1 dB | median 28.5, min 21.8 dB (-5.1 / -5.3 dB) |

At the continuous A1 sweep the daily volume does not depend on the step
count (each step adds its dwell as well as its column), so a lower start
costs volume only at a fixed cadence.  INS-10 stays above its 10 dB
threshold on this model, but the number is a difference against 900 km, not
a verdict, until the 7 dB drop on PHaRLAP is traced.

**For G2.**  Option A (keep h_sat to 1400 km, sweep from 0.33 MHz, score by
completeness) costs nothing at the antenna and gives, at a deep minimum, 65
per cent coverage and 0.89 mean completeness at 1400 km; option B (the goal
gets its own edge) means 0.164 MHz, a 913 m half-wave and a start nearly four
times below the 0.6 MHz bottom of the sky-noise reference; option C ("up to X km
where coverage is at least Y") gives X of about 1000 km for Y = 85 per cent
and just under 1200 km for Y = 75 per cent in 2019-20 conditions (O1).  Phase 5
(ISIS-1 at 900-1400 km) is still to come and tests the topside these numbers
rest on.
