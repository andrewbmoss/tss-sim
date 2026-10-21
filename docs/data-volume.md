# Data products and daily volume (`tss study data_volume`)

WP0 closure task 2, Phases 7.1-7.3 and 7.7 (2026-10-20).  INS-13 states the data
product and its volume; Rev G gave the raw-ionogram and trace-extraction
arithmetic only.  `src/tss/studies/data_volume.py` defines every product the
four observing modes of INS-27 make, sizes it from the requirements that set its
axes, and sums it over a day for reference duty cycles per orbit class.

**Scope (Andrew, 2026-10-20).**  The contact analysis, the link budgets and the
on-board processing budget (plan items 7.4-7.6) are outside this dissertation.
The work stops at how much data each mission-specific design must handle; the
downlink, the ground segment and the processor are sized by the mission.

    tss study data_volume --out OUT

writes `data_volume.csv`, `.json` and `.md` (12 cases: A1 and A2, 900 and
1400 km, three duty cycles).  Every byte-level choice is a field of
`ProductParams` and can be overridden under `data_volume:` in the config.

## 7.1 The products

| product | mode | contents | size, A1, 900 km | what sets it |
|---|---|---|---|---|
| L0 raw | calibration bursts | 16-bit digitiser samples at f_s = 4 K_f R_max / c, per receive channel | 157 kB/s; 668 kB per sweep | INS-06, INS-04 |
| L1a ionogram | ionogram (A1, A2) | 8-bit amplitude per step and 10 km bin over R_max, per channel; 64 B header; 1 B AGC per step | 16,695 B matrix; 16,918 B delivered | INS-03, INS-04, INS-13 |
| L1b traces | ionogram | up to 4 echoes per step (O, X, Z, one off-vertical or ground), 4 B each: virtual range, amplitude, mode and quality flags (INS-16); vector adds azimuth, zenith and polarization, 8 B | 2,608 B (2,784 B for A2's 85 steps x 8 B) | INS-02/03, INS-16, INS-18, SCI-27 |
| L2 profile | ionogram | Ne(h) and sigma at SCI-06's 20 km sampling from h_sat to hmF2, float32, plus NmF2, hmF2, foF2, fxF2, flags; goal: 10 km and H(h) with sigma | 384 B threshold; 1,232 B goal | SCI-06, SCI-18 |
| Doppler | Doppler | per sub-band dwell: range x Doppler map, 10 km bins x N bins at 8 bit (N = 20 threshold, 60 goal, 18 ms sub-sweeps), plus 8 echoes x 8 B (range, LOS velocity, sigma, amplitude, flag) | 2,228 B per 0.36 s dwell; 6,428 B per 1.08 s (goal) | INS-09, INS-15, SCI-15, INS-27 |
| Drift | drift (route one) | phase, amplitude and range of up to 4 echoes at 45 Hz, 5 B each, 2 s dwells | 1,864 B per dwell | TV Ch. 8 route one, INS-27 |
| Housekeeping | all | fixed rate | 50 B/s (4.3 MB/day) | - |

SCI-23's drift (route two) is the displacement of one depletion between two
nodes' L1b and L2 products, so it adds no product.  At 1400 km R_max is 1800 km:
180 bins, a 28,620 B matrix, 584 B / 2,032 B profiles and 3,728 B Doppler dwells.

**Product sets.**

* **Threshold:** L1b and L2 (threshold fidelity) for every ionogram; full L1a
  for one ionogram in ten, as a calibration and validation subset; the Doppler
  (20-bin) and drift products when those modes run; housekeeping; one L0 sweep
  per orbit.
* **Goal:** the same with L1a for every ionogram, L2 at 10 km with H(h), and
  60-bin Doppler maps.
* **Raw reference:** L0 for every sweep.  What the instrument makes before any
  on-board processing; shown for scale, not a delivery option.

## 7.3 Reference duty cycles

Reference cases for sizing, not requirements; each mission declares its own
(MIS-08, and the proposed MIS-11 below).  Doppler and drift dwells displace
ionograms; they are not added to a full ionogram day.

| duty cycle | ionogram | Doppler | drift | meaning |
|---|---|---|---|---|
| O1 | 100 % | 0 | 0 | global class, ionograms back to back at the sweep time |
| O1_6p7s | 100 % at 6.7 s | 0 | 0 | the same at INS-07's threshold cadence |
| O2 | 88.75 % | 7.5 % | 3.75 % | bubble class: ionograms outside 18-03 LT; inside it (9 h of 24, SCI-26, MIS-07) Doppler 20 %, drift 10 %, ionograms the rest |

## 7.2 Daily volume, MB/day

| config | h_sat | duty | ionograms/day | L1a (all) | L1b | L2 | Doppler | L0 cal | **threshold set** | **goal set** | raw reference |
|---|---|---|---|---|---|---|---|---|---|---|---|
| A1 | 900 | O1 | 20,276 | 343 | 52.9 | 7.8 | 0 | 9.3 | **109** | **435** | 13.5 GB |
| A1 | 900 | O1_6p7s | 12,896 | 218 | 33.6 | 5.0 | 0 | 9.3 | **74** | **281** | 13.5 GB |
| A1 | 900 | O2 | 17,995 | 304 | 46.9 | 6.9 | 40.1 | 9.3 | **141** | **429** | 13.5 GB |
| A2 | 900 | O1 | 37,928 | 1,021 | 105.6 | 14.6 | 0 | 94.6 | **321** | **1,272** | 256.5 GB |
| A2 | 900 | O1_6p7s | 12,896 | 347 | 35.9 | 5.0 | 0 | 94.6 | **174** | **498** | 256.5 GB |
| A2 | 900 | O2 | 33,661 | 906 | 93.7 | 12.9 | 40.1 | 94.6 | **339** | **1,182** | 256.5 GB |
| A1 | 1400 | O1 | 20,276 | 585 | 52.9 | 11.8 | 0 | 14.5 | **142** | **698** | 23.2 GB |
| A1 | 1400 | O1_6p7s | 12,896 | 372 | 33.6 | 7.5 | 0 | 14.5 | **97** | **451** | 23.2 GB |
| A1 | 1400 | O2 | 17,995 | 519 | 46.9 | 10.5 | 67.1 | 14.5 | **198** | **690** | 23.2 GB |
| A2 | 1400 | O1 | 37,928 | 1,747 | 105.6 | 22.1 | 0 | 146.7 | **453** | **2,080** | 439.8 GB |
| A2 | 1400 | O1_6p7s | 12,896 | 594 | 35.9 | 7.5 | 0 | 146.7 | **254** | **807** | 439.8 GB |
| A2 | 1400 | O2 | 33,661 | 1,550 | 93.7 | 19.7 | 67.1 | 146.7 | **490** | **1,932** | 439.8 GB |

(L2 and Doppler at threshold fidelity; drift is 3.0 MB/day in every O2 case;
housekeeping 4.3 MB/day everywhere.)

* **Rev G reproduced** (tests): 16,695 B per raw ionogram, 215 MB/day at 6.7 s
  and 339 MB/day at 4.26 s; trace extraction 2,544 B per ionogram, 33-52 MB/day
  (Rev G 32-51, which rounded the record to 2.5 kB); raw samples 13.5 GB/day at
  78.4 kHz.
* **What drives the volume** is the L1a ionogram and how often it is kept.
  Delivering traces and profiles for every ionogram and full ionograms for one in
  ten cuts the A1 day from 435 to 109 MB (a factor of 4); everything else is
  small beside it.
* **Back to back, the step count does not change the volume:** each step adds a
  column and its 26.8 ms dwell, so a start below 0.33 MHz costs volume only at a
  fixed cadence.  What A2 costs is cadence: 2.28 s sweeps back to back nearly
  double the ionograms per day (SCI-09's goal), and three channels triple each
  L1a.  At INS-07's 6.7 s A2 is 174 / 498 MB/day.
* **The goal altitude** (1400 km, Phase 4) raises R_max to 1800 km: L1a and the
  Doppler maps by 1.71, the threshold set at A1 / O1 from 109 to 142 MB/day,
  the goal set from 435 to 698 MB/day.
* A2's L0 calibration sweep is 6.8 MB (three channels at INS-06's 494 kHz; the
  schedule's peak rate is 35.3 MHz/s, below the 37.3 MHz/s INS-06 quotes for the
  step at exactly 20 MHz), so one per orbit is 95 MB/day; a mission may thin it.

**INS-19, three receive channels against one** (A1 sweep, 900 km, O1): raw
x3.00, L1a x2.97, L1b x1.98, threshold set x2.27, goal set x2.72.  "Roughly a
factor of two, not three" holds for the trace product; the delivered sets come
out higher because the L1a subset and the calibration sweeps keep all three
channels.

## 7.7 Draft requirement text (for the Rev G amendment, G3)

**INS-13, Data products and volume (redrafted).**  *Threshold:* the instrument
shall deliver, for every ionogram, the L1b echo traces (up to four echoes per
frequency step, each with virtual range, amplitude and mode and quality flags,
the INS-16 flag among them) and the L2 profile Ne(h) with its per-point
uncertainty at SCI-06's height sampling (SCI-18); full L1a ionograms (8-bit
amplitude per frequency step and 10 km range bin over INS-04's range, with time,
position, attitude, sweep settings and AGC) for at least one ionogram in ten;
the Doppler-mode and drift-mode products whenever those modes run; housekeeping;
and one L0 sweep per orbit for calibration.  *Goal:* L1a for every ionogram, L2
at 10 km with H(h) and its uncertainty, and 60-bin Doppler maps.  *Volume:* the
daily volume follows from the configuration, altitude and duty cycle the mission
declares (MIS-11); for reference, A1 at 900 km with ionograms back to back is
109 MB/day at threshold and 435 MB/day at goal (74 and 281 MB/day at INS-07's
6.7 s), the O2 bubble duty cycle 141 and 429 MB/day, the vector configuration
A2 321 and 1,272 MB/day back to back (174 and 498 MB/day at 6.7 s), and A1 at
SCI-02's 1400 km goal altitude 142 and 698 MB/day (tss-sim `data_volume`).
*Verification:* A (data-volume arithmetic, `tss study data_volume`), I (product
definitions).

**INS-28, On-board product generation (new).**  The instrument shall form the
L1b traces and the L2 profile with its uncertainty on board from every
ionogram, keeping pace with the ionogram cadence, so that the threshold product
set of INS-13 can be delivered without the L1a or L0 data of every sweep.  It
shall retain L0 data for commandable calibration bursts and L1a ionograms for a
commandable subset, and shall flag the profiles INS-16 identifies.  The
processor's throughput, memory and storage are sized by the mission against its
declared duty cycle (MIS-11); this requirement fixes what is produced, not the
hardware that produces it.  *Verification:* T (bench: products formed from
recorded sweeps at cadence), A.

**INS-19, Data volume in the vector configuration (updated).**  Three receive
channels triple the raw data rate and the L1a ionogram.  With on-board trace
extraction and arrival-angle reduction (INS-28) the L1b product grows by a factor
of 2.0, since the range and frequency axes are shared and each echo gains only
its arrival angle and polarization sense; the INS-13 threshold product set grows
by 2.3 and the goal set by 2.7 (A1 sweep, 900 km, O1), because the L1a subset
and the calibration sweeps keep all three channels.

**MIS-11, Data-handling declaration (proposed; this is gate G5).**  Each mission
shall declare its observing-mode duty cycle, the INS-13 product set it delivers
(threshold or goal) and the resulting daily volume from the INS-13 arithmetic,
and shall show that its downlink and on-board storage carry that volume.  For a
constellation, the declaration is per node and the ground segment is shared at
mission level.  Ground stations, contact plan, link budget and processor design
are mission designs outside this requirement set.  *Verification:* I
(declaration review).

**Adopted 2026-10-20.**  Andrew accepted MIS-11 as drafted (gate G5), and the
four texts above entered `TSS_Science_Requirements_RevG.docx` (sections 5.13,
6.1, 6.2, 6.5) and `TSS_SRM_RevG.xlsx` (Instrument & Mission; Derivations
section 24, which reproduces the A1 / O1 / 900 km threshold and goal sets,
108.6 and 434.6 MB/day, as live formulas) with the amendment that closed WP0.
The catalogue now holds INS-28 and MIS-11; `requirements.yaml` lists them as
bench and mission requirements.
