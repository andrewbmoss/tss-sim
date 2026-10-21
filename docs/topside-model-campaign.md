# The topside difference between PyIRI and Fortran IRI-2020

Campaign of 2026-10-17/18. Reproduce with

    python scripts/topside_campaign.py grid    OUT_DIR
    python scripts/topside_campaign.py isis    OUT_DIR --crc ../Dissertation/ISISValidation/isis1_ne.asc --t1 ../Dissertation/ISISValidation/t1_results.csv
    python scripts/topside_campaign.py analyse OUT_DIR
    python scripts/topside_campaign.py isis-epoch    OUT_DIR --crc ... --t1 ...   # real indices
    python scripts/topside_campaign.py analyse-epoch OUT_DIR

`grid` needs PyIRI and MATLAB with PHaRLAP (`scripts/topside_campaign_iri.m`); `analyse` writes
`summary.json` and four figures. Written up in Technical Volume section "Two Implementations of One
Topside Model" (`tv_literature.tex`, label `sec:twoiri`).

## Grid

Lat -80..80 step 10, lon -180..150 step 30, UT every 3 h, 15 Mar/Jun/Sep/Dec 2024, F10.7 70/120/180
(R12 from F10.7 = 63.75 + 0.728 R12 + 0.00089 R12^2, the PHaRLAP mex relation). 19,584 profiles per
code, 100-2000 km. Statistics cos(lat) weighted. IRI-2020 options: NeQuick topside, URSI foF2,
SHU-2015 hmF2, ABT-2009, F10.7 limited for topside.

| | F10.7 70 | 120 | 180 |
|---|---|---|---|
| NmF2 ratio PyIRI/IRI (median) | 1.016 | 1.026 | 1.011 |
| Ne ratio 500 km | 0.45 | 0.52 | 0.57 |
| Ne ratio 700 km | 0.41 | 0.44 | 0.46 |
| Ne ratio 900 km | 0.42 | 0.44 | 0.45 |
| fp p05 at 900 km, IRI / PyIRI (MHz) | 0.35 / 0.22 | 0.68 / 0.44 | 0.99 / 0.65 |

The ratio barely varies with geography, season or local time (within 0.03 at 900 km).

## Cause

Both codes rebuild exactly (to 1e-5) from their own NmF2, hmF2, foF2, M(3000)F2 with the NeQuick
semi-Epstein topside (g = 0.125, r = 100), B2bot from (dN/dh)max, and the Coisson et al. (2006) k.
The whole difference is the H0 law:

* IRI-2020 (irisub.for, itopn = 2): H0 = k' B2bot, k' = (k e^{2(k-1)} + 1)/(e^{2(k-1)} + 1).
  Same as NeQuick 2 / ITU-R P.531 as reimplemented in Orekit 13 (`NeQuickItu.computeH0`).
* PyIRI 0.1.7 (`sh_library.thickness_F2`, and `main_library.thickness`): x = (k B2bot - 150)/100,
  v = 0.041163 x^2 - 0.183981 x + 1.424472, H0 = k B2bot / v. This is the step of the original
  NeQuick topside, which there follows a seasonal k clipped to 2..8 (Coisson et al. 2008, JASTP 70,
  eqs. 4-8); the revision to the Coisson k dropped it (eqs. 9-10), and NeQuick G kept the original. With the Coisson k (median 2.35), k B2bot is
  54-76 km and v is 1.57-1.65 (minimum possible 1.22), so H0 is 0.63-0.66 of IRI-2020's.

Factor split at 900 km, F10.7 120: NmF2 1.028, hmF2 1.000, H0 inputs 0.998, H0 law 0.431.
PyIRI with IRI-2020's H0 law matches IRI-2020 to 0.6-2.5% at 500-900 km.

## ISIS-1 (5,118 CRC profiles, as T1), solar maximum, median model/measured dB

| band | N | IRI-2020 | PyIRI | IRI law, measured peak | PyIRI law, measured peak | p05 obs / IRI / PyIRI (MHz) |
|---|---|---|---|---|---|---|
| 500-600 | 971 | +1.09 | -2.32 | +0.93 | -2.44 | 1.51 / 2.22 / 1.45 |
| 600-700 | 1655 | +0.77 | -2.63 | +0.53 | -2.99 | 1.16 / 1.73 / 1.13 |
| 700-800 | 1181 | +0.55 | -3.01 | +0.34 | -3.38 | 1.05 / 1.30 / 0.86 |
| 800-900 | 1004 | +0.59 | -2.98 | +0.15 | -3.61 | 0.90 / 1.07 / 0.70 |

H0 fitted to each profile (5,078 fits): fit / IRI law median 0.93 (p10 0.70, p90 1.38); fit / PyIRI
law 1.44; 74% of profiles closer to the IRI law (90% above 60 deg latitude). Solar minimum: 61
profiles, fit / IRI law 0.79, evenly split.

## ISIS-1 with each epoch's own indices (2026-10-18)

`isis-epoch` reruns the comparison with real indices instead of the T1 solar level (the CRC record's
Rz12 rounded to 5, storm model off). Fortran IRI-2020 runs with R12 = -1, so it reads R12, IG12,
daily and 81-day F10.7 and 3-hourly ap for the epoch from its own `ig_rz.dat` and `apf107.dat`
(`scripts/topside_isis_epoch_iri.m`, density evaluated at h_sat directly). Configurations: binned
(control), epoch with every storm model off, epoch with IRI's default foF2 storm model on, epoch with
all storm models on, and the two IRI-2001 topsides. PyIRI (no storm model) gets the epoch's daily
F10.7, its 81-day mean, and the F10.7 of the epoch's R12. The thickness laws on the measured peak use
the epoch R12 and IRI's M(3000)F2. Checks: the binned configuration reproduces T1 to 1e-4; the
indices IRI reports equal `apf107.dat` exactly; a fresh end-to-end run reproduces `summary_epoch.json`.

Indices: epoch R12 - bucket is -0.9..+1.6 (5-95%). At solar maximum the daily F10.7 is 126-201 sfu
(5-95%, 12% above the 188 topside clamp) against 150-155 for the binned level; 9% of the profiles are
at Kp >= 4 (435 at solar maximum).

Solar maximum, median model/measured dB:

| band | IRI binned | IRI epoch, quiet | IRI epoch, foF2 storm | PyIRI binned | PyIRI daily F10.7 | PyIRI 81-day | p05 obs / IRI epoch / PyIRI daily (MHz) |
|---|---|---|---|---|---|---|---|
| 500-600 | +1.09 | +1.10 | +1.16 | -2.32 | -2.33 | -2.06 | 1.51 / 2.21 / 1.31 |
| 600-700 | +0.77 | +0.81 | +0.78 | -2.63 | -2.65 | -2.49 | 1.16 / 1.71 / 1.06 |
| 700-800 | +0.55 | +0.56 | +0.57 | -3.01 | -2.70 | -2.77 | 1.05 / 1.26 / 0.81 |
| 800-900 | +0.59 | +0.57 | +0.58 | -2.98 | -2.95 | -2.88 | 0.90 / 1.01 / 0.66 |

* The result does not depend on the indices: IRI-2020 moves by at most 0.07 dB, PyIRI by at most
  0.31 dB, and the fitted H0 is still 0.93 of the IRI law (0.91 at Kp < 3, 0.98 at Kp >= 4); 74% of
  profiles closer to it, 86% of the 435 at Kp >= 4.
* The storm model lowers IRI-2020's p05 at 800-900 km from 1.07 to 1.01 MHz (measured 0.90), and at
  Kp >= 4 its median from +0.58 to +0.46 dB.
* Solar minimum is where the indices matter. The mex derives IG12 from the user R12
  (IG12 = -12.35 + 1.468 R12 - 0.00268 R12^2), which for the bucket R12 15 is 9.1, while `ig_rz.dat`
  gives 17.3 for those epochs. The binned run's density there was 16% low (median); with epoch indices
  IRI-2020 goes from -0.39 to +0.32 dB (500-600 km, 23 profiles) and -1.27 to -0.56 dB (600-700 km,
  31 profiles). Too few profiles to settle anything, but it says the band-edge study's user-R12
  runs at solar minimum may be biased low (conservative for the lower edge): check against epoch
  runs.

## ISIS-1 at 900-1400 km (WP0 Phase 5, 2026-10-20)

`scripts/isis_high_altitude.py run|analyse` repeats the epoch comparison for the CRC profiles taken
with the spacecraft at 900-1400 km: 5,233 profiles (quality 0-7, Rz12 recorded), 4,659 at solar
maximum (Rz12 >= 100), 509 intermediate and 65 at solar minimum (Rz12 <= 30, all 1975-77, 63 of
them poleward of 60 deg). Same six IRI-2020 configurations (`topside_isis_epoch_iri.m`, 85 s, no
failures). Outputs in `data/reference/isis_high_altitude`.

Solar maximum, median model/measured dB, and the 5th percentiles (MHz):

| band | n | NeQuick epoch | IRI-2001 | IRI-2001 corr. | p05 obs | ensemble p05 | ens/obs | realised coverage | H0 fit / IRI law (median) |
|---|---|---|---|---|---|---|---|---|---|
| 900-1000 | 980 | +0.67 | +4.08 | +1.14 | 0.733 | 0.902 | 1.23 | 89.5 % | 0.96 |
| 1000-1100 | 947 | +0.67 | +5.25 | +1.55 | 0.601 | 0.810 | 1.35 | 87.1 % | 0.97 |
| 1100-1200 | 947 | +0.57 | +5.77 | +1.60 | 0.510 | 0.712 | 1.40 | 86.6 % | 0.97 |
| 1200-1300 | 900 | +0.44 | +6.52 | +1.64 | 0.460 | 0.672 | 1.46 | 86.4 % | 0.99 |
| 1300-1400 | 885 | +0.26 | +7.11 | +1.49 | 0.465 | 0.620 | 1.34 | 87.3 % | 0.96 |

The ensemble edge is INS-01a's construction (the minimum of the three topsides' 5th
percentiles; NeQuick gives it everywhere); realised coverage is the fraction of measured
fN_sat at or above it, 95 % if the low tail were right.

* **NeQuick's median holds** above 900 km (+0.3 to +0.7 dB, as at 700-900 km), the storm models move
  it by under 0.05 dB, and the fitted H0 sits on IRI-2020's law (median 0.96-0.99, p10-p90
  0.7-1.7). IRI-2001 is 4-7 dB high and grows with altitude; IRI-2001 corrected 1.1-1.6 dB.
* **Its low tail does not.** The modelled 5th percentile is 23-46 % above the measured one, against
  12 % at 800-900 km in the 2026-10-18 run and the 19 % Rev G reserved INS-01a's tolerance for.
  The edge INS-01a's method would set realises 86-90 % coverage instead of 95 %.
* **The overstatement is mostly polar.** At |lat| < 30 deg the ratio is 1.24 (88 % realised), at
  30-60 deg 1.12 (91 %), at 60-70 deg 1.28 (91 %), 70-80 deg 1.48 (78 %) and 80-90 deg 2.05 (64 %);
  the polar rows run 1.39-2.13 by band, the non-polar 1.13-1.29. Night (18-06 LT) is 1.32 against
  1.19 by day. The polar records are mostly CRC quality index 6-7 (73 % of index 7 are polar, 0.5 %
  of index 4), so latitude and scaling quality are confounded here: non-polar index 6-7 records
  give 1.25, non-polar index 4-5 1.19, polar index 6-7 1.93 (only 40 polar records are index 4-5).
  The polar-cap topside above 1000 km being less dense than IRI is plausible physics (IRI has no
  polar-wind depletion of H+), but this archive cannot separate it from the scaling.
* **Solar minimum** (63 polar, 2 non-polar records, 1975-77): measured p05 0.146 MHz against an
  ensemble 0.331 MHz (2.3 times), realised coverage 73 %, NeQuick median -0.66 dB. A polar-cap,
  shallow-minimum sample, too small to say anything about the deep minima or low latitudes.
* **Intermediate activity** (509): ratio 0.79-0.99, realised 95-98 %: the tail is right or low.

## Consequence

`l0: iri2020` (PyIRI) must not be used for band-edge or topside-shape results until its H0 is
corrected; `matlab_iri` is the reference. Raise the H0 step with the PyIRI developers.

Since 2026-10-18 campaigns choose the implementation with the `iri` tier's `model` option,
which defaults to PHaRLAP (`l0: {name: iri, f107: 120}`); `model: pyiri` remains for quick
runs that need no MATLAB.  The `baseline` campaign uses `model: pharlap`.
