# The topside difference between PyIRI and Fortran IRI-2020

Campaign of 2026-10-17/18. Reproduce with

    python scripts/topside_campaign.py grid    OUT_DIR
    python scripts/topside_campaign.py isis    OUT_DIR --crc ../Dissertation/ISISValidation/isis1_ne.asc --t1 ../Dissertation/ISISValidation/t1_results.csv
    python scripts/topside_campaign.py analyse OUT_DIR

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
* PyIRI 0.1.7 (`sh_library.thickness_F2`): x = (k B2bot - 150)/100,
  v = 0.041163 x^2 - 0.183981 x + 1.424472, H0 = k B2bot / v. This is the NeQuick G (Galileo) step,
  which there follows a seasonal k clipped to 2..8. With the Coisson k (median 2.35), k B2bot is
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

## Consequence

`l0: iri2020` (PyIRI) must not be used for band-edge or topside-shape results until its H0 is
corrected; `matlab_iri` is the reference. Raise the H0 step with the PyIRI developers.

Since 2026-10-18 campaigns choose the implementation with the `iri` tier's `model` option,
which defaults to PHaRLAP (`l0: {name: iri, f107: 120}`); `model: pyiri` remains for quick
runs that need no MATLAB.  The `baseline` campaign uses `model: pharlap`.
