# PyIRI against Fortran IRI-2020

The `l0: iri2020` tier is IRI-2020 through PyIRI 0.1.7 (Forsythe et al. 2023,
*Space Weather*; Servan-Schreiber et al. 2026, *Space Weather* 24,
e2025SW004853), an independent Python implementation of the International
Reference Ionosphere.  This note records how far it sits from the Fortran
IRI-2020 that PHaRLAP's `iri2020` mex runs (the `matlab_iri` tier), on the
quantities this work rests on.

Reproduce with

    python scripts/pyiri_validation.py OUT_DIR --f107 120

## Setup

* 576 points: latitude -80 to 80 in 20 degree steps, longitude -180 to 135 in
  45 degree steps, UT 0, 6, 12 and 18 h, on 2024-03-15 and 2024-06-15.
* Altitude 100 to 2000 km in 10 km steps.
* Both models with IRI-2020's default options and the NeQuick topside.  PyIRI:
  URSI foF2, SHU-2015 hmF2, ABT-2009 bottomside, F10.7 = 120.  PHaRLAP:
  `iri2020` with the NeQuick topside at a user R12 of 71.1, which the Davies
  relation of `tss.l0_ionosphere.base.r12_to_f107` maps to F10.7 = 120.  The
  pairing is not exact, because the two codes interpolate their coefficient
  maps in different solar indices.

## Results (2026-10-17)

| Quantity | PyIRI against IRI-2020 |
|---|---|
| NmF2 ratio | median 1.026; 10th to 90th percentile 1.009 to 1.054 |
| hmF2 difference | median -1.7 km; 10th to 90th percentile -6.9 to +4.3 km (10 km grid) |
| Ne ratio at 500 km | median 0.51; 0.46 to 0.57 |
| Ne ratio at 700 km | median 0.44; 0.41 to 0.46 |
| Ne ratio at 900 km | median 0.44; 0.42 to 0.46 |
| Ne / NmF2 at 900 km | PyIRI 0.020, IRI-2020 0.046 (medians) |
| 5th percentile of the plasma frequency at 500 / 700 / 900 km | PyIRI 1.07 / 0.57 / 0.39 MHz; IRI-2020 1.52 / 0.87 / 0.59 MHz |

The peak agrees.  The topside does not: PyIRI's density from 500 to 900 km is
about 0.44 to 0.51 of IRI-2020's, and the ratio is nearly constant from 700 km
up, so the difference is in the topside formulation near the peak rather than
in the scale height far above it.  Both are labelled NeQuick topsides.  PyIRI
builds an Epstein layer with the NeQuick-2 thickness law (`thickness_F2`,
`epstein_function_top_array` in PyIRI 0.1.7); where the factor of about two
comes from in IRI-2020's implementation has not yet been traced.

## What it means here

* The plasma frequency at the spacecraft is the quantity INS-01a was derived
  from, and PyIRI puts it about 35 per cent lower at 900 km in this sample.  A
  band-edge study run on PyIRI would therefore give a lower edge than the
  PHaRLAP IRI-2020 ensemble did.  The ISIS-1 comparison (Technical Volume
  Chapter 4) found the models overstate the low tail by about 19 per cent, so
  measurement sits between the two; which model is closer needs the archive,
  not a model-to-model comparison.
* A profile scored against SCI-05 and SCI-07 is scored against its own truth,
  so the topside difference does not bias those scores, but it changes the
  ionosphere the instrument is tested in, and the trade study inherits it.
* The `iri2016` Fortran wrapper is installed but not built (its first call
  tries to compile and needs ninja), so it could not be a third opinion.

Until the difference is traced, a result made on `iri2020` should name the
tier, and a band-edge or topside-shape result should be checked on
`matlab_iri` as well.
