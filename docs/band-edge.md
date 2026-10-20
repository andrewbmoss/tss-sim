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
