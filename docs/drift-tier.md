# WP5 Phase 3, part 1: curved-Earth facets, the sheet wrap and the drift tier

2026-11-25.  Plan: Dissertation `Plans\WP5_plan.md`, P3; decisions (a) to (h)
of 2026-11-25 09:01 (project note `claude/wp5-p3-decisions-2026-11-25.md`).
Phase 2's sequence path is in `docs/drift-sequence.md`.  Not yet done in
Phase 3: the diffuse halo, the Cramer-Rao harness, the study
`wp5_route_one`, Poole's comparator on sequences.

## 1. Rays on a sphere (`l1_propagation/sphere.py`)

In a medium whose index depends on the radius alone, every ray stays in the
plane through the Earth's centre, its start and its end, and Bouguer's law
holds along it.  So each scatterer's ray is a two-dimensional problem in its
own plane, solved by Earth flattening about the platform's radius:
`x = R psi`, `d = R ln(R/r)`, `mu'^2 = mu^2 (r/R)^2`, Snell's invariant by
Newton's method.  The phase path is exact; the group path is
`int (r/R)^2 dd / sqrt(mu'^2 - q^2)` (because `ds = (r/R) ds'` and
`mu = mu' R/r`), with no frequency difference.

Checks (`tests/unit/test_l1_sphere.py`):

| Check | Measured |
|---|---|
| Phase and group path against Bouguer's law integrated in radius (no flattening), 3.6, 5, 10 MHz, 18 points to 40 km off nadir | within 1e-9 m |
| Flattening in free space (index `(r/R)^2`, not 1) against the chord | 4e-10 m |
| Group path against `d(fP)/df` | 7e-4 m (the difference's truncation) |
| Strata and point rising together: phase path change | `-mu_sat w` to 1e-4 |
| Precision of the local frame on a sphere of 1e9 to 1e13 m | heights to 1e-9 m |

The frame is local (east, north, up) with the centre at `(0, 0, -r0)`;
radii and angles are formed from offsets, so a radius of 1e15 m is the flat
limit without loss of precision.  174 000 points take about 1 s cold, 0.5 s
warm-started.

## 2. Plates in three dimensions (decision (b))

`facet.plate_surface`: a rectangle `length_m` along `azimuth_deg` by
`width_m` across (default: the length, so a square facet), tilted by
`tilt_deg` along and `cross_tilt_deg` across; quarter-wavelength samples; the
centre as a zero-weight reference ray.  `plate_model` is `"2d"` or `"3d"`;
`"auto"` is 2-D on the straight track (M3 is scored there, decision (c)) and
3-D on the orbit.  Straight 3-D is the same engine on a sphere of 1e15 m;
orbit 2-D is Phase 2's flattened plate.  The four combinations split the
plate model's effect from the sphere's.

Checks (`tests/unit/test_l1_facet_3d.py`):

* The engine in the flat limit against `sheet_echo.echo_stratified` (a 30 degree
  dip, 10 degree rotation, curved, moving sheet, oblique track): within
  0.01 rad and 1 % in magnitude at nine times over 2 s.
* A tilted plate against the sheet it equals (dip = tilt, no rotation, no
  curvature, no taper), two independent point layouts: within 0.01 rad in
  the main lobe.
* The 3-D orbit plate against Phase 2's 2-D one: group delay within 5 cm, the
  Doppler at overflight within 1e-3, the phase different by a constant whose
  spread over the lobe is under 0.01 rad (the cross-track integral).
* A facet 3 km to the side of the track: the 3-D group path is the
  spherical chord to its centre to 1e-9.

**`|E|` of a 3-D plate is not the 2-D plate's** (the cross-track integral), so
`amplitude_db`, and the SNR a study takes from it, change meaning on the orbit.
The default width (a square) is a choice for Andrew to confirm.

## 3. The sheet wrap (P3 first item, D4's second half)

`sheet_echo.echo_track` reproduces `echo_free_space` and `echo_stratified`
exactly and adds the shortest ray's phase path.  On the straight track a sheet
(`kind: sheet`: `height_m`, `heading_deg`, `dip_deg`, `rotation_deg`, sizes,
velocities) is that model; on the orbit it is the 3-D engine with the sheet
laid out in the magnetic frame at its reference point, dip and declination
from IGRF (`field_model="igrf"`, provenance in `describe()`; `dip_deg` and
`declination_deg` override), curvature from the dipole formula.  Rows are
flagged `facet:<id>;sheet`, so L2 keys fading to the sheet's id unchanged.
The phase delay is taken in the cycle centred on -pi/2 from the shortest
ray's (the stationary phase of a two-dimensional minimum), so it does not jump
by a cycle during a pass.

**Exit met** (`tests/unit/test_l1_facet_sheet.py`): a noise-free sequence
reproduces the sheet model's phase to **2.0e-5 rad** (0.01 asked) and its
Doppler to **3.4e-7** (1e-3 asked), measured at 4208661 with the geometry-time
correction of section 4 (1.2e-10 and 1e-10 before it); `field_aligned_verify`
unchanged.  Rows are flagged `facet:<id>;sheet`, not `sheet:<id>` as the plan
wrote, so that L2's fading and facet handling apply unchanged: a deviation to
confirm.  IGRF is called with the spherical latitude of `tss.platform.orbits`
where it expects a geodetic one (up to 0.19 degrees at mid-latitudes).  The exit
track crosses the sheet's normal plane square on; along the normal's azimuth
the platform's own Doppler is about 100 Hz, which 45 Hz sweeps sample only
modulo 45 Hz (the echo set carries it exactly).

## 4. Geometry time and receive time (an L2 correction)

An echo received at `t` was reflected at `t - tau/2`, so the echo set L1
evaluates for the sweep's centre `T/2` describes what arrives at
`T/2 + tau/2`.  Phase 2's L2 turned each row's Doppler about `T/2`, which put
its phase half a delay early against its amplitude; route one read that as a
drift bias of `v^2 / (mu_sat c)`, **0.2 m/s at any range** (24.80 for 25 m/s
before, 24.997 after).  `SequenceContext.geometry_delay_fraction` (0.5 in a
sequence) fixes it, and the drift tier stamps each value with
`t_n + t_ref - tau/2`.  Phase 2's exit tests compare at the geometry time
now.  Their errors grew but stay far inside the criterion, because the record
centre is now about half a delay from the time the echo set describes and the
fixed Doppler is extrapolated over it (second order): straight 5.4e-4 rad over
the window (8.4e-5 over the main lobe; 1.3e-9 before), Doppler 5.9e-5; orbit
(3-D plate since part 1) 1.9e-4 rad.  This changes Phase 2's pushed
behaviour, for Andrew to accept.

## 5. The drift tier (`l3_signal/drift.py`, `l3/drift`)

Route one at signal level, one channel: detect (sum of spectra over the
sequence, peaks within 25 dB of the strongest), gate out the first delay
(decision (f); `blanked_fraction` 0.10 at 340 km to 0.16 at 530 km of group
path), one complex value per sweep at the gate's centre converted at the
frequency transmitted for it (`f_eff`, D8 condition 2), the platform Doppler
for a local-down wave normal removed, then `estimate_route_one` (reused).
Products in the contract's names (`docs/drift-mode.md`); `drift_sigma_mps`
is NaN until the Cramer-Rao harness and `not_done` says so.  The N(h) prior
is an explicit input (`OnboardPrior`); `run_sequence` has no ionogram chain to
supply one, so tests pass the truth profile and the product says
`prior_source: truth` (advisor, 2026-11-25).

Noise-free (`tests/unit/test_l3_drift.py`; 5 MHz, `mu_f` 0.5, 2 km plate):

| Case | Truth | Drift tier |
|---|---|---|
| Straight, overflight | 25 m/s | 24.997 |
| Straight, facet 2 km along track | -60 m/s | -60.35 |
| Orbit, 3-D plate 300 km below 600 km, exponential profile (not the L0 profile) | 50 m/s | 49.996 |
| Extended plate (80 km, tapered) | declined | `extended` |

All within D5's 1 m/s.  The 0.35 m/s at 2 km along track is the power-envelope
fit, not the chain: the same estimator on the echo set's own field gives the
same value, and on the plate model's field without the spreading term gives
0.004; a 3e-4 change in the envelope's shape moves the four-parameter
`sinc^2` fit's centre by 3 ms when the peak falls between sweeps.  Noise
(sigma of `t_pk` 23 to 54 ms in the budget) will swamp it; it is reported, not
fixed.

### The effective speed: a correction to Phase 2

Phase 2 said route one would read the drift about 4 % wrong on the orbit
without `v_eff^2 = v^2 P'/P'_flat`.  Measured: noise-free, to first order, it
does not.  Route one's drift is
`(c / 2 f mu) [<f_D> + S (t_bar - t_pk)]`, with `t_bar` the power-weighted
mean time of the Doppler samples; `v` enters only through the last term,
which vanishes when the lobe is symmetric.  Noise-free on the orbit the drift
with and without `v_eff` differs by 0.0004 m/s.  With the envelope made
asymmetric (a per-sweep gain) it differs by 0.06 m/s (a +/-50 % ramp), 0.06
(a 1.3 Hz fade) and 0.28 (a factor of 2 step), while those gains move the
drift itself by up to 28 m/s: second order against the fading error.  What
the curvature moves to first order is the **slope ratio**, the compact-facet
test, by `(v_eff/v)^2`: 0.965 with `v_eff`, 0.922 without (0.956).  So
decision (d) (the facet's true height from the onboard profile) matters to
the classification first and to the drift only at second order.

## 6. What is not yet known or done

* **Noise.**  SNR per sweep is the budget's definition: the per-sweep
  complex value at the envelope's peak over its noise, which is the value
  coherent over the gate (checked: 10.07 dB measured for 10 nominal).
  Straight track, 25 m/s, +/-2.5 envelope widths (209 sweeps), 200 seeds,
  the beat given so the estimator is isolated (robust sd = 1.4826 MAD):

  | SNR/sweep | formed | slope ratio median (p16, p84) | drift median, robust sd | sd `t_pk` | budget (overflight, Doppler; sigma_t0) |
  |---|---|---|---|---|---|
  | 20 dB | 200/200 | 0.986 (0.96, 1.01) | 25.13, 1.31 m/s | 11.0 ms | 1.07, 1.11 m/s; 9.4 ms |
  | 15 dB | 194/200 | 0.990 (0.91, 1.05) | 25.13, 2.35 m/s | 19.9 ms | |
  | 10 dB | 92/200 | 0.884 (0.49, 1.20) | 24.46, 6.67 m/s | 34.6 ms | 3.43, 3.51 m/s; 30 ms |

  The same estimator on an ideal field (a `sinc^2` envelope of width W, the
  ephemeris Doppler, white noise; no L1, L2 or tier) gives at 10 dB 87/200
  formed, slope ratio 0.909, drift sd 6.13 m/s, `t_pk` 28.5 ms; on the
  tier's own noise-free field with white noise added, 88/200, 0.881, 6.19.
  So the 10 dB behaviour is the estimator's, not the chain's.  A first run
  of 40 seeds (1000 to 1039) gave 11/40 formed and was an unlucky draw; the
  same 40 seeds give the same result again.

  The right comparator is the budget's overflight and Doppler terms (RSS 1.5
  m/s at 20 dB, 4.9 at 10 dB), not its full RSS (10.6 at 10 dB), which adds
  the oscillator and range-scale terms the simulation does not include.  An
  earlier first look (20 dB, drifts 12 to 38 m/s) was a bug: the envelope fit
  (scipy `least_squares`, absolute tolerances) stopped early on powers of
  1e-12 W.  The field is now scaled to a unit peak before the estimator
  (`power_peak_db` keeps the scale), and a test pins `t_pk` to the
  `route_one` Monte Carlo within a factor of 1.6.

  **At 10 dB route one's compact-facet test declines more than half of the
  compact facets** (108 of 200 declined; in the 40-seed run every decline
  was on the slope, none on the contrast): the free, power-weighted
  Doppler slope has a median 0.88 of the ephemeris slope and a spread
  (0.49 to 1.20) wider than the tolerance 0.8 to 1.2.  A likely cause, not
  tested: the Doppler history keeps samples down to 0.1 of the (noisy)
  peak power, which at 10 dB is the noise level, and their random Doppler
  pulls the weighted slope towards zero.  No earlier study ran the slope
  test with noise: the parameter-level `route_one` study sizes `t_pk` and
  the Doppler only, and `c3_verify` is noise-free.  So **the budget's 10 dB operating point overstates availability** by about
  a factor of two, and the formed drifts scatter by 6.7 m/s against the
  budget's overflight and Doppler RSS of 4.9.  `t_pk` (34.6 ms) is near the
  budget's 30 ms, so most of the excess appears to be in the zero crossing
  (not separated: `t_pk` and `t_zc` share the noise).  The estimator is
  pinned by `wp5_baseline` and is not changed here; a slope test with a
  higher power cut, or the ephemeris slope as a prior, is for the study.
* **Detection at 10 dB.**  The part-1 summed-spectrum detector found 0 of 5
  at 10 dB: its statistic is the spectrum summed over all 91 sweeps, of which
  the facet's lobe fills about 15, against a 13 dB threshold.  Replaced by
  the windowed detector (section 8).
* Contract products not formed yet: per-sweep `range_m`, `power_db` and
  `doppler_hz`; `snr_db`, `resolution_mps`, `v_unambiguous_mps`,
  `range_debiased_m`, `los_rate_sigma_mps`, `tau_c_s`, `doppler_spread_hz`,
  `drift_sigma_mps` (NaN).
* The facet tier keeps state between calls (the latched plate azimuth, the
  sheet's field direction, the surfaces and the Newton warm start), so a study
  must build a tier per scene.
* The halo, the Cramer-Rao harness, the study `wp5_route_one`, Poole's
  comparator.  The onboard prior is section 7; its full N(h) (L4's
  retrieval rather than an assumed scale height) is later work.
* Noise so far is white, added to noise-free records or L2's thermal noise.
  The negative-beat half of an L2 drift record with thermal noise on was
  checked flat (-103.0 dB in both halves, 0 to 1050 km, no notch: the
  drift-mode records carry no leakage high-pass), but sky noise, clutter and
  quantisation in drift mode are untested.
* Noise-free on the Chapman orbit scene with the truth prior the slope ratio
  is 0.916 (0.965 on the exponential scene of part 1): 0.084 from the
  decline edge and unexplained.

## 7. The onboard prior (part 2, Option A, decisions (i) and (k))

`onboard_prior_from_ionogram(ionogram, platform)`: the WP4 tracker's `f_N`
at the spacecraft (`flags["tracker"]["priors"]["f_n_hz"]`), accepted only when
below the lowest detected trace frequency (L4's rule), and an exponential
profile below the spacecraft of an **assumed** scale height (200 km
default, flagged `scale_height_source: assumed`).  No valid `f_N` gives no
prior, so no drift and `NO_PRIOR` (counted as an M3 failure, decision (k)).
The report travels in the drift product's `flags["prior"]`.
`Pipeline.run_prior_ionogram` makes the ionogram from the drift scene's own
L0, through a given L1/L2/L3 (here `analytic`, the baseline instrument at
600 km, `baseline` L3 with the tracker, the A1 sweep), on random streams of
its own; `run_configuration` and `run_sequence` are unchanged.

Measured on the Chapman L0 (f_N 1.8224 MHz at 600 km, topside scale height
100 km), a 3-D plate at 430 km drifting 50 m/s, 5 MHz, noise-free sequence
(`tests/unit/test_l3_drift_prior.py`):

| Prior | f_N, MHz (method) | drift, m/s | facet depth, km | slope ratio |
|---|---|---|---|---|
| truth (the L0 profile) | 1.8224 | 49.982 | 169.9 | 0.916 |
| onboard, clutter off, H 200 km | 1.8203 (trace_start) | 49.992 | 189.8 | 0.918 |
| onboard, clutter off, H 100 km | 1.8203 | 49.992 | 169.0 | 0.916 |
| onboard, clutter on, H 200 km | 1.8658 (resonance) | 50.187 | 188.5 | 0.917 |

`f_N` sets `mu_sat` and so scales the drift: the clutter-on prior pass
took a resonance 2.4 % high and moved the drift 0.4 %.  The scale height
moves the inferred depth (20 km for a factor of two) and `v_eff`, not the
drift.  The clutter-on row is one realisation; the spread of the tracker's
`f_N` and the NO_PRIOR rate over seeds are not yet measured.  The
tracker's O and X labels swap near 5 MHz on this scene (one seed has the O
trace in the X column), which matters to the range gate (section 8) and
not here.

Fixed on the way: `effective_speed_ratio` halved an evanescent bracket
until it propagated and could cut below the answer (108 km returned for
170 km); it now bisects to the reflection level first, and an echo whose
`P'` the prior cannot reach is flagged `PRIOR_REFLECTS_ABOVE`.

## 8. The windowed detector (part 2, decision (j))

Per sweep, the samples after `2 max_range / c` are Hann windowed and
transformed without padding (one bin, 1.46 km, is a range cell); the power
is summed over `window_sweeps` consecutive sweeps, sliding by one.  The
noise reference is the negative-beat half (median over ln 2, about 1.4e5
values per sequence).  Under noise a cell is Gamma(L) distributed, and the
threshold is its quantile for `pfa` per cell (1e-4, WP4's value).  The
strongest `max_echoes` (2) range peaks, four bins apart and within
`dynamic_db` of the strongest, are echoes; the beat is refined on the
window that peaked.  Window: 21 sweeps (0.47 s), the half-power width of a
4 km facet's envelope at 530 km (0.463 s) and half that of a 2 km facet's
(0.926 s); the facet's length is not known onboard.  `range_gate_from_ionogram` bounds the search
by the trace's virtual range at the frequency (largest over the O and X
columns, because the labels can swap) plus 20 km.

Straight track, 25 m/s, 209 sweeps, SNR per sweep as in section 6:

| Case | pfa per cell | L = 11 | L = 21 | L = 41 |
|---|---|---|---|---|
| Noise only, cells over / cells (200 seq.) | 1e-4 | 1.01e-4 | 1.07e-4 | 1.15e-4 |
| Noise only, echoes declared per sequence | 1e-4 | 1.98 | 1.91 | 1.55 |
| Noise only, drifts formed on them (200 seq.) | 1e-4 | 8 | 15 | 12 |
| Noise only, cells over / cells (100 seq.) | 1e-7 | 0 | 2.3e-7 | 2.5e-7 |
| Noise only, echoes declared per sequence | 1e-7 | 0.00 | 0.01 | 0.02 |
| 10 dB, facet detected (50 seeds) | 1e-4 | 50 | 50 | 50 |
| 10 dB, facet detected (50 seeds) | 1e-7 | 50 | 50 | 50 |
| 7 dB, facet detected (50 seeds) | 1e-7 | 46 | 50 | 50 |
| 4 dB, facet detected (50 seeds) | 1e-7 | 9 | 20 | 34 |

End to end at 10 dB, L = 21, 200 seeds, nothing given: the facet is
detected in 200, a drift formed in 96 and declined in 104 (route one's
slope test, section 6).  At 1e-4 per cell the second echo slot holds a
false echo in 196 of 200 sequences and 8 drifts form on them; at 1e-7, one
false echo and none formed.

**At 1e-4 per cell the stated rate is met, but a sequence holds 1.2e5 to
1.4e5 cells, so nearly every sequence declares a false echo and about 1 in
13 forms a false drift.**  1e-4 is WP4's rate for an ionogram cell, where
the tracker then integrates along the trace; here nothing follows the
detector but route one.  A rate stated per sequence (1e-7 per cell is about
0.01 false echoes per sequence) costs about 1 dB of threshold and no
detections at 10 dB.  The default stays 1e-4 as decided; the choice is
put back to Andrew.

The range gate on the Chapman orbit scene (section 7): the trace gives
359 km at 5 MHz, gate 379 km; the facet is at 215.65 km.  The segment
grows from 1438 to 1861 samples (+1.1 dB) and the range bins fall from 705
to 319.  At 10 dB (30 seeds) the facet's excess over threshold rises from
2.96 to 4.45 dB and false echoes fall from 25 to 14 per 30 sequences.
The gated drift equals the ungated one.  The gate is tested on this one
consistent scene only; in the exponential straight-track scenes the
ionogram's L0 is not the facet's medium, and the gate is not used.
