# WP5 Phase 3, part 1: curved-Earth facets, the sheet wrap and the drift tier

2026-11-25.  Plan: Dissertation `Plans\WP5_plan.md`, P3; decisions (a) to (h)
of 2026-11-25 09:01 (project note `claude/wp5-p3-decisions-2026-11-25.md`).
Phase 2's sequence path is in `docs/drift-sequence.md`.  Part 2 (the same
night; decisions (i) to (l), Option A): noise (section 6), the onboard prior
(7), the windowed detector and range gate (8), the Cramer-Rao bound (9).
Part 3 (decisions (m) to (u), 2026-11-27): the false-alarm rate stated per
sequence (10), the joint estimators (11), the fading outliers (12), the power
cut (13), the diffuse halo (14), Poole's velocity as the comparator (15).
Not yet done in Phase 3: the study `wp5_route_one`.

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
`prior_source: truth` (advisor, 2026-11-25).  (As built in part 1.  Since
part 2 the default detector is the windowed one of section 8, the onboard
prior is section 7, and `drift_sigma_mps` stays NaN: the bound of section
9 is not the estimator's error.)

Noise-free (`tests/unit/test_l3_drift.py`; 5 MHz, `mu_f` 0.5, 2 km plate):

| Case | Truth | Part 1 (23a4270) | HEAD, summed detector | HEAD, windowed (default) |
|---|---|---|---|---|
| Straight, overflight | 25 m/s | 24.997 | 25.029 | 25.066 |
| Straight, facet 2 km along track | -60 m/s | -60.35 | -59.901 | -59.826 |
| Orbit, 3-D plate 300 km below 600 km, exponential profile (not the L0 profile) | 50 m/s | 49.996 | 49.999 | 50.015 |
| Extended plate (80 km, tapered) | declined | `extended` | `extended` | `extended` |

All within D5's 1 m/s.  **Corrected 2026-11-25 night.**  Part 1 put the
0.35 m/s at 2 km along track on the envelope fit's sensitivity to shape (a
3e-4 change in shape moving the four-parameter `sinc^2` fit's centre by
3 ms).  Most of it was the fit's absolute tolerances on physical-scale
powers (section 6): the same estimator on the plate model's own field,
which is near unit scale, gave 0.004, and with the field scaled to a unit
peak the tier gives 0.10 (summed beat) or 0.17 (windowed beat).  What is
left moves with the beat the detector reports (it sets the gate's first
sample): the noise-free `t_pk` shifts by about 0.5 ms between the two
detectors, 0.06 m/s.  The parameter-level studies are not affected:
`c3_verify` hands the estimator fields of peak power 0.98 to 1.0, and
scaling them to a unit peak changes no drift at 1e-4 m/s (21 cases
checked; the other studies were not checked).

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
  (scipy `least_squares`, absolute tolerances) stopped early on
  physical-scale powers (a peak of about -66 dB in the tier's units on
  these scenes).  The field is now scaled to a unit peak before the estimator
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
  `drift_sigma_mps` (NaN).  The last is formed by the joint estimators of
  section 11.
* The facet tier keeps state between calls (the latched plate azimuth, the
  sheet's field direction, the surfaces and the Newton warm start), so a study
  must build a tier per scene.
* The study `wp5_route_one`.  (The harness is section 9, the halo section 14,
  Poole's comparator section 15.)  The onboard prior is section 7; its full N(h) (L4's
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
drift.  Over 20 seeds of the prior pass (platform at 600 km on this
Chapman scene): clutter off, NO_PRIOR 0/20 and the tracker's `f_N` -0.11 %
in every seed (`trace_start`); clutter on, NO_PRIOR 0/20 and +2.38 % in
every seed (`resonance`): a systematic error of the method on this
scene, not noise (0.4 % of drift through `mu_sat` at 5 MHz).  Other
scenes and other `f_N` are not measured.  The tracker's O and X labels
swap near 5 MHz on this scene, which matters to the range gate (section
8) and not here.

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
detections at 10 dB.  The default stayed 1e-4 as decided and the choice was
put back to Andrew, who set the rate per sequence (section 10).

The range gate takes the larger of the O and X columns' virtual ranges and
is a bound only when the O trace is among the columns that bracket the
frequency (`columns_used` in its report): over 20 prior-pass seeds at
5 MHz, both columns bracketed in 10 (clutter off) and 14 (on), only column
1 in 10 and 5, only column 0 in 0 and 1.  The single-column cases on this
scene carry the O trace under the X label (the values match), but a true X
trace alone would give a gate tighter than the O bound by the O-X
separation (about 11 km at 5 MHz here, so inside the 20 km margin on this
scene) less the margin.

On the Chapman orbit scene (section 7): the trace gives 359 km at 5 MHz,
gate 379 km; the facet is at 215.65 km.  The segment
grows from 1438 to 1861 samples (+1.1 dB) and the range bins fall from 705
to 319.  At 10 dB (30 seeds) the facet's excess over threshold rises from
2.96 to 4.45 dB and false echoes fall from 25 to 14 per 30 sequences.
The gated drift equals the ungated one.  The gate is tested on this one
consistent scene only; in the exponential straight-track scenes the
ionogram's L0 is not the facet's medium, and the gate is not used.

## 9. The Cramer-Rao bound (part 2, `l3_signal/novel/drift_crb.py`)

The per-sweep samples of one compact facet, `y = A sinc(x) e^{j phi} g + w`,
with route one's envelope (`x = 2 X (t - t0) / W`) and Doppler law
(`phi = phi0 + 2 pi S [kappa v (t - t0) - (t - t0)^2 / 2]`, `kappa v =
t_zc - t0`).  `fisher_route_one` gives the deterministic matrix (`g = 1`,
parameters `A, phi0, t0, v, W`) and the zero-mean Gaussian one (Rayleigh
`g` of correlation `exp(-|tau|/tau_c)`, parameters `P, t0, v, W`).  SNR per
sweep `A^2`, the budget's definition.

Checks (`tests/unit/test_l3_drift_crb.py`): the Doppler-sample model of the
degeneracy note gives 0.4240 m/s with `t0` known and a singular matrix
with it free; with `t0` known the drift bound is Stein's FDOA accuracy as
Moss (2014, eqs. 2.25 to 2.28) gives it, for one noisy channel, to 1e-6;
for a flat envelope Stein's form is the budget's `sqrt(6) / (2 pi T
sqrt(rho N))`; a flat envelope (phase only) leaves `t0` and `v`
inseparable (condition 1e16); both matrices match finite differences of
the full model (the phase terms left to `phi0` kept in) to 1e-4 and 1e-3.
Moss (2014) is Andrew's NPS thesis; only its Section 2.3 (Stein's TDOA and
FDOA accuracies) is used.

The straight scene of section 6 (4.988 MHz, P' 529.8 km, mu_sat 0.954,
W 0.926 s, 209 sweeps):

| SNR/sweep | CRB `sigma_v`, `t0` free (known) | CRB `sigma_t0` | corr(t0, v) | estimator `sigma_v` | estimator `sigma_t_pk` | (CRB / est.)^2 | CRB with fading, tau_c 0.5 s (0.1, 2.0) |
|---|---|---|---|---|---|---|---|
| 20 dB | 0.735 (0.108) m/s | 6.36 ms | -0.989 | 1.31 m/s | 11.0 ms | 0.31 | 4.30 (10.7, 2.40) m/s |
| 10 dB | 2.33 (0.342) m/s | 20.1 ms | -0.989 | 6.67 m/s (92/200 formed) | 34.6 ms | 0.12 | 6.76 (15.2, 4.55) m/s |

The bound limits the RMS error, so the estimator's RMS over the formed
drifts (a rerun with the beat from the windowed detector, so the robust
figures differ from the table in the second digit): 20 dB, RMS 1.19 m/s,
bias +0.11, robust sd 1.29, `t_pk` RMS 9.7 ms (200/200 formed); 10 dB, RMS
6.42 m/s, bias +0.10, robust sd 6.86, `t_pk` RMS 31.6 ms (93/200 formed).
In RMS the estimator is 1.6 times the bound at 20 dB and 2.8 times at
10 dB ((CRB / RMS)^2 0.38 and 0.13).

`W` free or known gives the same bound (the envelope's width and centre
decouple for a window centred on the peak).  **The estimator is not
efficient**: 1.6 to 1.8 times the bound at 20 dB and 2.8 to 2.9 times at
10 dB, against the plan's exit of 20 % (or the shortfall explained, TV
risk 3).  The
shortfall is in both halves: `t_pk` comes from a fit to the power alone
(the budget's Monte Carlo of that fit gives 9.4 and 30 ms against the
bound's 6.4 and 20), and the zero crossing from a power-weighted mean of
phase steps.  The pilot of the motion-updated formulation reported 0.09
to 0.46 of the bound (TV Section 8.22); this is 0.12 to 0.31.  The bound
is not the budget: the budget adds the overflight and Doppler terms in
quadrature (4.9 m/s at 10 dB), while the bound treats them jointly
(their correlation is -0.989).  Not in the bound: the oscillator, the range
scale, the ephemeris, a Rician `K`; the 10 dB estimator figure is over the
formed drifts only.  With the pinned estimator `drift_sigma_mps` in the
product stays NaN: the bound is not the estimator's error (section 11 has the
joint estimators' own sigma).

## 10. The false-alarm rate stated per sequence (part 3, decisions (m) and (q))

Andrew, 2026-11-27: the detector's rate is stated per sequence, 1e-2 by
default (`pfa_per_sequence`), and the per-cell rate is that over the number
of (window, range bin) cells in the sequence.  This is a union bound on the
chance that any cell of a noise-only sequence crosses the threshold, and it
is conservative because overlapping windows share sweeps.  A range gate
leaves fewer cells, so the rate per cell rises and the rate per sequence
stays.  `pfa`, per cell as WP4 states it, still overrides it when given.  The
report carries `pfa_per_cell` and `pfa_per_sequence` (None when `pfa` was
given).

It differs from WP4's per-cell rate because WP4's tracker integrates along
the trace after the detector, so a false cell rarely reaches a product.  In
drift mode nothing follows the detector but route one.

Measured at the default (the straight scene of section 6, 209 sweeps, 200
noise-only sequences per window; the noise-only seeds are 20000 to 20199):

| Window L | cells | rate per cell | threshold over the mean | sequences with any cell over | cells over | echoes declared | drifts formed |
|---|---|---|---|---|---|---|---|
| 11 | 140295 | 7.1e-8 | 38.0 | 2/200 | 3 | 2 | 0 |
| 21 (default) | 133245 | 7.5e-8 | 54.5 | 2/200 | 4 | 2 | 0 |
| 41 | 119145 | 8.4e-8 | 83.7 | 1/200 | 2 | 1 | 0 |

At 1e-4 per cell, L = 21 declared 1.91 echoes per noise-only sequence and
formed 15 false drifts in 200 (section 8).

Detection and formation with the facet present (nothing given; 200 sequences
at 10 dB per sweep, seeds from 30000; 50 at 7 and 4 dB, seeds from 40000 and
50000):

| SNR/sweep | L | facet detected | false echoes | drifts formed on the facet | on a false echo |
|---|---|---|---|---|---|
| 10 dB | 11 / 21 / 41 | 200 / 200 / 200 of 200 | 1 / 2 / 0 | 71 / 80 / 77 | 0 / 1 / 0 |
| 7 dB | 11 / 21 / 41 | 47 / 50 / 50 of 50 | 1 / 0 / 0 | 8 / 9 / 9 | 0 / 0 / 0 |
| 4 dB | 11 / 21 / 41 | 7 / 20 / 34 of 50 | 0 / 1 / 1 | 0 / 2 / 2 | 0 / 0 / 0 |

The new threshold does not change the facet's estimate: on the same 200 seeds
at 10 dB with L = 21, both rates form the same 80 drifts with identical
values; at 1e-4 per cell there are 197 false echoes and 6 false drifts, at
the default 2 and 1.  The 80 of 200 formed on these seeds stands beside 96 of
200 on section 8's seeds and 92 of 200 with the beat given (section 6): at
10 dB route one's slope test forms about 40 to 48 % of compact facets on
this scene.

Tests (`tests/unit/test_l3_drift.py`): noise-only sequences rarely declare
any cell; the per-cell rate follows the number of cells, gated and not; `pfa`
overrides; the facet is found at 10 dB with the second echo slot empty.

## 11. The joint estimators (part 3 step 2, decisions (o), (r), (s) and (t))

`l3_signal/novel/route_one_joint.py`; the tier's `estimator` is `"pinned"`
(default), `"joint"` or `"joint_fading"`.  The deterministic fit was decision
(r); the fading-aware fit was added by Andrew on 2026-11-27.

* **Model** (section 9's): `y_n = c s_n g_n + w_n`, `s_n` the `sinc`
  envelope and route one's Doppler law with the ephemeris slope.
* **`joint`**: `g_n = 1`.  Least squares on the complex values, `c` solved in
  closed form, a search in `t0`, `v`, `W` from the pinned estimate.  Sigma
  from the Jacobian of the projected residual.
* **`joint_fading`**: `g` zero-mean complex Gaussian with correlation
  `exp(-|dt|/tau_c)`, so `y ~ CN(0, P D R D^H + sigma^2 I)`; `P`, `t0`, `v`,
  `W`, `tau_c` fitted, `sigma^2` given.  `R^-1` is tridiagonal, so the
  likelihood costs `O(N)` (it matches the dense form to 1e-9).  Sigma from the
  observed information over the parameters off their bounds.
* **Gate** (decision (s)): the pinned free-slope and contrast tests decide, so
  all three form a drift on the same echoes and availability is unchanged.
  The pinned answer is kept (`drift_pinned_mps`, `t_pk_pinned_s`).  A joint
  fit that fails leaves the pinned drift, flagged `JOINT_FAILED_PINNED_USED`.
* **Noise for the fading fit**: the windowed detector's negative-beat power
  per bin over `sum(hann^2)` is the sample variance; times the number of
  gated samples is the variance of the per-sweep value; then the unit-peak
  scaling.  The tier's value is 1.00 +/- 0.02 of the injected one (5
  sequences at 20 dB).  The summed detector has no noise reference, so
  `joint_fading` then falls back to the pinned drift.
* **`tau_c_fit_s`** is the fitted time of the exponential correlation, with
  `TAU_FIT_LOWER_BOUND` above half the dwell.  `tau_c_s` is left for the
  lag-one definition of `docs/drift-mode.md` section 4, not formed yet
  (Andrew, 2026-11-27).

### The noise floor of the fading fit

The envelope is not exactly `sinc^2`: on the noise-free `c3_verify` fields the
deterministic fit leaves a residual of 1.6e-4 to 2.7e-3 of the peak power.  A
fading likelihood told that less than that is noise explains the misfit as
fading and moves `t0` and the width:

| noise given, of the peak | untilted (case A) over 1 m/s | tilted (case B) over 1 % of the bias |
|---|---|---|
| 1e-6 | 45/56 (up to 241 m/s) | 30/32 |
| 1e-4 | 0/56 | 1/32 (30 %) |
| 1e-3 | 0/56 | 0/32 |
| 3e-3 (`MODEL_ERROR_FLOOR`) | 0/56 (max 0.04 m/s) | 0/32 (max 0.16 %) |
| 1e-2 | 0/56 | 0/32 |

The fit therefore does not take `sigma^2` below 3e-3 of the largest sample
power (accepted by Andrew, 2026-11-27).  Under fading the floor binds in up to 21 % of fits at 20 dB and never
at 10 dB; the RMS errors below change by 0.01 of the bound or less from no
floor to this one (1e-2 raises one case from 1.07 to 1.16).

### The search

The fading likelihood has several minima that differ in width and time.  The
search is a grid (`t0` over the window in 0.1 s steps, `v` over +/-300 m/s,
six widths, three `tau_c`) and a local fit from the six best distinct cells at
each width and from the pinned estimate.  How often each setting missed the
lowest likelihood that any setting or a fit started at the truth reached (100
seeds per case, formed echoes only):

| Setting | 20 dB, tau_c 0.5 s | 10 dB, 0.5 s | 20 dB, 2 s | 10 dB, 2 s | s per fit |
|---|---|---|---|---|---|
| pinned estimate alone | 68 % | 50 % | 38 % | 29 % | 0.03 |
| six best cells, any width | 46 % | 24 % | 35 % | 13 % | 1.2 |
| best 1 per width | 18 % | 8 % | 10 % | 2 % | 1.2 |
| best 2 per width | 8 % | 0 | 3 % | 0 | 1.3 |
| best 6 per width (default) | 0 | 0 | 0 | 0 | 1.8 |
| 6 per width, +/-600 m/s | 0 | 0 | 0 | 0 | 3.0 |

A grid twice as fine in `t0`, `v` and width with seven `tau_c` and four starts
per width gave the same errors as the default on 24 seeds per case, at about
20 s per fit.  (Times: 14 processes at once.)

### Noise-free, through the tier

| Case | Truth | pinned | joint | joint_fading |
|---|---|---|---|---|
| Straight, overflight | 25 m/s | 25.066 | 25.090 | 25.100 |
| Straight, facet 2 km along track | -60 m/s | -59.826 | -59.743 | -59.575 |
| Orbit, 3-D plate | 50 m/s | 50.015 | 50.010 | 49.976 |
| Extended plate (80 km) | declined | declined | declined | declined |

At parameter level both fits read all 56 untilted `c3_verify` facets within
1 m/s and all 32 tilted ones within 1 % of the predicted bias.

### Parameter level against the bounds

Section 9's scene (209 sweeps, drift 25 m/s), 200 seeds per case from 10000,
echoes the pinned test accepts.  RMS drift error over the matching bound
(deterministic without fading, Rayleigh-fading with it):

| SNR, tau_c | formed | bound, m/s | pinned | joint | joint_fading |
|---|---|---|---|---|---|
| 20 dB, none | 200 | 0.735 | 1.60 | 0.96 | 1.12 |
| 15 dB, none | 194 | 1.307 | 1.79 | 0.97 | 1.06 |
| 10 dB, none | 81 | 2.325 | 2.46 | 0.91 | 0.98 |
| 20 dB, 0.5 s | 122 | 4.298 | 5.22 | 6.50 | 1.08 |
| 10 dB, 0.5 s | 72 | 6.759 | 3.83 | 6.50 | 1.78 |
| 20 dB, 2 s | 176 | 2.403 | 8.80 | 9.48 | 1.37 |
| 10 dB, 2 s | 85 | 4.554 | 4.32 | 4.22 | 2.90 |

The fading fit by its outliers (errors beyond three bounds):

| SNR, tau_c | RMS | robust sd (1.4826 MAD) | outliers | RMS without them | sigma / robust sd | within 1, 2 sigma | tau_c fitted / true, median (16-84 %) |
|---|---|---|---|---|---|---|---|
| 20 dB, 0.5 s | 1.08 | 0.98 | 0.8 % | 1.05 | 1.06 | 0.69, 0.94 | 1.04 (0.77-1.59) |
| 10 dB, 0.5 s | 1.78 | 1.11 | 5.6 % | 1.25 | 0.81 | 0.58, 0.88 | 1.20 (0.76-2.00) |
| 20 dB, 2 s | 1.37 | 1.08 | 4.0 % | 1.06 | 1.04 | 0.69, 0.92 | 0.95 (0.56-1.71) |
| 10 dB, 2 s | 2.90 | 1.03 | 8.2 % | 1.03 | 1.02 | 0.62, 0.86 | 1.14 (0.56-2.35) |

In none of the 200-seed cases did a fit started at the truth reach a lower
likelihood than the search's answer, so these outliers are not search
failures: the likelihood's lowest minimum found lies far from the truth (a
threshold effect).  The cloud prototype (not merged) reported 1.14 and 1.60
at 0.5 s and 1.13 to 3.22 at 2 s on 100 seeds, with its grid placed about the
true `t0` and drift.

`joint`'s sigma matches its spread without fading (1.05, 1.03, 1.09 of the
standard deviation at 20, 15, 10 dB).  **Under fading it is 0.05 to 0.11 of
the spread**, and its RMS error is worse than the pinned estimator's in three
of the four cases: it must not be quoted on a fading echo.

### Signal level, no fading

The straight scene through L1, L2 and the tier with white noise, 200 seeds per
SNR from 60000:

| SNR | formed | pinned | joint | joint_fading |
|---|---|---|---|---|
| 20 dB | 200/200 | 1.24 m/s (1.69) | 0.81 (1.10), sigma/sd 0.95 | 1.00 (1.37), sigma/sd 1.03 |
| 10 dB | 90/200 | 6.58 m/s (2.83) | 2.58 (1.11), sigma/sd 0.90 | 2.83 (1.22), sigma/sd 0.96 |

(In brackets: over the deterministic bound.)  No fit failed.

**Run time** (one process, 209 sweeps): the tier takes 0.02 s per sequence
with `pinned`, 0.03 s with `joint` and about 1.5 s with `joint_fading`, per
echo.  Longer windows cost more.

### The exit

Within 20 % of the matching bound, or the shortfall explained:

* **`joint`, without fading: met** (0.91 to 0.97 at parameter level, 1.10 and
  1.11 at signal level).
* **`joint_fading`, under fading: met at 20 dB, tau_c 0.5 s (1.08).**  Not met
  in RMS elsewhere (1.37, 1.78, 2.90): outliers in 4 to 8 % of formed echoes;
  without them 1.03 to 1.25.  The availability (the pinned gate, 36 to 88 %
  of echoes under fading) is unchanged.
* **`joint_fading`, without fading: met at parameter level (0.98 to 1.12);
  1.22 and 1.37 at signal level.**  The excess over the parameter level at
  20 dB is probably the envelope model's error; not tested.

The default stays `"pinned"` for now.  Andrew, 2026-11-27, amending decision
(t): switch to `"joint_fading"` when appropriate, not to `"joint"`, which
under fading states a sigma far below its spread.

Not done: fading at signal level (L2's specular fading through the tier), a
Rician gain or another correlation shape, noise on the orbit scenes, the
several-echo case, and the power-cut experiment (n).

## 12. The fading outliers (decision 2B, 2026-11-27)

Andrew chose to keep working on section 11's shortfall inside the 8 January
2027 timebox.  Parameter level, section 9's scene.  Two seed sets, pooled
below: the 200 per case of section 11 (from 10000) and fresh ones from 30000
(600 per fading case, 200 per case without fading).  Outliers are errors
beyond three fading bounds.  Echoes are those the pinned test accepts; "kept"
is the fraction of them that still give a drift.

### The ambiguity check

`estimate_route_one_joint_fading` reports `nll_margin`: how far the lowest
other maximum of the likelihood lies above the best one in negative log
likelihood, "other" meaning more than three of the best fit's sigmas away in
drift.  It sees the maxima its starts reach.  The tier's `ambiguity_margin`
(off by default) declines an echo whose margin is smaller: no drift, flag
`DECLINED: ambiguous`; `ambiguity_nll_margin`, `drift_alt_mps` and
`drift_pinned_mps` stay.

Each cell: kept; outliers among the kept; RMS error of the kept over the
fading bound.

| margin | 20 dB, 0.5 s (523) | 10 dB, 0.5 s (305) | 20 dB, 2 s (719) | 10 dB, 2 s (305) |
|---|---|---|---|---|
| off | 100 %; 1.1 %; 1.16 | 100 %; 5.6 %; 1.58 | 100 %; 3.3 %; 1.36 | 100 %; 6.6 %; 2.11 |
| 0.5 | 100; 1.1; 1.16 | 97.0; 4.1; 1.48 | 100; 3.3; 1.36 | 99.3; 6.3; 2.10 |
| 1 | 100; 1.1; 1.16 | 95.7; 3.8; 1.47 | 100; 3.3; 1.36 | 98.0; 5.7; 1.60 |
| 2 | 100; 1.1; 1.16 | 91.1; 2.9; 1.33 | 99.4; 3.4; 1.35 | 94.1; 5.6; 1.61 |
| 3 | 99.8; 1.0; 1.14 | 87.5; 2.6; 1.34 | 99.4; 3.4; 1.35 | 92.1; 5.3; 1.60 |
| 5 | 99.2; 1.0; 1.14 | 72.5; 2.7; 1.34 | 98.9; 3.2; 1.35 | 77.4; 4.2; 1.34 |
| 8 | 98.5; 1.0; 1.14 | 46.2; 3.5; 1.47 | 95.7; 2.8; 1.27 | 49.2; 5.3; 1.39 |
| 12 | 94.5; 1.0; 1.13 | 18.7; 7.0; 1.82 | 89.2; 2.2; 1.23 | 21.0; 6.2; 1.56 |

RMS over the bound without the outliers is 1.05 (20 dB, 0.5 s), 1.12, 1.12
and 1.04 with the check off.

Without fading (400, 391 and 154 echoes at 20, 15 and 10 dB) no margin up to
8 declines any echo; 12 declines 4.5 % at 10 dB.

Outliers declined, of all outliers, at margin 2: 0/6, 9/17, 0/24, 4/20.

* **No margin brings all four fading cases within 20 % of the bound.**
* At 10 dB the check works: many misses have a competing maximum (the second
  maximum is within three bounds of the truth in 11 of 17 at 0.5 s).
* At 20 dB it does little: the likelihood has one clear maximum in the wrong
  place (1 of 24 outliers at 2 s has another near the truth), so no test of
  the likelihood's maxima can catch them.
* The 10 dB, 2 s drop from 2.11 to 1.60 at margin 1 is mostly one echo of the
  first seed set, 105 m/s in error.  On the fresh seeds alone that case goes
  from 1.71 to 1.65 at margins 1 and 2.
* The two seed sets cost similar availability at margin 2 (6 to 9 % at 10 dB)
  but differ in what they gain; the counts of outliers are small.

### Does a fade that hides part of the envelope cause the misses?

Tested with the true gain, which only the simulation knows.  Median over
outliers against the rest, and the chance that an outlier ranks as more faded
than another echo (AUC; 0.5 is no information), pooled:

| Case | envelope power received, over its mean: outliers, others | AUC | weaker half's power: outliers, others | AUC |
|---|---|---|---|---|
| 20 dB, 0.5 s | 0.85, 0.90 | 0.57 | 0.56, 0.64 | 0.53 |
| 10 dB, 0.5 s | 0.72, 1.11 | 0.71 | 0.52, 0.77 | 0.67 |
| 20 dB, 2 s | 0.33, 0.82 | 0.80 | 0.22, 0.60 | 0.78 |
| 10 dB, 2 s | 0.47, 1.36 | 0.82 | 0.32, 1.00 | 0.83 |

**Yes for slow fading (2 s), weakly at 0.5 s.**  The shift of the received
envelope's centroid says little (AUC 0.53 to 0.62).

### Cheaper ideas tried

Flags compared at equal cost, dropping the 10 % of echoes each ranks most
suspect (outliers caught; RMS over the bound of the rest), pooled:

| Flag | 20 dB, 0.5 s | 10 dB, 0.5 s | 20 dB, 2 s | 10 dB, 2 s |
|---|---|---|---|---|
| ambiguity margin, smallest | 2/6; 1.13 | 9/17; 1.34 | 8/24; 1.26 | 5/20; 1.61 |
| stated sigma, largest | 2/6; 1.11 | 3/17; 1.48 | 12/24; 1.19 | 7/20; 2.04 |
| power seen in the fit's own lobe, lowest | 0/6; 1.16 | 3/17; 1.51 | 6/24; 1.34 | 3/20; 1.62 |
| the same, weaker half | 0/6; 1.17 | 4/17; 1.51 | 7/24; 1.33 | 4/20; 1.61 |
| quiet samples in the lobe, most | 0/6; 1.17 | 5/17; 1.42 | 6/24; 1.34 | 5/20; 1.59 |

* The fade shows in the truth but not usefully in the data: the three
  measures formed from the fit and the samples rank outliers poorly (AUC 0.26
  to 0.73), probably because the fit shapes its envelope to what it received.
* **The stated sigma is the better flag at 20 dB** (AUC 0.86 and 0.87) and
  poor at 10 dB, 0.5 s (0.59).  Section 11's rule of twice the median sigma
  is too blunt to show it.  A flag combining the margin (10 dB) with sigma
  against the sigma expected at the fitted parameters (20 dB) is not tested.
* **Averaging the maxima** (the drift weighted by each maximum's likelihood,
  no echo declined) gains little: 1.16, 1.47, 1.36, 2.02.  Not built.

### Decided (Andrew, 2026-11-27)

* `ambiguity_margin` stays **off**, which is the default.
* The outlier work stops for now.  M3 reports the outlier fraction beside the
  RMS error.
* The combined margin and sigma flag is left untested.

### Not done

Through the tier with fading; other facet scenes; the combined flag above.

## 13. The power cut of the Doppler history (decision (n))

Route one's Doppler history keeps the sweep pairs whose power is at least
`main_lobe` = 0.1 of the peak.  It feeds the free-slope test, which decides
whether a drift is formed, and the pinned zero crossing.  Section 6 guessed
that at 10 dB the cut sits at the noise and pulls the slope down.  This
compares 0.1 with 0.25 and 0.5.  **An experiment: nothing is adopted and the
default stays 0.1.**

`scripts/wp5_power_cut.py results/wp5_power_cut` at commit 35832cb: section 9's
scene, 200 echoes per case from seed 40000, parameter level.  `joint_fading` is
fitted once per echo behind the same test.  About 3 minutes.

Availability (drift formed), and RMS drift error in m/s of the formed echoes,
pinned / `joint_fading`:

| Case | cut 0.1 | cut 0.25 | cut 0.5 |
|---|---|---|---|
| 15 dB | 98.5 %; 2.52 / 1.47 | 95.5 %; 2.55 / 1.47 | 47.0 %; 4.84 / 1.29 |
| 10 dB | 40.0 %; 7.49 / 2.57 | 53.0 %; 6.86 / 2.43 | 11.0 %; 18.17 / 2.14 |
| 15 dB, tau_c 0.5 s | 62.0 %; 21.8 / 6.46 | 40.5 %; 23.5 / 5.26 | 15.0 %; 25.0 / 4.20 |
| 10 dB, tau_c 0.5 s | 40.5 %; 21.7 / 8.58 | 27.0 %; 23.3 / 10.67 | 8.0 %; 34.8 / 6.70 |
| 15 dB, tau_c 2 s | 74.0 %; 24.6 / 4.63 | 56.0 %; 21.1 / 4.85 | 25.0 %; 22.4 / 4.46 |
| 10 dB, tau_c 2 s | 41.0 %; 24.2 / 12.77 | 31.5 %; 22.5 / 8.49 | 10.5 %; 36.2 / 5.00 |

Noise-free, `c3_verify`:

| cut | untilted formed, worst error | tilted formed | extended declined | plate lengths formed (case D) | slope ratio, untilted |
|---|---|---|---|---|---|
| 0.1 | 56/56, 0.030 m/s | 32/32 | 12/12 | 0.5, 1, 2 km | 0.942-0.996 |
| 0.25 | 56/56, 0.026 | 32/32 | 12/12 | 0.5, 1, 2, **4** km | 0.953-0.997 |
| 0.5 | 56/56, 0.024 | 32/32 | 12/12 | 0.5, 1, 2, **4** km | 0.961-0.997 |

* **0.25 forms more drifts only at 10 dB without fading** (40 to 53 %, pinned
  error 7.49 to 6.86 m/s).  In every other case it forms fewer: 3 points fewer
  at 15 dB, 9 to 22 points fewer under fading.
* **0.5 forms far fewer everywhere**, and the pinned drift is worse.
* **Either higher cut lets the 4 km plate pass as compact**, which 0.1
  declines: the compact-facet test gets weaker.
* The `joint_fading` error changes mostly through which echoes pass the gate:
  on the echoes formed at every cut it is the same by construction (the fit
  does not use the cut).  At 10 dB those common sets are 4 to 7 echoes, too
  few to compare.
* The cut stays at 0.1.

## 14. The diffuse halo (decision (e))

**Model** (`l1_propagation/facet.py`; option A, Andrew, 2026-11-27).  A plate
may carry `halo: {k, ratio, n_points, seed}`: 200 weak point scatterers, equal
in all directions, uniform along the track over `ratio` facet lengths, random
phases, moving with the medium.  `K` is the facet's power over the halo's mean
power **at overflight**; the halo's power then stays about constant through the
pass while the facet's rises and falls.  `K = 0` leaves the facet out and keeps
the halo at the power the facet would have had.  Checked in
`tests/unit/test_l1_facet_halo.py`: mean power `1/K` through the pass; the
field correlation `1/e` at the uniform patch's time (`fading.halo_decorrelation_s`);
the drift's Doppler; the facet's specular fading leaves it alone.

**Measured** with `scripts/wp5_halo.py results/wp5_halo` at commit b1a37f2: the
straight scene of section 6 (2 km facet, 5 MHz, drift 25 m/s, 209 sweeps)
through L1, L2 and the tier, 48 halo draws per case, white noise at 20 and
10 dB per sweep of the facet's echo, or none.  About 4 minutes.

Drifts formed, of 48 (no noise / 20 dB / 10 dB):

| Halo | formed |
|---|---|
| none | 48 / 48 / 26 |
| 2 facet lengths, K 10 | 47 / 47 / 39 |
| 2, K 3 | 37 / 38 / 41 |
| 2, K 1 | 29 / 34 / 33 |
| 1, K 3 | 25 / 29 / 32 |
| 4, K 3 | 41 / 41 / 41 |
| 2, **K 0 (no facet)** | **24 / 27 / 25** |

RMS drift error of the formed, m/s, pinned / joint / joint_fading:

| Halo | no noise | 20 dB | 10 dB |
|---|---|---|---|
| none | 0.07 / 0.11 / 0.08 | 1.10 / 0.71 / 1.08 | 5.65 / 2.81 / 3.04 |
| 2, K 10 | 17.0 / 14.8 / 53.5 | 17.1 / 14.8 / 49.4 | 21.1 / 23.1 / 46.1 |
| 2, K 3 | 32.7 / 36.4 / 66.9 | 32.6 / 36.2 / 65.7 | 33.7 / 36.0 / 69.4 |
| 2, K 1 | 60.4 / 64.1 / 80.0 | 62.8 / 65.3 / 73.2 | 65.3 / 69.1 / 81.6 |
| 1, K 3 | 44.2 / 47.6 / 67.2 | 48.7 / 49.3 / 65.0 | 53.0 / 63.8 / 80.6 |
| 4, K 3 | 62.7 / 53.9 / 47.7 | 62.5 / 53.1 / 42.7 | 56.0 / 56.3 / 56.9 |
| 2, K 0 | 100 / 103 / 81 | 100 / 116 / 77 | 105 / 133 / 90 |

An echo at the facet's range was found in every draw.  Robust spreads (1.4826
MAD, in `summary.txt`) are of the same order, except `pinned` and `joint` at
4 lengths, where they are about half the RMS (outliers).

* **The halo, not the noise, sets the drift error.**  At K = 10 (the halo 10 dB
  under the facet at overflight) the noise-free error goes from 0.1 to
  15-17 m/s with `pinned` and `joint`; at K = 3 it is 33 to 63 m/s.
* **`joint_fading` is the worst of the three under a halo** with the facet
  present, except at 4 lengths: its model has fading on the facet but no halo, and a halo of this
  kind is added power, not a gain on the facet.
* **`K = 0` should be declined and is not**: a drift is formed on half the
  draws, wrong by about 100 m/s.
* Option A is the harsher form of the halo: its power is received through the
  whole 4.6 s window, the facet's mostly within its 0.9 s lobe.  A halo that
  rose and fell with the pass (option B) was not run.

### Why a random patch passes the compact-facet test

Each point of a randomly phased patch sweeps in Doppler at the ephemeris slope,
with its own zero crossing.  Their sum keeps that slope, so route one's
free-slope test reads about 1; only the contrast test (a peak above the window's
edges) can decline it, and a random envelope has such a peak often.  Checked
without L1, L2 or the tier, from route one's own Doppler law (400 draws, no
noise):

| Patch length | 0.5 km | 2 km | 4 km | 8 km | 32 km |
|---|---|---|---|---|---|
| drift formed | 0.8 % | 18.5 % | 30.5 % | 41 % | 54.5 % |
| slope ratio within 0.8 to 1.2 | 100 % | 100 % | 100 % | 100 % | 93.8 % |
| RMS drift error of the formed, m/s | 10 | 105 | 123 | 133 | 140 |

Through the chain at K = 0, 12 draws formed 5, with slope ratios 0.93 to 1.02
and a median contrast of 2.7 against the threshold of 4.  A facet plus a patch
built the same way gives errors of the same size as the chain's (pinned, 300
draws: 17.5 m/s at 2 lengths and K 10, 50 at K 3, 84 at K 1, 49 at 1 length and
K 3, 37 at 4 lengths and K 3).

This **contradicts an expectation** carried since Phase 1
(`docs/decorrelation.md`, "that route one declines a randomly phased patch is
expected, not yet run").  The extended reflector route one does decline is a
specular one (the tapered 80 km plate), whose Doppler does not sweep.  The
rough extended reflector is not declined.  The Technical Volume's statement
that a wholly random echo is an extended reflector rather than a compact facet
holds as a description of the scene, but route one's tests do not detect it.
Telling the two apart is open; nothing is changed in the estimator.

### Not done

Orbit scenes with the halo (the code runs; only checked for amplitude); a
cross-track extent; scatterers with a narrow scattering angle; the halo
together with the facet's specular fading.

## 15. Poole's velocity as the comparator (TV Section 8.21, item 10)

Poole (1985) forms a group range from the phase step across a change of start
frequency and a velocity from the phase step between sweeps.  A drift sequence
has one start frequency, so **only the velocity half is formed; the range half
cannot be formed from this data** (option A, Andrew, 2026-11-27).
`poole_range.velocity_series` takes the principal value between successive
sweeps (+/-675 m/s of phase-range rate at 5 MHz and 45 Hz).  `-V* / mu_sat` is
the line-of-sight rate in route one's sense: platform, drift and tilt together.

`scripts/wp5_poole.py results/wp5_poole` at commit f06e2fa, on the sequences of
sections 6 and 14, 48 noise (or halo) draws per case.  About 4 minutes.

**Pair by pair**, against the line-of-sight rate L1 put into the facet's echo,
over the main lobe, median RMS in m/s: Poole / the tier's phase steps, then
Poole against the tier.

| Scene | no noise | 20 dB | 10 dB |
|---|---|---|---|
| overhead | 0.51 / 0.22; 0.31 | 41 / 35; 20 | 231 / 223; 159 |
| 2 km along | 0.51 / 0.22; 0.31 | 42 / 36; 20 | 233 / 217; 154 |
| halo, 2 lengths, K 3 | 28.0 / 28.3; 0.72 | 47 / 45; 19 | 135 / 119; 62 |
| 80 km plate | 1.20 / 0.63; 0.58 | not detected | not detected |

**As a drift** (true 25 m/s), RMS error in m/s.  "Peak" is Poole's rate at the
strongest pair; "mean" is its power-weighted mean over the lobe.  The joint
columns and `pinned` form a drift only where the compact test passes (count in
the last column); the Poole columns here are over all 48.

| Scene, level | Poole peak | Poole mean | pinned | joint | joint_fading | route one formed |
|---|---|---|---|---|---|---|
| overhead, none | 1.31 | 0.36 | 0.07 | 0.11 | 0.08 | 48 |
| overhead, 20 dB | 23.6 | 1.26 | 1.10 | 0.71 | 1.08 | 48 |
| overhead, 10 dB | 77.6 | 9.44 | 5.65 | 2.81 | 3.04 | 26 |
| 2 km along, none | 0.91 | 0.28 | 0.07 | 0.21 | 0.12 | 48 |
| 2 km along, 20 dB | 26.5 | 1.72 | 1.25 | 0.88 | 1.07 | 48 |
| 2 km along, 10 dB | 80.7 | 10.1 | 8.32 | 2.19 | 2.28 | 22 |
| 80 km plate, none | 0.03 | 0.09 | declined | declined | declined | 0 |
| halo, none | 70.5 | 29.6 | 32.7 | 36.4 | 66.9 | 37 |
| halo, 20 dB | 76.9 | 29.2 | 32.6 | 36.2 | 65.7 | 38 |
| halo, 10 dB | 78.8 | 25.3 | 33.7 | 36.0 | 69.4 | 41 |

On the echoes route one formed at 10 dB, Poole's mean is 10.4 m/s overhead and
9.6 m/s 2 km along.

* **Pair by pair Poole's rate is the tier's own observable**: without noise the
  two agree to 0.3 m/s.  With noise Poole's is 4 to 17 % further from the truth,
  probably from its window over the whole record against the tier's gate (not
  tested).
* **One pair is no drift estimate under noise** (24 to 81 m/s).
* **Averaged over the lobe it comes close to the drift**, because the platform's
  part cancels when the whole lobe is in the window.  That average is the first
  term of route one's pinned drift (section 5).  It is behind `joint` at 20 dB
  (1.3 to 1.7 against 0.7 to 0.9 m/s) and well behind at 10 dB.
* **What Poole's method does not give**: whether the echo is a compact facet (it
  returns a number for the 80 km plate; right here only because that plate is
  flat and overhead, so the line-of-sight rate is the drift; a tilted layer was
  not run), the overflight time that corrects a lobe not centred in the window
  (not run: every window here held the whole lobe), and an uncertainty.
* **Under the halo none of them works**; Poole's mean is no worse than route
  one's.
* **Item 10, velocity half**: in these cases Poole's observable adds no
  information beyond the phase steps route one already uses.  The range half
  stays open: it needs a second start frequency in the drift mode.
* The 80 km plate was not detected with noise, which is set from the compact
  facet's echo; its noisy cases are not measured.

