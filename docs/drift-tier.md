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
reproduces the sheet model's phase to **1.2e-10 rad** (0.01 asked) and its
Doppler to **1e-10** (1e-3 asked); `field_aligned_verify` unchanged.  The exit
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
now; their measured errors are unchanged.

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
| Orbit, 3-D plate 300 km below 600 km, L0-style exponential profile | 50 m/s | 49.996 |
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
without `v_eff^2 = v^2 P'/P'_flat`.  Measured: it does not.  Route one's
drift is the Doppler at the power peak converted to a velocity (the zero
crossing is taken with the ephemeris slope, so `v^2` cancels to first order);
with and without `v_eff` the orbit drift differs by 0.0004 m/s.  What the
curvature moves is the **slope ratio**, the compact-facet test, by
`(v_eff/v)^2`: 0.965 with `v_eff`, 0.922 without (0.956).  So decision (d)
(the facet's true height from the onboard profile) matters to the
classification, not to the drift.
