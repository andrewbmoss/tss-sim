# The drift mode and what one node delivers

WP5, claim C3.  This is the output contract for the drift-mode L3 tier
(`l3/drift`, planned for WP5 Phase 3) and for the ionogram-mode INS-09
correction (Phase 5).  It says, for each class of echo, what a single
spacecraft with **one receive channel** delivers, how well, under which
conditions, and which flags it sets when it cannot deliver.  The code that
implements it does not exist yet; `DriftProduct` in `interfaces.py` is written
from this document in Phase 2.

Status: **complete** (2026-11-24).  Part 1 (morning): the products common to
every echo, and compact facets.  Part 2 (afternoon): field-aligned sheets, the
extended layer and the ionogram mode, with the results of the two screens run
for it (the wideband check D8 and the null scan C4, section 10).

Every accuracy here is **parameter level, to be replaced at signal level in
Phase 3**.  Each figure names its source: a study in this repository, a TV
table, or the screening script
`Dissertation\Plans\wp5_screening\single_node_gaps.py` (written for this
contract, output `single_node_gaps.csv`, parts A to G).

## Decisions this contract rests on

- C3 is the three-route separation of plasma drift from platform motion
  (Proposal Rev E, TV Ch. 8).  Up and down sweep pairing is not the
  contribution; it is kept as the ionogram-mode INS-09 bias check.
- Route one is the best solution available with one receive channel (T1,
  2026-11-14).  Arrival-direction measurement, which needs a second channel,
  belongs to WP11.
- M3 is route one alone at signal level (2026-11-24).
- The several-frequency lever is a standard facet product, set per facet when
  the facet's reflection band allows it (2026-11-24).
- **SCI-15 is scored against the plasma's own rate** (2026-11-24, evening):
  the plasma's velocity along the echo's true wave normal, in the Earth-fixed
  frame.  The platform's projection on a reflector tilt, which survives the
  removal of the ephemeris term for a local-down wave normal, counts as error.
  The combined rate `los_rate_mps` is reported, but it scores SCI-15 only when
  a class product removes or bounds the tilt term (section 4, "What SCI-15 is
  scored against").
- **The unambiguous interval is limited to +/-1 km/s, and claimed only where
  sampling supports it** (2026-11-24, evening): `f mu_sat` <= 5.3 MHz at 71 Hz,
  3.4 MHz at 45 Hz.  The +/-1.5 km/s goal is not claimed.  Resolving the wrap
  with drift-mode sub-bands at different frequencies is recorded as a possible
  future improvement of the plasma-velocity mode, not adopted.
- **`tau_c_s` and `facet_size_m` are reportable data products** (2026-11-24,
  evening).
- Reference case, as in TV Ch. 8: plasma frequency at the spacecraft
  `f_N,sat` = 1.5 MHz, so `mu_sat` = 0.87, 0.95, 0.99 at 3, 5, 10 MHz;
  `v_s` = 7.6 km/s; 10 dB per sweep; 45 Hz sweep repetition, 71 Hz cap
  (INS-27); ranges are apparent ranges `P'`.

## 1. Inputs

| Input | Source | Used for |
|---|---|---|
| De-chirped record per sweep, drift-mode sub-band | L2 via `Pipeline.run_sequence` (Phase 2) | Everything |
| Sweep start times | sweep law (`sweep_law.drift_mode`, Phase 2) | Phase steps, envelope timing |
| Platform position and velocity per sweep, **Earth-fixed (corotating) frame** | Ephemeris (`PlatformState`) | Removing the platform Doppler; route one's `v_s` and timing |
| Plasma frequency at the spacecraft `f_N,sat` | WP4 tracker (resonance at the trace floor) | `mu_sat` for every velocity |
| Noise power per range cell | L3 detector noise reference (WP4 P2) | SNR, decorrelation time |
| Ionogram-mode trace near the drift-mode frequency | WP4 tracker | Range prior, O/X label, class prior |
| Magnetic field direction | IGRF | Field-aligned sheets (part 2) |

Attitude is not needed for the Doppler of one channel, because the platform
term depends on the wave normal and not on the antenna axis.  It enters only
through the antenna gain.

The velocity must be the spacecraft's velocity relative to the corotating
Earth, because the plasma corotates.  An inertial velocity leaves up to about
0.46 km/s of Earth rotation in the platform term.

## 2. Processing common to every echo

1. For each sweep `n` and drift-mode frequency `f`, predict the platform
   Doppler for a wave normal along local down:
   `f_D,plat = (2 f / c) mu_sat (v_ECEF . k0)`, with `k0` the unit vector from
   the spacecraft toward the sub-satellite point.  This is the
   ephemeris-predicted platform term of SCI-15 and INS-09.  Remove it from the
   record by a phase ramp before any Doppler processing, so that the residual
   stays inside the unambiguous interval (TV Table 7.2, WP5 row).
2. Range-gate each detected echo with the existing transform
   (`tss.l3_signal.transform`) and take one complex value per sweep.
3. Doppler per sweep pair from the phase step (`doppler_history_from_field`),
   main lobe only (steps across a pattern null are dropped).
4. Power per sweep, and its envelope.
5. Classify the echo (section 3).
6. Form the products of section 4 for every echo, then the class products
   (section 5 onward).

## 3. Echo classification

| Label | Rule | Source of the rule |
|---|---|---|
| `compact_facet` | Free-fit Doppler slope within +/-20 % of the ephemeris slope `(2f/c) v_s^2 / P'`, **and** power contrast (peak over window edges) >= 4 | `route_one_estimator.estimate_route_one`; `c3_verify` (switch at the Fresnel radius, 3.9 km at 5 MHz and 500 km) |
| `extended` | Slope ratio outside +/-20 % (a layer's Doppler history does not sweep; ratio <= 6e-4 for an 80 km plate) | `c3_verify` case C |
| `field_aligned_candidate` | Tested before `extended`.  The mean residual Doppler after the platform term is predicted for a field-aligned normal (section 6, step 2) lies inside the drift prior (+/-250 m/s), and the residual Doppler history is flat (slope ratio inside +/-20 % of zero against the ephemeris facet slope), while the local-down prediction leaves an offset larger than the drift prior.  On a near-meridional track at dip above a few degrees the two predictions differ by `v_s sin I mu_f/mu_sat`, which is 1.3 km/s times `mu_f/mu_sat` at a dip of 10 degrees, so the test is decisive there.  Near the dip equator the two coincide and the echo is classed by the other rules | Section 6; `field_aligned_verify` geometry.  Thresholds provisional, to be set in Phase 3 |
| `unresolved` | Power envelope narrower than about 4 sweeps, or SNR below the detector threshold, or contrast < 4 with a sweeping slope | `route_one_tilt_map` (unresolved envelope); detector P_fa 1e-4 |

The label is an output, with the slope ratio and the contrast that decided it.
The rules run in the order `compact_facet`, `field_aligned_candidate`,
`extended`, `unresolved`.  The order matters because a field-aligned sheet on a
near-meridional track is long along track, so it fails route one's slope test
and would otherwise be classed as an extended layer.

## 4. Products for every echo

Accuracy terms for the combined rate (part A of the screen): the Doppler term is
the tone-frequency bound over the dwell, `sqrt(6) / (2 pi T sqrt(rho f_rep T))`,
converted with `c / (2 f mu_sat)`; the oscillator term is `c sigma_y / (2 mu_sat)`
(INS-26); the ephemeris term is 0.1 m/s of velocity knowledge
(`route_one.ephemeris_term_mps`).  These assume the echo stays coherent over the
dwell (`tau_c` >= `T`).

| Product | Units | Definition | Accuracy (parameter level) and conditions | Source |
|---|---|---|---|---|
| `t_s` | s | Sweep times of the sequence | Exact | Sweep law |
| `range_m` | m | Apparent range `P'` of the echo per sweep (range gate) | Resolution `c / (2 B_sub)`: 1.0 km at the recommended `B_sub` of 150 kHz (section 10).  Fixed in Phase 2 | Section 10; Phase 2 |
| `power_db`, `snr_db` | dB | Echo power and SNR per sweep and at the peak | SNR against the detector's noise reference | WP4 P2 |
| `platform_doppler_hz` | Hz | `f_D,plat` removed in step 1 | Ephemeris velocity to 0.1 m/s gives 0.1 m/s | `route_one.ephemeris_term_mps` |
| `doppler_hz` | Hz | Residual Doppler per sweep pair | As `los_rate_mps` | Step 3 |
| `los_rate_mps` | m/s | **Combined line-of-sight rate**, `c (f_D - f_D,plat) / (2 f mu_sat)`, averaged over the dwell.  For a reflector tilted by `alpha` along track it is `w cos(alpha) + v_s (mu_f / mu_sat) sin(alpha)`: drift plus the platform's projection on the tilt, which one channel cannot remove | Dwell 1.15 s (INS-15 goal): **8.7-10.0 m/s** with the oscillator at the INS-26 threshold, **0.5-1.0 m/s** at the INS-26 goal, 3-20 MHz.  Dwell 0.3 s: 8.7-11.9 m/s (threshold), 1.0-6.5 m/s (goal).  Dwell 0.05 s (INS-15 threshold): 15-95 m/s, Doppler term dominant, worst at 3 MHz | Screen part A; the oscillator term is `route_one`'s; matches the Doppler and oscillator terms of TV Table `c3r1budget` |
| `los_rate_sigma_mps` | m/s | 1-sigma of `los_rate_mps` from the terms above | | |
| `resolution_mps` | m/s | Separation at which two scatterers in one range cell are resolved in Doppler, `c / (2 f mu_sat T)` | Dwell 1.15 s: 50, 27, 13 m/s at 3, 5, 10 MHz.  The SCI-15 goal of 50 m/s between scatterers is met at 1.15 s from 3 MHz up, which is how INS-15's goal was derived | Arithmetic; INS-15 |
| `v_unambiguous_mps` | m/s | `c f_rep / (4 f mu_sat)` | 45 Hz: 1298, 707, 341, 169 m/s at 3, 5, 10, 20 MHz.  71 Hz: 2048, 1116, 538, 267 m/s.  **The +/-1 km/s interval of INS-09 and SCI-15 is claimed only where sampling gives it: `f mu_sat` <= 5.3 MHz at 71 Hz, 3.4 MHz at 45 Hz** (decision of 2026-11-24).  Above that the product reports the interval it has and sets `ALIAS_RISK`.  The +/-1.5 km/s goal is not claimed | Screen part B |
| `range_debiased_m` | m | `range_m` corrected for the measured rate in the drift-mode sweep, bias `mu_sat Rdot f / K_sub` | Follows `los_rate_sigma_mps` | `doppler.apparent_range_bias_m` |
| `doppler_spread_hz` | Hz | Spectral width of the residual after the mean Doppler (and, for a facet, the ephemeris slope) is removed | To be established in Phase 3.  Used to flag several scatterers in one cell or a rough reflector | Phase 3 |
| `tau_c_s` | s | **Data product.**  **The echo's own decorrelation time**, from `|R(1)| / (R(0) - noise)` of the complex echo (the modulus removes the constant Doppler, so no Doppler estimate is needed) | 10 dB per sweep, 45 Hz, at the decorrelation times carried since G1 (2026-11-24; `Plans\wp5_screening\decorrelation_update.py`), Rayleigh fading: tau_c 0.1 s within a factor 0.69-1.38 (16-84 %) over a 1.15 s dwell and 0.79-1.27 over 3 s; 0.5 s within 0.37-2.4 and 0.55-1.66; 1.0 s a lower bound over 1.15 s and 0.45-2.8 over 3 s.  A phase-wandering specular echo (the medium's own decorrelation in `l2_instrument.fading`) is measured a little better (0.82-1.38, 0.66-2.3 over 1.15 s; 0.88-1.19, 0.78-1.53, 0.65-2.1 over 3 s).  Rule: an estimate when `tau_c` <= about half the dwell, otherwise a lower bound.  The time measured is the echo's, from the medium and the platform's motion together (`docs/decorrelation.md`); the estimator does not separate them | Screen part G |
| `echo_class` | label | Section 3 | | |

`tau_c_s` matters beyond the echo.  The topside decorrelation time gates route
one (WP5 decision D3) and nothing measures it.  A sounder that reports it for
every echo measures it in flight, at the frequencies and heights it sounds.  It
is a reportable data product (decision of 2026-11-24), with `TAU_LOWER_BOUND`
set where it is only a bound.

### What SCI-15 is scored against

The L5 truth for SCI-15 is the plasma's own velocity along the echo's true wave
normal at the spacecraft, in the Earth-fixed frame (decision of 2026-11-24).
The derivation is TV section 8.10.4 (`sec:c3sci15`, figure `fig:c3sci15`).  In
short, with `x` along track and `z` up, spacecraft velocity
`V_s = (v_s, 0, v_z)`, plasma velocity `U = (u, 0, w)` and wave normal at the
spacecraft `k = (sin th, 0, -cos th)`, `th` = `theta_sat`:

    f_D    = (2 f mu_sat / c) k . (V_s - U)                 (eikonal)
    V_meas = c (f_D - f_D,plat) / (2 f mu_sat)
           = -k . U  +  (k - k0) . V_s                       (k0 = local down)
    V_true = -k . U = w cos(th) - u sin(th)                  (SCI-15 truth)
    eps    = v_s sin(th) + v_z (1 - cos th)  ~  v_s sin(th)
           = v_s (mu_f / mu_sat) sin(alpha)                   (Snell: mu sin(th) conserved)

So `los_rate_mps` = `V_true + eps`, and against this truth `eps` is an error
of 132.6 m/s per degree of `theta_sat`.  The goal (50 m/s) is used up at
0.38 degrees and the threshold (100 m/s) at 0.75 degrees.  A product that uses
a predicted or estimated wave normal `k_p` in place of `k0` leaves
`v_s (sin th - sin th_p)`, so SCI-15 is met to the extent the along-track wave
normal is known.  A bias `b` left in a product counts with its random part as
`sqrt(sigma^2 + b^2)`.  The consequences:

| Product | Scores SCI-15? |
|---|---|
| `los_rate_mps` of any echo | Only for a reflector whose tilt is known to be negligible; otherwise the tilt term is in the error.  Its 8.7-10.0 m/s at 1.15 s is a floor, not the SCI-15 result |
| `drift_mps` (route one, compact facet) | Yes, with the tilt bias of `tilt_allowance_deg` counted as error (`TILT_UNFLAGGED`) |
| `drift_mf_mps` (several-frequency lever) | Yes; the tilt term is estimated and removed |
| `drift_fa_mps` (field-aligned sheet) | Yes within `NEAR_MERIDIAN` / `OFF_MERIDIAN`, with `phi_error_mps` counted; refused beyond |
| Extended layer | No product; SCI-15 is not met by one node with one channel (section 7) |

The requirement text (`requirements.yaml`, SCI-15 notes) says "the
line-of-sight velocity of an isolated echo, in the Earth-fixed frame after the
ephemeris-predicted platform term is removed".  The decision reads "velocity"
as the plasma's, and the notes now say so.

`mu_sat` sensitivity (part C): 1 % of error in `f_N,sat` scales every rate by
1.3 % at 2 MHz, 0.33 % at 3 MHz, 0.10 % at 5 MHz and 0.02 % at 10 MHz: 0.3 m/s
of 100 m/s at 3 MHz, 1.3 m/s at 2 MHz.  It matters only below about
2 `f_N,sat`, which is where INS-09 stops applying.

### Flags for every echo

| Flag | Set when | Consequence |
|---|---|---|
| `EVANESCENT` | `f` <= `f_N,sat` | No echo; products NaN |
| `NEAR_PLASMA_FREQUENCY` | `f` < 2 `f_N,sat` | Rate carries the `mu_sat` sensitivity above; outside INS-09's range |
| `SHORT_DWELL` | Dwell under 0.05 s (INS-15 threshold) | Rate reported; accuracy from the Doppler term |
| `DECORRELATED` | `tau_c_s` < dwell | The coherent Doppler bound does not hold; the rate is incoherently averaged and its sigma is widened |
| `TAU_LOWER_BOUND` | `tau_c_s` > half the dwell | `tau_c_s` is a lower bound |
| `ALIAS_RISK` | Predicted residual rate (drift prior plus the tilt allowance) beyond `v_unambiguous_mps`, or `f mu_sat` above the limit of the +/-1 km/s claim | Rate may be wrapped; outside the +/-1 km/s claim |
| `SPREAD` | `doppler_spread_hz` above the single-scatterer value | Several scatterers or a rough reflector; class products are not formed |
| `TILT_DEGENERATE` | Always, unless a class product removes it | `los_rate_mps` is drift plus `v_s (mu_f/mu_sat) sin(alpha)`; it is not the plasma's own rate and does not score SCI-15 |

## 5. Compact facets: route one

Inputs are those of section 1 plus the facet's range gate.  Processing is
`estimate_route_one`: the overflight time `t_pk` from the four-parameter
`sinc^2` fit of the power envelope; the zero crossing `t_zc` with the slope
imported from the ephemeris; drift `w = v_s^2 (t_zc - t_pk) / (mu_sat P')`.

| Product | Units | Definition | Accuracy (parameter level) and conditions | Source |
|---|---|---|---|---|
| `drift_mps` | m/s | Route-one drift of the facet | 2 km facet, 10 dB per sweep, 45 Hz, 3-10 MHz, P' 500-1000 km: **10.4-12.2 m/s** without fading, **14.5-33.7 m/s** at tau_c 0.5 s (oscillator at the INS-26 threshold; 5.1-7.0 m/s at the goal, no fading).  A 0.5 km facet at tau_c 0.5 s: 59 m/s.  Noise-free: within 0.027 m/s (56/56).  **Biased by facet tilt and not flagged** unless the several-frequency lever runs | `route_one` study, TV Table `c3r1budget`; `c3_verify` |
| `drift_sigma_mps` | m/s | Budget RSS | | `route_one_budget` |
| `t_pk_s`, `t_zc_s` | s | Overflight time and Doppler zero crossing | `sigma_t0` 0.023-0.054 s without fading, 0.091-0.434 s at tau_c 0.5 s | `route_one` |
| `slope_ratio`, `peak_contrast` | - | The classification statistics | | Section 3 |
| `along_track_m` | m | Along-track position where the spacecraft crosses the **facet normal**, `v_s t_pk` | **175-410 m** without fading, **0.7-3.3 km** at tau_c 0.5 s.  It is the facet's position only for an untilted facet; a tilt of `alpha` moves it by `h tan(alpha)` | Screen part D from the `route_one` budget; `c3_verify` finding |
| `facet_range_m` | m | Apparent range at the power peak | As `range_m`; true height needs the WP6 inversion | Phase 2; WP6 |
| `facet_size_m` | m | **Data product.**  Along-track size of the specular region, `L = 0.443 lambda P' / (W v_s)` from the fitted envelope width `W` | **3-9 %** (1-sigma) without fading for W 0.44-2.9 s; **14-51 %** with fading at tau_c 0.1-1.0 s, biased low by 3-30 % (fading narrows the fitted envelope).  It measures the coherent specular region, which a curved or rough facet makes smaller than the facet | Screen part E; `route_one.envelope_width_s` |
| `tilt_allowance_deg` | deg | Facet tilt at which the unflagged bias `v_s (mu_f/mu_sat) sin(alpha)` plus the random budget reaches 50 or 100 m/s | Goal / threshold: about 0.3 / 0.7 deg near the spacecraft (`mu_f` = `mu_sat`), 0.6 / 1.3 deg at `mu_f` 0.5, 1.5 / 3.3 deg near reflection (`mu_f` 0.2); 0.06-0.35 deg at 3 MHz with fading | `route_one_tilt_map`, TV Table `c3tiltallow` |
| `arrival_angle_deg` | deg | Along-track arrival angle, by-product for a facet of known orientation | 0.02-0.45 deg | `route_one` study note (2026-10-06) |
| `tau_c_s`, `doppler_spread_hz` | | As for every echo, with the envelope divided out | | Section 4 |

### Several-frequency lever (standard product when the band allows)

Route one at frequencies inside one facet's reflection band returns
`v(f) = w + B g(f)`, with `B = v_s sin(alpha)` and `g = mu_f / mu_sat`, which is
linear in the drift and the tilt term (TV §8.9.3).  Three frequencies spread
from 1.02 `f_Nf` to 0.995 `f_Nf sqrt(1 + delta)`, where `delta` is the facet's
density enhancement (Dyson 1969: 5-70 %), separate them.

| Product | Units | Accuracy (parameter level), three frequencies, `f_Nf` 2.5-4 MHz | Source |
|---|---|---|---|
| `drift_mf_mps` | m/s | delta >= 20 %: **18-28 m/s** without fading, **46-69 m/s** at tau_c 0.5 s.  delta 10 %: 47-48 and 118-121 m/s | Screen part F; TV Table `c3multifreq` |
| `tilt_term_mps` | m/s | delta >= 20 %: 35-80 m/s without fading, 86-199 m/s at tau_c 0.5 s.  delta 10 %: 160-179 and 400-448 m/s | Screen part F |
| `facet_tilt_deg` | deg | delta >= 20 %: **0.26-0.60 deg** without fading, 0.65-1.5 deg at tau_c 0.5 s.  delta 10 %: 1.2-1.4 and 3.0-3.4 deg | Screen part F |
| `along_track_corrected_m` | m | `along_track_m` minus `h tan(alpha)` with the measured tilt, so the facet's position rather than its normal's crossing | Follows `facet_tilt_deg` | Derived |

The lever is the only single-channel route to a facet drift without the tilt
bias, and the only single-channel measurement of a facet's along-track tilt.
Its conditions:

- the background is horizontally stratified (the assumption `c3_verify` made);
- the facet reflects over the whole band used and persists while the drift
  mode visits the three sub-bands;
- the drift mode schedules about three sub-bands per facet, which shares the
  sweep repetition between them and raises the per-frequency error;
- errors are taken as equal and independent between frequencies;
- it does nothing for an extended layer.

### Facet flags and refusals

| Flag or refusal | Set when | Source |
|---|---|---|
| Declined: `extended` | Slope ratio outside +/-20 % | `estimate_route_one` |
| Declined: no peak in window | Contrast < 4 | `estimate_route_one` |
| `UNRESOLVED_ENVELOPE` | Envelope narrower than about 4 sweeps (shallow facets at high frequency) | `route_one_tilt_map` |
| `TILT_UNFLAGGED` | Set whenever the lever did not run: `drift_mps` may carry a tilt bias up to the `tilt_allowance_deg` scale | `c3_verify` finding |
| `MF_APPLIED` / `MF_NOT_AVAILABLE` | Lever run, or not (band narrower than three usable frequencies, facet lost between sub-bands, or sweep repetition insufficient) | Section above |
| `FACET_ALIAS` | Tilt beyond the facet Nyquist limit at 45 Hz: about 0.8 deg at 20 MHz, 1.6 deg at 10 MHz, 3.2 deg at 5 MHz (the Doppler history's extremes alias before a steady echo would) | `route_one_tilt_map` |
| `DECORRELATED` | `tau_c_s` shorter than the envelope width `W` | Section 4 |

### What one node cannot give for a facet

- The cross-track position and the cross-track tilt: one channel has no
  cross-track observable, and a facet left or right of the track at the same
  range is indistinguishable.  That is WP11's.
- Several facets in one range cell closer than `resolution_mps` in Doppler.
- A facet larger than about the Fresnel radius behaves as an extended
  reflector and is declined.

## 6. Field-aligned sheets

A sheet that contains the magnetic field has its normal perpendicular to B.  In
coordinates along track, across track and up, with the track at heading `psi`
from magnetic north and dip `I`, the along-track tilt of the normal is

    tan(alpha_along) = (tan I - (n_c / n_v) sin(psi)) / cos(psi)

(TV eq. `c3fieldtilt`).  On the magnetic meridian (`psi` = 0) the tilt is the
dip, whatever the sheet's rotation `phi` about B, so the field model supplies
the orientation that one channel cannot measure.  Off the meridian the unknown
rotation enters with weight `sin(psi)`.  A sheet on a near-meridional track is
long along track, so the estimator is the plane form with the tilt known, not
route one's peak timing.

### Inputs

Those of section 1, plus: the IGRF field direction (dip `I`, declination) along
the track; the track heading, from which `psi` = heading minus declination; and
the field line the echo comes from, which IGRF and the measured range identify.

### Processing

1. Steps 1 to 4 of section 2.
2. Predict the platform Doppler for the specular normal of a field-aligned
   sheet on the echo's field line, with `phi` = 0, curvature of the field line
   included (`field_aligned_estimator`, the "assumed zero" treatment).  Remove
   it by a phase ramp.  The predicted term is large (`v_s sin I mu_f/mu_sat`, about
   1.3 km/s times `mu_f/mu_sat` at a dip of 10 degrees), so it has to come off
   before the drift mode samples the Doppler (TV §8.17).
3. Classify (section 3).  For a `field_aligned_candidate`, read the drift as
   the mean of the residual, converted with `mu_sat`.
4. Form the error term for `phi` from the heading and the dip (below) and set
   the flags.

### Products

| Product | Units | Definition | Accuracy (parameter level) and conditions | Source |
|---|---|---|---|---|
| `drift_fa_mps` | m/s | Drift along the ray direction at the spacecraft, `s . V`, with the sheet's orientation from IGRF and `phi` = 0 | Noise-free model error, `phi` given (an ideal cross-track arrival angle, which one channel does not have): **0.38 m/s** in free space, 0.08 m/s in the plasma, every case.  With `phi` assumed zero, which is what one channel delivers: on the meridian a second-order error `v_s sin I (1 - cos phi)`, 0 at dip 0, **5.0 and 20.0 m/s** at dip 10 degrees and `phi` 5 and 10 degrees, 14.5 and 57.7 m/s at dip 30 degrees; off the meridian a **first-order error `v_s sin(psi) sin(phi)`**, -167 m/s at `psi` 15 degrees and `phi` 5 degrees, -465 m/s at 45 degrees.  Inside the plasma both are smaller by about `mu_f/mu_sat` (-87 to -90 m/s at `psi` 15 degrees; 2.6 m/s on the meridian at dip 10 degrees).  The random part is that of `los_rate_mps` (8.7-10.0 m/s at a 1.15 s dwell) | `field_aligned_verify` (seed 20261009, 90 free-space and 8 plasma trials at 5 MHz); TV Table `c3sheet` |
| `phi_error_mps` | m/s | Model-error term for an assumed spread of `phi`, `v_s (mu_f/mu_sat) [sin(psi) sin(phi_max) + sin I (1 - cos phi_max)]` | Arithmetic from the closed form above.  `phi_max` is a configuration value (5 degrees until the spread of real rotations about B is known) | Closed form, TV Table `c3sheet` |
| `heading_from_meridian_deg` | deg | `psi` | From the ephemeris and IGRF | IGRF |
| `alpha_along_deg` | deg | Along-track tilt of the sheet's normal predicted from IGRF | Set by the IGRF dip error and the departure of real irregularities from field alignment (0.1-1 degree of aspect spread in the literature), neither modelled | Note `route-one-tilt-map-and-irregularity-orientation` |
| `los_rate_mps`, `tau_c_s`, ... | | As for every echo | | Section 4 |

The heading at which the first-order term reaches the SCI-15 goal (50 m/s, free
space) for `phi_max` = 5 degrees is about 4.3 degrees, and the threshold
(100 m/s) about 8.6 degrees.  The cross-track arrival-angle accuracy that would
remove the term is 0.38 degrees at `psi` 90 degrees, 0.54 at 45, 1.5 at 15,
4-9 at 5 and 11-35 at 2 (note `field-aligned-sheet-separation-result`).  That
measurement needs a second channel and is WP11's.

### Flags and refusals

| Flag or refusal | Set when | Consequence |
|---|---|---|
| `NEAR_MERIDIAN` | `phi_error_mps` <= 50 m/s (about `psi` <= 4 degrees at `phi_max` 5 degrees, free space) | `drift_fa_mps` reported as a drift |
| `OFF_MERIDIAN` | 50 < `phi_error_mps` <= 100 m/s | `drift_fa_mps` reported with the term in its sigma |
| Refused: `ROTATION_UNKNOWN` | `phi_error_mps` > 100 m/s | No drift.  `los_rate_mps` only, with `TILT_DEGENERATE` |
| `SECOND_ORDER_DIP` | Dip above about 20 degrees on the meridian, where `v_s sin I (1 - cos phi_max)` is a large part of the budget | Sigma widened by that term |
| `ALIGNMENT_UNMODELLED` | Always, until the aspect spread and IGRF error are modelled | Reader warned that the orientation error is not in the sigma |
| `ALIAS_RISK` | The residual after step 2 still exceeds `v_unambiguous_mps` | As section 4 |

### What one node cannot give for a sheet

- The rotation `phi` about B.  It cannot be fitted from the Doppler history: the
  fit returned -1.4 degrees against a true 5 degrees, with drift errors of -216
  and -599 m/s at `psi` 15 and 45 degrees (TV Table `c3sheet`).  Across B the
  sheet is a narrow strip, so route one's degeneracy recurs in that plane.
- Hence the drift on tracks far from the meridian.  Route one on its own
  separates only on east-west tracks, where it crosses the strip, and there it
  returns the tilt problem (-265 m/s at `phi` 2 degrees).
- The orbit decides how often a pass is near-meridional (heading minus
  declination), so this product's coverage is a mission-configuration result
  (WP10), not a processing one.

### What Phases 2 and 3 must verify

1. `field_aligned_verify` reproduced through L1-L3 at signal level: drift
   within 1 m/s of the closed form for `phi` = 0, and the closed-form error
   reproduced for `phi` not zero.
2. The classification rule: a sheet on a near-meridional track classed
   `field_aligned_candidate` and not `extended`; a compact facet and an
   extended plate not misclassed as sheets.
3. Noise, fading and the oscillator on top, with sigmas within 20 % of the
   figures above.
4. Not in Phase 3 (recorded as limits): the spread of `phi`, aspect spread,
   IGRF error, finite sheet length along B, magnetoionic propagation.

## 7. The extended layer

For an extended tilted layer no processing of one receive channel separates
the drift from the tilt.  A structure drifting at `(u, w)` under a spacecraft
moving at `v_s` gives the same ranges and Doppler shifts, at every frequency
and time, as a static structure stretched by `u/v_s` and sheared by `w/v_s`
(TV §8.10.2).  The difference in range is `h (w/v_s)^2 / 2`, 11 m at 500 km for
50 m/s, and it is constant.  What differs is the direction of the wave normal
at the spacecraft, by `w/v_s` (0.38 degrees for 50 m/s), and one channel senses
a wave normal only through a reflector-fixed event (route one, facets only) or
through the magnetoionic indices (below).

### Inputs

Those of section 1, plus the ordinary and extraordinary traces of the
ionogram-mode sounding nearest in time (WP4 tracker), and IGRF.

### Processing

1. Steps 1 to 4 of section 2.
2. Classify `extended` (section 3).
3. Form the combined rate (section 4).
4. In the ionogram mode, fit the ordinary trace with a free profile, predict the
   extraordinary trace with the field direction and gyrofrequency, and fit the
   angle between the density gradient and B that best explains the mismatch
   (`ox_grouppath.py` in the screening set; WP6 owns the inversion).
5. If the drift mode visited several frequencies, form the Doppler ratio
   between them and between modes (TV §8.10.1).
6. Set the tilt flag.

### Products

| Product | Units | Definition | Accuracy (parameter level) and conditions | Source |
|---|---|---|---|---|
| `los_rate_mps` | m/s | **The combined rate only**: `w cos(alpha) + v_s (mu_f/mu_sat) sin(alpha)` | As section 4: 8.7-10.0 m/s (INS-26 threshold) or 0.5-1.0 m/s (goal) at 1.15 s, 3-20 MHz | Screen part A |
| `tilt_ox_deg` | deg | Meridional component of the layer's tilt from the consistency of the O and X traces of one ionogram | Bound 0.06-0.20 degrees at 300 m per sample, 0.19-0.67 at 1 km per sample (dip 5-70 degrees, dense topside).  **Bias 0.05-0.29 degrees per 100 m** of range error common to the X samples.  With the present tracker (X predicted from O to 1-2 km) the bias is 0.5-6 degrees.  Exact only for strata tilted as a whole | TV §8.10.1, Table `c3oxcheck`; `ox_grouppath.py` |
| `tilt_ox_sigma_deg` | deg | Bound plus the bias for the tracker's differential range error | | As above |
| `tilt_ratio_deg` | deg | The same angle from the Doppler ratio between modes and frequencies | 0.5-0.7 degrees at a combined rate of 200 m/s, 1.8-2.7 degrees at 50 m/s (four frequencies 2-4 MHz, both modes, 30 dB per sweep over 45 sweeps).  Information proportional to the rate | TV §8.10.1 |
| `tilt_flag` | label | `TILT_LARGE` when `tilt_ox_deg` or `tilt_ratio_deg` exceeds three times its sigma, else `TILT_NOT_DETECTED` | Detects meridional tilts of several degrees, the size a travelling disturbance produces (3-5 degrees, Reinisch et al. 2018).  Does not detect the 0.38-0.76 degrees that matter for SCI-15 | TV §8.10.1 |
| `tilt_term_bound_mps` | m/s | `v_s (mu_f/mu_sat) sin(alpha_max)` for the tilt the flag cannot exclude | 1 degree is 133 m/s at `mu_f/mu_sat` = 1 | Arithmetic |

### The bounded negative

No single-channel method separates the drift of an extended layer at the
accuracy of SCI-15 (TV §8.10.2, Table `c3onechannel`, fifteen methods).  The
node reports the combined rate to a few m/s, the meridional tilt as a flag for
tilts of several degrees, and no drift.  The two-channel methods that do
separate it (dipole plus nadir stub, 0.14 degrees at 30 dB over 45 sweeps;
polarimetry on the crossed pair, 0.04-0.2 degrees within 20 degrees of dip
below 5 MHz) are WP11's.  Route two (a second node) is the WP5 route for the
layer.  The one-channel attitude option, the null scan (C4), was screened on
2026-11-24 and is recorded in section 10.

### Flags and refusals

| Flag or refusal | Set when | Consequence |
|---|---|---|
| Refused: `DRIFT_NOT_SEPARABLE` | Always for `extended` | No `drift_mps`.  `los_rate_mps` with `TILT_DEGENERATE` |
| `TILT_LARGE` / `TILT_NOT_DETECTED` | As `tilt_flag` | `TILT_NOT_DETECTED` does not mean untilted |
| `OX_UNAVAILABLE` | No X trace (below the X cutoff or not tracked), or dip below 5 degrees or above 70 degrees where the bound was not computed | No `tilt_ox_deg` |
| `OX_RANGE_LIMITED` | The tracker's differential O/X range error puts the bias above the bound | Sigma carries the bias |

### What one node cannot give for a layer

- The drift, at any accuracy useful to SCI-15.
- The cross-track tilt (nothing senses it on one channel).
- The along-track tilt better than several degrees with the present tracker.

### What Phases 2 and 3 must verify

1. The extended plate of `c3_verify` case C classed `extended` and refused, with
   `los_rate_mps` equal to the combined rate to 0.1 m/s (noise-free).
2. `tilt_ox_deg` through the forward model is WP6's (field inversion) and
   WP4's (differential O/X range).  Phase 3 checks only that the flag is
   raised for a 5 degree meridional tilt and not for a horizontal layer, given
   a tracker error of 100 m.

## 8. The ionogram mode (INS-09)

The ionogram mode does not deliver SCI-15 (Requirements Rev G; `requirements.yaml`
SCI-15 notes).  Its C3 job is INS-09's ionogram part: the residual
apparent-range bias, `mu_sat Rdot f / K_f`, must be under one range cell after
WP5 processing.  The range cell is 1.0 km above 3 MHz and up to 9.1 km at
0.33 MHz for the A1 stepped sweep.  The correction is computed from the
ephemeris velocity.  Up and down sweep pairing is kept as a bias check only,
because it is not the contribution and because the two visits to a frequency
are separated in time (Phase 5, `novel/sweep_pairing`).

### Inputs

The A1 up-sweep ionogram (and the down-sweep when one is scheduled); the
platform velocity in the Earth-fixed frame at each **step's** time; `f_N,sat`;
the WP4 traces.

### Processing

1. For each step `i` (centre frequency `f_i`, rate `K_i`, centre time `t_i`),
   predict the platform range rate for a wave normal along local down,
   `Rdot_plat = -v_ECEF(t_i) . k0`, and the bias `b_i = direction * mu_sat(f_i)
   Rdot_plat f_i / K_i`.
2. Subtract `b_i` from every range in step `i`.  The geometry is taken at the
   step's time, not once per sweep (section 10, D8: a sweep lasts 4.26 s).
3. If a down-sweep is paired, estimate the bias per trace sample as
   `b_hat = (R_up - R_down) / 2` and compare it with the prediction.  Report the
   difference.  Do not use it as a correction.

### Products

| Product | Units | Definition | Accuracy (parameter level) and conditions | Source |
|---|---|---|---|---|
| `range_corrected_m` | m | Apparent range per step, platform bias removed | Ephemeris term: 0.1 m/s gives `mu f/K` x 0.1 m/s, **0.36 m at 20 MHz**, 0.09 m at 5 MHz | `route_one.ephemeris_term_mps`; arithmetic |
| `bias_per_mps` | m per m/s | `mu_sat f_i / K_i` | 3.57 m at 20 MHz, 1.79 at 10, 0.85 at 5 (`mu` 0.95), 0.46 at 3 (`mu` 0.87); 0.536 `mu_sat` m below 3 MHz, where `K_i` is proportional to `f_i` | A1 law (`sweep_law`) |
| `residual_bias_bound_m` | m | `bias_per_mps` x (drift prior + `v_s (mu_f/mu_sat) sin(alpha_max)`) | The residual the ephemeris cannot remove.  One cell (1 km) is reached at a combined rate of 280 m/s at 20 MHz, 1.2 km/s at 5 MHz.  A 1 degree tilt gives 475 m at 20 MHz, 5 degrees 2.4 km | `requirements.yaml` INS-09 notes; arithmetic |
| `pair_bias_m` | m | `b_hat`, when paired | `sigma_r / sqrt(2)` for a per-trace range error `sigma_r` | Arithmetic |
| `pair_rate_mps` | m/s | `b_hat K_i / (mu_sat f_i)`: the combined rate the pair implies | `sigma_r K / (sqrt(2) mu f)`: with `sigma_r` = 100 m (assumed; WP4 to confirm), **20 m/s at 20 MHz**, 83 m/s at 5 MHz, 152 m/s at 3 MHz.  The up and down visits to a frequency are separated by about 0 at the top of the band and up to two sweep periods at the bottom, so the layer and the geometry change between them | `sweep_pairing` docstring; arithmetic |

### Flags and refusals

| Flag | Set when | Consequence |
|---|---|---|
| `BIAS_RESIDUAL_UNKNOWN` | Always | The drift and tilt part of the bias is not removed.  `residual_bias_bound_m` states its size for the configured prior |
| `CELL_EXCEEDED` | `residual_bias_bound_m` > one range cell at that step | INS-09 not met at that frequency for that prior |
| `PAIR_UNAVAILABLE` | No down-sweep paired | No check |
| `BIAS_CHECK_FAIL` | `abs(pair_bias_m - b_i)` > 3 sigma | Reported for the reader, and as a candidate tilt or drift at the top of the band, where `pair_rate_mps` is best |

### Limits

- The correction removes the platform term for a wave normal along local down.
  For a tilted reflector the platform's projection on the tilt is in the
  residual, as for every echo (`TILT_DEGENERATE`).
- The pairing check is useful only near the top of the band, where the bias is
  largest and the two visits nearly coincide.

### What Phase 5 must verify

1. Noise-free: the corrected range equals the injected group range to 1 m on
   the A1 sweep for platform range rates of +/-1 km/s, with the geometry
   evaluated per step.
2. The pairing estimate recovers the injected combined rate for a horizontal
   layer, and its error grows as the closed form says as the visits separate.
3. INS-09 scored as residual bias in cells against the configured prior.

## 9. What Phases 2 and 3 must verify

1. Noise-free, through L1-L3: `los_rate_mps` equals the injected rate (with
   `w cos(alpha) + v_s (mu_f/mu_sat) sin(alpha)` for a tilted plate) to
   0.1 m/s; `platform_doppler_hz` equals `doppler.doppler_shift_hz` with
   `mu_sat` to 1e-3.
2. `drift_mps` within 1 m/s on the `c3_verify` grid; the tilted-facet bias
   reproduced and `TILT_UNFLAGGED` set; the extended plate declined.
3. With noise, fading and the oscillator: `los_rate_sigma_mps` and
   `drift_sigma_mps` within 20 % of the bound and of the figures above.
4. `tau_c_s` against the injected fading for tau_c 0.1, 0.5 and 2 s, with the
   platform-motion decorrelation of a rough reflector included.
5. `facet_size_m` against the facet size given to `plate_echo`.
6. The several-frequency lever through the forward model (TV §8.21): three
   sub-bands on one tilted facet, drift and tilt recovered within 1 m/s and
   0.01 degrees noise-free, and the screen's sigmas reproduced with noise.
7. `ALIAS_RISK` and `FACET_ALIAS` raised where the tables above say.
8. The class checks of sections 6, 7 and 8.
9. The wideband check of section 10 at signal level: the L2 drift-mode echo's
   sweep-to-sweep phase follows the phase path (`mu_sat Rdot`), not the group
   delay, and the processor converts with the effective frequency of the gated
   window.

M3 is route one alone, on one channel, at signal level (decision of
2026-11-24).  The arrival-direction test is WP11's.

## 10. Screens of 2026-11-24: the wideband check (D8) and the null scan (C4)

Both are stand-alone screening scripts in `Dissertation\Plans\wp5_screening\`
(`wideband_check.py`, `null_scan.py`; outputs `.txt` and `.csv` beside them),
run on this machine with the `tss` environment.  Parameter level.

### D8: is the narrowband model good enough, and how wide can a sub-band be?

The exact echo `s(t - tau(t))`, with the phase following the phase path
(`mu_sat Rdot`) and the envelope the group path (`Rdot / mu_sat`), was
compared through de-chirp with the narrowband model (constant delay per sweep
or step, one Doppler at the step frequency).  The stationary-phase form of the
exact echo agrees with a brute-force frequency-domain sum to 0.05 degrees.
Drift mode: 45 Hz, 1.15 s, window gated after the echo delay, `P'` 500 and
1000 km, combined rates -550 to +1000 m/s, errors scaled to 550 m/s (a
platform residual of 300 m/s plus a 250 m/s drift).

| Result | Value |
|---|---|
| Narrowband model error, phase taken along the phase path | Free space: under 2e-4 m/s at every width to 1 MHz and every frequency 2-20 MHz.  Reference plasma (`f_N,sat` 1.5 MHz): 0.02 m/s at 150 kHz and 0.23 m/s at 500 kHz at 3 MHz (1.0 m/s at 1 MHz); under 0.09 m/s to 1 MHz at 5 and 10 MHz.  `mu_sat` 0.5 (below the drift-mode floor): 0.18 m/s at 50 kHz, 0.73 m/s at 100 kHz.  The error comes from the curvature of `f mu(f)` across the sub-band |
| **Widest sub-band with the model error under 0.5 m/s** | **500 kHz at 3 MHz, at least 1 MHz from 5 MHz up** (reference plasma); unlimited in practice in free space; 50 kHz at `mu_sat` 0.5 |
| Phase error per sweep (model minus exact) | Drift of 3e-3 degrees per sweep at 150 kHz and 0.03 at 500 kHz (reference plasma, 550 m/s); a constant offset of up to 0.9 degrees, which does not affect the rate |
| **Condition 1: L2 carries the phase path** | The stepped construction in `baseline.py` (rows without a dispersion trace) forms the phase as `phi(t) - phi(t - tau)` with `tau` the group delay.  If a drift-mode sequence updates `tau` per sweep from L1, the sweep-to-sweep phase moves at the group-path rate, wrong by `(1/mu_sat^2 - 1)` of the rate: **183 m/s at 3 MHz, 54 at 5, 12.7 at 10** (550 m/s), 1.65 km/s at `mu_sat` 0.5, at any width.  Zero in free space.  To be fixed or confirmed in Phase 2 |
| **Condition 2: the processor converts at the effective frequency** | The gated window starts one echo delay into the sweep, so the range-gated phase belongs to the transmitted frequency at the window centre, `f_eff = f_start + K (u_c - tau)`, which depends on range.  Converting at the sub-band centre instead biases the rate by about `0.15 (B_sub/f)/mu_sat^2` of it at 1000 km.  That limits the sub-band to 10 kHz at 2-3 MHz, 20 at 5, 50 at 10 and 100 at 20 MHz.  Converting at `f_eff` removes it |
| Range walk (arithmetic) | Over 1.15 s at 550 m/s the echo moves 630 m (more by `1/mu_sat` in group range), one range cell at `B_sub` of about 240 kHz.  Wider sub-bands need range-walk compensation (keystone) or shorter dwells |
| Ionogram mode, A1 sweep, +/-1 km/s | Narrowband with the geometry per step against exact: 0.14 m (free space), 0.43 m (`f_N,sat` 1.5 MHz).  The bias formula `mu f Rdot / K` holds to 0.61 m.  **The full-sweep ionogram numbers do not change** (cell 1 km) |
| Ionogram mode, geometry frozen per sweep | The analytic tier takes one `PlatformState` per call, so every row of a 4.26 s sweep is at one epoch.  Against per-step geometry that is `Rdot t_i`: 1.27 km at 300 m/s and 4.25 km at 1 km/s at the top of the band.  Small for the platform's own radial rate (tens of m/s), several cells for the combined rates of a tilted layer |

**Recommendation for the Phase 2 sweep law:** `B_sub` = 150 kHz (one A1 step,
1.0 km cell, range walk under one cell over 1.15 s); the narrowband model is
then within 0.02 m/s of the exact echo down to 3 MHz.  Phase 2 must (1) make the
sweep-to-sweep phase follow the phase path, and (2) convert phase rate to
velocity at `f_eff`.  Phase 5 must evaluate the ionogram geometry per step, or
score INS-09 against truth at the sweep epoch.

### C4: the null scan on one channel

The spacecraft is rolled so one dipole of the crossed pair points to nadir in
the orbit plane and the other lies along track.  It transmits on the horizontal
dipole and receives on the near-nadir one (still crossed, so INS-23 isolation
is kept), and pitch dither or a slow spin scans the receive null through the
echo's wave normal.  The voltage goes as `sin(beta - alpha)`, so the null's
position gives the along-track wave-normal angle `alpha`, and the drift error
is `v_s sigma_alpha` with no plasma factor.  Heritage: ISIS 1 (James 1980) and
DE 1 (Sonwalkar and Inan 1986) used spin modulation of a dipole for arrival
and wave-normal directions.  Free-space transverse field; the plasma's
longitudinal component and the O/X polarisation are not modelled.

| Result | Value |
|---|---|
| Bound without fading | `1/sqrt(2 SNR N)` rad, almost independent of the dither amplitude: at 30 dB per sweep, 29 m/s over 1.15 s, 16 m/s over 10 s, 4 m/s over 40 s.  At 20 dB, 93, 51 and 13 m/s |
| With fading (AR(1), `tau_c` 0.1, 0.5, 1.0 s, re-run 2026-11-24 evening with `decorrelation_update.py`; the original 0.1-2 s screen is in `null_scan.txt`), ML estimator | Pitch +/-2 to +/-5 degrees at 0.05 Hz, 40 s dwell: **23-32 m/s at 20 dB, 8-10 m/s at 30 dB**.  Pitch +/-5 degrees at 0.2 Hz, 10 s: 49-59 m/s at 20 dB, 17-24 at 30 dB.  1.15 s dwell (1 Hz): 125-167 m/s at 20 dB, 48-59 at 30 dB.  Spin 36 degrees/s over 20 s: 58-62 m/s at 20 dB, 21-28 at 30 dB.  The fading of an extended layer can come from the medium or the platform's motion; the screen does not distinguish them.  A decorrelation time wrong by a factor of 2 costs under 1.5x; the residual Doppler must be given to the estimator (ignored: 535-639 m/s) |
| Loss of echo | Mean echo power on the receive dipole 24 dB (+/-5 degrees) to 32 dB (+/-2 degrees) below an aligned dipole; 9 dB for a spin.  The echo must be found in the ionogram mode first |
| Attitude (moments of inertia assumed: 0.6-2.2 kg m^2 about pitch, bus plus two 5 m dipoles at 20-100 g/m) | +/-2 degrees at 0.05 Hz: 0.63 degrees/s, **4 mN m (2-8), 13 mN m s (6-25)**, within a typical 6U wheel (2-10 mN m, 20-100 mN m s, assumed).  +/-5 degrees at 0.05 Hz: 10 mN m (5-19), 33 mN m s (16-61), at or above the top of that range.  1 Hz dithers need 0.8-15 N m: not feasible.  Spin: 0.06-3.5 N m s and 7-105 h of magnetorquer spin-up; the whole mission spins |
| Verdict | **Reaches the SCI-15 goal (50 m/s) and threshold (100 m/s) for an extended layer from 20 dB per sweep**, with a slow pitch dither of about +/-2 degrees at 0.05 Hz over 40 s dwells, within a 6U wheel's capacity.  It does not reach them at the 1.15 s dwell without a 1 Hz dither, which the attitude system cannot provide.  Costs: 300 km of track per dwell (the tilt is averaged over it), the echo 24-32 dB down during the scan, dwells taken from the ionogram cadence, pitch knowledge and the dipole's effective axis to about 0.1 degrees (0.1 degrees is 13 m/s), boom dynamics near 0.1 Hz not modelled, magnetoionic polarisation not modelled |

The null scan is the only single-channel method found that measures the
wave-normal direction of an extended layer, so it qualifies the bounded
negative of section 7: no **processing** of one channel separates drift from
tilt, but one channel and the attitude system can, at a cost.  It is an
arrival-direction method, and Andrew's decision of 2026-11-14 puts
arrival-direction work in WP11.  Where C4 goes is recorded as open question 5.

## 11. Open questions for Andrew

Decided on the evening of 2026-11-24 (recorded under "Decisions this contract
rests on"): 1, SCI-15 is scored against the plasma's own rate; 2, the
unambiguous interval is limited to +/-1 km/s where sampling supports it, with
sub-bands a possible future improvement; 3 and 4, `tau_c_s` and
`facet_size_m` are reportable data products.  Still open:

5. **Where the null scan (C4) goes.**  It is arrival-direction work (WP11 by
   the decision of 2026-11-14) but needs only one channel and the attitude
   system.  Options: keep it in WP11 with the two-channel methods; or make it a
   WP5 Phase 4 item (a signal-level test with the magnetoionic polarisation and
   an attitude-knowledge error).  Recommended: WP11, with the screen recorded in
   TV §8.10.2, because its open questions (boom dynamics, polarisation,
   effective axis) are WP11's calibration questions.
6. **D8 condition 1.**  Confirm in Phase 2 that a drift-mode sequence carries
   the sweep-to-sweep phase along the phase path (section 10).  If it does not,
   the simulator overstates the plasma Doppler by `1/mu_sat^2`.
7. **INS-13 product set.**  `tau_c_s` and `facet_size_m` are now data
   products.  The INS-13 product sets and the data-volume estimate
   (`docs/data-volume.md`) do not list them.  Adding them is a requirements
   amendment, which needs Andrew's approval; they add a few numbers per echo.

## Possible future improvement (not adopted)

In the plasma-velocity mode, two or more drift-mode sub-bands at different
frequencies have different unambiguous intervals, `c f_rep / (4 f mu_sat)`, and
together resolve a wrapped rate by the usual staggered-frequency argument.
That would extend the +/-1 km/s interval (and the +/-1.5 km/s goal) across the
band above `f mu_sat` of 5.3 MHz.  The several-frequency lever already
schedules several sub-bands for facets, so the cost is mainly in the
scheduler.  Recorded 2026-11-24; not part of the WP5 baseline.

## Sources

- Studies: `c3_verify`, `route_one`, `route_one_tilt_map` (this repository);
  `tss.l3_signal.novel.route_one_estimator`; `tss.l1_propagation.doppler`.
- TV Ch. 8: Tables `c3r1budget`, `c3tiltallow`, `c3multifreq`, `c3onechannel`;
  §8.9.3, §8.10.2.
- Project notes: `route-one-study-and-ch8-plasma-values` (2026-10-06),
  `c3-verify-route-one-tilt-finding` (2026-10-08),
  `route-one-tilt-map-and-irregularity-orientation` (2026-10-08),
  `wp5-single-vehicle-trade-study` (2026-11-14),
  `wp5-decisions-and-p0-start-2026-11-24`.
- Screening: `Dissertation\Plans\wp5_screening\single_node_gaps.py`,
  `multifreq_lever.py`, `wideband_check.py` and `null_scan.py`;
  `ox_grouppath.py` (O/X consistency).
- Studies for part 2: `field_aligned_verify`
  (`l3_signal/novel/field_aligned_estimator.py`), TV Tables `c3sheet`,
  `c3oxcheck`, `c3onechannel`; `requirements.yaml` INS-09 and SCI-15 notes.
- James (1980), JGR 85, doi:10.1029/JA085iA07p03367; Sonwalkar and Inan
  (1986), JGR 91, 154-164, doi:10.1029/JA091iA01p00154.
- Dyson (1969), JGR 74, 6291 (facet enhancement 5-70 %).
