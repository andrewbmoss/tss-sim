# The drift mode and what one node delivers

WP5, claim C3.  This is the output contract for the drift-mode L3 tier
(`l3/drift`, planned for WP5 Phase 3) and for the ionogram-mode INS-09
correction (Phase 5).  It says, for each class of echo, what a single
spacecraft with **one receive channel** delivers, how well, under which
conditions, and which flags it sets when it cannot deliver.  The code that
implements it does not exist yet; `DriftProduct` in `interfaces.py` is written
from this document in Phase 2.

Status: **part 1 of 2** (2026-11-24): the products common to every echo, and
compact facets.  Part 2 adds field-aligned sheets, the extended layer and the
ionogram mode.

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
| `field_aligned_candidate` | To be defined in part 2 (geometry predicted from IGRF and the track's angle to the magnetic meridian) | Part 2 |
| `unresolved` | Power envelope narrower than about 4 sweeps, or SNR below the detector threshold, or contrast < 4 with a sweeping slope | `route_one_tilt_map` (unresolved envelope); detector P_fa 1e-4 |

The label is an output, with the slope ratio and the contrast that decided it.

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
| `range_m` | m | Apparent range `P'` of the echo per sweep (range gate) | Resolution `c / (2 B_sub)`; set when D8 fixes the sub-band width `B_sub`.  To be established in Phase 2 | Phase 2 |
| `power_db`, `snr_db` | dB | Echo power and SNR per sweep and at the peak | SNR against the detector's noise reference | WP4 P2 |
| `platform_doppler_hz` | Hz | `f_D,plat` removed in step 1 | Ephemeris velocity to 0.1 m/s gives 0.1 m/s | `route_one.ephemeris_term_mps` |
| `doppler_hz` | Hz | Residual Doppler per sweep pair | As `los_rate_mps` | Step 3 |
| `los_rate_mps` | m/s | **Combined line-of-sight rate**, `c (f_D - f_D,plat) / (2 f mu_sat)`, averaged over the dwell.  For a reflector tilted by `alpha` along track it is `w cos(alpha) + v_s (mu_f / mu_sat) sin(alpha)`: drift plus the platform's projection on the tilt, which one channel cannot remove | Dwell 1.15 s (INS-15 goal): **8.7-10.0 m/s** with the oscillator at the INS-26 threshold, **0.5-1.0 m/s** at the INS-26 goal, 3-20 MHz.  Dwell 0.3 s: 8.7-11.9 m/s (threshold), 1.0-6.5 m/s (goal).  Dwell 0.05 s (INS-15 threshold): 15-95 m/s, Doppler term dominant, worst at 3 MHz | Screen part A; the oscillator term is `route_one`'s; matches the Doppler and oscillator terms of TV Table `c3r1budget` |
| `los_rate_sigma_mps` | m/s | 1-sigma of `los_rate_mps` from the terms above | | |
| `resolution_mps` | m/s | Separation at which two scatterers in one range cell are resolved in Doppler, `c / (2 f mu_sat T)` | Dwell 1.15 s: 50, 27, 13 m/s at 3, 5, 10 MHz.  The SCI-15 goal of 50 m/s between scatterers is met at 1.15 s from 3 MHz up, which is how INS-15's goal was derived | Arithmetic; INS-15 |
| `v_unambiguous_mps` | m/s | `c f_rep / (4 f mu_sat)` | 45 Hz: 1298, 707, 341, 169 m/s at 3, 5, 10, 20 MHz.  71 Hz: 2048, 1116, 538, 267 m/s.  The +/-1 km/s of INS-09 and SCI-15 is met by sampling alone only where `f mu_sat` <= 5.3 MHz at 71 Hz; above that, ambiguity has to be resolved (open question 2) | Screen part B |
| `range_debiased_m` | m | `range_m` corrected for the measured rate in the drift-mode sweep, bias `mu_sat Rdot f / K_sub` | Follows `los_rate_sigma_mps` | `doppler.apparent_range_bias_m` |
| `doppler_spread_hz` | Hz | Spectral width of the residual after the mean Doppler (and, for a facet, the ephemeris slope) is removed | To be established in Phase 3.  Used to flag several scatterers in one cell or a rough reflector | Phase 3 |
| `tau_c_s` | s | **The echo's own decorrelation time**, from `|R(1)| / (R(0) - noise)` of the complex echo (the modulus removes the constant Doppler, so no Doppler estimate is needed) | 10 dB per sweep, 45 Hz: tau_c 0.1 s within a factor 0.69-1.41 (16-84 %) over a 1.15 s dwell and 0.77-1.24 over 3 s; tau_c 0.5 s within 0.43-2.9 (1.15 s) and 0.58-1.8 (3 s); tau_c 2 s gives only a lower bound.  Rule: an estimate when `tau_c` <= about half the dwell, otherwise a lower bound.  Not modelled: decorrelation from the platform moving over a rough reflector, which biases `tau_c` short and has to be separated at signal level | Screen part G |
| `echo_class` | label | Section 3 | | |

`tau_c_s` matters beyond the echo.  The topside decorrelation time gates route
one (WP5 decision D3) and nothing measures it.  A sounder that reports it for
every echo measures it in flight, at the frequencies and heights it sounds.

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
| `ALIAS_RISK` | Predicted residual rate (drift prior plus the tilt allowance) beyond `v_unambiguous_mps` | Rate may be wrapped; resolve by open question 2 |
| `SPREAD` | `doppler_spread_hz` above the single-scatterer value | Several scatterers or a rough reflector; class products are not formed |
| `TILT_DEGENERATE` | Always, unless a class product removes it | `los_rate_mps` is drift plus `v_s (mu_f/mu_sat) sin(alpha)`; it is not the plasma's own rate |

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
| `facet_size_m` | m | Along-track size of the specular region, `L = 0.443 lambda P' / (W v_s)` from the fitted envelope width `W` | **3-9 %** (1-sigma) without fading for W 0.44-2.9 s; **14-48 %** with fading at tau_c 0.1-0.5 s, biased low by 2-27 % (fading narrows the fitted envelope).  It measures the coherent specular region, which a curved or rough facet makes smaller than the facet | Screen part E; `route_one.envelope_width_s` |
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

Part 2.

## 7. The extended layer

Part 2.

## 8. The ionogram mode (INS-09)

Part 2.

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
6. The several-frequency lever through the forward model (TV §8.21).
7. `ALIAS_RISK` and `FACET_ALIAS` raised where the tables above say.

## 10. Open questions for Andrew

1. **What SCI-15 is scored against.**  `requirements.yaml` scores SCI-15 as
   the line-of-sight velocity of an isolated echo after the
   ephemeris-predicted platform term is removed.  Read literally, that is
   `los_rate_mps`, which meets threshold and goal (8.7-10 m/s at 1.15 s).  But
   for a tilted reflector it contains `v_s (mu_f/mu_sat) sin(alpha)`, which is
   platform motion, not plasma motion.  Should the L5 truth be the apparent
   rate (what the requirement literally says) or the plasma's own rate along
   the wave normal (what the science needs)?  The second is met only for
   facets, through route one and the lever.
2. **Velocity ambiguity above `f mu_sat` of about 5.3 MHz.**  The drift mode
   cannot sample +/-1 km/s unambiguously there even at 71 Hz.  Options: two
   sub-bands at different frequencies resolve the wrap by their different
   unambiguous intervals (the several-frequency lever schedules several
   sub-bands anyway); or restrict the +/-1 km/s claim to the lower
   frequencies.
3. **Adopt `tau_c_s` as a data product.**  It would measure the topside
   decorrelation time in flight.  It adds a few numbers per echo to the INS-13
   product.
4. **Adopt `facet_size_m`.**  It comes free from the envelope fit.

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
- Screening: `Dissertation\Plans\wp5_screening\single_node_gaps.py` and
  `multifreq_lever.py`.
- Dyson (1969), JGR 74, 6291 (facet enhancement 5-70 %).
