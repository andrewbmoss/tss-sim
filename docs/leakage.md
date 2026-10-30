# Transmit-to-receive leakage (WP3)

Module: `tss.l2_instrument.leakage`.  Tests: `tests/unit/test_leakage.py`.
Plan: project note `claude/wp3-plan.md`.  Reference echo: `docs/reference-echo.md`.

**Status: WP3 closed 2026-10-29.**  Claim C4 closes in simulation.  The
transmitter's AM noise, not phase noise or isolation, sets the residual; the
budget closes with a passive RF canceller of at least 26 dB whose weight moves
continuously (INS-31), AM of at most -140 dBc/Hz (INS-30) and 35 dB of isolation
(INS-23), stated as leakage-borne noise no greater than k T_sys (INS-29).  The
two hardware allocations are measured by the WP9 conducted bench (study
`wp9_bench`; WP9 Bench Test Specification Rev A).  Phases 1 to 7 follow.

## Phase 1 (2026-10-27): the model

Decisions: G1, parametric isolation brackets (30, 40, 50, 60 dB, with ripple and
resonances) as the scored set and one moment-method anchoring run; G2, coupling
paths of 1, 3, 10 and 30 m (3.3 to 100 ns), with the baseline left at 3 ns
(Open Item 13).

### What it computes

`leakage_signal(sweep, fs, model, rng, ...)` returns the leakage at the receiver
input as de-chirped complex samples in root-watts, with the phase convention of
the L2 baseline's echoes, `phi_ref(t) - phi_rx(t)`:

    x(t) = sqrt(P_tx) conj(H(f_e)) exp(j [Dphi(t) + e_ref(t) - e_tx(t - tau)])

* `Dphi = phi(t) - phi(t - tau)` of the transmitted law, computed in closed form
  inside a step (full precision for a nanosecond delay) and as the exact
  difference across a step boundary, so the transition is in the samples.
* `H(f)`: the isolation profile.  A flat or tabulated base, a sinusoidal ripple
  in dB, parallel resonances (which turn the coupling phase), and an optional
  retuning transient at each step boundary.
* `tau(t)`: the coupling delay, with optional first-order Gauss-Markov jitter.
* `e_tx`, `e_ref`: transmit and reference phase errors.  `"common"` is INS-21;
  `None` is an independent reference and loses the range-correlation credit.
  Phase 1 carries the single-tone FM of Rev G section 5.9; Phase 2 adds the
  oscillator, spurs, AM noise and sweep nonlinearity.

`leakage_over_echo_db(instrument, sweep, model)` gives leakage over the Phase 0
echoes at each step with the model's profile; with a flat model it is the
Phase 0 reference exactly.

### Against the analytics (exit criterion 0.5 dB)

| Form | Cases | Largest error |
|---|---|---|
| Sideband after de-chirp, exact Bessel `J1/J0` of `(2 Df / f_m) sin(pi f_m tau)` | 4 paths x 3 offsets (3.74, 20, 39.2 kHz) | 0.008 dB |
| Suppression of a transmit spur, `4 sin^2(pi f_m tau)` | same 12 | 0.011 dB |
| Rev G's small-angle sideband, `(pi Df tau)^2` | same 12 | 0.008 dB |
| Independent reference: the spur survives whole | 10 m, 20 kHz | within 0.5 dB |
| Same, inside one A1 step at 5 MHz | 10 m | within 0.5 dB |
| Isolation ripple: AM sidebands at `K / period`, `I1(z) / I0(z)` | 0.5, 1, 3 dB | < 0.001 dB |
| Delay jitter: `2 pi f sigma_tau` rms phase | 2 ps at 20 MHz | within 0.5 dB |
| Leakage power over the Phase 0 echo at 1050 km | 4 isolations x 3 steps | within 0.5 dB |

### The moment-method anchor (G1(b))

`scripts/matlab/wp3_isolation_mom.m`, output in
`data/reference/wp3_isolation_mom` (README there), loaded by
`mom_anchor_model()`.  Two 5 m dipoles, coil Q 150, feeds 0.366 m apart, free
space, no bus.

* Below a few MHz both antennas are dominated by their coil loss and the pair
  behaves as two coupled resonators.  The isolation is then set by the
  coupling coefficient `k = |Z12| / |X|` times the loaded Q, not by antenna
  efficiency, and it is flat with frequency.
* Crossed dipoles with the receive feed 3 and 8 cm off the axis: k = 4.7e-4,
  29.0 dB from 0.33 to 3 MHz, rising to 35 dB at 11 MHz and 49 dB at 20 MHz.
* Parallel dipoles: k = 0.20, 0.3 dB with both ports matched together, 23.5 dB
  with each matched alone.
* On this model INS-23's 40 dB needs k below about 1.3e-4 at Q = 150.  That
  bears on Committee Review 5's point that 40 dB is assumed, not shown.  It is a
  free-space model of one assumed geometry, and the bus is not in it; WP3
  Phase 5 uses it as one profile among the brackets.

### The anchor on the bus (2026-10-27)

`scripts/matlab/wp3_isolation_mom_bus.m` adds a closed conducting 6U box and
compares mounts; see the README in `data/reference/wp3_isolation_mom`.  The
isolation obeys `coupled_resonator_isolation_db(k, Q)`, `-20 log10(kQ/2)`:
36.6 dB for crossed dipoles with a 3 and 8 cm feed offset on the bus
(k = 2.0e-4), 0.4 dB parallel.  The symmetric crossed mounts are limited by
their imperfections: a 0.1 deg angle error gives 36 to 40 dB, 1 deg 14 to
21 dB.  `mom_anchor_model()` now defaults to the bus mount.

## Phase 2 (2026-10-27): transmitter impairments

Module: `tss.l2_instrument.impairments`.  Tests: `tests/unit/test_impairments.py`
(19).  Configuration: the `transmit_impairments` block of
`config/instrument/baseline_6u.yaml`, all off by default.

`TransmitImpairments` switches each term on or off:

| Term | Shared with the reference under INS-21 | Prediction tested | Largest error |
|---|---|---|---|
| Reference-oscillator noise (Allan deviation and noise type, scaled to the carrier) | yes | `L(f) 4 sin^2(pi f tau)` at a 7 ms echo delay, and `2 L(f)` with an independent reference | 0.20 dB |
| Synthesiser phase-noise floor | yes | `L 4 sin^2(pi f tau)`, 3 and 30 m | 0.06 dB |
| INS-21 off | - | the credit lost: 50.43 dB measured against 50.47 dB predicted (10 m, 15-25 kHz) | 0.10 dB |
| AM noise | no | `L_AM` whatever the path or INS-21 | 0.008 dB |
| PM spur / AM spur | yes / no | exact Bessel sideband after de-chirp / unchanged | within 0.5 dB |
| Sweep nonlinearity (INS-25) | yes (ideal law if independent) | in-step residual phase in closed form; credit above 60 dB at 30 m | 1e-9 rad |
| Doppler from `EchoSet.doppler_hz` (mu_sat in L1) | - | beat moves by `-f_d` on up- and down-sweeps | 0.5 Hz |

One `ImpairmentRealization` per record holds the draws; the echoes in the L2
baseline and the leakage take the same one, as one oscillator would.  Noise is
synthesised in the frequency domain on a grid padded by the longest delay, so a
delay is an exact phase ramp; the oscillator term is the running integral of
the delayed `f_c y`.  The de-chirped phase of a nonlinear law is computed from
each step's local cycles, so it keeps full precision for nanosecond paths.

With the baseline part (flicker FM, 2e-9 at 1 s) at 10 MHz, the residual at a
7 ms echo delay is about -115 dBc/Hz from 300 to 700 Hz whether or not the
reference is shared: range correlation helps only below about 1/(6 tau),
24 Hz at 1050 km.  On the leakage path it is worth 50 dB and more.

## Phase 3 (2026-10-28): the leakage in the samples and the receive chain

Decisions (2026-10-28): leakage off by default until Phase 4 scores it; when
on, the MoM bus profile (`bus_crossed_ends_offset`); the high-pass corner
follows each step; INS-11 scored on the noise-referenced rule (G5b), Rev G's
rule reported beside it.

Module: `tss.l2_instrument.receive_chain`.  Study: `leakage_receive_chain`
(outputs in `data/reference/wp3_receive_chain`).  Tests:
`tests/unit/test_receive_chain.py` (12).  Configuration: the `leakage`,
`front_end` and `leakage_highpass` blocks of `baseline_6u.yaml`.

* **In the samples.**  With `leakage.enabled`, the L2 baseline adds
  `leakage_signal(...)` at the receiver input, before the noise, drawn from the
  record's impairment realization, so echoes and leakage share one oscillator.
* **Front end.**  A Rapp limiter on the composite's envelope, input P1dB
  +1 dBm (a level-7 mixer, assumed), optional AM/PM.  Applying it to the
  de-chirped envelope is checked against a passband simulation of a real tanh
  nonlinearity followed by the de-chirp: the leakage, the echo and the echo's
  image agree to 0.1 dB.  A weak echo beside the leakage loses
  `|a|^2` and gains an image `|b/a|^2` at the negative beat
  (`small_signal_gain_db`, `compression_image_dbc`), checked to 0.05 dB.
* **Leakage high-pass.**  Fourth-order Butterworth in state-variable sections,
  trapezoidal, corner `0.5 x 2 K_i R_min / c` per step (205 Hz at 0.33 MHz to
  1.87 kHz above 3 MHz in A1; to 11.8 kHz in A2).  It matches the analog
  response to 0.05 dB; a 100 km echo loses 0.017 dB.  Its states carry across a
  change of corner, and a 5 % change shakes loose less than -100 dBc of an
  absorbed near-DC leakage.

### What the study finds (A1, MoM bus isolation)

| | Transmitter starts with the sweep (cold) | Transmitter already on (primed) |
|---|---|---|
| INS-11, noise rule, least margin | 38.6 dB | 92 dB |
| Rule asks at most (noise / Rev G) | 60.7 / 96.8 dB | same |
| Longest settling | 25.1 ms, first step (0.94 of a dwell) | 1.65 ms |
| Steps settling longer than the 7 ms echo transition | 1 (the first) | 0 |
| Filtered leakage peak against full scale | +50 dB (the converter clips in step 1) | -43 dB |

* Without the high-pass every sample clips: INS-11 needs it.
* In steady state the filter removes the carrier by about 200 dB; the rules
  ask 61 dB (noise) or 97 dB (Rev G) at most.  What decides the first step is
  the leakage switching on: if the transmitter starts with the sweep, step 1
  (0.33 MHz) is lost for 25 ms.  Running the transmitter before the sweep, or
  blanking step 1, removes it.
* Front end: with the bus isolation the leakage is -6.7 dBm below 3 MHz,
  7.6 dB under the assumed +1 dBm P1dB, short of WP9-2's 10 dB in 83 of 159
  steps.  An echo there loses 0.11 dB and gains an image 41 dB down at the
  negative beat: a strong echo 60 dB above the noise would put a ghost 19 dB
  above it.  At flat 40 dB isolation the margin is 11 dB and the image 55 dB
  down.

## Phase 4 (2026-10-28): the residual and its score

Decisions (2026-10-28): INS-24 scored against the weakest admitted echo, noise
+ 10 dB (G5a), with the R_max reading beside it; scored without a canceller,
an ideal 30 dB bound beside it (G3); the default transmitter carries its
allocations at threshold; the transmitter runs before the sweep (primed).
Study: `leakage_residual` (`data/reference/wp3_residual`, pinned to its own
configuration).

* `dechirped_residual(...)`: per step, a Hann transform over the dwell in the
  echo band (100 km to R_max).  The **pedestal** is the largest of up to four
  sub-band means; a **line** is a bin above `ln(N / 1e-3)` times its part's mean,
  confirmed at the same offset in a neighbouring step (without the confirmation
  about one step per sweep carried a false line); the **level** is the larger.
* `score_ins23_ins24(...)`: INS-23 is the least isolation over the sweep;
  INS-24 the worst-step margin of the level below the weakest echo, with the
  R_max margin, the canceller bound and the number of failing steps.
* `leakage_signal(primed=True)` and `LeakageHighpass.steady_state` start the
  record with the filter settled on the carrier.
* The L2 baseline puts `l2_metrics` on every record and the per-step residual in
  `flags["leakage_per_step"]`; L5 maps INS-11, INS-14, INS-23, INS-24 and INS-26.
  Leakage is on by default from Phase 4.

**Finding.**  The pedestal is the transmitter's AM noise carried on the
leakage, `P_leak + L_AM + 10 log10 B_cell`, reproduced to 0.5 dB.  Range
correlation leaves the phase noise of a 3 ns path about 36 dB below the cell
noise but does nothing for AM.  The cell noise falls toward the band bottom with
the dipole's efficiency, so at -133 dBc/Hz the pedestal stood 9 dB over the
noise at 5 MHz and 27 dB at 0.338 MHz: INS-24 failed by 26.7 dB.  A flat 40 dB
isolation still left -13.5 dB; a cold start took INS-11 to 29.9 dB.

## Phase 5 (2026-10-28): the leakage budget and the canceller

Study: `leakage_budget` (`data/reference/wp3_leakage_budget`).  Tests:
`tests/unit/test_leakage_budget.py`.

`RFCanceller` (in this module) subtracts a weighted sample of the PA output
ahead of the front end, so it removes the leakage's AM with its carrier.  The
residual factor is `g = 1 - w [H(f_k)/H(f)] exp(-j 2 pi (f - f_k) dtau)`, with:

* calibration errors of amplitude, phase and delay: RMS per calibration,
  correlated over `error_correlation_steps`, or a fixed bias
  (`error_draw="worst"`);
* `retune`: `track` (a C2 spline through the per-step calibrations), `step`
  (switched at each boundary) or `sweep` (one setting at `calibration_hz`);
* `added_noise_dbm_hz` for an active vector modulator; `ideal_depth_db` for a
  bound.

The study runs the full L2 chain at two AM levels (-133 and -153 dBc/Hz) with
matched draws, which splits the residual exactly into the AM part `X` and the
rest `Y` per unit leakage; the grid over AM, isolation, power and canceller
then follows by scaling, `ped = P_leak (X 10^((AM+133)/10) + Y)`.  Sixteen
grid points re-simulated directly agree to 0.05 dB.

| Canceller (A1, bus isolation, 1 W) | AM that meets INS-24 (dBc/Hz) | INS-24 at -140 (dB) |
|---|---|---|
| none | -159.7 | -9.7 |
| ideal 26 dB (INS-31 threshold) | -133.9 | +16.1 |
| 0.1 dB / 1 deg, tracked (33.6 dB) | -133.8 | +16.2 |
| 0.3 dB / 2 deg, tracked | -140.4 | +9.6 |
| 0.5 dB / 3 deg, tracked | -144.0 | +6.0 |
| one setting at 1 MHz | -138.5 (A2 -148.1) | +11.5 |
| switched at each step | none suffices | -0.4 |
| independent errors per step | none suffices | -0.2 |
| active vector modulator (-158 dBm/Hz) | none suffices | +1.6 |

* Dynamics matter as much as depth: boundary jumps of 0.001 / 0.003 / 0.01 dB
  RMS leave INS-24 at +12.8 / +3.5 / -6.9 dB, so the weight must move
  continuously (steps no larger than 0.001 dB or 0.0066 deg, -79 dB of the
  leakage) and its errors must be smooth across the band.
* The canceller must be passive.  Lowering the power does not close the budget
  (0.1 W still needs -150 dBc/Hz without a canceller); isolation buys it dB
  for dB.
* Behind the canceller the front end is 38.6 dB clear of a level-7 mixer's
  P1dB and the echo image vanishes.
* One settling dwell before the first step is equivalent to priming.

Rev G gained INS-29 (leakage-borne noise <= k T_sys), INS-30 (AM <= -140 dBc/Hz,
goal -150) and INS-31 (passive canceller); all three are scored
(`leakage_noise_over_system_noise_db`, `transmit_am_noise_dbc_hz`,
`rf_canceller_min_depth_db`).  The default `baseline_6u.yaml` carries AM at
-140 dBc/Hz and an ideal 26 dB canceller (the `canceller` block).

INS-23 was relaxed to 35 dB on 2026-10-29 (`INS23_THRESHOLD_DB`;
`BASELINE_ISOLATION_DB` stays 40 as the analyses' bracket): with INS-31 the
leakage is at least 61 dB below the transmitter, INS-29 holds by 4.7 dB and
INS-22 by 21.9 dB.  Smoke campaign at the default: INS-24 +16.1 dB, INS-29
-6.1 dB, INS-10 25.9 dB over the scored echoes (each mode's cutoff step,
reported separately under G5c, at 18.9 dB).

## Phase 6 (2026-10-29): the WP9 bench predictions

Study: `wp9_bench` (`data/reference/wp9_bench`).  Tests:
`tests/unit/test_wp9_bench.py`.  The conducted bench emulates INS-31 by extra
attenuation in a flat leakage path from the PA output, which lowers the AM
with the carrier as a canceller does: E0 39 dB (no canceller; a level-7 mixer at
WP9-2's 10 dB margin), E1 61 dB (35 + 26), E2 69 dB (35 + 34), each through the
RG-58 delay set (0.66 to 19.8 m).  The study gives, per step, the INS-29 limit
in dBm/Hz, the predicted leakage-borne noise, the rise of an assumed
-165 dBm/Hz bench floor, the front-end margin and image, the stray-coupling
limits and the AM measurement floor.

* Behind E1 at -140 dBc/Hz: INS-29 -4.5 dB at 0.338 MHz, a floor rise of
  1.0 dB; at the INS-29 limit the floor would rise at least 2.4 dB.  E1 passes
  for AM up to -135.5 dBc/Hz.
* The delay set changes the residual by under 0.6 dB at -140 dBc/Hz, up to
  3.4 dB above 7 MHz at -150 on 19.8 m, where phase noise begins to show.
* AD9959 crosstalk (> 65 dB) limits a synthesiser-channel echo to about
  -86 dBm behind E1, so weak echoes with leakage present come from the DRFM.
* AM measurement floor by in-phase homodyne detection: about -159 dBc/Hz in one
  chain, -174 by a two-chain cross-spectrum of 1000 averages.

## Phase 7 (2026-10-29): close-out

* `leakage_budget` scores INS-10 without each mode's cutoff-step echoes, as L5
  has since G5c; the reference run is regenerated (ideal instrument 24.9 dB, no
  echo below 10 dB at any power in the grid).
* The Technical Volume's section 5.5 is rewritten around these results, with
  figures from `data/reference/wp3_isolation_mom` and
  `data/reference/wp3_leakage_budget`; Rev G and the SRM point the status of
  INS-21 to INS-24 and INS-29 to INS-31 at the bench tests.

## Not modelled, and where it goes

* The plasma sheath, booms, panels, harness and matching-network coupling in
  the isolation (Open Item 13).
* AM noise that is not white across an echo band; the canceller as a circuit
  rather than a weight with errors.  Both are bench items (WP9-5, WP9-8).
* `tss.l3_signal.leakage_suppress` (WP5) is no longer needed for C4: the
  suppression is analog (INS-31).
