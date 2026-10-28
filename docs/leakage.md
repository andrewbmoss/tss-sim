# Transmit-to-receive leakage (WP3)

Module: `tss.l2_instrument.leakage`.  Tests: `tests/unit/test_leakage.py`.
Plan: project note `claude/wp3-plan.md`.  Reference echo: `docs/reference-echo.md`.

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

### Not yet

* The leakage is not added to the L2 baseline's samples; receiver compression
  and the leakage high-pass are Phase 3.  When it is, it takes the record's
  realization (`leakage_signal(..., realization=...)`).
* `dechirped_residual` and `score_ins23_ins24` are Phase 4.
* The canceller (`tss.l3_signal.leakage_suppress`) is WP5.
