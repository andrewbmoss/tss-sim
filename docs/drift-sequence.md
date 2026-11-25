# The drift-mode sequence path (WP5 Phase 2)

2026-11-25.  Plan: project note `claude/wp5-plan-draft.md`, Phase 2.  Decisions
(Andrew, 2026-11-25 03:56): `phase_delay_s` on the echo set; the orbit and the
L0 profile now; the specular fading gain, oscillator continuity and the layer
echo in this phase; a 150 kHz chirp at 45 Hz.

## What it is

```
Pipeline.run_sequence(configuration, sweep=, t_s=, platforms=, alt_m=, epoch=, truth_at=)
    L0 once (at truth_at, the facet's sub-point)
    L1 per sweep, at the sweep's centre t_s[n] + T/2          -> EchoSet with phase_delay_s
    L2 receive_sequence over the whole sequence               -> ReceivedSequence
```

* **Sweep law** (`tss.sweep_law.drift_mode`, `drift_mode_times`): one linear
  up-chirp over 150 kHz centred on `f_c`, repeated at 45 Hz with no dead time
  (71 Hz an option).  6.75 MHz/s, about 1 km cells, 2101 samples per record at
  the INS-06 rate.
* **Interfaces** (`tss.interfaces`): `EchoSet.phase_delay_s` (optional; the
  two-way phase path over c), `PlatformState.t_s` (optional; seconds into the
  sequence), `ReceivedSequence`, `DriftProduct` (defined; filled in Phase 3),
  protocol `L3DriftProcessing`.
* **L1 `facet`** (`tss.l1_propagation.facet`): the plate model per sweep (D4),
  rows flagged `facet:<id>`, plus the layer echo (`layer:O`) when asked.
  General entry `plate_echo.echo_track` beside the unedited `echo_free_space`
  and `echo_stratified`, which it reproduces exactly on a straight level track.
* **L2** (`BaselineInstrument.receive_sequence`): per sweep, with the Doppler
  applied inside the sweep about its centre, the first group delay blanked,
  one impairment draw over the whole sequence (oscillator continuity), and the
  specular fading gain keyed to each facet's id.

## Conventions

* **Phase.**  L2's phase is `+2 pi f tau`: a closing echo's phase falls from
  sweep to sweep and its Doppler `-2 pi f_d t` turns it the same way inside a
  sweep.  The plate field is `exp(-j 2 k0 P)`, so `2 pi f tau_p = -arg(E)`.
* **Phase path (D8 condition 1).**  The beat and envelope follow the group delay;
  the carrier adds `2 pi f_k (tau_p - tau_g)`, the stationary-phase result for
  `d(f tau_p)/df = tau_g` over the sub-band.  Rows without the field are as
  before.  Condition 2 (convert at the gate's effective frequency) is L3's,
  Phase 3.
* **Evaluation at the sweep centre.**  L1 is evaluated at `t_s[n] + T/2` and L2
  holds the Doppler fixed about that time, so the phase is exact at the centre
  and the error of a fixed Doppler is second order.  With L1 at the sweep start
  the error near the first nulls of the envelope reached 0.02 rad.
* **Blanking.**  In a repeated chirp the first `tau_g` of each record holds the
  previous sweep's echo at another beat.  It is not modelled; L2 leaves it empty
  (`flags["first_delay_blanked"]`) and L3 must gate it out.  At 500 to 1000 km
  that is 15 to 30 % of each 22 ms record, a cost to the per-sweep SNR that
  Phase 3 must count.
* **Doppler in the drift mode is always on**, whatever
  `transmit_impairments.doppler` says (the ionogram campaigns keep it off).

## Exit criterion (WP5 plan, P2)

Noise-free, `quantise=False`, 5 MHz, reference plasma, a 2 km facet with
`mu_f` 0.5 drifting at 25 m/s, 91 sweeps across overflight
(`tests/unit/test_run_sequence.py`):

| Check | Criterion | Measured |
|---|---|---|
| L2 phase against `echo_stratified` called directly, whole window to the first nulls | 0.01 rad | 1.3e-9 rad (3.3e-10 over the main lobe) |
| Echo-set Doppler at overflight against `doppler_shift_hz(25 m/s, f, mu_sat)` | 1e-3 | within 1e-3 (central difference) |
| Doppler from L2 phase steps about overflight | 1e-3 | within 1e-3 |
| Same chain on the orbit in the L0 profile, against the tier itself at the record centre (chain self-consistency, not an independent check: no spherical reference exists) | 0.01 rad | 6.6e-9 rad, all 61 sweeps |
| Smoke campaign and `wp4_campaign` unchanged | byte-identical | identical except run time and the requirements path (see the commit) |

## Orbit and the L0 profile (new, unverified)

* **Earth flattening** about `R = R_E + ref_alt_m`: `x = R psi`,
  `d = R ln(R/r)`, `mu'^2 = mu^2 (r/R)^2`.  Conformal, and exact for rays and
  phase paths (Bouguer's law becomes Snell's).  The group path is
  `d(f P)/df` by a central difference in frequency, because the plate model's
  group integral assumes `mu_g = 1/mu`.
* **Checks** (`tests/unit/test_l1_facet_orbit.py`): free-space range equals the
  spherical chord to 2e-6; overflight Doppler is the drift Doppler; the
  tabulated L0 profile reproduces the exponential to under 0.01 cycle and its
  second difference across sweeps to 1e-3 rad.
* **Finding for Phase 3.**  Near overflight the range's second derivative is
  `v^2 (r_f / r_p) / D`, not the flat `v^2 / D`: 4.3 % less for a facet
  300 km below a 600 km orbit (`r_f / r_p = 0.957`; the tier matches the
  formula to 0.5 % on a point-like 50 m plate).  In a plasma the factor is
  `P' / P'_flat` with `P'_flat = int (r_p / r)^2 dz / mu`.  A 2 km plate,
  part mirror, has a smaller slope than a point; that is the envelope physics
  route one already handles.
* **Corrected 2026-11-25 (Phase 3, measured).**  This note said route one
  would read the drift about 4 % wrong without the effective speed.  It does
  not: `estimate_route_one` takes the zero crossing with the ephemeris slope
  `S`, so its drift `v^2 (t_zc - t_pk) / (mu_sat P')` is the Doppler at the
  power peak converted to a velocity, and `v^2` cancels to first order.  On the
  orbit at signal level the drift with and without `v_eff` differs by 0.0004
  m/s (`tests/unit/test_l3_drift.py`).  What the curvature moves is the slope
  ratio, the compact-facet test of `docs/drift-mode.md` section 3, by
  `(v_eff / v)^2` (0.956 here: 0.965 with `v_eff` against 0.922 without).  The
  drift tier uses `v_eff` for that reason.
* **Limits.**  The plate is two-dimensional: the facet is taken in the track
  plane and the cross-track offset is ignored (largest seen in
  `FacetTier.max_cross_track_m`; under 50 m over +/-1 s in the tests).  The
  medium's drift is mapped with the scale at the platform, along-track drift
  with the scale at the facet.

## Fading and the oscillator

* **Specular gain only** (`BaselineInstrument(fading=FadingModel(...))`, with
  `k_factor=inf`): one stream per facet id, sorted, the drawn `tau_i` recorded
  in each record's `flags["fading"]`.  A finite K is refused: a diffuse halo
  wider than the facet does not share its envelope, so it cannot be a gain on
  it (Phase 3 needs its own scatterer field).
* **Oscillator continuity**: one `ImpairmentRealization` over the whole
  sequence, its carrier the sawtooth `f0 + K mod(t, T)`, cut per record by
  `RealizationWindow`; `flags["oscillator_continuous"]`.  The leakage, when on,
  reads the same window; that path runs but has not been timed or scored over a
  sequence.

## Not done here (deviations from the plan)

* `sheet_echo` is not wrapped (plan D4 named both); plates only.
* No drift tier, no Cramer-Rao harness, no diffuse halo: Phase 3.
* No study or campaign file drives `run_sequence` yet; Phase 3's
  `wp5_route_one` will.
