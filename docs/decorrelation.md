# Echo decorrelation time and the fading model (WP5 Phase 1)

2026-11-24.  Gate G1 decided by Andrew the same evening.  Module:
`src/tss/l2_instrument/fading.py` (standalone; Phase 2 wires it in).  Project
note: `claude/wp5-p1-decorrelation.md`.

## Definition

`tau_c` is the **field** autocorrelation 1/e time, `|rho(tau_c)| = 1/e`.  That
is the definition of the AR(1) gain `route_one` has used since October
(`r = exp(-dt / tau)`).  Every number below is in it.  A Gaussian Doppler
spectrum of standard deviation `sigma_f` has `tau_c = 1 / (sqrt(2) pi sigma_f)`.

## What the literature gives

No measurement of echo decorrelation from a topside sounder in orbit was
found.  Every number comes from a stationary receiver on the ground, mostly
looking at the bottomside, and bounds only the intrinsic change of the
irregularities.

| Source | Geometry | tau_c (field 1/e) |
|---|---|---|
| ITU-R F.1487 (2000), Doppler spread table | ground, oblique | mid-latitude quiet / moderate / disturbed 4.5 / 0.9 / 0.45 s; low-latitude 0.9 / 0.30 / 0.045 s; high-latitude 0.9 / 0.045 / 0.015 s |
| Lowell DPS (Digisonde) practice | ground, vertical | coherent integration 0.16 to 2.56 s routine; coherence for tens of seconds "not rare" |
| Ganguly et al. (2001) | review | tens of ms to several s; bottomside often >= 100 ms |
| Briggs and Phillips (1950), via NBS TN 133 | ground, vertical | angular spread about 5 degrees (unverified; secondary source) |

The topside is not known to be the more benign.  In 99 coincident Alouette
and ground soundings irregularities appeared on the topside alone in 46, on
both in 23, on neither in 30 and on the bottomside alone in none (Dyson 1967a,
via Jackson 1988, NSSDC 88-10).  Topside irregularities are hundreds of metres
to many kilometres across the field and normally a few percent in amplitude
(Calvert and Warnock); auroral and polar fine structure is about 1.3 km thick
at about 25 % (Dyson 1969).  Transionospheric fading seen by ISIS-II and e-POP
RRI is mostly Faraday fading and is not a source for `tau_c`.

## G1 decisions (2026-11-24)

1. Start from the ground-based literature, held at **low confidence and high
   variance**.  Each echo's intrinsic time is drawn log-normally about the
   carried value with `sigma(ln tau_i) = ln 3`
   (`LOW_CONFIDENCE_LOG_SIGMA`).
2. **Tie the diffuse decorrelation to the facet.**
3. Carry the intrinsic time at **0.1, 0.5 and 1.0 s** (`CARRY_TAU_I_S`).
4. The **0.1 s fading case of the route-one Monte Carlo is an extended
   reflector** and is no longer run there (`route_one.EXTENDED_FADING_CASES`).

## The tie to the facet

A randomly phased patch of along-track extent `L` at apparent range `P'` has a
Doppler spectrum uniform over `B = 2 v_s L / (lambda P')`, in route one's
convention, so its field correlation is `sinc(pi B tau)`.  It falls to 1/e at

    tau_d = x_e lambda P' / (2 pi v_s L),    x_e = 2.199

and has its first null at `lambda P' / (2 v_s L)`, as the envelope of a
coherent facet of size `L` does.  Against the half-power width `W` of a facet
of size `D`:

    tau_d / W = 0.79 D / L        (any frequency and range)

A diffuse halo wider than about 0.8 D decorrelates faster than the facet's
envelope.  Patch extent for a given time (7.6 km/s):

| f, P' | 0.5 s | 0.1 s | one sweep at 45 Hz |
|---|---|---|---|
| 3 MHz, 500 km | 4.6 km | 23 km | 104 km |
| 5 MHz, 500 km | 2.8 km | 14 km | 63 km |
| 5 MHz, 1000 km | 5.5 km | 28 km | 125 km |
| 10 MHz, 500 km | 1.4 km | 6.9 km | 31 km |

At 5 MHz and 500 km the D_length scenes switch from compact to extended
between 2 and 4 km, so the 0.1 s case (14 km) is extended and the 0.5 s case
(2.8 km) is compact.

The halo also gives a coherence bandwidth: the two-way excess delay
`x^2 / (c P')` of a uniform patch has standard deviation
`L^2 / (6 sqrt 5 c P')`, and `B_c = 1 / (2 pi sigma_tau)`.  It is an order of
magnitude only.

## The model

    g = sqrt(K / (K+1)) s + sqrt(1 / (K+1)) d

* `s`, specular: 1 when steady, otherwise a unit-power complex Gaussian
  process with the per-echo intrinsic time `tau_i`.
* `d`, diffuse: a unit-power complex Gaussian process with `tau_d`, set from a
  halo of `CARRY_HALO_RATIOS` = 1, 2 or 4 facets (`FadingModel.from_halo`).
  White when shorter than the sweep step.
* `kind`: `"exp"` (AR(1), any increasing grid) or `"gauss"` (Gaussian kernel,
  uniform grid).
* Across frequency: `exp(-|df| / B_c)` per part, by an AR(1) across the
  sorted frequencies; no `B_c` is one draw for every frequency.
* K carried: inf, 10, 3, 1, 0.  `K = inf` with `tau_i = 0.5 s` and no spread is
  the route-one budget's model; `K = 0` is Rayleigh.
* Streams: one per echo (`echo_streams`); inside an echo, the intrinsic-time
  draw, the specular part and the diffuse part each take a spawned child.

Tests: `tests/unit/test_fading.py` (ACF within 0.05 of the closed form at lags
1 to 50, unit power, Rician moments, stream independence, frequency
correlation, the halo relations).

## Budget-point facets (same evening)

The six budget points cannot all hold a facet in the c3_verify reference
profile (f_N,sat 1.5 MHz, exponential, H = 200 km): a facet below the
reflection level is seen at no more than 527 km at 3 MHz, 750 km at 5 MHz and
1034 km at 10 MHz.  `wp5_campaign` adds the two that can (`E_budget`: 5 and
10 MHz at 500 km, `mu_f` 0.554 and 0.870); the other four are listed in
`BUDGET_POINTS_UNREACHABLE`.  The analytic budget in TV Table `c3r1budget`
treats `P'` and `mu_sat` as independent, so it pairs four of its ranges with
a plasma frequency at the spacecraft that this profile cannot deliver.

## Not yet done

* Wiring into the drift-mode sequence path (Phase 2), and the `tau_c_s` data
  product's estimator.
* TV Table `c3nodeproducts` quotes the `tau_c_s` product's accuracy at
  0.1 s, and the C4 null-scan screen ran 0.1 to 2 s.  Both are left as they
  are; they measure or tolerate a time, not a compact-facet budget.
