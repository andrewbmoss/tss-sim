# WP3 isolation anchor: moment method (G1(b))

Produced 2026-10-27 with MATLAB R2026a, Antenna Toolbox.  Two 5 m dipoles, radius 5 mm (the
`baseline_6u` antenna).  Each antenna's matching coil is a series loss |X|/Q with Q = 150, as
`tss.l2_instrument.antenna.ShortDipole`.  The MoM self-impedance agrees with the Balanis
short-dipole form to 3 % in reactance from 0.33 to 20 MHz.

## Files

| File | Script | What |
|---|---|---|
| `isolation_mom.csv` | `wp3_isolation_mom.m` | Free space, feeds 0.366 m apart: `parallel` and `crossed_offset` (rx feed displaced 3, 8 cm). Rerun at a 5 cm dipole mesh; the first run's default mesh was 0.6 dB off. |
| `isolation_mom_bus.csv` | `wp3_isolation_mom_bus.m` | On a closed conducting 6U box (0.100 x 0.2263 x 0.366 m, 1.5 cm mesh), full MoM. Four mounts, with the coupling coefficient `k_<mount>`. |
| `isolation_mom_bus_sensitivity.csv` | same | Each imperfection alone on the symmetric crossed mounts: deployment angle error 0 to 3 deg, rx feed displaced along its own axis 0 to 10 cm. |

Columns `iso_<mount>_db` are 1/MAG (both ports conjugate-matched together, the other antenna present:
the least isolation the pair allows); `_self_db` matches each port to its own self-impedance;
`_lossless_db` omits the coil loss (not physical here); `eta_tx_db` is the transmit efficiency.

Bus mounts (bus from z = 0 to 0.366 m; x the 10 cm side, y the 22.6 cm side; feeds 5 cm off the face):
`crossed_ends_offset` (tx along x below, rx along y above, rx feed at x = 3, y = 8 cm),
`crossed_ends_centred` (rx feed on the axis), `parallel_ends`, `crossed_colocated` (both below,
tx along x at z = -5 cm, rx along y at z = -8 cm: a turnstile pair).

## What it shows

1. **The isolation is the coupled-resonator law.**  Below a few MHz the coil loss dominates both
   antennas, and the isolation is `-20 log10(k Q / 2)` (exactly `kappa + sqrt(kappa^2 - 1)`,
   `kappa = 1 + 2/(kQ)^2`; `tss.l2_instrument.leakage.coupled_resonator_isolation_db`), flat with
   frequency.  It does not grow as efficiency falls.  INS-23's 40 dB needs k < 1.3e-4 at Q = 150.
2. **The bus helps a little.**  Crossed with the (3, 8) cm offset: k = 3.2e-4 in free space at the same
   feed spacing, 2.0e-4 on the bus, so 36.6 dB from 0.33 to 4 MHz, 41 dB at 11 MHz, 53 dB at 20 MHz.
   Parallel: k = 0.16, 0.4 dB (21.7 dB self-matched).
3. **A symmetric crossed pair is limited by its imperfections, and angle dominates.**  The symmetric
   mounts sit at the solver's floor (k ~ 5e-6, 67 to 69 dB), which is not a physical number.  Each
   imperfection alone, isolation at 0.33 MHz:

   | Deployment angle error | 0.1 deg | 0.3 deg | 1 deg | 3 deg |
   |---|---|---|---|---|
   | Ends (feeds 0.47 m apart) | 39.7 dB | 30.7 dB | 20.6 dB | 11.5 dB |
   | Co-located (3 cm apart) | 36.1 dB | 26.8 dB | 13.9 dB | 8.2 dB |

   A feed displaced along its own dipole's axis costs little at the ends (10 cm: 61.9 dB) and more
   co-located (3 cm: 46.4 dB).  A displacement across both axes, as in `crossed_ends_offset`, is what
   costs (36.6 dB).
4. **Reading.**  40 dB at the bottom of the band needs the two dipoles orthogonal to about 0.1 deg
   (4 mm at the tip of a 2.5 m arm) and their feeds centred to a few cm.  Deployable dipoles are not
   usually that straight, so on this model INS-23 is not met by the antenna geometry alone at the
   bottom of the band, and the margin has to come from somewhere else (a lower-Q match at an
   efficiency cost, a hybrid or cancellation stage, or a relaxed INS-23 with a tighter INS-22).
   That is WP3 Phase 5's trade.

Mesh: 5 cm on the end mounts (coupling converged to 2 % against 10 cm); 3 cm co-located, which
converges only to about 3 dB.  Not modelled: the plasma sheath, deployment-boom structure, solar
panels, the harness, and the matching networks' own coupling.

Loaded by `tss.l2_instrument.leakage.mom_anchor_model()` (default: `bus_crossed_ends_offset`).
