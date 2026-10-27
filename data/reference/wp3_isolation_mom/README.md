# WP3 isolation anchor: moment method (G1(b))

Produced 2026-10-27 by `scripts/matlab/wp3_isolation_mom.m` (MATLAB R2026a, Antenna Toolbox).
Two 5 m dipoles (radius 5 mm, the `baseline_6u` antenna), feeds 0.366 m apart, free space.
Each antenna's matching coil is a series loss |X|/Q with Q = 150, as `tss.l2_instrument.antenna.ShortDipole`.
The bus is not in the model.

| Column | Meaning |
|---|---|
| `iso_<geom>_db` | 1 / maximum available gain: both ports conjugate-matched simultaneously, the other antenna present. The least isolation the pair allows. |
| `iso_<geom>_self_db` | 1 / transducer gain with each port matched to its own self-impedance. |
| `*_lossless_db` | The same without the coil loss (not physical for this antenna; kept for the record). |
| `eta_tx_db` | Transmit radiation efficiency, R_rad / (R_rad + |X|/Q), from the MoM self-impedance. |

Geometries: `parallel` (both along x: the worst case) and `crossed_offset` (transmit along x,
receive along y, receive feed displaced 3 cm and 8 cm from the axis: the anchor).

What it shows: below a few MHz the coil loss dominates both antennas and the pair behaves as two
coupled resonators, so the isolation is set by the coupling coefficient k = |Z12| / |X| times the
loaded Q, not by the antenna efficiency, and is flat with frequency. Crossed-offset: k = 4.7e-4,
29.0 dB from 0.33 to 3 MHz, 35 dB at 11 MHz, 49 dB at 20 MHz. Parallel: k = 0.20, 0.3 dB (simultaneous
match) or 23.5 dB (self-matched). On this model INS-23's 40 dB needs k below about 1.3e-4 at Q = 150.
The MoM self-impedance agrees with the Balanis short-dipole form to 3 % in reactance (0.33 to 20 MHz).

Loaded by `tss.l2_instrument.leakage.mom_anchor_model()`.
