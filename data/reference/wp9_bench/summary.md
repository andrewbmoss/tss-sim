# wp9_bench

WP3 Phase 6: WP9 conducted-bench predictions and pass limits (INS-29 behind emulated canceller paths)

Status: **implemented**.  Work packages: WP3, WP9.  Claims: C4.  Requirements: INS-11, INS-24, INS-29, INS-30, INS-31.

WP9 conducted bench: predictions and pass limits (WP3 Phase 6).

Paths from the PA output at 1 W; delay 1.98 m RG-58 unless stated; worst of 3 seeds; bench floor assumed -165 dBm/Hz.

## INS-29 behind each bench path (worst step)

| sweep | path | AM (dBc/Hz) | INS-29 (dB re k T_sys) | at (MHz) | least floor rise (dB) |
|---|---|---|---|---|---|
| A1 | E0 (39 dB) | -140 | +17.5 | 0.338 | 16.27 |
| A1 | E0 (39 dB) | -150 | +7.5 | 0.338 | 7.11 |
| A1 | E1 (61 dB) | -140 | -4.5 | 0.338 | 1.01 |
| A1 | E1 (61 dB) | -150 | -14.5 | 0.338 | 0.11 |
| A1 | E2 (69 dB) | -140 | -12.5 | 0.338 | 0.18 |
| A1 | E2 (69 dB) | -150 | -22.5 | 0.338 | 0.02 |
| A2 | E0 (39 dB) | -140 | +17.9 | 0.355 | 16.10 |
| A2 | E0 (39 dB) | -150 | +7.9 | 0.355 | 6.97 |
| A2 | E1 (61 dB) | -140 | -4.1 | 0.355 | 0.97 |
| A2 | E1 (61 dB) | -150 | -14.1 | 0.355 | 0.11 |
| A2 | E2 (69 dB) | -140 | -12.1 | 0.355 | 0.17 |
| A2 | E2 (69 dB) | -150 | -22.1 | 0.355 | 0.02 |

## Checks the bench makes on the model

- Delay sensitivity, E0, 0.66 m to 19.8 m: A1 0.50 dB at -140, 3.31 dB at -150 (AM binds, not phase noise).
- Linear scaling E0 to E1 (22 dB): largest departure 0.000 dB (A1).
- AM at which E1 just meets INS-29: A1 -135.5 dBc/Hz at 0.338 MHz; A2 -135.9.
- Floor rise at the INS-29 limit, least over the band: 2.42 dB (k T_sys down to -166.3 dBm/Hz).

## Front end

| path | leakage (dBm) | level-7: margin, image | level-17: margin, image |
|---|---|---|---|
| E0 | -9.0 | 10.0 dB, -50.7 dBc | 20.0 dB, -90.7 dBc |
| E1 | -31.0 | 32.0 dB, -138.7 dBc | 42.0 dB, -178.7 dBc |
| E2 | -39.0 | 40.0 dB, -170.7 dBc | 50.0 dB, -212.2 dBc |
| INS-23 35 dB, no canceller | -5.0 | 6.0 dB, -34.8 dBc | 16.0 dB, -74.7 dBc |

## Stray coupling and echo levels (WP9-1)

- Weakest admitted echo on the bench -137.5 dBm; strong echo -87.5 dBm.
- Transmit-chain stray: >= 79 dB from the PA output (10 dB under E2); its noise needs >= 65 dB.
- A synthesiser-channel echo (AD9959, > 65 dB channel isolation): weakest usable echo -64 / -86 / -94 dBm behind E0 / E1 / E2, i.e. SNR 84 / 62 / 54 dB on the bench.  Weak echoes with leakage present need the DRFM.

## AM measurement (WP9-5) and canceller (WP9-8)

- AM floor, one chain, level-7 mixer at -9 dBm: -159.0 dBc/Hz; two-chain cross-spectrum, 1000 averages: -174.0; level-17 mixer: -169.0.
- INS-31 weight step limit -78.8 dB of the leakage; added noise <= -170 dBm/Hz, which bounds the weight's own noise under drive at -165 dBc/Hz of the leakage at 35 dB isolation.

---

Run 2026-10-30T00:26:52+00:00.  Provenance is in `manifest.json` next to this file.  Rerunning `tss study wp9_bench` from the same commit and the same configuration reproduces it.
