# leakage_budget

WP3 Phase 5: leakage budget; AM noise x RF canceller x isolation x power x front end x start-up

Status: **implemented**.  Work packages: WP3.  Claims: C4.  Requirements: INS-10, INS-11, INS-14, INS-22, INS-23, INS-24.

Leakage budget, A1 unless stated.  INS-24 margin is the worst step of the worst of 3 seeds against the weakest admitted echo; it must be at least 10 dB.  'AM needed' is the AM density at which the worst step just meets it, from the characterisation.

Model checked against direct simulation at 16 grid points: largest error 0.05 dB.  INS-10 prediction checked against the smoke campaign at 3 points: largest error 1.48 dB.

## Canceller against AM (MoM bus isolation, 1 W)

| canceller | nominal depth (dB) | AM needed (dBc/Hz) | INS-24 at -133 (dB) | at -140 | at -150 | INS-10 at -140 (dB) |
|---|---|---|---|---|---|---|
| none |  | -159.7 | -16.7 | -9.7 | 0.3 | 16.5 |
| ideal20 | 20.0 | -139.9 | 3.1 | 10.1 | 20.1 | 19.2 |
| ideal26_ins31 | 26.0 | -133.9 | 9.1 | 16.1 | 26.1 | 19.2 |
| ideal30 | 30.0 | -129.9 | 13.1 | 20.1 | 30.1 | 19.2 |
| ideal40 | 40.0 | -119.9 | 23.1 | 30.1 | 40.1 | 19.2 |
| rf_good | 33.6 | -133.8 | 9.2 | 16.2 | 26.2 | 18.7 |
| rf_moderate | 26.1 | -140.4 | 2.6 | 9.6 | 19.6 | 18.7 |
| rf_poor | 22.0 | -144.0 | -1.0 | 6.0 | 16.0 | 18.7 |
| rf_good_bias_only | 33.5 | -126.4 | 16.6 | 17.7 | 18.0 | 18.7 |
| rf_good_uncorrelated | 33.6 | none suffices | -0.3 | -0.2 | -0.2 | 18.7 |
| rf_good_stepped | 33.6 | none suffices | -2.8 | -0.4 | 0.4 | 18.7 |
| rf_good_stepped_uncorrelated | 33.6 | none suffices | -33.3 | -33.3 | -33.3 | 18.7 |
| rf_one_setting | 33.5 | -138.5 | 4.5 | 11.5 | 21.5 | 18.7 |
| jump_0.001dB | 78.8 | -112.8 | 12.8 | 12.8 | 12.8 | 18.7 |
| jump_0.003dB | 69.2 | none suffices | 3.5 | 3.5 | 3.5 | 18.7 |
| jump_0.01dB | 58.8 | none suffices | -6.9 | -6.9 | -6.9 | 18.7 |
| jump_0.03dB | 49.2 | none suffices | -16.4 | -16.4 | -16.4 | 18.7 |
| rf_good_active | 33.6 | none suffices | 1.1 | 1.6 | 1.7 | 18.4 |

## Configuration A2 (MoM bus isolation, 1 W)

| canceller | A2: AM needed (dBc/Hz) | A2: INS-24 at -140 (dB) |
|---|---|---|
| none | -160.1 | -10.1 |
| ideal26_ins31 | -134.3 | 15.7 |
| rf_good | -123.0 | 27.0 |
| rf_moderate | -132.1 | 17.9 |
| rf_one_setting | -148.1 | 1.9 |
| rf_good_stepped | none suffices | 5.7 |
| rf_good_active | none suffices | 1.7 |

## Isolation (1 W)

| isolation | AM needed, no canceller | AM needed, rf_good | AM needed, rf_one_setting |
|---|---|---|---|
| mom_bus | -159.7 | -133.8 | -138.5 |
| flat20 | -176.8 | -150.4 | -144.4 |
| flat30 | -166.5 | -140.4 | -134.4 |
| flat40 | -156.5 | -130.4 | -124.4 |
| flat50 | -146.5 | -120.4 | -114.4 |
| flat60 | -136.5 | -110.4 | -104.4 |

## Transmit power (MoM bus)

| power (W) | INS-24 at -133, no canceller (dB) | AM needed, no canceller | AM needed, rf_good | INS-10 min at -133, no canceller (dB) | INS-10 min, rf_good at its AM limit (dB) | echoes below 10 dB |
|---|---|---|---|---|---|---|
| 1 | -16.7 | -159.7 | -133.8 | 11.8 | 18.7 | 0 of 44 |
| 0.5 | -13.7 | -156.7 | -130.8 | 11.1 | 15.7 | 0 of 44 |
| 0.25 | -10.7 | -153.7 | -127.8 | 10.0 | 12.7 | 0 of 44 |
| 0.1 | -6.7 | -149.7 | -123.8 | 7.7 | 8.7 | 2 of 44 |

## Front end (MoM bus, 1 W)

| canceller | AM (dBc/Hz) | P1dB +1 dBm: margin (dB), image (dBc) | P1dB +7 dBm: margin (dB), image (dBc) | P1dB +11 dBm: margin (dB), image (dBc) |
|---|---|---|---|---|
| none | -140 | 7.6, -41.2 | 13.6, -65.1 | 17.6, -81.1 |
| rf_good | -140 | 38.6, -165.1 | 44.6, -189.1 | 48.6, -205.4 |
| rf_one_setting | -140 | 9.0, -46.6 | 15.0, -70.6 | 19.0, -86.6 |

## Start-up (MoM bus, -133 dBc/Hz)

| start-up | canceller | INS-11 (dB) | residual peak over clip (dB) | INS-24 (dB) | INS-24 step 1 (dB) |
|---|---|---|---|---|---|
| primed | none | 60.0 | -40.2 | -16.7 | -16.7 |
| cold | none | 29.9 | 30.1 | -43.9 | -43.9 |
| settling_step | none | 60.0 | -40.2 | -16.9 | -16.9 |
| primed | rf_good | 60.0 | -74.1 | 24.0 | 24.5 |
| cold | rf_good | 60.0 | -11.3 | -2.5 | -2.5 |
| settling_step | rf_good | 60.0 | -74.2 | 23.9 | 24.3 |

---

Run 2026-10-29T02:40:00+00:00.  Provenance is in `manifest.json` next to this file.  Rerunning `tss study leakage_budget` from the same commit and the same configuration reproduces it.
