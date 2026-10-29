# leakage_residual

WP3 Phase 4: residual leakage in the echo band and the INS-24 margin, per step, A1 and A2

Status: **implemented**.  Work packages: WP3.  Claims: C4.  Requirements: INS-10, INS-11, INS-23, INS-24.

Residual leakage in each step's echo band, scored against the weakest echo INS-10 admits (INS-24, 10 dB), with no canceller.

* Default, A1: margin -17.6 dB at 0.338 MHz; 105 of 159 steps fail.  R_max reading -32.2 dB.  Pedestal up to 27.6 dB above the noise in a cell.
* Behind an ideal 30 dB canceller: margin 12.3 dB, 0 steps failing.  Canceller needed: 27.6 dB.
* Without AM noise: margin 24.0 dB.  AM sets the pedestal; the AM density that would meet INS-24 at every step is -160.6 dBc/Hz.
* On a flat 40 dB isolation: margin -14.3 dB.
* INS-11 verified at 60 dB; least rejection margin on the noise rule 23.4 dB.  INS-23: 36.6 dB at 1.39 MHz.

## A1, default (allocations at threshold, primed)

| f (MHz) | iso (dB) | pedestal (dBW) | over noise (dB) | weakest echo (dBW) | margin (dB) | margin, R_max (dB) | margin, canceller (dB) | AM needed (dBc/Hz) |
|---|---|---|---|---|---|---|---|---|
| 0.338 | 36.6 | -151.2 | 27.6 | -168.8 | -17.6 | -32.2 | 12.3 | -160.6 |
| 0.500 | 36.6 | -151.9 | 24.8 | -166.7 | -14.8 | -24.7 | 15.1 | -157.8 |
| 0.989 | 36.6 | -152.6 | 19.7 | -162.4 | -9.7 | -12.2 | 20.1 | -152.7 |
| 1.263 | 36.6 | -151.6 | 19.1 | -160.7 | -9.1 | -9.0 | 20.7 | -152.1 |
| 1.959 | 36.6 | -152.0 | 15.6 | -157.6 | -5.6 | -0.9 | 24.2 | -148.6 |
| 3.039 | 36.7 | -151.9 | 12.5 | -154.4 | -2.5 | 6.5 | 27.3 | -145.5 |
| 4.988 | 37.0 | -152.0 | 9.0 | -151.0 | 1.0 | 14.8 | 30.8 | -142.0 |
| 9.938 | 40.1 | -155.4 | 2.4 | -147.7 | 7.6 | 26.6 | 37.6 | -135.4 |
| 19.957 | 52.7 | -167.7 | -9.2 | -148.4 | 19.2 | 39.4 | 49.2 | -123.8 |

## A2, default (allocations at threshold, primed)

| f (MHz) | iso (dB) | pedestal (dBW) | over noise (dB) | weakest echo (dBW) | margin (dB) | margin, R_max (dB) | margin, canceller (dB) | AM needed (dBc/Hz) |
|---|---|---|---|---|---|---|---|---|
| 0.338 | 36.6 | -152.4 | 26.4 | -168.8 | -16.4 | -31.0 | 13.4 | -159.4 |
| 0.500 | 36.6 | -151.0 | 25.7 | -166.7 | -15.7 | -25.7 | 14.1 | -158.7 |
| 0.989 | 36.6 | -151.9 | 20.4 | -162.4 | -10.4 | -12.8 | 19.4 | -153.4 |
| 1.263 | 36.6 | -152.2 | 18.5 | -160.7 | -8.5 | -8.4 | 21.3 | -151.5 |
| 1.959 | 36.6 | -151.8 | 15.8 | -157.6 | -5.8 | -1.2 | 24.0 | -148.8 |
| 3.039 | 36.7 | -152.1 | 12.4 | -154.4 | -2.4 | 6.7 | 27.4 | -145.4 |
| 4.951 | 37.0 | -152.3 | 8.8 | -151.1 | 1.2 | 15.0 | 31.1 | -141.8 |
| 9.802 | 40.0 | -155.3 | 2.5 | -147.8 | 7.5 | 26.4 | 37.5 | -135.5 |
| 19.940 | 52.6 | nan | nan | -148.4 | nan | nan | nan | nan |

---

Run 2026-10-29T01:15:58+00:00.  Provenance is in `manifest.json` next to this file.  Rerunning `tss study leakage_residual` from the same commit and the same configuration reproduces it.
