# leakage_receive_chain

WP3 Phase 3: leakage through the front end and the leakage high-pass, per step, A1 and A2

Status: **implemented**.  Work packages: WP3.  Claims: C4.  Requirements: INS-11, INS-23, INS-24.

Leakage on the moment-method isolation of the dipoles on the 6U bus, front end 1.0 dBm P1dB (Rapp), high-pass order 4 at 0.5 of each step's 100 km beat.

* INS-11, noise-referenced rule (G5b): met at every step of A1, least margin 39 dB; the rule asks at most 60.7 dB.  Rev G's rule asks at most 96.8 dB.
* Settling, cold: longest 25.07 ms (0.94 of a dwell); first step 25.07 ms; steps settling longer than the 7.0 ms echo transition: 1.
* Front end: least margin below P1dB 7.6 dB; steps under WP9-2's 10 dB: 83; worst echo desensitisation -0.113 dB; worst image -41.2 dBc.
* Echo loss at 100 km through the high-pass: at most 0.017 dB.

## A1, MoM bus isolation, filter cold

| f (MHz) | corner (Hz) | iso (dB) | L in (dBm) | L/N_fs (dB) | need, noise rule (dB) | need, Rev G rule (dB) | achieved (dB) | settle (ms) | settle / dwell | FE margin (dB) | desens (dB) | image (dBc) | loss at 100 km (dB) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.338 | 205 | 36.6 | -6.7 | 110.7 | 60.7 | 96.8 | 99 | 25.07 | 0.936 | 7.6 | -0.113 | -41.2 | 0.017 |
| 0.500 | 303 | 36.6 | -6.7 | 108.6 | 58.6 | 90.0 | 223 | 0.00 | 0.000 | 7.6 | -0.113 | -41.2 | 0.017 |
| 0.989 | 601 | 36.6 | -6.7 | 104.2 | 54.2 | 78.2 | 219 | 1.17 | 0.044 | 7.6 | -0.113 | -41.2 | 0.017 |
| 1.263 | 767 | 36.6 | -6.7 | 102.6 | 52.6 | 73.9 | 216 | 0.94 | 0.035 | 7.6 | -0.113 | -41.2 | 0.017 |
| 1.959 | 1189 | 36.6 | -6.7 | 99.4 | 49.4 | 66.3 | 211 | 0.00 | 0.000 | 7.6 | -0.113 | -41.2 | 0.017 |
| 3.039 | 1845 | 36.7 | -6.7 | 96.3 | 46.3 | 58.7 | 206 | 0.00 | 0.000 | 7.7 | -0.110 | -41.4 | 0.017 |
| 4.988 | 1867 | 37.0 | -7.1 | 92.5 | 42.5 | 50.2 | 204 | 0.00 | 0.000 | 8.0 | -0.095 | -42.8 | 0.017 |
| 9.938 | 1867 | 40.1 | -10.2 | 86.1 | 36.1 | 38.6 | 200 | 0.00 | 0.000 | 11.1 | -0.023 | -55.2 | 0.017 |
| 19.957 | 1079 | 52.7 | -22.7 | 74.3 | 24.3 | 25.7 | 192 | 0.00 | 0.000 | 23.7 | -0.000 | -105.3 | 0.017 |

## A2, MoM bus isolation, filter cold

| f (MHz) | corner (Hz) | iso (dB) | L in (dBm) | L/N_fs (dB) | need, noise rule (dB) | need, Rev G rule (dB) | achieved (dB) | settle (ms) | settle / dwell | FE margin (dB) | desens (dB) | image (dBc) | loss at 100 km (dB) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0.338 | 205 | 36.6 | -6.7 | 102.7 | 52.7 | 96.8 | 99 | 22.60 | 0.843 | 7.6 | -0.113 | -41.2 | 0.017 |
| 0.500 | 303 | 36.6 | -6.7 | 100.6 | 50.6 | 90.0 | 223 | 0.00 | 0.000 | 7.6 | -0.113 | -41.2 | 0.017 |
| 0.989 | 601 | 36.6 | -6.7 | 96.2 | 46.2 | 78.2 | 219 | 0.08 | 0.003 | 7.6 | -0.113 | -41.2 | 0.017 |
| 1.263 | 767 | 36.6 | -6.7 | 94.6 | 44.6 | 73.9 | 216 | 0.06 | 0.002 | 7.6 | -0.113 | -41.2 | 0.017 |
| 1.959 | 1189 | 36.6 | -6.7 | 91.4 | 41.4 | 66.3 | 211 | 0.00 | 0.000 | 7.6 | -0.113 | -41.2 | 0.017 |
| 3.039 | 1845 | 36.7 | -6.7 | 88.3 | 38.3 | 58.7 | 206 | 0.00 | 0.000 | 7.7 | -0.110 | -41.4 | 0.017 |
| 4.951 | 3006 | 37.0 | -7.1 | 84.6 | 34.6 | 50.4 | 205 | 0.00 | 0.000 | 8.0 | -0.095 | -42.7 | 0.017 |
| 9.802 | 5951 | 40.0 | -10.2 | 78.1 | 28.1 | 38.8 | 202 | 0.00 | 0.000 | 11.0 | -0.024 | -54.7 | 0.017 |
| 19.940 | 1494 | 52.6 | -22.7 | 66.3 | 16.3 | 25.7 | 192 | 0.00 | 0.000 | 23.6 | -0.000 | -105.2 | 0.017 |

---

Run 2026-10-29T00:24:34+00:00.  Provenance is in `manifest.json` next to this file.  Rerunning `tss study leakage_receive_chain` from the same commit and the same configuration reproduces it.
