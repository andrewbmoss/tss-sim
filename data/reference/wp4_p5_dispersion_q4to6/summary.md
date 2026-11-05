# wp4_dispersion

WP4 Phase 3: dispersion and step-transition treatments compared at echo level and through the chain

Status: **implemented**.  Work packages: WP4.  Claims: C1, C2.  Requirements: INS-05, INS-07, INS-10, SCI-05, SCI-07, SCI-09.

# wp4_dispersion

18 scenes x O, X x 3 seeds; echo-level, clutter off.

| Treatment | Loss dB med (p90) | Width km med (p90) | Edge err km med | Detected | Cell km max | Loss near f_N | Loss near cusp | Steps | Sweep s |
|---|---|---|---|---|---|---|---|---|---|
| D0 | 3.8 (9.0) | 3.5 (10.5) | -5.4 | 0.981 | 8.0 | 8.0 | 7.7 | 159 | 4.26 |
| D1x2 | 6.3 (8.9) | 4.3 (8.9) | -7.6 | 0.983 | 15.9 | 7.9 | 8.2 | 159 | 4.26 |
| D1x4 | 8.3 (11.9) | 9.5 (20.8) | -10.9 | 0.981 | 31.9 | 11.3 | 11.4 | 159 | 4.26 |
| D2 | 0.1 (5.0) | 1.5 (4.2) | -2.6 | 0.987 | 8.0 | 1.5 | 1.6 | 159 | 4.26 |
| D2w | 0.1 (2.0) | 1.5 (3.2) | -2.6 | 0.990 | 8.0 | 0.4 | 0.7 | 159 | 4.26 |
| D3 | 0.3 (4.4) | 4.5 (7.5) | -5.4 | 0.983 | 8.0 | 4.5 | 4.5 | 237 | 6.35 |
| T1 | 3.0 (7.9) | 3.4 (8.8) | -3.7 | 0.981 | 10.8 | 6.2 | 7.3 | 159 | 4.26 |
| T2 | 3.8 (9.0) | 3.6 (10.5) | -5.4 | 0.980 | 8.0 | 8.0 | 7.7 | 159 | 4.26 |

Loss by band (median dB):

| Treatment | 0-0.5 MHz | 0.5-2 MHz | 2-5 MHz | 5-10 MHz | 10-20 MHz |
|---|---|---|---|---|---|
| D0 | nan | 3.1 | 4.8 | 3.1 | 4.0 |
| D1x2 | nan | 5.3 | 6.6 | 6.1 | 6.2 |
| D1x4 | nan | 7.4 | 9.0 | 7.8 | 8.3 |
| D2 | nan | 0.4 | 0.2 | 0.1 | 0.3 |
| D2w | nan | 0.3 | 0.2 | 0.1 | 0.2 |
| D3 | nan | 3.2 | 0.3 | 0.2 | 4.3 |
| T1 | nan | 2.4 | 3.8 | 2.5 | 3.3 |
| T2 | nan | 3.1 | 4.8 | 3.1 | 4.0 |

Chain (wp4_campaign, clutter off):

| Treatment | SCI-05 med | SCI-07 km med | INS-10 dB med (worst) | INS-05 worst km | O complete | SCI-05 verdicts | run s |
|---|---|---|---|---|---|---|---|
| D0 | 0.037 | 4.1 | 14.7 (-14.1) | 13.6 | 0.82 | {'fail': 3, 'goal': 12, 'threshold': 3} | 3.0 |
| D1x2 | 0.056 | 10.8 | 11.9 (-8.8) | 27.3 | 0.38 | {'fail': 5, 'goal': 12, 'threshold': 1} | 3.0 |
| D1x4 | 0.099 | 13.3 | 15.9 (9.8) | 54.5 | 0.20 | {'fail': 6, 'goal': 9, 'threshold': 3} | 3.0 |
| D2 | 0.030 | 2.8 | 20.4 (-4.6) | 13.6 | 0.97 | {'goal': 17, 'threshold': 1} | 6.8 |
| D2w | 0.031 | 5.1 | 22.2 (9.8) | 13.6 | 0.95 | {'goal': 18} | 11.7 |
| D3 | 0.041 | 5.9 | 17.9 (11.7) | 13.6 | 0.84 | {'fail': 3, 'goal': 15} | 3.5 |
| T1 | 0.033 | 2.4 | 16.3 (-9.9) | 18.5 | 0.77 | {'fail': 4, 'goal': 12, 'threshold': 2} | 2.1 |
| T2 | 0.037 | 4.1 | 14.7 (-14.1) | 13.6 | 0.80 | {'fail': 4, 'goal': 12, 'threshold': 2} | 2.2 |

---

Run 2026-11-03T21:04:02+00:00.  Provenance is in `manifest.json` next to this file.  Rerunning `tss study wp4_dispersion` from the same commit and the same configuration reproduces it.
