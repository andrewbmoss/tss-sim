# wp4_dispersion

WP4 Phase 3: dispersion and step-transition treatments compared at echo level and through the chain

Status: **implemented**.  Work packages: WP4.  Claims: C1, C2.  Requirements: INS-05, INS-07, INS-10, SCI-05, SCI-07, SCI-09.

# wp4_dispersion

18 scenes x O, X x 3 seeds; echo-level, clutter off.

| Treatment | Loss dB med (p90) | Width km med (p90) | Edge err km med | Detected | Cell km max | Loss near f_N | Loss near cusp | Steps | Sweep s |
|---|---|---|---|---|---|---|---|---|---|
| D0 | 4.1 (9.2) | 4.0 (10.5) | -5.9 | 0.981 | 5.7 | 8.2 | 8.2 | 159 | 4.26 |
| D1x2 | 6.3 (9.3) | 4.0 (9.4) | -8.2 | 0.987 | 11.3 | 8.1 | 8.5 | 159 | 4.26 |
| D1x4 | 8.4 (12.0) | 9.5 (20.7) | -11.4 | 0.982 | 22.7 | 11.3 | 11.5 | 159 | 4.26 |
| D2 | 0.1 (5.4) | 1.5 (4.3) | -2.8 | 0.989 | 5.7 | 1.1 | 3.3 | 159 | 4.26 |
| D2w | 0.1 (2.0) | 1.5 (3.0) | -2.8 | 0.992 | 5.7 | 0.3 | 0.6 | 159 | 4.26 |
| D3 | 0.3 (4.4) | 4.5 (7.5) | -5.3 | 0.975 | 5.7 | 5.0 | 5.6 | 237 | 6.35 |
| T1 | 3.2 (8.2) | 3.4 (9.5) | -3.8 | 0.982 | 7.7 | 6.5 | 7.6 | 159 | 4.26 |
| T2 | 4.1 (9.2) | 4.0 (10.5) | -6.0 | 0.981 | 5.7 | 8.2 | 8.2 | 159 | 4.26 |

Loss by band (median dB):

| Treatment | 0-0.5 MHz | 0.5-2 MHz | 2-5 MHz | 5-10 MHz | 10-20 MHz |
|---|---|---|---|---|---|
| D0 | nan | 3.6 | 4.9 | 3.5 | 5.0 |
| D1x2 | nan | 6.0 | 6.5 | 6.1 | 6.4 |
| D1x4 | nan | 7.8 | 9.0 | 7.9 | 9.2 |
| D2 | nan | 0.4 | 0.2 | 0.1 | 0.2 |
| D2w | nan | 0.3 | 0.1 | 0.1 | 0.1 |
| D3 | nan | 3.6 | 0.3 | 0.2 | 5.1 |
| T1 | nan | 3.1 | 3.7 | 2.7 | 4.3 |
| T2 | nan | 3.6 | 4.9 | 3.5 | 5.0 |

Chain (wp4_campaign, clutter off):

| Treatment | SCI-05 med | SCI-07 km med | INS-10 dB med (worst) | INS-05 worst km | O complete | SCI-05 verdicts | run s |
|---|---|---|---|---|---|---|---|
| D0 | 0.177 | 18.9 | 14.7 (-14.1) | 13.6 | 0.12 | {'fail': 9, 'goal': 3, 'threshold': 6} | 0.9 |
| D1x2 | 0.200 | 20.8 | 11.9 (-8.8) | 27.3 | 0.10 | {'fail': 10, 'goal': 3, 'threshold': 5} | 0.9 |
| D1x4 | 0.208 | 22.5 | 15.9 (9.8) | 54.5 | 0.14 | {'fail': 9, 'goal': 3, 'threshold': 6} | 1.0 |
| D2 | 0.147 | 17.4 | 20.2 (-4.6) | 13.6 | 0.14 | {'fail': 9, 'goal': 6, 'threshold': 3} | 3.3 |
| D2w | 0.146 | 17.4 | 20.2 (9.2) | 13.6 | 0.16 | {'fail': 9, 'goal': 6, 'threshold': 3} | 5.4 |
| D3 | 0.172 | 18.2 | 17.9 (11.7) | 13.6 | 0.05 | {'fail': 9, 'goal': 3, 'threshold': 6} | 1.4 |
| T1 | 0.153 | 18.9 | 16.3 (-9.9) | 18.5 | 0.16 | {'fail': 9, 'goal': 4, 'threshold': 5} | 0.9 |
| T2 | 0.177 | 18.9 | 14.7 (-14.1) | 13.6 | 0.12 | {'fail': 9, 'goal': 3, 'threshold': 6} | 1.2 |

---

Run 2026-10-31T13:14:26+00:00.  Provenance is in `manifest.json` next to this file.  Rerunning `tss study wp4_dispersion` from the same commit and the same configuration reproduces it.
