# reference_echo

INS-10/INS-24 reference echo, leakage ratios and the echo band, A1 and A2

Status: **implemented**.  Work packages: WP3.  Claims: C1, C4.  Requirements: INS-10, INS-11, INS-22, INS-23, INS-24.

The reference INS-24 is measured against, computed from the baseline instrument
(1 W, 5 m dipole with the Balanis reactance, 3 dB noise figure, Franklin and
MacLean sky, Hann window, 26.8 ms dwell, 40 dB isolation).  A1 sweep; A2 has the
same levels per cell because the dwell is the same, and a wider echo band
(to 247 kHz against 39.2 kHz).

| f (MHz) | step rate (MHz/s) | echo band (kHz) | G (dBi) | N cell (dBm) | weakest echo (dBm) | free-space echo at 1050 km (dBm) | SNR at 1050 km (dB) | L/weakest (dB) | L/echo(1050) (dB) | L/N cell (dB) |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.338 | 0.616 | 0.41-4.3 | -44.0 | -148.8 | -138.8 | -153.4 | -4.6 | 128.8 | 143.4 | 138.8 |
| 0.500 | 0.910 | 0.61-6.4 | -38.9 | -146.7 | -136.7 | -146.6 | 0.1 | 126.7 | 136.6 | 136.7 |
| 0.989 | 1.801 | 1.20-12.6 | -30.0 | -142.4 | -132.4 | -134.8 | 7.6 | 122.4 | 124.8 | 132.4 |
| 1.959 | 3.566 | 2.38-25.0 | -21.1 | -137.6 | -127.6 | -122.9 | 14.6 | 117.6 | 112.9 | 127.6 |
| 3.039 | 5.532 | 3.69-38.7 | -15.4 | -134.4 | -124.4 | -115.4 | 19.0 | 114.4 | 105.4 | 124.4 |
| 4.988 | 5.597 | 3.73-39.2 | -9.2 | -131.0 | -121.0 | -107.2 | 23.8 | 111.0 | 97.2 | 121.0 |
| 9.938 | 5.597 | 3.73-39.2 | -1.9 | -127.7 | -117.7 | -98.7 | 29.0 | 107.7 | 88.7 | 117.7 |
| 19.957 | 3.235 | 2.16-22.7 | 1.3 | -128.4 | -118.4 | -98.3 | 30.1 | 108.4 | 88.3 | 118.4 |

The old figure.  SRM §16 takes leakage over echo as 122 dB less the isolation.
The inherited link budget (`linkbudget.SounderConfig` defaults: 300 km, the
legacy -6.8 dBi antenna) puts its 5 MHz echo 122.5 dB
below the transmitter, which is that 122 dB.  On the instrument, at INS-04's
1050 km, the free-space echo at 5 MHz is
137.2 dB below
the transmitter, so leakage over echo is 97.2 dB
at 40 dB isolation, not 82.  Against the weakest echo INS-10 admits it is
111.0 dB at 5 MHz and
107.5 to 128.8 dB
over the band.  All of these move one for one with the isolation.

What that does to the two derivations built on the old figure (A1, 10 m coupling
path, 10 dB margin, single-tone small-angle as in Rev G §5.9; INS-11 at 60 dB
with 10 dB headroom):

| f (MHz) | INS-11 HPF, Rev G rule (dB) | INS-11 HPF, noise rule (dB) | INS-22 spur, 10 m, vs echo at 1050 km (dBc) | INS-22 spur, 10 m, vs weakest echo (dBc) |
|---|---|---|---|---|
| 0.338 | 93.4 | 57.3 | -92.6 | -77.9 |
| 0.500 | 86.6 | 55.2 | -89.2 | -79.2 |
| 0.989 | 74.8 | 50.9 | -83.2 | -80.8 |
| 1.959 | 62.9 | 46.1 | -77.3 | -81.9 |
| 3.039 | 55.4 | 43.0 | -73.6 | -82.6 |
| 4.988 | 47.2 | 39.6 | -65.5 | -79.3 |
| 9.938 | 38.7 | 36.3 | -57.0 | -76.0 |
| 19.957 | 38.3 | 37.0 | -51.8 | -72.0 |

The "Rev G rule" is SRM §22, leakage over echo less 60 dB plus 10 dB.  The
"noise rule" refers the digitiser's 60 dB to the noise in the full digitiser
band (78.4 kHz, 31.5 dB above one
cell), which is what the 60 dB has to span.  Which rule INS-11 uses, and which
echo INS-22 and INS-24 refer to, are WP3 decisions (plan, G5); no requirement
document is changed by this study.

Free-space SNR at 1050 km falls below INS-10's 10 dB below
1.26 MHz on the 5 m antenna
(minimum -4.6 dB at the band edge), because only
transmit efficiency costs SNR and the receive efficiency stops cancelling where
the sky no longer dominates.  INS-10 is worded "at R_max across the band"; at
the bottom of the band the echoes come from far nearer than 1050 km.

Chain check.  A synthetic echo at 1050 km set to a predicted 20 dB, through the
L2 baseline (noise, receiver, quantiser) and measured by L3 exactly as INS-10 is
scored:

| f (MHz) | predicted SNR (dB) | measured, mean of 4 (dB) | loss (dB) |
|---|---|---|---|
| 0.355 | 20.0 | 19.34 | +0.66 |
| 0.989 | 20.0 | 19.11 | +0.89 |
| 2.500 | 20.0 | 19.01 | +0.99 |
| 4.988 | 20.0 | 19.93 | +0.07 |
| 9.938 | 20.0 | 19.32 | +0.68 |
| 18.038 | 20.0 | 20.25 | -0.25 |

Below 3 MHz the loss is mostly the step-transition skirt: for the first
tau = 7.0 ms of each 26.8 ms dwell the echo beats at the previous step's rate,
which costs about 0.9 dB with the Hann window.  Above 3 MHz the rate does not
change, and what is left is the scatter of four noise draws (about +/-0.7 dB).
The prediction and the chain agree to within 1 dB everywhere, so the reference
is the level the chain measures.

INS-22 with the coupling path.  Against the weakest admitted echo, A1 at 40 dB:
-62.6 dBc at 1 m, -82.6 dBc
at 10 m and -92.2 dBc at 30 m (worst step); A2 at 10 m
-92.6 dBc.  A longer coupling path suppresses less,
because the suppression grows as (pi f_b tau)^2.  Rev G states -70 dBc.

---

Run 2026-10-21T02:44:33+00:00.  Provenance is in `manifest.json` next to this file.  Rerunning `tss study reference_echo` from the same commit and the same configuration reproduces it.
