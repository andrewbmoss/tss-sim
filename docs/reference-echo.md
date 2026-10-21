# The reference echo, and the INS-10 minimum on PHaRLAP (WP3 Phase 0)

WP3 Phase 0, 2026-10-20.  Plan: project note `claude/wp3-plan.md`.  Code:
`tss.l2_instrument.reference_echo`, study `reference_echo`,
`scripts/ins10_trace.py`.  Outputs: `data/reference/wp3_reference_echo`.

Phase 0 had four jobs: trace the 7 dB INS-10 drop on PHaRLAP, resolve the
`baseline_6u` -20 dB, confirm that the link budget carries the Balanis
reactance, and fix "the weakest echo INS-10 admits" as a tested number per
band and per configuration, because INS-24 (and so INS-22 and INS-11) is
measured against it.

## 1. The -20 dB in `baseline_6u.yaml`

Already resolved on 2026-10-17 and tested.  `transmit.efficiency_db` is 0 dB.
The matching-coil loss lives only in the antenna's radiation efficiency, and
that efficiency follows the sounding frequency
(`test_review_fixes.py::test_the_matching_loss_is_counted_once_in_the_antenna_efficiency`
and the two tests after it).  Nothing further to do.

## 2. Does the link budget carry the Balanis reactance?

The L2 baseline does.  The `link_budget` study and the 122 dB figure do not.

* **The L2 baseline.**  It builds `ShortDipole` from the instrument file.  The
  5 m dipole has -10.9 dB efficiency at 5 MHz, which is Rev G §5.8
  (`test_l2_antenna.py`, and `test_reference_echo.py` checks it again).
* **The `link_budget` study.**  It ran `linkbudget.SounderConfig()` with its
  defaults.  That is the inherited TISC point design:
  * the legacy antenna curve (-10 dBi at 1 MHz, -6.8 dBi at 5 MHz);
  * a reflection 300 km away;
  * a Chapman layer that removes most of the band as evanescent or
    penetrating.

  That is where the SRM's INS-10 status, "median 38.5 dB, minimum 33.1 dB",
  comes from.  The study now also reports the instrument
  (`instrument_snr_at_rmax.csv` and `instrument_*` metrics).  It keeps the
  inherited numbers for the attribution table and labels them as inherited.
* **The 122 dB.**  SRM Derivations §16 row 179 takes "leakage over echo =
  122 - isolation" from the M&S plan, and INS-22 (§16) and the INS-11
  high-pass rejection (§22) are built on it.  The inherited budget's echo at
  5 MHz, 300 km and -6.8 dBi is 122.5 dB below the transmitter
  (`reference_echo.legacy_echo_db`), so the 122 dB is that budget.  It is not
  at INS-04's 1050 km and not on the instrument's antenna.

## 3. The reference

`reference_levels(instrument, sweep)` uses the same terms the L2 baseline
puts into the samples:
* transmit power;
* the antenna gain at both ends;
* free-space spreading over 2R;
* the -6 dB polarization loss;
* `T_ant = eta T_sky + (1 - eta) T_phys` plus the receiver's `(F - 1) T0`;
* a cell bandwidth of `ENBW / T_d`, which is 1.5 / 26.8 ms = 56 Hz with Hann.

It returns, per step:

| Quantity | Definition |
|---|---|
| `noise_cell_dbw` | `k T_sys ENBW / T_d` at the receiver input |
| `weakest_echo_dbw` | `noise_cell_dbw + 10 dB`: the weakest echo INS-10 admits.  This is the reading of INS-24 implemented here. |
| `echo_rmax_dbw` | the free-space echo from 1050 km (no absorption): the echo the 122 dB stood in for, now at the right range on the right antenna |
| `leakage_dbw` | transmit power less INS-23's isolation |
| ratios | leakage over each of the three |
| echo band | `2 K_i R / c` for R from 100 km to 1050 km at the step's rate: where INS-24 applies |

A1, 40 dB isolation (dBm at the receiver input):

| f (MHz) | step rate (MHz/s) | echo band (kHz) | G (dBi) | noise per cell | weakest echo | echo at 1050 km | SNR at 1050 km (dB) | L / weakest (dB) | L / echo at 1050 km (dB) | L / noise per cell (dB) |
|---|---|---|---|---|---|---|---|---|---|---|
| 0.338 | 0.616 | 0.41-4.3 | -44.0 | -148.8 | -138.8 | -153.4 | -4.6 | 128.8 | 143.4 | 138.8 |
| 0.500 | 0.910 | 0.61-6.4 | -38.9 | -146.7 | -136.7 | -146.6 | 0.1 | 126.7 | 136.6 | 136.7 |
| 0.989 | 1.801 | 1.20-12.6 | -30.0 | -142.4 | -132.4 | -134.8 | 7.6 | 122.4 | 124.8 | 132.4 |
| 1.959 | 3.566 | 2.38-25.0 | -21.1 | -137.6 | -127.6 | -122.9 | 14.6 | 117.6 | 112.9 | 127.6 |
| 3.039 | 5.532 | 3.69-38.7 | -15.4 | -134.4 | -124.4 | -115.4 | 19.0 | 114.4 | 105.4 | 124.4 |
| 4.988 | 5.597 | 3.73-39.2 | -9.2 | -131.0 | -121.0 | -107.2 | 23.8 | 111.0 | 97.2 | 121.0 |
| 9.938 | 5.597 | 3.73-39.2 | -1.9 | -127.7 | -117.7 | -98.7 | 29.0 | 107.7 | 88.7 | 117.7 |
| 19.957 | 3.235 | 2.16-22.7 | 1.3 | -128.4 | -118.4 | -98.3 | 30.1 | 108.4 | 88.3 | 118.4 |

Every leakage ratio moves one for one with the isolation.  A2 has the same
levels per cell, because the dwell is the same, but a wider echo band: to
247 kHz at its peak rate against A1's 39.2 kHz.

**Tested.**  `tests/unit/test_reference_echo.py`:
* every term against an independent calculation;
* the 122 dB reproduced from the inherited budget;
* the SRM D188 spur rule reproduced (-50.3 dBc);
* a synthetic echo at 1050 km put through the L2 baseline and measured by L3
  exactly as INS-10 is scored.

The study runs that last check at six frequencies, with a predicted 20 dB:

| f (MHz) | 0.355 | 0.989 | 2.50 | 4.99 | 9.94 | 18.0 |
|---|---|---|---|---|---|---|
| loss (dB) | +0.66 | +0.89 | +0.99 | +0.07 | +0.68 | -0.25 |

* **Below 3 MHz** the loss is mostly the step-transition skirt.  In the first
  7 ms of each 26.8 ms dwell the echo beats at the previous step's rate, which
  costs about 0.9 dB with the Hann window.
* **Above 3 MHz** what is left is the scatter of four noise draws.

The reference is the level the chain measures, to within 1 dB.

## 4. What the reference does to the derivations built on 122 dB

At 40 dB isolation, leakage over the echo is 97.2 dB at 5 MHz on the
instrument at 1050 km, not 82 dB.  Against the weakest echo INS-10 admits it
is 111.0 dB at 5 MHz, and 107.5 to 128.8 dB across the band.  Carried through
the two Rev G rules (A1, 10 m path, 10 dB margin, single-tone small-angle):

| f (MHz) | INS-11 HPF, Rev G rule (dB) | INS-11 HPF, noise rule (dB) | INS-22 spur vs echo at 1050 km (dBc) | INS-22 spur vs weakest echo (dBc) |
|---|---|---|---|---|
| 0.338 | 93.4 | 57.3 | -92.6 | -77.9 |
| 0.500 | 86.6 | 55.2 | -89.2 | -79.2 |
| 0.989 | 74.8 | 50.9 | -83.2 | -80.8 |
| 1.959 | 62.9 | 46.1 | -77.3 | -81.9 |
| 3.039 | 55.4 | 43.0 | -73.6 | -82.6 |
| 4.988 | 47.2 | 39.6 | -65.5 | -79.3 |
| 9.938 | 38.7 | 36.3 | -57.0 | -76.0 |
| 19.957 | 38.3 | 37.0 | -51.8 | -72.0 |

* **The Rev G rule** is SRM §22: leakage over echo, less 60 dB, plus 10 dB.
  It gives 32 dB from 82 dB.
* **The noise rule** refers the digitiser's 60 dB to the noise in the whole
  digitiser band (78.4 kHz, 31.5 dB above one cell).  That is what the 60 dB
  actually has to span.
* **INS-22 against the weakest echo**, worst step:
  * -62.6 dBc at a 1 m path;
  * -82.6 dBc at 10 m;
  * -92.2 dBc at 30 m;
  * A2 at 10 m: -92.6 dBc.

  A longer path suppresses less, because the suppression grows as
  `(pi f_b tau)^2`.

**Rev G states -70 dBc and 32 dB.**  Both were derived from the 122 dB.  They
are not changed here (WP3 plan, G5).  The numbers above are for the decision.

## 5. INS-10 "at R_max across the band"

INS-10 is worded "SNR >= 10 dB per range cell at R_max across the band".  On
the 5 m antenna the free-space SNR at 1050 km is:
* 23.8 dB at 5 MHz;
* 27.8 dB median over the 159 steps;
* below 10 dB under 1.26 MHz, falling to -4.6 dB at 0.338 MHz.

Only transmit efficiency costs SNR while the sky dominates.  Receive
efficiency stops cancelling at the bottom of the band, where `eta T_sky` falls
to the receiver's own noise.  At those frequencies the real echoes come from
far nearer than 1050 km: the first PHaRLAP echo at 800 km is at 3.04 MHz and
111 km.  So the literal wording fails at the band bottom on any antenna short
enough to fly.

## 6. The INS-10 minimum on PHaRLAP: the cusp, not the topside

`tss run baseline` scored INS-10 at 7.0 dB at 800 km on PHaRLAP's IRI-2020,
against 18.7 dB on PyIRI.  The question left open (RESUME, 2026-10-18) was
whether the denser IRI-2020 topside is responsible.  It is not.

**Which echo.**  On both codes, and in all three configurations, the minimum
is the X echo at 10.988 MHz.  That is the last step below fxF2 (A1 steps are
150 kHz apart there).

| | PHaRLAP | PyIRI |
|---|---|---|
| foF2 at the sub-satellite point (0 N, 0 E, 12 UT) | 10.637 MHz | 10.674 MHz (+0.35 %) |
| hmF2 | 405 km | 400 km |
| fxF2 (from the field strength at hmF2) and the 10.988 MHz step's distance below it | 10.9985 MHz, 10 kHz | 11.036 MHz, 48 kHz |
| 10.988 MHz X echo: virtual range | 873 km | 610 km |
| its absorption (two-way) | 6.2 dB | 3.7 dB |
| its SNR | 7.0 dB | 18.7 dB |
| the next X echo down (10.838 MHz) | 656 km, 17.0 dB | 563 km, 24.4 dB |
| median SNR over all truth echoes | 30.8 dB | 30.9 dB |

The medians agree to 0.1 dB, so the topside does not move the sensitivity.
On PHaRLAP, foF2 is 0.35 % lower, so the last step sits closer to the X
cutoff.  There:
* the group delay diverges, putting the echo at 873 km;
* the deviative absorption rises;
* the delay changes fastest within the dwell, so in-dwell dispersion smears
  the echo across cells.

The path gain is 4 dB below the neighbouring echo's.  The chain loses about
6 dB more to the smear.

**Proof: move the step, not the topside.**  `scripts/ins10_trace.py cusp`
scales the density by 0.97 to 1.00, which moves foF2 by 1.5 %, one step.  It
uses one fixed noise draw.  At 800 km:

| | PHaRLAP | PyIRI |
|---|---|---|
| INS-10 minimum over the 13 scalings | 5.6 to 13.2 dB | 12.2 to 18.9 dB |
| minimum without each mode's last step below its cutoff | 11.9 to 18.6 dB | 16.2 to 20.9 dB |
| median | 29.2 to 31.0 dB | 30.5 to 31.1 dB |

At 500 km on PHaRLAP the same sweep gives a minimum of 7.5 to 19.2 dB.
Without the last step it is 21.5 to 27.1 dB.

The minimum is a sawtooth in foF2, with the period of one step.  Two more
things move it by several dB:
* **The noise draw.**  The same PHaRLAP profile gives 7.0 dB with the
  campaign's generator and 10.8 dB with the fixed one, because an echo near
  10 dB SNR is measured as a peak over noise.
* **O/X interference in one channel.**  The same steps sometimes reappear as
  the minimum at mid-band (5.74, 7.24 and 7.99 MHz, about 10 dB below their
  neighbours, with the O and X readings identical).  These are steps where the
  O and X echoes fall in the same range cell and add with opposing phase in a
  single receive channel.  That is the A1 configuration's own limitation.
  A2's polarization separation is meant to remove it, but it is not yet
  modelled, because the vector channels are still placeholders.

**What it means.**  Scored as the minimum over every truth echo, INS-10 is
decided by three things:
* where the last step lands below each cutoff, which is a matter of luck at
  the 0.1 % level in foF2;
* the noise draw on that one echo;
* O/X fading.

It is not decided by sensitivity.  Any instrument has an echo arbitrarily
close to the cusp, where the group range diverges, so no fixed SNR holds
there.  The 2026-10-17 decision to score at every truth echo stands; how to
treat the cusp is the decision this raises (WP3 plan, G5; project note
`claude/wp3-p0-reference-echo.md`).
