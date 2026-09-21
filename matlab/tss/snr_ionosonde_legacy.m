function [snrR, pwrR, nR] = snr_ionosonde_legacy(Kf, R, Rmax, f, Pt)
%   Layer:        L2
%   WP:           WP8
%   Claims:       C1
%   Status:       implemented
%   Requirements: -
%
%   The 2014 model, kept verbatim as a regression case.
%
%   Inputs
%     Kf, R, Rmax, f, Pt   as the original: Hz/s, m, m, Hz, W.
%
%   Outputs
%     snrR, pwrR, nR       dB, dBW, dBW.
%
%   Example
%     [s, p, n] = snr_ionosonde_legacy(5.6e6, 300e3, 500e3, 5e6, 1.0);
%
%   ARCHIVE.  Do not use this for results.  It is on the path so that the
%   difference between the old and the new model can be demonstrated rather
%   than silently absorbed, and TestLinkBudget pins its output so that a
%   change to it is caught.
%
%   Migrated into tss-sim from the author's LinkBudget/matlab tree.  The
%   body below is unchanged; only this header was added.  Original help
%   text follows.
%
%SNR_IONOSONDE_LEGACY  The 2014 model, kept verbatim as a regression case.
%
%   This is J. Pratt's original snr_ionosonde.m, renamed and otherwise
%   unchanged.  It is retained so the difference between the old and new models
%   can be demonstrated rather than silently absorbed.  Do not use it for
%   results: see SNR_IONOSONDE_V2 for what is wrong with it.

c = 299792458;
Tref = 290;
k = 1.38e-23;

Tx = f;
BW = 2*Kf*Rmax/c;
Lt = -3;  Gt = -6;
EIRP = 10*log10(Pt) + Lt + Gt;
lam = c./Tx;
Ls = 20*log10(lam./(4*pi*2*R));
La = -.3;  Lf = -5;  Lp = -6.51;
Gr = Gt;  Ta = 250;  Lr = -3;
Tr = (10^(-Lr/10)-1)*Tref;
Gl = 20;  Fl = 3;
Tl = (10^(Fl/10)-1)*Tref;
Len = 15;

nA = 10*log10(k*Ta.*BW) + Len;
nF = 10*log10(10.^((nA+Lr)/10) + 10^(Lr/10)*k*Tr.*BW);
nR = 10*log10(10.^((nF+Gl)/10) + 10^(Gl/10)*k*Tl.*BW);

pwrA = EIRP + Gr + Ls + La + Lf + Lp;
pwrF = pwrA + Lr;
pwrR = pwrF + Gl;
snrR = pwrR - nR;
end
