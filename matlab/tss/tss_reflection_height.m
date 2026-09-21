function [h_refl, status] = tss_reflection_height(f_hz, cfg, tol)
%   Layer:        L1
%   WP:           WP2
%   Claims:       C1
%   Status:       implemented
%   Requirements: INS-01
%
%   Altitude at which a wave of frequency f_hz reflects.
%
%   Inputs
%     f_hz   sounding frequency, Hz.  Scalar.
%     cfg    configuration struct.
%     tol    bisection tolerance in metres.  Default 1e-3.
%
%   Outputs
%     h_refl   reflection altitude, m, or NaN when there is none.
%     status   'ok', 'evanescent' or 'penetrates'.
%
%   Example
%     [h, st] = tss_reflection_height(5e6, tss_config());
%
%   The 'evanescent' return is the physical floor on the sounding band and is
%   the quantity the minimum sounding frequency study measures over IRI-2020
%   profiles.  See TSS_STUDY_IRI_PAIRS and the L0 CSV that TSS_IRI_PROFILES
%   writes.
%
%   Migrated into tss-sim from the author's LinkBudget/matlab tree.  The
%   body below is unchanged; only this header was added.  Original help
%   text follows.
%
%TSS_REFLECTION_HEIGHT  Altitude at which a wave of frequency f_hz reflects.
%
%   [h_refl, status] = TSS_REFLECTION_HEIGHT(f_hz, cfg) searches downward from
%   the spacecraft for the altitude where the local plasma frequency equals the
%   sounding frequency, locating it by bisection.
%
%   status is one of:
%     'ok'          h_refl is the reflection altitude in metres
%     'evanescent'  the plasma frequency at the spacecraft already exceeds
%                   f_hz, so the wave cannot propagate away from the antenna at
%                   all.  This is the physical floor on the sounding band --
%                   the open item the science requirements flag as needing an
%                   IRI/NeQuick study to fix the minimum sounding frequency.
%     'penetrates'  f_hz exceeds foF2; there is no topside echo.
%
%   The bisection matters: the absorption integrand rises steeply just below the
%   reflection point, so a reflection height quantized to a coarse search grid
%   puts visible jitter on absorption versus frequency.

if nargin < 3, tol = 1e-3; end
c = tss_constants();
fn = @(h) c.plasma_const * sqrt(tss_ne(h, cfg)) - f_hz;

h_floor = cfg.hmf2;
h_refl  = NaN;

if fn(cfg.h_sat) >= 0
    status = 'evanescent';  return
end
if fn(h_floor) < 0
    status = 'penetrates';  return
end

lo = h_floor;  hi = cfg.h_sat;          % fn(lo) >= 0, fn(hi) < 0
n_iter = ceil(log2(max(cfg.h_sat - h_floor, 1.0) / tol)) + 4;
for k = 1:n_iter
    mid = 0.5*(lo + hi);
    if fn(mid) >= 0, lo = mid; else, hi = mid; end
end
h_refl = 0.5*(lo + hi);
status = 'ok';
end
