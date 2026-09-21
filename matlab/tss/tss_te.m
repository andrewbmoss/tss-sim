function te = tss_te(h, cfg)
%   Layer:        L0
%   WP:           WP1
%   Claims:       C1
%   Status:       implemented
%   Requirements: -
%
%   Electron temperature, K, at altitude h (m).
%
%   Inputs
%     h     altitude, m.  Scalar or array.
%     cfg   configuration struct; uses hmf2, te_peak and te_top.
%
%   Outputs
%     te    electron temperature, K, the same size as h.
%
%   Example
%     te = tss_te(600e3, tss_config());
%
%   Migrated into tss-sim from the author's LinkBudget/matlab tree.  The
%   body below is unchanged; only this header was added.  Original help
%   text follows.
%
%TSS_TE  Electron temperature, K, at altitude h (m).  Linear from hmF2 to 1000 km.
frac = min(max((h - cfg.hmf2) ./ (1000e3 - cfg.hmf2), 0.0), 1.5);
te = cfg.te_peak + (cfg.te_top - cfg.te_peak) .* frac;
end
