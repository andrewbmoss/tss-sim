function ne = tss_ne(h, cfg)
%   Layer:        L0
%   WP:           WP1
%   Claims:       C1, C2
%   Status:       implemented
%   Requirements: -
%
%   alpha-Chapman topside electron density, m^-3, at altitude h (m).
%
%   Inputs
%     h     altitude above the ellipsoid, m.  Scalar or array.
%     cfg   configuration struct; uses nmf2, hmf2 and plasma_scale_ht.
%
%   Outputs
%     ne    electron density, m^-3, the same size as h.
%
%   Example
%     ne = tss_ne(600e3, tss_config());
%
%   This is the analytic L0 stand-in.  The IRI-2020 L0 truth used by the band
%   edge study comes from TSS_IRI_PROFILES instead, and nothing downstream
%   needs to change to swap one for the other.
%
%   Migrated into tss-sim from the author's LinkBudget/matlab tree.  The
%   body below is unchanged; only this header was added.  Original help
%   text follows.
%
%TSS_NE  alpha-Chapman topside electron density, m^-3, at altitude h (m).
%
%   cfg.plasma_scale_ht is the PLASMA scale height Hp -- the quantity Meehan
%   (2017) measures as 150 +/- 30 km at night and Stankov et al. (2003) show
%   ranging 78-244 km across profiler models.  The alpha-Chapman topside decays
%   asymptotically with twice the neutral scale height, so H_chapman = Hp/2.
%
%   Replace this function with an IRI-2020 or NeQuick interpolant to run the
%   budget against a modelled profile; nothing else needs to change.

h_chap = cfg.plasma_scale_ht / 2.0;
z = (h - cfg.hmf2) ./ h_chap;
ne = cfg.nmf2 .* exp(0.5 .* (1.0 - z - exp(-z)));
end
