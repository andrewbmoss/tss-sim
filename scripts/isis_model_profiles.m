function isis_model_profiles(points_csv, out_mat, pharlap_root)
%ISIS_MODEL_PROFILES  IRI-2020 profiles at the ISIS-1 records (WP4 Phase 6).
%
%   ISIS_MODEL_PROFILES(POINTS_CSV, OUT_MAT) evaluates PHaRLAP's iri2020 mex
%   from 100 to 1200 km in 5 km steps at the time and place of every row of
%   POINTS_CSV, on the epoch's own indices (R12 = -1), once per topside:
%
%     NeQuick            IRI-2020's default topside, foF2 storm model on
%     IRI2001            IRI-2001 topside, foF2 storm model on
%     IRI2001corr        IRI-2001 corrected topside, foF2 storm model on
%
%   These are the epoch configurations of topside_isis_epoch_iri.m (WP0),
%   evaluated over the whole height range instead of at the spacecraft only.
%
%   POINTS_CSV  headerless: year, month, day, hour_utc, minute_utc, lat_deg,
%               lon_deg (one row per record, in isis1_retrieval's order)
%   OUT_MAT     ne (n x 221 x 3, single, m^-3, NaN where IRI gives -1),
%               alt_km (1 x 221), nmf2, hmf2 (n x 3), conf_names, elapsed_s,
%               n_failed.
%   Andrew B. Moss -- PhD dissertation, Utah State University.

if nargin < 3 || isempty(pharlap_root)
    pharlap_root = 'C:\Users\Andrew\PhD\Chartier Verification\pharlap_4.7.4';
end
addpath(fullfile(pharlap_root, 'src', 'matlab'), fullfile(pharlap_root, 'dat'));
addpath(fullfile(pharlap_root, 'mex'));
setenv('DIR_MODELS_REF_DAT', fullfile(pharlap_root, 'dat'));

base = struct();
base.iri_messages     = 'off';
base.foF2_coeffs      = 'URSI';
base.Ne_B0B1_model    = 'ABT-2009';
base.hmF2_model       = 'Shubin-COSMIC';
base.Ne_tops_limited  = 'f10.7 limited';
base.Ne_profile_calc  = 'standard';
base.Ne_topside_model = 'NeQuick';
base.foF2_storm       = 'on';
base.hmF2_storm       = 'off';
base.foE_storm        = 'off';
base.topside_storm    = 'off';

conf_names = {'NeQuick', 'IRI2001', 'IRI2001corr'};
opts = repmat(base, 1, 3);
opts(2).Ne_topside_model = 'IRI-2001';
opts(3).Ne_topside_model = 'IRI-2001 corrected';

h0 = 100; dh = 5; nh = 221;
alt_km = h0 + dh * (0:nh-1);
pts = readmatrix(points_csv);
n = size(pts, 1);  nc = numel(conf_names);
ne = nan(n, nh, nc, 'single');
nmf2 = nan(n, nc); hmf2 = nan(n, nc);
n_failed = zeros(1, nc);
t0 = tic;
for c = 1:nc
    for k = 1:n
        UT = pts(k, 1:5);
        try
            [iono, e] = iri2020(pts(k,6), pts(k,7), -1, UT, h0, dh, nh, opts(c));
            v = iono(1, :);  v(v <= 0) = NaN;
            ne(k, :, c) = single(v);
            if e(1) > 0, nmf2(k, c) = e(1); end
            if e(2) > 0, hmf2(k, c) = e(2); end
        catch ME
            n_failed(c) = n_failed(c) + 1;
            if n_failed(c) <= 3
                fprintf(2, '%s record %d failed: %s\n', conf_names{c}, k, ME.message);
            end
        end
    end
    fprintf('%-12s done  %7.1f s  %d failed\n', conf_names{c}, toc(t0), n_failed(c));
end
elapsed_s = toc(t0);
save(out_mat, 'ne', 'alt_km', 'nmf2', 'hmf2', 'conf_names', 'elapsed_s', 'n_failed', '-v7');
fprintf('isis_model_profiles: %d records x %d topsides in %.1f s\n', n, nc, elapsed_s);
end
