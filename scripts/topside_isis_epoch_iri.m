function topside_isis_epoch_iri(points_csv, out_mat, pharlap_root)
%TOPSIDE_ISIS_EPOCH_IRI  Fortran IRI-2020 at the ISIS-1 points with real indices.
%
%   TOPSIDE_ISIS_EPOCH_IRI(POINTS_CSV, OUT_MAT) evaluates PHaRLAP's iri2020 mex
%   at the satellite height of every row of POINTS_CSV, once per configuration
%   below, and saves the density there and the full iono_extra vector.
%
%   POINTS_CSV  headerless, one row per ISIS-1 profile:
%               year, month, day, hour_utc, minute_utc, lat_deg, lon_deg,
%               h_sat_km, R12_bucket
%   OUT_MAT     ne_sat (n x n_conf, m^-3, NaN for IRI's -1), extra
%               (n x 100 x n_conf, single, NaN for -1), conf_names,
%               conf_options, elapsed_s, n_failed.
%
%   Configurations:
%     binned_NeQuick          user R12 = the 5-unit bucket of T1 (control: it
%                             must reproduce Dissertation/ISISValidation)
%     epoch_NeQuick_quiet     R12 = -1, so R12, IG12, F10.7 daily and 81-day
%                             and ap come from IRI's own index files for the
%                             epoch; every storm model off
%     epoch_NeQuick           R12 = -1 with IRI's epoch default, the foF2
%                             storm model on
%     epoch_NeQuick_allstorm  R12 = -1, foF2, hmF2, topside and foE storm
%                             models all on
%     epoch_IRI2001           R12 = -1, foF2 storm on, IRI-2001 topside
%     epoch_IRI2001corr       R12 = -1, foF2 storm on, IRI-2001 corrected
%
%   Every other option is IRI-2020's default, stated explicitly as in
%   topside_campaign_iri.m.  Called by scripts/topside_campaign.py isis-epoch.
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
base.foF2_storm       = 'off';
base.hmF2_storm       = 'off';
base.foE_storm        = 'off';
base.topside_storm    = 'off';

conf_names = {'binned_NeQuick', 'epoch_NeQuick_quiet', 'epoch_NeQuick', ...
              'epoch_NeQuick_allstorm', 'epoch_IRI2001', 'epoch_IRI2001corr'};
opts = repmat(base, 1, numel(conf_names));
epoch = [false true true true true true];
opts(3).foF2_storm = 'on';
opts(4).foF2_storm = 'on'; opts(4).hmF2_storm = 'on';
opts(4).foE_storm  = 'on'; opts(4).topside_storm = 'on';
opts(5).foF2_storm = 'on'; opts(5).Ne_topside_model = 'IRI-2001';
opts(6).foF2_storm = 'on'; opts(6).Ne_topside_model = 'IRI-2001 corrected';
conf_options = opts;

pts = readmatrix(points_csv);
n = size(pts, 1);  nc = numel(conf_names);
ne_sat = nan(n, nc);
extra  = nan(n, 100, nc, 'single');
n_failed = zeros(1, nc);
t0 = tic;
for c = 1:nc
    for k = 1:n
        UT = pts(k, 1:5);
        if epoch(c), R12 = -1; else, R12 = pts(k, 9); end
        try
            [iono, e] = iri2020(pts(k,6), pts(k,7), R12, UT, pts(k,8), 1, 2, opts(c));
            v = iono(1, 1);
            if v > 0, ne_sat(k, c) = v; end
            e = e(:).';  e(e == -1) = NaN;
            m = min(100, numel(e));
            extra(k, 1:m, c) = single(e(1:m));
        catch ME
            n_failed(c) = n_failed(c) + 1;
            if n_failed(c) <= 3
                fprintf(2, '%s point %d failed: %s\n', conf_names{c}, k, ME.message);
            end
        end
    end
    fprintf('%-24s done  %7.1f s  %d failed\n', conf_names{c}, toc(t0), n_failed(c));
end
elapsed_s = toc(t0);
save(out_mat, 'ne_sat', 'extra', 'conf_names', 'conf_options', 'elapsed_s', 'n_failed', '-v7');
fprintf('topside_isis_epoch_iri: %d points x %d configurations in %.1f s\n', n, nc, elapsed_s);
end
