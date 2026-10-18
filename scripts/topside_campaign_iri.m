function topside_campaign_iri(points_csv, alt_csv, out_mat, pharlap_root)
%TOPSIDE_CAMPAIGN_IRI  Fortran IRI-2020 half of the PyIRI topside campaign.
%
%   TOPSIDE_CAMPAIGN_IRI(POINTS_CSV, ALT_CSV, OUT_MAT) evaluates PHaRLAP's
%   iri2020 mex at every row of POINTS_CSV and saves the electron density
%   profiles and the full iono_extra vector to OUT_MAT.
%
%   POINTS_CSV  headerless, one row per profile:
%               year, month, day, hour_utc, minute_utc, lat_deg, lon_deg, R12
%               R12 is per row, so one batch spans several solar levels.
%   ALT_CSV     one column, km, uniformly spaced, 65 to 2000 km.
%   OUT_MAT     written with variables ne (n x n_alt, m^-3, single, NaN for
%               IRI's -1 sentinel), extra (n x 100, iono_extra, NaN for -1),
%               alt_km, iri_options, elapsed_s, n_failed.
%
%   The options are IRI-2020's defaults, stated explicitly so that the record
%   does not depend on what the wrapper happens to default to: NeQuick
%   topside, URSI foF2, Shubin-COSMIC (SHU-2015) hmF2, ABT-2009 B0 and B1,
%   F10.7 limited to 188 for the topside.  They match scripts/pyiri_validation
%   and Dissertation/MinFreqStudy/matlab/iri_batch.m.
%
%   Called by scripts/topside_campaign.py.  Andrew B. Moss -- PhD
%   dissertation, Utah State University.

if nargin < 4 || isempty(pharlap_root)
    pharlap_root = getenv('PHARLAP_ROOT');
end
if isempty(pharlap_root)
    pharlap_root = 'C:\Users\Andrew\PhD\Chartier Verification\pharlap_4.7.4';
end
addpath(fullfile(pharlap_root, 'src', 'matlab'), fullfile(pharlap_root, 'dat'));
addpath(fullfile(pharlap_root, 'mex'));
setenv('DIR_MODELS_REF_DAT', fullfile(pharlap_root, 'dat'));

pts    = readmatrix(points_csv);
alt_km = readmatrix(alt_csv);  alt_km = alt_km(:).';
na = numel(alt_km);  ht_step = alt_km(2) - alt_km(1);
if any(abs(diff(alt_km) - ht_step) > 1e-6)
    error('topside_campaign_iri:grid', 'altitude grid must be uniform');
end

iri_options = struct();
iri_options.iri_messages     = 'off';
iri_options.Ne_topside_model = 'NeQuick';
iri_options.foF2_coeffs      = 'URSI';
iri_options.Ne_B0B1_model    = 'ABT-2009';
iri_options.hmF2_model       = 'Shubin-COSMIC';
iri_options.Ne_tops_limited  = 'f10.7 limited';
iri_options.Ne_profile_calc  = 'standard';

n = size(pts, 1);
ne = nan(n, na, 'single');
extra = nan(n, 100);
n_failed = 0;
t0 = tic;
for k = 1:n
    UT = pts(k, 1:5);
    try
        [iono, e] = iri2020(pts(k,6), pts(k,7), pts(k,8), UT, alt_km(1), ...
                            ht_step, na, iri_options);
        row = iono(1, 1:na);  row(row < 0) = NaN;
        ne(k, :) = single(row);
        e = e(:).';  e(e == -1) = NaN;
        m = min(100, numel(e));
        extra(k, 1:m) = e(1:m);
    catch ME
        n_failed = n_failed + 1;
        if n_failed <= 5
            fprintf(2, 'point %d failed: %s\n', k, ME.message);
        end
    end
    if mod(k, 2000) == 0
        fprintf('%7d / %7d  %6.1f s\n', k, n, toc(t0));
    end
end
elapsed_s = toc(t0);
save(out_mat, 'ne', 'extra', 'alt_km', 'iri_options', 'elapsed_s', 'n_failed', '-v7');
fprintf('topside_campaign_iri: %d profiles in %.1f s, %d failed\n', n, elapsed_s, n_failed);
end
