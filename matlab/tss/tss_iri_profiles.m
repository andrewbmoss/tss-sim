function info = tss_iri_profiles(outdir, points, alt_km, varargin)
%   Layer:        L0
%   WP:           WP1
%   Claims:       C1, C2
%   Status:       implemented
%   Requirements: INS-01
%
%   Batch IRI-2020 profile generation through PHaRLAP, written as CSV for L0.
%
%   Inputs
%     outdir     directory for the output files.  Created if absent.
%     points     N-by-7 numeric [year month day hour_utc minute_utc lat_deg
%                lon_deg], or a path to a headerless CSV of the same columns
%                (the layout the author's points.csv already uses).
%     alt_km     altitude grid, km, uniformly spaced, 2 to 1000 samples.  May
%                also be a path to a CSV holding one column.
%
%     'Topside'      'NeQuick' (default), 'IRI-2001' or 'IRI-2001 corrected'.
%                    These are the only topside models the PHaRLAP iri2020
%                    wrapper exposes.  IRI's Gulyaeva option is not reachable.
%     'R12'          user specified yearly smoothed sunspot number, 1 to 200.
%                    Default 105.  See WHY R12 below.
%     'IriOptions'   struct, or a JSON string, of extra iri_options fields
%                    merged over the explicit defaults set here.
%     'PharlapRoot'  installation root.  Default: TSS_PHARLAP_PATH, then
%                    PHARLAP_ROOT, then a short list of known locations.
%     'Study'        study name used in the file names.  Default 'iri'.
%     'Tag'          discriminator used in the file names.  Default is the
%                    topside model and R12, for example 'nequick_r12-105'.
%     'Progress'     print a progress line every N points.  Default 2000.
%     'Manifest'     write manifest.json beside the outputs.  Default true.
%
%   Outputs
%     info       struct with fields ne_csv, extra_csv, alt_csv, sidecar,
%                manifest, n_points, n_alt, n_failed, elapsed_s, pharlap_root,
%                iri_options.
%
%   Example
%     pts = [2015 6 21 12 0 41.7 -111.8
%            2015 6 21 12 0 -20.0  130.0];
%     info = tss_iri_profiles(tss_repo_root('data','bridge','matlab'), ...
%                             pts, 100:10:1000, 'Topside', 'NeQuick', 'R12', 105);
%
%   THE FILE CONTRACT WITH THE PYTHON L0 TIER
%
%   Three CSV files and a JSON sidecar.  All CSV files carry a header row.
%
%   <study>_profiles_<tag>.csv
%     point_id, year, month, day, hour_utc, minute_utc, lat_deg, lon_deg,
%     then one column per altitude sample named ne_m3_i001 .. ne_m3_iNNN in
%     grid order, electron density in m^-3.
%
%     The density columns are named by index and not by altitude.  An altitude
%     in a column name would be friendlier to read and it would also be a
%     second, informal copy of the grid that could drift out of step with the
%     real one.  The grid lives in exactly one place, the alt_km file, and the
%     sidecar repeats it for a reader that has only the sidecar.
%
%     The epoch is carried as five numeric columns rather than an ISO string so
%     that the whole table is numeric and can be written in one pass.  At 1.2
%     million profiles that is the difference between seconds and minutes.  The
%     columns are UTC.
%
%   <study>_extra_<tag>.csv
%     point_id, the same epoch and location columns, then nmf2_m3, hmf2_km,
%     solar_zenith_deg, dip_lat_deg, modip_deg, spread_f_prob, rz12,
%     f107_daily, kp, ap.  These come from IRI's iono_extra and are more
%     accurate than reading them off the sampled altitude grid.
%
%   <study>_alt_km_<tag>.csv
%     one column, alt_km, the authoritative grid.
%
%   IRI's -1 sentinel, which means the height is outside 65 to 2000 km or no
%   valid value was produced, is converted to NaN everywhere so that it cannot
%   be mistaken for a real density or a real scalar.  A point that raises an
%   error is recorded as an all-NaN row and the batch continues; the count of
%   such points is in info.n_failed and in the sidecar.
%
%   WHY R12 RATHER THAN AN EPOCH
%
%   Calling iri2020 with R12 = -1 makes it read R12 and IG12 from ig_rz.dat by
%   date.  That file has a finite span and the wrapper errors if the epoch falls
%   outside it, which would kill a long batch partway through.  A user specified
%   R12 removes the dependence on the index file and gives direct reproducible
%   control over solar activity, which is what a solar minimum and maximum study
%   wants.  The documented cost is that the foF2 storm model is forced off.
%   Storm sensitivity is a separate run against real epochs.
%
%   The wrapper derives the other indices as
%       F10.7 = 63.7 + 0.728*R12 + 0.00089*R12^2
%   so R12 = 10 gives F10.7 = 71 and R12 = 120 gives 164.  IRI clamps F10.7 to
%   188 for the topside calculation unless Ne_tops_limited is overridden, so an
%   R12 above about 140 buys no further topside response.  That clamp is set
%   explicitly below rather than left to the default, so the choice is visible
%   in the record.
%
%   Generalised from the author's MinFreqStudy/matlab/iri_batch.m.  The IRI call
%   sequence, the option block, the -1 handling and the per-point error handling
%   are that file's; what is new is that the points come in as an argument
%   rather than through a fixed points.csv, and that the output is the CSV
%   contract above rather than two headerless matrices.

p = inputParser();
p.FunctionName = 'tss_iri_profiles';
p.addRequired('outdir',  @(x) ischar(x) || isstring(x));
p.addRequired('points',  @(x) isnumeric(x) || ischar(x) || isstring(x) || istable(x));
p.addRequired('alt_km',  @(x) isnumeric(x) || ischar(x) || isstring(x));
p.addParameter('Topside', 'NeQuick', @(x) ischar(x) || isstring(x));
p.addParameter('R12', 105, @(x) isnumeric(x) && isscalar(x));
p.addParameter('IriOptions', struct(), @(x) isstruct(x) || ischar(x) || isstring(x));
p.addParameter('PharlapRoot', '', @(x) ischar(x) || isstring(x));
p.addParameter('Study', 'iri', @(x) ischar(x) || isstring(x));
p.addParameter('Tag', '', @(x) ischar(x) || isstring(x));
p.addParameter('Progress', 2000, @(x) isnumeric(x) && isscalar(x));
p.addParameter('Manifest', true, @(x) islogical(x) || isnumeric(x));
p.parse(outdir, points, alt_km, varargin{:});
a = p.Results;

outdir  = char(a.outdir);
topside = char(a.Topside);
R12     = double(a.R12);

pts = load_points(a.points);
alt = load_alt(a.alt_km);

n  = size(pts, 1);
na = numel(alt);

ht_start = alt(1);
ht_step  = alt(2) - alt(1);
if any(abs(diff(alt) - ht_step) > 1e-6)
    error('tss_iri_profiles:grid', 'altitude grid must be uniformly spaced');
end
if na < 2 || na > 1000
    error('tss_iri_profiles:grid', ...
        'num_hts must be between 2 and 1000 (got %d)', na);
end
if ht_start < 65 || alt(end) > 2000
    warning('tss_iri_profiles:range', ...
        'IRI returns no profile below 65 km or above 2000 km; grid is %.0f-%.0f km', ...
        ht_start, alt(end));
end

pharlap_root = ensure_pharlap(char(a.PharlapRoot));

iri_options = default_iri_options(topside);
iri_options = merge_options(iri_options, a.IriOptions);

if ~isfolder(outdir)
    [ok, msg] = mkdir(outdir);
    if ~ok, error('tss_iri_profiles:mkdir', 'could not create %s: %s', outdir, msg); end
end

tag = char(a.Tag);
if isempty(tag)
    tag = sprintf('%s_r12-%g', lower(strrep(strrep(topside, ' ', '-'), '.', '')), R12);
end
study = char(a.Study);
stem  = @(what, ext) fullfile(outdir, sprintf('%s_%s_%s.%s', study, what, tag, ext));

fprintf('tss_iri_profiles: %d points x %d altitudes (%.0f-%.0f km step %.0f)\n', ...
        n, na, ht_start, alt(end), ht_step);
fprintf('                  topside = %s, R12 = %g (F10.7 ~ %.1f)\n', ...
        topside, R12, 63.7 + 0.728*R12 + 0.00089*R12^2);

ne    = nan(n, na);
extra = nan(n, 10);
nfail = 0;
progress_every = max(round(a.Progress), 0);

t0 = tic;
for k = 1:n
    UT = [pts(k,1) pts(k,2) pts(k,3) pts(k,4) pts(k,5)];   % yr mo dy hh mm
    la = pts(k,6);
    lo = pts(k,7);

    try
        [iono, iono_extra] = iri2020(la, lo, R12, UT, ht_start, ht_step, na, ...
                                     iri_options);

        row = iono(1, 1:na);            % electron density, m^-3
        row(row < 0) = NaN;             % IRI's -1 sentinel is not a density
        ne(k,:) = row;

        e = iono_extra;
        extra(k,:) = [pick(e,1), pick(e,2), pick(e,23), pick(e,26), pick(e,27), ...
                      pick(e,48), pick(e,33), pick(e,41), pick(e,83), pick(e,51)];
    catch ME
        % One bad point must not lose a half-hour batch.  Record it and move on.
        nfail = nfail + 1;
        if nfail <= 5
            fprintf(2, '  point %d (lat %.1f lon %.1f) failed: %s\n', ...
                    k, la, lo, ME.message);
        end
    end

    if progress_every > 0 && mod(k, progress_every) == 0
        el = toc(t0);
        fprintf('  %7d / %7d   %5.1f s elapsed, %5.1f s remaining\n', ...
                k, n, el, el*(n-k)/k);
    end
end
elapsed = toc(t0);

% ---- write the CSV contract -------------------------------------------
info = struct();
info.ne_csv    = stem('profiles', 'csv');
info.extra_csv = stem('extra', 'csv');
info.alt_csv   = stem('alt_km', 'csv');
info.sidecar   = stem('profiles', 'json');

ne_names = arrayfun(@(i) sprintf('ne_m3_i%03d', i), 1:na, 'UniformOutput', false);
meta_names = {'point_id', 'year', 'month', 'day', 'hour_utc', 'minute_utc', ...
              'lat_deg', 'lon_deg'};

write_numeric_csv(info.ne_csv, [meta_names, ne_names], ...
    [(1:n).', pts], ne, '%.6e');

extra_names = {'nmf2_m3', 'hmf2_km', 'solar_zenith_deg', 'dip_lat_deg', ...
               'modip_deg', 'spread_f_prob', 'rz12', 'f107_daily', 'kp', 'ap'};
write_numeric_csv(info.extra_csv, [meta_names, extra_names], ...
    [(1:n).', pts], extra, '%.6g');

write_numeric_csv(info.alt_csv, {'alt_km'}, zeros(na, 0), alt(:), '%.6g');

info.n_points     = n;
info.n_alt        = na;
info.n_failed     = nfail;
info.elapsed_s    = elapsed;
info.pharlap_root = pharlap_root;
info.iri_options  = iri_options;
info.topside      = topside;
info.r12          = R12;
info.alt_km       = alt(:).';

% ---- sidecar, per docs/matlab-python-bridge.md ------------------------
side = struct();
side.producer      = 'matlab';
side.generator     = 'tss_iri_profiles';
side.tss_version   = tss_version();
side.created_utc   = char(datetime('now', 'TimeZone', 'UTC', ...
                                   'Format', 'uuuu-MM-dd''T''HH:mm:ss'));
side.study         = study;
side.tag           = tag;
side.topside       = topside;
side.r12           = R12;
side.n_points      = n;
side.n_alt         = na;
side.n_failed      = nfail;
side.alt_km        = alt(:).';
side.iri_options   = iri_options;
side.pharlap_root  = pharlap_root;
side.columns       = struct('profiles', {[meta_names, ne_names]}, ...
                            'extra',    {[meta_names, extra_names]}, ...
                            'alt_km',   {{'alt_km'}});
side.sentinel_note = ['IRI -1 means no valid value and is written as NaN.  ' ...
                      'A point that raised an error is an all-NaN row.'];
write_json(info.sidecar, side);

% ---- manifest ---------------------------------------------------------
info.manifest = '';
if a.Manifest
    cfg = struct('topside', topside, 'r12', R12, 'n_points', n, ...
                 'alt_km_start', alt(1), 'alt_km_step', ht_step, 'n_alt', na, ...
                 'iri_options', iri_options);
    [~, info.manifest] = tss_manifest(sprintf('%s_profiles', study), ...
        'Configuration', tag, ...
        'Config', cfg, ...
        'Tiers', struct('l0', 'matlab_iri'), ...
        'Outputs', {info.ne_csv, info.extra_csv, info.alt_csv, info.sidecar}, ...
        'Notes', sprintf('IRI-2020 through PHaRLAP at %s', pharlap_root), ...
        'Directory', outdir, ...
        'Filename', sprintf('%s_%s_manifest.json', study, tag));
end

fprintf('tss_iri_profiles: wrote %s in %.1f s', info.ne_csv, elapsed);
if nfail > 0
    fprintf(' (%d of %d points failed and are NaN)', nfail, n);
end
fprintf('\n');
end

% =======================================================================
function pts = load_points(src)
%LOAD_POINTS  Accept a matrix, a table or the author's headerless points.csv.
if istable(src)
    pts = table2array(src);
elseif ischar(src) || isstring(src)
    pts = readmatrix(char(src));
else
    pts = double(src);
end
if size(pts, 2) < 7
    error('tss_iri_profiles:points', ...
        ['points needs 7 columns [year month day hour_utc minute_utc ' ...
         'lat_deg lon_deg], got %d'], size(pts, 2));
end
pts = pts(:, 1:7);
end

% -----------------------------------------------------------------------
function alt = load_alt(src)
if ischar(src) || isstring(src)
    alt = readmatrix(char(src));
else
    alt = double(src);
end
alt = alt(:).';
end

% -----------------------------------------------------------------------
function root = ensure_pharlap(given)
%ENSURE_PHARLAP  Put PHaRLAP on the path and prove iri2020 is callable.
%
%   This repeats what startup_tss does, unconditionally and on purpose.  A
%   `matlab -batch` session inherits no interactive path and no interactive
%   environment, so a generator that assumes startup_tss already ran fails with
%   "Undefined function 'iri2020' for input arguments of type 'struct'" in
%   exactly the case where a long batch was being launched unattended.
root = given;
if isempty(root), root = getenv('TSS_PHARLAP_PATH'); end
if isempty(root), root = getenv('PHARLAP_ROOT'); end
if isempty(root)
    guesses = { ...
        'C:\Users\Andrew\PhD\Chartier Verification\pharlap_4.7.4', ...
        fullfile(getenv('USERPROFILE'), 'pharlap'), ...
        '/usr/local/pharlap'};
    for g = 1:numel(guesses)
        if isfolder(fullfile(guesses{g}, 'mex'))
            root = guesses{g};
            break
        end
    end
end
if isempty(root) || ~isfolder(fullfile(root, 'mex'))
    error('tss_iri_profiles:pharlap', ...
        ['Cannot locate PHaRLAP.  Pass PharlapRoot, or set TSS_PHARLAP_PATH, ' ...
         'to the folder containing mex/ and dat/.']);
end

addpath(fullfile(root, 'src', 'matlab'), fullfile(root, 'dat'));
addpath(fullfile(root, 'mex'));   % added last => highest precedence
setenv('DIR_MODELS_REF_DAT', fullfile(root, 'dat'));

% Verify by DOING, not by asking.  exist(name,'file') collapses MEX and .m into
% one return code, and PHaRLAP ships an iri2020.m help stub in src/matlab beside
% the mex, so a check for exist(...) == 3 rejects a working install.  One real
% call settles it, and also proves DIR_MODELS_REF_DAT is right, which no exist()
% check can do.
try
    [~, ~] = iri2020(0, 0, 100, [2015 1 1 12 0], 300, 10, 2, ...
                     struct('iri_messages', 'off'));
catch ME
    error('tss_iri_profiles:iri2020', ...
        ['iri2020 is present but not callable after adding %s to the path.\n' ...
         '  which(iri2020) = %s\n  DIR_MODELS_REF_DAT = %s\n  platform = %s\n' ...
         '  underlying error: %s'], ...
        fullfile(root, 'mex'), which('iri2020'), ...
        getenv('DIR_MODELS_REF_DAT'), computer, ME.message);
end
fprintf('tss_iri_profiles: PHaRLAP at %s (iri2020 verified callable)\n', root);
end

% -----------------------------------------------------------------------
function o = default_iri_options(topside)
%DEFAULT_IRI_OPTIONS  Everything that affects the answer, set explicitly.
%   A study whose result depends on an unstated default is not reproducible.
o = struct();
o.iri_messages     = 'off';
o.Ne_topside_model = topside;          % the study variable
o.foF2_coeffs      = 'URSI';
o.Ne_B0B1_model    = 'ABT-2009';
o.hmF2_model       = 'Shubin-COSMIC';
o.Ne_tops_limited  = 'f10.7 limited';  % clamps F10.7 to 188 for topside
o.Ne_profile_calc  = 'standard';
o.Ni_model         = 'RBV-2010 & TTS-2005';
o.Te_topside_model = 'TBT-2012';
o.D_model          = 'IRI-1990';
o.plasmasphere     = 'Ozhogin';
o.plasmapause      = 'without';
end

% -----------------------------------------------------------------------
function o = merge_options(o, extra)
if isempty(extra), return; end
if ischar(extra) || isstring(extra)
    if strlength(string(extra)) == 0, return; end
    extra = jsondecode(char(extra));
end
f = fieldnames(extra);
for k = 1:numel(f)
    o.(f{k}) = extra.(f{k});
end
end

% -----------------------------------------------------------------------
function v = pick(e, i)
%PICK  Read iono_extra(i), mapping IRI's -1 sentinel to NaN.
if numel(e) >= i && e(i) ~= -1
    v = e(i);
else
    v = NaN;
end
end

% -----------------------------------------------------------------------
function write_numeric_csv(path, names, meta, data, datafmt)
%WRITE_NUMERIC_CSV  Header row then one fprintf pass over the whole matrix.
%   meta is written with %g, data with datafmt.  Both may be empty.
fid = fopen(path, 'w');
if fid < 0
    error('tss_iri_profiles:write', 'could not open %s for writing', path);
end
cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>

fprintf(fid, '%s\n', strjoin(names, ','));

nmeta = size(meta, 2);
ndata = size(data, 2);
parts = cell(1, nmeta + ndata);
for k = 1:nmeta,  parts{k} = '%.10g'; end
for k = 1:ndata,  parts{nmeta + k} = datafmt; end
fmt = [strjoin(parts, ',') '\n'];

block = [meta, data];
fprintf(fid, fmt, block.');    % fprintf walks column major, so transpose
end

% -----------------------------------------------------------------------
function write_json(path, s)
fid = fopen(path, 'w');
if fid < 0
    error('tss_iri_profiles:write', 'could not open %s for writing', path);
end
cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
try
    txt = jsonencode(s, 'PrettyPrint', true);
catch
    txt = jsonencode(s);
end
fwrite(fid, txt, 'char');
end
