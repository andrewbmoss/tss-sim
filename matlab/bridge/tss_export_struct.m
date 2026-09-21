function info = tss_export_struct(s, matpath, varargin)
%   Layer:        bridge
%   WP:           WP0
%   Claims:       -
%   Status:       implemented
%   Requirements: -
%
%   Write a struct to .mat v7.3 with a JSON sidecar describing every field.
%
%   Inputs
%     s          scalar struct.  Each field becomes one top level variable in
%                the .mat file, which is what the Python reader expects.
%     matpath    output path.  A missing or different extension is replaced
%                with .mat.  The parent directory is created if absent.
%     'Producer'   default 'matlab'.
%     'Study'      study or work package name, recorded in the sidecar.
%     'Tag'        discriminator, recorded in the sidecar.
%     'Inputs'     cellstr of input paths, checksummed into the sidecar.
%     'Config'     struct of settings, hashed into the sidecar.
%     'Notes'      free text.
%
%   Outputs
%     info       struct with fields mat, sidecar, variables (the same records
%                written to the sidecar) and bytes.
%
%   Example
%     rays = struct('freq_hz', f, 'group_delay_s', tg, 'mode', 'O');
%     info = tss_export_struct(rays, ...
%              tss_repo_root('data','bridge','matlab','v2_rays_vertical.mat'), ...
%              'Study', 'raytrace_v2', 'Tag', 'vertical');
%
%   WHY A SIDECAR
%
%   A .mat v7.3 file is HDF5, so Python can open it with h5py and discover the
%   variable names.  What it cannot discover is what they mean: whether a 3-by-N
%   array is three components of N samples or N rows of a 3-vector, which
%   quantity is in which column, or what units anything is in.  The sidecar
%   carries the field names, classes, sizes and complexity so that the reader
%   checks rather than guesses, and so that a shape change shows up as a
%   mismatch against a recorded size instead of as a silent transpose.
%
%   THE TRANSPOSE, ONCE
%
%   HDF5 is row major and MATLAB is column major, so a v7.3 file presents each
%   array with its dimensions reversed relative to MATLAB's own view.  The
%   sidecar records size_matlab, which is the size as MATLAB sees it and as
%   SIZE reports it.  tss.bridge.read_mat on the Python side transposes back to
%   that.  Do not transpose a second time in the caller.
%
%   v7.3 is required and not merely preferred.  It is the only format that is
%   HDF5 underneath, so it is the only one h5py can read without a scipy.io
%   fallback, and it is the only one without the 2 GB limit.  A ray set from a
%   long run will exceed 2 GB.

p = inputParser();
p.FunctionName = 'tss_export_struct';
p.addRequired('s', @(x) isstruct(x) && isscalar(x));
p.addRequired('matpath', @(x) ischar(x) || isstring(x));
p.addParameter('Producer', 'matlab', @(x) ischar(x) || isstring(x));
p.addParameter('Study', '', @(x) ischar(x) || isstring(x));
p.addParameter('Tag', '', @(x) ischar(x) || isstring(x));
p.addParameter('Inputs', {}, @iscell);
p.addParameter('Config', struct(), @isstruct);
p.addParameter('Notes', '', @(x) ischar(x) || isstring(x));
p.parse(s, matpath, varargin{:});
a = p.Results;

matpath = char(matpath);
[outdir, base, ext] = fileparts(matpath);
if ~strcmpi(ext, '.mat')
    matpath = fullfile(outdir, [base '.mat']);
end
if ~isempty(outdir) && ~isfolder(outdir)
    [ok, msg] = mkdir(outdir);
    if ~ok
        error('tss_export_struct:mkdir', 'could not create %s: %s', outdir, msg);
    end
end

names = fieldnames(s);
if isempty(names)
    error('tss_export_struct:empty', ...
        'refusing to write an empty struct to %s', matpath);
end
bad = names(~cellfun(@isvarname, names));
if ~isempty(bad)
    error('tss_export_struct:fieldName', ...
        'field name is not a valid variable name and cannot be a .mat variable: %s', ...
        strjoin(bad, ', '));
end

save(matpath, '-struct', 's', '-v7.3');

% ---- describe every field ---------------------------------------------
vars = cell(1, numel(names));
for k = 1:numel(names)
    v = s.(names{k});
    rec = struct();
    rec.name         = names{k};
    rec.class        = class(v);
    rec.size_matlab  = double(size(v));
    rec.numel        = double(numel(v));
    rec.complex      = ~isreal(v);
    rec.sparse       = issparse(v);
    if ischar(v) && (isempty(v) || size(v, 1) <= 1)
        rec.value = v;               % short char fields are worth carrying
    end
    vars{k} = rec;
end

finfo = dir(matpath);

side = struct();
side.producer     = char(a.Producer);
side.generator    = 'tss_export_struct';
side.tss_version  = tss_version();
side.created_utc  = char(datetime('now', 'TimeZone', 'UTC', ...
                                  'Format', 'uuuu-MM-dd''T''HH:mm:ss'));
side.file         = [base '.mat'];
side.format       = 'mat-v7.3';
side.study        = char(a.Study);
side.tag          = char(a.Tag);
side.notes        = char(a.Notes);
side.bytes        = finfo(1).bytes;
side.struct_expanded = true;   % each field is a top level variable
side.orientation_note = ['size_matlab is the size MATLAB reports.  A v7.3 ' ...
    'file presents arrays with dimensions reversed; the Python reader ' ...
    'transposes back to size_matlab.  Do not transpose again.'];
side.variables = vars;

side.inputs = {};
for k = 1:numel(a.Inputs)
    pth = char(a.Inputs{k});
    rec = struct('path', pth, 'sha256', '', 'bytes', '');
    if isfile(pth)
        di = dir(pth);
        rec.bytes  = sprintf('%d', di(1).bytes);
        rec.sha256 = tss_sha256(pth);
    end
    side.inputs{end+1} = rec; %#ok<AGROW>
end

if ~isempty(fieldnames(a.Config))
    side.config = a.Config;
    try
        side.config_sha256 = tss_sha256( ...
            unicode2native(jsonencode(a.Config), 'UTF-8'));
    catch
        side.config_sha256 = '';
    end
end

sidecar = fullfile(outdir, [base '.json']);
fid = fopen(sidecar, 'w');
if fid < 0
    error('tss_export_struct:write', 'could not open %s for writing', sidecar);
end
cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
try
    txt = jsonencode(side, 'PrettyPrint', true);
catch
    txt = jsonencode(side);
end
fwrite(fid, txt, 'char');

info = struct('mat', matpath, 'sidecar', sidecar, ...
              'variables', {vars}, 'bytes', finfo(1).bytes);
end
