function [data, info] = tss_import_npz(npzpath, varargin)
%   Layer:        bridge
%   WP:           WP0
%   Claims:       -
%   Status:       implemented
%   Requirements: -
%
%   Read a NumPy .npz written by the Python side.
%
%   Inputs
%     npzpath      path to the .npz file.
%     'Names'      cellstr of array names to read.  Default: all of them.
%     'TempDir'    where to unzip.  Default: a fresh tempname directory,
%                  deleted when the function returns.
%     'Fallback'   when true (the default), a .npz that cannot be opened or
%                  parsed falls back to a sibling .mat or .csv of the same
%                  stem, with a warning.  Set false to make a broken .npz an
%                  error.
%
%   Outputs
%     data         struct, one field per array.  A name that is not a valid
%                  MATLAB field name is made valid and the original is kept in
%                  info.names.
%     info         struct: source ('npz', 'mat' or 'csv'), file, names (a
%                  2-column cell of original name and field name), dtypes,
%                  shapes, fortran_order flags.
%
%   Example
%     d = tss_import_npz(tss_repo_root('tests','fixtures','bridge', ...
%                                      'link_budget_inputs.npz'));
%     f = d.freq_hz;
%
%   WHICH FILE IS AUTHORITATIVE
%
%   The .npz is authoritative.  It is what the Python side writes and what the
%   checked-in fixtures are.  The fallback to a sibling .mat or .csv exists
%   only for a session where UNZIP is unavailable or a file is truncated, and
%   it warns every time it fires so that it cannot be mistaken for the normal
%   path.  If the fallback file and the .npz ever disagree, the .npz is right
%   and the other one is stale; regenerate it rather than reconciling by hand.
%
%   WHAT THIS READER HANDLES
%
%   The subset of the .npy format this repository uses, which is the subset
%   docs/matlab-python-bridge.md commits the Python side to writing:
%
%     handled      numeric arrays of any shape, real and complex, in any of
%                  bool, int8 to int64, uint8 to uint64, float32, float64,
%                  complex64 and complex128; C order and Fortran order; .npy
%                  format versions 1.0 and 2.0; little and big endian.
%     not handled  object arrays, pickled objects, structured and record
%                  dtypes, datetime64 and timedelta64, float16.  These raise a
%                  clear error naming the dtype rather than returning
%                  something plausible.
%
%   Always save with numpy.savez and explicit keyword names, never positional.
%   Positional names come out as arr_0, arr_1 and so on, which tells a reader
%   nothing and breaks the moment an argument is inserted.

p = inputParser();
p.FunctionName = 'tss_import_npz';
p.addRequired('npzpath', @(x) ischar(x) || isstring(x));
p.addParameter('Names', {}, @(x) iscell(x) || ischar(x) || isstring(x));
p.addParameter('TempDir', '', @(x) ischar(x) || isstring(x));
p.addParameter('Fallback', true, @(x) islogical(x) || isnumeric(x));
p.parse(npzpath, varargin{:});
a = p.Results;

npzpath = char(a.npzpath);
wanted  = normalise_names(a.Names);

data = struct();
info = struct('source', '', 'file', npzpath, 'names', {cell(0, 2)}, ...
              'dtypes', {{}}, 'shapes', {{}}, 'fortran_order', {{}});

if ~isfile(npzpath)
    if a.Fallback
        [data, info] = try_fallback(npzpath, 'the .npz file does not exist');
        return
    end
    error('tss_import_npz:noFile', 'no such file: %s', npzpath);
end

try
    [data, info] = read_npz(npzpath, wanted, char(a.TempDir));
catch ME
    if a.Fallback
        [data, info] = try_fallback(npzpath, ME.message);
        return
    end
    rethrow(ME);
end
end

% =======================================================================
function [data, info] = read_npz(npzpath, wanted, tmpdir)
%READ_NPZ  Unzip to a temporary directory and parse each member as .npy.
%   A .npz is an ordinary zip archive of .npy members, so UNZIP opens it.

owned = isempty(tmpdir);
if owned
    tmpdir = [tempname() '_npz'];
end
if ~isfolder(tmpdir)
    [ok, msg] = mkdir(tmpdir);
    if ~ok
        error('tss_import_npz:mkdir', 'could not create %s: %s', tmpdir, msg);
    end
end
if owned
    cleaner = onCleanup(@() rmdir(tmpdir, 's')); %#ok<NASGU>
end

files = unzip(npzpath, tmpdir);
if isempty(files)
    error('tss_import_npz:empty', '%s contained no members', npzpath);
end

data = struct();
names   = cell(0, 2);
dtypes  = {};
shapes  = {};
forder  = {};

for k = 1:numel(files)
    [~, stem, ext] = fileparts(files{k});
    if ~strcmpi(ext, '.npy')
        continue                      % savez writes nothing else
    end
    if ~isempty(wanted) && ~any(strcmp(stem, wanted))
        continue
    end
    [v, meta] = read_npy(files{k});
    field = matlab.lang.makeValidName(stem);
    data.(field)  = v;
    names(end+1, :) = {stem, field};      %#ok<AGROW>
    dtypes{end+1}   = meta.descr;         %#ok<AGROW>
    shapes{end+1}   = meta.shape;         %#ok<AGROW>
    forder{end+1}   = meta.fortran_order; %#ok<AGROW>
end

if ~isempty(wanted)
    missing = setdiff(wanted, names(:, 1));
    if ~isempty(missing)
        error('tss_import_npz:missingName', ...
            '%s does not contain: %s', npzpath, strjoin(missing, ', '));
    end
end

info = struct('source', 'npz', 'file', npzpath, 'names', {names}, ...
              'dtypes', {dtypes}, 'shapes', {shapes}, 'fortran_order', {forder});
end

% =======================================================================
function [v, meta] = read_npy(pth)
%READ_NPY  Parse one .npy member.

fid = fopen(pth, 'r', 'l');
if fid < 0
    error('tss_import_npz:open', 'could not open %s', pth);
end
closer = onCleanup(@() fclose_if_open(fid));

magic = fread(fid, 6, 'uint8=>uint8').';
if numel(magic) < 6 || ~isequal(double(magic), [147 78 85 77 80 89])
    error('tss_import_npz:magic', '%s is not a .npy file', pth);
end
major = fread(fid, 1, 'uint8');
minor = fread(fid, 1, 'uint8'); %#ok<NASGU>
if major == 1
    hlen = fread(fid, 1, 'uint16');
    offset = 10 + hlen;
elseif major == 2 || major == 3
    hlen = fread(fid, 1, 'uint32');
    offset = 12 + hlen;
else
    error('tss_import_npz:version', ...
        '%s is .npy format version %d, which this reader does not handle', ...
        pth, major);
end
hdr = char(fread(fid, hlen, 'uint8=>uint8').');
clear closer                    % the cleanup object closes fid

meta = parse_header(hdr, pth);

[cls, bytes_each, is_complex] = dtype_to_class(meta.descr, pth);
machinefmt = 'l';
if meta.descr(1) == '>'
    machinefmt = 'b';
end

n = prod(max(meta.shape, 0));
if isempty(meta.shape)
    n = 1;                                   % a 0-d array is one scalar
end
nread = n;
if is_complex
    nread = 2 * n;
end

fid = fopen(pth, 'r', machinefmt);
if fid < 0
    error('tss_import_npz:open', 'could not open %s', pth);
end
closer2 = onCleanup(@() fclose_if_open(fid)); %#ok<NASGU>
if fseek(fid, offset, 'bof') ~= 0
    error('tss_import_npz:seek', 'could not seek past the header of %s', pth);
end

raw = fread(fid, nread, sprintf('%s=>%s', cls, cls));
if numel(raw) ~= nread
    error('tss_import_npz:truncated', ...
        '%s holds %d values of %s but the header declares %d', ...
        pth, numel(raw), cls, nread);
end

if is_complex
    v = complex(raw(1:2:end), raw(2:2:end));
else
    v = raw;
end
if strcmp(meta.descr(2:end), 'b1')
    v = logical(v);
end

v = reshape_like_numpy(v, meta.shape, meta.fortran_order);

meta.bytes_each = bytes_each;
end

% =======================================================================
function meta = parse_header(hdr, pth)
%PARSE_HEADER  Pull descr, fortran_order and shape out of the Python dict.
descr = regexp(hdr, '''descr''\s*:\s*''([^'']+)''', 'tokens', 'once');
if isempty(descr)
    error('tss_import_npz:header', ...
        ['could not read a simple descr from the header of %s.  A structured ' ...
         'or record dtype is not handled; see the WHAT THIS READER HANDLES ' ...
         'block in tss_import_npz.'], pth);
end
meta.descr = descr{1};

meta.fortran_order = ~isempty(regexp(hdr, '''fortran_order''\s*:\s*True', 'once'));

shape_txt = regexp(hdr, '''shape''\s*:\s*\(([^)]*)\)', 'tokens', 'once');
if isempty(shape_txt)
    error('tss_import_npz:header', 'could not read a shape from %s', pth);
end
parts = strtrim(strsplit(strtrim(shape_txt{1}), ','));
parts = parts(~cellfun(@isempty, parts));
meta.shape = cellfun(@str2double, parts);
if any(isnan(meta.shape))
    error('tss_import_npz:header', 'unreadable shape in %s: (%s)', ...
        pth, shape_txt{1});
end
end

% =======================================================================
function [cls, nbytes, is_complex] = dtype_to_class(descr, pth)
%DTYPE_TO_CLASS  Map a NumPy dtype string onto a MATLAB fread class.
is_complex = false;
kind = descr(2:end);
switch kind
    case 'b1', cls = 'uint8';   nbytes = 1;
    case 'i1', cls = 'int8';    nbytes = 1;
    case 'i2', cls = 'int16';   nbytes = 2;
    case 'i4', cls = 'int32';   nbytes = 4;
    case 'i8', cls = 'int64';   nbytes = 8;
    case 'u1', cls = 'uint8';   nbytes = 1;
    case 'u2', cls = 'uint16';  nbytes = 2;
    case 'u4', cls = 'uint32';  nbytes = 4;
    case 'u8', cls = 'uint64';  nbytes = 8;
    case 'f4', cls = 'single';  nbytes = 4;
    case 'f8', cls = 'double';  nbytes = 8;
    case 'c8', cls = 'single';  nbytes = 8;  is_complex = true;
    case 'c16', cls = 'double'; nbytes = 16; is_complex = true;
    otherwise
        error('tss_import_npz:dtype', ...
            ['dtype ''%s'' in %s is not handled.  Object arrays, pickled ' ...
             'objects, structured dtypes, datetime64 and float16 are out of ' ...
             'scope on purpose; see the header comment.  Write the array as ' ...
             'a plain numeric dtype on the Python side.'], descr, pth);
end
end

% =======================================================================
function v = reshape_like_numpy(v, shape, fortran_order)
%RESHAPE_LIKE_NUMPY  Undo C ordering so MATLAB sees the array NumPy saw.
%   NumPy defaults to C order, which is row major.  MATLAB is column major.
%   Reshaping to the reversed shape and then reversing the dimensions puts the
%   elements where MATLAB expects them.  A Fortran-ordered .npy needs neither
%   step.
if isempty(shape)
    return                                  % 0-d: already a scalar
end
if numel(shape) == 1
    v = reshape(v, [shape(1), 1]);          % 1-d comes back as a column
    return
end
if fortran_order
    v = reshape(v, shape);
else
    v = permute(reshape(v, fliplr(shape)), numel(shape):-1:1);
end
end

% =======================================================================
function [data, info] = try_fallback(npzpath, why)
%TRY_FALLBACK  Read the sibling .mat or .csv, loudly.
[d, stem] = fileparts(npzpath);
matpath = fullfile(d, [stem '.mat']);
csvpath = fullfile(d, [stem '.csv']);

if isfile(matpath)
    warning('tss_import_npz:fallbackMat', ...
        ['could not use %s (%s).  Falling back to %s.  The .npz is the ' ...
         'authoritative file; if the two disagree, regenerate the .mat.'], ...
        npzpath, why, matpath);
    data = load(matpath);
    info = struct('source', 'mat', 'file', matpath, 'names', {cell(0, 2)}, ...
                  'dtypes', {{}}, 'shapes', {{}}, 'fortran_order', {{}});
    f = fieldnames(data);
    names = cell(numel(f), 2);
    for k = 1:numel(f)
        names(k, :) = {f{k}, f{k}};
    end
    info.names = names;
    return
end

if isfile(csvpath)
    warning('tss_import_npz:fallbackCsv', ...
        ['could not use %s (%s).  Falling back to %s.  The .npz is the ' ...
         'authoritative file; if the two disagree, regenerate the .csv.'], ...
        npzpath, why, csvpath);
    t = tss_read_csv_table(csvpath);
    data = struct();
    names = cell(width(t), 2);
    vn = t.Properties.VariableNames;
    for k = 1:numel(vn)
        field = matlab.lang.makeValidName(vn{k});
        data.(field) = t.(vn{k});
        names(k, :) = {vn{k}, field};
    end
    info = struct('source', 'csv', 'file', csvpath, 'names', {names}, ...
                  'dtypes', {{}}, 'shapes', {{}}, 'fortran_order', {{}});
    return
end

error('tss_import_npz:noFallback', ...
    ['could not use %s (%s), and neither %s nor %s exists.  Regenerate the ' ...
     'bridge file on the Python side.'], npzpath, why, matpath, csvpath);
end

% =======================================================================
function c = normalise_names(x)
%NORMALISE_NAMES  Accept a char, a string array or a cellstr.
if isempty(x)
    c = {};
elseif ischar(x)
    c = {x};
elseif isstring(x)
    c = cellstr(x(:).');
else
    c = cellfun(@char, x(:).', 'UniformOutput', false);
end
end

% =======================================================================
function fclose_if_open(fid)
%FCLOSE_IF_OPEN  Close fid, tolerating a fid that is already closed.
try
    if ~isempty(fopen(fid))
        fclose(fid);
    end
catch
    % already closed, or never opened; nothing to do
end
end
