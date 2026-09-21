function [t, info] = tss_read_csv_table(path, varargin)
%   Layer:        bridge
%   WP:           WP0
%   Claims:       -
%   Status:       implemented
%   Requirements: -
%
%   Read a CSV written by TSS_WRITE_CSV_TABLE, or by the Python side, as a
%   table, keeping the column names exactly as they appear in the file.
%
%   Inputs
%     path        path to the file.  The extension is forced to .csv, so a
%                 stem or a .json sidecar path is accepted and resolves to the
%                 same file the writer produced.
%     'AsStruct'  return a scalar struct of column vectors instead of a table.
%                 Default false.
%     'Columns'   cellstr of column names that must be present.  A missing one
%                 is an error naming every column that is missing and every
%                 column the file actually has.  Default {}, which checks
%                 nothing.
%     'Sidecar'   read the JSON sidecar beside the file when there is one, and
%                 cross-check its recorded row count and column list.
%                 Default true.
%
%   Outputs
%     t           table, or scalar struct when 'AsStruct' is true.
%     info        struct: csv, sidecar, has_sidecar, n_rows, n_cols, columns
%                 (the names as the header spells them), matlab_names (the
%                 names the table actually carries, which differ only if the
%                 release refused to preserve one), classes.
%
%   Example
%     [t, info] = tss_read_csv_table(fullfile(outdir, 'linkbudget_snr.csv'));
%     snr = t.snr_db;
%     s   = tss_read_csv_table(csv, 'AsStruct', true, ...
%                              'Columns', {'f_mhz', 'snr_db'});
%
%   THE COLUMN NAMING CONVENTION
%
%   This is the reading half of the convention documented in the header of
%   TSS_WRITE_CSV_TABLE, and that header is the authority.  In short: lower
%   snake case, the same identifier the code uses, with the SI unit as a
%   suffix.  f_mhz, snr_db, alt_km, ne_m3, sys_temp_k, noise_bw_hz.
%
%   The names are preserved and not made into valid MATLAB identifiers where
%   the release allows it.  READTABLE's default behaviour is to modify a header
%   it does not like, which would silently rename a column between writing and
%   reading and defeat the point of having one convention on both sides.
%   VariableNamingRule = 'preserve' is requested where it exists; on a release
%   that does not have it the header is read separately and put back, and a
%   name that cannot be a MATLAB variable name at all is left as the release
%   made it and reported in info.matlab_names.  Every name the convention
%   permits is a valid MATLAB variable name, so that case means the file does
%   not follow the convention.
%
%   A MISSING SIDECAR IS A WARNING AND NOT AN ERROR
%
%   docs/matlab-python-bridge.md says a bridge file without a sidecar is
%   treated as untrusted and a study that reads one warns.  That is what
%   happens here.  It is not an error because a CSV handed over by a
%   collaborator, or written by a spreadsheet, is still worth reading; it is a
%   warning because such a file carries no record of what produced it.
%
%   The sidecar is cross-checked and not trusted.  When it disagrees with the
%   file about the row count or the column list, the file wins and the
%   disagreement is a warning, because the file is the data and the sidecar is
%   a description of it that can go stale.
%
%   See also TSS_WRITE_CSV_TABLE, TSS_IMPORT_NPZ.

p = inputParser();
p.FunctionName = 'tss_read_csv_table';
p.addRequired('path', @(x) ischar(x) || isstring(x));
p.addParameter('AsStruct', false, @(x) islogical(x) || isnumeric(x));
p.addParameter('Columns', {}, @(x) iscell(x) || isstring(x) || ischar(x));
p.addParameter('Sidecar', true, @(x) islogical(x) || isnumeric(x));
p.parse(path, varargin{:});
a = p.Results;

path = char(a.path);
[outdir, base, ext] = fileparts(path);
if ~strcmpi(ext, '.csv')
    path = fullfile(outdir, [base '.csv']);
end
if ~isfile(path)
    error('tss_read_csv_table:missing', ...
        ['cannot read %s.  If a study should have produced it, run that ' ...
         'study first; the file is not created on demand.'], path);
end

header = read_header(path);
tbl    = read_table(path);

if numel(header) == width(tbl)
    keep = cellfun(@isvarname, header);
    if all(keep)
        tbl.Properties.VariableNames = header;
    else
        warning('tss_read_csv_table:naming', ...
            ['these header names are not valid MATLAB variable names and ' ...
             'could not be preserved: %s'], strjoin(header(~keep), ', '));
    end
else
    warning('tss_read_csv_table:headerWidth', ...
        ['%s has %d header names but readtable found %d columns.  The ' ...
         'names readtable chose are being used.'], ...
        path, numel(header), width(tbl));
    header = tbl.Properties.VariableNames;
end

tbl = strings_to_cellstr(tbl);
check_names(tbl.Properties.VariableNames);

info = struct();
info.csv          = path;
info.sidecar      = '';
info.has_sidecar  = false;
info.n_rows       = height(tbl);
info.n_cols       = width(tbl);
info.columns      = header(:).';
info.matlab_names = tbl.Properties.VariableNames;
info.classes      = cellfun(@(n) class(tbl.(n)), info.matlab_names, ...
                            'UniformOutput', false);

if a.Sidecar
    info = read_sidecar(info, outdir, base, tbl);
end

want = to_cellstr(a.Columns);
if ~isempty(want)
    missing = want(~ismember(want, info.matlab_names));
    if ~isempty(missing)
        error('tss_read_csv_table:columns', ...
            '%s is missing required column(s): %s.  It has: %s.', ...
            path, strjoin(missing, ', '), strjoin(info.matlab_names, ', '));
    end
end

if a.AsStruct
    t = struct();
    for k = 1:numel(info.matlab_names)
        t.(info.matlab_names{k}) = tbl.(info.matlab_names{k});
    end
else
    t = tbl;
end
end

% =======================================================================
function names = read_header(path)
%READ_HEADER  The first line, split on commas, as the file spells it.
%   Splitting on commas without a quote aware parser is correct here because
%   the convention forbids a space or a comma in a column name.  A file whose
%   header needs quoting does not follow the convention, and the width check in
%   the caller catches it.
fid = fopen(path, 'r');
if fid < 0
    error('tss_read_csv_table:open', 'could not open %s for reading', path);
end
cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
line = fgetl(fid);
if ~ischar(line)
    error('tss_read_csv_table:empty', '%s is empty', path);
end
parts = strsplit(line, ',', 'CollapseDelimiters', false);
names = cell(1, numel(parts));
for k = 1:numel(parts)
    names{k} = strtrim(strrep(parts{k}, '"', ''));
end
end

% =======================================================================
function tbl = read_table(path)
%READ_TABLE  READTABLE with name preservation where the release offers it.
try
    opts = detectImportOptions(path, 'FileType', 'text', 'Delimiter', ',');
catch
    opts = [];
end

if ~isempty(opts)
    try
        opts.VariableNamingRule = 'preserve';   % R2019b and later
    catch
        % Older release.  The header is put back by the caller instead.
    end
    tbl = readtable(path, opts);
else
    tbl = readtable(path, 'FileType', 'text', 'Delimiter', ',', ...
                    'ReadVariableNames', true);
end
end

% =======================================================================
function tbl = strings_to_cellstr(tbl)
%STRINGS_TO_CELLSTR  Text columns come back as cellstr on every release.
%   Which of string and cellstr READTABLE returns depends on the release and on
%   the import options.  Fixing it here means a caller comparing a status
%   column can use STRCMP without first asking what class it got.
names = tbl.Properties.VariableNames;
for k = 1:numel(names)
    v = tbl.(names{k});
    if isstring(v)
        tbl.(names{k}) = cellstr(v);
    end
end
end

% =======================================================================
function info = read_sidecar(info, outdir, base, tbl)
%READ_SIDECAR  Read the JSON sidecar and cross-check it against the file.
sidecar = fullfile(outdir, [base '.json']);
if ~isfile(sidecar)
    warning('tss_read_csv_table:noSidecar', ...
        ['%s has no sidecar %s.  Per docs/matlab-python-bridge.md a bridge ' ...
         'file without a sidecar carries no record of what produced it and ' ...
         'is treated as untrusted.'], info.csv, [base '.json']);
    return
end

try
    side = jsondecode(fileread(sidecar));
catch ME
    warning('tss_read_csv_table:sidecarUnreadable', ...
        'could not parse %s: %s', sidecar, ME.message);
    return
end

info.sidecar      = sidecar;
info.has_sidecar  = true;
info.sidecar_data = side;

if isfield(side, 'n_rows') && isnumeric(side.n_rows) && ...
        double(side.n_rows) ~= height(tbl)
    warning('tss_read_csv_table:rowCount', ...
        ['%s records %d rows but the file has %d.  The file wins; the ' ...
         'sidecar is stale.'], sidecar, double(side.n_rows), height(tbl));
end

if isfield(side, 'columns')
    recorded = to_cellstr(side.columns);
    actual   = info.matlab_names;
    same = numel(recorded) == numel(actual);
    if same
        same = all(strcmp(recorded(:).', actual(:).'));
    end
    if ~same
        warning('tss_read_csv_table:columnList', ...
            '%s records columns %s but the file has %s.  The file wins.', ...
            sidecar, strjoin(recorded, ', '), strjoin(actual, ', '));
    end
end
end

% =======================================================================
function c = to_cellstr(x)
if isempty(x)
    c = {};
elseif ischar(x)
    c = {x};
elseif isstring(x)
    c = cellstr(x(:).');
elseif iscell(x)
    c = cellfun(@char, x(:).', 'UniformOutput', false);
else
    c = {};
end
end

% =======================================================================
function check_names(names)
%CHECK_NAMES  Warn, do not fail, on a name that breaks the convention.
%   The same rule and the same wording as TSS_WRITE_CSV_TABLE, applied on the
%   way in so that a file from elsewhere is reported once rather than causing
%   confusion later.
bad = {};
for k = 1:numel(names)
    n = names{k};
    if ~strcmp(n, lower(n)) || contains(n, ' ')
        bad{end+1} = n; %#ok<AGROW>
    end
end
if ~isempty(bad)
    warning('tss_read_csv_table:naming', ...
        ['these column names are not lower snake case: %s.  The convention ' ...
         'is in the header of tss_write_csv_table.'], strjoin(bad, ', '));
end
end
