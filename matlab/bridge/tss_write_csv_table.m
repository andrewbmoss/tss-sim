function info = tss_write_csv_table(path, t, varargin)
%   Layer:        bridge
%   WP:           WP0
%   Claims:       -
%   Status:       implemented
%   Requirements: -
%
%   Write a table or a struct of columns as CSV, the always-works handoff.
%
%   Inputs
%     path       output path.  The extension is forced to .csv.  The parent
%                directory is created if absent.
%     t          a table, or a scalar struct whose fields are equal length
%                column vectors, or a numeric matrix with 'Names' supplied.
%     'Names'    column names, required when t is a matrix and otherwise used
%                to rename in order.
%     'Sidecar'  write a JSON sidecar beside the file.  Default true.
%     'Study'    study name, recorded in the sidecar.
%     'Tag'      discriminator, recorded in the sidecar.
%     'Notes'    free text, recorded in the sidecar.
%
%   Outputs
%     info       struct: csv, sidecar, n_rows, n_cols, columns.
%
%   Example
%     info = tss_write_csv_table(fullfile(outdir, 'linkbudget_snr.csv'), ...
%                                struct('f_mhz', f, 'snr_db', snr), ...
%                                'Study', 'linkbudget', 'Tag', 'baseline_6u');
%
%   THE COLUMN NAMING CONVENTION
%
%   One rule, applied without exception: a column is named with the same
%   identifier the code uses, in lower snake case, carrying its unit as a
%   suffix.  f_mhz, snr_db, alt_km, ne_m3, group_delay_s, lat_deg, sys_temp_k,
%   noise_bw_hz.  Never a bare alt or freq, never a capital, never a space,
%   never a unit in the header as a separate row or in parentheses.
%
%   The suffixes are the ones src/tss/units.py lists, which is the authority:
%     _m metres, _km kilometres, _s seconds, _hz hertz, _mhz megahertz,
%     _k kelvin, _w watts, _dbm dBm, _db decibels, _dbi dBi, _dbw dBW,
%     _deg degrees, _rad radians, _mps metres per second, _m3 per cubic metre.
%
%   The point is that a column in a CSV and a field in a Python dataclass are
%   obviously the same quantity, so a reader never has to work out whether
%   snr_db in one place is the same thing as snr in another.  A dimensionless
%   column takes no suffix.  A count takes n_, an identifier takes _id, and a
%   flag or status column holds a short lower case word and not a number.
%
%   Rows keep the order they were produced in.  Where there is a natural key,
%   for example frequency, the file is sorted by it ascending, so that two
%   files can be compared line by line.
%
%   Missing numeric values are written as NaN, which pandas.read_csv reads as
%   NaN without any converter argument.  An empty field would be read as NaN
%   too, but it is indistinguishable from a truncated line, so NaN is written
%   explicitly.
%
%   See also TSS_READ_CSV_TABLE.

p = inputParser();
p.FunctionName = 'tss_write_csv_table';
p.addRequired('path', @(x) ischar(x) || isstring(x));
p.addRequired('t', @(x) istable(x) || isstruct(x) || isnumeric(x));
p.addParameter('Names', {}, @(x) iscell(x) || isstring(x));
p.addParameter('Sidecar', true, @(x) islogical(x) || isnumeric(x));
p.addParameter('Study', '', @(x) ischar(x) || isstring(x));
p.addParameter('Tag', '', @(x) ischar(x) || isstring(x));
p.addParameter('Notes', '', @(x) ischar(x) || isstring(x));
p.parse(path, t, varargin{:});
a = p.Results;

path = char(a.path);
[outdir, base, ext] = fileparts(path);
if ~strcmpi(ext, '.csv')
    path = fullfile(outdir, [base '.csv']);
end
if ~isempty(outdir) && ~isfolder(outdir)
    [ok, msg] = mkdir(outdir);
    if ~ok
        error('tss_write_csv_table:mkdir', 'could not create %s: %s', outdir, msg);
    end
end

names = {};
if ~isempty(a.Names)
    if isstring(a.Names)
        names = cellstr(a.Names(:).');
    else
        names = cellfun(@char, a.Names(:).', 'UniformOutput', false);
    end
end

tbl = to_table(a.t, names);
check_names(tbl.Properties.VariableNames);

writetable(tbl, path, 'FileType', 'text', 'Delimiter', ',', ...
           'WriteVariableNames', true);

info = struct();
info.csv     = path;
info.n_rows  = height(tbl);
info.n_cols  = width(tbl);
info.columns = tbl.Properties.VariableNames;
info.sidecar = '';

if a.Sidecar
    side = struct();
    side.producer    = 'matlab';
    side.generator   = 'tss_write_csv_table';
    side.tss_version = tss_version();
    side.created_utc = char(datetime('now', 'TimeZone', 'UTC', ...
                                     'Format', 'uuuu-MM-dd''T''HH:mm:ss'));
    side.file        = [base '.csv'];
    side.format      = 'csv';
    side.study       = char(a.Study);
    side.tag         = char(a.Tag);
    side.notes       = char(a.Notes);
    side.n_rows      = height(tbl);
    side.columns     = info.columns;
    side.classes     = cellfun(@(n) class(tbl.(n)), info.columns, ...
                               'UniformOutput', false);
    side.convention  = ['lower snake case with an SI unit suffix; see the ' ...
                        'header comment of tss_write_csv_table'];
    info.sidecar = fullfile(outdir, [base '.json']);
    fid = fopen(info.sidecar, 'w');
    if fid < 0
        error('tss_write_csv_table:write', ...
            'could not open %s for writing', info.sidecar);
    end
    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
    try
        txt = jsonencode(side, 'PrettyPrint', true);
    catch
        txt = jsonencode(side);
    end
    fwrite(fid, txt, 'char');
end
end

% =======================================================================
function tbl = to_table(t, names)
%TO_TABLE  Accept a table, a struct of columns or a matrix plus names.
if istable(t)
    tbl = t;
    if ~isempty(names)
        if numel(names) ~= width(tbl)
            error('tss_write_csv_table:names', ...
                'Names has %d entries but the table has %d columns', ...
                numel(names), width(tbl));
        end
        tbl.Properties.VariableNames = names;
    end
    return
end

if isstruct(t)
    if ~isscalar(t)
        t = columns_from_struct_array(t);
    end
    f = fieldnames(t);
    if isempty(f)
        error('tss_write_csv_table:empty', 'nothing to write');
    end
    n = numel(t.(f{1}));
    cols = cell(1, numel(f));
    for k = 1:numel(f)
        v = t.(f{k});
        if numel(v) ~= n
            error('tss_write_csv_table:ragged', ...
                'field %s has %d elements, expected %d', f{k}, numel(v), n);
        end
        cols{k} = v(:);
    end
    tbl = table(cols{:}, 'VariableNames', f(:).');
    if ~isempty(names)
        tbl.Properties.VariableNames = names;
    end
    return
end

% numeric matrix
if isempty(names)
    error('tss_write_csv_table:names', ...
        'Names is required when writing a numeric matrix');
end
if numel(names) ~= size(t, 2)
    error('tss_write_csv_table:names', ...
        'Names has %d entries but the matrix has %d columns', ...
        numel(names), size(t, 2));
end
tbl = array2table(t, 'VariableNames', names);
end

% =======================================================================
function s = columns_from_struct_array(sa)
%COLUMNS_FROM_STRUCT_ARRAY  Turn a struct array of scalars into columns.
%   This is the shape tss_link_budget results arrive in when a study sweeps a
%   band: one struct per frequency.  Char fields such as status become a
%   cellstr column, which writetable quotes only when it has to.
f = fieldnames(sa);
s = struct();
for k = 1:numel(f)
    vals = {sa.(f{k})};
    if all(cellfun(@(v) isnumeric(v) && isscalar(v), vals))
        s.(f{k}) = cell2mat(vals(:));
    elseif all(cellfun(@(v) islogical(v) && isscalar(v), vals))
        s.(f{k}) = cell2mat(vals(:));
    else
        s.(f{k}) = cellfun(@char_or_empty, vals(:), 'UniformOutput', false);
    end
end
end

% =======================================================================
function c = char_or_empty(v)
if ischar(v)
    c = v;
elseif isstring(v) && isscalar(v)
    c = char(v);
elseif isnumeric(v) && isscalar(v)
    c = sprintf('%.17g', v);
else
    error('tss_write_csv_table:cellValue', ...
        'a column holds a %s, which cannot be written as a CSV cell', class(v));
end
end

% =======================================================================
function check_names(names)
%CHECK_NAMES  Warn, do not fail, on a name that breaks the convention.
%   A warning and not an error: a one-off comparison table is not worth
%   blocking, but a name that will end up in the dissertation should be
%   noticed while it is still cheap to change.
bad = {};
for k = 1:numel(names)
    n = names{k};
    if ~strcmp(n, lower(n)) || contains(n, ' ')
        bad{end+1} = n; %#ok<AGROW>
    end
end
if ~isempty(bad)
    warning('tss_write_csv_table:naming', ...
        ['these column names are not lower snake case: %s.  The convention ' ...
         'is in the header of tss_write_csv_table.'], strjoin(bad, ', '));
end
end
