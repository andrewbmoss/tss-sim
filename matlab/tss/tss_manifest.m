function [manifest, outpath] = tss_manifest(study, varargin)
%   Layer:        support
%   WP:           WP0
%   Claims:       -
%   Status:       implemented
%   Requirements: -
%
%   Build, and optionally write, a run manifest matching the Python one.
%
%   Inputs
%     study          char, the study or generator name.  Required.
%     'Configuration'  char, the instrument or platform configuration scored.
%     'Config'         struct of resolved settings.  Hashed into the manifest.
%     'Tiers'          struct naming the tier implementation per layer.
%     'Seed'           numeric RNG seed, or [] for none.
%     'Inputs'         cellstr of paths, or a cell of {path, role} pairs, or a
%                      struct array with fields path and role.  Each existing
%                      file is checksummed.
%     'Outputs'        cellstr of paths written by the run.
%     'Notes'          char, free text.
%     'Directory'      where to write manifest.json.  Omit to build only.
%     'Filename'       default 'manifest.json'.
%
%   Outputs
%     manifest       struct with the same top level fields the Python
%                    tss.provenance.Manifest writes: study, configuration,
%                    config, tiers, seed, inputs, outputs, notes, created_utc,
%                    git, environment, config_sha256.
%     outpath        path written, or '' when Directory was not given.
%
%   Example
%     m = tss_manifest('link_budget', 'Config', cfg, ...
%                      'Outputs', {csvfile}, 'Directory', outdir);
%
%   WHERE THIS DIFFERS FROM THE PYTHON MANIFEST, AND IT DOES
%
%   The field names and the nesting are the same, so the two can be read by the
%   same consumer.  Two values will not match byte for byte and pretending
%   otherwise would be worse than saying so.
%
%   First, environment.  Python records a python version; this records a matlab
%   version and an environment.producer of 'matlab'.  A reader that wants to
%   know which side wrote a manifest should look there.
%
%   Second, config_sha256.  Python hashes json.dumps(config, sort_keys=True).
%   This hashes a canonical JSON built here with sorted keys, and MATLAB and
%   Python do not format every float identically, nor does MATLAB distinguish an
%   integer from a float the way Python does.  So the two hashes agree on simple
%   configurations and are not guaranteed to agree in general.  The field
%   config_sha256_method records which side computed it.  Do not use a hash
%   equality between a MATLAB manifest and a Python manifest as a test; compare
%   the config objects themselves.

p = inputParser();
p.FunctionName = 'tss_manifest';
p.addRequired('study', @(x) ischar(x) || isstring(x));
p.addParameter('Configuration', '', @(x) ischar(x) || isstring(x));
p.addParameter('Config', struct(), @isstruct);
p.addParameter('Tiers', struct(), @isstruct);
p.addParameter('Seed', [], @(x) isempty(x) || isnumeric(x));
p.addParameter('Inputs', {}, @(x) iscell(x) || isstruct(x));
p.addParameter('Outputs', {}, @(x) iscell(x) || ischar(x) || isstring(x));
p.addParameter('Notes', '', @(x) ischar(x) || isstring(x));
p.addParameter('Directory', '', @(x) ischar(x) || isstring(x));
p.addParameter('Filename', 'manifest.json', @(x) ischar(x) || isstring(x));
p.parse(study, varargin{:});
a = p.Results;

manifest = struct();
manifest.study         = char(a.study);
manifest.configuration = char(a.Configuration);
manifest.config        = a.Config;
manifest.tiers         = a.Tiers;

if isempty(a.Seed)
    manifest.seed = NaN;          % jsonencode writes NaN as null, matching None
else
    manifest.seed = double(a.Seed);
end

manifest.inputs      = build_inputs(a.Inputs);
manifest.outputs     = to_cellstr(a.Outputs);
manifest.notes       = char(a.Notes);
manifest.created_utc = utc_stamp();
manifest.git         = git_state(tss_repo_root());
manifest.environment = environment_record();

[manifest.config_sha256, manifest.config_sha256_method] = config_hash(a.Config);

outpath = '';
if ~isempty(char(a.Directory))
    d = char(a.Directory);
    if ~isfolder(d)
        [ok, msg] = mkdir(d);
        if ~ok
            error('tss_manifest:mkdir', 'could not create %s: %s', d, msg);
        end
    end
    outpath = fullfile(d, char(a.Filename));
    txt = encode_pretty(manifest);
    fid = fopen(outpath, 'w');
    if fid < 0
        error('tss_manifest:write', 'could not open %s for writing', outpath);
    end
    cleaner = onCleanup(@() fclose(fid));
    fwrite(fid, txt, 'char');
end
end

% -----------------------------------------------------------------------
function items = build_inputs(spec)
%BUILD_INPUTS  Normalise the Inputs argument and checksum each existing file.
items = {};
if isempty(spec), return; end

if isstruct(spec)
    for k = 1:numel(spec)
        role = '';
        if isfield(spec(k), 'role'), role = char(spec(k).role); end
        items{end+1} = one_input(char(spec(k).path), role); %#ok<AGROW>
    end
    return
end

for k = 1:numel(spec)
    e = spec{k};
    if iscell(e) && numel(e) >= 2
        items{end+1} = one_input(char(e{1}), char(e{2})); %#ok<AGROW>
    else
        items{end+1} = one_input(char(e), ''); %#ok<AGROW>
    end
end
end

% -----------------------------------------------------------------------
function rec = one_input(pth, role)
%ONE_INPUT  One inputs entry, with the same keys the Python side writes.
rec = struct('path', pth, 'role', role, 'sha256', '', 'bytes', '');
if isfile(pth)
    d = dir(pth);
    rec.bytes  = sprintf('%d', d(1).bytes);
    rec.sha256 = tss_sha256(pth);
end
end

% -----------------------------------------------------------------------
function c = to_cellstr(x)
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

% -----------------------------------------------------------------------
function s = utc_stamp()
%UTC_STAMP  ISO 8601 to the second, in UTC, matching Python's isoformat.
t = datetime('now', 'TimeZone', 'UTC');
s = sprintf('%04d-%02d-%02dT%02d:%02d:%02d+00:00', ...
            year(t), month(t), day(t), hour(t), minute(t), floor(second(t)));
end

% -----------------------------------------------------------------------
function g = git_state(repo)
%GIT_STATE  Commit, branch and dirty flag, or an honest note saying why not.
g = struct('commit', NaN, 'branch', NaN, 'dirty', NaN, 'dirty_files', {{}});

[st, commit] = run_git(repo, 'rev-parse HEAD');
if st ~= 0
    g.note = 'not a git working tree, or git is not on PATH';
    return
end
g.commit = strtrim(commit);

[~, branch] = run_git(repo, 'rev-parse --abbrev-ref HEAD');
g.branch = strtrim(branch);

[st, status] = run_git(repo, 'status --porcelain');
if st == 0
    lines = strsplit(strtrim(status), newline);
    lines = lines(~cellfun(@isempty, lines));
    g.dirty = ~isempty(lines);
    files = cell(1, numel(lines));
    for k = 1:numel(lines)
        L = lines{k};
        if numel(L) > 3, files{k} = L(4:end); else, files{k} = L; end
    end
    if numel(files) > 50, files = files(1:50); end
    g.dirty_files = files;
end
end

% -----------------------------------------------------------------------
function [status, out] = run_git(repo, args)
%RUN_GIT  Run git in repo without changing pwd.
cmd = sprintf('git -C "%s" %s', repo, args);
[status, out] = system(cmd);
end

% -----------------------------------------------------------------------
function e = environment_record()
%ENVIRONMENT_RECORD  Who and what produced this run.
e = struct();
e.producer = 'matlab';
e.matlab   = version();
e.platform = sprintf('%s %s', computer('arch'), computer());
e.host     = local_hostname();
e.user     = local_username();
end

% -----------------------------------------------------------------------
function h = local_hostname()
h = getenv('COMPUTERNAME');
if isempty(h), h = getenv('HOSTNAME'); end
if isempty(h)
    [st, out] = system('hostname');
    if st == 0, h = strtrim(out); else, h = 'unknown'; end
end
end

% -----------------------------------------------------------------------
function u = local_username()
u = getenv('USERNAME');
if isempty(u), u = getenv('USER'); end
if isempty(u), u = 'unknown'; end
end

% -----------------------------------------------------------------------
function [h, method] = config_hash(cfg)
%CONFIG_HASH  SHA-256 over a canonical JSON rendering with sorted keys.
method = 'matlab canonical json, sorted keys, sha256';
try
    blob = canonical_json(cfg);
    h = tss_sha256(unicode2native(blob, 'UTF-8'));
catch ME
    h = '';
    method = sprintf('unavailable: %s', ME.message);
end
end

% -----------------------------------------------------------------------
function s = canonical_json(v)
%CANONICAL_JSON  Deterministic JSON: object keys sorted, no insignificant space.
if isstruct(v) && numel(v) == 1
    f = sort(fieldnames(v));
    parts = cell(1, numel(f));
    for k = 1:numel(f)
        parts{k} = sprintf('%s: %s', jsonencode(f{k}), canonical_json(v.(f{k})));
    end
    s = ['{' strjoin(parts, ', ') '}'];
elseif isstruct(v)
    parts = arrayfun(@(x) canonical_json(x), v, 'UniformOutput', false);
    s = ['[' strjoin(parts, ', ') ']'];
elseif iscell(v)
    parts = cellfun(@(x) canonical_json(x), v(:).', 'UniformOutput', false);
    s = ['[' strjoin(parts, ', ') ']'];
else
    s = jsonencode(v);
end
end

% -----------------------------------------------------------------------
function txt = encode_pretty(s)
%ENCODE_PRETTY  Indented JSON where the release supports it, compact otherwise.
try
    txt = jsonencode(s, 'PrettyPrint', true);
catch
    txt = jsonencode(s);
end
end
