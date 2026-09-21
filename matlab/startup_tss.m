function info = startup_tss(pharlap_path)
%   Layer:        support
%   WP:           WP0
%   Claims:       -
%   Status:       implemented
%   Requirements: -
%
%   Put the tss_ functions, the bridge and the studies on the MATLAB path.
%
%   Inputs
%     pharlap_path   optional PHaRLAP installation root.  When omitted the
%                    search order is the TSS_PHARLAP_PATH environment
%                    variable, then PHARLAP_ROOT (the variable the author's
%                    iri_batch.m used), then the known install location, then
%                    a bounded search of the usual parents.
%
%   Outputs
%     info           struct describing what was found:
%                      repo_root, matlab_root, added {cellstr}
%                      pharlap_root, pharlap_ok, pharlap_reason
%                      iri2020_ok, raytrace_3d_on_path
%                      data_dir, bridge_dir, results_dir, version
%
%   Example
%     run('C:\path\to\tss-sim\matlab\startup_tss.m')
%     info = startup_tss;                       % once matlab/ is on the path
%     info = startup_tss('C:\opt\pharlap_4.7.4');
%
%   This function never calls RESTOREDEFAULTPATH.  The author's MATLAB tooling
%   depends on the interactive path, and a startup script that resets it breaks
%   more than it fixes.  Where a clean path is genuinely needed, save the
%   current one with PATH and restore it explicitly in a cleanup object.
%
%   It also does not change pwd, does not call CLEAR, and does not SAVEPATH.
%   The path it adds lasts for the session and no longer.

narginchk(0, 1);
if nargin < 1, pharlap_path = ''; end

here       = fileparts(mfilename('fullpath'));   % <repo>/matlab
repo_root  = fileparts(here);

info = struct();
info.repo_root   = repo_root;
info.matlab_root = here;
info.added       = {};

% ---- our own code ------------------------------------------------------
% Explicit list rather than genpath(here).  genpath would also pick up
% matlab/tests, and a test class on the path in every session is a good way to
% shadow something by accident.  tss_run_tests adds the tests folder itself.
for d = {'tss', 'bridge', 'studies'}
    p = fullfile(here, d{1});
    if isfolder(p)
        addpath(p);
        info.added{end+1} = p; %#ok<AGROW>
    else
        warning('startup_tss:missingFolder', 'expected folder is absent: %s', p);
    end
end

info.version = tss_version();

% ---- data and results directories, relative to the repo root -----------
info.data_dir    = fullfile(repo_root, 'data');
info.bridge_dir  = fullfile(repo_root, 'data', 'bridge', 'matlab');
info.results_dir = fullfile(repo_root, 'results');
for d = {info.bridge_dir, info.results_dir}
    if ~isfolder(d{1})
        [ok, msg] = mkdir(d{1});
        if ~ok
            warning('startup_tss:mkdir', 'could not create %s: %s', d{1}, msg);
        end
    end
end

% ---- PHaRLAP -----------------------------------------------------------
[info.pharlap_root, searched] = locate_pharlap(pharlap_path);

info.pharlap_ok          = false;
info.iri2020_ok          = false;
info.raytrace_3d_on_path = false;
info.pharlap_reason      = '';

if isempty(info.pharlap_root)
    info.pharlap_reason = sprintf( ...
        'no candidate contained a mex/ folder.  Looked at: %s', ...
        strjoin(searched, '; '));
else
    % Add the whole tree, then put mex first so the compiled entry points win
    % over the .m help stubs that ship beside them in src/matlab.
    addpath(genpath(info.pharlap_root));
    mexdir = fullfile(info.pharlap_root, 'mex');
    if isfolder(mexdir)
        addpath(mexdir, '-begin');
    end
    setenv('DIR_MODELS_REF_DAT', fullfile(info.pharlap_root, 'dat'));
    info.added{end+1} = info.pharlap_root;

    % WHY A TRIAL CALL AND NOT exist().
    %
    % An earlier version of this guard was exist('iri2020', 'file') and it
    % wrongly rejected a working PHaRLAP install.  Two reasons, and both are
    % worth keeping in front of whoever edits this next.
    %
    % First, exist() answers a question about the MATLAB path, not about the
    % disk.  Asking it before addpath tells you only that the path has not been
    % set up yet, which you already knew.  So the addpath has to come first,
    % unconditionally, and it does above.
    %
    % Second, even after the addpath, exist() collapses a MEX file and a .m
    % file into codes that do not distinguish the working case from the broken
    % one.  PHaRLAP ships an iri2020.m help stub in src/matlab alongside the
    % mex binary, so a test for exist(...) == 3 fails on a perfectly good
    % install where the stub happens to be found first, and a test for
    % exist(...) > 0 passes on an install whose mex never built.
    %
    % One real call settles both.  It also proves DIR_MODELS_REF_DAT points at
    % a dat/ folder that actually contains ig_rz.dat and apf107.dat, which no
    % exist() check can establish at all.  The arguments below are the same
    % ones the author's iri_batch.m uses for its own verification call, so this
    % is a known-good signature and not a guess.
    try
        [~, ~] = iri2020(0, 0, 100, [2015 1 1 12 0], 300, 10, 2, ...
                         struct('iri_messages', 'off'));
        info.iri2020_ok = true;
        info.pharlap_ok = true;
    catch ME
        info.pharlap_reason = sprintf( ...
            ['iri2020 is on the path but not callable.  which = %s, ' ...
             'DIR_MODELS_REF_DAT = %s, platform = %s, error = %s'], ...
            which('iri2020'), getenv('DIR_MODELS_REF_DAT'), computer, ME.message);
    end

    % raytrace_3d is reported as visible or not, and nothing stronger.  Its
    % argument list is not fixed here (see tss_raytrace_3d, which is a stub
    % until WP2), so there is no trial call that could be made honestly.  This
    % line answers a path question and says so.
    info.raytrace_3d_on_path = ~isempty(which('raytrace_3d'));
end

% ---- banner ------------------------------------------------------------
fprintf('\n  tss-sim MATLAB side, version %s\n', info.version);
fprintf('  repo root     %s\n', info.repo_root);
fprintf('  on the path   %s\n', strjoin(cellfun(@(p) local_rel(p, repo_root), ...
                                                info.added, 'UniformOutput', false), ', '));
fprintf('  bridge out    %s\n', local_rel(info.bridge_dir, repo_root));
fprintf('  results out   %s\n', local_rel(info.results_dir, repo_root));

if info.pharlap_ok
    fprintf('  PHaRLAP       %s\n', info.pharlap_root);
    fprintf('                iri2020 called successfully, so the IRI tiers will run.\n');
    if info.raytrace_3d_on_path
        fprintf('                raytrace_3d is visible on the path (not called, see tss_raytrace_3d).\n');
    else
        fprintf('                raytrace_3d is NOT visible; the ray-traced L1 tier will not run.\n');
    end
else
    fprintf('  PHaRLAP       NOT USABLE\n');
    fprintf('                %s\n', wrap_reason(info.pharlap_reason));
    fprintf('                Set TSS_PHARLAP_PATH or pass the path to startup_tss.\n');
    fprintf('                Everything that does not need IRI or ray tracing still works:\n');
    fprintf('                the link budget, the bridge, and their tests.\n');
end
fprintf('\n');
end

% -----------------------------------------------------------------------
function [root, searched] = locate_pharlap(given)
%LOCATE_PHARLAP  Env var, then the known install, then a bounded search.
searched = {};
candidates = {};

if ~isempty(given)
    candidates{end+1} = given;
end
for v = {'TSS_PHARLAP_PATH', 'PHARLAP_ROOT'}
    e = getenv(v{1});
    if ~isempty(e), candidates{end+1} = e; end %#ok<AGROW>
end

% The author's install.  Kept as a literal because it is the machine the
% dissertation runs on, not because anything should depend on it.
candidates{end+1} = 'C:\Users\Andrew\PhD\Chartier Verification\pharlap_4.7.4';
candidates{end+1} = fullfile(getenv('USERPROFILE'), 'pharlap');
candidates{end+1} = '/usr/local/pharlap';

% Bounded search: one level of pharlap* under a short list of parents.  This
% is deliberately not a recursive scan of the whole drive.
parents = {getenv('USERPROFILE'), getenv('HOME'), ...
           fullfile(getenv('USERPROFILE'), 'PhD'), ...
           'C:\opt', '/opt', '/usr/local'};
for k = 1:numel(parents)
    if isempty(parents{k}) || ~isfolder(parents{k}), continue; end
    d = dir(fullfile(parents{k}, 'pharlap*'));
    for j = 1:numel(d)
        if d(j).isdir
            candidates{end+1} = fullfile(d(j).folder, d(j).name); %#ok<AGROW>
        end
    end
end

root = '';
for k = 1:numel(candidates)
    cand = candidates{k};
    if isempty(cand), continue; end
    searched{end+1} = cand; %#ok<AGROW>
    if isfolder(fullfile(cand, 'mex'))
        root = cand;
        return
    end
end
end

% -----------------------------------------------------------------------
function s = local_rel(p, root)
%LOCAL_REL  Shorten a path for the banner when it sits under the repo root.
if strncmp(p, root, numel(root))
    s = p(numel(root)+2:end);
    if isempty(s), s = '.'; end
else
    s = p;
end
end

% -----------------------------------------------------------------------
function s = wrap_reason(msg)
%WRAP_REASON  Keep the banner to one screen width.
if numel(msg) > 180
    s = [msg(1:177) '...'];
else
    s = msg;
end
end
