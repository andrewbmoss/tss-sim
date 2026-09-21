function v = tss_version()
%   Layer:        support
%   WP:           WP0
%   Claims:       -
%   Status:       implemented
%   Requirements: -
%
%   Return the tss package version string.
%
%   Inputs
%     (none)
%
%   Outputs
%     v   character vector, for example '0.1.0'.
%
%   Example
%     fprintf('tss %s\n', tss_version());
%
%   The version is read from src/tss/__init__.py rather than held here, so
%   there is one place to bump it and the two language halves cannot drift
%   apart.  If that file cannot be read the function returns '0.0.0+unknown'
%   and warns once, because a manifest that records a made-up version is worse
%   than one that records an honest unknown.

persistent cached
if ~isempty(cached)
    v = cached;
    return
end

v = '0.0.0+unknown';
initfile = tss_repo_root('src', 'tss', '__init__.py');

if isfile(initfile)
    try
        txt = fileread(initfile);
        tok = regexp(txt, '__version__\s*=\s*["'']([^"'']+)["'']', 'tokens', 'once');
        if ~isempty(tok)
            v = tok{1};
        else
            warning('tss_version:noMatch', ...
                'no __version__ assignment found in %s', initfile);
        end
    catch ME
        warning('tss_version:unreadable', ...
            'could not read %s: %s', initfile, ME.message);
    end
else
    warning('tss_version:missing', ...
        'cannot find %s, reporting %s', initfile, v);
end

cached = v;
end
