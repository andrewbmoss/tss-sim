function h = tss_sha256(x)
%   Layer:        support
%   WP:           WP0
%   Claims:       -
%   Status:       implemented
%   Requirements: -
%
%   SHA-256 of a file or of a byte vector, as lower case hex.
%
%   Inputs
%     x   a path to an existing file (char or string), or a uint8 vector of
%         bytes.  A char vector is always treated as a path, never as text;
%         hash text by converting it first with unicode2native(txt,'UTF-8').
%
%   Outputs
%     h   64 character lower case hex digest, or '' if no digest engine is
%         available and the failure was non-fatal.
%
%   Example
%     h = tss_sha256(fullfile(tss_repo_root(), 'README.md'));
%     h = tss_sha256(unicode2native('hello', 'UTF-8'));
%
%   The digest comes from the JVM, which is present in every normal MATLAB
%   session.  When MATLAB was started with -nojvm the function falls back to
%   sha256sum or certutil on the system path.  If neither is available it
%   errors rather than returning a wrong or empty digest silently, because a
%   manifest that records a blank checksum beside a real path invites the
%   reader to assume the file was verified.

if isstring(x), x = char(x); end

if ischar(x)
    if ~isfile(x)
        error('tss_sha256:noFile', 'not a file: %s', x);
    end
    h = digest_file(x);
elseif isa(x, 'uint8')
    h = digest_bytes(x(:));
else
    error('tss_sha256:type', ...
        'expected a file path or a uint8 vector, got %s', class(x));
end
end

% -----------------------------------------------------------------------
function h = digest_file(pth)
if usejava('jvm')
    md  = java.security.MessageDigest.getInstance('SHA-256');
    fid = fopen(pth, 'r');
    if fid < 0
        error('tss_sha256:open', 'could not open %s', pth);
    end
    cleaner = onCleanup(@() fclose(fid)); %#ok<NASGU>
    while true
        block = fread(fid, 1048576, '*uint8');
        if isempty(block), break; end
        md.update(typecast(block, 'int8'));
    end
    h = hexify(typecast(md.digest(), 'uint8'));
else
    h = digest_by_shell(pth);
end
end

% -----------------------------------------------------------------------
function h = digest_bytes(bytes)
if usejava('jvm')
    md = java.security.MessageDigest.getInstance('SHA-256');
    if ~isempty(bytes)
        md.update(typecast(bytes, 'int8'));
    end
    h = hexify(typecast(md.digest(), 'uint8'));
else
    tmp = [tempname() '.bin'];
    fid = fopen(tmp, 'w');
    fwrite(fid, bytes, 'uint8');
    fclose(fid);
    c = onCleanup(@() delete(tmp)); %#ok<NASGU>
    h = digest_by_shell(tmp);
end
end

% -----------------------------------------------------------------------
function h = digest_by_shell(pth)
%DIGEST_BY_SHELL  Last resort for a -nojvm session.
if ispc
    [st, out] = system(sprintf('certutil -hashfile "%s" SHA256', pth));
    if st == 0
        lines = strsplit(strtrim(out), newline);
        for k = 1:numel(lines)
            t = lower(strtrim(strrep(lines{k}, ' ', '')));
            if numel(t) == 64 && all(isstrprop(t, 'xdigit'))
                h = t;
                return
            end
        end
    end
else
    [st, out] = system(sprintf('sha256sum "%s"', pth));
    if st == 0
        tok = regexp(strtrim(out), '^([0-9a-f]{64})', 'tokens', 'once');
        if ~isempty(tok)
            h = tok{1};
            return
        end
    end
    [st, out] = system(sprintf('shasum -a 256 "%s"', pth));
    if st == 0
        tok = regexp(strtrim(out), '^([0-9a-f]{64})', 'tokens', 'once');
        if ~isempty(tok)
            h = tok{1};
            return
        end
    end
end
error('tss_sha256:noEngine', ...
    ['no SHA-256 engine: MATLAB has no JVM and neither sha256sum nor ' ...
     'certutil answered.  Cannot checksum %s.'], pth);
end

% -----------------------------------------------------------------------
function h = hexify(bytes)
h = lower(reshape(sprintf('%02x', bytes), 1, []));
end
