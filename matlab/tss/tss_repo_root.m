function root = tss_repo_root(varargin)
%   Layer:        support
%   WP:           WP0
%   Claims:       -
%   Status:       implemented
%   Requirements: -
%
%   Locate the tss-sim repository root from this file, never from pwd.
%
%   Inputs
%     varargin   optional path parts appended to the root, in the style of
%                FULLFILE.  tss_repo_root('data','bridge','matlab') returns
%                the bridge output directory.
%
%   Outputs
%     root       absolute path.
%
%   Example
%     r   = tss_repo_root();
%     cfg = tss_repo_root('config', 'default.yaml');
%
%   pwd is not used.  A function that locates the repository from the working
%   directory works when it is called from the root and fails when it is called
%   from anywhere else, which is the worst kind of failure because it depends on
%   how the caller happened to be sitting.  This file lives at
%   <repo>/matlab/tss/, so three fileparts up is the root and that holds
%   wherever it is called from.

here = fileparts(mfilename('fullpath'));    % <repo>/matlab/tss
root = fileparts(fileparts(here));          % <repo>

if ~isempty(varargin)
    root = fullfile(root, varargin{:});
end
end
