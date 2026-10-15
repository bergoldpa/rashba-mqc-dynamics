function P = get_case_parameters(start_dir)
% GET_CASE_PARAMETERS  Parameters of the test case of the current folder.
%
%   P = get_case_parameters() searches upwards from the current folder
%   (pwd) for the file case_parameters.m (located in the case folder, e.g.
%   case2/) and runs it. It returns the struct P defined there.
%
%   P = get_case_parameters(start_dir) does the same, starting from the
%   folder start_dir instead of the current folder (e.g. a case folder).
%
%   If the start folder also contains a local_parameters.m, that file is
%   run afterwards and can override single values (e.g. P.n_b = 5000).
%
%   The result is cached per folder, because the simulations call this
%   function very often; "clear all" resets the cache.

persistent cacheP cacheDir
if nargin < 1
    start_dir = pwd;
end
here = start_dir;
if ~isempty(cacheP) && strcmp(cacheDir, here)
    P = cacheP;
    return
end

% find the case folder
d = here;
while ~isfile(fullfile(d, 'case_parameters.m'))
    parent = fileparts(d);
    if strcmp(parent, d)
        error('get_case_parameters: no case_parameters.m found above %s.', here);
    end
    d = parent;
end

% read case parameters and optional local overrides
P = struct();
run(fullfile(d, 'case_parameters.m'));
if isfile(fullfile(here, 'local_parameters.m'))
    run(fullfile(here, 'local_parameters.m'));
end

cacheP = P;
cacheDir = here;
end
