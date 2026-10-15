% SETUP  Add the project code to the MATLAB path.
%
%   Run once from anywhere (e.g. "run /path/to/project/setup.m" or open
%   the project folder and type "setup"). It adds core/ and all its
%   subfolders to the MATLAB path for the current session. To keep the path
%   for future sessions, run "savepath" afterwards.

project_root = fileparts(mfilename('fullpath'));
addpath(genpath(fullfile(project_root, 'core')));
fprintf('Added %s to the MATLAB path.\n', fullfile(project_root, 'core'));
clear project_root
