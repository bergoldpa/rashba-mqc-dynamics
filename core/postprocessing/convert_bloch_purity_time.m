% CONVERT_BLOCH_PURITY_TIME  Convert the time axis of Bloch vector / purity data.
%
%   Reads <method>_n1.dat, _n2.dat, _n3.dat, _purity.dat (time in atomic
%   units) from the current folder and writes <method>_n1_N.dat ... with the
%   time in ps or ns (P.time_scale in case_parameters.m) for the LaTeX plots.
%
%   Usage:   set `method` below ('S' = SOFT, 'K' = koopmons, 'E' = Ehrenfest)
%            and run the script in the folder containing the .dat files.

% Clear the workspace, but keep run_case_file / run_case_method if this
% script is called from run_case.m
clearvars -global
clearvars -except run_case_file run_case_method
close all;

method = 'E';   % 'S' = SOFT, 'K' = koopmons, 'E' = Ehrenfest
if exist('run_case_method','var')
    method = run_case_method;    % set by run_case.m
end
% conversion factor from atomic time units to ps or ns
P = get_case_parameters();
if P.time_scale == "NANO"
    tfac = 2.4188843265857e-17*1e9;   % atomic time units -> ns
else
    tfac = 2.4188843265857e-17*1e12;  % atomic time units -> ps
end

% file prefix of the selected method
switch method
    case 'S'
        prefix = 'soft';
    case 'K'
        prefix = 'koopmons';
    case 'E'
        prefix = 'ehrenfest';
end

% read, convert the time column, write <prefix>_<quantity>_N.dat
for quantity = {'n1','n2','n3','purity'}
    A = readmatrix([prefix '_' quantity{1} '.dat'], 'FileType','text');
    out = [A(:,1)*tfac A(:,2)];
    writematrix(out, [prefix '_' quantity{1} '_N.dat'], 'Delimiter',' ', 'FileType','text');
end
