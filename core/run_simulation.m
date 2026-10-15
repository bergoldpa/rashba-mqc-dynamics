function run_simulation(n, method)
% RUN_SIMULATION  Run one simulation of test case n in its own folder.
%
%   run_simulation(n, method) changes to the run folder of caseN, runs the
%   simulation there (results are written to that folder) and returns to
%   the previous folder. method is one of
%     "soft"             caseN/soft               soft_simulation
%     "ehrenfest"        caseN/ehrenfest/n500     ehrenfest_simulation
%     "ehrenfest_n5000"  caseN/ehrenfest/n5000    ehrenfest_simulation (case 3)
%     "koopmons"         caseN/koopmons/n500_a05  koopmons_simulation
%     "compare_n500_vs_n5000"  caseN/ehrenfest/n500_vs_n5000
%                              compare_n500_vs_n5000 (case 3, after both runs)

switch string(method)
    case "soft"
        sub = 'soft';                   script = 'soft_simulation';
    case "ehrenfest"
        sub = fullfile('ehrenfest','n500');     script = 'ehrenfest_simulation';
    case "ehrenfest_n5000"
        sub = fullfile('ehrenfest','n5000');    script = 'ehrenfest_simulation';
    case "koopmons"
        sub = fullfile('koopmons','n500_a05');  script = 'koopmons_simulation';
    case "compare_n500_vs_n5000"
        sub = fullfile('ehrenfest','n500_vs_n5000'); script = 'compare_n500_vs_n5000';
    otherwise
        error('run_simulation: unknown method "%s".', method);
end

root = fileparts(fileparts(mfilename('fullpath')));   % CODE folder
folder = fullfile(root, sprintf('case%d', n), sub);
if ~isfolder(folder)
    error('run_simulation: folder %s not found.', folder);
end

old = cd(folder);
fprintf('=== case %d: %s (%s)\n', n, script, folder);
log_file = '';
if string(method) == "compare_n500_vs_n5000"   % keep the result (console is cleared later)
    log_file = fullfile(folder, 'comparison_result.txt');
    if isfile(log_file), delete(log_file); end
    diary(log_file);
end
try
    run_isolated(script);
catch err
    if ~isempty(log_file), diary off; end
    cd(old);
    rethrow(err);
end
if ~isempty(log_file), diary off; end
cd(old);
end

% -------------------------------------------------------------------------
function run_isolated(script__)
    % Separate workspace for the simulation script (it starts with "clear all").
    eval(script__);
end
