function run_case(n, parts)
% RUN_CASE  Run the complete postprocessing of test case n.
%
%   run_case(n) uses the existing simulation results of caseN
%   (soft/soft_simulation.mat, ehrenfest/n500/ehrenfest_simulation.mat,
%   koopmons/n500_a05/koopmons_simulation.mat) and
%     1. creates the phase-space snapshots (SOFT, Ehrenfest, koopmons),
%     2. computes Bloch vector and purity,
%     3. computes the spin-momentum correlations,
%     4. converts the time axes to ps or ns (P.time_scale), and
%     5. copies all files needed by the LaTeX figures into caseN/figures/*/input/.
%
%   run_case(n, parts) runs only the selected parts, e.g.
%   run_case(1, "phase_space") or run_case(1, ["bloch_purity","spin_correlation"]).
%   Possible parts: "phase_space", "bloch_purity", "spin_correlation".
%
%   The simulations themselves are not started (see run_simulation.m).

if nargin < 2
    parts = ["phase_space", "bloch_purity", "spin_correlation"];
end
parts = string(parts);

root = fileparts(fileparts(mfilename('fullpath')));   % CODE folder
C = fullfile(root, sprintf('case%d', n));
if ~isfolder(C)
    error('run_case: folder %s not found.', C);
end
post = fullfile(C, 'postprocessing');
figs = fullfile(C, 'figures');

E = "../../ehrenfest/n500/ehrenfest_simulation.mat";
K = "../../koopmons/n500_a05/koopmons_simulation.mat";

%% 1. Phase space (SOFT first: it creates the colour scale WIGNER_VAL_REF.mat)
if any(parts == "phase_space")
    d = fullfile(post, 'phase_space');
    run_in(d, 'plot_phase_space_soft', "", "");
    run_in(d, 'plot_phase_space_ehrenfest_koopmons', E, "");
    run_in(d, 'plot_phase_space_ehrenfest_koopmons', K, "");
    copy_files(d, 'ps*_*.png', fullfile(figs, 'phase_space', 'input'));
end

%% 2. Bloch vector and purity
if any(parts == "bloch_purity")
    d = fullfile(post, 'bloch_purity');
    run_in(d, 'plot_bloch_purity_soft', "", "");
    run_in(d, 'plot_bloch_purity', E, "");
    run_in(d, 'plot_bloch_purity', K, "");
    for m = ["S", "E", "K"]
        run_in(fullfile(d, 'output'), 'convert_bloch_purity_time', "", m);
    end
    copy_files(fullfile(d, 'output'), '*_N.dat', fullfile(figs, 'bloch_purity', 'input'));
end

%% 3. Spin-momentum correlation
if any(parts == "spin_correlation")
    d = fullfile(post, 'spin_correlation');
    run_in(d, 'compute_spin_correlation_soft', "", "");
    run_in(d, 'compute_spin_correlation', E, "");
    run_in(d, 'compute_spin_correlation', K, "");
    for m = ["S", "E", "K"]
        run_in(fullfile(d, 'output'), 'convert_spin_correlation_time', "", m);
    end
    copy_files(fullfile(d, 'output'), '*_N.dat', fullfile(figs, 'spin_correlation', 'input'));
end

fprintf('run_case(%d) finished (%s). Compile the .tex files in %s.\n', n, strjoin(parts, ', '), figs);
end

% -------------------------------------------------------------------------
function run_in(folder, script, file, method)
    % Run a postprocessing script in the given folder, passing the result
    % file and/or the method to it. The working folder is restored afterwards.
    if ~isfolder(folder)
        mkdir(folder);
    end
    old = cd(folder);
    fprintf('--- %s  (%s)\n', script, folder);
    try
        run_isolated(script, file, method);
    catch err
        cd(old);
        rethrow(err);
    end
    cd(old);
end

% -------------------------------------------------------------------------
function run_isolated(script__, run_case_file, run_case_method)
    % Separate workspace for the script; empty arguments are removed so the
    % script uses its own default.
    if run_case_file == ""
        clear run_case_file
    end
    if run_case_method == ""
        clear run_case_method
    end
    eval(script__);
end

% -------------------------------------------------------------------------
function copy_files(src, pattern, dst)
    % Copy all files matching pattern from src to dst.
    if ~isfolder(dst)
        mkdir(dst);
    end
    f = dir(fullfile(src, pattern));
    for i = 1:numel(f)
        copyfile(fullfile(src, f(i).name), dst);
    end
end
