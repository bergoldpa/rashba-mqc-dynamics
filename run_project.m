% RUN_PROJECT  Run simulations and postprocessing from one place.
%
%   Comment in (or out) the lines you want to run, choose the test cases in
%   `cases`, and press Run. Each simulation runs in its own folder and the
%   results are written there; the postprocessing reads them from there.

setup                               % add core/ to the MATLAB path

cases = [1 2 3 4 5];        % test cases, e.g. [1 2 3 4 5]

for n = cases
    %% Simulations
    run_simulation(n, "soft");
    run_simulation(n, "ehrenfest");
    run_simulation(n, "koopmons");
    if n == 3                       % case 3 only: N = 5000 run and N = 500 vs 5000 comparison
        run_simulation(n, "ehrenfest_n5000");
        run_simulation(n, "compare_n500_vs_n5000");
    end

    %% Postprocessing (all parts, or select single parts)
    run_case(n);
    % run_case(n, "phase_space");
    % run_case(n, "bloch_purity");
    % run_case(n, "spin_correlation");
end

%% Other tools
run(fullfile(fileparts(mfilename('fullpath')), 'potential_energy_surfaces', 'momentum_surfaces.m'))
%                               % energy surfaces of all cases -> potential_energy_surfaces/figures_momentum
