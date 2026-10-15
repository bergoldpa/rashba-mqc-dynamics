% MAX_ENERGY_DRIFT  Maximum relative energy drift of a finished simulation.
%
%       drift = max_{0<=t<=T} |H(t)-H(0)| / |H(0)|
%
%   Reads the result file(s) of ehrenfest_simulation.m or
%   koopmons_simulation.m.
%
%   Usage:   run from the folder containing the result file, or list
%            several result files (with paths) in `files` below.

clear; clc

% Result file(s) to analyse - choose one line (add more files to the list if needed):
% files = {'koopmons_simulation'};
files = {'ehrenfest_simulation'};

drift_all = zeros(1,numel(files));

for i = 1:numel(files)
    S = load(files{i}, 'E_tot_val', 'l', 'L', 'dt');

    % The simulations save before "l = l-1", so the number
    % of stored time steps (t = 0, dt, ..., (n-1)*dt) is l-1
    n = S.l - 1;
    H = S.E_tot_val(1:n);
    t = (0:n-1)*S.dt;

    H0 = H(1);
    [drift, idx] = max(abs(H - H0) / abs(H0));
    drift_all(i) = drift;

    fprintf('%s\n', files{i});
    fprintf('  H(0)             = %.6e\n', H0);
    fprintf('  T                = %.6e  (%d of %d steps)\n', t(end), n, S.L);
    fprintf('  max|H-H0|/|H0|   = %.4e  (at t = %.4e)\n\n', drift, t(idx));
    if n < S.L
        warning('%s: run stopped early (energy tolerance exceeded).', files{i});
    end
end

fprintf('Overall max|H-H0|/|H0| = %.4e\n', max(drift_all));
