% COMPUTE_SPIN_CORRELATION  Spin-momentum correlation of Ehrenfest or koopmon results.
%
%   For k = 1,2,3 computes the covariance U_k = <p sigma_k> - <p><sigma_k>
%   and the correlation coefficient C_k = U_k / (std(p) * std(sigma_k)),
%   with averages over all particles.
%
%   Usage:   run from caseN/postprocessing/<topic>; select the Ehrenfest or
%            koopmon result in `string_name` below. The result is read
%            directly from the simulation folder of the case (not modified).
%            The method (K = koopmons, E = Ehrenfest) follows from the file.
%
%   Output (folder output/):  <method>_sc1.dat, _sc2.dat, _sc3.dat (time, C_k).

% Clear the workspace, but keep run_case_file / run_case_method if this
% script is called from run_case.m
clearvars -global
clearvars -except run_case_file run_case_method
close all;

%% Load simulation result (choose one)
% string_name = "../../ehrenfest/n500/ehrenfest_simulation.mat";
string_name = "../../koopmons/n500_a05/koopmons_simulation.mat";
if exist('run_case_file','var')
    string_name = run_case_file;    % set by run_case.m
end
load(string_name)
if ~isfolder('output')
    mkdir('output');     % results are written to output/
end

%% Method
% method follows from the loaded file: K = koopmons, E = Ehrenfest
if contains(string_name, "koopmons")
    plot_title = "K";
else
    plot_title = "E";
end
if plot_title == "K"
    prefix = 'koopmons';     % prefix of the output files
else
    prefix = 'ehrenfest';
end

t_steps = length(q_val(1,:));   % number of stored time steps

% Pauli matrices
sigma0 = [1,0;0,1];
sigma1 = [0,1;1,0];
sigma2 = [0,-1i;1i,0];
sigma3 = [1,0;0,-1];

%% Expectation values, covariances and correlation coefficients per time step
Ek = zeros(3,t_steps);
Ekk = zeros(3,t_steps);
Ep = zeros(1,t_steps);
Ekp = zeros(3,t_steps);
Epp = zeros(1,t_steps);
Uk = zeros(3,t_steps);
Ck = zeros(3,t_steps);
sig = {sigma1, sigma2, sigma3};
for k = 1:t_steps
    % E_k = <sigma_k>, E_kk = <sigma_k^2>, E_kp = <p sigma_k>
    for c = 1:3
        sum_k = 0;
        sum_kk = 0;
        sum_kp = 0;
        for n = 1:n_b
            sum_k = sum_k+trace(rho_val(:,:,n,k)*sig{c});
            sum_kk = sum_kk+trace(rho_val(:,:,n,k)*sig{c}*sig{c});
            sum_kp = sum_kp+p_val(n,k)*trace(rho_val(:,:,n,k)*sig{c});
        end
        Ek(c,k) = sum_k/n_b;
        Ekk(c,k) = sum_kk/n_b;
        Ekp(c,k) = sum_kp/n_b;
    end

    % E_p = <p>, E_pp = <p^2>
    sum_p = 0;
    sum_pp = 0;
    for n = 1:n_b
        sum_p = sum_p+p_val(n,k)*trace(rho_val(:,:,n,k)*sigma0);
        sum_pp = sum_pp+p_val(n,k)^2*trace(rho_val(:,:,n,k)*sigma0);
    end
    Ep(1,k) = sum_p/n_b;
    Epp(1,k) = sum_pp/n_b;

    % covariance U_k = <p sigma_k> - <p><sigma_k> and correlation coefficient C_k
    sigma_p = sqrt(Epp(1,k)-Ep(1,k)^2);
    for c = 1:3
        Uk(c,k) = Ekp(c,k)-Ek(c,k)*Ep(1,k);
        sigma_c = sqrt(Ekk(c,k)-Ek(c,k)^2);
        Ck(c,k) = corr_coeff(Uk(c,k),sigma_p,sigma_c);
    end
end

time = dt*(0:t_steps-1);    % first stored state is at t = 0

%% Plot covariances and correlation coefficients
figure
for c = 1:3
    subplot(3,2,2*c-1:2*c);
    plot(time,Uk(c,:))
    ylim([-1,1])
    grid on
    title(sprintf('sigma%d', c))
end

figure
for c = 1:3
    subplot(3,2,2*c-1:2*c);
    plot(time,Ck(c,:))
    ylim([-1,1])
    grid on
    title(sprintf('Csigma%d', c))
end

%% Export correlation coefficients (time, C_k) to output/<prefix>_sc1.dat ... _sc3.dat
for c = 1:3
    out = [time' real(Ck(c,:))'];
    writematrix(out, fullfile("output", sprintf('%s_sc%d.dat', prefix, c)), 'Delimiter',' ', 'FileType','text');
end

% -------------------------------------------------------------------------
function C = corr_coeff(U,sigma_a,sigma_b)
    % Correlation coefficient U/(sigma_a*sigma_b). It is undefined (NaN) if
    % one of the standard deviations vanishes, e.g. at t = 0 when the spin
    % points exactly along z.
    if abs(sigma_a*sigma_b) < 1e-12
        C = NaN;
    else
        C = U/(sigma_a*sigma_b);
    end
end
