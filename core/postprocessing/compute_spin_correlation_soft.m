% COMPUTE_SPIN_CORRELATION_SOFT  Spin-momentum correlation of the exact (SOFT) solution.
%
%   For k = 1,2,3 computes the covariance U_k = <p sigma_k> - <p><sigma_k>
%   and the correlation coefficient C_k = U_k / (std(p) * std(sigma_k)),
%   with quantum expectation values of the two-component wave function.
%   The momentum operator is applied in Fourier space.
%
%   Usage:   run from caseN/postprocessing/<topic>. The result is read
%            directly from caseN/soft/soft_simulation.mat (not modified).
%
%   Output (folder output/):  soft_sc1.dat, soft_sc2.dat, soft_sc3.dat (time, C_k).

% Clear the workspace, but keep run_case_file / run_case_method if this
% script is called from run_case.m
clearvars -global
clearvars -except run_case_file run_case_method
close all;

%% Load simulation result
string_name = "../../soft/soft_simulation.mat";
load(string_name)
if ~isfolder('output')
    mkdir('output');     % results are written to output/
end
prefix = 'soft';    % prefix of the output files

% Pauli matrices
sigma0 = [1,0;0,1];
sigma1 = [0,1;1,0];
sigma2 = [0,-1i;1i,0];
sigma3 = [1,0;0,-1];

%% Expectation values, covariances and correlation coefficients per time step
Ek = zeros(3,n_snap);
Ekk = zeros(3,n_snap);
Ep = zeros(1,n_snap);
Ekp = zeros(3,n_snap);
Epp = zeros(1,n_snap);
Uk = zeros(3,n_snap);
Ck = zeros(3,n_snap);
sig = {sigma1, sigma2, sigma3};
for k = 1:n_snap
    % wave function v and its derivative dv = -i d/dx v (via FFT), i.e. p applied to v
    v1 = psi1_val(:,k).';
    v2 = psi2_val(:,k).';
    v = [v1;v2];

    kvec = dp*[0:np/2-1,-np/2:-1];
    ftpsi1 = fft(v1)*dx/sqrt(2*pi);
    ftpsi2 = fft(v2)*dx/sqrt(2*pi);
    ftpsi1 = kvec.*ftpsi1;
    ftpsi2 = kvec.*ftpsi2;
    dv1 = ifft(ftpsi1)*dp/sqrt(2*pi)*np;
    dv2 = ifft(ftpsi2)*dp/sqrt(2*pi)*np;
    dv = [dv1; dv2];

    x_d = x(1:end);
    c1 = psi1_val(1:end,k)';    % complex conjugate of the components (rows)
    c2 = psi2_val(1:end,k)';

    % E_k = <sigma_k>, E_kk = <sigma_k^2>, E_kp = <p sigma_k>
    for c = 1:3
        Ek(c,k) = real(trapz(x_d,braket(c1,c2,sig{c}*v)));
        Ekk(c,k) = real(trapz(x_d,braket(c1,c2,sig{c}*sig{c}*v)));
        Ekp(c,k) = real(trapz(x_d,braket(c1,c2,sig{c}*dv)));
    end

    % E_p = <p>, E_pp = <p^2>
    Ep(1,k) = real(trapz(x_d,braket(c1,c2,sigma0*dv)));
    fEpp_pos = conj(dv1).*dv1+conj(dv2).*dv2;
    Epp(1,k) = real(trapz(x_d,fEpp_pos));

    % covariance U_k = <p sigma_k> - <p><sigma_k> and correlation coefficient C_k
    sigma_p = sqrt(Epp(1,k)-Ep(1,k)^2);
    for c = 1:3
        Uk(c,k) = Ekp(c,k)-Ek(c,k)*Ep(1,k);
        sigma_c = sqrt(Ekk(c,k)-Ek(c,k)^2);
        Ck(c,k) = corr_coeff(Uk(c,k),sigma_p,sigma_c);
    end
end

time = tau*(0:n_snap-1);    % snapshot l is at time (l-1)*tau

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
function f = braket(c1,c2,w)
    % Integrand psi^dagger * w for the conjugated components c1, c2 and the
    % two-component vector w (rows).
    f = c1.*w(1,:)+c2.*w(2,:);
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
