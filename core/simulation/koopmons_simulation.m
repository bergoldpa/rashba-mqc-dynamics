% KOOPMONS_SIMULATION  Mixed quantum-classical dynamics with koopmons.
%
%   Propagates n_b koopmons (q,p), each carrying a 2x2 electronic density
%   matrix rho. In contrast to Ehrenfest dynamics, the equations of motion
%   contain quantum-classical correlation terms. These are computed with
%   Gaussian kernels of width sigma_kernel_x / sigma_kernel_mu, integrated
%   with the trapezoidal rule on a uniform (x,mu) grid. Time integration
%   uses a fourth-order Runge-Kutta scheme. The run stops at t_max or when
%   the relative error of the total energy exceeds tol_energy.
%
%   Usage:   cd into a run folder (e.g. case1/koopmons/n500_a05) and run
%            koopmons_simulation. Parameters: koopmons_input.m and the
%            case_parameters.m of the case.
%
%   Output (written to the current folder):
%     koopmons_simulation.mat    complete workspace (trajectories, energies, ...)
%     time_seconds.dat           wall-clock run time
%     plot_energy.fig            energy components over time
%     plot_energy_error.fig      relative energy error over time

clear all; close all; clc

%% Parameters
koopmons_input

tic

%% Initialisation
L = ceil(t_max/dt)+1; % total number of time steps
q_val = zeros(n_b,L); % classical distribution (position space)
p_val = zeros(n_b,L); % classical distribution (momentum space)
rho_val = zeros(2,2,n_b,L); % quantum density matrix
E_tot_val = zeros(1,L); % total energy
E_kin_val = zeros(1,L); % kinetic energy
E_mfd_val = zeros(1,L); % mean field energy
E_dif_val = zeros(1,L); % correlation energy

[w,q,p,rho] = fun_initialize();
q = rescale1*q;
p = rescale2*p;
pos_var = var(q);
mom_var = var(p);
[E_tot, E_kin, E_mfd, E_dif] = fun_energy(w,q/rescale1,p/rescale2,rho);

%% Time propagation
n_steps = L-1;                  % number of time steps
fprintf('  0 %%   relative energy error %.2e\n', 0);
l = 1;
E_tot0 = E_tot;

while l < L+1 && abs((E_tot-E_tot0)/E_tot0) < tol_energy
    % store current state
    q_val(:,l) = q;
    p_val(:,l) = p;
    rho_val(:,:,:,l) = rho;
    E_tot_val(l) = E_tot;
    E_kin_val(l) = E_kin;
    E_mfd_val(l) = E_mfd;
    E_dif_val(l) = E_dif;

    % one RK4 step (computed in unscaled variables)
    [q,p,rho] = fun_propagate(w,q/rescale1,p/rescale2,rho);
    q = rescale1*q;
    p = rescale2*p;

    l = l+1;
    [E_tot, E_kin, E_mfd, E_dif] = fun_energy(w,q/rescale1,p/rescale2,rho);

    % progress (every 10 %) and relative energy error |E(t)-E(0)|/|E(0)|
    if mod(l-1, round(n_steps/10)) == 0
        fprintf('%3.0f %%   relative energy error %.2e\n', 100*(l-1)/n_steps, abs((E_tot-E_tot0)/E_tot0));
    end
end
toc
time_seconds = toc;

% keep only the computed time steps (the run stops early if the relative
% energy error exceeds tol_energy)
n_done = l-1;
if n_done < L
    warning('Run stopped early after %d of %d time steps (energy tolerance exceeded).', n_done, L);
end
q_val = q_val(:,1:n_done);
p_val = p_val(:,1:n_done);
rho_val = rho_val(:,:,:,1:n_done);
E_tot_val = E_tot_val(1:n_done);
E_kin_val = E_kin_val(1:n_done);
E_mfd_val = E_mfd_val(1:n_done);
E_dif_val = E_dif_val(1:n_done);

%% Save results
writematrix(time_seconds, 'time_seconds.dat', 'FileType','text');
save(mfilename,'-v7.3')

%% Plot energy components
l = l-1;    % number of stored time steps
figure
plot((0:l-1)*dt,E_tot_val(1:l),'-k');
hold on
plot((0:l-1)*dt,E_kin_val(1:l),'-b');
plot((0:l-1)*dt,E_mfd_val(1:l),'-r');
plot((0:l-1)*dt,E_dif_val(1:l),'-g');
xlabel('time')
legend('tot','kin','mfd','cor','Location','northwest')
str = strcat({'Energy, '}, num2str(n_b),...
    {' Koompons, varp = '}, num2str(round(pos_var,2)),...
    {', varm = '}, num2str(round(mom_var,2)),...
    {', alpha = '}, num2str(alpha));
sgtitle(str)
savefig('plot_energy.fig')
close

%% Plot relative energy error |E(t)-E(0)|/|E(0)|
rel_err = abs((E_tot_val(1:l)-E_tot_val(1))/E_tot_val(1));
if P.time_scale == "NANO"
    t_plot = (0:l-1)*dt*2.4188843265857e-17*1e9;
    t_label = 'time (ns)';
else
    t_plot = (0:l-1)*dt*2.4188843265857e-17*1e12;
    t_label = 'time (ps)';
end
figure
semilogy(t_plot(2:end),rel_err(2:end),'-k')    % t = 0 omitted (error is zero there)
grid on
xlabel(t_label)
ylabel('|E(t)-E(0)| / |E(0)|')
title(['Relative energy error, koopmons, N = ' num2str(n_b)])
savefig('plot_energy_error.fig')
close

%% ========================================================================
%  Local functions
%
%  Kernel naming: fun_K1a2b(qk,x,pk,mu) = K_a(x-qk) * K_b(mu-pk), where
%  K_0 is a normalised Gaussian and K_1, K_2 are its first and second
%  derivatives (a: derivative order in x, b: derivative order in mu).
%  fun_J1a2bij is the kernel-weighted ensemble average of rho_ij.
%  ========================================================================

% -------------------------------------------------------------------------
function kernel = fun_K0(sigma_kernel,y)
    % Normalised one-dimensional Gaussian kernel with standard deviation sigma_kernel.
    kernel = exp(-.5*(y/sigma_kernel).^2)/sqrt(2*pi)/sigma_kernel;
end

% -------------------------------------------------------------------------
function kernel1 = fun_K1(sigma_kernel,y)
    % First derivative of the Gaussian kernel fun_K0.
    kernel1 = -y/sigma_kernel^2.*fun_K0(sigma_kernel,y);
end

% -------------------------------------------------------------------------
function kernel2 = fun_K2(sigma_kernel,y)
    % Second derivative of the Gaussian kernel fun_K0.
    kernel2 = -1/sigma_kernel^2*(fun_K0(sigma_kernel,y)...
        +y.*fun_K1(sigma_kernel,y));
end

% -------------------------------------------------------------------------
function K1020 = fun_K1020(qk,x,pk,mu)
    % K_0(x-qk) * K_0(mu-pk)
    global sigma_kernel_x sigma_kernel_mu
    K1020 = fun_K0(sigma_kernel_x,x-qk).*fun_K0(sigma_kernel_mu,mu-pk);
end

% -------------------------------------------------------------------------
function K1021 = fun_K1021(qk,x,pk,mu)
    % K_0(x-qk) * K_1(mu-pk)
    global sigma_kernel_x sigma_kernel_mu
    K1021 = fun_K0(sigma_kernel_x,x-qk).*fun_K1(sigma_kernel_mu,mu-pk);
end

% -------------------------------------------------------------------------
function K1022 = fun_K1022(qk,x,pk,mu)
    % K_0(x-qk) * K_2(mu-pk)
    global sigma_kernel_x sigma_kernel_mu
    K1022 = fun_K0(sigma_kernel_x,x-qk).*fun_K2(sigma_kernel_mu,mu-pk);
end

% -------------------------------------------------------------------------
function K1120 = fun_K1120(qk,x,pk,mu)
    % K_1(x-qk) * K_0(mu-pk)
    global sigma_kernel_x sigma_kernel_mu
    K1120 = fun_K1(sigma_kernel_x,x-qk).*fun_K0(sigma_kernel_mu,mu-pk);
end

% -------------------------------------------------------------------------
function K1121 = fun_K1121(qk,x,pk,mu)
    % K_1(x-qk) * K_1(mu-pk)
    global sigma_kernel_x sigma_kernel_mu
    K1121 = fun_K1(sigma_kernel_x,x-qk).*fun_K1(sigma_kernel_mu,mu-pk);
end

% -------------------------------------------------------------------------
function K1220 = fun_K1220(qk,x,pk,mu)
    % K_2(x-qk) * K_0(mu-pk)
    global sigma_kernel_x sigma_kernel_mu
    K1220 = fun_K2(sigma_kernel_x,x-qk).*fun_K0(sigma_kernel_mu,mu-pk);
end

% -------------------------------------------------------------------------
function J1020ij = fun_J1020ij(i,j,q,p,w,rho,x,mu)
    % Kernel-weighted average of rho_ij on the grid, numerator kernel K1020.
    % (The small constant in the denominator avoids division by zero.)
    global n_b
    sum_n_J = 0;
    sum_d_J = 1e-12;
    for k = 1:n_b
        sum_n_J = sum_n_J+w(k)*rho(i,j,k).*fun_K1020(q(k),x,p(k),mu);
        sum_d_J = sum_d_J+w(k).*fun_K1020(q(k),x,p(k),mu);
    end
    J1020ij = sum_n_J./sum_d_J;
end

% -------------------------------------------------------------------------
function J1021ij = fun_J1021ij(i,j,q,p,w,rho,x,mu)
    % Kernel-weighted average of rho_ij on the grid, numerator kernel K1021.
    global n_b
    sum_n_J = 0;
    sum_d_J = 1e-12;
    for k = 1:n_b
        sum_n_J = sum_n_J+w(k)*rho(i,j,k).*fun_K1021(q(k),x,p(k),mu);
        sum_d_J = sum_d_J+w(k).*fun_K1020(q(k),x,p(k),mu);
    end
    J1021ij = sum_n_J./sum_d_J;
end

% -------------------------------------------------------------------------
function J1120ij = fun_J1120ij(i,j,q,p,w,rho,x,mu)
    % Kernel-weighted average of rho_ij on the grid, numerator kernel K1120.
    global n_b
    sum_n_J = 0;
    sum_d_J = 1e-12;
    for k = 1:n_b
        sum_n_J = sum_n_J+w(k)*rho(i,j,k).*fun_K1120(q(k),x,p(k),mu);
        sum_d_J = sum_d_J+w(k).*fun_K1020(q(k),x,p(k),mu);
    end
    J1120ij = sum_n_J./sum_d_J;
end

% -------------------------------------------------------------------------
function grid = fun_grid_x(q)
    % Uniform grid in x covering all particles plus nx kernel widths.
    global dx nx sigma_kernel_x
    q_min = min(q)-nx*sigma_kernel_x;
    q_max = max(q)+nx*sigma_kernel_x;
    grid = q_min:dx:q_max;
end

% -------------------------------------------------------------------------
function grid = fun_grid_mu(p)
    % Uniform grid in mu covering all particles plus nmu kernel widths.
    global dmu nmu sigma_kernel_mu
    p_min = min(p)-nmu*sigma_kernel_mu;
    p_max = max(p)+nmu*sigma_kernel_mu;
    grid = p_min:dmu:p_max;
end

% -------------------------------------------------------------------------
function [w,q,p,rho] = fun_initialize()
    % Sample the initial ensemble (quasi-random Sobol sampling): equal
    % weights w, Gaussian positions q and momenta p, identical density matrices rho.
    global n_b q_init p_init rho_init sigma_init_q sigma_init_p
    % quasi-random Sobol points, mapped to normal distributions
    sob_z = sobolset(2,'Skip',10,'Leap',10);
    u_z = net(sob_z,n_b);
    norm_q = norminv(u_z(:,1));     % inverse CDF transform
    norm_p = norminv(u_z(:,2));     % inverse CDF transform
    w = zeros(1,n_b);
    q = zeros(1,n_b);
    p = zeros(1,n_b);
    rho = zeros(2,2,n_b);
    for k = 1:n_b
        w(k) = 1/n_b;
        q(k) = q_init+sigma_init_q*norm_q(k);
        p(k) = p_init+sigma_init_p*norm_p(k);
        rho(:,:,k) = rho_init;
    end
end

% -------------------------------------------------------------------------
function [E_tot, E_kin, E_mfd, E_dif] = fun_energy(w,q,p,rho)
    % Ensemble-averaged total, kinetic, mean-field and correlation energy.
    global n_b M hbar lambda_x lambda_mu

    % integration grid in (x,mu)
    grid_x = fun_grid_x(q);
    grid_mu = fun_grid_mu(p);
    [grid_X, grid_Mu] = meshgrid(grid_x,grid_mu);

    % kernel-averaged density matrices on the grid (entry 21 by hermiticity)
    J1021_11_grid = fun_J1021ij(1,1,q,p,w,rho,grid_X,grid_Mu);
    J1021_12_grid = fun_J1021ij(1,2,q,p,w,rho,grid_X,grid_Mu);
    J1021_22_grid = fun_J1021ij(2,2,q,p,w,rho,grid_X,grid_Mu);
    J1021_21_grid = conj(J1021_12_grid);

    J1120_11_grid = fun_J1120ij(1,1,q,p,w,rho,grid_X,grid_Mu);
    J1120_12_grid = fun_J1120ij(1,2,q,p,w,rho,grid_X,grid_Mu);
    J1120_22_grid = fun_J1120ij(2,2,q,p,w,rho,grid_X,grid_Mu);
    J1120_21_grid = conj(J1120_12_grid);

    % derivatives of the electronic Hamiltonian on the grid
    [DHe11_x_grid, DHe11_mu_grid] = fun_DHe11_full(grid_X,grid_Mu);
    [DHe12_x_grid, DHe12_mu_grid] = fun_DHe12_full(grid_X,grid_Mu);
    [DHe22_x_grid, DHe22_mu_grid] = fun_DHe22_full(grid_X,grid_Mu);
    DHe21_x_grid = conj(DHe12_x_grid);
    DHe21_mu_grid = conj(DHe12_mu_grid);

    % commutators of the Hamiltonian derivatives with the averaged density matrices
    [com1021_11_x_grid,com1021_12_x_grid,com1021_21_x_grid,com1021_22_x_grid] =...
        fun_commutator(DHe11_x_grid,DHe12_x_grid,DHe21_x_grid,DHe22_x_grid,...
        J1021_11_grid,J1021_12_grid,J1021_21_grid,J1021_22_grid);

    [com1120_11_mu_grid,com1120_12_mu_grid,com1120_21_mu_grid,com1120_22_mu_grid] =...
        fun_commutator(DHe11_mu_grid,DHe12_mu_grid,DHe21_mu_grid,DHe22_mu_grid,...
        J1120_11_grid,J1120_12_grid,J1120_21_grid,J1120_22_grid);

    Ek_tot = zeros(1,n_b);
    Ek_kin = zeros(1,n_b);
    Ek_mfd = zeros(1,n_b);
    Ek_cor = zeros(1,n_b);
    for k = 1:n_b
        % kernel centred at particle k, truncated by the window
        K1020k_grid = fun_K1020(q(k),grid_X,p(k),grid_Mu);
        windowk_grid = fun_window(q(k),grid_X,lambda_x,p(k),grid_Mu,lambda_mu);
        K1020k_grid = K1020k_grid.*windowk_grid;

        % correlation energy of particle k
        trace_x_grid = 1i*imag(fun_trace(rho(1,1,k),rho(1,2,k),rho(2,1,k),rho(2,2,k),...
            com1021_11_x_grid,com1021_12_x_grid,com1021_21_x_grid,com1021_22_x_grid));
        trace_mu_grid = 1i*imag(fun_trace(rho(1,1,k),rho(1,2,k),rho(2,1,k),rho(2,2,k),...
            com1120_11_mu_grid,com1120_12_mu_grid,com1120_21_mu_grid,com1120_22_mu_grid));

        int_x_grid = K1020k_grid.*trace_x_grid;
        int_mu_grid = K1020k_grid.*trace_mu_grid;

        int_xk = trapz(grid_mu,trapz(grid_x,int_x_grid,2));
        int_muk = trapz(grid_mu,trapz(grid_x,int_mu_grid,2));

        % kinetic, mean-field and correlation energy of particle k
        [He_fullk,~,~] = fun_He_full(q(k),p(k));
        Ek_kin(k) = p(k)^2/2/M;
        Ek_mfd(k) = real(trace(rho(:,:,k)'*He_fullk));
        Ek_cor(k) = 1i*hbar/2*int_xk-1i*hbar/2*int_muk;
        Ek_tot(k) = Ek_kin(k)+Ek_mfd(k)+Ek_cor(k);
    end
    E_tot = real(sum(w.*Ek_tot));
    E_kin = real(sum(w.*Ek_kin));
    E_mfd = real(sum(w.*Ek_mfd));
    E_dif = real(sum(w.*Ek_cor));
end

% -------------------------------------------------------------------------
function [qdot,pdot,rhodot] = fun_dot(w,q,p,rho)
    % Right-hand side of the koopmon equations of motion.
    global n_b M hbar lambda_x lambda_mu
    qdot = zeros(1,n_b);
    pdot = zeros(1,n_b);
    rhodot = zeros(2,2,n_b);

    % integration grid in (x,mu)
    grid_x = fun_grid_x(q);
    grid_mu = fun_grid_mu(p);
    [grid_X, grid_Mu] = meshgrid(grid_x,grid_mu);

    % kernel-averaged density matrices on the grid (entry 21 by hermiticity)
    J1020_11_grid = fun_J1020ij(1,1,q,p,w,rho,grid_X,grid_Mu);
    J1020_12_grid = fun_J1020ij(1,2,q,p,w,rho,grid_X,grid_Mu);
    J1020_22_grid = fun_J1020ij(2,2,q,p,w,rho,grid_X,grid_Mu);
    J1020_21_grid = conj(J1020_12_grid);

    J1021_11_grid = fun_J1021ij(1,1,q,p,w,rho,grid_X,grid_Mu);
    J1021_12_grid = fun_J1021ij(1,2,q,p,w,rho,grid_X,grid_Mu);
    J1021_22_grid = fun_J1021ij(2,2,q,p,w,rho,grid_X,grid_Mu);
    J1021_21_grid = conj(J1021_12_grid);

    J1120_11_grid = fun_J1120ij(1,1,q,p,w,rho,grid_X,grid_Mu);
    J1120_12_grid = fun_J1120ij(1,2,q,p,w,rho,grid_X,grid_Mu);
    J1120_22_grid = fun_J1120ij(2,2,q,p,w,rho,grid_X,grid_Mu);
    J1120_21_grid = conj(J1120_12_grid);

    % derivatives of the electronic Hamiltonian on the grid
    [DHe11_x_grid, DHe11_mu_grid] = fun_DHe11_full(grid_X,grid_Mu);
    [DHe12_x_grid, DHe12_mu_grid] = fun_DHe12_full(grid_X,grid_Mu);
    [DHe22_x_grid, DHe22_mu_grid] = fun_DHe22_full(grid_X,grid_Mu);
    DHe21_x_grid = conj(DHe12_x_grid);
    DHe21_mu_grid = conj(DHe12_mu_grid);

    % commutators of the Hamiltonian derivatives with the averaged density matrices
    [com1020_11_x_grid,com1020_12_x_grid,com1020_21_x_grid,com1020_22_x_grid] =...
        fun_commutator(DHe11_x_grid,DHe12_x_grid,DHe21_x_grid,DHe22_x_grid,...
        J1020_11_grid,J1020_12_grid,J1020_21_grid,J1020_22_grid);

    [com1021_11_x_grid,com1021_12_x_grid,com1021_21_x_grid,com1021_22_x_grid] =...
        fun_commutator(DHe11_x_grid,DHe12_x_grid,DHe21_x_grid,DHe22_x_grid,...
        J1021_11_grid,J1021_12_grid,J1021_21_grid,J1021_22_grid);

    [com1120_11_mu_grid,com1120_12_mu_grid,com1120_21_mu_grid,com1120_22_mu_grid] =...
        fun_commutator(DHe11_mu_grid,DHe12_mu_grid,DHe21_mu_grid,DHe22_mu_grid,...
        J1120_11_grid,J1120_12_grid,J1120_21_grid,J1120_22_grid);

    [com1020_11_mu_grid,com1020_12_mu_grid,com1020_21_mu_grid,com1020_22_mu_grid] =...
        fun_commutator(DHe11_mu_grid,DHe12_mu_grid,DHe21_mu_grid,DHe22_mu_grid,...
        J1020_11_grid,J1020_12_grid,J1020_21_grid,J1020_22_grid);

    % not needed in the equations of motion (would pair two momentum derivatives):
    % [com1021_11_mu_grid,com1021_12_mu_grid,com1021_21_mu_grid,com1021_22_mu_grid] =...
    %     fun_commutator(DHe11_mu_grid,DHe12_mu_grid,DHe21_mu_grid,DHe22_mu_grid,...
    %     J1021_11_grid,J1021_12_grid,J1021_21_grid,J1021_22_grid);

    % terms quadratic in the averaged density matrices
    [comJJ_11_x_grid,comJJ_12_x_grid,comJJ_21_x_grid,comJJ_22_x_grid] =...
        fun_commutator(J1021_11_grid,J1021_12_grid,J1021_21_grid,J1021_22_grid,...
        J1020_11_grid,J1020_12_grid,J1020_21_grid,J1020_22_grid);

    traceJJ_x_grid = 1i*imag(fun_trace(DHe11_x_grid,DHe12_x_grid,DHe21_x_grid,DHe22_x_grid,...
        comJJ_11_x_grid,comJJ_12_x_grid,comJJ_21_x_grid,comJJ_22_x_grid));

    [comJJ_11_mu_grid,comJJ_12_mu_grid,comJJ_21_mu_grid,comJJ_22_mu_grid] =...
        fun_commutator(J1120_11_grid,J1120_12_grid,J1120_21_grid,J1120_22_grid,...
        J1020_11_grid,J1020_12_grid,J1020_21_grid,J1020_22_grid);

    traceJJ_mu_grid = 1i*imag(fun_trace(DHe11_mu_grid,DHe12_mu_grid,DHe21_mu_grid,DHe22_mu_grid,...
        comJJ_11_mu_grid,comJJ_12_mu_grid,comJJ_21_mu_grid,comJJ_22_mu_grid));

    for k = 1:n_b
        % kernels centred at particle k, truncated by the window
        K1020k_grid = fun_K1020(q(k),grid_X,p(k),grid_Mu);
        K1021k_grid = fun_K1021(q(k),grid_X,p(k),grid_Mu);
        K1022k_grid = fun_K1022(q(k),grid_X,p(k),grid_Mu);
        K1120k_grid = fun_K1120(q(k),grid_X,p(k),grid_Mu);
        K1121k_grid = fun_K1121(q(k),grid_X,p(k),grid_Mu);
        K1220k_grid = fun_K1220(q(k),grid_X,p(k),grid_Mu);

        windowk_grid = fun_window(q(k),grid_X,lambda_x,p(k),grid_Mu,lambda_mu);
        K1020k_grid = K1020k_grid.*windowk_grid;
        K1021k_grid = K1021k_grid.*windowk_grid;
        K1022k_grid = K1022k_grid.*windowk_grid;
        K1120k_grid = K1120k_grid.*windowk_grid;
        K1121k_grid = K1121k_grid.*windowk_grid;
        K1220k_grid = K1220k_grid.*windowk_grid;

        % traces with the density matrix of particle k
        trace1020_x_grid = 1i*imag(fun_trace(rho(1,1,k),rho(1,2,k),rho(2,1,k),rho(2,2,k),...
            com1020_11_x_grid,com1020_12_x_grid,com1020_21_x_grid,com1020_22_x_grid));
        trace1021_x_grid = 1i*imag(fun_trace(rho(1,1,k),rho(1,2,k),rho(2,1,k),rho(2,2,k),...
            com1021_11_x_grid,com1021_12_x_grid,com1021_21_x_grid,com1021_22_x_grid));

        trace1020_mu_grid = 1i*imag(fun_trace(rho(1,1,k),rho(1,2,k),rho(2,1,k),rho(2,2,k),...
            com1020_11_mu_grid,com1020_12_mu_grid,com1020_21_mu_grid,com1020_22_mu_grid));
        trace1120_mu_grid = 1i*imag(fun_trace(rho(1,1,k),rho(1,2,k),rho(2,1,k),rho(2,2,k),...
            com1120_11_mu_grid,com1120_12_mu_grid,com1120_21_mu_grid,com1120_22_mu_grid));

        % momentum equation: correlation integral
        int_grid_pk_x1 = K1120k_grid.*trace1021_x_grid;
        int_grid_pk_x2 = K1121k_grid.*trace1020_x_grid;
        int_grid_pk_x3 = K1120k_grid.*traceJJ_x_grid;
        int_grid_pk_mu1 = K1120k_grid.*trace1120_mu_grid;
        int_grid_pk_mu2 = K1220k_grid.*trace1020_mu_grid;
        int_grid_pk_mu3 = K1120k_grid.*traceJJ_mu_grid;

        int_grid_pk_x = int_grid_pk_x1-int_grid_pk_x2-int_grid_pk_x3;
        int_grid_pk_mu = -int_grid_pk_mu1+int_grid_pk_mu2+int_grid_pk_mu3;

        int_pk_x = trapz(grid_mu,trapz(grid_x,int_grid_pk_x,2));
        int_pk_mu = trapz(grid_mu,trapz(grid_x,int_grid_pk_mu,2));
        int_pk = int_pk_x+int_pk_mu;

        % momentum equation
        [He_fullk, DHe_xk, DHe_muk] = fun_He_full(q(k),p(k));

        pdot(k) = -trace(rho(:,:,k)'*DHe_xk)+1i*hbar/2*int_pk;

        % position equation: correlation integral
        int_grid_qk_x1 = K1021k_grid.*trace1021_x_grid;
        int_grid_qk_x2 = K1022k_grid.*trace1020_x_grid;
        int_grid_qk_x3 = K1021k_grid.*traceJJ_x_grid;
        int_grid_qk_mu1 = K1021k_grid.*trace1120_mu_grid;
        int_grid_qk_mu2 = K1121k_grid.*trace1020_mu_grid;
        int_grid_qk_mu3 = K1021k_grid.*traceJJ_mu_grid;

        int_grid_qk_x = -int_grid_qk_x1+int_grid_qk_x2+int_grid_qk_x3;
        int_grid_qk_mu = int_grid_qk_mu1-int_grid_qk_mu2-int_grid_qk_mu3;

        int_qk_x = trapz(grid_mu,trapz(grid_x,int_grid_qk_x,2));
        int_qk_mu = trapz(grid_mu,trapz(grid_x,int_grid_qk_mu,2));
        int_qk = int_qk_x+int_qk_mu;

        % position equation
        qdot(k) = p(k)/M+trace(rho(:,:,k)'*DHe_muk)+1i*hbar/2*int_qk;

        % density matrix equation: correlation correction to the effective Hamiltonian
        int_grid_rhokC1_11 = -K1020k_grid.*com1120_11_mu_grid;
        int_grid_rhokC1_12 = -K1020k_grid.*com1120_12_mu_grid;
        int_grid_rhokC1_21 = -K1020k_grid.*com1120_21_mu_grid;
        int_grid_rhokC1_22 = -K1020k_grid.*com1120_22_mu_grid;

        int_grid_rhokC2_11 = K1020k_grid.*com1021_11_x_grid;
        int_grid_rhokC2_12 = K1020k_grid.*com1021_12_x_grid;
        int_grid_rhokC2_21 = K1020k_grid.*com1021_21_x_grid;
        int_grid_rhokC2_22 = K1020k_grid.*com1021_22_x_grid;

        int_grid_rhokC3_11 = K1120k_grid.*com1020_11_mu_grid;
        int_grid_rhokC3_12 = K1120k_grid.*com1020_12_mu_grid;
        int_grid_rhokC3_21 = K1120k_grid.*com1020_21_mu_grid;
        int_grid_rhokC3_22 = K1120k_grid.*com1020_22_mu_grid;

        int_grid_rhokC4_11 = -K1021k_grid.*com1020_11_x_grid;
        int_grid_rhokC4_12 = -K1021k_grid.*com1020_12_x_grid;
        int_grid_rhokC4_21 = -K1021k_grid.*com1020_21_x_grid;
        int_grid_rhokC4_22 = -K1021k_grid.*com1020_22_x_grid;

        int_grid_rhok11 = int_grid_rhokC2_11+int_grid_rhokC4_11...
            +int_grid_rhokC1_11+int_grid_rhokC3_11;
        int_grid_rhok12 = int_grid_rhokC2_12+int_grid_rhokC4_12...
            +int_grid_rhokC1_12+int_grid_rhokC3_12;
        int_grid_rhok21 = int_grid_rhokC2_21+int_grid_rhokC4_21...
            +int_grid_rhokC1_21+int_grid_rhokC3_21;
        int_grid_rhok22 = int_grid_rhokC2_22+int_grid_rhokC4_22...
            +int_grid_rhokC1_22+int_grid_rhokC3_22;

        int_rhok_11 = trapz(grid_mu,trapz(grid_x,int_grid_rhok11,2));
        int_rhok_12 = trapz(grid_mu,trapz(grid_x,int_grid_rhok12,2));
        int_rhok_21 = trapz(grid_mu,trapz(grid_x,int_grid_rhok21,2));
        int_rhok_22 = trapz(grid_mu,trapz(grid_x,int_grid_rhok22,2));

        % density matrix equation (von Neumann equation with effective Hamiltonian)
        com_rhok = He_fullk+1i*hbar/2*...
            [int_rhok_11, int_rhok_12; int_rhok_21, int_rhok_22];
        rhodot(:,:,k) = -1i/hbar*(com_rhok*rho(:,:,k)-rho(:,:,k)*com_rhok);
    end
    % remove round-off imaginary parts
    pdot = real(pdot);
    qdot = real(qdot);
end

% -------------------------------------------------------------------------
function [q,p,rho] = fun_propagate(w,q,p,rho)
    % One step of the classical fourth-order Runge-Kutta method.
    global dt
    q_old = q;
    p_old = p;
    rho_old = rho;

    % stage 1
    [qdot,pdot,rhodot] = fun_dot(w,q,p,rho);
    dq1 = dt*qdot;
    q = q_old+dq1/2;
    dp1 = dt*pdot;
    p = p_old+dp1/2;
    drho1 = dt*rhodot;
    rho = rho_old+drho1/2;

    % stage 2
    [qdot,pdot,rhodot] = fun_dot(w,q,p,rho);
    dq2 = dt*qdot;
    q = q_old+dq2/2;
    dp2 = dt*pdot;
    p = p_old+dp2/2;
    drho2 = dt*rhodot;
    rho = rho_old+drho2/2;

    % stage 3
    [qdot,pdot,rhodot] = fun_dot(w,q,p,rho);
    dq3 = dt*qdot;
    q = q_old+dq3;
    dp3 = dt*pdot;
    p = p_old+dp3;
    drho3 = dt*rhodot;
    rho = rho_old+drho3;

    % stage 4 and update
    [qdot,pdot,rhodot] = fun_dot(w,q,p,rho);
    dq4 = dt*qdot;
    q = q_old+dq1/6+dq2/3+dq3/3+dq4/6;
    dp4 = dt*pdot;
    p = p_old+dp1/6+dp2/3+dp3/3+dp4/6;
    drho4 = dt*rhodot;
    rho = rho_old+drho1/6+drho2/3+drho3/3+drho4/6;
end

% -------------------------------------------------------------------------
function trace = fun_trace(a11,a12,a21,a22,b11,b12,b21,b22)
    % tr(A^dagger * B) for 2x2 matrices given entrywise (works on whole grids).
    trace = conj(a11).*b11+conj(a21).*b21+conj(a12).*b12+conj(a22).*b22;
end

% -------------------------------------------------------------------------
function [com11,com12,com21,com22] = fun_commutator(a11,a12,a21,a22,b11,b12,b21,b22)
    % Commutator [A,B] = A*B - B*A for 2x2 matrices given entrywise.
    com11 = a12.*b21-b12.*a21;
    com12 = a11.*b12+a12.*b22-b11.*a12-b12.*a22;
    com21 = a21.*b11+a22.*b21-b21.*a11-b22.*a21;
    com22 = a21.*b12-b21.*a12;
end

% -------------------------------------------------------------------------
function window = fun_window(qk,x,lambda_x,pk,mu,lambda_mu)
    % Indicator of the box |x-qk| < lambda_x, |mu-pk| < lambda_mu (truncates the kernels).
    window = (abs(x-qk)<lambda_x).*(abs(mu-pk)<lambda_mu);
end

% -------------------------------------------------------------------------
function [He_full, DHe_x, DHe_mu] = fun_He_full(x,mu)
    % Electronic Hamiltonian He(x,mu) and its derivatives with respect to
    % position x and momentum mu (MomC model).
    global model
    if (model == "MomC")
        % momentum coupling: Rashba (a_R) and Dresselhaus (a_D) terms,
        % constant field B0 and optional harmonic trap
        global M
        [a_R,a_D,omega,B0,harmonic] = parameters_MomC();
        sigma0 = eye(2);
        sigma1 = [0,1;1,0];
        sigma2 = [0,-1i;1i,0];
        He_full = harmonic*.5*M*omega^2*x.^2*sigma0+a_R*mu*sigma2+a_D*mu*sigma1+B0*sigma1;
        DHe_x = harmonic*M*omega^2*x*sigma0;
        DHe_mu = a_R*sigma2+a_D*sigma1;
    else
        error('Only the model "MomC" is implemented (model = "%s").', model);
    end
end

% -------------------------------------------------------------------------
function [DHe11_x, DHe11_mu] = fun_DHe11_full(x,mu)
    % Derivatives of the entry He(1,1) with respect to x and mu (on a grid).
    global model
    if (model == "MomC")
        global M
        [~,~,omega,~,harmonic] = parameters_MomC();
        DHe11_x = harmonic*M*omega^2*x;
        DHe11_mu = 0*mu;
    else
        error('Only the model "MomC" is implemented (model = "%s").', model);
    end
end

% -------------------------------------------------------------------------
function [DHe12_x, DHe12_mu] = fun_DHe12_full(x,mu)
    % Derivatives of the entry He(1,2) with respect to x and mu (on a grid).
    global model
    if (model == "MomC")
        [a_R,a_D,~,~] = parameters_MomC();
        DHe12_x = 0*x;
        DHe12_mu = a_R*(-1i)+a_D;
    else
        error('Only the model "MomC" is implemented (model = "%s").', model);
    end
end

% -------------------------------------------------------------------------
function [DHe22_x, DHe22_mu] = fun_DHe22_full(x,mu)
    % Derivatives of the entry He(2,2) with respect to x and mu (on a grid).
    global model
    if (model == "MomC")
        global M
        [~,~,omega,~,harmonic] = parameters_MomC();
        DHe22_x = harmonic*M*omega^2*x;
        DHe22_mu = 0*mu;
    else
        error('Only the model "MomC" is implemented (model = "%s").', model);
    end
end

% -------------------------------------------------------------------------
% Model parameters (MomC)
% -------------------------------------------------------------------------
function [a_R,a_D,omega,B0,harmonic] = parameters_MomC()
    % Values are taken from case_parameters.m of the current case. They are
    % read once per run and then kept in memory ("clear all" at the start of
    % the script resets them), because this function is called very often.
    persistent par
    if isempty(par)
        P = get_case_parameters();
        par = struct('a_R',P.a_R,'a_D',P.a_D,'omega',P.omega,'B0',P.B0,'harmonic',P.harmonic);
    end
    a_R = par.a_R;
    a_D = par.a_D;
    omega = par.omega;
    B0 = par.B0;
    harmonic = par.harmonic;    % 1 = harmonic potential on, 0 = off
end
