% EHRENFEST_SIMULATION  Mixed quantum-classical Ehrenfest (mean-field) dynamics.
%
%   Propagates an ensemble of n_b classical trajectories (q,p), each carrying
%   a 2x2 electronic density matrix rho, with a fourth-order Runge-Kutta
%   scheme. The run stops at t_max or when the relative error of the total
%   energy exceeds tol_energy.
%
%   Usage:   cd into a run folder (e.g. case1/ehrenfest/n500) and run
%            ehrenfest_simulation. Parameters: ehrenfest_input.m and the
%            case_parameters.m of the case.
%
%   Output (written to the current folder):
%     ehrenfest_simulation.mat   complete workspace (trajectories, energies, ...)
%     time_seconds.dat           wall-clock run time
%     plot_energy.fig            energy components over time
%     plot_energy_error.fig      relative energy error over time

clear all; close all; clc;

%% Parameters
ehrenfest_input

tic
%% Initialisation
L = ceil(t_max/dt)+1;           % total number of time steps
q_val = zeros(n_b,L);           % trajectories, position
p_val = zeros(n_b,L);           % trajectories, momentum
rho_val = zeros(2,2,n_b,L);     % electronic density matrices
E_tot_val = zeros(1,L);         % total energy
E_kin_val = zeros(1,L);         % kinetic energy
E_mfd_val = zeros(1,L);         % mean-field energy

[w,q,p,rho] = fun_initialize();
q = rescale1*q;
p = rescale2*p;
pos_var = var(q);
mom_var = var(p);
[E_tot, E_kin, E_mfd] = fun_energy(w,q/rescale1,p/rescale2,rho);

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

    % one RK4 step (computed in unscaled variables)
    [q,p,rho] = fun_propagate(q/rescale1,p/rescale2,rho);
    q = rescale1*q;
    p = rescale2*p;

    l = l+1;
    [E_tot, E_kin, E_mfd] = fun_energy(w,q/rescale1,p/rescale2,rho);

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
xlabel('time')
legend('tot','kin','mfd','Location','northwest')
str = strcat({'Energy, Ehrenfest, N = '}, num2str(n_b),...
    {', varp = '}, num2str(round(pos_var,2)),...
    {', varm = '}, num2str(round(mom_var,2)));
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
title(['Relative energy error, Ehrenfest, N = ' num2str(n_b)])
savefig('plot_energy_error.fig')
close


%% ========================================================================
%  Local functions
%  ========================================================================

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
function [E_tot, E_kin, E_mfd] = fun_energy(w,q,p,rho)
    % Ensemble-averaged total, kinetic, mean-field and correlation energy.
    global n_b M

    Ek_tot = zeros(1,n_b);
    Ek_kin = zeros(1,n_b);
    Ek_mfd = zeros(1,n_b);
    for k = 1:n_b
        [He_fullk,~,~] = fun_He_full(q(k),p(k));

        Ek_kin(k) = p(k)^2/2/M;
        Ek_mfd(k) = trace(rho(:,:,k)'*He_fullk);
        Ek_tot(k) = Ek_kin(k)+Ek_mfd(k);
    end
    E_tot = real(sum(w.*Ek_tot));
    E_kin = real(sum(w.*Ek_kin));
    E_mfd = real(sum(w.*Ek_mfd));
end

% -------------------------------------------------------------------------
function [qdot,pdot,rhodot] = fun_dot(q,p,rho)
    % Right-hand side of the Ehrenfest equations of motion.
    global n_b M hbar
    qdot = zeros(1,n_b);
    pdot = zeros(1,n_b);
    rhodot = zeros(2,2,n_b);

    for k = 1:n_b
        [He_fullk, DHe_xk, DHe_muk] = fun_He_full(q(k),p(k));

        qdot(k) = p(k)/M+trace(rho(:,:,k)'*DHe_muk);
        pdot(k) = -trace(rho(:,:,k)'*DHe_xk);
        rhodot(:,:,k) = -1i/hbar*(He_fullk*rho(:,:,k)-rho(:,:,k)*He_fullk);
    end
end

% -------------------------------------------------------------------------
function [q,p,rho] = fun_propagate(q,p,rho)
    % One step of the classical fourth-order Runge-Kutta method.
    global dt
    q_old = q;
    p_old = p;
    rho_old = rho;

    % stage 1
    [qdot,pdot,rhodot] = fun_dot(q,p,rho);
    dq1 = dt*qdot;
    q = q_old+dq1/2;
    dp1 = dt*pdot;
    p = p_old+dp1/2;
    drho1 = dt*rhodot;
    rho = rho_old+drho1/2;

    % stage 2
    [qdot,pdot,rhodot] = fun_dot(q,p,rho);
    dq2 = dt*qdot;
    q = q_old+dq2/2;
    dp2 = dt*pdot;
    p = p_old+dp2/2;
    drho2 = dt*rhodot;
    rho = rho_old+drho2/2;

    % stage 3
    [qdot,pdot,rhodot] = fun_dot(q,p,rho);
    dq3 = dt*qdot;
    q = q_old+dq3;
    dp3 = dt*pdot;
    p = p_old+dp3;
    drho3 = dt*rhodot;
    rho = rho_old+drho3;

    % stage 4 and update
    [qdot,pdot,rhodot] = fun_dot(q,p,rho);
    dq4 = dt*qdot;
    q = q_old+dq1/6+dq2/3+dq3/3+dq4/6;
    dp4 = dt*pdot;
    p = p_old+dp1/6+dp2/3+dp3/3+dp4/6;
    drho4 = dt*rhodot;
    rho = rho_old+drho1/6+drho2/3+drho3/3+drho4/6;
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
