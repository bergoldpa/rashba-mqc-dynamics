% SOFT_SIMULATION  Exact quantum reference solution (split-operator Fourier transform).
%
%   Solves the time-dependent Schroedinger equation for a two-component
%   wave function psi = (psi{1}, psi{2}) with Hamiltonian
%
%       H = p^2/(2m) + a_R*p*sigma_y + a_D*p*sigma_x + V(x),
%       V(x) = harmonic * m*omega^2*x^2/2 * I + B0*sigma_x,
%
%   using a symmetric Strang splitting in every time step (half kinetic
%   step - potential step - half kinetic step): kinetic and momentum-coupling
%   terms are applied in Fourier space, the potential in position space.
%   Snapshots are stored at t = 0, tau, ..., t_final (n_snap = t_steps+1).
%   Positions are rescaled by rescale = 1/sqrt(m*omega).
%
%   Usage:   cd into the SOFT folder of a case (e.g. case1/soft) and run
%            soft_simulation. Parameters: case_parameters.m of the case.
%            The norm of the initial state is printed (warning if it is not 1).
%
%   Output:  soft_simulation.mat (complete workspace) in the current folder.
%            Densities are not stored; they follow from psi1_val, psi2_val.

clear all; close all; clc
clc

%% Model and physical constants
global model
model = "MomC";
P = get_case_parameters();  % case-specific values (case_parameters.m)

global m hbar
m = P.M;        % effective electron mass (units of the electron mass)
hbar = 1;       % reduced Planck constant

%% Simulation box and grids
global a b nx dx np dp
a = P.soft_a;   % box [a, b] before rescaling
b = -a;

[~,~,omega,~] = parameters_MomC();
global rescale
rescale = 1/sqrt(m*omega);
a = a/rescale;
b = b/rescale;

L = b-a;
nx = P.soft_nx;  % number of grid points in position space (case_parameters.m)
dx = L/nx; % grid spacing x
x = a+(1:nx)*dx;

np = nx; % number of grid points in momentum space
dp = 2*pi/L; % grid spacing p
% momentum grid in FFT ordering
p = zeros(1,np);
for j = 1:np
    if (j < (np/2+2))
        p(j) = (j-1)*dp;
    else
        p(j) = (j-1-np)*dp;
    end
end

global KE ME
KE = p.^2/2/m/rescale^2;
ME = p/rescale;

%% Initial wave packet (Gaussian)
gz = @(q0,p0,gamma,hbar,x) (pi*hbar)^(-.25)*imag(gamma)^.25*...
    exp(1i/hbar*(gamma*(x-q0).^2/2+p0*(x-q0)));

q0 = 0;
p0 = P.p_init;
[~,~,omega,~] = parameters_MomC();
gamma0 = 1i*m*omega;

% initial spin along sigma3
psi{1} = 1*gz(q0,p0,gamma0,hbar,rescale*x)*sqrt(rescale);
psi{2} = 0*gz(q0,p0,gamma0,hbar,rescale*x)*sqrt(rescale);

% normalisation check of the initial state
abs0 = abs(psi{1}.^2)+abs(psi{2}.^2);
norm0 = trapz(x,abs0);
fprintf('Initial norm: %.12f\n', norm0);
if abs(norm0-1) > 1e-6
    warning('Initial wave packet is not normalised (norm = %.6g).', norm0);
end


%% Time stepping and storage
global t_steps tau
t_final = P.t_max;
t_steps = P.soft_t_steps;   % number of time steps (case_parameters.m)
tau = t_final/t_steps;

% snapshots at t = 0, tau, ..., t_final: snapshot l is at time (l-1)*tau
n_snap = t_steps+1;
rho_val = zeros(2,2,n_snap);
psi1_val = zeros(length(x),n_snap);
psi2_val = zeros(length(x),n_snap);

%% Propagators (the Hamiltonian does not depend on time: computed once)
[M0,MC] = mix(x);                 % potential, full step tau
KIN = kinetic_propagators();        % kinetic and momentum couplings, half step tau/2

%% Time propagation
tic;
for l = 1:n_snap
    % one symmetric split-operator step from (l-2)*tau to (l-1)*tau
    if (l > 1)
        psi = trotter(psi,M0,MC,KIN);
    end

    % reduced electronic density matrix rho = int psi psi^dagger dx
    rho_11 = abs(psi{1}).^2;
    rho_22 = abs(psi{2}).^2;
    rho_12 = psi{1}.*conj(psi{2});
    rho_21 = psi{2}.*conj(psi{1});

    rho_11_int = trapz(x,rho_11);
    rho_22_int = trapz(x,rho_22);
    rho_12_int = trapz(x,rho_12);
    rho_21_int = trapz(x,rho_21);
    rho = [rho_11_int rho_12_int; rho_21_int rho_22_int];
    rho_val(:,:,l) = rho;

    % store the wave function
    psi1_val(:,l) = psi{1};
    psi2_val(:,l) = psi{2};

    % progress (every 10 %)
    if mod(l-1, round(t_steps/10)) == 0
        fprintf('%3.0f %%\n', 100*(l-1)/t_steps);
    end
end
toc;

% norm check of the final state
norm_end = trapz(x,abs(psi{1}).^2+abs(psi{2}).^2);
fprintf('Final norm: %.12f\n', norm_end);
if abs(norm_end-1) > 1e-6
    warning('Norm is not conserved (final norm = %.6g).', norm_end);
end

time = (0:t_steps)*tau;   % times of the snapshots

%% Save results
save(mfilename,'-v7.3')
close

%% ========================================================================
%  Local functions
%  ========================================================================

% -------------------------------------------------------------------------
function [M0,Mc] = mix(x)
    % Potential propagator exp(-i*He(x)*tau) at every grid point, split into
    % diagonal entries M0{1}, M0{2} and the off-diagonal entry Mc
    % (He is real symmetric, so the propagator is symmetric: M(2,1) = M(1,2)).
    global tau nx

    M0{1} = zeros(1,nx);
    M0{2} = zeros(1,nx);
    Mc = zeros(1,nx);
    for j = 1:nx
        He = potential(x(j));
        [V,D] = eig(He);
        M = V*[exp(-1i*D(1,1)*tau),0;0,exp(-1i*D(2,2)*tau)]*V';   % exp(-i*He*tau) = V*exp(-i*D*tau)*V'
        M0{1}(j) = M(1,1);
        M0{2}(j) = M(2,2);
        Mc(j) = M(1,2);
    end
end

% -------------------------------------------------------------------------
function [psi] = trotter(psi,M0,MC,KIN)
    % One symmetric (Strang) split-operator step of length tau:
    % half kinetic step - full potential step - half kinetic step.
    psi = kinetic_half_step(psi,KIN);
    psi = potential_step(psi,M0,MC);
    psi = kinetic_half_step(psi,KIN);
end

% -------------------------------------------------------------------------
function KIN = kinetic_propagators()
    % Propagators of the kinetic part T for a half step tau/2 in momentum
    % space: p^2/(2m) and the momentum couplings a_D*p*sigma_x and
    % a_R*p*sigma_y. The two couplings do not commute; they are applied in
    % the symmetric order sigma_x (tau/4) - sigma_y (tau/2) - sigma_x (tau/4).
    global tau KE ME
    sigma1 = [0,1;1,0];
    sigma2 = [0,-1i;1i,0];
    [KIN.V1,D1] = eig(sigma1);
    [KIN.V2,D2] = eig(sigma2);

    l1_x = D1(1,1);
    l2_x = D1(2,2);
    l1_y = D2(1,1);
    l2_y = D2(2,2);

    [a_R,a_D,~,~] = parameters_MomC();
    KIN.KP = exp(-1i*KE*tau/2);                     % p^2/(2m), tau/2
    KIN.MP_1_x = exp(-l1_x*1i*a_D*ME*tau/4);        % a_D*p*sigma_x, tau/4
    KIN.MP_2_x = exp(-l2_x*1i*a_D*ME*tau/4);
    KIN.MP_1_y = exp(-l1_y*1i*a_R*ME*tau/2);        % a_R*p*sigma_y, tau/2
    KIN.MP_2_y = exp(-l2_y*1i*a_R*ME*tau/2);
end

% -------------------------------------------------------------------------
function [psi] = kinetic_half_step(psi,KIN)
    % Apply exp(-i*T*tau/2) in Fourier space with the precomputed
    % propagators KIN (see kinetic_propagators).
    global dx dp np

    % to momentum space
    ftpsi{1} = fft(psi{1})*dx/sqrt(2*pi);
    ftpsi{2} = fft(psi{2})*dx/sqrt(2*pi);

    ftpsi{1} = KIN.KP.*ftpsi{1}; % apply p^2
    ftpsi{2} = KIN.KP.*ftpsi{2}; % apply p^2

    ftpsi = apply_spin_term(ftpsi,KIN.V1,KIN.MP_1_x,KIN.MP_2_x);    % p*sigma_x, tau/4
    ftpsi = apply_spin_term(ftpsi,KIN.V2,KIN.MP_1_y,KIN.MP_2_y);    % p*sigma_y, tau/2
    ftpsi = apply_spin_term(ftpsi,KIN.V1,KIN.MP_1_x,KIN.MP_2_x);    % p*sigma_x, tau/4

    % back to position space
    psi{1} = ifft(ftpsi{1})*dp/sqrt(2*pi)*np;
    psi{2} = ifft(ftpsi{2})*dp/sqrt(2*pi)*np;
end

% -------------------------------------------------------------------------
function ftpsi = apply_spin_term(ftpsi,V,MP_1,MP_2)
    % Apply exp(-i*c*p*sigma) in momentum space: transform to the eigenbasis
    % V of sigma, multiply by the phases MP_1, MP_2, transform back.
    v = [ftpsi{1};ftpsi{2}];
    v = V'*v;
    ftpsi{1} = MP_1.*v(1,:);
    ftpsi{2} = MP_2.*v(2,:);
    v = [ftpsi{1};ftpsi{2}];
    v = V*v;
    ftpsi{1} = v(1,:);
    ftpsi{2} = v(2,:);
end

% -------------------------------------------------------------------------
function [psi] = potential_step(psi,M0,MC)
    % Apply the potential propagator exp(-i*He(x)*tau) in position space.
    for s = 1:2
        psiA{s} = psi{s}.*MC;
        psi{s} = psi{s}.*M0{s};
    end
    psi{1} = psi{1}+psiA{2};
    psi{2} = psi{2}+psiA{1};
end

% -------------------------------------------------------------------------
function He = potential(x)
    % Electronic Hamiltonian He(x) (2x2) of the MomC model.

    %% Model and physical constants
    global model
    if (model == "MomC")
        global rescale m
        [~,~,omega,B0,harmonic] = parameters_MomC();
        sigma0 = eye(2);
        sigma1 = [0,1;1,0];
        He = harmonic*.5*m*omega^2*(rescale*x).^2*sigma0+B0*sigma1;
    end
end

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
