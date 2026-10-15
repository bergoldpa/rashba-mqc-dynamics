% CASE_PARAMETERS  Parameters of test case 1.
%
%   Read by core/get_case_parameters.m. Only these values differ between the
%   test cases; all code is shared and lives in CODE/core.

%% Model (MomC: momentum-coupled two-level system), in atomic units
%   Material InSb: hbar*alpha_R = 3e-10 eV*cm, g = -51, B = 140 mT, hbar*omega = 13 mueV.
%   The values below are these inputs in atomic units, rounded as in the
%   computations for the paper.
P.M        = 0.014;              % effective mass (in units of the electron mass)
P.a_R      = 0.0021;             % Rashba coupling
P.a_D      = 0.0053;             % Dresselhaus coupling (not part of Hamiltonian (1.10) in the paper)
P.omega    = 4.7746e-07;         % trap frequency (also sets the length scale)
P.B0       = -7.5913e-06;        % Zeeman term g*mu_B*B/2 along sigma_x
P.harmonic = 0;                  % 1 = harmonic potential on, 0 = off

%% Initial condition and scaling
P.p_init   = 0;                 % initial mean momentum (p0 in the SOFT code)
P.rescale  = 1/sqrt(P.M*P.omega);   % length scale; rescale2 = rescale, rescale1 = 1/rescale
P.gamma0   = 1i*P.M*P.omega;        % width parameter of the initial Gaussian (oscillator ground state)
P.n_b      = 500;               % number of trajectories (override per folder with local_parameters.m)

%% Time integration
P.dt       = 1e5/4;             % time step (Ehrenfest and koopmons)
P.t_max    = 1e7;               % final time (t_final in the SOFT code)

%% Numerical settings
P.alpha      = 0.5;             % koopmon kernel width (rescaled units)
P.soft_t_steps = 1e3;           % SOFT: number of time steps
P.soft_nx    = 2^12;            % SOFT: number of grid points

%% SOFT simulation box
P.soft_a   = -3e5;              % box [a, -a] before rescaling

%% Postprocessing
P.xlim       = [-10,10];        % phase-space plots: position range
P.ylim       = [-1,1];          % phase-space plots: momentum range
P.ps_xstep   = 5;               % phase-space plots: tick/grid spacing in position
P.ps_ystep   = 0.2;             % phase-space plots: tick/grid spacing in momentum
P.time_final = .9e7;            % last time shown in the phase-space plots
P.time_scale = "PICO";          % time axis of the exported .dat files: "PICO" (ps) or "NANO" (ns)

%% Momentum surfaces (potential_energy_surfaces/momentum_surfaces.m)
P.pes_energy_unit = "meV";         % energy unit of the exported surfaces: "eV", "meV" or "mueV"
