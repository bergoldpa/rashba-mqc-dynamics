% CASE_PARAMETERS  Parameters of test case 5.
%
%   Read by core/get_case_parameters.m. Only these values differ between the
%   test cases; all code is shared and lives in CODE/core.

%% Model (MomC: momentum-coupled two-level system), in atomic units
%   Material GaAs: hbar*alpha_R = 0.68e-11 eV*cm, g = -0.44, B = -4.2 mT, hbar*omega = 0.182 mueV.
%   The values below are these inputs in atomic units, rounded as in the
%   computations for the paper.
P.M        = 0.067;              % effective mass (in units of the electron mass)
P.a_R      = 4.7159e-05;         % Rashba coupling
P.a_D      = 0;                  % Dresselhaus coupling
P.omega    = 6.6844e-09;         % trap frequency (also sets the length scale)
P.B0       = 1.9659e-09;         % Zeeman term g*mu_B*B/2 along sigma_x
P.harmonic = 1;                  % 1 = harmonic potential on, 0 = off

%% Initial condition and scaling
P.p_init   = 1e-4;              % initial mean momentum (p0 in the SOFT code)
P.rescale  = 1/sqrt(P.M*P.omega);   % length scale; rescale2 = rescale, rescale1 = 1/rescale
P.gamma0   = 1i*P.M*P.omega;        % width parameter of the initial Gaussian (oscillator ground state)
P.n_b      = 500;               % number of trajectories (override per folder with local_parameters.m)

%% Time integration
P.dt       = 9e7/4;             % time step (Ehrenfest and koopmons)
P.t_max    = 9e9;               % final time (t_final in the SOFT code)

%% Numerical settings
P.alpha      = 0.5;             % koopmon kernel width (rescaled units)
P.soft_t_steps = 1e3;           % SOFT: number of time steps
P.soft_nx    = 2^12;            % SOFT: number of grid points

%% SOFT simulation box
P.soft_a   = -5e5;              % box [a, -a] before rescaling

%% Postprocessing
P.xlim       = [-17.5,15];      % phase-space plots: position range
P.ylim       = [-.6,.6];        % phase-space plots: momentum range
P.ps_xstep   = 5;               % phase-space plots: tick/grid spacing in position
P.ps_ystep   = 0.2;             % phase-space plots: tick/grid spacing in momentum
P.time_final = 8.9e9;           % last time shown in the phase-space plots
P.time_scale = "NANO";          % time axis of the exported .dat files: "PICO" (ps) or "NANO" (ns)

%% Momentum surfaces (potential_energy_surfaces/momentum_surfaces.m)
P.pes_energy_unit = "mueV";        % energy unit of the exported surfaces: "eV", "meV" or "mueV"
