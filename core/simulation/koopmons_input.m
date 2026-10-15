% KOOPMONS_INPUT  Input parameters for the koopmon simulation.
%
%   Called at the start of koopmons_simulation.m. All values that differ
%   between test cases are read from case_parameters.m (via
%   get_case_parameters); everything else is defined here and shared by all
%   cases. The parameters are passed to the local functions as globals.

% Case-specific values (case_parameters.m, optionally local_parameters.m)
P = get_case_parameters();

%% Physical constants and model
global M hbar model
M = P.M;            % effective electron mass (units of the electron mass)
hbar = 1;           % reduced Planck constant
model = "MomC";     % Hamiltonian (only "MomC" is implemented)

%% Initial condition
global n_b q_init p_init rho_init sigma_init_q sigma_init_p
n_b = P.n_b;            % number of koopmons
q_init = 0;             % initial mean position
p_init = P.p_init;      % initial mean momentum
rho_init = [1,0;0,0];   % initial electronic density matrix

% Rescaling of position (q -> rescale1*q) and momentum (p -> rescale2*p)
global rescale1 rescale2
rescale1 = 1/P.rescale;
rescale2 = P.rescale;

% Width of the initial Gaussian wave packet, derived from gamma0
gamma0 = P.gamma0;
sigma0 = sqrt(real(1i/2/gamma0));
Delta0 = sqrt(2)*sigma0;
Delta = Delta0;

sigma_init_q = Delta/sqrt(2);       % standard deviation of the initial distribution (position)
sigma_init_p = 1/2/sigma_init_q;    % standard deviation of the initial distribution (momentum)

%% Kernel parameters
global sigma_kernel_x sigma_kernel_mu lambda_x lambda_mu
alpha = P.alpha/rescale1;           % kernel width (P.alpha in rescaled units, case_parameters.m)
sigma_kernel_x = alpha/sqrt(2);     % standard deviation of the Gaussian kernel (position)
sigma_kernel_mu = alpha/sqrt(2);    % standard deviation of the Gaussian kernel (momentum)
lambda_x = 4*sigma_kernel_x;        % kernel cut-off (position)
lambda_mu = 4*sigma_kernel_mu;      % kernel cut-off (momentum)

%% Integration grid and time stepping
global dx dmu nx nmu dt t_max
dx = sigma_kernel_x/2;      % grid spacing of the trapezoidal rule (position)
dmu = sigma_kernel_mu/2;    % grid spacing of the trapezoidal rule (momentum)
nx = 2;     % grid covers [min(q)-nx*sigma_kernel_x, max(q)+nx*sigma_kernel_x]
nmu = 2;    % grid covers [min(p)-nmu*sigma_kernel_mu, max(p)+nmu*sigma_kernel_mu]
dt = P.dt;          % time step
t_max = P.t_max;    % final time

%% Tolerances
global tol_energy
tol_energy = .1;    % stop if the relative energy error exceeds this value
