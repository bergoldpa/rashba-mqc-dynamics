% MOMENTUM_SURFACES  Energy surfaces of the MomC model as functions of momentum.
%
%   For every test case (case1, case2, ...) the parameters are read
%   automatically from its case_parameters.m: mass M, couplings a_R and
%   a_D, field B0, the momentum scaling rescale and the energy unit of the
%   output (P.pes_energy_unit). The script diagonalises
%
%       H(p) = p^2/(2M) + a_R*p*sigma_y + B0*sigma_x
%
%   on a momentum grid and plots the two eigenvalues against the momentum
%   in eV/c. The Dresselhaus term a_D*p*sigma_x is deliberately switched off
%   (factor 0): the figures show the surfaces without it. This differs from
%   the simulations, which include a_D (only relevant for case 1, the only
%   case with a_D ~= 0).
%
%   Usage:   run from anywhere (core must be on the MATLAB path).
%
%   Output:  figures_momentum/MOM<case>_1.dat (lower surface) and
%            figures_momentum/MOM<case>_2.dat (upper surface), columns:
%            momentum (eV/c), energy (P.pes_energy_unit). These files are
%            the input of the LaTeX figures in figures_momentum/.

clear all; close all; clc

%% Model
global model M rescale P_case
model = "MomC";

%% Folders
here = fileparts(mfilename('fullpath'));    % potential_energy_surfaces
root = fileparts(here);                     % CODE
out_dir = fullfile(here, 'figures_momentum');

%% Unit conversion
h_a0_in_eV_c = 3.7289e+03;                  % 1 au of momentum in eV/c
E_unit.eV   = 27.211386;                    % 1 au of energy in eV
E_unit.meV  = 27.211386*1e+3;               % ... in meV
E_unit.mueV = 27.211386*1e+6;               % ... in micro-eV

%% Momentum grid (rescaled units)
p = -20:.01:20;

%% Loop over all test cases
cases = dir(fullfile(root, 'case*'));
for c = 1:numel(cases)
    case_dir = fullfile(cases(c).folder, cases(c).name);
    if ~cases(c).isdir || ~isfile(fullfile(case_dir, 'case_parameters.m'))
        continue
    end
    n = sscanf(cases(c).name, 'case%d');

    % parameters of this case
    P_case = get_case_parameters(case_dir);
    M = P_case.M;
    rescale = P_case.rescale;
    unit = char(P_case.pes_energy_unit);

    % eigenvalues on the momentum grid
    eig1 = zeros(1,length(p));
    eig2 = zeros(1,length(p));
    for k = 1:length(p)
        E = eig(fun_He_mom(p(k)));      % eigenvalues in ascending order
        eig1(k) = E(1);
        eig2(k) = E(2);
    end

    % physical units
    p_value = 1/rescale*h_a0_in_eV_c;       % eV/c per rescaled momentum unit
    p_phys = p*p_value;
    E_fac = E_unit.(unit);

    % plot
    figure
    plot(p_phys,eig1*E_fac,'-k')
    hold on
    plot(p_phys,eig2*E_fac,'-k')
    grid on
    xlabel('Momentum (eV/c)')
    ylabel(['Energy (' unit ')'])
    title(sprintf('Eigenvalues of the Hamiltonian, case %d', n))

    % export
    out1 = [p_phys' (eig1*E_fac)'];
    out2 = [p_phys' (eig2*E_fac)'];
    writematrix(out1, fullfile(out_dir, sprintf('MOM%d_1.dat', n)), 'Delimiter',' ', 'FileType','text');
    writematrix(out2, fullfile(out_dir, sprintf('MOM%d_2.dat', n)), 'Delimiter',' ', 'FileType','text');
    fprintf('case %d: MOM%d_1.dat, MOM%d_2.dat written (energy in %s)\n', n, n, n, unit);
end


%% ========================================================================
%  Local functions
%  ========================================================================

% -------------------------------------------------------------------------
function He = fun_He_mom(p)
    % Hamiltonian of the MomC model at momentum p (in rescaled units).
    global model rescale
    p_sca = p/rescale;

    if (model == "MomC")
        global M
        [a_R,a_D,~,B0] = parameters_MomC();
        sigma0 = eye(2);
        sigma1 = [0,1;1,0];
        sigma2 = [0,-1i;1i,0];
        % Dresselhaus term deliberately switched off (factor 0), see header
        He = p_sca^2/(2*M)*sigma0+a_R*p_sca*sigma2+0*a_D*p_sca*sigma1+B0*sigma1;
    end
end

% -------------------------------------------------------------------------
function [a_R,a_D,omega,B0] = parameters_MomC()
    % Coupling parameters of the current case (from case_parameters.m).
    global P_case
    a_R = P_case.a_R;
    a_D = P_case.a_D;
    omega = P_case.omega;
    B0 = P_case.B0;
end
