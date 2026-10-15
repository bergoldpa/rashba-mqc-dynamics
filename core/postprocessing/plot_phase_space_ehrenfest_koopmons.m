% PLOT_PHASE_SPACE_EHRENFEST_KOOPMONS  Phase-space snapshots of Ehrenfest or koopmon results.
%
%   Plots the particle positions together with a kernel-smoothed phase-space
%   density at the 8 snapshot times defined in phase_space_plot_settings.m.
%   The colour scale is taken from WIGNER_VAL_REF.mat, which is created by
%   plot_phase_space_soft.m - run that script first.
%
%   Usage:   run from caseN/postprocessing/<topic>; select the Ehrenfest or
%            koopmon result in `string_name` below. The result is read
%            directly from the simulation folder of the case (not modified).
%            The method (K = koopmons, E = Ehrenfest) follows from the file;
%            the smoothing width is set by alpha_opt below.
%
%   Output:  ps1_<method>.png ... ps8_<method>.png (300 dpi) in the current folder.

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

load("WIGNER_VAL_REF.mat",...
    "numColors","num_blue","ticks_all","ticks_disp","rescale")

tau = dt;   % time step name used in phase_space_plot_settings

%% Plot settings (axes, snapshot times, colour map)
loop = 1;
phase_space_plot_settings

SOFT_yn = 0;
loop = 2;
phase_space_plot_settings

%% Method and smoothing
% method follows from the loaded file: K = koopmons, E = Ehrenfest
if contains(string_name, "koopmons")
    plot_title = "K";
else
    plot_title = "E";
end
alpha_opt = 0.2;    % width Delta of the smoothed phase-space density (rescaled units; paper: Delta = 0.2)

%% Plot grid in physical units
global q_scale p_scale
% conversion from atomic units to micrometre (position) and eV/c (momentum)
a0_in_m = 5.291772e-11;     % 1 au of length in m
h_a0_in_eV_c = 3.7289e+03;  % 1 au of momentum in eV/c

q_scale = rescale*a0_in_m*1e6;          % position axis in micrometres

p_scale = 1/rescale*h_a0_in_eV_c;       % momentum axis in eV/c

fig_q_min = xlim_var(1);
fig_q_max = xlim_var(2);
fig_p_min = ylim_var(1);
fig_p_max = ylim_var(2);
fig_dx = (fig_q_max-fig_q_min)/100;
fig_dmu = (fig_p_max-fig_p_min)/100;
fig_grid_x = fig_q_min:fig_dx:fig_q_max;
fig_grid_mu = fig_p_min:fig_dmu:fig_p_max;
[fig_X,fig_Mu] = meshgrid(fig_grid_x/q_scale,fig_grid_mu/p_scale);

%% Snapshots 1-ANZ (one PDF each, ANZ = 8 from phase_space_plot_settings)

for k = 1:ANZ
    plot_wigner(k,q_val,p_val,n_b,w,alpha_opt,fig_X,fig_Mu)
    str = sprintf('ps%d_%s.png', k, plot_title);    % e.g. ps3_K.png
    exportgraphics(gcf,str,'Resolution',300,'BackgroundColor','white');
    close
end

%% ========================================================================
%  Local functions
%  ========================================================================

% -------------------------------------------------------------------------
function plot_wigner(k,q_val,p_val,n_b,w,alpha_opt,fig_X,fig_Mu)
    % Plot snapshot k: particles and kernel-smoothed density in phase space.
    global xlim_var ylim_var xticks_var yticks_var times_all ticks_all ticks_disp colorMapMatrix q_scale p_scale

    time_val = times_all(k);
    q = q_val(:,time_val);
    p = p_val(:,time_val);

    % kernel-smoothed phase-space density on the plot grid
    fig_psi = 0;
    for j = 1:n_b
        fig_psi = fig_psi...
            +w(j)*fun_K0(alpha_opt,(fig_X-q(j))*1)...
            .*fun_K0(alpha_opt,(fig_Mu-p(j))*1);
    end

    % particles
    % separate (undocked) window of fixed size: since R2025a figures are
    % docked by default, and the automatic ticks depend on the figure size
    figure('WindowStyle','normal','Units','pixels','Position',[0,0,700,600]);

    hold on
    scatter(q_val(:,time_val)*q_scale,p_val(:,time_val)*p_scale,...
    'MarkerEdgeColor',[0 0 0],...
    'MarkerEdgeAlpha',1,...
    'MarkerFaceColor',[0 0 0],...
    'MarkerFaceAlpha',1,...
    'LineWidth',.1,...
    'SizeData',5);

    % smoothed density
    contourf(fig_X*q_scale,fig_Mu*p_scale,fig_psi,ticks_all,LineColor='none',FaceAlpha=.75);

    grid on
    set(gca,'layer','top');

    hold on

    % axes and colour scale
    clim([ticks_disp(1),ticks_disp(end)])
    xlim(xlim_var)
    ylim(ylim_var)
    xticks(xticks_var)
    yticks(yticks_var)

    set(gca,'YTickLabel',[]);
    if k < 8
        set(gca,'XTickLabel',[]);
    end

    colormap(colorMapMatrix);
    fontsize(gcf,scale = 2)
    set(gcf,'position',[0,0,700,600])
end

% -------------------------------------------------------------------------
function kernel = fun_K0(sigma_kernel,y)
    % Normalised one-dimensional Gaussian kernel with standard deviation sigma_kernel.
    kernel = exp(-.5*(y/sigma_kernel).^2)/sqrt(2*pi)/sigma_kernel;
end
