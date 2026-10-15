% PLOT_BLOCH_PURITY_SOFT  Bloch vector and purity of the exact (SOFT) solution.
%
%   Computes the Bloch vector n = (n_x, n_y, n_z) of the reduced electronic
%   density matrix and the purity (1+|n|^2)/2 over time.
%
%   Usage:   run from caseN/postprocessing/<topic>. The result is read
%            directly from caseN/soft/soft_simulation.mat (not modified).
%
%   Output (folder output/):  <prefix>_bloch_*.fig, <prefix>_purity.fig and
%            soft_n1.dat, soft_n2.dat, soft_n3.dat, soft_purity.dat (time, value).

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

dt = tau;   % SOFT time step
prefix = 'soft';    % prefix of all output files

%% Bloch vector n = (2 Re rho_12, 2 Im rho_21, rho_11 - rho_22) and its length
n = zeros(3,n_snap);
norm_n = zeros(1,n_snap);
for l = 1:n_snap
    rho_l = rho_val(:,:,l);
    n(1,l) = 2*real(rho_l(1,2));
    n(2,l) = 2*imag(rho_l(2,1));
    n(3,l) = real(rho_l(1,1)-rho_l(2,2));
    norm_n(l) = sqrt(n(1,l)^2+n(2,l)^2+n(3,l)^2);
end

%% Plot: all components and purity
t_plot = (0:(l-1))*dt;              % time axis (atomic units)
labels = {'n_x','n_y','n_z'};
figure
for c = 1:3
    subplot(4,2,2*c-1:2*c);
    plot(t_plot,n(c,:))
    xlabel('time')
    ylabel(labels{c})
    grid on
    ylim([-1.05,1.05])
end

subplot(4,2,7:8);
plot(t_plot,.5*(1+norm_n.^2))
xlabel('time')
ylabel('Purity')
grid on
ylim([0.45,1.05])
sgtitle('Components of the Bloch Sphere')
saveas(gcf, fullfile("output", [prefix '_bloch_components.fig']));
close

%% Plot and export the components n_x, n_y, n_z (files <prefix>_n1.dat ... _n3.dat)
for c = 1:3
    figure
    plot(t_plot,n(c,:))
    xlabel('time')
    ylabel(labels{c})
    grid on
    ylim([-1.05,1.05])
    saveas(gcf, fullfile("output", sprintf('%s_bloch_n%d.fig', prefix, c)));
    close

    out = [t_plot' n(c,:)'];
    writematrix(out, fullfile("output", sprintf('%s_n%d.dat', prefix, c)), 'Delimiter',' ', 'FileType','text');
end

%% Plot and export purity (file <prefix>_purity.dat)
figure
plot(t_plot*2.4188843265857e-17*1e12,.5*(1+norm_n.^2))   % time in ps
xlabel('time (in ps)')
ylabel('Purity')
grid on
ylim([0.45,1.05])
saveas(gcf, fullfile("output", [prefix '_purity.fig']));
close

out = [t_plot' (.5*(1+norm_n.^2))'];
writematrix(out, fullfile("output", [prefix '_purity.dat']), 'Delimiter',' ', 'FileType','text');

%% Plot: Bloch vector trajectory on the sphere
figure
x_plot = n(1,:);
y_plot = n(2,:);
z_plot = n(3,:);
time_plot = (0:length(x_plot)-1)*dt;
surface([x_plot;x_plot],[y_plot;y_plot],[z_plot;z_plot],[time_plot;time_plot],...
    'EdgeColor', 'flat', 'FaceColor', 'none', 'linewidth', 4);
xlabel('n_x')
ylabel('n_y')
zlabel('n_z')

hold on
[X,Y,Z] = sphere;
surf(X,Y,Z, 'FaceAlpha', .05, 'EdgeColor', 'white');
axis equal
grid on
colorbar
h = colorbar;
ylabel(h,'time')
view([30 15])

title('Time Evolution of the Bloch Vector')
saveas(gcf, fullfile("output", [prefix '_bloch_full.fig']));
close
