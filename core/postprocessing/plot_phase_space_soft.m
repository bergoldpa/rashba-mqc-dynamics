% PLOT_PHASE_SPACE_SOFT  Wigner function of the exact (SOFT) solution.
%
%   Computes the Wigner function W(q,p) of the two-component wave function
%   (sum of both components) at the snapshot times defined in
%   phase_space_plot_settings.m and plots it.
%
%   Usage:   run from caseN/postprocessing/<topic>. The result is read
%            directly from caseN/soft/soft_simulation.mat (not modified).
%
%   Output:  ps1_S.png ... ps8_S.png (300 dpi) and WIGNER_VAL_REF.mat (colour scale and
%            settings, needed by plot_phase_space_ehrenfest_koopmons.m; a few MB).
%            The Wigner functions stay in the workspace (wigner_fft_val) but
%            are not saved - rerun this script to get them back exactly.

% Clear the workspace, but keep run_case_file / run_case_method if this
% script is called from run_case.m
clearvars -global
clearvars -except run_case_file run_case_method
close all;
string_name = "../../soft/soft_simulation.mat";
load(string_name)

global rescale

%% Plot settings, part 1 (axes and snapshot times)
loop = 1;
phase_space_plot_settings

%% Wigner function at every snapshot time
l_times = length(times_all);
wig_Nx = nx;
wig_Np = nx;
wigner_fft_val = zeros(wig_Np,wig_Nx,l_times);
for k = 1:l_times
    time_val = times_all(k);
    y1 = psi1_val(:,time_val).';
    y2 = psi2_val(:,time_val).';

    % norm of the position density
    density_pos = abs(y1).^2+abs(y2).^2;
    norm_pos = trapz(x,density_pos);

    % norm of the momentum density (p is in FFT order, so it is sorted with
    % fftshift before integrating)
    p = 2*pi*hbar/(b-a)*[0:nx/2-1 -nx/2:-1];
    xi1 = exp(-1i/hbar*p*a)*(b-a)/nx/sqrt(2*pi*hbar).*fft(y1);
    xi2 = exp(-1i/hbar*p*a)*(b-a)/nx/sqrt(2*pi*hbar).*fft(y2);
    density_mom = abs(xi1).^2+abs(xi2).^2;
    norm_mom = trapz(fftshift(p),fftshift(density_mom));

    fprintf('snapshot %d of %d: norm in position %.6f, in momentum %.6f\n', k, l_times, norm_pos, norm_mom);

    % Wigner function of both components
    W_qp1 = mywigner_fft(y1,x,wig_Np);
    W_qp2 = mywigner_fft(y2,x,wig_Np);
    W_qp = W_qp1+W_qp2;
    wigner_fft_val(:,:,k) = W_qp;
end

%% Plot settings, part 2 (colour scale from the Wigner functions) and plots
SOFT_yn = 1;
loop = 2;
phase_space_plot_settings
for k = 1:ANZ
    str = ['ps' num2str(k) '_S.png'];
    plot_wigner_fft(k,wigner_fft_val,x,wig_Np);
    exportgraphics(gcf,str,'Resolution',300,'BackgroundColor','white');
    close
end

%% Save colour scale and settings for the Ehrenfest/koopmon plots
% The large arrays (Wigner functions and wave functions) are not saved: they
% are recomputed exactly by running this script on soft_simulation.mat.
save("WIGNER_VAL_REF.mat",'-regexp','^(?!(wigner_fft_val|W_qp|W_qp1|W_qp2|psi1_val|psi2_val)$).')

%% ========================================================================
%  Local functions
%  ========================================================================

% -------------------------------------------------------------------------
function W_qp = mywigner_fft(psi,x,wig_Np)
    % Wigner function W(q,p) of the wave function psi on the grid x,
    % computed with an FFT along the anti-diagonals of psi(x)*conj(psi(x')).
    global hbar nx dx

    s = wig_Np/nx;
    L = nx*dx;
    p = -pi/s/L*[0:wig_Np/2-1 -wig_Np/2:-1]';

    rho_x = psi.'*conj(psi);
    rho_rot = antidiagonals_reshape(rho_x);
    rho_rot_ft = fft(rho_rot,wig_Np,1);
    W_qp = exp(2i/hbar*p*x(1))*dx/pi/hbar.*rho_rot_ft;

    [~,idp] = sort(p);
    W_qp = real(W_qp(idp,:));
end

% -------------------------------------------------------------------------
function plot_wigner_fft(k,wigner_fft_val,x,wig_Np)
    % Plot the Wigner function of snapshot k (physical units, fixed colour scale).
    global nx dx xlim_var ylim_var xticks_var yticks_var ticks_all ticks_disp colorMapMatrix rescale times_disp

    s = wig_Np/nx;
    L = nx*dx;

    W_qp = wigner_fft_val(:,:,k);
    q_w = x;
    p_w = (-wig_Np/2:wig_Np/2-1)*pi/L/s;

    % conversion from atomic units to micrometre (position) and eV/c (momentum)
    a0_in_m = 5.291772e-11;     % 1 au of length in m
    h_a0_in_eV_c = 3.7289e+03;  % 1 au of momentum in eV/c

    q_scale = rescale*a0_in_m*1e6;          % position axis in micrometres

    p_scale = 1/rescale*h_a0_in_eV_c;       % momentum axis in eV/c

    % only the part of the grid inside the plot window (plus 2 grid points on
    % each side) is drawn: the picture is identical, but the vector PDF is
    % much smaller and faster to create
    qv = q_w*q_scale;
    pv = p_w*p_scale;
    iq = window_indices(qv, xlim_var);
    ip = window_indices(pv, ylim_var);
    [fig_q,fig_p] = meshgrid(qv(iq),pv(ip));

    % separate (undocked) window of fixed size: since R2025a figures are
    % docked by default, and the automatic ticks depend on the figure size
    figure('WindowStyle','normal','Units','pixels','Position',[0,0,700,600]);

    % filled contour plot of the Wigner function
    % (values below the lowest colour level are shown in the lowest colour)
    contourf(fig_q,fig_p,max(W_qp(ip,iq),ticks_all(1)),ticks_all,LineColor = 'none')
    grid on
    set(gca,'layer','top');

    hold on
    clim([ticks_disp(1),ticks_disp(end)])
    xlim(xlim_var)
    ylim(ylim_var)
    xticks(xticks_var)
    yticks(yticks_var)

    colormap(colorMapMatrix);
    CB = colorbar;
    set(CB, 'Ticks', ticks_disp);
    fontsize(gcf,scale = 2)
    set(gcf,'position',[0,0,700,600])
    % print the snapshot time
    P = get_case_parameters();
    if P.time_scale == "NANO"
        fprintf('%.2f ns\n', round(times_disp(k)*2.4188843265857e-17*1e9,2));
    else
        fprintf('%.2f ps\n', round(times_disp(k)*2.4188843265857e-17*1e12,2));
    end

    if k < 8
        set(gca,'XTickLabel',[]);
    end
end

% -------------------------------------------------------------------------
function idx = window_indices(v, lims)
    % Indices of the grid points v inside [lims(1), lims(2)], extended by
    % 2 points on each side (all points if none lies inside).
    idx = find(v >= lims(1) & v <= lims(2));
    if isempty(idx)
        idx = 1:numel(v);
        return
    end
    idx = max(idx(1)-2,1):min(idx(end)+2,numel(v));
end

% -------------------------------------------------------------------------
function R = antidiagonals_reshape(A)
    % Rearrange the anti-diagonals of the square matrix A into the columns of R
    % (centred), as needed for the FFT in mywigner_fft.
    n = size(A,1);
    R = zeros(n,n);
    if mod(n,2) == 1
        for k = 1:n
            middle = floor(n/2)+1;
            if middle+(k-1) <= n
                range = middle-(k-1):middle+(k-1);
            else
                range = range(2:end-1);
            end
            antidiag = diag(flipud(A),-n+2*(k-1)+1);
            R(range,k) = antidiag;
        end
    else
        for k = 1:n
            middle = floor(n/2);
            if middle-(k-1) >= 1
                range = middle-(k-1):middle+(k-1);
            elseif middle-(k-1) == 0
            elseif middle-(k-1) < 0
                range = range(2:end-1);
            end
            antidiag = diag(flipud(A),-n+2*(k-1)+1);
            R(range,k) = antidiag;
        end
    end
end
