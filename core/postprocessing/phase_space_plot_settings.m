% PHASE_SPACE_PLOT_SETTINGS  Shared settings of the phase-space plots.
%
%   Script called by plot_phase_space_soft.m and
%   plot_phase_space_ehrenfest_koopmons.m (twice, controlled by `loop`):
%     loop = 1   axis ranges (from case_parameters.m) and the 8 snapshot times
%     loop = 2   colour map; if SOFT_yn = 1 the colour levels are computed
%                from the SOFT Wigner functions, otherwise they are loaded
%                from WIGNER_VAL_REF.mat by the calling script.

global xlim_var ylim_var xticks_var yticks_var times_all times_disp ticks_all ticks_disp colorMapMatrix

if loop == 1
    P = get_case_parameters();  % axis ranges and final time from case_parameters.m
    xlim_var = P.xlim;
    ylim_var = P.ylim;

    % fixed tick positions (and thus grid lines), independent of figure and
    % font size: all multiples of P.ps_xstep / P.ps_ystep inside the axis
    % ranges (rounded to remove floating-point errors such as -3*0.2)
    xticks_var = round((ceil(xlim_var(1)/P.ps_xstep-1e-9):floor(xlim_var(2)/P.ps_xstep+1e-9))*P.ps_xstep, 10);
    yticks_var = round((ceil(ylim_var(1)/P.ps_ystep-1e-9):floor(ylim_var(2)/P.ps_ystep+1e-9))*P.ps_ystep, 10);

    time_final = P.time_final;
    ANZ = 8;    % number of snapshots
    % snapshot times and the corresponding time-step indices
    times_disp = floor(linspace(0,time_final,ANZ));
    times_all = [];
    for k = 1:ANZ
        times_all = [times_all, floor(times_disp(k)/tau)+1];
    end
end
if loop == 2
    ind = 6;            % controls the number of colours in the colour bar

    if (SOFT_yn == 1)
        % colour levels: rounded range of the Wigner function, split such that
        % zero lies on a level boundary
        W_max = max(max(max(wigner_fft_val(:,:,:))));
        W_min = min(min(min(wigner_fft_val(:,:,:))));

        d1_min = ceil(W_min*10)-1;
        z_min = d1_min/10;
        d1_max = floor(W_max*10)+1;
        z_max = d1_max/10;
        if isfield(P,'ps_clim')         % fixed colour range from case_parameters.m
            d1_min = round(P.ps_clim(1)*10);
            d1_max = round(P.ps_clim(2)*10);
            z_min = d1_min/10;
            z_max = d1_max/10;
        end

        r = -(d1_max-d1_min)/d1_min;
        gcd_get = gcd(d1_max-d1_min,d1_min);

        m = (d1_max-d1_min)/gcd_get;
        l = -d1_min/gcd_get;

        k = 1:20;
        list_k = mod(k*m+floor(m/2),l);
        k_new = k(list_k == 0);
        numColors = (k_new(ind)*m+floor(m/2))/l;
        num_blue = floor(numColors/r);

        ticks_all = linspace(z_min,z_max,numColors+1);
        ticks_disp = z_min:.1:z_max;
    end

    % blue (negative) - white - red (positive) colour map
    redPalette = linspace(.9,0, numColors-num_blue-1)';
    redColors = [ones(numColors-num_blue-1,1),redPalette,redPalette];
    whiteColor = [1, 1, 1];
    bluePalette = linspace(0, .9, num_blue)';
    blueColors = [bluePalette, bluePalette, ones(num_blue,1)];
    colorMapMatrix = [blueColors; whiteColor; redColors];
end
