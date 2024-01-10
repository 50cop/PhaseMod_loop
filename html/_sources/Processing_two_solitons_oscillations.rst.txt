.. _Processing_two_solitons_oscillations:



Processing_two_solitons_oscillations.m
======================================

.. contents:: Table of Contents


Generating the matrix containing each runs
""""""""""""""""""""""""""""""""""""""""""

::

    %% Generating the matrix containing each runs
    lst = who('data_*');
    for dd = 1:(length(lst))
        matrice_STD(dd, :, :) = eval([num2str(cell2mat(lst(dd))) '.z']);
    end
    
    lst = who('dat_sol*');
    for dd = 1:(length(lst))
        matrice_STD_sol(dd, :, :) = eval([num2str(cell2mat(lst(dd))) '.z']);
    end
    
    t = (data_01.x - mean(data_01.x)).*1e9;
    z = data_01.y;
    
    box_STD = dat_box.z;
    t_box   = (dat_box.x - mean(dat_box.x)).*1e9;
    
::

    %% Array of the initial separations 
    
    sep = [356 318 275 237 200 162 119 81]; % ps (To be adjusted to this specific file (numbers here correspond to exp of 16/10!)
    
    %% Plot all runs
    
    % figure,
    % for aa = 1:32
    %     subplot(5, 8, aa),
    %     imagesc(t, z, squeeze(matrice_STD(aa, :, :))*1e3), axis('xy'), colormap(turbo)
    %     caxis([0 200])
    %     xlim(.5*[-1 1])
    %     ylim([0 4000])
    %     if mod(aa, 8) == 1; ylabel('z (km)'); end
    % %     if aa <= 8; title(['Initial separation = ' num2str(sep(aa)) ' ps']); end
    % end
    % 
    % for bb = 1:8
    %     subplot(5, 8, aa + bb),
    %     area(t, squeeze(matrice_STD(bb, 1, :))*1e3), axis('xy'), colormap(turbo)
    %     xlim(.5*[-1 1])
    %     ylim([0 100])
    %     xlabel('t (ns)')
    %     if bb == 1; ylabel('Power (mW)'); end
    % end
    % 
    % set(gcf, 'Color', 'w')
    
    % % print(gcf, 'Generated figures/All_runs', '-dpng')




Fits of the initial conditions of the 5 runs
""""""""""""""""""""""""""""""""""""""""""""
    

::
    
    %% Fits of the initial conditions of the 5 runs
    for dd = 1:size(matrice_STD_sol, 1)
        data = squeeze(matrice_STD_sol(dd, :, :));
        ind_t = find((t > -.5).*(t < .5));
        t_cropped = t(ind_t);
        init_param   = [-0.372119815668203	0.0846938775510204
                        -0.369815668202765	0.0829446064139942
                        -0.365207373271889	0.0797376093294461
                        -0.369815668202765	0.0838192419825073
                        -0.372119815668203	0.0785714285714286
                       ];
    
    %     f1 = figure();
    %     plot(t, data(1, :)), hold on
    %     xlim([.5*[-1 1]])
    %     ylim([0 0.1])
    %     title(num2str(dd))
    %     init_param(dd, :) = ginput(1);
    %     close(f1)
    
        fit_array{dd, 1} = cell(numel(z), 1);
        sechEqn     = 'b*sech((x - a)/c).^2';
        startPoints = [init_param(dd, :) 15e-3];
        fit_array{dd, 1} = fit(t_cropped', data(1, ind_t)', sechEqn, 'Start', startPoints);
    end
    


Fit of the soliton at each sequence and extraction of the relevant fit parameters
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""
    
::

    %% Fit of the soliton at each sequence
    for dd = 5:size(matrice_STD_sol, 1)
        data = squeeze(matrice_STD_sol(dd, :, :));
        for ff = 2:(numel(z))
            startPoints = [fit_array{dd, ff - 1}.a fit_array{dd, ff - 1}.b fit_array{dd, ff - 1}.c];
            if (dd == 5 && ff >=153 && ff <=204) || (dd == 5 && ff == 31) % Manually input the starting parameters for the fit for the 5th soliton between sequences 154 and 204 
                figure(1),
                plot(t, data(ff, :)), hold on
                ylim([-0.001 0.140])
                xlim(.5*[-1 1])
                init_param2 = ginput(1);
                startPoints = [init_param2 15e-3];
                fit_array{dd, ff} = fit(t_cropped', data(ff, ind_t)', sechEqn, 'Start', startPoints, 'Lower', [-0.5 0.01 1e-3]);
                plot(t_cropped, fit_array{dd, ff}(t_cropped), 'linewidth', 2), hold off
            else
                fit_array{dd, ff} = fit(t_cropped', data(ff, ind_t)', sechEqn, 'Start', startPoints, 'Lower', [-0.5 0.01 1e-3]);
            end
            clc, disp([dd,ff])
        end
    end
    
    %% Extract relevant fit parameters
    for dd = 1:size(matrice_STD_sol, 1)
        for ff = 1:numel(z)
            trajectory(dd, ff) = fit_array{dd, ff}.a;
            pow(dd, ff)        = fit_array{dd, ff}.b;
            width(dd, ff)      = fit_array{dd, ff}.c;
        end
    end


Truncating to an integer number of oscillations
"""""""""""""""""""""""""""""""""""""""""""""""  

::

    %% Truncating to an integer number of oscillations
    for dd = 1:size(matrice_STD_sol, 1)
        figure(1),
        plot(z, trajectory(dd, :)), hold on
        yline(trajectory(dd, 1)), hold off
        ylim(0.5*[-1 1])
        xlim([8600 9600])
        xmax(dd, :) = ginput(1);
    end
    xmax = 5*round(xmax(:, 1)/5);


Fourier filtering of the trajectory
"""""""""""""""""""""""""""""""""""

::

    %% Fourier filtering of the trajectory
    cutoff = 0.01;
    Fs     = 1/diff(z(1:2));
    npts   = ((xmax - z(1))/5 + 1);
    sig    = zeros(size(matrice_STD_sol, 1), max(npts));
    
    for dd = 1:size(matrice_STD_sol, 1)
    
        dF          = Fs/npts(dd);
        freq        = ifftshift(-Fs/2:dF:Fs/2-dF);
        sig_FT      = fft(trajectory(dd, 1:npts(dd)) - mean(trajectory(dd, 1:npts(dd))));
        sig_FT_filt = sig_FT.*exp(-(freq./0.01).^1e10);
        sig_buff    = ifft(sig_FT_filt);
        sig(dd, 1:npts(dd))  = real(sig_buff) - mean(real(sig_buff));
    %     % Create figure of a Fourier transform of a trajectory + filtering
    %         figure, semilogy(fftshift(freq), abs(fftshift(sig_FT)), 'displayName', 'Raw spectrum', 'LineWidth', 2), grid on, hold on
    %         semilogy(fftshift(freq), abs(fftshift(sig_FT_filt)), 'displayName', 'Filtered spectrum', 'LineWidth', 2)
    %         title('Fourier transform of a soliton trajectory (trapping potential)')
    %         legend()
    %         xlim([-.03 .03])
    %         ylim([1e-2 1e3])
    %         xlabel('Wavenumber (km^{-1})')
    %         ylabel('Amplitude (dB)')
    %         set(gca, 'Fontsize', 24)
    %         set(gcf, 'Color', 'w')
        %     pause()
        %     print(gcf, 'Generated figures/Trajectory_Fourier_tranform_filtering', '-dpng')
    end
    clear cutoff sig_FT sig_FT_filt sig_buff

Generating the ST diagrams with jitter correction inferred from the filtered trajectories
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

::

    %% Generating the ST diagrams with jitter correction inferred from the filtered trajectories
    matrice_STD_corr = matrice_STD;
    correction = (trajectory(5, 1:max(npts)) - mean(trajectory(5, 1:max(npts)))) - sig(5, :);
    corr_pts          = round(correction*1e-9*160e9);
    for dd = 1:size(matrice_STD, 1)
        for cc = 1:max(npts)
            matrice_STD_corr(dd, cc, :) = circshift(squeeze(matrice_STD(dd, cc, :)), - corr_pts(cc));
        end
        clc, dd
    end
    correction = [correction zeros(1, (length(z) - max(npts)))];
    
    matrice_STD_sol_corr = matrice_STD_sol;
    % correction = (trajectory(5, 1:max(npts)) - mean(trajectory(5, 1:max(npts)))) - sig(5, :);
    corr_pts          = round(correction*1e-9*160e9);
    for dd = 1:size(matrice_STD_sol, 1)
        for cc = 1:max(npts)
            matrice_STD_sol_corr(dd, cc, :) = circshift(squeeze(matrice_STD_sol(dd, cc, :)), - corr_pts(cc));
        end
        clc, dd
    end
    
    % Case of the large square pulse
    % box_STD      = data_box.z; t_box = (data_box.x - mean(data_box.x)).*1e9; z_box = data_box.y;
    % box_STD_corr = zeros(size(box_STD));
    % for cc = 1:max(npts)
    %     box_STD_corr(cc, :) = circshift(squeeze(data_box.z(cc, :)), -corr_pts(cc));
    % end


Visualising in loop the ST diagrams of each run with and without correction
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

::

    %% Visualising in loop the ST diagrams of each run with and without correction
    % xl = [0 5000];
    % for aa = 1:1
    %     figure(1),
    %     if aa == 1
    %         pause
    %     end
    %     for dd = 1:size(matrice_STD, 1)
    %         subplot(2, 1, 1)
    %         pcolor(z, t, 1e3*squeeze(matrice_STD(dd, :, :))'), colormap(turbo), shading interp
    %         hcb = colorbar;
    %         set(get(hcb, 'Title'), 'String', 'Power (mW)');
    %         xlim(xl)
    %         ylim(.5*[-1 1])
    %         caxis([0 100])
    %         ylabel('Delay (ns)')
    %         %         title(num2str(aa))
    %         set(gca, 'FontSize', 20)
    %         subplot(2, 1, 2)
    %         pcolor(z, t, 1e3*squeeze(matrice_STD_corr(dd, :, :))'), colormap(turbo), shading interp
    %         hcb = colorbar;
    %         set(get(hcb, 'Title'), 'String', 'Power (mW)');
    %         xlim(xl)
    %         ylim(.5*[-1 1])
    %         caxis([0 100])
    %         xlabel('Distance (km)')
    %         ylabel('Delay (ns)')
    %         set(gca, 'FontSize', 20)
    %         set(gcf, 'Color', 'w')
    %         pause(.1)
    %         %         if dd == 1 % Creates an animated GIF
    %         %             pause()
    %         %             gif(['Generated figures/Anim_soliton_oscillation_correction.gif'], 'DelayTime', .2)
    %         %         else
    %         %             gif
    %         %         end
    %     end
    % end


Visualising in loop the ST diagrams of each run (single soliton) with and without correction
""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""


::

    %% Visualising in loop the ST diagrams of each run (single soliton) with and without correction
    % xl = [0 5000];
    % for aa = 1:10
    %     figure(1),
    %     if aa == 1
    %         pause
    %     end
    %     for dd = 1:size(matrice_STD_sol, 1)
    %         subplot(2, 1, 1)
    %         pcolor(z, t, 1e3*squeeze(matrice_STD_sol(dd, :, :))'), colormap(turbo), shading interp
    %         hcb = colorbar;
    %         set(get(hcb, 'Title'), 'String', 'Power (mW)');
    %         xlim(xl)
    %         ylim(.5*[-1 1])
    %         caxis([0 100])
    %         ylabel('Delay (ns)')
    %         %         title(num2str(aa))
    %         set(gca, 'FontSize', 20)
    %         subplot(2, 1, 2)
    %         pcolor(z, t, 1e3*squeeze(matrice_STD_sol_corr(dd, :, :))'), colormap(turbo), shading interp
    %         hcb = colorbar;
    %         set(get(hcb, 'Title'), 'String', 'Power (mW)');
    %         xlim(xl)
    %         ylim(.5*[-1 1])
    %         caxis([0 100])
    %         xlabel('Distance (km)')
    %         ylabel('Delay (ns)')
    %         set(gca, 'FontSize', 20)
    %         set(gcf, 'Color', 'w')
    %         pause(.1)
    %         %         if dd == 1 % Creates an animated GIF
    %         %             pause()
    %         %             gif(['Generated figures/Anim_soliton_oscillation_correction.gif'], 'DelayTime', .2)
    %         %         else
    %         %             gif
    %         %         end
    %     end
    % end
    
    
Plot all runs
"""""""""""""
    
::

    %% Plot all runs
    % figure,
    % for aa = 1:32
    %     subplot(5, 8, aa),
    %     imagesc(t, z, squeeze(matrice_STD_corr(aa, :, :))*1e3), axis('xy'), colormap(turbo)
    %     caxis([0 200])
    %     xlim(.5*[-1 1])
    %     ylim([0 8500])
    %     if mod(aa, 8) == 1; ylabel('z (km)'); end
    %     if aa <= 8; title(['Initial separation = ' num2str(sep(aa)) ' ps']); end
    % end
    % 
    % for bb = 1:8
    %     subplot(5, 8, aa + bb),
    %     area(t, squeeze(matrice_STD_corr(bb, 1, :))*1e3), axis('xy'), colormap(turbo)
    %     xlim(.5*[-1 1])
    %     ylim([0 100])
    %     xlabel('t (ns)')
    %     if bb == 1; ylabel('Power (mW)'); end
    % end
    % 
    % set(gcf, 'Color', 'w')
    
    % % print(gcf, 'Generated figures/All_runs', '-dpng')
    


Calculation and superposition of the particle trajectories for one specific realisation
"""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""""

::

    %% Calculation and superposition of the particle trajectories for one specific realisation
    
    lst_style  = {'-', ':'};
    lst_width  = [2 3];
    
    figure(1),
    imagesc(t, z, (squeeze(matrice_STD_corr(30, :, :)))), axis('xy'), colormap(gray) % Plots runs with normalised amplitudes
    xlim(.4*[-1 1])
    ylim([0 9000])
    xlabel('t (ns)')
    ylabel('z (km)'); end
    set(gca, 'Fontsize', 20)
    set(gcf, 'Color', 'w')
    
    t_init = -170e-12;
    P0     = 43e-3;
    period = 760;
    Hamiltonian_2_solitons_trapped
    
    figure(1),
    for qq = 1:length(eta)
        % subplot(1, 3, 1)
        hold on
        plot(q_values_phys(:, (2*qq - 1)), t_values_phys, 'linewidth', lst_width(qq), 'linestyle', lst_style(qq), 'color', col_vect(qq))
        title(['Initial separation = ' num2str(abs(t_init - pos_init(2))*1e12) ' ps'])
    end
    hold off
    
    % set(gcf, 'units', 'normalized', 'position', [1.0445    0.0563    0.2758    0.8569])
    % set(gcf, 'OuterPosition', [2563         200        1777        1066])
    % print(gcf, 'Generated figures/Soliton_pair', '-dpng')
