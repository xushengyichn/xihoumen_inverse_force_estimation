function [result]=Cal_aero_damping_ratio()
    
input.ncycle = 1;%计算气动阻尼时n个周期算一次阻尼比
ncycle = input.ncycle;
input.wind_dir = "F:\test\result_wind_10min";
 % input.wind_dir = "Z:\Drive\Backup\SHENGYI_HP\F\test\result_wind_10min";
[Wind_Data] = read_wind_data(start_time, end_time, wind_dir);
 %% 9 calculate aerodynamic damping ratio
    
    t = Acc_Data.mergedData.Time;
    f_keep = [Freq * 0.9, Freq * 1.1];



    Fa = p_filt_m;
    dis = x_k_k(1:nmodes, :);
    vel = x_k_k(nmodes + 1:2*nmodes, :);
    
    Fa_filtered = bandpass(Fa', [min(Freq) * 0.9, max(Freq) * 1.1], fs)';
    vel_filtered = bandpass(vel', [min(Freq) * 0.9, max(Freq) * 1.1], fs)';
    dis_filtered = bandpass(dis', [min(Freq) * 0.9, max(Freq) * 1.1], fs)';

    % 找到峰值设定保存频率成分的变量
    top_freqs = cell(1, nmodes);
    for k1 = 1:nmodes
        [top_freqs{k1}, ~, ~] = extractSignificantFrequencies(fs, Fa_filtered(k1, :));
        [top_freqs_vel{k1}, ~, ~] = extractSignificantFrequencies(fs, vel_filtered(k1, :));
    end
    
    ifq_interpolated_allmodes = cell(1, nmodes);
    for k1=1:nmodes
        for k2 = 1:length(top_freqs{k1})
            f_keep_temp = [top_freqs{k1}(k2) * 0.9, top_freqs{k1}(k2) * 1.1];
            [ifq_interpolated_allmodes{k1}{k2}] = instfreq_samelength(fs, Fa(k1,:), f_keep_temp, t);
        end
    end

    tic
    bandwidth = 0.01; %根据需要调整带宽
    filtered_Fa = cell(1,nmodes);
    for k1 = 1:nmodes
        frequencies = top_freqs{k1};
        Fa_current = Fa_filtered(k1, :);
        filtered_Fa{k1} =filterSignals(Fa_current, frequencies, fs, bandwidth);
    end

    filtered_vel = cell(1,nmodes);
    for k1 = 1:nmodes
        frequencies = top_freqs{k1};
        vel_current = vel_filtered(k1, :);
        filtered_vel{k1} =filterSignals(vel_current, frequencies, fs, bandwidth);
    end

    filtered_dis = cell(1,nmodes);
    for k1 = 1:nmodes
        frequencies = top_freqs{k1};
        dis_current = dis_filtered(k1, :);
        filtered_dis{k1} =filterSignals(dis_current, frequencies, fs, bandwidth);
    end
    toc

    peaks_locs_cell = cell(nmodes, 1); % 初始化一个cell数组以保存peaks和locs的结构

    
    % 寻找不同模态不同频率力信号的周期
    for k1 = 1:nmodes
        filtered_Fa_current = filtered_Fa{k1}; % 获取当前模式的滤波数据
        peaks_locs_struct = struct(); % 初始化一个结构体以保存peaks和locs       
        for k2 = 1:length(filtered_Fa_current)
            f_temp = top_freqs{k1}(k2);
            T_temp = 1/f_temp;
            d = T_temp * 50*0.9;
            pp = 0;
            [peaks, locs] = findpeaks(filtered_Fa_current{k2}, ...
                                      'MinPeakDistance', d, ...
                                      'MinPeakProminence', pp);
            % 将peaks和locs保存为结构体的字段
            peaks_locs_struct(k2).peaks = peaks;
            peaks_locs_struct(k2).locs = locs;
            % [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
            % plot(t,filtered_Fa_current{k2})
            % hold on 
            % scatter(t(locs),peaks)
            % hold off
            % str = sprintf('mode freq %.2f Hz, mode sel %d',Freq(k1),modesel(k1));
            % title(str)
        end 
        peaks_locs_cell{k1} = peaks_locs_struct; % 将结构体保存在cell数组中
    end


    amp_cell = cell(size(filtered_Fa)); % 初始化用于保存amp数组的cell
    zeta_all_cell = cell(size(filtered_Fa)); % 初始化用于保存zeta_all数组的cell

    for i = 1:length(filtered_Fa)
        force_mode_signals = filtered_Fa{i}; % 获取当前模式下的力信号
        vel_filtered_mode_signal = filtered_vel{i}; % 获取当前模式下的速度信号
        dis_filtered_mode_signal = filtered_dis{i}; % 获取当前模式下的位移信号
        ifq_interpolated_mode = ifq_interpolated_allmodes{i}; % 获取当前模式下的ifq_interpolated
        peaks_locs_cell_mode =peaks_locs_cell{i};

        amp_mode = {}; % 初始化用于保存当前模式amp的cell
        zeta_all_mode = {}; % 初始化用于保存当前模式zeta_all的cell
        
        for j = 1:length(force_mode_signals) % 遍历当前模式下的每个信号
            Fa_temp = force_mode_signals{j};
            vel_temp = vel_filtered_mode_signal{j};
            freq_temp = ifq_interpolated_mode{j};
            dis_temp = dis_filtered_mode_signal{j};
            locs = peaks_locs_cell_mode(j).locs;
            [result] = compute_dynamics_parameters(ncycle, t, Fa_temp, vel_temp, freq_temp, dis_temp , locs, Wind_Data);
            
            
            amp_mode{j} = result.amp; % 将此信号的amp数组保存到当前模式的cell中
            zeta_all_mode{j} = result.zeta_all; % 将此信号的zeta_all数组保存到当前模式的cell中

            [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
            scatter(result.amp*max(mode_deck(:,i)), result.zeta_all)
            titlestr= sprintf('mode freq %.2f Hz, mode sel %d, freq %.2f Hz',Freq(i),modesel(i),top_freqs{i}(j));
            title(titlestr)
            

        end
        
        amp_cell{i} = amp_mode; % 将当前模式的amp cell保存到总cell中
        zeta_all_cell{i} = zeta_all_mode; % 将当前模式的zeta_all cell保存到总cell中
    end












    if 0
        
            % Colorize scatter plot based on wind_U
        scatter(amp_filt_kalman, zeta_aero_filt_kalman, [], wind_color, 'filled');

        hold on;

        reference_amp = [0 600];
        reference_zeta = [0 0];
        plot(reference_amp, reference_zeta);

        title("amplitude dependent aerodynamic damping ratio");
        ylim([-0.5, 0.5]);

        % Choose a color map (for example, 'jet') and display the color scale
        colormap('jet');
        colorbar;

        [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);

        % Colorize scatter plot based on wind_U
        scatter(amp_filt_kalman * max(mode_deck), zeta_aero_filt_kalman, [], wind_color, 'filled');

        hold on;

        reference_amp = [0 0.12];
        reference_zeta = [0 0];
        plot(reference_amp, reference_zeta);

        title("amplitude dependent aerodynamic damping ratio");
        ylim([-0.05, 0.05]);

        % Choose a color map (for example, 'jet') and display the color scale
        colormap('jet');
        colorbar;

        [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);

        % Colorize scatter plot based on wind_U
        scatter(amp_filt_kalman * max(mode_deck), work, [], wind_color, 'filled');
        hold on;
        title("amplitude dependent work");
        % Choose a color map (for example, 'jet') and display the color scale
        colormap('jet');
        colorbar;

                % 画图来验证峰值检测的准确性
        [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
        plot(t, dis);
        hold on;
        plot(t(locs), peaks, 'ro'); % 红色的圆圈表示检测到的峰值
        hold off;
        title('Peak detection');
        xlabel('Time');
        ylabel('Displacement');

        [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
        [p, fd, td] = pspectrum(dis, t, 'spectrogram', 'FrequencyResolution', 0.005);
        instfreq(p, fd, td);

        [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);

    end
end