%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: ShengyiXu xushengyichn@outlook.com
%Date: 2023-10-09 22:23:15
%LastEditors: ShengyiXu xushengyichn@outlook.com
%LastEditTime: 2023-10-12 22:10:13
%FilePath: \Exercises-for-Techniques-for-estimation-in-dynamics-systemsf:\git\xihoumen_inverse_force_estimation\20231005 first version\Main.m
%Description: TODO:加上更多模态，不要只留下单一模态，看看能不能起到滤波的作用
%
%Copyright (c) 2023 by ${git_name_email}, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [result_Main] = Main(input,varargin)

    if nargin == 0
        clc; clear; close all;
        addpath(genpath("F:\git\ssm_tools\"))
        addpath(genpath("F:\git\Function_shengyi_package\"))
        addpath(genpath("F:\git\xihoumen_inverse_force_estimation\FEM_model\"))

        subStreamNumberDefault = 2132;

        params = Init_fun();
        input.ON = params.ON;
        input.OFF = params.OFF;
        %% 0 绘图参数
        input.fig_bool = params.ON;
        input.num_figs_in_row = 6; %每一行显示几个图
        input.figPos = params.figPosSmall; %图的大小，参数基于InitScript.m中的设置
        %设置图片间隔
        input.gap_between_images = [0, 0];
        input.figureIdx = 0;
        
        input.displayText = params.ON;

        n = 4;
        [result] = viv2013(n, params.OFF);
        input.start_time = result.startDate;
        input.end_time = result.endDate;

        % input.start_time = datetime('2013-02-06 01:30:00', 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
        % input.end_time = datetime('2013-02-06 01:45:59', 'InputFormat', 'yyyy-MM-dd HH:mm:ss');

        input.wind_dir = "F:\test\result_wind_10min";
        input.acc_dir = "F:\test\result";
        input.ncycle = 1;%计算气动阻尼时n个周期算一次阻尼比
        input.lambada = 1e-1;
        input.sigma_p = 100000;
        input.omega_0_variation =1;
        Main(input)
        return;
    end
    p = inputParser;
      addParameter(p, 'showtext', true, @islogical)
      addParameter(p, 'shouldFilterYn', false, @islogical)
      addParameter(p, 'shouldFilterp_filt_m', false, @islogical)
      
      parse(p, varargin{:});
    showtext = p.Results.showtext;
    shouldFilterYn = p.Results.shouldFilterYn;
    shouldFilterp_filt_m = p.Results.shouldFilterp_filt_m;
    %% 1 读取数据
    start_time = input.start_time;
    end_time = input.end_time;
    wind_dir = input.wind_dir;
    acc_dir = input.acc_dir;
    ncycle = input.ncycle;



    fig_bool = input.fig_bool;
    num_figs_in_row = input.num_figs_in_row;
    figPos = input.figPos;
    gap_between_images = input.gap_between_images;
    figureIdx = input.figureIdx;
    ON = input.ON;
    OFF = input.OFF;
   
    [Wind_Data] = read_wind_data(start_time, end_time, wind_dir);
    [Acc_Data] = read_acceleration_data(start_time, end_time, acc_dir);
   
    [uniqueTimestamps, ia, ~] = unique(Acc_Data.mergedData.Time, 'stable');

    % Find duplicates by checking the difference in lengths
    if showtext
        if length(uniqueTimestamps) < height(Acc_Data.mergedData)
            disp('There are duplicate timestamps:');
    
            % Get duplicated indices
            duplicatedIndices = setdiff(1:height(Acc_Data.mergedData), ia);
    
            % Display duplicates
            for i = 1:length(duplicatedIndices)
                disp(Acc_Data.mergedData.Time(duplicatedIndices(i)));
            end
    
        else
            disp('No duplicates found.');
        end
    end


    % select which wind profile to use

    %% 判断涡振模态
    % acc_1 = Acc_Data.mergedData.AC3_1;
    % t = Acc_Data.mergedData.Time;
    % % [p, fd, td] = pspectrum(acc_1, t, 'spectrogram', 'FrequencyResolution', 0.005);
    % % instfreq(p, fd, td);
    % 
    % % fs = 50;
    % % [f, magnitude] = fft_transform(fs, acc_1);
    % % figure
    % % plot(f, magnitude)


    %% 2 有限元模型
    % 读入ANSYS梁桥模型质量刚度矩阵  MCK矩阵 Import MCK matrix from ANSYS
    % 将ANSYS中的稀疏矩阵处理为完全矩阵 Handling sparse matrices in ANSYS as full matrices

    modesel= [2,3,5,6,7,9,15,21,23,29,33,39,44,45];
    % modesel = [23];
    nmodes = length(modesel); ns = nmodes * 2;
    Result = ImportMK(nmodes, 'KMatrix.matrix', 'MMatrix.matrix', 'nodeondeck.txt', 'KMatrix.mapping', 'nodegap.txt', 'modesel', modesel,'showtext',showtext);
    mode_deck = Result.mode_deck; mode_deck_re = Result.mode_deck_re; node_loc = Result.node_loc;
    Freq = Result.Freq;
    MM_eq = Result.MM_eq; KK_eq = Result.KK_eq;

    mode_vec = Result.mode_vec;
    nodeondeck = Result.nodeondeck;
    Mapping_data = Result.Mapping;
    zeta = ones(size(modesel)) * 0.3/100;
    omega = diag(2 * pi * Freq);
    CC_eq = 2 .* MM_eq .* omega .* zeta;

    phi = mode_vec; %模态向量 每一列是一个模态
    C = CC_eq; K = KK_eq; M = MM_eq;

    Gamma = C; % 对应文献中表述的符号
    omega2 = K;

    %% 3 传感器布置
    % accelerometer location
    % loc_acc= [578+1650/4*3;578+1650/2;578+1650/4];
    loc_acc = [990.5; 1403; 1815.5];
    loc_vel = [];
    loc_dis = [];

    timeDifferences = diff(Acc_Data.mergedData.Time); % Returns a duration array
    dt = seconds(timeDifferences(1)); % Converts the first duration value to seconds and assigns to dt
    fs = 1 / dt;
    acc_names = ["Main span 1/4", "Main span 1/2", "Main span 3/4"];
    yn(1, :) = Acc_Data.mergedData.AC2_1 / 1000 * 9.8;
    yn(2, :) = Acc_Data.mergedData.AC2_3 / 1000 * 9.8;
    yn(3, :) = Acc_Data.mergedData.AC3_1 / 1000 * 9.8;
    yn(4, :) = Acc_Data.mergedData.AC3_3 / 1000 * 9.8;
    yn(5, :) = Acc_Data.mergedData.AC4_1 / 1000 * 9.8;
    yn(6, :) = Acc_Data.mergedData.AC4_3 / 1000 * 9.8;

    if shouldFilterYn == true
        [f, magnitude] = fft_transform(fs, yn(3,:));
        [~, idx] = max(magnitude);
        f_keep = [f(idx) * 0.9, f(idx) * 1.1];
        yn = bandpass(yn', f_keep, fs)';
    end


    [S_a, S_v, S_d, n_sensors] = sensor_selection(loc_acc, loc_vel, loc_dis, node_loc, phi, nodeondeck, Mapping_data);
    % establish continuous time matrices
    [A_c, B_c, G_c, J_c] = ssmod_c_mode(nmodes, omega2, Gamma, phi, S_a, S_v, S_d);

    
    maxvalue = max(max(abs(yn)));
    
    % establish discrete time matrices
    
    [A_d, B_d, G_d, J_d, ~] = ssmod_c2d(A_c, B_c, G_c, J_c, dt);

    %% 4 反算模态力
    Q = 10 ^ (-8) * eye(ns);
    R = 10 ^ (-6) * eye(n_sensors);
    S = zeros(ns, n_sensors);
    x0 = zeros(ns, 1);

    Q_xd = Q;
    np_m = nmodes;
    B_c_m = [zeros(nmodes, np_m); ...
                 eye(np_m, np_m)];
    J_c_m = [S_a * phi];
    B_d_m = A_c \ (A_d - eye(size(A_d))) * B_c_m;
    J_d_m = J_c_m;

    lambdas_m = [input.lambada] * ones(1, np_m);
    sigma_ps_m = [input.sigma_p] * ones(1, np_m);
    omega_0 = 2 * pi * Freq*input.omega_0_variation;

    [F_c_m, L_c_m, H_c_m, sigma_w_m12] = ssmod_quasiperiod_coninue(lambdas_m, sigma_ps_m, omega_0, np_m);

    [~, ~, ~, ~, Fad_m, ~, Had_m, ~, Qad_m] = ssmod_lfm_aug(A_c, B_c_m, G_c, J_c_m, F_c_m, H_c_m, L_c_m, Q_xd, sigma_w_m12, dt);
    A_a_m = Fad_m;
    G_a_m = Had_m;
    Q_a_m = Qad_m;
    R_a_m = R;
    yn_a = yn;
    N = length(Acc_Data.mergedData.Time);
    NN = N;
    xa_history = zeros(ns + np_m * (2), NN);
    pa_history = zeros(ns + np_m * (2), NN);

    x_ak = zeros(ns + np_m * (2), 1);
    P_ak = 10 ^ (1) * eye(ns + np_m * (2));


    [x_k_k, x_k_kmin, P_k_k, P_k_kmin, result] = KalmanFilterNoInput(A_a_m, G_a_m, Q_a_m, R_a_m, yn_a, x_ak, P_ak, 'debugstate', true,'showtext',showtext);
    % [x_k_k, P_k_k] = RTSFixedInterval(A_a_m, x_k_k, x_k_kmin, P_k_k, P_k_kmin);
    xa_history = x_k_k;
    pa_history = P_k_k;
    x_filt_original = xa_history(1:ns, :);
    H_d_m = H_c_m;
    p_filt_m = H_d_m * xa_history(ns + 1:end, :);

    
    Pp_filt_m = H_d_m * pa_history(ns + 1:end, :);
    [f, magnitude] = fft_transform(fs, x_k_k(1, :));
    %% 5 fft and bandpass filter for the estimated modal force
    for k1 = 1:nmodes
            if shouldFilterp_filt_m == true
                [f, magnitude] = fft_transform(fs, x_k_k(k1,:));
                [~, idx] = max(magnitude);
                f_keep_temp = [f(idx) * 0.8, f(idx) * 1.2];
                % p_filt_m = fft_filter(fs, p_filt_m, f_keep_temp);
                p_filt_m = bandpass(p_filt_m', f_keep_temp, fs)';
            end
        [f_p_filt_m(k1, :), magnitude_filt_m(k1, :)] = fft_transform(fs, p_filt_m(k1, :));
    end

    %% 6 virtual sensoring
    loc_acc_v = [990.5; 1403; 1815.5];
    % loc_acc_v = [578+1650/4*3;578+1650/2;578+1650/4];
    loc_vel_v = [];
    loc_dis_v = [];
    [S_a_v, S_v_v, S_d_v, n_sensors_v] = sensor_selection(loc_acc_v, loc_vel_v, loc_dis_v, node_loc, phi, nodeondeck, Mapping_data);

    G_c_v = [S_d_v * phi - S_a_v * phi * omega2, S_v_v * phi - S_a_v * phi * Gamma];
    J_c_v = [S_a_v * phi];

    h_hat = G_c_v * x_filt_original + J_c_v * p_filt_m;

    %% 7 重构涡振响应
    p_reconstruct = p_filt_m;
    % p_reconstruct([1 2 3 ],:)=0;
    [~, yn_reconstruct, ~] = CalResponse(A_d, B_d, G_d, J_d, p_reconstruct, 0, 0, N, x0, ns, n_sensors);

    % fft
    [f_origin, magnitude_origin] = fft_transform(1 / dt, yn(3, :));
    [f_re, magnitude_re] = fft_transform(1 / dt, yn_reconstruct(3, :));

    %% 8 marginal likelihood
    logL = result.logL;
    logSk = result.logSk;
    logek = result.logek;

    %% 9 calculate aerodynamic damping ratio
    t = Acc_Data.mergedData.Time;
    f_keep = [Freq * 0.9, Freq * 1.1];

    for k1=1:nmodes
        [ifq_interpolated(k1, :)] = instfreq_samelength(fs, x_k_k(k1,:), f_keep(k1,:), t);
    end

    Fa = p_filt_m;
    dis = x_k_k(1:nmodes, :);
    vel = x_k_k(nmodes + 1:2*nmodes, :);
    
    Fa_filtered = bandpass(Fa', [min(Freq) * 0.9, max(Freq) * 1.1], fs)';

    % 找到峰值设定保存频率成分的变量
    top_freqs = cell(1, nmodes);
    for k1 = 1:nmodes
        [top_freqs{k1}, ~, ~] = extractSignificantFrequencies(fs, Fa_filtered(k1, :));
    end

    

    filtered_Fa = cell(1,nmodes);
    bandwidth = 0.01; %根据需要调整带宽
    for k1 = 1:nmodes
        frequencies = top_freqs{k1};
        Fa_current = Fa_filtered(k1, :);
        filtered_Fa{k1} =filterSignals(Fa_current, frequencies, fs);
    end
    

    peaks_locs_cell = cell(nmodes, 1); % 初始化一个cell数组以保存peaks和locs的结构

    
    % 寻找不同模态不同频率力信号的周期
    for k1 = 1:nmodes
        filtered_Fa_current = filtered_Fa{k1}; % 获取当前模式的滤波数据
        peaks_locs_struct = struct(); % 初始化一个结构体以保存peaks和locs       
        for k2 = 1:size(filtered_Fa_current, 1)
            f_temp = top_freqs{k1}(k2);
            T_temp = 1/f_temp;
            d = T_temp * 50;
            pp = 0;
            [peaks, locs] = findpeaks(filtered_Fa_current(k2, :), ...
                                      'MinPeakDistance', d, ...
                                      'MinPeakProminence', pp);
            % 将peaks和locs保存为结构体的字段
            peaks_locs_struct(k2).peaks = peaks;
            peaks_locs_struct(k2).locs = locs;
            % [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
            % plot(t,filtered_Fa_current(k2, :))
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
        each_mode_signals = filtered_Fa{i}; % 获取当前模式下的所有信号
        
        amp_mode = {}; % 初始化用于保存当前模式amp的cell
        zeta_all_mode = {}; % 初始化用于保存当前模式zeta_all的cell
        
        for j = 1:size(each_mode_signals, 1) % 遍历当前模式下的每个信号
            current_signal = each_mode_signals(j, :);
            
            k2 = 1;
            locs = peaks_locs_cell{i}(j).locs; % 获取当前信号的locs
            Fa_temp = filtered_Fa{i}(j, :); % 获取当前信号的Fa
            freq_temp = top_freqs{i}(j); % 获取当前信号的频率
            vel_temp = vel(i, :); % 获取当前信号的速度
            vel_temp_filt = bandpass(vel_temp, [freq_temp * 0.9, freq_temp * 1.1], fs); % 对速度信号进行滤波
            for k1 = 1:ncycle:(length(locs) - ncycle + 1)
                % Ensure we don't exceed array bounds
                if k1 + ncycle > length(locs)
                    break; % Exit the loop if we're going to exceed the array bounds
                end
            
                % This difference is a duration over ncycle periods
                dt_duration = t(locs(k1 + ncycle)) - t(locs(k1));
                dt = seconds(dt_duration); % Convert duration to seconds
            
                % Create an array of time intervals in seconds over ncycle periods
                timeIntervals = seconds(t(locs(k1):locs(k1 + ncycle)) - t(locs(k1)));
            
                % Work calculation over ncycle periods
                work(k2) = trapz(timeIntervals, Fa_temp(locs(k1):locs(k1 + ncycle)) .* vel_temp_filt(locs(k1):locs(k1 + ncycle)));
            
                % Average Amplitude calculation over ncycle periods
                sum_amplitudes = 0;
            
                for cycle_offset = 0:(ncycle - 1)
                    current_dis = dis(locs(k1 + cycle_offset):locs(k1 + cycle_offset + 1));
                    current_amp = (max(current_dis) - min(current_dis)) / 2;
                    sum_amplitudes = sum_amplitudes + current_amp;
                end
            
                amp(k2) = sum_amplitudes / ncycle;
            
                % Frequency calculation using the start and end of the ncycle periods
                f(k2) = (ifq1_interpolated(locs(k1)) + ifq1_interpolated(locs(k1 + ncycle))) / 2;
                omega(k2) = 2 * pi * f(k2);
                c(k2) = work(k2) / pi / amp(k2) ^ 2 / omega(k2) / ncycle;
                zeta_all(k2) = -c(k2) / 2 / omega(k2); %气动力做正功，移动到等号左边就为负阻尼
            
                % Time stamp calculation over ncycle periods
                timestamp_cycle(k2) = mean(t(locs(k1):locs(k1 + ncycle)));
                wind_color(k2) = interp1(Wind_Data.resultsTable_UA4.Time_Start, Wind_Data.resultsTable_UA4.U, timestamp_cycle(k2));
    
                
                % amp(k2) = sum_amplitudes / ncycle;
                % zeta_all(k2) = -c(k2) / 2 / omega(k2); %气动力做正功，移动到等号左边就为负阻尼
                
                
                k2 = k2 + 1; % Increment the continuous index
            end
            
            amp_mode{j} = amp; % 将此信号的amp数组保存到当前模式的cell中
            zeta_all_mode{j} = zeta_all; % 将此信号的zeta_all数组保存到当前模式的cell中
        end
        
        amp_cell{i} = amp_mode; % 将当前模式的amp cell保存到总cell中
        zeta_all_cell{i} = zeta_all_mode; % 将当前模式的zeta_all cell保存到总cell中
    end

    a=1;

    % k2 = 1; % Continuous index for output arrays
    
    % for k1 = 1:ncycle:(length(locs) - ncycle + 1)
    %     % Ensure we don't exceed array bounds
    %     if k1 + ncycle > length(locs)
    %         break; % Exit the loop if we're going to exceed the array bounds
    %     end
    
    %     % This difference is a duration over ncycle periods
    %     dt_duration = t(locs(k1 + ncycle)) - t(locs(k1));
    %     dt = seconds(dt_duration); % Convert duration to seconds
    
    %     % Create an array of time intervals in seconds over ncycle periods
    %     timeIntervals = seconds(t(locs(k1):locs(k1 + ncycle)) - t(locs(k1)));
    
    %     % Work calculation over ncycle periods
    %     work(k2) = trapz(timeIntervals, Fa(locs(k1):locs(k1 + ncycle)) .* vel(locs(k1):locs(k1 + ncycle)));
    
    %     % Average Amplitude calculation over ncycle periods
    %     sum_amplitudes = 0;
    
    %     for cycle_offset = 0:(ncycle - 1)
    %         current_dis = dis(locs(k1 + cycle_offset):locs(k1 + cycle_offset + 1));
    %         current_amp = (max(current_dis) - min(current_dis)) / 2;
    %         sum_amplitudes = sum_amplitudes + current_amp;
    %     end
    
    %     amp(k2) = sum_amplitudes / ncycle;
    
    %     % Frequency calculation using the start and end of the ncycle periods
    %     f(k2) = (ifq1_interpolated(locs(k1)) + ifq1_interpolated(locs(k1 + ncycle))) / 2;
    %     omega(k2) = 2 * pi * f(k2);
    %     c(k2) = work(k2) / pi / amp(k2) ^ 2 / omega(k2) / ncycle;
    %     zeta_all(k2) = -c(k2) / 2 / omega(k2); %气动力做正功，移动到等号左边就为负阻尼
    
    %     % Time stamp calculation over ncycle periods
    %     timestamp_cycle(k2) = mean(t(locs(k1):locs(k1 + ncycle)));
    %     wind_color(k2) = interp1(Wind_Data.resultsTable_UA4.Time_Start, Wind_Data.resultsTable_UA4.U, timestamp_cycle(k2));
    
    %     k2 = k2 + 1; % Increment the continuous index
    % end
    
    % amp_filt_kalman = amp;
    
    % zeta_aero_filt_kalman = zeta_all; %气动阻尼= 总阻尼-机械阻尼


    % zeta_filt_kalman = zeta_all;



    if fig_bool == ON

        for k1 = 1:length(acc_names)
            [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);

            % subplot(1, nmodes, k1)
            plot(t, yn(2 * k1 - 1, :), 'Color', 'r')
            hold on
            plot(t, yn(2 * k1, :), 'Color', 'b')

            xlabel('time (s)')
            ylabel('acc (m/s^2)')
            set(gca, 'FontSize', 12)
            legend('left','right', 'Location', 'northwest')
            title([acc_names(k1)]);
            ylim([-maxvalue, maxvalue])
        end

        for k1 = 1:length(acc_names)
            [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);

            % subplot(1, nmodes, k1)
            plot(t, h_hat(2 * k1 - 1, :), 'Color', 'r')
            hold on
            plot(t, h_hat(2 * k1, :), 'Color', 'b')

            xlabel('time (s)')
            ylabel('acc (m/s^2)')
            set(gca, 'FontSize', 12)
            legend('left','right', 'Location', 'northwest')
            title([acc_names(k1)] + "filter");
            ylim([-maxvalue, maxvalue])
        end

        for k1 = 1:length(acc_names)
            [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);

            % subplot(1, nmodes, k1)
            plot(t, yn_reconstruct(2 * k1 - 1, :), 'Color', 'r')
            hold on
            plot(t, yn_reconstruct(2 * k1, :), 'Color', 'b')

            xlabel('time (s)')
            ylabel('acc (m/s^2)')
            set(gca, 'FontSize', 12)
            legend('left','right', 'Location', 'northwest')
            title([acc_names(k1)] + "recalculate");
            ylim([-maxvalue, maxvalue])
        end

        for k1 = 1:nmodes
            [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);

            % subplot(1, nmodes, k1)
            plot(t, p_filt_m(k1, :), 'Color', 'r')
            hold on
            %         plot(t, p_m_real(k1, :), 'Color', 'b', 'LineStyle', '--')
            xlabel('time (s)')
            ylabel('Modal force ')
            set(gca, 'FontSize', 12)
            legend('filtered modal force', 'Location', 'northwest')
            title(['mode ', num2str(modesel(k1))]);
            % ylim([-1e3, 1e3])
        end

        set(hFigure, 'name', 'modal force estimation', 'Numbertitle', 'off');

        for k1 = 1:nmodes
            [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
            % subplot(1, nmodes, k1)
            plot(f_p_filt_m(k1, :), magnitude_filt_m(k1, :), 'Color', 'r')
            hold on
            %         plot(f_p_real_m(k1, :), magnitude_real_m(k1, :), 'Color', 'b', 'LineStyle', '--')
            xlabel('frequency (Hz)')
            ylabel('magnitude (N)')
            set(gca, 'FontSize', 12)
            legend('filtered modal force frequency', 'Location', 'northwest')
            title(['mode ', num2str(modesel(k1))]);
            xlim([0, 0.5])
            % ylim([0, 50])
        end

        set(hFigure, 'name', 'filtered modal force frequency', 'Numbertitle', 'off');
        % figureIdx=0;
        [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
        plot(f_re, magnitude_re)
        hold on
        plot(f_origin, magnitude_origin)
        legend("cal", "measure")
        xlim([0, 0.5])

        [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
        plot(f_re, log(magnitude_re))
        hold on
        plot(f_origin, log(magnitude_origin))
        legend("cal", "measure")
        xlim([0, 0.5])

        [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
        plot(t2, ifq1)
        title("inst frequency")

        [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
        plot(f, magnitude)
        title("frequency of estimated vibration")

        instfreq(p, fd, td);

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
    end

    result_Main.logL = logL;
    result_Main.logSk = logSk;
    result_Main.logek = logek;
    

end



