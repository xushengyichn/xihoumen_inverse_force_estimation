%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: ShengyiXu xushengyichn@outlook.com
%Date: 2023-10-09 22:23:15
%LastEditors: ShengyiXu xushengyichn@outlook.com
%LastEditTime: 2023-10-17 11:05:15
%FilePath: \Exercises-for-Techniques-for-estimation-in-dynamics-systemsf:\git\xihoumen_inverse_force_estimation\20231005 first version\KalmanMain.m
%Description: TODO:加上更多模态，不要只留下单一模态，看看能不能起到滤波的作用
%
%Copyright (c) 2023 by ${git_name_email}, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [result_Main] = KalmanMain(input,varargin)

    if nargin == 0
        clc; clear; close all;
        addpath(genpath("F:\git\ssm_tools\"))
        addpath(genpath("F:\git\Function_shengyi_package\"))
        addpath(genpath("F:\git\xihoumen_inverse_force_estimation\FEM_model\"))
        addpath(genpath("C:\Users\xushengyi\Documents\Github\"))
        

        subStreamNumberDefault = 2132;

        params = Init_fun();
        % input.ON = params.ON;
        % input.OFF = params.OFF;
        %% 0 绘图参数
        
        % input.num_figs_in_row = 12; %每一行显示几个图
        % input.figPos = params.figPosSmall; %图的大小，参数基于InitScript.m中的设置
        %设置图片间隔
        % input.gap_between_images = [0, 0];
        % input.figureIdx = 0;
        % input.fig_bool = params.ON;
        % input.displayText = params.ON;

        n = 4;
        [result] = viv2013(n, params.OFF);
        input.start_time = result.startDate;
        input.end_time = result.endDate;

        % input.start_time = datetime('2013-02-06 01:30:00', 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
        % input.end_time = datetime('2013-02-06 01:45:59', 'InputFormat', 'yyyy-MM-dd HH:mm:ss');

        
        input.acc_dir = "F:\test\result";
       
        % input.acc_dir = "Z:\Drive\Backup\SHENGYI_HP\F\test\result";
        
        input.lambda = 1e-1;
        input.sigma_p = 10000;
        input.omega_0_variation =1;
        input.Q_value =10 ^ (-8);
        input.R_value = 10 ^ (-6);
        
        input.modesel= [2,3,5,6,7,9,15,21,23,29,33,39,44,45];
        % KalmanMain(input,'showtext', false,'showplot',true,'shouldFilterYn', true,'shouldFilterp_filt_m', true);
        KalmanMain(input,'showtext', false,'showplot',false)
        return;
    end
    p = inputParser;
      addParameter(p, 'showtext', true, @islogical)
      addParameter(p, 'shouldFilterYn', false, @islogical)
      addParameter(p, 'shouldFilterp_filt_m', false, @islogical)
      addParameter(p, 'showplot', true, @islogical)
      addParameter(p, 'num_figs_in_row', 12, @isnumeric)
      addParameter(p, 'figPos', [100,100,400,300], @isnumeric)
      addParameter(p, 'gap_between_images', [0,0], @isnumeric)
      addParameter(p, 'figureIdx', 0, @isnumeric)
      addParameter(p, 'f_keep', 0.33*[0.9,1.1], @isnumeric)
      addParameter(p,'filterstyle','nofilter',@ischar);% fft use the function by myself, bandpass use the function in matlab. fft is much faster than bandpass but may not be accurate
      parse(p, varargin{:});
    filterstyle = p.Results.filterstyle;
    showtext = p.Results.showtext;
    shouldFilterYn = p.Results.shouldFilterYn;
    shouldFilterp_filt_m = p.Results.shouldFilterp_filt_m;
    showplot = p.Results.showplot;
    f_keep = p.Results.f_keep;
    num_figs_in_row = p.Results.num_figs_in_row;
    figPos = p.Results.figPos;
    gap_between_images = p.Results.gap_between_images;
    figureIdx = p.Results.figureIdx;

    %% 1 读取数据
    start_time = input.start_time;
    end_time = input.end_time;
    
    acc_dir = input.acc_dir;


    modesel= input.modesel;

    lambda = input.lambda;
    sigma_p = input.sigma_p;
    omega_0_variation = input.omega_0_variation;
    Q_value = input.Q_value;
    R_value = input.R_value;

    fig_bool = showplot;
    ON = true;
    OFF = false;
    
    [Acc_Data] = read_acceleration_data(start_time, end_time, acc_dir);

    switch filterstyle
        case 'fft'
        timeDifferences = diff(Acc_Data.mergedData.Time); % Returns a duration array
        dt = seconds(timeDifferences(1)); % Converts the first duration value to seconds and assigns to dt
        fs = 1 / dt;
        Acc_Data.mergedData.AC2_1 = fft_filter(fs, Acc_Data.mergedData.AC2_1, f_keep)';
        Acc_Data.mergedData.AC2_2 = fft_filter(fs, Acc_Data.mergedData.AC2_2, f_keep)'; 
        Acc_Data.mergedData.AC2_3 = fft_filter(fs, Acc_Data.mergedData.AC2_3, f_keep)';
        Acc_Data.mergedData.AC3_1 = fft_filter(fs, Acc_Data.mergedData.AC3_1, f_keep)';
        Acc_Data.mergedData.AC3_2 = fft_filter(fs, Acc_Data.mergedData.AC3_2, f_keep)';
        Acc_Data.mergedData.AC3_3 = fft_filter(fs, Acc_Data.mergedData.AC3_3, f_keep)';
        Acc_Data.mergedData.AC4_1 = fft_filter(fs, Acc_Data.mergedData.AC4_1, f_keep)';
        Acc_Data.mergedData.AC4_2 = fft_filter(fs, Acc_Data.mergedData.AC4_2, f_keep)';
        Acc_Data.mergedData.AC4_3 = fft_filter(fs, Acc_Data.mergedData.AC4_3, f_keep)';


        
        case 'bandpass'
        timeDifferences = diff(Acc_Data.mergedData.Time); % Returns a duration array
        dt = seconds(timeDifferences(1)); % Converts the first duration value to seconds and assigns to dt
        fs = 1 / dt;
        Acc_Data.mergedData.AC2_1 = bandpass(Acc_Data.mergedData.AC2_1', f_keep, fs)';
        Acc_Data.mergedData.AC2_2 = bandpass(Acc_Data.mergedData.AC2_2', f_keep, fs)';
        Acc_Data.mergedData.AC2_3 = bandpass(Acc_Data.mergedData.AC2_3', f_keep, fs)';
        Acc_Data.mergedData.AC3_1 = bandpass(Acc_Data.mergedData.AC3_1', f_keep, fs)';
        Acc_Data.mergedData.AC3_2 = bandpass(Acc_Data.mergedData.AC3_2', f_keep, fs)';
        Acc_Data.mergedData.AC3_3 = bandpass(Acc_Data.mergedData.AC3_3', f_keep, fs)';

        case 'nofilter'
        

        otherwise
        error('Invalid filter style. Choose either "fft" or "bandpass".')
    end
   
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

    % modesel= [2,3,5,6,7,9,15,21,23,29,33,39,44,45];
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
    % zeta = ones(size(modesel)) * 0/100;
    if showtext
        disp("Damping ratio of the structure is set as "+num2str(zeta));
    end
    omega = diag(2 * pi * Freq);
    CC_eq = 2 .* MM_eq .* omega .* zeta;

    phi = mode_vec; %模态向量 每一列是一个模态
    C = CC_eq; K = KK_eq; M = MM_eq;

    Gamma = C; % 对应文献中表述的符号
    omega2 = K;

    %% 3 传感器布置
    % accelerometer location
    % loc_acc= [578+1650/4*3;578+1650/2;578+1650/4];
    loc_acc = [1403];
    % loc_acc = [990.5; 1403; 1815.5];
    loc_vel = [];
    loc_dis = [];

    timeDifferences = diff(Acc_Data.mergedData.Time); % Returns a duration array
    dt = seconds(timeDifferences(1)); % Converts the first duration value to seconds and assigns to dt
    fs = 1 / dt;
    acc_names = ["Main span 1/4", "Main span 1/2", "Main span 3/4"];
    % yn(1, :) = Acc_Data.mergedData.AC2_1 / 1000 * 9.8;
    % yn(2, :) = Acc_Data.mergedData.AC2_3 / 1000 * 9.8;
    % yn(3, :) = Acc_Data.mergedData.AC3_1 / 1000 * 9.8;
    % yn(4, :) = Acc_Data.mergedData.AC3_3 / 1000 * 9.8;
    % yn(5, :) = Acc_Data.mergedData.AC4_1 / 1000 * 9.8;
    % yn(6, :) = Acc_Data.mergedData.AC4_3 / 1000 * 9.8;

    yn(1, :) = Acc_Data.mergedData.AC3_1 / 1000 * 9.8;
    yn(2, :) = Acc_Data.mergedData.AC3_3 / 1000 * 9.8;

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
    Q = Q_value * eye(ns);
    R = R_value * eye(n_sensors);
    S = zeros(ns, n_sensors);
    x0 = zeros(ns, 1);

    Q_xd = Q;
    np_m = nmodes;
    B_c_m = [zeros(nmodes, np_m); ...
                 eye(np_m, np_m)];
    J_c_m = [S_a * phi];
    B_d_m = A_c \ (A_d - eye(size(A_d))) * B_c_m;
    J_d_m = J_c_m;

    lambdas_m = [lambda] * ones(1, np_m);
    sigma_ps_m = [sigma_p] * ones(1, np_m);
    omega_0 = 2 * pi * Freq*omega_0_variation;


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
    % [f, magnitude] = fft_transform(fs, x_k_k(1, :));
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
    % loc_acc_v = [990.5; 1403; 1815.5];
    loc_acc_v = [1403];
    % loc_acc_v = [578+1650/4*3;578+1650/2;578+1650/4];
    loc_vel_v = [1403];
    loc_dis_v = [1403];
    [S_a_v, S_v_v, S_d_v, n_sensors_v] = sensor_selection(loc_acc_v, loc_vel_v, loc_dis_v, node_loc, phi, nodeondeck, Mapping_data);

    G_c_v = [S_d_v * phi - S_a_v * phi * omega2, S_v_v * phi - S_a_v * phi * Gamma];
    J_c_v = [S_a_v * phi];

    h_hat = G_c_v * x_filt_original + J_c_v * p_filt_m;

    %% 7 重构涡振响应
    p_reconstruct = p_filt_m;
    % p_reconstruct([1 2 3 ],:)=0;
    [~, yn_reconstruct, ~] = CalResponse(A_d, B_d, G_d, J_d, p_reconstruct, 0, 0, N, x0, ns, n_sensors);

    % fft
    % [f_origin, magnitude_origin] = fft_transform(1 / dt, yn(3, :));
    % [f_re, magnitude_re] = fft_transform(1 / dt, yn_reconstruct(3, :));

    %% 8 marginal likelihood
    logL = result.logL;
    logSk = result.logSk;
    logek = result.logek;
    invSk = result.invSk;
   
    t = Acc_Data.mergedData.Time;
    
    %% output
    result_Main.logL = logL;
    result_Main.logSk = logSk;
    result_Main.logek = logek;
    result_Main.invSk = invSk;
    result_Main.t = t;
    result_Main.p_filt_m = p_filt_m;
    result_Main.x_k_k = x_k_k;
    result_Main.nmodes = nmodes;
    result_Main.Freq = Freq;
    result_Main.fs = fs;
    result_Main.mode_deck = mode_deck;
    result_Main.yn=yn;
    result_Main.h_hat = h_hat;
    result_Main.yn_reconstruct = yn_reconstruct;
    result_Main.x_filt_original=x_filt_original;
    result_Main.MM = MM_eq;
    result_Main.CC = CC_eq;
    result_Main.KK = KK_eq;


    if fig_bool == ON

        for k1 = 1:length(loc_acc)
            [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);

            % subplot(1, nmodes, k1)
            plot(t, yn(2 * k1 - 1, :), 'Color', 'r')
            hold on
            plot(t, yn(2 * k1, :), 'Color', 'b')

            xlabel('time (s)')
            ylabel('acc (m/s^2)')
            set(gca, 'FontSize', 12)
            legend('left','right', 'Location', 'northwest')
            title([num2str(loc_acc(k1))]);
            ylim([-maxvalue, maxvalue])
        end

        for k1 = 1:length(loc_acc)
            [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);

            % subplot(1, nmodes, k1)
            plot(t, h_hat(2 * k1 - 1, :), 'Color', 'r')
            hold on
            plot(t, h_hat(2 * k1, :), 'Color', 'b')

            xlabel('time (s)')
            ylabel('acc (m/s^2)')
            set(gca, 'FontSize', 12)
            legend('left','right', 'Location', 'northwest')
            title([num2str(loc_acc(k1))] + "filter");
            ylim([-maxvalue, maxvalue])
        end

        for k1 = 1:length(loc_acc)
            [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);

            % subplot(1, nmodes, k1)
            plot(t, yn_reconstruct(2 * k1 - 1, :), 'Color', 'r')
            hold on
            plot(t, yn_reconstruct(2 * k1, :), 'Color', 'b')

            xlabel('time (s)')
            ylabel('acc (m/s^2)')
            set(gca, 'FontSize', 12)
            legend('left','right', 'Location', 'northwest')
            title([num2str(loc_acc(k1))] + "recalculate");
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

        % set(hFigure, 'name', 'filtered modal force frequency', 'Numbertitle', 'off');
        % % figureIdx=0;
        % [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
        % plot(f_re, magnitude_re)
        % hold on
        % plot(f_origin, magnitude_origin)
        % legend("cal", "measure")
        % xlim([0, 0.5])
        % 
        % [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
        % plot(f_re, log(magnitude_re))
        % hold on
        % plot(f_origin, log(magnitude_origin))
        % legend("cal", "measure")
        % xlim([0, 0.5])

        % [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
        % plot(t2, ifq1)
        % title("inst frequency")

        % [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
        % plot(f, magnitude)
        % title("frequency of estimated vibration")

        % instfreq(p, fd, td);


        
    end


    

end



