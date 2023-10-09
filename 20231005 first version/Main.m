function [result] = Main(input)

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
        input.num_figs_in_row = 12; %每一行显示几个图
        input.figPos = params.figPosSmall; %图的大小，参数基于InitScript.m中的设置
        %设置图片间隔
        input.gap_between_images = [0, 0];
        input.figureIdx = 0;

        %% 1 读取数据
        n = 4;
        displayDates = 1;
        [result] = viv2013(n, displayDates);
        input.start_time = result.startDate;
        input.end_time = result.endDate - hours(1.6);
        % start_time = datetime('2013-04-03 16:24:00', 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
        % end_time = datetime('2013-04-08 15:30:59', 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
        input.wind_dir = "F:\test\result_wind_10min";
        input.acc_dir = "F:\test\result";
        Main(input)
        return;
    end

    start_time = input.start_time;
    end_time = input.end_time;
    wind_dir = input.wind_dir;
    acc_dir = input.acc_dir;
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

    % select which wind profile to use

    %% 判断涡振模态
    acc_1 = Acc_Data.mergedData.AC3_1;
    t = Acc_Data.mergedData.Time;
    [p, fd, td] = pspectrum(acc_1, t, 'spectrogram', 'FrequencyResolution', 0.005);
    instfreq(p, fd, td);

    fs = 50;
    [f, magnitude] = fft_transform(fs, acc_1);
    figure
    plot(f, magnitude)
    %% 2 有限元模型
    % 读入ANSYS梁桥模型质量刚度矩阵  MCK矩阵 Import MCK matrix from ANSYS
    % 将ANSYS中的稀疏矩阵处理为完全矩阵 Handling sparse matrices in ANSYS as full matrices

    % modesel= [2,3,5,6,7,9,15,21,23,29,33,39,44,45];
    modesel = [23];
    nmodes = length(modesel); ns = nmodes * 2;
    Result = ImportMK(nmodes, 'KMatrix.matrix', 'MMatrix.matrix', 'nodeondeck.txt', 'KMatrix.mapping', 'nodegap.txt', 'modesel', modesel);
    mode_deck = Result.mode_deck; mode_deck_re = Result.mode_deck_re; node_loc = Result.node_loc;
    Freq = Result.Freq;
    % Freq = 0.328194;
    MM_eq = Result.MM_eq; KK_eq = Result.KK_eq;

    mode_vec = Result.mode_vec;
    nodeondeck = Result.nodeondeck;
    Mapping_data = Result.Mapping;
    % CC_eq = 0.1 * MM_eq + 0.005 * KK_eq; %人为指定瑞利阻尼
    zeta = ones(size(modesel)) * 0.3/100;
    omega = diag(2 * pi * Freq);
    CC_eq = 2 .* MM_eq .* omega .* zeta;

    phi = mode_vec; %模态向量 每一列是一个模态
    C = CC_eq; K = KK_eq; M = MM_eq;

    Gamma = C; % 对应文献中表述的符号
    omega2 = K;

    %% 3 传感器布置
    %

    % accelerometer location
    % loc_acc= [578+1650/4*3;578+1650/2;578+1650/4];
    loc_acc = [990.5; 1403; 1815.5];

    loc_vel = [];
    loc_dis = [];
    % yn(1, :) = Acc_Data.mergedData.AC3_1/1000*9.8;
    % yn(2, :) = Acc_Data.mergedData.AC3_3/1000*9.8;

    yn(1, :) = Acc_Data.mergedData.AC2_1 / 1000 * 9.8;
    yn(2, :) = Acc_Data.mergedData.AC2_3 / 1000 * 9.8;
    yn(3, :) = Acc_Data.mergedData.AC3_1 / 1000 * 9.8;
    yn(4, :) = Acc_Data.mergedData.AC3_3 / 1000 * 9.8;
    yn(5, :) = Acc_Data.mergedData.AC4_1 / 1000 * 9.8;
    yn(6, :) = Acc_Data.mergedData.AC4_3 / 1000 * 9.8;
    timeDifferences = diff(Acc_Data.mergedData.Time); % Returns a duration array
    dt = seconds(timeDifferences(1)); % Converts the first duration value to seconds and assigns to dt

    [S_a, S_v, S_d, n_sensors] = sensor_selection(loc_acc, loc_vel, loc_dis, node_loc, phi, nodeondeck, Mapping_data);

    % establish continuous time matrices
    [A_c, B_c, G_c, J_c] = ssmod_c_mode(nmodes, omega2, Gamma, phi, S_a, S_v, S_d);

    acc_names = ["Main span 1/4", "Main span 1/2", "Main span 3/4"];

    maxvalue = max(max(abs(yn)));
    fs = 1 / dt;

    % establish discrete time matrices
    N = length(Acc_Data.mergedData.Time);
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

    lambdas_m = [1e-1] * ones(1, np_m);

    sigma_ps_m = [100000] * ones(1, np_m);

    omega_0 = 2 * pi * Freq;

    [F_c_m, L_c_m, H_c_m, sigma_w_m12] = ssmod_quasiperiod_coninue(lambdas_m, sigma_ps_m, omega_0, np_m);

    [~, ~, ~, ~, Fad_m, ~, Had_m, ~, Qad_m] = ssmod_lfm_aug(A_c, B_c_m, G_c, J_c_m, F_c_m, H_c_m, L_c_m, Q_xd, sigma_w_m12, dt);
    A_a_m = Fad_m;
    G_a_m = Had_m;
    Q_a_m = Qad_m;
    R_a_m = R;
    yn_a = yn;
    NN = N;
    xa_history = zeros(ns + np_m * (2), NN);
    pa_history = zeros(ns + np_m * (2), NN);

    x_ak = zeros(ns + np_m * (2), 1);
    P_ak = 10 ^ (1) * eye(ns + np_m * (2));

    % G_a=G_a_m; A_a=A_a_m; Q_a=Q_a_m;
    [x_k_k, x_k_kmin, P_k_k, P_k_kmin, result] = KalmanFilterNoInput(A_a_m, G_a_m, Q_a_m, R_a_m, yn_a, x_ak, P_ak, 'debugstate', true);
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

        if 0
            % p_filt_m(k1, :) = fft_filter(fs, p_filt_m(k1, :), [0.15 0.55]);
            % p_filt_m(k1, :) = bandpass(p_filt_m(k1, :), [0.3 0.55], fs);
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
    % p_reconstruct([1 2 3 4 5 6 7 8 10 11 12 13 14],:)=0;
    % p_reconstruct([1],:)=0;
    [~, yn_reconstruct, ~] = CalResponse(A_d, B_d, G_d, J_d, p_reconstruct, 0, 0, N, x0, ns, n_sensors);

    seq1 = yn(:);
    seq2 = yn_reconstruct(:);
    % mse = immse(seq1, seq2);

    % fft
    [f_origin, magnitude_origin] = fft_transform(1 / dt, yn(3, :));
    [f_re, magnitude_re] = fft_transform(1 / dt, yn_reconstruct(3, :));

    %% 8 marginal likelihood

    logL = result.logL;
    logSk = result.logSk;
    logek = result.logek;

    %% 9 calculate aerodynamic damping ratio
    t = Acc_Data.mergedData.Time;
    % t = 0:1/fs:1/fs*(length(t)-1);
    f_keep = [Freq * 0.9, Freq * 1.1];
    [x_k_k_filtered] = fft_filter(fs, x_k_k(1, :), f_keep);
    [p, fd, td] = pspectrum(x_k_k_filtered, t, 'spectrogram', 'FrequencyResolution', 0.005);
    [ifq1, t1] = instfreq(p, fd, td);

    % 首先，我们对于t1范围内的t值执行插值
    inside_indices = t >= min(t1) & t <= max(t1);
    t_inside = t(inside_indices);
    ifq1_interpolated_inside = interp1(t1, ifq1, t_inside, 'linear');

    % 然后，我们对于t1范围外的t值保留最近的ifq1值
    t_outside_left = t(t < min(t1));
    t_outside_right = t(t > max(t1));
    ifq1_extrapolated_left = repmat(ifq1(1), size(t_outside_left));
    ifq1_extrapolated_right = repmat(ifq1(end), size(t_outside_right));

    % 最后，我们将这些值合并到一个数组中
    ifq1_interpolated = [ifq1_extrapolated_left; ifq1_interpolated_inside; ifq1_extrapolated_right];

    Fa = p_filt_m;
    vel = x_k_k(2, :);
    dis = x_k_k(1, :);
    t = t;

    % 找到峰值
    d = 100;
    pp = 0;
    [peaks, locs] = findpeaks(dis, 'MinPeakDistance', d, 'MinPeakProminence', pp);

    ncycle = 10;

    k2 = 1; % Continuous index for output arrays

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
        work(k2) = trapz(timeIntervals, Fa(locs(k1):locs(k1 + ncycle)) .* vel(locs(k1):locs(k1 + ncycle)));

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

        k2 = k2 + 1; % Increment the continuous index
    end

    amp_filt_kalman = amp;
    zeta_filt_kalman = zeta_all;
    zeta_aero_filt_kalman = zeta_filt_kalman - zeta; %气动阻尼= 总阻尼-机械阻尼

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
            legend('filtered modal force', 'Location', 'northwest')
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
            legend('filtered modal force', 'Location', 'northwest')
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
            legend('filtered modal force', 'Location', 'northwest')
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
        plot(t1, ifq1)
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
        ylim([-0.25, 0.25]);

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

    result = 1;
    %% necessary functions
    function node = FindNodewithLocation(loc, node_loc, nodeondeck)
        %myFun - Description
        %
        % Syntax: result = findnode(input)
        %
        % Long description
        node_seq = zeros(1, length(loc));

        for k1 = 1:length(loc)
            node_seq(k1) = find(node_loc >= loc(k1), 1);
        end

        node = nodeondeck(node_seq, :);
    end

    function [S_a, S_v, S_d, n_sensors] = sensor_selection(loc_acc, loc_vel, loc_dis, node_loc, phi, nodeondeck, Mapping_data)

        if ~isempty(loc_acc)
            acc_node = FindNodewithLocation(loc_acc, node_loc, nodeondeck);
            acc_node_list = reshape(permute(acc_node, [2 1]), [], 1); % 交错重塑
            % acc_node_list=acc_node(:);
        else
            acc_node_list = [];
        end

        if ~isempty(loc_vel)
            vel_node = FindNodewithLocation(loc_vel, node_loc, nodeondeck);
            vel_node_list = reshape(permute(vel_node, [2 1]), [], 1);
        else
            vel_node_list = [];
        end

        if ~isempty(loc_dis)
            dis_node = FindNodewithLocation(loc_dis, node_loc, nodeondeck);
            dis_node_list = reshape(permute(dis_node, [2 1]), [], 1);
        else
            dis_node_list = [];
        end

        n_acc = length(acc_node_list);
        n_vel = length(vel_node_list);
        n_dis = length(dis_node_list);
        acc_matrix_seq = node2matrixseq(acc_node_list, Mapping_data);
        vel_matrix_seq = node2matrixseq(vel_node_list, Mapping_data);
        dis_matrix_seq = node2matrixseq(dis_node_list, Mapping_data);

        n_sensors = n_acc + n_vel + n_dis;

        S_a = zeros(n_sensors, size(phi, 1));
        S_v = zeros(n_sensors, size(phi, 1));
        S_d = zeros(n_sensors, size(phi, 1));

        for k1 = 1:n_sensors

            if k1 <= n_acc
                S_a(k1, acc_matrix_seq(k1)) = 1;
            elseif k1 <= n_acc + n_vel
                S_v(k1, vel_matrix_seq(k1 - n_acc)) = 1;
            else
                S_d(k1, dis_matrix_seq(k1 - n_acc - n_vel)) = 1;
            end

        end

    end

    function matrix_seq = node2matrixseq(node_list, KMmapping)
        UYNode = KMmapping.Node(KMmapping.DOF == 'UY');
        UYMatrix = KMmapping.MatrixEqn(KMmapping.DOF == 'UY');
        indices = [];

        for k1 = 1:length(node_list)
            indices(k1) = find(UYNode == node_list(k1));
        end

        if isempty(indices)
            matrix_seq = [];
        else
            matrix_seq = UYMatrix(indices);
        end

    end

    function [A_c, B_c, G_c, J_c] = ssmod_c_mode(nmodes, omega2, Gamma, phi, S_a, S_v, S_d)
        A_c = [zeros(nmodes), eye(nmodes); ...
                   -omega2, -Gamma];
        B_c = [zeros(nmodes, nmodes); ...
                   eye(nmodes, nmodes)];
        G_c = [S_d * phi - S_a * phi * omega2, ...
                   S_v * phi - S_a * phi * Gamma];
        J_c = [S_a * phi];
    end

    function [F_c, L_c, H_c, sigma_w12] = ssmod_quasiperiod_coninue(lambdas, sigma_ps, omega_0, np)
        % create arrays for latent force model
        F_c_array = zeros(2, 2, np);
        L_c_array = zeros(2, 2, np);
        H_c_array = zeros(1, 2, np);
        sigma_w = zeros(np, 1);

        for k1 = 1:np
            [F_c_array(:, :, k1), L_c_array(:, :, k1), H_c_array(:, :, k1), sigma_w(k1)] = ssmod_quasiperiod(lambdas(k1), sigma_ps(k1), omega_0(k1));
        end

        F_c = [];
        L_c = [];
        H_c = [];

        for k1 = 1:np
            F_c = blkdiag(F_c, F_c_array(:, :, k1));
            L_c = blkdiag(L_c, L_c_array(:, :, k1));
            H_c = blkdiag(H_c, H_c_array(:, :, k1));
        end

        sigma_w12 = [];

        for k1 = 1:np
            sigma_w12 = blkdiag(sigma_w12, sigma_w(k1) * eye(2));
        end

    end

    function [xn, yn, xn_true] = CalResponse(A_d, B_d, G_d, J_d, p, Q, R, N, x0, ns, n_sensors)
        w = sqrt(Q) * randn(ns, N);
        v = sqrt(R) * randn(n_sensors, N);
        xn = x0;
        xn_true = x0;

        for t1 = 1:N - 1
            xn(:, t1 + 1) = A_d * xn(:, t1) + B_d * p(:, t1) + w(:, t1);
            yn(:, t1) = G_d * xn(:, t1) + J_d * p(:, t1) + v(:, t1);
            xn_true(:, t1 + 1) = A_d * xn_true(:, t1) + B_d * p(:, t1);
        end

        yn(:, end + 1) = G_d * xn(:, end) + J_d * p(:, end) + v(:, end);
    end

end
