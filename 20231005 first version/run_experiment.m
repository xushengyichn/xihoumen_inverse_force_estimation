function [results_experiment] = run_experiment(input, varargin)

    if nargin == 0
        clc; clear; close all;
        addpath(genpath("F:\git\ssm_tools\"))
        addpath(genpath("F:\git\Function_shengyi_package\"))
        addpath(genpath("F:\git\xihoumen_inverse_force_estimation\FEM_model\"))
        addpath(genpath("F:\git\HHT-Tutorial\"))
        addpath(genpath("/Users/xushengyi/Documents/GitHub/Function_shengyi_package"))
        addpath(genpath("/Users/xushengyi/Documents/GitHub/ssm_tools_sy"))
        addpath(genpath("/Users/xushengyi/Documents/GitHub/xihoumen_inverse_force_estimation/FEM_model"))
        addpath(genpath("/Users/xushengyi/Documents/GitHub/HHT-Tutorial"))
        subStreamNumberDefault = 2132;

        params = Init_fun();

        %% 0 绘图参数
        fig_bool = params.OFF;
        num_figs_in_row = 12; %每一行显示几个图
        figPos = params.figPosSmall; %图的大小，参数基于InitScript.m中的设置
        %设置图片间隔
        gap_between_images = [0, 0];
        figureIdx = 0;

        %% 输入参数
        % input.num_figs_in_row = 12; %每一行显示几个图
        % input.figPos = figPos; %图的大小，参数基于InitScript.m中的设置
        %设置图片间隔
        % input.ON =ON;
        % input.OFF =OFF;
        % input.gap_between_images = [0, 0];
        % input.figureIdx = 0;
        n = 4;
        [result] = viv2013(n, params.OFF);
        startDate_global = result.startDate;
        endDate_global = result.endDate;
        input.start_time = startDate_global;
        input.end_time = endDate_global;
        % input.acc_dir = "/Users/xushengyi/Documents/xihoumendata/acc";
        input.acc_dir = "F:\test\result";
        % input.acc_dir = "Z:\Drive\Backup\SHENGYI_HP\F\test\result";
        % input.wind_dir = "/Users/xushengyi/Documents/xihoumendata/wind";
        input.wind_dir = "F:\test\result_wind_10min";
        % input.lambda = 10 ^ (-1);
        % input.sigma_p = 10000;
        % input.omega_0_variation =1;
        % input.Q_value =10 ^ (-8);
        % input.R_value = 10 ^ (-6);

        input.lambda = 10 ^ (-1.001242541230301);
        % input.sigma_p = 9.063096667830060e+04;
        % input.omega_0_variation =0.902647415472734;
        % input.Q_value =10 ^ (-1.016211706804576);
        % input.R_value = 10 ^ (-1.003148221125874);

        % input.lambda = 10 ^ (-4.934808796013671);
        % input.sigma_p = 6.895548550856822e+03;
        % input.omega_0_variation =1.097383030422062;
        % input.Q_value =10 ^ (-9.633948257379021);
        % input.R_value = 10 ^ (-2.415076745081128);

        % input.lambda = 10 ^ (-4.993486819657864);
        input.sigma_p = 1.441514803767596e+04;
        input.omega_0_variation = 0.904969462898074;
        input.Q_value = 10 ^ (-9.901777612793937);
        input.R_value = 10 ^ (-3.866588296864785);

        % modesel= [2,3,5,6,7,9,15,21,23,29,33,39,44,45];
        modesel = 23;
        input.modesel = modesel;

        results_experiment = run_experiment(input, 'showtext', false, 'showplot', false);
        return;
    end

    p = inputParser;
    addParameter(p, 'showtext', true, @islogical)
    addParameter(p, 'showplot', false, @islogical)
    addParameter(p, 'caldamp', true, @islogical)
    addParameter(p, 'caldamp_recalculated_v', false, @islogical)
    addParameter(p, 'num_figs_in_row', 12, @isnumeric)
    addParameter(p, 'figPos', [100, 100, 400, 300], @isnumeric)
    addParameter(p, 'gap_between_images', [0, 0], @isnumeric)
    addParameter(p, 'figureIdx', 0, @isnumeric)
    parse(p, varargin{:});
    showtext = p.Results.showtext;
    showplot = p.Results.showplot;
    num_figs_in_row = p.Results.num_figs_in_row;
    figPos = p.Results.figPos;
    gap_between_images = p.Results.gap_between_images;
    figureIdx = p.Results.figureIdx;
    caldamp = p.Results.caldamp;
    caldamp_recalculated_v = p.Results.caldamp_recalculated_v;

    [result_Main] = KalmanMain(input, 'showtext', showtext, 'showplot', showplot, 'filterstyle', 'fft', 'f_keep', 0.33 * [0.9, 1.1]);
    logL = result_Main.logL;
    logSk = result_Main.logSk;
    logek = result_Main.logek;
    invSk = result_Main.invSk;
    % display logL logSk logek
    if showtext
        dispstr = sprintf("logL = %f, logSk = %f, logek = %f", logL, logSk, logek);
        disp(dispstr)
    end

    mode_deck = result_Main.mode_deck;



    %% Recalcualte the modal displacement
    nmodes = result_Main.nmodes;
    Result = ImportMK(nmodes, 'KMatrix.matrix', 'MMatrix.matrix', 'nodeondeck.txt', 'KMatrix.mapping', 'nodegap.txt', 'modesel', [23], 'showtext', showtext);
    mode_deck = Result.mode_deck; mode_deck_re = Result.mode_deck_re; node_loc = Result.node_loc; nodeondeck = Result.nodeondeck;
    eig_vec = Result.eig_vec;
    Mapping_data = Result.Mapping;
    loc_acc = [1403];
    acc_node = FindNodewithLocation(loc_acc, node_loc, nodeondeck);
    acc_node_list = reshape(permute(acc_node, [2 1]), [], 1); % 交错重塑
    acc_matrix_seq = node2matrixseq(acc_node_list, Mapping_data);
    node_shape = mean(eig_vec(acc_matrix_seq));

    h_hat = result_Main.h_hat;
    MM = result_Main.MM;
    CC = result_Main.CC;
    KK = result_Main.KK;
    t_temp = result_Main.t;
    t_temp = (datenum(t_temp) - datenum('1970-01-01 00:00:00')) * 86400; % Convert to seconds since epocht = result_Main.t;
    t_temp = t_temp - t_temp(1); % Subtract the first element of t from all elements of t
    p_filt_m = result_Main.p_filt_m;
    [u udot u2dot] = NewmarkInt(t_temp, MM, CC, KK, p_filt_m, 1/2, 1/4, 0, 0);
    x_filt_original = result_Main.x_filt_original;

    input.u = u;
    input.udot = udot;
    input.u2dot = u2dot;

    %% Caldamping ratio
    fields = fieldnames(result_Main);

    for i = 1:numel(fields)
        input.(fields{i}) = result_Main.(fields{i});
    end

    if caldamp
        % input = result_Main;
        input.ncycle = 10;

        [result_Damping] = Cal_aero_damping_ratio(input, 'showplot', false, 'filterstyle', 'nofilter');

        %% read wind data

        % [result] = viv2013(n, OFF);
        start_time = input.start_time;
        end_time = input.end_time;
        wind_dir = input.wind_dir;

        % wind_dir = "F:\test\result_wind_10min";
        [result_wind] = read_wind_data(start_time, end_time, wind_dir);

        %% compare with ee
        yn = result_Main.yn;
        fs = 50;
        t = result_Main.t;
        disp_dir = acc2dsip(yn(1, :), 50);
        [ex, frex] = ee(disp_dir.disp, 1 / fs); %经验包络法求瞬时频率和瞬时振幅 % empirical envelope method to find instantaneous frequency and instantaneous amplitude

        omgx = frex * 2 * pi;

        for k1 = 1:length(ex) - 1
            epsx(k1) = log(ex(k1) / ex(k1 + 1)) / omgx(k1) * fs;
        end

        epsx = [epsx epsx(end)];

        % figure
        % plot(t,epsx)
        % grid

        %% compare with damping ratio calculated by acc

        amp_cell = result_Damping.amp_cell;
        t_cycle_mean_cell = result_Damping.t_cycle_mean_cell;
        amp_temp = amp_cell{1}{1};
        t_cycle_mean_temp = t_cycle_mean_cell{1}{1};
        m_cycle = input.ncycle; %cycles to be averaged
        zetam = zeros(1, length(t_cycle_mean_temp)); % Pre-allocate zetam with zeros

        for k1 = 1:length(t_cycle_mean_temp)

            if k1 <= m_cycle / 2 % Beginning boundary
                start_idx = 1;
                end_idx = start_idx + m_cycle;
            elseif k1 > length(t_cycle_mean_temp) - m_cycle / 2 % Ending boundary
                end_idx = length(t_cycle_mean_temp);
                start_idx = end_idx - m_cycle;
            else % Middle
                start_idx = k1 - floor(m_cycle / 2);
                end_idx = k1 + floor(m_cycle / 2);
            end

            deltam = log(amp_temp(start_idx) / amp_temp(end_idx));

            zetam(k1) = sqrt(deltam ^ 2 / (4 * m_cycle ^ 2 * pi ^ 2 + deltam ^ 2));

            if deltam > 0
                zetam(k1) = abs(zetam(k1));
            else
                zetam(k1) = -abs(zetam(k1));
            end

        end

    end

    if caldamp_recalculated_v

        [result_Damping_recalculated_v] = Cal_aero_damping_ratio(input, 'showplot', false, 'filterstyle', 'nofilter');
        
    end

    results_experiment = struct();
    results_experiment.result_Main = result_Main;
    
    results_experiment.h_hat = h_hat;
    results_experiment.nmodes = nmodes;
    results_experiment.yn_reconstruct = result_Main.yn_reconstruct;
    results_experiment.node_shape = node_shape;
    results_experiment.invSk = invSk;
    
    results_experiment.u = u;
    results_experiment.udot = udot;
    results_experiment.u2dot = u2dot;
    results_experiment.x_filt_original = x_filt_original;
    results_experiment.t_temp = t_temp;
    results_experiment.p_filt_m = p_filt_m;
    results_experiment.MM = MM;
    results_experiment.CC = CC;
    results_experiment.KK = KK;
    results_experiment.acc_matrix_seq = acc_matrix_seq;
    results_experiment.acc_node_list = acc_node_list;
    results_experiment.acc_node = acc_node;
    results_experiment.node_loc = node_loc;
    results_experiment.nodeondeck = nodeondeck;
    results_experiment.eig_vec = eig_vec;
    results_experiment.Mapping_data = Mapping_data;
    results_experiment.loc_acc = loc_acc;
    results_experiment.mode_deck = mode_deck;
    results_experiment.mode_deck_re = mode_deck_re;

    if caldamp
        results_experiment.result_wind = result_wind;
        results_experiment.work_cell = result_Damping.work_cell;
        results_experiment.zetam = zetam;
        results_experiment.result_Damping = result_Damping;
        results_experiment.amp_cell = result_Damping.amp_cell;
        results_experiment.zeta_all_cell = result_Damping.zeta_all_cell;
        results_experiment.top_freqs = result_Damping.top_freqs;
        results_experiment.t_cycle_mean_cell = result_Damping.t_cycle_mean_cell;
        results_experiment.disp_dir = disp_dir;
        results_experiment.ex = ex;
        results_experiment.epsx = epsx;
        results_experiment.yn = yn;
        results_experiment.t = t;
    end

    if caldamp_recalculated_v
        

    end    

    %% plot
    if showplot
        [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
        plot(t_temp, x_filt_original(1, :))
        hold on
        plot(t_temp, u)
        legend("filtered", "recalculate")
        title("u")

        [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
        plot(t_temp, x_filt_original(2, :))
        hold on
        plot(t_temp, udot)
        legend("filtered", "recalculate")
        title("udot")

        [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
        plot(t_temp, h_hat(1, :) / node_shape)
        hold on
        plot(t_temp, u2dot)
        legend("filtered", "recalculate")
        title("u2dot")

        t = result_Main.t;
        yn = result_Main.yn;
        h_hat = result_Main.h_hat;
        nmodes = result_Main.nmodes;
        amp_cell = result_Damping.amp_cell;
        zeta_all_cell = result_Damping.zeta_all_cell;
        top_freqs = result_Damping.top_freqs;
        t_cycle_mean_cell = result_Damping.t_cycle_mean_cell;
        yn_reconstruct = result_Main.yn_reconstruct;

        [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
        scatter(result_wind.resultsTable_UA6.Time_Start, result_wind.resultsTable_UA6.U);
        xlabel('Time (s)')
        ylabel('Wind speed (m/s^2)')
        title("Wind speed vs. Time")

        [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
        plot(t, yn(1, :));
        hold on
        plot(t, h_hat(1, :));
        plot(t, yn_reconstruct(1, :))
        xlabel('Time (s)')
        ylabel('Acceleration (m/s^2)')
        title("Acceleration vs. Time")
        legend("measure", "filtered", "reconstruct")

        [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
        % Result = ImportMK(nmodes, 'KMatrix.matrix', 'MMatrix.matrix', 'nodeondeck.txt', 'KMatrix.mapping', 'nodegap.txt', 'modesel', [23]);
        % mode_deck = Result.mode_deck; mode_deck_re = Result.mode_deck_re; node_loc = Result.node_loc;nodeondeck = Result.nodeondeck;
        % eig_vec=Result.eig_vec;
        % Mapping_data = Result.Mapping;
        % loc_acc = [1403];
        % acc_node=FindNodewithLocation(loc_acc, node_loc, nodeondeck);
        % acc_node_list = reshape(permute(acc_node, [2 1]), [], 1); % 交错重塑
        % acc_matrix_seq = node2matrixseq(acc_node_list, Mapping_data);
        % node_shape = mean(eig_vec(acc_matrix_seq));
        plot(t, disp_dir.disp);
        hold on
        plot(t, h_hat(3, :));

        xlabel('Time (s)')
        ylabel('Displacement (m)')
        title("Displacement vs. Time")
        legend("measure", "filtered")

        [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
        [p, fd, td] = pspectrum(h_hat(1, :), t, 'spectrogram', 'FrequencyResolution', 0.005);
        instfreq(p, fd, td);
        ylim([0.1, 0.5])

        for k1 = 1:nmodes

            for k2 = 1:length(top_freqs{k1})
                % 假设 t_cycle_mean_cell{k1}{k2} 是一个包含 datetime 对象的数组
                datetimeArray = t_cycle_mean_cell{k1}{k2};

                % 提取第一个 datetime 对象作为参考点
                referenceDatetime = datetimeArray(1);

                % 计算每个 datetime 对象相对于参考点的秒数
                secondsFromReference = seconds(datetimeArray - referenceDatetime);

                % 现在，secondsFromReference 包含相对于第一个时间戳的秒数

                [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
                scatter(amp_cell{k1}{k2} * max(mode_deck(:, k1)), zeta_all_cell{k1}{k2}, [], secondsFromReference, 'filled');
                % 设置 colormap
                colormap('jet')
                colorbar

                hold on
                plot([0, 0.15], [-0.003, -0.003])
                % scatter(ex,epsx,'green')
                str = "Mode : %d, Frequency : %.2f Hz";
                title(sprintf(str, modesel(k1), top_freqs{k1}(k2)));
                xlim([0.05, 0.12])
                ylim([-0.5, 0.5] / 100)
                xlabel("Amplitude(m)")
                ylabel("Damping ratio")
            end

        end

        for k1 = 1:nmodes

            for k2 = 1:length(top_freqs{k1})
                % 假设 t_cycle_mean_cell{k1}{k2} 是一个包含 datetime 对象的数组
                datetimeArray = t_cycle_mean_cell{k1}{k2};

                % 提取第一个 datetime 对象作为参考点
                referenceDatetime = datetimeArray(1);

                % 计算每个 datetime 对象相对于参考点的秒数
                secondsFromReference = seconds(datetimeArray - referenceDatetime);

                % 现在，secondsFromReference 包含相对于第一个时间戳的秒数

                [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
                scatter(datetimeArray, amp_cell{k1}{k2} * max(mode_deck(:, k1)), [], secondsFromReference, 'filled');
                % 设置 colormap
                colormap('jet')
                % colorbar

                % hold on
                % plot([0,0.15],[-0.003,-0.003])

                str = "Mode : %d, Frequency : %.2f Hz, Amp vs. Time";
                title(sprintf(str, modesel(k1), top_freqs{k1}(k2)));
                % xlim([0,0.15])
                % ylim([-25,25]/100)
                xlabel("Time(s)")
                ylabel("Amplitude")
            end

        end

        for k1 = 1:nmodes

            for k2 = 1:length(top_freqs{k1})
                % 假设 t_cycle_mean_cell{k1}{k2} 是一个包含 datetime 对象的数组
                datetimeArray = t_cycle_mean_cell{k1}{k2};

                % 提取第一个 datetime 对象作为参考点
                referenceDatetime = datetimeArray(1);

                % 计算每个 datetime 对象相对于参考点的秒数
                secondsFromReference = seconds(datetimeArray - referenceDatetime);

                % 现在，secondsFromReference 包含相对于第一个时间戳的秒数

                [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
                scatter(datetimeArray, zeta_all_cell{k1}{k2}, [], secondsFromReference, 'filled');
                % 设置 colormap
                colormap('jet')
                % colorbar

                % hold on
                % plot([datetimeArray(1),datetimeArray(end)],[-0.003,-0.003])

                str = "Mode : %d, Frequency : %.2f Hz";
                title(sprintf(str, modesel(k1), top_freqs{k1}(k2)));
                % xlim([0.05,0.12])
                hold on
                plot([datetimeArray(1), datetimeArray(end)], [-0.003, -0.003])

                ylim([-0.5, 0] / 100)
                xlabel("Time(s)")
                ylabel("Damping ratio")

            end

        end

        [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);

        data = disp_dir.disp;
        minData = min(data);
        maxData = max(data);
        normalizedData = 2 * ((data - minData) / (maxData - minData)) - 1;
        plot(t, normalizedData);

        xlabel('Time (s)')
        ylabel('Displacement (m)')
        title("Displacement vs. Time")
        legend("measure", "filtered")
        hold on

        for k1 = 1:nmodes

            for k2 = 1:length(top_freqs{k1})
                % 假设 t_cycle_mean_cell{k1}{k2} 是一个包含 datetime 对象的数组
                datetimeArray = t_cycle_mean_cell{k1}{k2};

                % 提取第一个 datetime 对象作为参考点
                referenceDatetime = datetimeArray(1);

                % 计算每个 datetime 对象相对于参考点的秒数
                secondsFromReference = seconds(datetimeArray - referenceDatetime);

                % 现在，secondsFromReference 包含相对于第一个时间戳的秒数
                data = zeta_all_cell{k1}{k2};
                minData = min(data);
                maxData = max(data);
                normalizedData = 2 * ((data - minData) / (maxData - minData)) - 1;
                scatter(datetimeArray, (normalizedData - mean(normalizedData(end / 4:end * 3/4))) * 1000, [], secondsFromReference, 'filled');
                % 设置 colormap
                colormap('jet')
                % colorbar

                % hold on
                % plot([datetimeArray(1),datetimeArray(end)],[-0.003,-0.003])

                str = "Mode : %d, Frequency : %.2f Hz";
                title(sprintf(str, modesel(k1), top_freqs{k1}(k2)));
                % xlim([0.05,0.12])
                ylim([-1, 1])
                xlabel("Time(s)")
                ylabel("Damping ratio")

            end

        end

        hold off
        [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
        plot(t, frex)
        xlabel('Time (s)')
        ylabel('Frequency (Hz)')
        title("Frequency vs. Time calculated by ee")

        [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
        plot(t, ex)
        xlabel('Time (s)')
        ylabel('Amplitude (m)')
        title("Amplitude vs. Time calculated by ee")
        [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
        scatter(ex, epsx, 'green')
        xlim([0.05, 0.12])
        ylim([-0.5, 0.5] / 100)
        xlabel('Amplitude (m)')
        ylabel('Damping ratio')
        title("Damping ratio vs. Amplitude calculated by ee")

    end

end

%% functions
function logL = fitnessFunction(params, external_params)
    input.lambda = 10 ^ params(1);
    input.sigma_p = params(2);
    input.omega_0_variation = params(3);
    input.Q_value = 10 ^ params(4);
    input.R_value = 10 ^ params(5);

    input.modesel = external_params.modesel;
    %  figPos = external_params.figPos;
    % ON = external_params.ON;
    % OFF = external_params.OFF;

    % input.num_figs_in_row = 12; %每一行显示几个图
    % input.figPos = figPos; %图的大小，参数基于InitScript.m中的设置
    %设置图片间隔
    % input.ON =ON;
    % input.OFF =OFF;
    % input.gap_between_images = [0, 0];
    % input.figureIdx = 0;
    n = 4;
    [result] = viv2013(n, false);
    input.start_time = result.startDate;
    input.end_time = result.endDate;
    % input.acc_dir = "F:\test\result";
    input.acc_dir = "Z:\Drive\Backup\SHENGYI_HP\F\test\result";

    result_Main = KalmanMain(input, 'showtext', false, 'showplot', false, 'filterstyle', 'fft', 'f_keep', 0.33 * [0.9, 1.1]);
    logL = -result_Main.logL; % 因为 ga 试图最小化函数，所以取负数
end

function [amp, fre] = ee(data, dt)
    %EE Summary of this function goes here
    %   Detailed explanation goes here
    [normalizeddata, amp] = splinenormalizeep(data);
    frecarrier = gradient(normalizeddata, dt);
    [normalizedfrecarrier, ampfrecarrier] = splinenormalizeep(frecarrier);
    fre = ampfrecarrier / 2 / pi;
end
