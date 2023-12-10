clc; clear; close all;
run('CommonCommand.m');
% computer_name = getenv('COMPUTERNAME');
% 
% if strcmp(computer_name, 'SHENGYI_HP')
%     addpath(genpath("F:\git\ssm_tools\"))
%     addpath(genpath("F:\git\Function_shengyi_package\"))
%     addpath(genpath("F:\git\xihoumen_inverse_force_estimation\FEM_model\"))
%     addpath(genpath("F:\git\HHT-Tutorial\"))
% elseif strcmp(computer_name, 'mac')
%     addpath(genpath("/Users/xushengyi/Documents/GitHub/Function_shengyi_package"))
%     addpath(genpath("/Users/xushengyi/Documents/GitHub/ssm_tools_sy"))
%     addpath(genpath("/Users/xushengyi/Documents/GitHub/xihoumen_inverse_force_estimation/FEM_model"))
%     addpath(genpath("/Users/xushengyi/Documents/GitHub/HHT-Tutorial"))
% elseif strcmp(computer_name, 'ROG-SHENGYI')
%     addpath(genpath("D:\Users\xushe\Documents\GitHub\Function_shengyi_package"))
%     addpath(genpath("D:\Users\xushe\Documents\GitHub\ssm_tools"))
%     addpath(genpath("D:\git\xihoumen_inverse_force_estimation\FEM_model"))
%     addpath(genpath("D:\Users\xushe\Documents\GitHub\xihoumen_data_extract"))
%     addpath(genpath("D:\Users\xushe\Documents\GitHub\HHT-Tutorial\"))
% elseif strcmp(computer_name, 'NTNU08916')
%     addpath(genpath("C:\Users\shengyix\Documents\GitHub\Function_shengyi_package"))
%     addpath(genpath("C:\Users\shengyix\Documents\GitHub\ssm_tools_sy"))
%     addpath(genpath("C:\Users\shengyix\Documents\GitHub\xihoumen_inverse_force_estimation\FEM_model"))
%     addpath(genpath("C:\Users\shengyix\Documents\GitHub\xihoumen_data_extract"))
% elseif strcmp(computer_name, 'WIN-IUMUERP66UT')
%     addpath(genpath("C:\Users\xushengyi\Documents\Github\Function_shengyi_package"))
%     addpath(genpath("C:\Users\xushengyi\Documents\Github\Function_shengyi_package"))
%     addpath(genpath("C:\Users\xushengyi\Documents\Github\xihoumen_inverse_force_estimation\FEM_model"))
%     addpath(genpath("C:\Users\xushengyi\Documents\Github\xihoumen_data_extract"))
% else
%     error("Please add path first.")
% end
% 
% subStreamNumberDefault = 2132;
% 
% run("InitScript.m")
% 
% %% 0 绘图参数
% fig_bool = ON;
% num_figs_in_row = 6; %每一行显示几个图
% figPos = figPosX3Large; %图的大小，参数基于InitScript.m中的设置
% %设置图片间隔
% gap_between_images = [0, 0];
% figureIdx = 0;

%% 输入参数
% input.num_figs_in_row = 12; %每一行显示几个图
% input.figPos = figPos; %图的大小，参数基于InitScript.m中的设置
%设置图片间隔
% input.ON =ON;
% input.OFF =OFF;
% input.gap_between_images = [0, 0];
% input.figureIdx = 0;
% n = 13;
% n = 4;
n = 4;
[result] = viv2013(n, OFF);
startDate_global = result.startDate;
endDate_global = result.endDate + hours(1);
input.start_time = startDate_global;
input.end_time = endDate_global;

% computer_name = getenv('COMPUTERNAME');
% 
% if strcmp(computer_name, 'SHENGYI_HP')
%     input.acc_dir = "F:\test\result";
%     input.wind_dir = "F:\test\result_wind_10min";
% elseif strcmp(computer_name, 'mac')
%     input.acc_dir = "/Users/xushengyi/Documents/xihoumendata/acc";
%     input.wind_dir = "/Users/xushengyi/Documents/xihoumendata/wind";
% elseif strcmp(computer_name, 'ROG-SHENGYI')
%     input.acc_dir = "D:\xihoumendata\acc";
%     input.wind_dir = "D:\xihoumendata\wind";
% elseif strcmp(computer_name, 'ketizu')
%     input.acc_dir = "Z:\Drive\Backup\SHENGYI_HP\F\test\result";
%     input.wind_dir = "Z:\Drive\Backup\SHENGYI_HP\F\test\result_wind_10min";
% elseif strcmp(computer_name, 'NTNU08916')
%     input.acc_dir = "C:\Users\shengyix\Documents\xihoumendata\acc";
%     input.wind_dir = "C:\Users\shengyix\Documents\xihoumendata\wind";
% elseif strcmp(computer_name, 'WIN-IUMUERP66UT')
%     input.acc_dir = "Z:\Drive\Backup\SHENGYI_HP\F\test\result";
%     input.wind_dir = "Z:\Drive\Backup\SHENGYI_HP\F\test\result_wind_10min";
% else
%     error("Please add data folder first.")
% end


input.omega_0_variation = 1;
input.Q_value = 10 ^ (-4);


input.sigma_buff = 100;
input.sigma_noise = 1.0000e-04;

input.lambda_VIV = 10 ^ (-4.987547778158018);
input.sigma_p_VIV = 1e+04;
input.omega_0_variation_VIV =1;

modesel= [2,3,5,6,7,9,15,21,23,29,33,39,44,45];
% modesel = [2,3,23];

VIV_mode_seq = find(modesel ==23);
% modesel = 23;
input.modesel = modesel;
input.VIV_mode_seq = VIV_mode_seq;
nVIV = length(VIV_mode_seq);
input.nVIV = nVIV ;

%% Apply Kalman Filter
% [result_Main] = KalmanMain(input, 'showtext', true, 'showplot', false, 'filterstyle', 'fft', 'f_keep', 0.33 * [0.9, 1.1]);
[result_Main] = KalmanMain(input, 'showtext', false, 'showplot', false, 'filterstyle', 'nofilter');
logL = result_Main.logL;
logSk = result_Main.logSk;
logek = result_Main.logek;
% display logL logSk logek
dispstr = sprintf("logL = %f, logSk = %f, logek = %f", logL, logSk, logek);
disp(dispstr)

mode_deck = result_Main.mode_deck;

if 0
    %% optimization logL to get the maximum with changing lambda sigma_p omega_0_variation Q_value R_value
    % 在调用 ga 函数之前，您可以这样设置 external_params：
    external_params.modesel = [2,3,5,6,7,9,15,21,23,29,33,39,44,45];
    % external_params.modesel = [23];
    % 定义参数的范围
    lb = [-5, 10, 0.9, -10, -10]; % 这里的值是假设的，请根据您的情况进行修改
    ub = [-1, 1e5, 1.1, 1, 1]; % 这里的值也是假设的

    % 定义整数和连续变量
    IntCon = []; % 如果没有整数变量，否则提供整数变量的索引

    options = optimoptions('ga', 'MaxGenerations', 100, 'Display', 'iter', 'UseParallel', true);
    [x, fval] = ga(@(params) fitnessFunction(params, external_params), 5, [], [], [], [], lb, ub, [], IntCon, options);
    % 保存结果
    save('optimization_results.mat', 'x', 'fval');

end

%% Recalcualte the modal displacement
nmodes = result_Main.nmodes;
Result = ImportMK(nmodes, 'KMatrix.matrix', 'MMatrix.matrix', 'nodeondeck.txt', 'KMatrix.mapping', 'nodegap.txt', 'modesel', modesel);
mode_deck = Result.mode_deck; mode_deck_re = Result.mode_deck_re; node_loc = Result.node_loc; nodeondeck = Result.nodeondeck;
KMmapping = Result.Mapping;
nodegap = Result.nodegap;
mode_vec = Result.mode_vec;
loc_acc = [1403];
% loc_acc = [990.5; 1403; 1815.5];
node_shape = FindModeShapewithLocation(loc_acc, node_loc, nodeondeck, KMmapping, nodegap, mode_vec);

h_hat = result_Main.h_hat;
MM = result_Main.MM;
CC = result_Main.CC;
KK = result_Main.KK;
t_temp = result_Main.t;
t_temp = (datenum(t_temp) - datenum('1970-01-01 00:00:00')) * 86400; % Convert to seconds since epocht = result_Main.t;
t_temp = t_temp - t_temp(1); % Subtract the first element of t from all elements of t
p_filt_m = result_Main.p_filt_m;
u0 = zeros(nVIV, 1);
udot0 = zeros(nVIV, 1);
[u udot u2dot] = NewmarkInt(t_temp, MM(VIV_mode_seq,VIV_mode_seq), CC(VIV_mode_seq,VIV_mode_seq), KK(VIV_mode_seq,VIV_mode_seq), p_filt_m, 1/2, 1/4, u0, udot0);

input.u = u;
input.udot = udot;
input.u2dot = u2dot;

%% replace the force
yn = result_Main.yn;
yn_reconstruct = result_Main.yn_reconstruct;
% disp_dir = acc2dsip(yn(1, :), 50);
% F_direct = MM.*yn(1,:)/node_shape+CC.*disp_dir.vel/node_shape + KK.*disp_dir.disp/node_shape;
%
% input.p_filt_m=F_direct;

%% Caldamping ratio
fields = fieldnames(result_Main);

for i = 1:numel(fields)
    input.(fields{i}) = result_Main.(fields{i});
end

% input = result_Main;
input.ncycle = 10;

[result_Damping] = Cal_aero_damping_ratio(input, 'showplot', false, 'filterstyle', 'fft');
% [result_Damping] = Cal_aero_damping_ratio(input, 'showplot', false, 'filterstyle', 'nofilter');

amp_cell = result_Damping.amp_cell;
zeta_all_cell = result_Damping.zeta_all_cell;
top_freqs = result_Damping.top_freqs;
t_cycle_mean_cell = result_Damping.t_cycle_mean_cell;
%% read wind data

% [result] = viv2013(n, OFF);
start_time = input.start_time;
end_time = input.end_time;
wind_dir = input.wind_dir;
[result_wind] = read_wind_data(start_time, end_time, wind_dir);

t = result_Main.t;
F_filter = p_filt_m;

%% plot
if fig_bool
    % 定义总子图数量
    total_plots = 16; % 或任何你需要的子图数量
    current_plot = 1;
    num_figs_in_row = [];
    figWidthFactor = 1.5;
    figPosition = [1080*2.5,100];
    % figPosition = [100,100];
    newfigure = false;
    holdon = true;
    
    %% modal force
    for k1 = 1:nVIV
        if k1 == 1
            create_subplot(@plot, total_plots, current_plot, {t, F_filter(k1,:)}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor,'figPosition',figPosition,'newfigure',true,'firstfigure',true);
            hold on
        else
            create_subplot(@plot, total_plots, current_plot, {t, F_filter(k1,:)}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor,'figPosition',figPosition,'newfigure',newfigure,'holdon',holdon);
        end
        legend("Force from direct integration")
        title("modal force for "+"mode"+modesel(VIV_mode_seq(k1)));
        current_plot = current_plot + 1;
        ylim([-100,100])
    end
    
    save temp.mat t F_filter
    %% installation of the sensors and the rms of the sensors
    % 三个传感器位置的振型大小
    loc_acc = result_Main.loc_acc;
    loc_acc_shape = FindModeShapewithLocation(loc_acc, node_loc, nodeondeck, KMmapping, nodegap, mode_vec);
    max_loc_acc_shape = max(abs(loc_acc_shape(:,VIV_mode_seq)));
    % 三个传感器位置的振动rms大小，并基于振型大小进行缩放
    yn_rms = rms(yn, 2);
    max_acc = max(yn_rms);
    scale_factor = max_acc / max_loc_acc_shape;
    yn_rms_scaled = yn_rms / scale_factor;
    create_subplot(@plot, total_plots, current_plot, {node_loc, mode_deck(:,VIV_mode_seq)}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure);
    hold on
    create_subplot(@scatter, total_plots, current_plot, {loc_acc,  loc_acc_shape(:,VIV_mode_seq) }, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure);
    create_subplot(@scatter, total_plots, current_plot, {loc_acc,  yn_rms_scaled([1,3,5],:) }, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure);
    legend("modal shape","modal shape at the sensors location","scaled rms value of the sensors")
    title("Mode shape and the locaiton of the sensors")
    current_plot = current_plot + 1;
    xlabel('Time (s)')
    ylabel('Dimensionless displacement')

    %% reconstructed data comparison (reconstructed vs kalman filter vitural sensoring)
    create_subplot(@plot, total_plots, current_plot, {t, node_shape(VIV_mode_seq)*u, t, h_hat(15, :)}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure);
    legend("Kalman filter reconstruct", 'filtered from kalman filter')
    title("reconstructed displacement vs kalman filter vitural sensoring")
    current_plot = current_plot + 1;
    xlabel('Time (s)')
    ylabel('Displacement (m)')

    create_subplot(@plot, total_plots, current_plot, {t, node_shape(VIV_mode_seq)*udot, t, h_hat(9, :)}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure);
    legend("Kalman filter reconstruct", 'filtered from kalman filter')
    title("reconstructed velocity vs kalman filter vitural sensoring")
    current_plot = current_plot + 1;
    xlabel('Time (s)')
    ylabel('Velocity (m/s)')
    
    create_subplot(@plot, total_plots, current_plot, {t, node_shape(VIV_mode_seq)*u2dot, t, h_hat(3, :)}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure);
    legend("Kalman filter reconstruct", 'filtered from kalman filter')
    title("reconstructed acceleration vs kalman filter vitural sensoring")
    current_plot = current_plot + 1;
    xlabel('Time (s)')
    ylabel('Acceleration (m/s^2)')
    
    %% wind speed during the period
    beta_deg_mean_UA5 = result_wind.resultsTable_UA5.beta_deg_mean;
    beta_deg_mean_UA6 = result_wind.resultsTable_UA6.beta_deg_mean;
    % judge if the wind direction of both UA5 and UA6 is from 45°-225°
    for k1 = 1:length(beta_deg_mean_UA5)

        if and(beta_deg_mean_UA5(k1) > 45, beta_deg_mean_UA5(k1) < 225) && and(beta_deg_mean_UA6(k1) > 45, beta_deg_mean_UA6(k1) < 225)
            U_sel(k1) = result_wind.resultsTable_UA5.U(k1);
            % judge if the wind direction of both UA5 and UA6 is from 225°-360° or 0°-45°
        elseif or(and(beta_deg_mean_UA5(k1) > 225, beta_deg_mean_UA5(k1) < 360), and(beta_deg_mean_UA5(k1) > 0, beta_deg_mean_UA5(k1) < 45)) && or(and(beta_deg_mean_UA6(k1) > 225, beta_deg_mean_UA6(k1) < 360), and(beta_deg_mean_UA6(k1) > 0, beta_deg_mean_UA6(k1) < 45))
            U_sel(k1) = result_wind.resultsTable_UA6.U(k1);
        else

            if abs(result_wind.resultsTable_UA5.alpha_deg_mean(k1)) < abs(result_wind.resultsTable_UA6.alpha_deg_mean)
                U_sel(k1) = result_wind.resultsTable_UA5.U(k1);
            else
                U_sel(k1) = result_wind.resultsTable_UA6.U(k1);
            end

            disp("该时间点两个风速仪风向不一致，取风攻角较小的值")
        end

    end

    create_subplot(@scatter, total_plots, current_plot, {result_wind.resultsTable_UA6.Time_Start, U_sel}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure);
    xlabel('Time (s)')
    ylabel('Wind speed (m/s)')
    title("Wind speed vs. Time")
    current_plot = current_plot + 1;

    %% measured, filtered and recalcuated acceleration
    create_subplot(@plot, total_plots, current_plot, {t, yn(1, :), t, h_hat(1, :), t, yn_reconstruct(1, :)}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure);
    xlabel('Time (s)')
    ylabel('Acceleration (m/s^2)')
    title("Acceleration vs. Time")
    legend("measure", "kalman filter", "reconstruct using force from kalman filter")
    current_plot = current_plot + 1;

    create_subplot(@plot, total_plots, current_plot, {t, yn(3, :), t, h_hat(3, :), t, yn_reconstruct(3, :)}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure);
    xlabel('Time (s)')
    ylabel('Acceleration (m/s^2)')
    title("Acceleration vs. Time")
    legend("measure", "kalman filter", "reconstruct using force from kalman filter")
    current_plot = current_plot + 1;

    create_subplot(@plot, total_plots, current_plot, {t, yn(5, :), t, h_hat(5, :), t, yn_reconstruct(5, :)}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure);
    xlabel('Time (s)')
    ylabel('Acceleration (m/s^2)')
    title("Acceleration vs. Time")
    legend("measure", "kalman filter", "reconstruct using force from kalman filter")
    current_plot = current_plot + 1;

    %% Damping ratio calculation
     for k1 = 1:nVIV
    
            for k2 = 1:length(top_freqs{k1})
                % 假设 t_cycle_mean_cell{k1}{k2} 是一个包含 datetime 对象的数组
                datetimeArray = t_cycle_mean_cell{k1}{k2};
    
                % 提取第一个 datetime 对象作为参考点
                referenceDatetime = datetimeArray(1);
    
                % 计算每个 datetime 对象相对于参考点的秒数
                secondsFromReference = seconds(datetimeArray - referenceDatetime);
    
                % 现在，secondsFromReference 包含相对于第一个时间戳的秒数
    
                create_subplot(@scatter, total_plots, current_plot, {amp_cell{k1}{k2} * max(mode_deck(:, VIV_mode_seq(k1))), zeta_all_cell{k1}{k2}, [], secondsFromReference, 'filled'}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure);
                current_plot = current_plot + 1;
                % scatter(amp_cell{k1}{k2} * max(mode_deck(:, k1)), zeta_all_cell{k1}{k2}, [], secondsFromReference, 'filled');
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

end


if 0


   

    



    %% instant frequency of the filtered acceleration
    [p, fd, td] = pspectrum(h_hat(1, :), t, 'spectrogram', 'FrequencyResolution', 0.005);
    create_subplot(@instfreq, total_plots, current_plot, {p, fd, td}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure);
    ylim([0.1, 0.5])
    current_plot = current_plot + 1;
    %% amplitude dependent damping ratio

    for k1 = 1:nmodes

        for k2 = 1:length(top_freqs{k1})
            % 假设 t_cycle_mean_cell{k1}{k2} 是一个包含 datetime 对象的数组
            datetimeArray = t_cycle_mean_cell{k1}{k2};

            % 提取第一个 datetime 对象作为参考点
            referenceDatetime = datetimeArray(1);

            % 计算每个 datetime 对象相对于参考点的秒数
            secondsFromReference = seconds(datetimeArray - referenceDatetime);

            % 现在，secondsFromReference 包含相对于第一个时间戳的秒数

            create_subplot(@scatter, total_plots, current_plot, {amp_cell{k1}{k2} * max(mode_deck(:, k1)), zeta_all_cell{k1}{k2}, [], secondsFromReference, 'filled'}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure);
            current_plot = current_plot + 1;
            % scatter(amp_cell{k1}{k2} * max(mode_deck(:, k1)), zeta_all_cell{k1}{k2}, [], secondsFromReference, 'filled');
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

    %% amp vs. time
    for k1 = 1:nmodes

        for k2 = 1:length(top_freqs{k1})
            % 假设 t_cycle_mean_cell{k1}{k2} 是一个包含 datetime 对象的数组
            datetimeArray = t_cycle_mean_cell{k1}{k2};

            % 提取第一个 datetime 对象作为参考点
            referenceDatetime = datetimeArray(1);

            % 计算每个 datetime 对象相对于参考点的秒数
            secondsFromReference = seconds(datetimeArray - referenceDatetime);

            % 现在，secondsFromReference 包含相对于第一个时间戳的秒数

            create_subplot(@scatter, total_plots, current_plot, {datetimeArray, amp_cell{k1}{k2} * max(mode_deck(:, k1)), [], secondsFromReference, 'filled'}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure);
            current_plot = current_plot + 1;

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

    %% damping ratio vs time
    for k1 = 1:nmodes

        for k2 = 1:length(top_freqs{k1})
            % 假设 t_cycle_mean_cell{k1}{k2} 是一个包含 datetime 对象的数组
            datetimeArray = t_cycle_mean_cell{k1}{k2};

            % 提取第一个 datetime 对象作为参考点
            referenceDatetime = datetimeArray(1);

            % 计算每个 datetime 对象相对于参考点的秒数
            secondsFromReference = seconds(datetimeArray - referenceDatetime);

            % 现在，secondsFromReference 包含相对于第一个时间戳的秒数

            create_subplot(@scatter, total_plots, current_plot, {datetimeArray, zeta_all_cell{k1}{k2}, [], secondsFromReference, 'filled'}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure);
            current_plot = current_plot + 1;
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

    %% damping ratio and displacement with time

    data = disp_dir.disp;
    minData = min(data);
    maxData = max(data);
    normalizedData = 2 * ((data - minData) / (maxData - minData)) - 1;
    create_subplot(@plot, total_plots, current_plot, {t, normalizedData}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure);

    % plot(t, normalizedData);

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

            create_subplot(@scatter, total_plots, current_plot, {datetimeArray, (normalizedData - mean(normalizedData(end / 4:end * 3/4))) * 1000, [], secondsFromReference, 'filled'}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure);

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

    current_plot = current_plot + 1;

    %%
    for k1 = 1:nmodes

        for k2 = 1:length(top_freqs{k1})
            % 假设 t_cycle_mean_cell{k1}{k2} 是一个包含 datetime 对象的数组
            datetimeArray = t_cycle_mean_cell{k1}{k2};

            % 提取第一个 datetime 对象作为参考点
            referenceDatetime = datetimeArray(1);

            % 计算每个 datetime 对象相对于参考点的秒数
            secondsFromReference = seconds(datetimeArray - referenceDatetime);

            % 现在，secondsFromReference 包含相对于第一个时间戳的秒数

            create_subplot(@scatter, total_plots, current_plot, {amp_cell{k1}{k2} * max(mode_deck(:, k1)), zeta_all_cell{k1}{k2}, 'red'}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure);

            % 设置 colormap
            colormap('jet')
            colorbar

            hold on

            create_subplot(@scatter, total_plots, current_plot, {amp_cell_direct{k1}{k2} * max(mode_deck(:, k1)), zeta_all_cell_direct{k1}{k2}, 'green'}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure);

            % plot([0, 0.15], [-0.003, -0.003])
            % scatter(ex,epsx,'green')

            create_subplot(@scatter, total_plots, current_plot, {amp_temp * max(mode_deck(:, k1)), zetam - 0.3/100, 'blue'}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure);

            str = "Mode : %d, Frequency : %.2f Hz";
            title(sprintf(str, modesel(k1), top_freqs{k1}(k2)));
            xlim([0.05, 0.12])
            ylim([-0.5, 0.5] / 100)
            xlabel("Amplitude(m)")
            ylabel("Damping ratio")
            legend("use force kalman", 'use force direct integral', "use acc")
        end

    end

    current_plot = current_plot + 1;
end

if 0
    %% reconstructed data comparison (direct integral and kalman filter)
    [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
    plot(t, u_direct)
    hold on
    plot(t, u)
    plot(t, h_hat(5, :) / node_shape)
    legend("direct integral reconstruct", "filter reconstruct", 'filtered')
    title("reconstructed displacement")

    [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
    plot(t, udot_direct)
    hold on
    plot(t, udot)
    plot(t, h_hat(3, :) / node_shape)
    legend("direct integral reconstruct", "filter reconstruct", 'filtered')
    title("reconstructed velocity")

    [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
    plot(t, u2dot_direct)
    hold on
    plot(t, u2dot)
    plot(t, h_hat(1, :) / node_shape)
    legend("direct integral reconstruct", "filter reconstruct", 'filtered')
    title("reconstructed acceleration")

    %% filtered and recalcuated displacement
    [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
    plot(t_temp, x_filt_original(1, :))
    hold on
    plot(t_temp, u)
    legend("filtered", "recalculate")
    title("u")

    %% filtered and recalcuated velocity
    [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
    plot(t_temp, x_filt_original(2, :))
    hold on
    plot(t_temp, udot)
    legend("filtered", "recalculate")
    title("udot")

    %% filtered and recalcuated acceleration
    [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
    plot(t_temp, h_hat(1, :) / node_shape)
    hold on
    plot(t_temp, u2dot)
    legend("filtered", "recalculate")
    title("u2dot")

    %% filtered force
    [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
    plot(t_temp, p_filt_m)
    hold on
    legend("filtered")
    title("filtered viv force")
    %% wind speed during the period
    [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
    beta_deg_mean_UA5 = result_wind.resultsTable_UA5.beta_deg_mean;
    beta_deg_mean_UA6 = result_wind.resultsTable_UA6.beta_deg_mean;
    % judge if the wind direction of both UA5 and UA6 is from 45°-225°
    for k1 = 1:length(beta_deg_mean_UA5)

        if and(beta_deg_mean_UA5(k1) > 45, beta_deg_mean_UA5(k1) < 225) && and(beta_deg_mean_UA6(k1) > 45, beta_deg_mean_UA6(k1) < 225)
            U_sel(k1) = result_wind.resultsTable_UA5.U(k1);
            % judge if the wind direction of both UA5 and UA6 is from 225°-360° or 0°-45°
        elseif or(and(beta_deg_mean_UA5(k1) > 225, beta_deg_mean_UA5(k1) < 360), and(beta_deg_mean_UA5(k1) > 0, beta_deg_mean_UA5(k1) < 45)) && or(and(beta_deg_mean_UA6(k1) > 225, beta_deg_mean_UA6(k1) < 360), and(beta_deg_mean_UA6(k1) > 0, beta_deg_mean_UA6(k1) < 45))
            U_sel(k1) = result_wind.resultsTable_UA6.U(k1);
        else

            if abs(result_wind.resultsTable_UA5.alpha_deg_mean(k1)) < abs(result_wind.resultsTable_UA6.alpha_deg_mean)
                U_sel(k1) = result_wind.resultsTable_UA5.U(k1);
            else
                U_sel(k1) = result_wind.resultsTable_UA6.U(k1);
            end

            disp("该时间点两个风速仪风向不一致，取风攻角较小的值")
        end

    end

    scatter(result_wind.resultsTable_UA6.Time_Start, U_sel);
    % scatter(result_wind.resultsTable_UA6.Time_Start, result_wind.resultsTable_UA6.U);
    % hold on
    % scatter(result_wind.resultsTable_UA5.Time_Start, result_wind.resultsTable_UA5.U);
    xlabel('Time (s)')
    ylabel('Wind speed (m/s)')
    title("Wind speed vs. Time")

    %% measured, filtered and recalcuated acceleration
    [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
    plot(t, yn(1, :));
    hold on
    plot(t, h_hat(1, :));
    plot(t, yn_reconstruct(1, :))
    xlabel('Time (s)')
    ylabel('Acceleration (m/s^2)')
    title("Acceleration vs. Time")
    legend("measure", "filtered", "reconstruct")

    %% measured(direct intergal), filtered displacement
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
    legend("direct", "filtered")

    %% instant frequency of the filtered acceleration
    [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
    [p, fd, td] = pspectrum(h_hat(1, :), t, 'spectrogram', 'FrequencyResolution', 0.005);
    instfreq(p, fd, td);
    ylim([0.1, 0.5])

    %% amplitude dependent damping ratio
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

    %% amp vs. time
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

    %% damping ratio vs time
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

    %% damping ratio and displacement with time
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

    %% calculated by ee instant frequency
    % [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
    % plot(t, frex)
    % xlabel('Time (s)')
    % ylabel('Frequency (Hz)')
    % title("Frequency vs. Time calculated by ee")

    %%
    % [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
    % plot(t, ex)
    % xlabel('Time (s)')
    % ylabel('Amplitude (m)')
    % title("Amplitude vs. Time calculated by ee")

    %%
    % [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
    % scatter(ex, epsx, 'green')
    % xlim([0.05, 0.12])
    % ylim([-0.5, 0.5] / 100)
    % xlabel('Amplitude (m)')
    % ylabel('Damping ratio')
    % title("Damping ratio vs. Amplitude calculated by ee")

    %%
    % for k1 = 1:nmodes
    %
    % for k2 = 1:length(top_freqs{k1})
    %     % 假设 t_cycle_mean_cell{k1}{k2} 是一个包含 datetime 对象的数组
    %     datetimeArray = t_cycle_mean_cell{k1}{k2};
    %
    %     % 提取第一个 datetime 对象作为参考点
    %     referenceDatetime = datetimeArray(1);
    %
    %     % 计算每个 datetime 对象相对于参考点的秒数
    %     secondsFromReference = seconds(datetimeArray - referenceDatetime);
    %
    %     % 现在，secondsFromReference 包含相对于第一个时间戳的秒数
    %
    %     [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
    %     % scatter(amp_cell{k1}{k2} * max(mode_deck(:, k1)), zeta_all_cell{k1}{k2}, [], secondsFromReference, 'filled');
    %     scatter(amp_cell{k1}{k2} * max(mode_deck(:, k1)),zeta_all_cell{k1}{k2},'red');
    %     % 设置 colormap
    %     colormap('jet')
    %     colorbar
    %
    %     hold on
    %     % plot([0, 0.15], [-0.003, -0.003])
    %     % scatter(ex,epsx,'green')
    %     scatter(amp_temp* max(mode_deck(:, k1)),zetam-0.3/100,'blue')
    %     scatter(ex, epsx-0.3/100, 'green')
    %     str = "Mode : %d, Frequency : %.2f Hz";
    %     title(sprintf(str, modesel(k1), top_freqs{k1}(k2)));
    %     xlim([0.05, 0.12])
    %     ylim([-0.5, 0.5] / 100)
    %     xlabel("Amplitude(m)")
    %     ylabel("Damping ratio")
    %     legend("use force","use acc","use ee")
    % end
    %
    % end

    %%
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
            % scatter(amp_cell{k1}{k2} * max(mode_deck(:, k1)), zeta_all_cell{k1}{k2}, [], secondsFromReference, 'filled');
            scatter(amp_cell{k1}{k2} * max(mode_deck(:, k1)), zeta_all_cell{k1}{k2}, 'red');
            % 设置 colormap
            colormap('jet')
            colorbar

            hold on
            scatter(amp_cell_direct{k1}{k2} * max(mode_deck(:, k1)), zeta_all_cell_direct{k1}{k2}, 'green');

            % plot([0, 0.15], [-0.003, -0.003])
            % scatter(ex,epsx,'green')
            scatter(amp_temp * max(mode_deck(:, k1)), zetam - 0.3/100, 'blue')
            str = "Mode : %d, Frequency : %.2f Hz";
            title(sprintf(str, modesel(k1), top_freqs{k1}(k2)));
            xlim([0.05, 0.12])
            ylim([-0.5, 0.5] / 100)
            xlabel("Amplitude(m)")
            ylabel("Damping ratio")
            legend("use force kalman", 'use force direct integral', "use acc")
        end

    end

end

% %% functions
% function logL = fitnessFunction(params, external_params)
%     input.lambda = 10 ^ params(1);
%     input.sigma_p = params(2);
%     input.omega_0_variation = params(3);
%     input.Q_value = 10 ^ params(4);
%     input.R_value = 10 ^ params(5);

%     input.modesel = external_params.modesel;
%     %  figPos = external_params.figPos;
%     % ON = external_params.ON;
%     % OFF = external_params.OFF;

%     % input.num_figs_in_row = 12; %每一行显示几个图
%     % input.figPos = figPos; %图的大小，参数基于InitScript.m中的设置
%     %设置图片间隔
%     % input.ON =ON;
%     % input.OFF =OFF;
%     % input.gap_between_images = [0, 0];
%     % input.figureIdx = 0;
%     n = 4;
%     [result] = viv2013(n, false);
%     input.start_time = result.startDate;
%     input.end_time = result.endDate;
%     % input.acc_dir = "F:\test\result";
%     input.acc_dir = "Z:\Drive\Backup\SHENGYI_HP\F\test\result";

%     result_Main = KalmanMain(input, 'showtext', false, 'showplot', false, 'filterstyle', 'fft', 'f_keep', 0.33 * [0.9, 1.1]);
%     logL = -result_Main.logL; % 因为 ga 试图最小化函数，所以取负数
% end

%% functions
function target = fitnessFunction(params, external_params)
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
    fields = fieldnames(result_Main);

    for i = 1:numel(fields)
        input.(fields{i}) = result_Main.(fields{i});
    end

    input.ncycle = 10;
    [result_Damping] = Cal_aero_damping_ratio(input, 'showplot', false, 'filterstyle', 'nofilter');
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

    zeta_all_cell = result_Damping.zeta_all_cell;
    zeta1 = zeta_all_cell{1}{1};
    zeta2 = zetam - 0.3/100;

    target = norm(zeta1 - zeta2);
    logL = target;

end

function [amp, fre] = ee(data, dt)
    %EE Summary of this function goes here
    %   Detailed explanation goes here
    [normalizeddata, amp] = splinenormalizeep(data);
    frecarrier = gradient(normalizeddata, dt);
    [normalizedfrecarrier, ampfrecarrier] = splinenormalizeep(frecarrier);
    fre = ampfrecarrier / 2 / pi;
end
