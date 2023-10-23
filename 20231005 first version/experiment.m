%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn xushengyichn@outlook.com
%Date: 2023-10-20 12:45:02
%LastEditors: xushengyichn xushengyichn@outlook.com
%LastEditTime: 2023-10-20 15:15:00
%FilePath: /ssm_tools_sy/Users/xushengyi/Documents/GitHub/xihoumen_inverse_force_estimation/20231005 first version/experiment.m
%Description: 由于选择不同参数会对阻尼比计算结果产生影响，因此需要进行参数选择实验
%
%Copyright (c) 2023 by ${git_name_email}, All Rights Reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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

run("InitScript.m")
%% 0 绘图参数
fig_bool = ON;
num_figs_in_row = 12; %每一行显示几个图
figPos = figPosSmall; %图的大小，参数基于InitScript.m中的设置
%设置图片间隔
gap_between_images = [0, 0];
figureIdx = 0;

n = 4;
[result] = viv2013(n, OFF);
startDate_global = result.startDate;
endDate_global = result.endDate;
input.start_time = startDate_global;
input.end_time = endDate_global;
% input.acc_dir = "/Users/xushengyi/Documents/xihoumendata/acc";
% input.wind_dir = "/Users/xushengyi/Documents/xihoumendata/wind";
input.acc_dir = "F:\test\result";
input.wind_dir = "F:\test\result_wind_10min";


experiment_names = {'exp1', 'exp2', 'exp3','exp4'}; % 定义实验的名称
% 定义每个实验的参数
parameters = {
    {10^(-1.001242541230301), 1.441514803767596e+04, 0.904969462898074, 10^(-9.901777612793937), 10^(-3.866588296864785)},
    {10^(-2.001242541230301), 1.441514803767596e+04, 0.904969462898074, 10^(-9.901777612793937), 10^(-3.866588296864785)},
    {10^(-3.001242541230301), 1.441514803767596e+04, 0.904969462898074, 10^(-9.901777612793937), 10^(-3.866588296864785)},
    {10^(-4.001242541230301), 1.441514803767596e+04, 0.904969462898074, 10^(-9.901777612793937), 10^(-3.866588296864785)},
    % {10 ^ (-1.001242541230301), 9.063096667830060e+04, 0.902647415472734, 10 ^ (-1.016211706804576), 10 ^ (-1.003148221125874)},
    % {10 ^ (-4.934808796013671), 6.895548550856822e+03, 1.097383030422062, 10 ^ (-9.633948257379021), 10 ^ (-2.415076745081128)},
    % {10 ^ (-4.993486819657864), 1.441514803767596e+04, 0.904969462898074, 10^(-9.901777612793937), 10^(-3.866588296864785)},
    % 在这里添加其他实验的参数
};


% parameters = {
%     {10^(-1.001242541230301), 1.441514803767596e+02, 0.904969462898074, 10^(-9.901777612793937), 10^(-3.866588296864785)},
%     {10^(-1.001242541230301), 1.441514803767596e+03, 0.904969462898074, 10^(-9.901777612793937), 10^(-3.866588296864785)},
%     {10^(-1.001242541230301), 1.441514803767596e+04, 0.904969462898074, 10^(-9.901777612793937), 10^(-3.866588296864785)},
%     {10^(-1.001242541230301), 1.441514803767596e+05, 0.904969462898074, 10^(-9.901777612793937), 10^(-3.866588296864785)},
%     % 在这里添加其他实验的参数
% };


% parameters = {
%     {10^(-1.001242541230301), 1.441514803767596e+04, 0.9, 10^(-9.901777612793937), 10^(-3.866588296864785)},
%     {10^(-1.001242541230301), 1.441514803767596e+04, 1, 10^(-9.901777612793937), 10^(-3.866588296864785)},
%     {10^(-1.001242541230301), 1.441514803767596e+04, 1.1, 10^(-9.901777612793937), 10^(-3.866588296864785)},
%     {10^(-1.001242541230301), 1.441514803767596e+04, 1.2, 10^(-9.901777612793937), 10^(-3.866588296864785)},
%     % 在这里添加其他实验的参数
% };


% parameters = {
%     {10^(-1.001242541230301), 1.441514803767596e+04, 0.904969462898074, 10^(-9), 10^(-3.866588296864785)},
%     {10^(-1.001242541230301), 1.441514803767596e+04, 0.904969462898074, 10^(-8), 10^(-3.866588296864785)},
%     {10^(-1.001242541230301), 1.441514803767596e+04, 0.904969462898074, 10^(-7), 10^(-3.866588296864785)},
%     {10^(-1.001242541230301), 1.441514803767596e+04, 0.904969462898074, 10^(-6), 10^(-3.866588296864785)},
%     % 在这里添加其他实验的参数
% };

% parameters = {
%     {10^(-1.001242541230301), 1.441514803767596e+04, 0.904969462898074, 10^(-9.901777612793937), 10^(-3)},
%     {10^(-1.001242541230301), 1.441514803767596e+04, 0.904969462898074, 10^(-9.901777612793937), 10^(-4)},
%     {10^(-1.001242541230301), 1.441514803767596e+04, 0.904969462898074, 10^(-9.901777612793937), 10^(-5)},
%     {10^(-1.001242541230301), 1.441514803767596e+04, 0.904969462898074, 10^(-9.901777612793937), 10^(-6)},
%     % 在这里添加其他实验的参数
% };

for i = 1:length(experiment_names)
    exp_name = experiment_names{i};
    
    % 设置input的参数
    input.lambda = parameters{i}{1};
    input.sigma_p = parameters{i}{2};
    input.omega_0_variation = parameters{i}{3};
    input.Q_value = parameters{i}{4};
    input.R_value = parameters{i}{5};
    
    modesel = 23;
    input.modesel = modesel;
    
    % 运行实验并保存结果
    results_experiment.(exp_name) = run_experiment(input, 'showtext', false, 'showplot', false);
end




% results_experiment = run_experiment(input, 'showtext', false, 'showplot', false);

% result_Main = results_experiment.result_Main;
% result_Damping = results_experiment.result_Damping;
% result_wind = results_experiment.result_wind;
% disp_dir = results_experiment.disp_dir;
% ex = results_experiment.ex;
% epsx = results_experiment.epsx;
% yn = results_experiment.yn;
% h_hat = results_experiment.h_hat;
% nmodes = results_experiment.nmodes;
% amp_cell = results_experiment.result_Damping.amp_cell;
% zeta_all_cell = results_experiment.result_Damping.zeta_all_cell;
% top_freqs = results_experiment.result_Damping.top_freqs;
% t_cycle_mean_cell = results_experiment.result_Damping.t_cycle_mean_cell;
% yn_reconstruct = results_experiment.result_Main.yn_reconstruct;
% node_shape = results_experiment.node_shape;
% t = results_experiment.t;
% u = results_experiment.u;
% udot = results_experiment.udot;
% u2dot = results_experiment.u2dot;
% x_filt_original = results_experiment.x_filt_original;
% t_temp = results_experiment.t_temp;
% p_filt_m = results_experiment.p_filt_m;
% MM = results_experiment.MM;
% CC = results_experiment.CC;
% KK = results_experiment.KK;
% acc_matrix_seq = results_experiment.acc_matrix_seq;
% acc_node_list = results_experiment.acc_node_list;
% acc_node = results_experiment.acc_node;
% node_loc = results_experiment.node_loc;
% nodeondeck = results_experiment.nodeondeck;
% eig_vec = results_experiment.eig_vec;
% Mapping_data = results_experiment.Mapping_data;
% loc_acc = results_experiment.loc_acc;
% mode_deck = results_experiment.mode_deck;
% mode_deck_re = results_experiment.mode_deck_re;

if fig_bool
[figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
for t1 = 1:length(experiment_names)
    exp_name = experiment_names{t1};
    amp_cell = results_experiment.(exp_name).result_Damping.amp_cell;
    zeta_all_cell = results_experiment.(exp_name).result_Damping.zeta_all_cell;
    t_cycle_mean_cell = results_experiment.(exp_name).result_Damping.t_cycle_mean_cell;
    mode_deck = results_experiment.(exp_name).mode_deck;
    top_freqs = results_experiment.(exp_name).result_Damping.top_freqs;
    zetam = results_experiment.(exp_name).zetam;
    k1 = 1;
    k2 = 1;

    % 假设 t_cycle_mean_cell{k1}{k2} 是一个包含 datetime 对象的数组
    datetimeArray = t_cycle_mean_cell{k1}{k2};

    % 提取第一个 datetime 对象作为参考点
    referenceDatetime = datetimeArray(1);

    % 计算每个 datetime 对象相对于参考点的秒数
    secondsFromReference = seconds(datetimeArray - referenceDatetime);

    % 现在，secondsFromReference 包含相对于第一个时间戳的秒数

    scatter(amp_cell{k1}{k2} * max(mode_deck(:, k1)), zeta_all_cell{k1}{k2}, []);
    

    hold on
   
    % % scatter(ex,epsx,'green')
    str = "Using Force: Mode : %d, Frequency : %.2f Hz";
    title(sprintf(str, modesel(k1), top_freqs{k1}(k2)));
    xlim([0.05, 0.12])
    ylim([-0.5, 0.5] / 100)
    xlabel("Amplitude(m)")
    ylabel("Damping ratio")
    
end
legend(experiment_names)
plot([0, 0.15], [-0.003, -0.003])

[figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
for t1 = 1:length(experiment_names)
    exp_name = experiment_names{t1};
    amp_cell = results_experiment.(exp_name).result_Damping.amp_cell;
    zeta_all_cell = results_experiment.(exp_name).result_Damping.zeta_all_cell;
    t_cycle_mean_cell = results_experiment.(exp_name).result_Damping.t_cycle_mean_cell;
    mode_deck = results_experiment.(exp_name).mode_deck;
    top_freqs = results_experiment.(exp_name).result_Damping.top_freqs;
    zetam = results_experiment.(exp_name).zetam;
    k1 = 1;
    k2 = 1;

    % 假设 t_cycle_mean_cell{k1}{k2} 是一个包含 datetime 对象的数组
    datetimeArray = t_cycle_mean_cell{k1}{k2};

    % 提取第一个 datetime 对象作为参考点
    referenceDatetime = datetimeArray(1);

    % 计算每个 datetime 对象相对于参考点的秒数
    secondsFromReference = seconds(datetimeArray - referenceDatetime);

    % 现在，secondsFromReference 包含相对于第一个时间戳的秒数

    scatter(amp_cell{k1}{k2} * max(mode_deck(:, k1)), zetam-0.3/100, []);
    

    hold on
   
    % % scatter(ex,epsx,'green')
    str = "Using Acceleration: Mode : %d, Frequency : %.2f Hz";
    title(sprintf(str, modesel(k1), top_freqs{k1}(k2)));
    xlim([0.05, 0.12])
    ylim([-0.5, 0.5] / 100)
    xlabel("Amplitude(m)")
    ylabel("Damping ratio")
    
end
legend(experiment_names)
plot([0, 0.15], [-0.003, -0.003])


[figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
for t1 = 1:length(experiment_names)
    exp_name = experiment_names{t1};
    h_hat = results_experiment.(exp_name).h_hat;
    t = results_experiment.(exp_name).t;
    plot(t, h_hat(1, :));
    hold on
    xlabel("Time(s)")
    ylabel("Acceleration(m/s^2)")
    title("Acceleration comparison")
end
legend(experiment_names)

[figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
for t1 = 1:length(experiment_names)
    exp_name = experiment_names{t1};
    h_hat = results_experiment.(exp_name).h_hat;
    t = results_experiment.(exp_name).t;
    plot(t, h_hat(3, :));
    hold on
    xlabel("Time(s)")
    ylabel("Velocity(m/s)")
    title("Velocity comparison")
end
legend(experiment_names)

[figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
for t1 = 1:length(experiment_names)
    exp_name = experiment_names{t1};
    h_hat = results_experiment.(exp_name).h_hat;
    t = results_experiment.(exp_name).t;
    plot(t, h_hat(5, :));
    hold on
    xlabel("Time(s)")
    ylabel("Displacement(m)")
    title("Displacement comparison")
end
legend(experiment_names)

[figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
for t1 = 1:length(experiment_names)
    exp_name = experiment_names{t1};
    p_filt_m = results_experiment.(exp_name).p_filt_m;
    t = results_experiment.(exp_name).t;
    plot(t, p_filt_m(1, :));
    hold on
    xlabel("Time(s)")
    ylabel("Modal Force")
    title("Modal Force comparison")
end
legend(experiment_names)

[figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
for t1 = 1:length(experiment_names)
    exp_name = experiment_names{t1};
    yn_reconstruct = results_experiment.(exp_name).yn_reconstruct;
    t = results_experiment.(exp_name).t;
    plot(t, yn_reconstruct(1, :));
    hold on
    xlabel("Time(s)")
    ylabel("Acceleration(m/s^2)")
    title("Recalculated acceleration")
end
legend(experiment_names)


[figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
for t1 = 1:length(experiment_names)
    exp_name = experiment_names{t1};
    h_hat = results_experiment.(exp_name).h_hat;
    udot = results_experiment.(exp_name).udot;
    t = results_experiment.(exp_name).t;
    t_seconds = seconds(t - t(1));
    p_filt_m = results_experiment.(exp_name).p_filt_m;
    v= h_hat(3, :);
    v= udot;
    f = p_filt_m(1, :);

    power = f.*v;
    total_work_done = trapz(t_seconds, power);
    disp(['Total work done by the force = ', num2str(total_work_done), ' experiment name: ', exp_name])

    cumulated_work_done = cumtrapz(t_seconds, power);
    plot(t, cumulated_work_done);
    hold on
    xlabel("Time(s)")
    ylabel("Cumulated work done")
    title("Cumulated work done comparison")
    grid on;
end

% [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
% for t1 = 1:length(experiment_names)
%     exp_name = experiment_names{t1};
%     workcell = results_experiment.(exp_name).result_Damping.work_cell;
%     workcell_temp = workcell{1}{1}
% 
% end

[figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
for t1 = 1:length(experiment_names)
    exp_name = experiment_names{t1};
    h_hat = results_experiment.(exp_name).h_hat;
    t = results_experiment.(exp_name).t;
    t_seconds = seconds(t - t(1));
    p_filt_m = results_experiment.(exp_name).p_filt_m;
    v= h_hat(3, :);
    v= udot;
    f = p_filt_m(1, :);

    power = f.*v;
    
    plot(t, power);
    hold on
    xlabel("Time(s)")
    ylabel("Power")
    title("Power comparison")
    grid on;
end
end
% TODO: 为什么累计的功都是上升的，因为气动力在这里基本都做正功，阻尼力做负功
save results_experiment
% 
% offset = [0,0,-4.5/1000,1.5/1000];
% [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
% for t1 = 1:length(experiment_names)
%     exp_name = experiment_names{t1};
%     amp_cell = results_experiment.(exp_name).result_Damping.amp_cell;
%     zeta_all_cell = results_experiment.(exp_name).result_Damping.zeta_all_cell;
%     t_cycle_mean_cell = results_experiment.(exp_name).result_Damping.t_cycle_mean_cell;
%     mode_deck = results_experiment.(exp_name).mode_deck;
%     top_freqs = results_experiment.(exp_name).result_Damping.top_freqs;
%     k1 = 1;
%     k2 = 1;
% 
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
%     scatter(amp_cell{k1}{k2} * max(mode_deck(:, k1)), zeta_all_cell{k1}{k2}+offset(t1), []);
% 
% 
%     hold on
% 
%     % scatter(ex,epsx,'green')
%     str = "Mode : %d, Frequency : %.2f Hz";
%     title(sprintf(str, modesel(k1), top_freqs{k1}(k2)));
%     xlim([0.05, 0.12])
%     ylim([-0.5, 0.5] / 100)
%     xlabel("Amplitude(m)")
%     ylabel("Damping ratio")
% 
% end
% legend(experiment_names)
% plot([0, 0.15], [-0.003, -0.003])
