%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn xushengyichn@outlook.com
%Date: 2023-10-20 12:45:02
%LastEditors: ShengyiXu xushengyichn@outlook.com
%LastEditTime: 2023-10-23 23:30:40
%FilePath: \Exercises-for-Techniques-for-estimation-in-dynamics-systemsf:\git\xihoumen_inverse_force_estimation\20231005 first version\tune_hyperparameter.m
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
addpath(genpath("C:\Users\xushengyi\Documents\Github\ssm_tools"))
addpath(genpath("C:\Users\xushengyi\Documents\Github\Function_shengyi_package\"))
addpath(genpath("C:\Users\xushengyi\Documents\Github\xihoumen_inverse_force_estimation\FEM_model\"))

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
% input.acc_dir = "F:\test\result";
% input.wind_dir = "F:\test\result_wind_10min";
input.acc_dir = "Z:\Drive\Backup\SHENGYI_HP\F\test\result";
input.wind_dir = "Z:\Drive\Backup\SHENGYI_HP\F\test\result_wind_10min";


% % 定义每个实验的参数
% parameters = {
%     {10^(-1.001242541230301), 1.441514803767596e+04, 0.904969462898074, 10^(-9.901777612793937), 10^(-3.866588296864785)},
%     {10^(-2.001242541230301), 1.441514803767596e+04, 0.904969462898074, 10^(-9.901777612793937), 10^(-3.866588296864785)},
%     {10^(-3.001242541230301), 1.441514803767596e+04, 0.904969462898074, 10^(-9.901777612793937), 10^(-3.866588296864785)},
%     {10^(-4.001242541230301), 1.441514803767596e+04, 0.904969462898074, 10^(-9.901777612793937), 10^(-3.866588296864785)},
%     % {10 ^ (-1.001242541230301), 9.063096667830060e+04, 0.902647415472734, 10 ^ (-1.016211706804576), 10 ^ (-1.003148221125874)},
%     % {10 ^ (-4.934808796013671), 6.895548550856822e+03, 1.097383030422062, 10 ^ (-9.633948257379021), 10 ^ (-2.415076745081128)},
%     % {10 ^ (-4.993486819657864), 1.441514803767596e+04, 0.904969462898074, 10^(-9.901777612793937), 10^(-3.866588296864785)},
%     % 在这里添加其他实验的参数
% };

n1=10;
n2=30;
lambda_list = logspace(-5,-1,n1);
sigma_p_list = linspace(1e2,1e5,n2);
% sigma_p_list = logspace(2,5,n2);

[X, Y] = meshgrid(lambda_list, sigma_p_list);
combinations = [reshape(X, [], 1), reshape(Y, [], 1)];

lambda_vals = combinations(:, 1);
sigma_p_vals = combinations(:, 2);
numIterations = size(combinations,1);

if isempty(gcp('nocreate'))
    parpool();
end

b = ProgressBar(numIterations, ...
    'IsParallel', true, ...
    'WorkerDirectory', pwd(), ...
    'Title', 'Parallel 2' ...
    );
b.setup([], [], []);


parfor i = 1:numIterations
% for i = 1:numIterations
    % 创建input的局部副本
    local_input = input;
    
    % 设置local_input的参数
    local_input.lambda = lambda_vals(i);
    local_input.sigma_p = sigma_p_vals(i);
    local_input.omega_0_variation = 0.904969462898074;
    local_input.Q_value = 10 ^ (-9.901777612793937);
    local_input.R_value = 10 ^ (-3.866588296864785);
    
    modesel = 23;
    local_input.modesel = modesel;
    
    % 运行实验并保存结果
    % tic
    results_experiment = run_experiment(local_input, 'showtext', false, 'showplot', false,'caldamp',false);
    % toc
    % logL(i) = results_experiment.result_Main.logL;
    logSk(i) = results_experiment.result_Main.logSk;
    % logek(i) = results_experiment.result_Main.logek;
    node_shape = results_experiment.node_shape;
    h_hat = results_experiment.h_hat;
    u2dot = results_experiment.u2dot;
    invSk = results_experiment.invSk;
    u2dot_real = u2dot.*ones(2,1)*node_shape;
    logek_new = 0;
    for k2 = 1:length(u2dot)
        ek_temp = u2dot_real(:,k2)- h_hat(1:2,k2);
        logek_temp = -0.5*ek_temp'*invSk*ek_temp;
        logek_new =logek_new + logek_temp;
    end
    logek(i) = logek_new;
    logL(i) = logSk(i) + logek(i);
    % USE THIS FUNCTION AND NOT THE STEP() METHOD OF THE OBJECT!!!
    updateParallel([], pwd);
end



b.release();

Z_1 = reshape(logL, n2, n1);
Z_2 = reshape(logSk, n2, n1);
Z_3 = reshape(logek, n2, n1);

%% plot
if fig_bool == ON

    [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);

    
    contourf(X, Y, Z_1);  % 绘制等高线图
    set(gca, 'XScale', 'log');
    xlabel('lambdas');
    ylabel('sigma_ps');
    colorbar;  % 添加颜色栏
    title('logL');

[figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
    contourf(X, Y, Z_2);  % 绘制等高线图
    set(gca, 'XScale', 'log');
    xlabel('lambdas');
    ylabel('sigma_ps');
    colorbar;  % 添加颜色栏
    title('logSk');


    [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
    contourf(X, Y, Z_3);  % 绘制等高线图
    set(gca, 'XScale', 'log');
    xlabel('lambdas');
    ylabel('sigma_ps');
    colorbar;  % 添加颜色栏
    title('logek');

    % 画 logL 的三维图
[figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
surf(X, Y, Z_1);  % 创建三维表面图
set(gca, 'XScale', 'log');
xlabel('lambdas');
ylabel('sigma_ps');
zlabel('logL');
colorbar;  % 添加颜色栏
title('3D plot of logL');

% 画 logSk 的三维图
[figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
surf(X, Y, Z_2);  % 创建三维表面图
set(gca, 'XScale', 'log');
xlabel('lambdas');
ylabel('sigma_ps');
zlabel('logSk');
colorbar;  % 添加颜色栏
title('3D plot of logSk');

% 画 logek 的三维图
[figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
surf(X, Y, Z_3);  % 创建三维表面图
set(gca, 'XScale', 'log');
xlabel('lambdas');
ylabel('sigma_ps');
zlabel('logek');
colorbar;  % 添加颜色栏
title('3D plot of logek');


end