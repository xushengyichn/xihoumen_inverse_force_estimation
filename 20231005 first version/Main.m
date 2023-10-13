clc; clear; close all;
addpath(genpath("F:\git\ssm_tools\"))
addpath(genpath("F:\git\Function_shengyi_package\"))
addpath(genpath("F:\git\xihoumen_inverse_force_estimation\FEM_model\"))

subStreamNumberDefault = 2132;

run("InitScript.m")

%% 0 绘图参数
fig_bool = ON;
num_figs_in_row = 12; %每一行显示几个图
figPos = figPosSmall; %图的大小，参数基于InitScript.m中的设置
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
[result] = viv2013(n, OFF);
input.start_time = result.startDate;
input.end_time = result.endDate;
input.acc_dir = "F:\test\result";
% input.acc_dir = "Z:\Drive\Backup\SHENGYI_HP\F\test\result";

input.lambda = 10 ^ (-1);
input.sigma_p = 10000;
input.omega_0_variation =1;
input.Q_value =10 ^ (-8);
input.R_value = 10 ^ (-6);

input.lambda = 10 ^ (-4.934808796013671);
input.sigma_p = 6.895548550856822e+03;
input.omega_0_variation =1.097383030422062;
input.Q_value =10 ^ (-9.633948257379021);
input.R_value = 10 ^ (-2.415076745081128);

[result_Main] = KalmanMain(input,'showtext', false,'showfigure',false);
logL = result_Main.logL;



%% optimization logL to get the maximum with changing lambda sigma_p omega_0_variation Q_value R_value
% 在调用 ga 函数之前，您可以这样设置 external_params：
external_params.figPos = figPos;
external_params.ON = ON;
external_params.OFF = OFF;
% 定义参数的范围
lb = [-5, 10, 0.9, -10, -10]; % 这里的值是假设的，请根据您的情况进行修改
ub = [-1, 1e5, 1.1, -1, -1]; % 这里的值也是假设的

% 定义整数和连续变量
IntCon = []; % 如果没有整数变量，否则提供整数变量的索引

options = optimoptions('ga', 'MaxGenerations', 100, 'Display', 'iter', 'UseParallel', true);
[x, fval] = ga(@(params) fitnessFunction(params, external_params), 5, [], [], [], [], lb, ub, [], IntCon, options);
% 保存结果
save('optimization_results.mat', 'x', 'fval');

function logL = fitnessFunction(params,external_params)
    input.lambda = 10 ^ params(1);
    input.sigma_p = params(2);
    input.omega_0_variation = params(3);
    input.Q_value = 10 ^ params(4);
    input.R_value = 10 ^ params(5);

     figPos = external_params.figPos;
    ON = external_params.ON;
    OFF = external_params.OFF;

    input.num_figs_in_row = 12; %每一行显示几个图
    input.figPos = figPos; %图的大小，参数基于InitScript.m中的设置
    %设置图片间隔
    input.ON =ON;
    input.OFF =OFF;
    input.gap_between_images = [0, 0];
    input.figureIdx = 0;
    n = 4;
    [result] = viv2013(n, OFF);
    input.start_time = result.startDate;
    input.end_time = result.endDate;
    % input.acc_dir = "F:\test\result";
    input.acc_dir = "Z:\Drive\Backup\SHENGYI_HP\F\test\result";

    result_Main = KalmanMain(input, 'showtext', false, 'showfigure', false);
    logL = -result_Main.logL; % 因为 ga 试图最小化函数，所以取负数
end



