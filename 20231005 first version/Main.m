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
input.start_time = result.startDate-hours(0.5);
input.end_time = result.endDate+hours(0.5);
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

% modesel= [2,3,5,6,7,9,15,21,23,29,33,39,44,45];
modesel= 23;
input.modesel= modesel;

[result_Main] = KalmanMain(input,'showtext', true,'showplot',false);
logL = result_Main.logL;
mode_deck = result_Main.mode_deck;



if 0
    %% optimization logL to get the maximum with changing lambda sigma_p omega_0_variation Q_value R_value
    % 在调用 ga 函数之前，您可以这样设置 external_params：
    external_params.modelsel = [2,3,5,6,7,9,15,21,23,29,33,39,44,45];
    % 定义参数的范围
    lb = [-5, 10, 0.9, -10, -10]; % 这里的值是假设的，请根据您的情况进行修改
    ub = [-1, 1e5, 1.1, -1, -1]; % 这里的值也是假设的
    
    % 定义整数和连续变量
    IntCon = []; % 如果没有整数变量，否则提供整数变量的索引
    
    options = optimoptions('ga', 'MaxGenerations', 100, 'Display', 'iter', 'UseParallel', true);
    [x, fval] = ga(@(params) fitnessFunction(params, external_params), 5, [], [], [], [], lb, ub, [], IntCon, options);
    % 保存结果
    save('optimization_results.mat', 'x', 'fval');

end

%% Caldamping ratio
input = result_Main;
input.ncycle = 5;
tic
[result_Damping]=Cal_aero_damping_ratio(input,'showplot',false,'filterstyle','fft');
toc

%% read wind data

[result] = viv2013(n, OFF);
start_time = result.startDate;
end_time = result.endDate+hours(0.15);
acc_dir = "F:\test\result_wind_10min";
[result_wind] = read_wind_data(start_time,end_time,acc_dir);

%% plot 
t = result_Main.t;
yn = result_Main.yn;
h_hat = result_Main.h_hat;
nmodes = result_Main.nmodes;
amp_cell = result_Damping.amp_cell;
zeta_all_cell = result_Damping.zeta_all_cell;
top_freqs = result_Damping.top_freqs;
t_cycle_mean_cell = result_Damping.t_cycle_mean_cell;

[figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
plot(t,yn(1,:));
hold on
plot(t,h_hat(1,:));
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
title("Acceleration vs. Time")
legend("measure","filtered")



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
        scatter(amp_cell{k1}{k2}*max(mode_deck(:,k1)),zeta_all_cell{k1}{k2},[],secondsFromReference,'filled');
        % 设置 colormap
        colormap('jet')
        colorbar
        
        hold on
        plot([0,0.15],[-0.003,-0.003])

        str = "Mode : %d, Frequency : %.2f Hz";
        title(sprintf(str,modesel(k1),top_freqs{k1}(k2)));
        xlim([0,0.15])
        ylim([-5,5]/100)
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
        scatter(secondsFromReference,amp_cell{k1}{k2}*max(mode_deck(:,k1)),[],secondsFromReference,'filled');
        % 设置 colormap
        colormap('jet')
        colorbar
        
        hold on
        plot([0,0.15],[-0.003,-0.003])

        str = "Mode : %d, Frequency : %.2f Hz, Amp vs. Time";
        title(sprintf(str,modesel(k1),top_freqs{k1}(k2)));
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
        scatter(secondsFromReference,zeta_all_cell{k1}{k2},[],secondsFromReference,'filled');
        % 设置 colormap
        colormap('jet')
        colorbar
        
        hold on
        plot([0,0.15],[-0.003,-0.003])

        str = "Mode : %d, Frequency : %.2f Hz";
        title(sprintf(str,modesel(k1),top_freqs{k1}(k2)));
        % xlim([0,0.15])
        % ylim([-25,25]/100)
        xlabel("Time(s)")
        ylabel("Damping ratio")
        
    end
end

%% functions
function logL = fitnessFunction(params,external_params)
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
        [result] = viv2013(n, OFF);
        input.start_time = result.startDate;
        input.end_time = result.endDate;
        % input.acc_dir = "F:\test\result";
        input.acc_dir = "Z:\Drive\Backup\SHENGYI_HP\F\test\result";
    
        result_Main = KalmanMain(input, 'showtext', false, 'showfigure', false);
        logL = -result_Main.logL; % 因为 ga 试图最小化函数，所以取负数
    end