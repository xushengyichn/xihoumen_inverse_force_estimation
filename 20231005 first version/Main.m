clc; clear; close all;
addpath(genpath("F:\git\ssm_tools\"))
addpath(genpath("F:\git\Function_shengyi_package\"))
addpath(genpath("F:\git\xihoumen_inverse_force_estimation\FEM_model\"))
addpath(genpath("F:\git\HHT-Tutorial\"))

subStreamNumberDefault = 2132;

run("InitScript.m")

%% 0 绘图参数
fig_bool = ON;
num_figs_in_row = 3; %每一行显示几个图
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
startDate_global = result.startDate;
endDate_global= result.endDate;
input.start_time = startDate_global;
input.end_time = endDate_global;
input.acc_dir = "F:\test\result";
% input.acc_dir = "Z:\Drive\Backup\SHENGYI_HP\F\test\result";

% input.lambda = 10 ^ (-1);
% input.sigma_p = 10000;
% input.omega_0_variation =1;
% input.Q_value =10 ^ (-8);
% input.R_value = 10 ^ (-6);

% input.lambda = 10 ^ (-4.934808796013671);
% input.sigma_p = 6.895548550856822e+03;
% input.omega_0_variation =1.097383030422062;
% input.Q_value =10 ^ (-9.633948257379021);
% input.R_value = 10 ^ (-2.415076745081128);

input.lambda = 10 ^ (-4.993486819657864);
input.sigma_p = 1.441514803767596e+04;
input.omega_0_variation =0.904969462898074;
input.Q_value =10 ^ (-9.901777612793937);
input.R_value = 10 ^ (-3.866588296864785);

% modesel= [2,3,5,6,7,9,15,21,23,29,33,39,44,45];
modesel= 23;
input.modesel= modesel;

[result_Main] = KalmanMain(input,'showtext', true,'showplot',false,'filterstyle','fft','f_keep', 0.33*[0.9,1.1]);
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
    % external_params.modelsel = [2,3,5,6,7,9,15,21,23,29,33,39,44,45];
    external_params.modelsel = [23];
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
input.ncycle = 3;
tic
[result_Damping]=Cal_aero_damping_ratio(input,'showplot',false,'filterstyle','nofilter');
toc

%% read wind data

% [result] = viv2013(n, OFF);
start_time = startDate_global;
end_time = endDate_global;
acc_dir = "F:\test\result_wind_10min";
[result_wind] = read_wind_data(start_time,end_time,acc_dir);

%% compare with ee
yn = result_Main.yn;
fs = 50;
t = result_Main.t;
disp_dir=acc2dsip(yn(1,:),50);
[ex, frex] = ee(disp_dir.disp, 1 / fs); %经验包络法求瞬时频率和瞬时振幅 % empirical envelope method to find instantaneous frequency and instantaneous amplitude


omgx = frex*2*pi;
for k1 = 1:length(ex)-1
    epsx(k1)=log(ex(k1)/ex(k1+1))/omgx(k1)*fs;
end

epsx = [epsx epsx(end)];

% figure
% plot(t,epsx)
% grid

%% Recalcualte the modal displacement
nmodes = result_Main.nmodes;
Result = ImportMK(nmodes, 'KMatrix.matrix', 'MMatrix.matrix', 'nodeondeck.txt', 'KMatrix.mapping', 'nodegap.txt', 'modesel', [23]);
mode_deck = Result.mode_deck; mode_deck_re = Result.mode_deck_re; node_loc = Result.node_loc;nodeondeck = Result.nodeondeck;
eig_vec=Result.eig_vec;
Mapping_data = Result.Mapping;
loc_acc = [1403];
acc_node=FindNodewithLocation(loc_acc, node_loc, nodeondeck);
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
[u udot u2dot] = NewmarkInt(t_temp,MM,CC,KK,p_filt_m,1/2,1/4,0,0);
x_filt_original = result_Main.x_filt_original;

[figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
plot(t_temp,x_filt_original(1,:))
hold on
plot(t_temp,u)
legend("filtered","recalculate")
title("u")

[figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
plot(t_temp,x_filt_original(2,:))
hold on
plot(t_temp,udot)
legend("filtered","recalculate")
title("udot")

[figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
plot(t_temp,h_hat(1,:)/node_shape)
hold on
plot(t_temp,u2dot)
legend("filtered","recalculate")
title("u2dot")


%% plot 
t = result_Main.t;
yn = result_Main.yn;
h_hat = result_Main.h_hat;
nmodes = result_Main.nmodes;
amp_cell = result_Damping.amp_cell;
zeta_all_cell = result_Damping.zeta_all_cell;
top_freqs = result_Damping.top_freqs;
t_cycle_mean_cell = result_Damping.t_cycle_mean_cell;
yn_reconstruct =result_Main.yn_reconstruct;

[figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
scatter(result_wind.resultsTable_UA6.Time_Start,result_wind.resultsTable_UA6.U);
xlabel('Time (s)')
ylabel('Wind speed (m/s^2)')
title("Wind speed vs. Time")



[figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
plot(t,yn(1,:));
hold on
plot(t,h_hat(1,:));
plot(t,yn_reconstruct(1,:))
xlabel('Time (s)')
ylabel('Acceleration (m/s^2)')
title("Acceleration vs. Time")
legend("measure","filtered","reconstruct")


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
plot(t,disp_dir.disp);
hold on
plot(t,h_hat(3,:));

xlabel('Time (s)')
ylabel('Displacement (m)')
title("Displacement vs. Time")
legend("measure","filtered")


[figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
[p, fd, td] = pspectrum(h_hat(1,:), t, 'spectrogram', 'FrequencyResolution', 0.005);
instfreq(p, fd, td);
ylim([0.1,0.5])

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
        % scatter(ex,epsx,'green')
        str = "Mode : %d, Frequency : %.2f Hz";
        title(sprintf(str,modesel(k1),top_freqs{k1}(k2)));
        xlim([0.05,0.12])
        ylim([-0.5,0.5]/100)
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
        scatter(datetimeArray,amp_cell{k1}{k2}*max(mode_deck(:,k1)),[],secondsFromReference,'filled');
        % 设置 colormap
        colormap('jet')
        % colorbar
        
        % hold on
        % plot([0,0.15],[-0.003,-0.003])

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
        scatter(datetimeArray,zeta_all_cell{k1}{k2},[],secondsFromReference,'filled');
        % 设置 colormap
        colormap('jet')
        % colorbar
        
        % hold on
        % plot([datetimeArray(1),datetimeArray(end)],[-0.003,-0.003])

        str = "Mode : %d, Frequency : %.2f Hz";
        title(sprintf(str,modesel(k1),top_freqs{k1}(k2)));
        % xlim([0.05,0.12])
        ylim([-0.5,0.5]/100)
        xlabel("Time(s)")
        ylabel("Damping ratio")
        
    end
end



[figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);

data =disp_dir.disp;
minData = min(data);
maxData = max(data);
normalizedData = 2 * ((data - minData) / (maxData - minData)) - 1;
plot(t,normalizedData);


xlabel('Time (s)')
ylabel('Displacement (m)')
title("Displacement vs. Time")
legend("measure","filtered")
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
        scatter(datetimeArray,(normalizedData-mean(normalizedData(end/4:end*3/4)))*10,[],secondsFromReference,'filled');
        % 设置 colormap
        colormap('jet')
        % colorbar
        
        % hold on
        % plot([datetimeArray(1),datetimeArray(end)],[-0.003,-0.003])

        str = "Mode : %d, Frequency : %.2f Hz";
        title(sprintf(str,modesel(k1),top_freqs{k1}(k2)));
        % xlim([0.05,0.12])
        ylim([-1,1])
        xlabel("Time(s)")
        ylabel("Damping ratio")
        
    end
end

hold off
[figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
plot(t,frex)

[figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
plot(t,ex)
[figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
scatter(ex,epsx,'green')
ylim([-0.01,0.01])
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
    
        result_Main = KalmanMain(input,'showtext', false,'showplot',false,'filterstyle','fft','f_keep', 0.33*[0.9,1.1]);
        logL = -result_Main.logL; % 因为 ga 试图最小化函数，所以取负数
end

function [ amp,fre ] = ee( data,dt )
%EE Summary of this function goes here
%   Detailed explanation goes here
[normalizeddata,amp]=splinenormalizeep(data);
frecarrier=gradient(normalizeddata,dt);
[normalizedfrecarrier,ampfrecarrier]=splinenormalizeep(frecarrier);
fre=ampfrecarrier/2/pi;
end
