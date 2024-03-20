%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn xushengyichn@outlook.com
%Date: 2023-09-19 02:17:00
%LastEditors: ShengyiXu xushengyichn@outlook.com
%LastEditTime: 2024-01-20 13:19:03
%FilePath: \manuscriptf:\git\xihoumen_inverse_force_estimation\20231203 third edition\readdata.m
%Description:
%
%Copyright (c) 2023 by ${git_name_email}, All Rights Reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;close all
run('CommonCommand.m');

%% vibration (This part is replaced by the following code)
if 0
    startDate = [datetime(2013,1,14,0,0,0); datetime(2013,1,16,3,0,0); datetime(2013,1,22,2,0,0); ...
        datetime(2013,2,6,0,0,0); datetime(2013,3,6,18,0,0); datetime(2013,4,3,15,0,0); ...
        datetime(2013,4,8,15,0,0); datetime(2013,4,10,19,0,0); datetime(2013,5,19,11,0,0); ...
        datetime(2013,5,21,18,0,0); datetime(2013,6,10,4,0,0); datetime(2013,7,9,17,0,0); ...
        datetime(2013,7,11,18,0,0); datetime(2013,7,16,15,0,0); datetime(2013,7,21,16,0,0); ...
        datetime(2013,8,1,12,0,0); datetime(2013,8,2,22,0,0); datetime(2013,8,10,22,0,0); ...
        datetime(2013,8,16,18,0,0); datetime(2013,8,24,15,0,0); datetime(2013,8,28,14,0,0); ...
        datetime(2013,8,29,14,0,0)];
    
    endDate = [datetime(2013,1,14,2,0,0); datetime(2013,1,16,5,0,0); datetime(2013,1,22,4,0,0); ...
        datetime(2013,2,6,3,0,0); datetime(2013,3,6,20,0,0); datetime(2013,4,3,17,0,0); ...
        datetime(2013,4,8,16,0,0); datetime(2013,4,10,22,0,0); datetime(2013,5,19,13,0,0); ...
        datetime(2013,5,21,20,0,0); datetime(2013,6,10,7,0,0); datetime(2013,7,9,18,0,0); ...
        datetime(2013,7,11,19,0,0); datetime(2013,7,16,17,0,0); datetime(2013,7,21,20,0,0); ...
        datetime(2013,8,1,16,0,0); datetime(2013,8,3,0,0,0); datetime(2013,8,10,23,0,0); ...
        datetime(2013,8,16,19,0,0); datetime(2013,8,24,21,0,0); datetime(2013,8,28,17,0,0); ...
        datetime(2013,8,29,15,0,0)];
    
    % 创建案例编号数组
    caseNumber = (1:22)';
    
    % 使用table函数创建表格
    vivTable = table(caseNumber, startDate, endDate);
    
    % 显示表格
    disp(vivTable);
end
%% Read viv Data
% 创建导入选项对象
opts = detectImportOptions('vivData.csv');

% 设置日期时间格式
% 假设日期时间格式为 'MM/dd/yyyy HH:mm'，请根据您的实际情况进行调整
opts = setvaropts(opts, 'startDate', 'InputFormat', 'MM/dd/yyyy HH:mm');
opts = setvaropts(opts, 'endDate', 'InputFormat', 'MM/dd/yyyy HH:mm');
opts = setvaropts(opts, 'startDate_update', 'InputFormat', 'MM/dd/yyyy HH:mm');
opts = setvaropts(opts, 'endDate_update', 'InputFormat', 'MM/dd/yyyy HH:mm');

vivTable = readtable('vivData.csv',opts);

%% Readdata

% for k1 = 2
for k1 = 12
    acc_dir = input_data.acc_dir;
    wind_dir = input_data.wind_dir;
    start_time = vivTable.startDate(k1);
    end_time = vivTable.endDate(k1);
    % start_time_update = vivTable.update_startDate(k1);
    % end_time_update = vivTable.endDate_update(k1);
    start_time.Format = 'dd-MMM-yyyy HH:mm:ss';
    end_time.Format = 'dd-MMM-yyyy HH:mm:ss';
    disp(start_time);
    disp(end_time);
    
end
[Acc_Data] = read_acceleration_data(start_time, end_time, acc_dir);
[Wind_Data] = read_wind_data(start_time, end_time, wind_dir);
Acc_Data = Acc_Data.mergedData;
t_acc = Acc_Data.Time;
AC2_1 = Acc_Data.AC2_1;
AC2_3 = Acc_Data.AC2_3;
AC3_1 = Acc_Data.AC3_1;
AC3_3 = Acc_Data.AC3_3;
AC4_1 = Acc_Data.AC4_1;
AC4_3 = Acc_Data.AC4_3; 

[f2_1, magnitude2_1] = fft_transform(50,AC2_1);
[f2_3, magnitude2_3] = fft_transform(50,AC2_3);
[f3_1, magnitude3_1] = fft_transform(50,AC3_1);
[f3_3, magnitude3_3] = fft_transform(50,AC3_3);
[f4_1, magnitude4_1] = fft_transform(50,AC4_1);
[f4_3, magnitude4_3] = fft_transform(50,AC4_3);

t_wind = Wind_Data.resultsTable_UA1.Time_Start;
u1 = Wind_Data.resultsTable_UA1.U;
u2 = Wind_Data.resultsTable_UA2.U;
u3 = Wind_Data.resultsTable_UA3.U;
u4 = Wind_Data.resultsTable_UA4.U;
u5 = Wind_Data.resultsTable_UA5.U;
u6 = Wind_Data.resultsTable_UA6.U;

%% readupdate data
start_time_update = datetime(2013,7,6,2,20,0);
end_time_update = datetime(2013,7,6,3,45,0);
start_time_update.Format = 'MM/dd/yyyy HH:mm';
end_time_update.Format = 'MM/dd/yyyy HH:mm';
str = string(start_time_update)+","+string(end_time_update);
disp(str)
[Wind_Data_update] = read_wind_data(start_time_update, end_time_update, wind_dir);
[Acc_Data_update] = read_acceleration_data(start_time_update, end_time_update, acc_dir);
t_wind_update = Wind_Data_update.resultsTable_UA1.Time_Start;
u1_update = Wind_Data_update.resultsTable_UA1.U;
u2_update = Wind_Data_update.resultsTable_UA2.U;
u3_update = Wind_Data_update.resultsTable_UA3.U;
u4_update = Wind_Data_update.resultsTable_UA4.U;
u5_update = Wind_Data_update.resultsTable_UA5.U;
u6_update = Wind_Data_update.resultsTable_UA6.U;

t_acc_update = Acc_Data_update.mergedData.Time;
AC2_1_update = Acc_Data_update.mergedData.AC2_1;
AC2_3_update = Acc_Data_update.mergedData.AC2_3;
AC3_1_update = Acc_Data_update.mergedData.AC3_1;
AC3_3_update = Acc_Data_update.mergedData.AC3_3;
AC4_1_update = Acc_Data_update.mergedData.AC4_1;
AC4_3_update = Acc_Data_update.mergedData.AC4_3;


%% Plot
fig_bool = true;
fig_bool_update = true;
if fig_bool
    % 定义总子图数量
    total_plots = 18; % 或任何你需要的子图数量
    current_plot = 1;
    num_figs_in_row = [];
    figWidthFactor = 1.5;
    %     figPosition = [1080*2.5,100];
    figPosition = [100, 100];
    newfigure = true;
    firstfigure = true;
    holdon = false;

    create_subplot(@plot, total_plots, current_plot, {t_acc,AC2_1}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);
    firstfigure = false;
    current_plot = current_plot + 1;

    create_subplot(@plot, total_plots, current_plot, {t_acc,AC2_3}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);
    firstfigure = false;
    current_plot = current_plot + 1;

    create_subplot(@plot, total_plots, current_plot, {t_acc,AC3_1}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);
    firstfigure = false;
    current_plot = current_plot + 1;

    create_subplot(@plot, total_plots, current_plot, {t_acc,AC3_3}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);
    firstfigure = false;
    current_plot = current_plot + 1;

    create_subplot(@plot, total_plots, current_plot, {t_acc,AC4_1}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);
    firstfigure = false;
    current_plot = current_plot + 1;

    create_subplot(@plot, total_plots, current_plot, {t_acc,AC4_3}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);
    firstfigure = false;
    current_plot = current_plot + 1;

    create_subplot(@plot, total_plots, current_plot, {f2_1,magnitude2_1}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);
    xlim([0,0.5]);
    firstfigure = false;
    current_plot = current_plot + 1;

    create_subplot(@plot, total_plots, current_plot, {f2_3,magnitude2_3}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);
    xlim([0,0.5]);
    firstfigure = false;
    current_plot = current_plot + 1;

    create_subplot(@plot, total_plots, current_plot, {f3_1,magnitude3_1}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);
    xlim([0,0.5]);
    firstfigure = false;
    current_plot = current_plot + 1;

    create_subplot(@plot, total_plots, current_plot, {f3_3,magnitude3_3}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);
    xlim([0,0.5]);
    firstfigure = false;
    current_plot = current_plot + 1;

    create_subplot(@plot, total_plots, current_plot, {f4_1,magnitude4_1}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);
    xlim([0,0.5]);
    firstfigure = false;
    current_plot = current_plot + 1;

    create_subplot(@plot, total_plots, current_plot, {f4_3,magnitude4_3}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);
    xlim([0,0.5]);
    firstfigure = false;
    current_plot = current_plot + 1;

    create_subplot(@plot, total_plots, current_plot, {t_wind,u1}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', true, 'firstfigure', firstfigure, 'holdon', holdon);
    firstfigure = false;
    current_plot = current_plot + 1;

    create_subplot(@plot, total_plots, current_plot, {t_wind,u2}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', true, 'firstfigure', firstfigure, 'holdon', holdon);
    firstfigure = false;
    current_plot = current_plot + 1;

    create_subplot(@plot, total_plots, current_plot, {t_wind,u3}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', true, 'firstfigure', firstfigure, 'holdon', holdon);
    firstfigure = false;
    current_plot = current_plot + 1;

    create_subplot(@plot, total_plots, current_plot, {t_wind,u4}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', true, 'firstfigure', firstfigure, 'holdon', holdon);
    firstfigure = false;
    current_plot = current_plot + 1;

    create_subplot(@plot, total_plots, current_plot, {t_wind,u5}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', true, 'firstfigure', firstfigure, 'holdon', holdon);
    firstfigure = false;
    current_plot = current_plot + 1;

    create_subplot(@plot, total_plots, current_plot, {t_wind,u6}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', true, 'firstfigure', firstfigure, 'holdon', holdon);
    firstfigure = false;
    current_plot = current_plot + 1;


    % create_subplot(@plot, total_plots, current_plot, {t_wind_update,u1_update}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', true, 'firstfigure', firstfigure, 'holdon', holdon);
    % firstfigure = false;
    % current_plot = current_plot + 1;

    % create_subplot(@plot, total_plots, current_plot, {t_wind_update,u2_update}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', true, 'firstfigure', firstfigure, 'holdon', holdon);
    % firstfigure = false;
    % current_plot = current_plot + 1;

    % create_subplot(@plot, total_plots, current_plot, {t_wind_update,u3_update}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', true, 'firstfigure', firstfigure, 'holdon', holdon);
    % firstfigure = false;
    % current_plot = current_plot + 1;

    % create_subplot(@plot, total_plots, current_plot, {t_wind_update,u4_update}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', true, 'firstfigure', firstfigure, 'holdon', holdon);
    % firstfigure = false;
    % current_plot = current_plot + 1;

    % create_subplot(@plot, total_plots, current_plot, {t_wind_update,u5_update}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', true, 'firstfigure', firstfigure, 'holdon', holdon);
    % firstfigure = false;
    % current_plot = current_plot + 1;

    % create_subplot(@plot, total_plots, current_plot, {t_wind_update,u6_update}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', true, 'firstfigure', firstfigure, 'holdon', holdon);
    % firstfigure = false;
    % current_plot = current_plot + 1;


end


if fig_bool_update
    total_plots = 12; % 或任何你需要的子图数量
    current_plot = 1;
    num_figs_in_row = [];
    figWidthFactor = 1.5;
    %     figPosition = [1080*2.5,100];
    figPosition = [100, 100];
    newfigure = true;
    firstfigure = true;
    holdon = false;

    create_subplot(@plot, total_plots, current_plot, {t_acc_update,AC2_1_update}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);
    firstfigure = false;
    current_plot = current_plot + 1;

    create_subplot(@plot, total_plots, current_plot, {t_acc_update,AC2_3_update}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);
    firstfigure = false;
    current_plot = current_plot + 1;

    create_subplot(@plot, total_plots, current_plot, {t_acc_update,AC3_1_update}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);
    firstfigure = false;
    current_plot = current_plot + 1;

    create_subplot(@plot, total_plots, current_plot, {t_acc_update,AC3_3_update}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);
    firstfigure = false;
    current_plot = current_plot + 1;

    create_subplot(@plot, total_plots, current_plot, {t_acc_update,AC4_1_update}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);
    firstfigure = false;
    current_plot = current_plot + 1;

    create_subplot(@plot, total_plots, current_plot, {t_acc_update,AC4_3_update}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);
    firstfigure = false;
    current_plot = current_plot + 1;

    create_subplot(@plot, total_plots, current_plot, {t_wind_update,u1_update}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', true, 'firstfigure', firstfigure, 'holdon', holdon);
    firstfigure = false;
    current_plot = current_plot + 1;

    create_subplot(@plot, total_plots, current_plot, {t_wind_update,u2_update}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', true, 'firstfigure', firstfigure, 'holdon', holdon);
    firstfigure = false;
    current_plot = current_plot + 1;

    create_subplot(@plot, total_plots, current_plot, {t_wind_update,u3_update}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', true, 'firstfigure', firstfigure, 'holdon', holdon);
    firstfigure = false;
    current_plot = current_plot + 1;

    create_subplot(@plot, total_plots, current_plot, {t_wind_update,u4_update}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', true, 'firstfigure', firstfigure, 'holdon', holdon);
    firstfigure = false;
    current_plot = current_plot + 1;

    create_subplot(@plot, total_plots, current_plot, {t_wind_update,u5_update}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', true, 'firstfigure', firstfigure, 'holdon', holdon);
    firstfigure = false;
    current_plot = current_plot + 1;

    create_subplot(@plot, total_plots, current_plot, {t_wind_update,u6_update}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', true, 'firstfigure', firstfigure, 'holdon', holdon);
    firstfigure = false;
    current_plot = current_plot + 1;


end


