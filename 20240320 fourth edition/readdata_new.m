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
VIV_sels = [2;3;4;5;6;7;8;9;10;12;16;17;18;19];
% for k1 = 2
for k1 =  5
    acc_dir = input_data.acc_dir;
    wind_dir = input_data.wind_dir;
    dirName = input_data.wind_dir_all;
    start_time = vivTable.startDate(k1);
    end_time = vivTable.endDate(k1);
    start_time = datetime(2013,3,30,1,30,0);
    end_time = datetime(2013,3,30,2,0,0);
    % start_time_update = vivTable.update_startDate(k1);
    % end_time_update = vivTable.endDate_update(k1);
    start_time.Format = 'dd-MMM-yyyy HH:mm:ss';
    end_time.Format = 'dd-MMM-yyyy HH:mm:ss';
    disp(start_time);
    disp(end_time);
    
end
[Acc_Data] = read_acceleration_data(start_time, end_time, acc_dir);
[Wind_Data] = read_wind_data(start_time, end_time, wind_dir);
[Wind_Data_all] = read_wind_data_all(start_time,end_time,dirName);
Acc_Data = Acc_Data.mergedData;
t_acc = Acc_Data.Time;
AC2_1 = Acc_Data.AC2_1;
AC2_3 = Acc_Data.AC2_3;
AC3_1 = Acc_Data.AC3_1;
AC3_3 = Acc_Data.AC3_3;
AC4_1 = Acc_Data.AC4_1;
AC4_3 = Acc_Data.AC4_3; 


t_wind = Wind_Data.resultsTable_UA1.Time_Start;
u1 = Wind_Data.resultsTable_UA1.U;
u2 = Wind_Data.resultsTable_UA2.U;
u3 = Wind_Data.resultsTable_UA3.U;
u4 = Wind_Data.resultsTable_UA4.U;
u5 = Wind_Data.resultsTable_UA5.U;
u6 = Wind_Data.resultsTable_UA6.U;

aoa1 = Wind_Data.resultsTable_UA1.alpha_deg_mean;
aoa2 = Wind_Data.resultsTable_UA2.alpha_deg_mean;
aoa3 = Wind_Data.resultsTable_UA3.alpha_deg_mean;
aoa4 = Wind_Data.resultsTable_UA4.alpha_deg_mean;
aoa5 = Wind_Data.resultsTable_UA5.alpha_deg_mean;
aoa6 = Wind_Data.resultsTable_UA6.alpha_deg_mean;


%% Plot
total_plots = 5; % 或任何你需要的子图数量
current_plot = 1;
num_figs_in_row = [];
figWidthFactor = 1.5;
figPosition = [100, 100];
newfigure = true;
holdon = false;
firstfigure = true;


if fig_bool
 create_subplot(@plot, total_plots, current_plot, {t_wind,u1}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
 hold on
 create_subplot(@plot, total_plots, current_plot, {t_wind,u2}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
 create_subplot(@plot, total_plots, current_plot, {t_wind,u3}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
 create_subplot(@plot, total_plots, current_plot, {t_wind,u4}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
 create_subplot(@plot, total_plots, current_plot, {t_wind,u5}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
 create_subplot(@plot, total_plots, current_plot, {t_wind,u6}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
 title("U")
 current_plot = current_plot+1;
 legend

 create_subplot(@plot, total_plots, current_plot, {t_wind,aoa1}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
 hold on
 create_subplot(@plot, total_plots, current_plot, {t_wind,aoa2}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
 create_subplot(@plot, total_plots, current_plot, {t_wind,aoa3}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
 create_subplot(@plot, total_plots, current_plot, {t_wind,aoa4}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
 create_subplot(@plot, total_plots, current_plot, {t_wind,aoa5}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
 create_subplot(@plot, total_plots, current_plot, {t_wind,aoa6}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
 title("AOA")
 current_plot = current_plot+1;
 legend


  create_subplot(@plot, total_plots, current_plot, {Wind_Data_all.Time,Wind_Data_all.UA1_z}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
 hold on
 create_subplot(@plot, total_plots, current_plot, {Wind_Data_all.Time,Wind_Data_all.UA2_z}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
 create_subplot(@plot, total_plots, current_plot, {Wind_Data_all.Time,Wind_Data_all.UA3_z}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
 create_subplot(@plot, total_plots, current_plot, {Wind_Data_all.Time,Wind_Data_all.UA4_z}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
 create_subplot(@plot, total_plots, current_plot, {Wind_Data_all.Time,Wind_Data_all.UA5_z}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
 create_subplot(@plot, total_plots, current_plot, {Wind_Data_all.Time,Wind_Data_all.UA6_z}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
 title("z")
 current_plot = current_plot+1;
 legend

   create_subplot(@plot, total_plots, current_plot, {Wind_Data_all.Time,Wind_Data_all.UA1_x}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
 hold on
 create_subplot(@plot, total_plots, current_plot, {Wind_Data_all.Time,Wind_Data_all.UA2_z}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
 create_subplot(@plot, total_plots, current_plot, {Wind_Data_all.Time,Wind_Data_all.UA3_z}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
 create_subplot(@plot, total_plots, current_plot, {Wind_Data_all.Time,Wind_Data_all.UA4_z}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
 create_subplot(@plot, total_plots, current_plot, {Wind_Data_all.Time,Wind_Data_all.UA5_z}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
 create_subplot(@plot, total_plots, current_plot, {Wind_Data_all.Time,Wind_Data_all.UA6_z}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
 title("x")
 current_plot = current_plot+1;
 legend

    create_subplot(@plot, total_plots, current_plot, {Wind_Data_all.Time,Wind_Data_all.UA1_y}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
 hold on
 create_subplot(@plot, total_plots, current_plot, {Wind_Data_all.Time,Wind_Data_all.UA2_z}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
 create_subplot(@plot, total_plots, current_plot, {Wind_Data_all.Time,Wind_Data_all.UA3_z}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
 create_subplot(@plot, total_plots, current_plot, {Wind_Data_all.Time,Wind_Data_all.UA4_z}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
 create_subplot(@plot, total_plots, current_plot, {Wind_Data_all.Time,Wind_Data_all.UA5_z}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
 create_subplot(@plot, total_plots, current_plot, {Wind_Data_all.Time,Wind_Data_all.UA6_z}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
 title("y")
 current_plot = current_plot+1;
 legend
end