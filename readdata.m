%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn xushengyichn@outlook.com
%Date: 2023-09-19 02:17:00
%LastEditors: xushengyichn xushengyichn@outlook.com
%LastEditTime: 2023-09-19 02:24:45
%FilePath: \xihoumen_inverse_force_estimation\readdata.m
%Description: 
%
%Copyright (c) 2023 by ${git_name_email}, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;clear;close all


addpath(genpath("F:\git\ssm_tools"))
addpath(genpath("D:\Users\xushe\Documents\GitHub\ssm_tools"))
addpath(genpath("D:\Users\xushe\Documents\GitHub\Function_shengyi_package"))
addpath(genpath("F:\git\Function_shengyi_package"))
addpath(genpath("D:\OneDrive\NAS云同步\Drive\0博士研究生\3大论文\研究内容\研究内容 3：气动力模型及参数识别；\反算气动力"))
subStreamNumberDefault = 2132;
run('InitScript.m');
%% 0 绘图参数
fig_bool = ON;
num_figs_in_row = 5; %每一行显示几个图
figPos = figPosSmall; %图的大小，参数基于InitScript.m中的设置
%设置图片间隔
gap_between_images = [0, 0];
figureIdx = 0;

%% vibration
vibac2_name=["2013-02-06 00-vibac2.txt";"2013-02-06 01-vibac2.txt";"2013-02-06 02-vibac2.txt"];
vibac3_name=["2013-02-06 01-VIBac3.txt";"2013-02-06 01-VIBac3.txt";"2013-02-06 02-VIBac3.txt"];
vibac4_name=["2013-02-06 00-VIBac4.txt";"2013-02-06 01-VIBac4.txt";"2013-02-06 02-VIBac4.txt"];

vibac2_data=read_vib_data(vibac2_name);
vibac3_data=read_vib_data(vibac3_name);
vibac4_data=read_vib_data(vibac4_name);

dt = 0.02;
t = 0:dt:dt*(length(vibac2_data(:,1))-1);
startDate = datetime(2013,2,6,0,0,0);
t = startDate:seconds(dt):startDate+(length(vibac2_data(:,1))-1)*seconds(dt);
%% wind
% winds1_name = ["wind_property_result_s1.txt"];

winds1_name =["wind_property_result_s1.txt"];
winds3_name =["wind_property_result_s3.txt"];
winds4_name =["wind_property_result_s4.txt"];
winds1_data = read_wind_data(winds1_name);
winds3_data = read_wind_data(winds3_name);
winds4_data = read_wind_data(winds4_name);

% startDate = datetime(2013,2,6,0,0,0);
endDate = datetime(2013,2,6,3,0,0);
winds1_trim = winds1_data(winds1_data.DateTime >= startDate & winds1_data.DateTime <= endDate, :);
winds3_trim = winds3_data(winds3_data.DateTime >= startDate & winds3_data.DateTime <= endDate, :);
winds4_trim = winds4_data(winds4_data.DateTime >= startDate & winds4_data.DateTime <= endDate, :);

%% plot
[figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);

plot(t,vibac2_data(:,2), 'Color', 'b', 'LineStyle', '--')
hold on
plot(t,vibac2_data(:,3), 'Color', 'k', 'LineStyle', '--')
plot(t,vibac2_data(:,4), 'Color', 'r', 'LineStyle', '--')
set(gca, 'FontSize', 12)
xlabel('Time')
ylabel('acc')
% legend('filtered', 'true', 'Location', 'northwest')
title(['vib2']);
% xlim([0, 3])
% ylim([0, 50])



[figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);

plot(t,vibac3_data(:,2), 'Color', 'b', 'LineStyle', '--')
hold on
plot(t,vibac3_data(:,3), 'Color', 'k', 'LineStyle', '--')
plot(t,vibac3_data(:,4), 'Color', 'r', 'LineStyle', '--')
set(gca, 'FontSize', 12)
xlabel('Time')
ylabel('acc')
% legend('filtered', 'true', 'Location', 'northwest')
title(['vib3']);
% xlim([0, 3])
% ylim([0, 50])


[figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);

plot(t,vibac4_data(:,2), 'Color', 'b', 'LineStyle', '--')
hold on
plot(t,vibac4_data(:,3), 'Color', 'k', 'LineStyle', '--')
plot(t,vibac4_data(:,4), 'Color', 'r', 'LineStyle', '--')
set(gca, 'FontSize', 12)
xlabel('Time')
ylabel('acc')
% legend('filtered', 'true', 'Location', 'northwest')
title(['vib4']);
% xlim([0, 3])
% ylim([0, 50])

[figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);

plot(winds1_trim.DateTime,winds1_trim.u, 'Color', 'b', 'LineStyle', '--')
hold on
plot(winds3_trim.DateTime,winds3_trim.u, 'Color', 'k', 'LineStyle', '--')
plot(winds4_trim.DateTime,winds4_trim.u, 'Color', 'r', 'LineStyle', '--')
set(gca, 'FontSize', 12)
xlabel('Time')
ylabel('u')
% legend('filtered', 'true', 'Location', 'northwest')
title(['wind velocity']);
% xlim([0, 3])
% ylim([0, 50])


function [vibac_data]=read_vib_data(vibac_name)
    vibac_data=[];  
    for i=1:length(vibac_name)
            vibac_data_temp=importdata(vibac_name(i));
            vibac_data=[vibac_data;vibac_data_temp];
        end
end

function data = read_wind_data(wind_name)
    data = readtable(wind_name);
    data.DateTime = datetime(data.day, 'InputFormat', 'yyyy-MM-dd') + hours(data.hour) + minutes(data.x10min*10);
    data.DateTime.Format = 'yyyy-MM-dd HH:mm';
    data.day = [];
    data.hour = [];
    data.x10min = [];
    data = movevars(data, 'DateTime', 'Before', data.Properties.VariableNames{1});
end