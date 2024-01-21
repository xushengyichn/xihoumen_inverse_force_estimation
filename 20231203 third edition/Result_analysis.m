%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: ShengyiXu xushengyichn@outlook.com
%Date: 2024-01-21 12:45:33
%LastEditors: ShengyiXu xushengyichn@outlook.com
%LastEditTime: 2024-01-21 13:07:43
%FilePath: \manuscriptf:\git\xihoumen_inverse_force_estimation\20231203 third edition\Result_analysis.m
%Description:
%
%Copyright (c) 2024 by ${git_name_email}, All Rights Reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 在被调用的脚本中
if ~exist('skipClear', 'var')
    clc; clear; close all;
end

run('CommonCommand.m');

% 定义VIV_sel值的数组，例如：[1, 2, 3, ...]
VIV_sels = [2,3,4,5,6]; % 根据你的数据集进行修改
% 定义不同的标记样式
% 定义15种不同的标记样式和颜色
markers = {'o', '+', '*', '.', 'x', 's', 'd', '^', 'v', '>', '<', 'p', 'h'};
colors = {'r', 'g', 'b', 'c', 'm', 'y', 'k', [.5 .6 .7], [.8 .2 .6], [.2 .5 .8], [.3 .7 .9], [.4 .4 .4], [.6 .2 .2], [.7 .5 .3], [.1 .3 .5]};

opts = detectImportOptions('vivData.csv');
opts = setvaropts(opts, 'startDate', 'InputFormat', 'MM/dd/yyyy HH:mm');
opts = setvaropts(opts, 'endDate', 'InputFormat', 'MM/dd/yyyy HH:mm');
opts = setvaropts(opts, 'startDate_update', 'InputFormat', 'MM/dd/yyyy HH:mm');
opts = setvaropts(opts, 'endDate_update', 'InputFormat', 'MM/dd/yyyy HH:mm');

vivTable = readtable('vivData.csv',opts);

% 初始化图形参数
total_plots = 2; % 或任何你需要的子图数量
current_plot = 1;
num_figs_in_row = [];
figWidthFactor = 1.5;
figPosition = [100, 100];
newfigure = true;
holdon = false;
firstfigure = true;

% 在一个图上绘制所有VIV_sel值对应的数据
for i = 1:length(VIV_sels)
    VIV_sel = VIV_sels(i);
    start_time = vivTable.startDate(VIV_sel);
    end_time = vivTable.endDate(VIV_sel);

    start_time.Format = 'yyyy_MM_dd_HH_mm';
    filename = sprintf('result_%s_%s.mat', start_time, num2str(VIV_sel));
    load(filename);

    % 绘制散点图
    create_subplot(@scatter3, total_plots, current_plot, {amp_filter, U_sel,  zeta_filter,[],AoA_sel,markers{i}}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);
    hold on
    xlabel('Amp. (m)')
    ylabel('Wind speed (m/s)')
    title("Damping ratio （Amp,U,damping ratio,,AOA）")

    firstfigure = false;
    
end
current_plot = current_plot + 1;
% 创建透明平面
x = [min(amp_filter), max(amp_filter), max(amp_filter), min(amp_filter)];
y = [min(U_sel), min(U_sel), max(U_sel), max(U_sel)];
z = [-0.003, -0.003, -0.003, -0.003];
patch(x, y, z, 'blue', 'FaceAlpha', 0.3);


% 在一个图上绘制所有VIV_sel值对应的数据
for i = 1:length(VIV_sels)
    VIV_sel = VIV_sels(i);
    start_time = vivTable.startDate(VIV_sel);
    end_time = vivTable.endDate(VIV_sel);

    start_time.Format = 'yyyy_MM_dd_HH_mm';
    filename = sprintf('result_%s_%s.mat', start_time, num2str(VIV_sel));
    load(filename);

    % 绘制散点图
    create_subplot(@scatter3, total_plots, current_plot, {amp_filter, U_sel,  zeta_filter,[], colors{i}}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);
    hold on
    xlabel('Amp. (m)')
    ylabel('Wind speed (m/s)')
    title("Damping ratio （Amp,U,damping ratio,,AOA）")

    firstfigure = false;
    
end
current_plot = current_plot + 1;
% 创建透明平面
x = [min(amp_filter), max(amp_filter), max(amp_filter), min(amp_filter)];
y = [min(U_sel), min(U_sel), max(U_sel), max(U_sel)];
z = [-0.003, -0.003, -0.003, -0.003];
patch(x, y, z, 'blue', 'FaceAlpha', 0.3);

holdon=false;