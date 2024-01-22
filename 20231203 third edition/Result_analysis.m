%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: ShengyiXu xushengyichn@outlook.com
%Date: 2024-01-21 12:45:33
%LastEditors: ShengyiXu xushengyichn@outlook.com
%LastEditTime: 2024-01-22 23:01:41
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
% VIV_sels = [2,3,4,5,6]; % 根据你的数据集进行修改
VIV_sels = [2;3;4;5;6;7;8;9;10;11;12;16;17;18;19;22];

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



%% read data
% 在一个图上绘制所有VIV_sel值对应的数据
for i = 1:length(VIV_sels)
    modelupdate = true;
    VIV_sel = VIV_sels(i);
    start_time = vivTable.startDate(VIV_sel);
    end_time = vivTable.endDate(VIV_sel);

    start_time.Format = 'yyyy_MM_dd_HH_mm';
    if modelupdate == true
        str_modeupdate = "updatedmodel";
    else
        str_modeupdate = "nomodelupdate";
    end
    filename = sprintf('result_%s_%s_%s.mat', start_time, num2str(VIV_sel),str_modeupdate);
    % filename = sprintf('result_%s_%s.mat', start_time, num2str(VIV_sel));
    data_temp=load(filename);
    plotdata_modelupdate(i).amp_filter = data_temp.amp_filter;
    plotdata_modelupdate(i).zeta_filter = data_temp.zeta_filter;
    plotdata_modelupdate(i).AoA_sel_1 = data_temp.AoA_sel_1;
    plotdata_modelupdate(i).AoA_sel_2 = data_temp.AoA_sel_2;
    plotdata_modelupdate(i).AoA_sel_3 = data_temp.AoA_sel_3;
    plotdata_modelupdate(i).U_sel_1 = data_temp.U_sel_loc_1;
    plotdata_modelupdate(i).U_sel_2 = data_temp.U_sel_loc_2;
    plotdata_modelupdate(i).U_sel_3 = data_temp.U_sel_loc_3;
    plotdata_modelupdate(i).zeta_structure = data_temp.zeta;
    plotdata_modelupdate(i).VIV_mode_seq = data_temp.VIV_mode_seq;
    plotdata_modelupdate(i).beta_deg_mean_UA1a = data_temp.beta_deg_mean_UA1a;
    plotdata_modelupdate(i).beta_deg_mean_UA1b = data_temp.beta_deg_mean_UA1b;
    plotdata_modelupdate(i).beta_deg_mean_UA2a = data_temp.beta_deg_mean_UA2a;
    plotdata_modelupdate(i).beta_deg_mean_UA2b = data_temp.beta_deg_mean_UA2b;
    plotdata_modelupdate(i).beta_deg_mean_UA3a = data_temp.beta_deg_mean_UA3a;
    plotdata_modelupdate(i).beta_deg_mean_UA3b = data_temp.beta_deg_mean_UA3b;
    plotdata_modelupdate(i).t_cycle_mean_temp = data_temp.t_cycle_mean_temp;

end


% 在一个图上绘制所有VIV_sel值对应的数据
for i = 1:length(VIV_sels)
    modelupdate = false;
    VIV_sel = VIV_sels(i);
    start_time = vivTable.startDate(VIV_sel);
    end_time = vivTable.endDate(VIV_sel);

    start_time.Format = 'yyyy_MM_dd_HH_mm';
    if modelupdate == true
        str_modeupdate = "updatedmodel";
    else
        str_modeupdate = "nomodelupdate";
    end
    filename = sprintf('result_%s_%s_%s.mat', start_time, num2str(VIV_sel),str_modeupdate);
    % filename = sprintf('result_%s_%s.mat', start_time, num2str(VIV_sel));
    data_temp=load(filename);
    plotdata_nomodelupdate(i).amp_filter = data_temp.amp_filter;
    plotdata_nomodelupdate(i).zeta_filter = data_temp.zeta_filter;
    plotdata_nomodelupdate(i).AoA_sel_1 = data_temp.AoA_sel_1;
    plotdata_nomodelupdate(i).AoA_sel_2 = data_temp.AoA_sel_2;
    plotdata_nomodelupdate(i).AoA_sel_3 = data_temp.AoA_sel_3;
    plotdata_nomodelupdate(i).U_sel_1 = data_temp.U_sel_loc_1;
    plotdata_nomodelupdate(i).U_sel_2 = data_temp.U_sel_loc_2;
    plotdata_nomodelupdate(i).U_sel_3 = data_temp.U_sel_loc_3;
    plotdata_nomodelupdate(i).zeta_structure = data_temp.zeta;
    plotdata_nomodelupdate(i).VIV_mode_seq = data_temp.VIV_mode_seq;
    plotdata_nomodelupdate(i).beta_deg_mean_UA1a = data_temp.beta_deg_mean_UA1a;
    plotdata_nomodelupdate(i).beta_deg_mean_UA1b = data_temp.beta_deg_mean_UA1b;
    plotdata_nomodelupdate(i).beta_deg_mean_UA2a = data_temp.beta_deg_mean_UA2a;
    plotdata_nomodelupdate(i).beta_deg_mean_UA2b = data_temp.beta_deg_mean_UA2b;
    plotdata_nomodelupdate(i).beta_deg_mean_UA3a = data_temp.beta_deg_mean_UA3a;
    plotdata_nomodelupdate(i).beta_deg_mean_UA3b = data_temp.beta_deg_mean_UA3b;
    plotdata_nomodelupdate(i).t_cycle_mean_temp = data_temp.t_cycle_mean_temp;

end




%% print structural information
for k1 = 1:length(plotdata_modelupdate)
    zeta_temp = plotdata_modelupdate(k1).zeta_structure;
    VIV_mode_seq = plotdata_modelupdate(k1).VIV_mode_seq;
    zeta_VIV_mode_modelupdate(k1,:) = zeta_temp(VIV_mode_seq);
end

for k1 = 1:length(plotdata_nomodelupdate)
    zeta_temp = plotdata_nomodelupdate(k1).zeta_structure;
    VIV_mode_seq = plotdata_nomodelupdate(k1).VIV_mode_seq;
    zeta_VIV_mode_nomodelupdate(k1,:) = zeta_temp(VIV_mode_seq);
end

% create table 
zeta_compare = table(VIV_sels,zeta_VIV_mode_modelupdate,zeta_VIV_mode_nomodelupdate,'VariableNames',{'VIV_sel','zeta_VIV_mode_modelupdate','zeta_VIV_mode_nomodelupdate'});
% disp(zeta_compare)
%% print amplitude information
for k1 =1:length(plotdata_modelupdate)
    amp_temp =  plotdata_modelupdate(k1).amp_filter;
    amp_rms(k1,1) = rms(amp_temp);
end
    amp_rms_table = table(VIV_sels,amp_rms);
% disp(amp_rms_table)
table_all  = table(VIV_sels,zeta_VIV_mode_modelupdate,zeta_VIV_mode_nomodelupdate,amp_rms,'VariableNames',{'VIV_sel','zeta_VIV_mode_modelupdate','zeta_VIV_mode_nomodelupdate','amp_rms'});
disp(table_all)


%% select data
for i = 1:length(VIV_sels)
    zeta_filter_temp = plotdata_modelupdate(i).zeta_filter;
    zeta_struct_temp = table_all.zeta_VIV_mode_modelupdate(i);
    loc_temp = find(zeta_filter_temp <= -zeta_struct_temp);
    plotdata_modelupdate_sel_excite(i).amp_filter         = plotdata_modelupdate(i).amp_filter(loc_temp);
    plotdata_modelupdate_sel_excite(i).zeta_filter        = plotdata_modelupdate(i).zeta_filter(loc_temp);
    plotdata_modelupdate_sel_excite(i).AoA_sel_1          = plotdata_modelupdate(i).AoA_sel_1(loc_temp);
    plotdata_modelupdate_sel_excite(i).AoA_sel_2          = plotdata_modelupdate(i).AoA_sel_2(loc_temp);
    plotdata_modelupdate_sel_excite(i).AoA_sel_3          = plotdata_modelupdate(i).AoA_sel_3(loc_temp);
    plotdata_modelupdate_sel_excite(i).U_sel_1            = plotdata_modelupdate(i).U_sel_1(loc_temp);
    plotdata_modelupdate_sel_excite(i).U_sel_2            = plotdata_modelupdate(i).U_sel_2(loc_temp);
    plotdata_modelupdate_sel_excite(i).U_sel_3            = plotdata_modelupdate(i).U_sel_3(loc_temp);
    plotdata_modelupdate_sel_excite(i).zeta_structure     = plotdata_modelupdate(i).zeta_structure;
    plotdata_modelupdate_sel_excite(i).VIV_mode_seq       = plotdata_modelupdate(i).VIV_mode_seq; 
    plotdata_modelupdate_sel_excite(i).beta_deg_mean_UA1a = plotdata_modelupdate(i).beta_deg_mean_UA1a(loc_temp);
    plotdata_modelupdate_sel_excite(i).beta_deg_mean_UA1b = plotdata_modelupdate(i).beta_deg_mean_UA1b(loc_temp);
    plotdata_modelupdate_sel_excite(i).beta_deg_mean_UA2a = plotdata_modelupdate(i).beta_deg_mean_UA2a(loc_temp);
    plotdata_modelupdate_sel_excite(i).beta_deg_mean_UA2b = plotdata_modelupdate(i).beta_deg_mean_UA2b(loc_temp);
    plotdata_modelupdate_sel_excite(i).beta_deg_mean_UA3a = plotdata_modelupdate(i).beta_deg_mean_UA3a(loc_temp);
    plotdata_modelupdate_sel_excite(i).beta_deg_mean_UA3b = plotdata_modelupdate(i).beta_deg_mean_UA3b(loc_temp);
    plotdata_modelupdate_sel_excite(i).t_cycle_mean_temp  = plotdata_modelupdate(i).t_cycle_mean_temp(loc_temp);
end

for i = 1:length(VIV_sels)
    zeta_filter_temp = plotdata_nomodelupdate(i).zeta_filter;
    zeta_struct_temp = table_all.zeta_VIV_mode_nomodelupdate(i);
    loc_temp = find(zeta_filter_temp <= -zeta_struct_temp);
    plotdata_nomodelupdate_sel_excite(i).amp_filter         = plotdata_nomodelupdate(i).amp_filter(loc_temp);
    plotdata_nomodelupdate_sel_excite(i).zeta_filter        = plotdata_nomodelupdate(i).zeta_filter(loc_temp);
    plotdata_nomodelupdate_sel_excite(i).AoA_sel_1          = plotdata_nomodelupdate(i).AoA_sel_1(loc_temp);
    plotdata_nomodelupdate_sel_excite(i).AoA_sel_2          = plotdata_nomodelupdate(i).AoA_sel_2(loc_temp);
    plotdata_nomodelupdate_sel_excite(i).AoA_sel_3          = plotdata_nomodelupdate(i).AoA_sel_3(loc_temp);
    plotdata_nomodelupdate_sel_excite(i).U_sel_1            = plotdata_nomodelupdate(i).U_sel_1(loc_temp);
    plotdata_nomodelupdate_sel_excite(i).U_sel_2            = plotdata_nomodelupdate(i).U_sel_2(loc_temp);
    plotdata_nomodelupdate_sel_excite(i).U_sel_3            = plotdata_nomodelupdate(i).U_sel_3(loc_temp);
    plotdata_nomodelupdate_sel_excite(i).zeta_structure     = plotdata_nomodelupdate(i).zeta_structure;
    plotdata_nomodelupdate_sel_excite(i).VIV_mode_seq       = plotdata_nomodelupdate(i).VIV_mode_seq; 
    plotdata_nomodelupdate_sel_excite(i).beta_deg_mean_UA1a = plotdata_nomodelupdate(i).beta_deg_mean_UA1a(loc_temp);
    plotdata_nomodelupdate_sel_excite(i).beta_deg_mean_UA1b = plotdata_nomodelupdate(i).beta_deg_mean_UA1b(loc_temp);
    plotdata_nomodelupdate_sel_excite(i).beta_deg_mean_UA2a = plotdata_nomodelupdate(i).beta_deg_mean_UA2a(loc_temp);
    plotdata_nomodelupdate_sel_excite(i).beta_deg_mean_UA2b = plotdata_nomodelupdate(i).beta_deg_mean_UA2b(loc_temp);
    plotdata_nomodelupdate_sel_excite(i).beta_deg_mean_UA3a = plotdata_nomodelupdate(i).beta_deg_mean_UA3a(loc_temp);
    plotdata_nomodelupdate_sel_excite(i).beta_deg_mean_UA3b = plotdata_nomodelupdate(i).beta_deg_mean_UA3b(loc_temp);
    plotdata_nomodelupdate_sel_excite(i).t_cycle_mean_temp  = plotdata_nomodelupdate(i).t_cycle_mean_temp(loc_temp);
end


for i = 1:length(VIV_sels)
    zeta_filter_temp = plotdata_modelupdate(i).zeta_filter;
    zeta_struct_temp = table_all.zeta_VIV_mode_modelupdate(i);
    loc_temp = find(zeta_filter_temp > -zeta_struct_temp);
    plotdata_modelupdate_sel_drop(i).amp_filter         = plotdata_modelupdate(i).amp_filter(loc_temp);
    plotdata_modelupdate_sel_drop(i).zeta_filter        = plotdata_modelupdate(i).zeta_filter(loc_temp);
    plotdata_modelupdate_sel_drop(i).AoA_sel_1          = plotdata_modelupdate(i).AoA_sel_1(loc_temp);
    plotdata_modelupdate_sel_drop(i).AoA_sel_2          = plotdata_modelupdate(i).AoA_sel_2(loc_temp);
    plotdata_modelupdate_sel_drop(i).AoA_sel_3          = plotdata_modelupdate(i).AoA_sel_3(loc_temp);
    plotdata_modelupdate_sel_drop(i).U_sel_1            = plotdata_modelupdate(i).U_sel_1(loc_temp);
    plotdata_modelupdate_sel_drop(i).U_sel_2            = plotdata_modelupdate(i).U_sel_2(loc_temp);
    plotdata_modelupdate_sel_drop(i).U_sel_3            = plotdata_modelupdate(i).U_sel_3(loc_temp);
    plotdata_modelupdate_sel_drop(i).zeta_structure     = plotdata_modelupdate(i).zeta_structure;
    plotdata_modelupdate_sel_drop(i).VIV_mode_seq       = plotdata_modelupdate(i).VIV_mode_seq; 
    plotdata_modelupdate_sel_drop(i).beta_deg_mean_UA1a = plotdata_modelupdate(i).beta_deg_mean_UA1a(loc_temp);
    plotdata_modelupdate_sel_drop(i).beta_deg_mean_UA1b = plotdata_modelupdate(i).beta_deg_mean_UA1b(loc_temp);
    plotdata_modelupdate_sel_drop(i).beta_deg_mean_UA2a = plotdata_modelupdate(i).beta_deg_mean_UA2a(loc_temp);
    plotdata_modelupdate_sel_drop(i).beta_deg_mean_UA2b = plotdata_modelupdate(i).beta_deg_mean_UA2b(loc_temp);
    plotdata_modelupdate_sel_drop(i).beta_deg_mean_UA3a = plotdata_modelupdate(i).beta_deg_mean_UA3a(loc_temp);
    plotdata_modelupdate_sel_drop(i).beta_deg_mean_UA3b = plotdata_modelupdate(i).beta_deg_mean_UA3b(loc_temp);
    plotdata_modelupdate_sel_drop(i).t_cycle_mean_temp  = plotdata_modelupdate(i).t_cycle_mean_temp(loc_temp);
end

for i = 1:length(VIV_sels)
    zeta_filter_temp = plotdata_nomodelupdate(i).zeta_filter;
    zeta_struct_temp = table_all.zeta_VIV_mode_nomodelupdate(i);
    loc_temp = find(zeta_filter_temp > -zeta_struct_temp);
    plotdata_nomodelupdate_sel_drop(i).amp_filter         = plotdata_nomodelupdate(i).amp_filter(loc_temp);
    plotdata_nomodelupdate_sel_drop(i).zeta_filter        = plotdata_nomodelupdate(i).zeta_filter(loc_temp);
    plotdata_nomodelupdate_sel_drop(i).AoA_sel_1          = plotdata_nomodelupdate(i).AoA_sel_1(loc_temp);
    plotdata_nomodelupdate_sel_drop(i).AoA_sel_2          = plotdata_nomodelupdate(i).AoA_sel_2(loc_temp);
    plotdata_nomodelupdate_sel_drop(i).AoA_sel_3          = plotdata_nomodelupdate(i).AoA_sel_3(loc_temp);
    plotdata_nomodelupdate_sel_drop(i).U_sel_1            = plotdata_nomodelupdate(i).U_sel_1(loc_temp);
    plotdata_nomodelupdate_sel_drop(i).U_sel_2            = plotdata_nomodelupdate(i).U_sel_2(loc_temp);
    plotdata_nomodelupdate_sel_drop(i).U_sel_3            = plotdata_nomodelupdate(i).U_sel_3(loc_temp);
    plotdata_nomodelupdate_sel_drop(i).zeta_structure     = plotdata_nomodelupdate(i).zeta_structure;
    plotdata_nomodelupdate_sel_drop(i).VIV_mode_seq       = plotdata_nomodelupdate(i).VIV_mode_seq; 
    plotdata_nomodelupdate_sel_drop(i).beta_deg_mean_UA1a = plotdata_nomodelupdate(i).beta_deg_mean_UA1a(loc_temp);
    plotdata_nomodelupdate_sel_drop(i).beta_deg_mean_UA1b = plotdata_nomodelupdate(i).beta_deg_mean_UA1b(loc_temp);
    plotdata_nomodelupdate_sel_drop(i).beta_deg_mean_UA2a = plotdata_nomodelupdate(i).beta_deg_mean_UA2a(loc_temp);
    plotdata_nomodelupdate_sel_drop(i).beta_deg_mean_UA2b = plotdata_nomodelupdate(i).beta_deg_mean_UA2b(loc_temp);
    plotdata_nomodelupdate_sel_drop(i).beta_deg_mean_UA3a = plotdata_nomodelupdate(i).beta_deg_mean_UA3a(loc_temp);
    plotdata_nomodelupdate_sel_drop(i).beta_deg_mean_UA3b = plotdata_nomodelupdate(i).beta_deg_mean_UA3b(loc_temp);
    plotdata_nomodelupdate_sel_drop(i).t_cycle_mean_temp  = plotdata_nomodelupdate(i).t_cycle_mean_temp(loc_temp);
end

%% 初始化图形参数
total_plots = 16; % 或任何你需要的子图数量
current_plot = 1;
num_figs_in_row = [];
figWidthFactor = 1.5;
figPosition = [100, 100];
newfigure = true;
holdon = false;
firstfigure = true;

%% wind direction for different wind sensors
VIV_sel  =  [2,3,4,5,6,7,8,9,10,11,12,16,17,18,19,22];
plot_sel = find(ismember(VIV_sels, VIV_sel));
plotdata_sel = plotdata_modelupdate;
for k1 = 1:length(plot_sel)
x = plotdata_sel(k1).t_cycle_mean_temp;
y = plotdata_sel(k1).beta_deg_mean_UA1a;
create_subplot(@scatter, total_plots, current_plot, {x, y}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
hold on
y = plotdata_sel(k1).beta_deg_mean_UA1b;
create_subplot(@scatter, total_plots, current_plot, {x, y}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
y = plotdata_sel(k1).beta_deg_mean_UA2a;
create_subplot(@scatter, total_plots, current_plot, {x, y}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
y = plotdata_sel(k1).beta_deg_mean_UA2b;
create_subplot(@scatter, total_plots, current_plot, {x, y}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', false);firstfigure = false;
y = plotdata_sel(k1).beta_deg_mean_UA3a;
create_subplot(@scatter, total_plots, current_plot, {x, y}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', false);firstfigure = false;
y = plotdata_sel(k1).beta_deg_mean_UA3b;
create_subplot(@scatter, total_plots, current_plot, {x, y}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', false);firstfigure = false;
xlabel('Time')
ylabel('Wind direction (deg)')
legend('UA1a','UA1b','UA2a','UA2b','UA3a','UA3b')
current_plot = current_plot + 1;
title("VIV_sel:"+num2str(VIV_sel(k1)))
end

%% 初始化图形参数
total_plots = 16; % 或任何你需要的子图数量
current_plot = 1;
num_figs_in_row = [];
figWidthFactor = 1.5;
figPosition = [100, 100];
newfigure = true;
holdon = false;
firstfigure = true;
%% wind velocity for different wind sensors


plotdata_sel = plotdata_modelupdate;
plot_sel = find(ismember(VIV_sels, VIV_sel));
for k1 = 1:length(plot_sel)
x = plotdata_sel(k1).t_cycle_mean_temp;
y = plotdata_sel(k1).U_sel_1;
create_subplot(@scatter, total_plots, current_plot, {x, y}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
hold on
y = plotdata_sel(k1).U_sel_2;
create_subplot(@scatter, total_plots, current_plot, {x, y}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
y = plotdata_sel(k1).U_sel_3;
create_subplot(@scatter, total_plots, current_plot, {x, y}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;

xlabel('Time')
ylabel('Wind speed location (m/s)')
legend('location 1','location 2','location 3')
title("VIV_sel:"+num2str(VIV_sel(k1)))
current_plot = current_plot + 1;
end
%% 初始化图形参数
total_plots = 16; % 或任何你需要的子图数量
current_plot = 1;
num_figs_in_row = [];
figWidthFactor = 1.5;
figPosition = [100, 100];
newfigure = true;
holdon = false;
firstfigure = true;

%% AOA for different wind sensors

plotdata_sel = plotdata_modelupdate;
for k1 = 1:length(plot_sel)
plot_sel = find(ismember(VIV_sels, VIV_sel));

x = plotdata_sel(k1).t_cycle_mean_temp;
y = plotdata_sel(k1).AoA_sel_1;
create_subplot(@scatter, total_plots, current_plot, {x, y}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
hold on
y = plotdata_sel(k1).AoA_sel_2;
create_subplot(@scatter, total_plots, current_plot, {x, y}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
y = plotdata_sel(k1).AoA_sel_3;
create_subplot(@scatter, total_plots, current_plot, {x, y}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;

xlabel('Time')
ylabel('AOA (m/s)')
legend('location 1','location 2','location 3')
title("VIV_sel:"+num2str(VIV_sel(k1)))
ylim([-4,4])
current_plot = current_plot + 1;
end


%% new figure
total_plots = 3; % 或任何你需要的子图数量
current_plot = 1;
num_figs_in_row = [];
figWidthFactor = 1.5;
figPosition = [100, 100];
newfigure = true;
holdon = false;
firstfigure = true;

%% camparison between model updating or not
VIV_sel  = [9,12];
plot_sel = find(ismember(VIV_sels, VIV_sel));

for k1 = 1:length(plot_sel)
    plotdata_sel = plotdata_modelupdate;
    x= plotdata_sel(plot_sel(k1)).amp_filter;
    y= plotdata_sel(plot_sel(k1)).zeta_filter;
    create_subplot(@scatter, total_plots, current_plot, {x, y}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
    hold on
    plotdata_sel = plotdata_nomodelupdate;
    x= plotdata_sel(plot_sel(k1)).amp_filter;
    y= plotdata_sel(plot_sel(k1)).zeta_filter;
    create_subplot(@scatter, total_plots, current_plot, {x, y}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
    hold on
    ylim([-0.01,0.1/100])
    current_plot = current_plot + 1;
end

%% plot all the VIV cases for damping ratio with amplitude model updating
VIV_sel  = VIV_sels;
plot_sel = find(ismember(VIV_sels, VIV_sel));

for k1 = 1:length(plot_sel)
    plotdata_sel = plotdata_modelupdate;
    x= plotdata_sel(plot_sel(k1)).amp_filter;
    y= plotdata_sel(plot_sel(k1)).zeta_filter;
    create_subplot(@scatter, total_plots, current_plot, {x, y}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
    hold on
    ylim([-0.01,0.1/100])
    
end
xlabel("amplitude")
ylabel("damping ratio")
title("all the VIV cases (model updating)")
current_plot = current_plot + 1;
%% plot all the VIV cases for damping ratio with amplitude no model updating
VIV_sel  = VIV_sels;
plot_sel = find(ismember(VIV_sels, VIV_sel));

for k1 = 1:length(plot_sel)
    plotdata_sel = plotdata_nomodelupdate;
    x= plotdata_sel(plot_sel(k1)).amp_filter;
    y= plotdata_sel(plot_sel(k1)).zeta_filter;
    create_subplot(@scatter, total_plots, current_plot, {x, y}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
    hold on
    ylim([-0.01,0.1/100])
    
end

xlabel("amplitude")
ylabel("damping ratio")
title("all the VIV cases (no model updating)")
current_plot = current_plot + 1;

%% new figure
total_plots = 2; % 或任何你需要的子图数量
current_plot = 1;
num_figs_in_row = [];
figWidthFactor = 1.5;
figPosition = [100, 100];
newfigure = true;
holdon = false;
firstfigure = true;
%% different VIV case and use sensor 1 model updating
% VIV_sel  = [4,9,11];
VIV_sel  = [5,6,12,16];
plot_sel = find(ismember(VIV_sels, VIV_sel));
plotdata_sel = plotdata_modelupdate;
for k1 = 1:length(plot_sel)
    x = plotdata_sel(plot_sel(k1)).amp_filter;
    y = plotdata_sel(plot_sel(k1)).U_sel_1;
    z = plotdata_sel(plot_sel(k1)).zeta_filter;
    color_points = colors{k1};
    create_subplot(@scatter3, total_plots, current_plot, {x, y, z,[],color_points}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
    hold on
end
zlim([-0.02 0.02])
current_plot = current_plot + 1;
xlabel("amplitude")
ylabel("wind speed")
zlabel("damping ratio")
title("different VIV cases (model updating)")

%% different VIV case and use sensor 1 no model updating
VIV_sel  = VIV_sel;
plot_sel = find(ismember(VIV_sels, VIV_sel));
plotdata_sel = plotdata_nomodelupdate;
for k1 = 1:length(plot_sel)
    x = plotdata_sel(plot_sel(k1)).amp_filter;
    y = plotdata_sel(plot_sel(k1)).U_sel_1;
    z = plotdata_sel(plot_sel(k1)).zeta_filter;
    color_points = colors{k1};
    create_subplot(@scatter3, total_plots, current_plot, {x, y, z,[],color_points}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
    hold on
end
zlim([-0.02 0.02])
current_plot = current_plot + 1;
xlabel("amplitude")
ylabel("wind speed")
zlabel("damping ratio")
title("different VIV cases (no model updating)")


%% new figure
total_plots = 2; % 或任何你需要的子图数量
current_plot = 1;
num_figs_in_row = [];
figWidthFactor = 1.5;
figPosition = [100, 100];
newfigure = true;
holdon = false;
firstfigure = true;
%% different AOA case and use sensor 1 model updating
VIV_sel = VIV_sel;
plot_sel = find(ismember(VIV_sels, VIV_sel));

% 初始化图例标签数组
% legend_labels = cell(1, length(categories) - 1);
plotdata_sel = plotdata_modelupdate;
for k1 = 1:length(plot_sel)
    x = plotdata_sel(plot_sel(k1)).amp_filter;
    y = plotdata_sel(plot_sel(k1)).U_sel_1;
    z = plotdata_sel(plot_sel(k1)).zeta_filter;
    color_points = plotdata_sel(plot_sel(k1)).AoA_sel_1;

    % 分类 color_points
    categories = -6:1:10;
    [~, ~, categorized_color_points] = histcounts(color_points, categories);

    for k2 = 1:length(categories) - 1
        % 选择当前类别的数据点
        idx = categorized_color_points == k2;
        x_cat = x(idx);
        y_cat = y(idx);
        z_cat = z(idx);

        % 选择颜色
        color_index = mod(k2, length(colors)) + 1;
        plot_color = colors{color_index};

        % 绘制散点图并为每个类别添加图例
        create_subplot(@scatter3, total_plots, current_plot, {x_cat, y_cat, z_cat, [], plot_color}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;


        hold on;

        % 更新图例标签
        legend_labels{k2} = ['AOA ' num2str(categories(k2))];
    end
end

% 添加图例
legend(legend_labels);
zlim([-0.02 0.02])
xlabel("amplitude")
ylabel("wind speed")
zlabel("damping ratio")
title("different AOA cases (model updating)")
% 更新 current_plot
current_plot = current_plot + 1;

%% different AOA case and use sensor 1 no model updating
VIV_sel = VIV_sel;
plot_sel = find(ismember(VIV_sels, VIV_sel));

% 初始化图例标签数组
% legend_labels = cell(1, length(categories) - 1);
plotdata_sel = plotdata_nomodelupdate;
for k1 = 1:length(plot_sel)
    x = plotdata_sel(plot_sel(k1)).amp_filter;
    y = plotdata_sel(plot_sel(k1)).U_sel_1;
    z = plotdata_sel(plot_sel(k1)).zeta_filter;
    color_points = plotdata_sel(plot_sel(k1)).AoA_sel_1;

    % 分类 color_points
    categories = -6:1:10;
    [~, ~, categorized_color_points] = histcounts(color_points, categories);

    for k2 = 1:length(categories) - 1
        % 选择当前类别的数据点
        idx = categorized_color_points == k2;
        x_cat = x(idx);
        y_cat = y(idx);
        z_cat = z(idx);

        % 选择颜色
        color_index = mod(k2, length(colors)) + 1;
        plot_color = colors{color_index};

        % 绘制散点图并为每个类别添加图例
        create_subplot(@scatter3, total_plots, current_plot, {x_cat, y_cat, z_cat, [], plot_color}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;


        hold on;

        % 更新图例标签
        legend_labels{k2} = ['AOA ' num2str(categories(k2))];
    end
end

% 添加图例
legend(legend_labels);
zlim([-0.02 0.02])
xlabel("amplitude")
ylabel("wind speed")
zlabel("damping ratio")
title("different AOA cases (no model updating)")
% 更新 current_plot
current_plot = current_plot + 1;

%% new figure
total_plots = 2; % 或任何你需要的子图数量
current_plot = 1;
num_figs_in_row = [];
figWidthFactor = 1.5;
figPosition = [100, 100];
newfigure = true;
holdon = false;
firstfigure = true;
%% different VIV case and use sensor 1 model updating select VIV exciting case
VIV_sel  = VIV_sel;
plot_sel = find(ismember(VIV_sels, VIV_sel));
plotdata_sel = plotdata_modelupdate_sel_excite;
for k1 = 1:length(plot_sel)
    x = plotdata_sel(plot_sel(k1)).amp_filter;
    y = plotdata_sel(plot_sel(k1)).U_sel_1;
    z = plotdata_sel(plot_sel(k1)).zeta_filter;
    color_points = colors{k1};
    create_subplot(@scatter3, total_plots, current_plot, {x, y, z,[],color_points}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
    hold on
end
zlim([-0.01 0])
current_plot = current_plot + 1;
xlabel("amplitude")
ylabel("wind speed")
zlabel("damping ratio")
title("different VIV cases (model updating)")

%% different VIV case and use sensor 1 no model updating select VIV exciting case
VIV_sel  = VIV_sel;
plot_sel = find(ismember(VIV_sels, VIV_sel));
plotdata_sel = plotdata_nomodelupdate_sel_excite;
for k1 = 1:length(plot_sel)
    x = plotdata_sel(plot_sel(k1)).amp_filter;
    y = plotdata_sel(plot_sel(k1)).U_sel_1;
    z = plotdata_sel(plot_sel(k1)).zeta_filter;
    color_points = colors{k1};
    create_subplot(@scatter3, total_plots, current_plot, {x, y, z,[],color_points}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
    hold on
end
zlim([-0.01 0])
current_plot = current_plot + 1;
xlabel("amplitude")
ylabel("wind speed")
zlabel("damping ratio")
title("different VIV cases (no model updating)")


%% new figure
total_plots = 2; % 或任何你需要的子图数量
current_plot = 1;
num_figs_in_row = [];
figWidthFactor = 1.5;
figPosition = [100, 100];
newfigure = true;
holdon = false;
firstfigure = true;
%% different AOA case and use sensor 1 model updating select VIV exciting case
VIV_sel = VIV_sel;
plot_sel = find(ismember(VIV_sels, VIV_sel));

% 初始化图例标签数组
% legend_labels = cell(1, length(categories) - 1);
plotdata_sel = plotdata_modelupdate_sel_excite;
for k1 = 1:length(plot_sel)
    x = plotdata_sel(plot_sel(k1)).amp_filter;
    y = plotdata_sel(plot_sel(k1)).U_sel_1;
    z = plotdata_sel(plot_sel(k1)).zeta_filter;
    color_points = plotdata_sel(plot_sel(k1)).AoA_sel_1;

    % 分类 color_points
    categories = -6:1:10;
    [~, ~, categorized_color_points] = histcounts(color_points, categories);

    for k2 = 1:length(categories) - 1
        % 选择当前类别的数据点
        idx = categorized_color_points == k2;
        x_cat = x(idx);
        y_cat = y(idx);
        z_cat = z(idx);

        % 选择颜色
        color_index = mod(k2, length(colors)) + 1;
        plot_color = colors{color_index};

        % 绘制散点图并为每个类别添加图例
        create_subplot(@scatter3, total_plots, current_plot, {x_cat, y_cat, z_cat, [], plot_color}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;


        hold on;

        % 更新图例标签
        legend_labels{k2} = ['AOA ' num2str(categories(k2))];
    end
end

% 添加图例
legend(legend_labels);
zlim([-0.01 0])
xlabel("amplitude")
ylabel("wind speed")
zlabel("damping ratio")
title("different AOA cases (model updating)")
% 更新 current_plot
current_plot = current_plot + 1;

%% different AOA case and use sensor 1 no model updating select VIV exciting case
VIV_sel = VIV_sel;
plot_sel = find(ismember(VIV_sels, VIV_sel));

% 初始化图例标签数组
% legend_labels = cell(1, length(categories) - 1);
plotdata_sel = plotdata_nomodelupdate_sel_excite;
for k1 = 1:length(plot_sel)
    x = plotdata_sel(plot_sel(k1)).amp_filter;
    y = plotdata_sel(plot_sel(k1)).U_sel_1;
    z = plotdata_sel(plot_sel(k1)).zeta_filter;
    color_points = plotdata_sel(plot_sel(k1)).AoA_sel_1;

    % 分类 color_points
    categories = -6:1:10;
    [~, ~, categorized_color_points] = histcounts(color_points, categories);

    for k2 = 1:length(categories) - 1
        % 选择当前类别的数据点
        idx = categorized_color_points == k2;
        x_cat = x(idx);
        y_cat = y(idx);
        z_cat = z(idx);

        % 选择颜色
        color_index = mod(k2, length(colors)) + 1;
        plot_color = colors{color_index};

        % 绘制散点图并为每个类别添加图例
        create_subplot(@scatter3, total_plots, current_plot, {x_cat, y_cat, z_cat, [], plot_color}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;


        hold on;

        % 更新图例标签
        legend_labels{k2} = ['AOA ' num2str(categories(k2))];
    end
end

% 添加图例
legend(legend_labels);
zlim([-0.01 0])
xlabel("amplitude")
ylabel("wind speed")
zlabel("damping ratio")
title("different AOA cases (no model updating)")
% 更新 current_plot
current_plot = current_plot + 1;


%% new figure
total_plots = 2; % 或任何你需要的子图数量
current_plot = 1;
num_figs_in_row = [];
figWidthFactor = 1.5;
figPosition = [100, 100];
newfigure = true;
holdon = false;
firstfigure = true;
%% excite and drop plot model updating
VIV_sel = VIV_sel;
plot_sel = find(ismember(VIV_sels, VIV_sel));
plotdata_sel = plotdata_modelupdate_sel_excite;
for k1 = 1:length(plot_sel)
    x = plotdata_sel(plot_sel(k1)).amp_filter;
    y = plotdata_sel(plot_sel(k1)).U_sel_1;
    z = plotdata_sel(plot_sel(k1)).zeta_filter;
    color_points = colors{1};
    create_subplot(@scatter3, total_plots, current_plot, {x, y, z,[],color_points}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
    hold on
end
plotdata_sel = plotdata_modelupdate_sel_drop;
for k1 = 1:length(plot_sel)
    x = plotdata_sel(plot_sel(k1)).amp_filter;
    y = plotdata_sel(plot_sel(k1)).U_sel_1;
    z = plotdata_sel(plot_sel(k1)).zeta_filter;
    color_points = colors{2};
    create_subplot(@scatter3, total_plots, current_plot, {x, y, z,[],color_points}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
    hold on
end
zlim([-0.01 0])
xlabel("amplitude")
ylabel("wind speed")
zlabel("damping ratio")
title("different VIV cases (model updating)")


% 更新 current_plot
current_plot = current_plot + 1;

%% excite and drop plot no model updating
VIV_sel = VIV_sel;
plot_sel = find(ismember(VIV_sels, VIV_sel));
plotdata_sel = plotdata_nomodelupdate_sel_excite;
for k1 = 1:length(plot_sel)
    x = plotdata_sel(plot_sel(k1)).amp_filter;
    y = plotdata_sel(plot_sel(k1)).U_sel_1;
    z = plotdata_sel(plot_sel(k1)).zeta_filter;
    color_points = colors{1};
    create_subplot(@scatter3, total_plots, current_plot, {x, y, z,[],color_points}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
    hold on
end
plotdata_sel = plotdata_nomodelupdate_sel_drop;
for k1 = 1:length(plot_sel)
    x = plotdata_sel(plot_sel(k1)).amp_filter;
    y = plotdata_sel(plot_sel(k1)).U_sel_1;
    z = plotdata_sel(plot_sel(k1)).zeta_filter;
    color_points = colors{2};
    create_subplot(@scatter3, total_plots, current_plot, {x, y, z,[],color_points}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
    hold on
end
zlim([-0.01 0])
xlabel("amplitude")
ylabel("wind speed")
zlabel("damping ratio")
title("different VIV cases (no model updating)")


%% new figure
total_plots = 2; % 或任何你需要的子图数量
current_plot = 1;
num_figs_in_row = [];
figWidthFactor = 1.5;
figPosition = [100, 100];
newfigure = true;
holdon = false;
firstfigure = true;
%% excite and drop plot model updating (without plotting the damping ratio there is almost no differenct between results between model updating and no model updating) amp vs wind speed vs AOA
VIV_sel = VIV_sel;
plot_sel = find(ismember(VIV_sels, VIV_sel));
plotdata_sel = plotdata_modelupdate_sel_excite;
for k1 = 1:length(plot_sel)
    x = plotdata_sel(plot_sel(k1)).amp_filter;
    y = plotdata_sel(plot_sel(k1)).U_sel_1;
    z = plotdata_sel(plot_sel(k1)).AoA_sel_1;
    color_points = colors{1};
    create_subplot(@scatter3, total_plots, current_plot, {x, y, z,[],color_points}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
    hold on
end
plotdata_sel = plotdata_modelupdate_sel_drop;
for k1 = 1:length(plot_sel)
    x = plotdata_sel(plot_sel(k1)).amp_filter;
    y = plotdata_sel(plot_sel(k1)).U_sel_1;
    z = plotdata_sel(plot_sel(k1)).AoA_sel_1;
    color_points = colors{2};
    create_subplot(@scatter3, total_plots, current_plot, {x, y, z,[],color_points}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
    hold on
end
zlim([-3, 3])
xlabel("amplitude")
ylabel("wind speed")
zlabel("AoA")
title("different VIV cases (model updating)")


% 更新 current_plot
current_plot = current_plot + 1;

%% excite and drop plot model updating (without plotting the damping ratio there is almost no differenct between results between model updating and no model updating) amp vs wind speed vs beta
VIV_sel = VIV_sel;
plot_sel = find(ismember(VIV_sels, VIV_sel));
plotdata_sel = plotdata_modelupdate_sel_excite;
for k1 = 1:length(plot_sel)
    x = plotdata_sel(plot_sel(k1)).amp_filter;
    y = plotdata_sel(plot_sel(k1)).U_sel_1;
    z = plotdata_sel(plot_sel(k1)).beta_deg_mean_UA1a;
    color_points = colors{1};
    create_subplot(@scatter3, total_plots, current_plot, {x, y, z,[],color_points}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
    hold on
end
plotdata_sel = plotdata_modelupdate_sel_drop;
for k1 = 1:length(plot_sel)
    x = plotdata_sel(plot_sel(k1)).amp_filter;
    y = plotdata_sel(plot_sel(k1)).U_sel_1;
    z = plotdata_sel(plot_sel(k1)).beta_deg_mean_UA1a;
    color_points = colors{2};
    create_subplot(@scatter3, total_plots, current_plot, {x, y, z,[],color_points}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
    hold on
end
% zlim([-3, 3])
xlabel("amplitude")
ylabel("wind speed")
zlabel("wind direction")
title("different VIV cases (model updating)")


% 更新 current_plot
current_plot = current_plot + 1;



holdon = true;

return
%% test
total_plots = 8; % 或任何你需要的子图数量
current_plot = 1;
num_figs_in_row = [];
figWidthFactor = 1.5;
figPosition = [100, 100];
newfigure = true;
holdon = false;
firstfigure = true;
% figure
VIV_sel = 11;
plot_sel = find(ismember(VIV_sels, VIV_sel));
plotdata_sel = plotdata_modelupdate;
% 初始化图例标签数组
% legend_labels = cell(1, length(categories) - 1);

for k1 = 1:length(plot_sel)
    x = plotdata_sel(plot_sel(k1)).amp_filter;
    y = plotdata_sel(plot_sel(k1)).U_sel_1;
    z = plotdata_sel(plot_sel(k1)).zeta_filter;
    color_points = plotdata_sel(plot_sel(k1)).AoA_sel_1;

    % 分类 color_points
    categories = -6:1:10;
    [~, ~, categorized_color_points] = histcounts(color_points, categories);

    for k2 = 1:length(categories) - 1
        % 选择当前类别的数据点
        idx = categorized_color_points == k2;
        x_cat = x(idx);
        y_cat = y(idx);
        z_cat = z(idx);

        % 选择颜色
        color_index = mod(k2, length(colors)) + 1;
        plot_color = colors{color_index};

        % 绘制散点图并为每个类别添加图例
        create_subplot(@scatter3, total_plots, current_plot, {x_cat, y_cat, z_cat, [], plot_color}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;


        hold on;

        % 更新图例标签
        legend_labels{k2} = ['AOA ' num2str(categories(k2))];
    end
end
current_plot=current_plot+1
% 添加图例
legend(legend_labels);
zlim([-0.02 0.02])

hold off
% test wind


% VIV_sel  = 4;
plot_sel = find(VIV_sels == VIV_sel);

x = plotdata_sel(k1).t_cycle_mean_temp;
y = plotdata_sel(k1).beta_deg_mean_UA1a;
create_subplot(@scatter, total_plots, current_plot, {x, y}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
hold on
y = plotdata_sel(k1).beta_deg_mean_UA1b;
create_subplot(@scatter, total_plots, current_plot, {x, y}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
y = plotdata_sel(k1).beta_deg_mean_UA2a;
create_subplot(@scatter, total_plots, current_plot, {x, y}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
y = plotdata_sel(k1).beta_deg_mean_UA2b;
create_subplot(@scatter, total_plots, current_plot, {x, y}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', false);firstfigure = false;
y = plotdata_sel(k1).beta_deg_mean_UA3a;
create_subplot(@scatter, total_plots, current_plot, {x, y}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', false);firstfigure = false;
y = plotdata_sel(k1).beta_deg_mean_UA3b;
create_subplot(@scatter, total_plots, current_plot, {x, y}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', false);firstfigure = false;
xlabel('Time')
ylabel('Wind direction (deg)')
legend('UA1a','UA1b','UA2a','UA2b','UA3a','UA3b')
current_plot = current_plot + 1;

% test wind
% VIV_sel  = 4;
plot_sel = find(VIV_sels == VIV_sel);

x = plotdata_sel(k1).t_cycle_mean_temp;
y = plotdata_sel(k1).U_sel_1;
create_subplot(@scatter, total_plots, current_plot, {x, y}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
hold on
y = plotdata_sel(k1).U_sel_2;
create_subplot(@scatter, total_plots, current_plot, {x, y}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
y = plotdata_sel(k1).U_sel_3;
create_subplot(@scatter, total_plots, current_plot, {x, y}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;

xlabel('Time')
ylabel('Wind speed location (m/s)')
legend('location 1','location 2','location 3')

current_plot = current_plot + 1;


% test wind
% VIV_sel  = 4;
plot_sel = find(VIV_sels == VIV_sel);

x = plotdata_sel(k1).t_cycle_mean_temp;
y = plotdata_sel(k1).AoA_sel_1;
create_subplot(@scatter, total_plots, current_plot, {x, y}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
hold on
y = plotdata_sel(k1).AoA_sel_2;
create_subplot(@scatter, total_plots, current_plot, {x, y}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
y = plotdata_sel(k1).AoA_sel_3;
create_subplot(@scatter, total_plots, current_plot, {x, y}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;

xlabel('Time')
ylabel('AOA  (deg)')
legend('location 1','location 2','location 3')

current_plot = current_plot + 1;

% % time history
% start_time = vivTable.startDate(VIV_sel);
% end_time = vivTable.endDate(VIV_sel);
% [Acc_Data] = read_acceleration_data(start_time, end_time, input_data.acc_dir);
% [Wind_Data] = read_wind_data(start_time, end_time, input_data.wind_dir);
% Acc_Data = Acc_Data.mergedData;
% 
% t_acc = Acc_Data.Time;
% AC2_1 = Acc_Data.AC2_1;
% AC2_3 = Acc_Data.AC2_3;
% AC3_1 = Acc_Data.AC3_1;
% AC3_3 = Acc_Data.AC3_3;
% AC4_1 = Acc_Data.AC4_1;
% AC4_3 = Acc_Data.AC4_3; 
% 
% 
% create_subplot(@plot, total_plots, current_plot, {t_acc, AC2_1}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
% hold on
% create_subplot(@plot, total_plots, current_plot, {t_acc, AC2_3}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
% create_subplot(@plot, total_plots, current_plot, {t_acc, AC3_1}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
% create_subplot(@plot, total_plots, current_plot, {t_acc, AC3_3}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
% create_subplot(@plot, total_plots, current_plot, {t_acc, AC4_1}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
% create_subplot(@plot, total_plots, current_plot, {t_acc, AC4_3}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
% current_plot=current_plot+1;
% 
% t_wind = Wind_Data.resultsTable_UA1.Time_Start;
% u1 = Wind_Data.resultsTable_UA1.U;
% u2 = Wind_Data.resultsTable_UA2.U;
% u3 = Wind_Data.resultsTable_UA3.U;
% u4 = Wind_Data.resultsTable_UA4.U;
% u5 = Wind_Data.resultsTable_UA5.U;
% u6 = Wind_Data.resultsTable_UA6.U;
% 
% 
% create_subplot(@plot, total_plots, current_plot, {t_wind, u1}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
% hold on
% create_subplot(@plot, total_plots, current_plot, {t_wind, u2}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
% create_subplot(@plot, total_plots, current_plot, {t_wind, u3}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
% create_subplot(@plot, total_plots, current_plot, {t_wind, u4}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
% create_subplot(@plot, total_plots, current_plot, {t_wind, u5}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
% create_subplot(@plot, total_plots, current_plot, {t_wind, u6}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', firstfigure, 'holdon', holdon);firstfigure = false;
% 


