%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: ShengyiXu xushengyichn@outlook.com
%Date: 2023-11-01 11:06:32
%LastEditors: ShengyiXu xushengyichn@outlook.com
%LastEditTime: 2023-11-01 22:08:51
%FilePath: \Exercises-for-Techniques-for-estimation-in-dynamics-systemsf:\git\xihoumen_inverse_force_estimation\20231005 first version\choose_acc_location.m
%Description: 由于加速度传感器的布置位置可能会有一定的偏移，现在希望通过调整加速度计的位置，使不同位置的加速度计满足振型的规律
%
%Copyright (c) 2023 by ${git_name_email}, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



clc; clear; close all;
computer_name = getenv('COMPUTERNAME');
if strcmp(computer_name,'SHENGYI_HP')
    addpath(genpath("F:\git\ssm_tools\"))
    addpath(genpath("F:\git\Function_shengyi_package\"))
    addpath(genpath("F:\git\xihoumen_inverse_force_estimation\FEM_model\"))
    addpath(genpath("F:\git\HHT-Tutorial\"))
elseif strcmp(computer_name,'mac')
    addpath(genpath("/Users/xushengyi/Documents/GitHub/Function_shengyi_package"))
    addpath(genpath("/Users/xushengyi/Documents/GitHub/ssm_tools_sy"))
    addpath(genpath("/Users/xushengyi/Documents/GitHub/xihoumen_inverse_force_estimation/FEM_model"))
    addpath(genpath("/Users/xushengyi/Documents/GitHub/HHT-Tutorial"))
elseif strcmp(computer_name,'ROG-SHENGYI')
    addpath(genpath("D:\Users\xushe\Documents\GitHub\Function_shengyi_package"))
    addpath(genpath("D:\Users\xushe\Documents\GitHub\ssm_tools"))
    addpath(genpath("D:\git\xihoumen_inverse_force_estimation\FEM_model"))
    addpath(genpath("D:\Users\xushe\Documents\GitHub\xihoumen_data_extract"))
    addpath(genpath("D:\Users\xushe\Documents\GitHub\HHT-Tutorial\"))
elseif strcmp(computer_name,'NTNU08916')
    addpath(genpath("C:\Users\shengyix\Documents\GitHub\Function_shengyi_package"))
    addpath(genpath("C:\Users\shengyix\Documents\GitHub\ssm_tools_sy"))
    addpath(genpath("C:\Users\shengyix\Documents\GitHub\xihoumen_inverse_force_estimation\FEM_model"))
    addpath(genpath("C:\Users\shengyix\Documents\GitHub\xihoumen_data_extract"))
else
    error("Please add path first.")
end


subStreamNumberDefault = 2132;

run("InitScript.m")

%% 0 绘图参数
fig_bool = ON;
num_figs_in_row = 4; %每一行显示几个图
figPos = figPosSmall; %图的大小，参数基于InitScript.m中的设置
%设置图片间隔
gap_between_images = [0, 0];
figureIdx = 0;

modesel = 23;
nmodes = length(modesel); ns = nmodes * 2;
Result = ImportMK(nmodes, 'KMatrix.matrix', 'MMatrix.matrix', 'nodeondeck.txt', 'KMatrix.mapping', 'nodegap.txt', 'modesel', modesel,'showtext',true);

loc_acc= [578+1650/4;578+1650/2;578+1650/4*3];

node_loc = Result.node_loc;
nodeondeck = Result.nodeondeck;
eig_vec = Result.eig_vec;
nodegap = Result.nodegap;
mode_vec = Result.mode_vec;
KMmapping = Result.Mapping;
mode_deck = Result.mode_deck;

acc_node = FindNodewithLocation(loc_acc, node_loc, nodeondeck);

acc_node_left = acc_node.left;
acc_node_right = acc_node.right;



Nodeondeck_Info_left = Cal_Nodeondeck_Info(acc_node_left, KMmapping, nodegap, mode_vec,'showtext',true);
Nodeondeck_Info_right = Cal_Nodeondeck_Info(acc_node_right, KMmapping, nodegap, mode_vec,'showtext',true);


weight = acc_node.weights;
mode_deck_loc_two = [Nodeondeck_Info_left.mode_deck,Nodeondeck_Info_right.mode_deck];
mode_deck_loc = diag(mode_deck_loc_two *weight);


mode_deck_norm = mode_deck/mode_deck_loc(2);

mode_deck_loc_norm  = mode_deck_loc/mode_deck_loc(2);



n=18;
[result] = viv2013(n, OFF);
startDate_global = result.startDate;
endDate_global = result.endDate;
input.start_time = startDate_global;
input.end_time = endDate_global;


computer_name = getenv('COMPUTERNAME');
if strcmp(computer_name,'SHENGYI_HP')
    input.acc_dir = "F:\test\result";
    input.wind_dir = "F:\test\result_wind_10min";
elseif strcmp(computer_name,'mac')
    input.acc_dir = "/Users/xushengyi/Documents/xihoumendata/acc";
    input.wind_dir = "/Users/xushengyi/Documents/xihoumendata/wind";
elseif strcmp(computer_name,'ROG-SHENGYI')
    input.acc_dir = "D:\xihoumendata\acc";
    input.wind_dir = "D:\xihoumendata\wind";
elseif strcmp(computer_name,'ketizu')
    input.acc_dir = "Z:\Drive\Backup\SHENGYI_HP\F\test\result";
    input.wind_dir = "Z:\Drive\Backup\SHENGYI_HP\F\test\result_wind_10min";
elseif strcmp(computer_name,'NTNU08916')
     input.acc_dir = "C:\Users\shengyix\Documents\xihoumendata\acc";
    input.wind_dir = "C:\Users\shengyix\Documents\xihoumendata\wind";
else
    error("Please add data folder first.")
end

start_time = input.start_time;
end_time = input.end_time;

acc_dir = input.acc_dir;

[Acc_Data] = read_acceleration_data(start_time, end_time, acc_dir);


f_keep =  0.33 * [0.9, 1.1];


timeDifferences = diff(Acc_Data.mergedData.Time); % Returns a duration array
dt = seconds(timeDifferences(1)); % Converts the first duration value to seconds and assigns to dt
fs = 1 / dt;
Acc_Data.mergedData.AC2_1 = fft_filter(fs, Acc_Data.mergedData.AC2_1, f_keep)';
Acc_Data.mergedData.AC2_2 = fft_filter(fs, Acc_Data.mergedData.AC2_2, f_keep)'; 
Acc_Data.mergedData.AC2_3 = fft_filter(fs, Acc_Data.mergedData.AC2_3, f_keep)';
Acc_Data.mergedData.AC3_1 = fft_filter(fs, Acc_Data.mergedData.AC3_1, f_keep)';
Acc_Data.mergedData.AC3_2 = fft_filter(fs, Acc_Data.mergedData.AC3_2, f_keep)';
Acc_Data.mergedData.AC3_3 = fft_filter(fs, Acc_Data.mergedData.AC3_3, f_keep)';
Acc_Data.mergedData.AC4_1 = fft_filter(fs, Acc_Data.mergedData.AC4_1, f_keep)';
Acc_Data.mergedData.AC4_2 = fft_filter(fs, Acc_Data.mergedData.AC4_2, f_keep)';
Acc_Data.mergedData.AC4_3 = fft_filter(fs, Acc_Data.mergedData.AC4_3, f_keep)';

yn(1, :) = Acc_Data.mergedData.AC2_1 / 1000 * 9.8;
yn(2, :) = Acc_Data.mergedData.AC2_3 / 1000 * 9.8;
yn(3, :) = Acc_Data.mergedData.AC3_1 / 1000 * 9.8;
yn(4, :) = Acc_Data.mergedData.AC3_3 / 1000 * 9.8;
yn(5, :) = Acc_Data.mergedData.AC4_1 / 1000 * 9.8;
yn(6, :) = Acc_Data.mergedData.AC4_3 / 1000 * 9.8;

AC2 =( yn(1, :) +yn(2, :) )/2;
AC3 =( yn(3, :) +yn(4, :) )/2;
AC4 =( yn(5, :) +yn(6, :) )/2;

AC2_rms = rms(AC2-mean(AC2));
AC3_rms = rms(AC3-mean(AC3));
AC4_rms = rms(AC4-mean(AC4));

AC2_norm = AC2_rms/AC3_rms;
AC3_norm = AC3_rms/AC3_rms;
AC4_norm = AC4_rms/AC3_rms;

ACplot = [AC2_norm,AC3_norm,AC4_norm];


[figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
plot(node_loc,mode_deck_norm)
hold on
scatter(loc_acc,mode_deck_loc_norm)
scatter(loc_acc,ACplot)



[f, magnitude] = fft_transform(50,yn(3, :),'showtext',true);

[~, index] = max(magnitude);
peak_freq = f(index);

title("Frequency is "+num2str(peak_freq)+".Date:"+string(startDate_global) )