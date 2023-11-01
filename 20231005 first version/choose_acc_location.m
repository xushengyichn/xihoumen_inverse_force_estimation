%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: ShengyiXu xushengyichn@outlook.com
%Date: 2023-11-01 11:06:32
%LastEditors: ShengyiXu xushengyichn@outlook.com
%LastEditTime: 2023-11-01 11:09:14
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

[figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
plot(node_loc,mode_deck)
hold on
scatter(loc_acc,mode_deck_loc)



n=4;
[result] = viv2013(n, OFF);
startDate_global = result.startDate;
endDate_global = result.endDate+hours(1);
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