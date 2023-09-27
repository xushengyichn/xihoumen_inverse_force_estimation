%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: ShengyiXu xushengyichn@outlook.com
%Date: 2023-09-27 12:52:24
%LastEditors: ShengyiXu xushengyichn@outlook.com
%LastEditTime: 2023-09-27 13:01:11
%FilePath: \xihoumen_inverse_force_estimation\simulation signal\Main.m
%Description: 使用模拟荷载信号，测试反算气动力的准确性，以及参数选择的合理性
%
%Copyright (c) 2023 by ${git_name_email}, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear; close all
addpath(genpath("F:\git\ssm_tools"))
addpath(genpath("D:\Users\xushe\Documents\GitHub\ssm_tools"))
addpath(genpath("D:\Users\xushe\Documents\GitHub\Function_shengyi_package"))
addpath(genpath("F:\git\Function_shengyi_package"))
addpath(genpath("C:\Users\xushe\OneDrive\NAS云同步\Drive\0博士研究生\3大论文\研究内容\研究内容 3：气动力模型及参数识别；\反算气动力"))
addpath(genpath("D:\OneDrive\NAS云同步\Drive\0博士研究生\3大论文\研究内容\研究内容 3：气动力模型及参数识别；\反算气动力"))

subStreamNumberDefault = 2132;
run('InitScript.m');

%% 0 绘图参数
fig_bool = ON;
num_figs_in_row = 10; %每一行显示几个图
figPos = figPosSmall; %图的大小，参数基于InitScript.m中的设置
%设置图片间隔
gap_between_images = [0, 0];
figureIdx = 0;

f1=1;
dt=1/50;
t=0:dt:200;
lambda = 0.00001; sigma_p_2 = 100000000; sigma_p = sqrt(sigma_p_2);  

[P] = quasiperiod_force(lambda, sigma_p, f1, dt, t);
p = P;

rms(p)



if fig_bool == ON
[figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
plot(t,p)
xlabel('t (s)')
ylabel('Modal force')



end