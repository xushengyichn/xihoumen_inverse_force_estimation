
clc;clear;close all


addpath(genpath("F:\git\ssm_tools"))
addpath(genpath("D:\Users\xushe\Documents\GitHub\ssm_tools"))
addpath(genpath("D:\Users\xushe\Documents\GitHub\Function_shengyi_package"))
addpath(genpath("F:\git\Function_shengyi_package"))
addpath(genpath("D:\OneDrive\NAS云同步\Drive\0博士研究生\3大论文\研究内容\研究内容 3：气动力模型及参数识别；\反算气动力"))
addpath(genpath("C:\Users\xushe\OneDrive\NAS云同步\Drive\0博士研究生\3大论文\研究内容\研究内容 3：气动力模型及参数识别；\反算气动力"))
addpath(genpath("C:\Users\xushe\OneDrive\NAS云同步\Drive\0博士研究生\3大论文\研究内容\研究内容 3：气动力模型及参数识别；\反算气动力"))
addpath(genpath("G:\2013\allVIB（移动）"))
subStreamNumberDefault = 2132;
run('InitScript.m');
%% 0 绘图参数
fig_bool = ON;
num_figs_in_row = 5; %每一行显示几个图
figPos = figPosSmall; %图的大小，参数基于InitScript.m中的设置
%设置图片间隔
gap_between_images = [0, 0];
figureIdx = 0;

% 定义一个简单的测试信号
fs = 1000;  % 采样频率
t = 0:1/fs:1-1/fs;  % 时间向量
f1 = 50;  % 信号的频率
f2 = 60;
data1 = 2*sin(2*pi*f1*t);  % 生成信号
data2 = sin(2*pi*f2*t);  % 生成信号
data = data1+data2;
% 定义需要保留的频率范围
f_keep = [45 55];  % 保留 45 到 55 Hz 之间的频率

% 使用fft_filter函数
data_filtered = fft_filter(fs, data, f_keep);

[f, magnitude] =fft_transform(fs,data_filtered);

% 为了验证幅度是否正确, 计算原始信号和滤波信号的RMS(均方根)值
rms_original = sqrt(mean(data.^2));
rms_filtered = sqrt(mean(data_filtered.^2));

% 显示结果
fprintf('Original RMS: %.2f\n', rms_original);
fprintf('Filtered RMS: %.2f\n', rms_filtered);

% 可视化结果
figure;
subplot(2,1,1);
plot(t, data);
title('Original Signal');
subplot(2,1,2);
plot(t, data_filtered);
title('Filtered Signal');

figure
plot(f,magnitude)
title('fft');