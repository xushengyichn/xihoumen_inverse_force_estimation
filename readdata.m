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

%% vibration
% vibac2_name=["2013-02-06 00-vibac2.txt";"2013-02-06 01-vibac2.txt";"2013-02-06 02-vibac2.txt"];
% vibac3_name=["2013-02-06 01-VIBac3.txt";"2013-02-06 01-VIBac3.txt";"2013-02-06 02-VIBac3.txt"];
% vibac4_name=["2013-02-06 00-VIBac4.txt";"2013-02-06 01-VIBac4.txt";"2013-02-06 02-VIBac4.txt"];
% startDate = datetime(2013,2,6,0,0,0);
% endDate = datetime(2013,2,6,3,0,0);
% 创建开始和结束日期的数组
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

%% Readdata

% for k1 = 2
for k1 = 1:22
result = viv2013(k1);


vibac2_name=result.vibac2_name;
vibac3_name=result.vibac3_name;
vibac4_name=result.vibac4_name;
startDate = result.startDate;
endDate = result.endDate;


vibac2_data=read_vib_data(vibac2_name);
vibac3_data=read_vib_data(vibac3_name);
vibac4_data=read_vib_data(vibac4_name);

vibac2_data(:,2:4)=vibac2_data(:,2:4)/1000*9.8;%第一列是是时间，不要修改
vibac3_data(:,2:4)=vibac3_data(:,2:4)/1000*9.8;%第一列是是时间，不要修改
vibac4_data(:,2:4)=vibac4_data(:,2:4)/1000*9.8;%第一列是是时间，不要修改

maxvalue = max(max(abs([vibac2_data(:,2:4);vibac2_data(:,2:4);vibac4_data(:,2:4)])));

dt = 0.02;
t = 0:dt:dt*(length(vibac2_data(:,1))-1);

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

winds1_trim = winds1_data(winds1_data.DateTime >= startDate & winds1_data.DateTime <= endDate, :);
winds3_trim = winds3_data(winds3_data.DateTime >= startDate & winds3_data.DateTime <= endDate, :);
winds4_trim = winds4_data(winds4_data.DateTime >= startDate & winds4_data.DateTime <= endDate, :);

%% fftanalysis

[f, magnitude] = fft_transform(1/dt,vibac3_data(:,3));
[a,b]= max(magnitude);
f_dom(k1)=f(b);

%% filter
if 0
    % f_dom(k1)=0.32
    for k2 = 1:3
        vibac2_data(:,k2+1) = fft_filter(1/dt, vibac2_data(:,k2+1), [f_dom(k1)*0.95 f_dom(k1)*1.05]);
        vibac3_data(:,k2+1) = fft_filter(1/dt, vibac3_data(:,k2+1), [f_dom(k1)*0.95 f_dom(k1)*1.05]);
        vibac4_data(:,k2+1) = fft_filter(1/dt, vibac4_data(:,k2+1), [f_dom(k1)*0.95 f_dom(k1)*1.05]);
    end
end
%% plot

folderName = 'images';  % 为目标文件夹指定名称

% 确保文件夹存在
if ~exist(folderName, 'dir')
    mkdir(folderName);  % 如果不存在，则创建文件夹
end


[figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
hold on


plot(t,vibac2_data(:,3), 'Color', 'k', 'LineStyle', '--')
plot(t,vibac2_data(:,4), 'Color', 'r', 'LineStyle', '--')
plot(t,vibac2_data(:,2), 'Color', 'b', 'LineStyle', '--')
set(gca, 'FontSize', 12)
xlabel('Time')
ylabel('acc')
% legend('filtered', 'true', 'Location', 'northwest')
title(['vib2']);
% xlim([0, 3])
ylim([-maxvalue , maxvalue ])

set(gcf, 'unit', 'centimeters', 'position', [35 10 7.8 6]);
set(gca,'unit', 'centimeters', 'position', [1.2, 1.5, 6, 4]);
xlh = get(gca,'xlabel');
xlv = get(gca,'ylabel');
xlh.Position(2) = xlh.Position(2) + 0.007;
xlv.Position(1) = xlv.Position(1) + 0.005;
ax = gca; % current axes
ax.XRuler.TickLabelGapOffset = -3; % negative numbers move the ticklabels down (positive -> up)
ax.YRuler.TickLabelGapOffset = -1.5; % negative numbers move the ticklabels right (negative -> left)
fileName = char(startDate)+"_acc2.png";  % 为文件指定名称
fullFileName = fullfile(folderName, fileName);  % 创建完整的文件路径
print(gcf, fullFileName,'-r300','-dpng');
% exportgraphics(gcf, 'test2.png','Resolution',300);
% copygraphics(gcf)

[figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
hold on


plot(t,vibac3_data(:,3), 'Color', 'k', 'LineStyle', '--')
plot(t,vibac3_data(:,4), 'Color', 'r', 'LineStyle', '--')
plot(t,vibac3_data(:,2), 'Color', 'b', 'LineStyle', '--')
set(gca, 'FontSize', 12)
xlabel('Time')
ylabel('acc')
% legend('filtered', 'true', 'Location', 'northwest')
title(['vib3']);
% xlim([0, 3])
ylim([-maxvalue , maxvalue ])

set(gcf, 'unit', 'centimeters', 'position', [35 10 7.8 6]);
set(gca,'unit', 'centimeters', 'position', [1.2, 1.5, 6, 4]);
xlh = get(gca,'xlabel');
xlv = get(gca,'ylabel');
xlh.Position(2) = xlh.Position(2) + 0.007;
xlv.Position(1) = xlv.Position(1) + 0.005;
ax = gca; % current axes
ax.XRuler.TickLabelGapOffset = -3; % negative numbers move the ticklabels down (positive -> up)
ax.YRuler.TickLabelGapOffset = -1.5; % negative numbers move the ticklabels right (negative -> left)
fileName = char(startDate)+"_acc3.png";  % 为文件指定名称
fullFileName = fullfile(folderName, fileName);  % 创建完整的文件路径
print(gcf, fullFileName,'-r300','-dpng');
% exportgraphics(gcf, 'test2.png','Resolution',300);
% copygraphics(gcf)


[figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
hold on

plot(t,vibac4_data(:,3), 'Color', 'k', 'LineStyle', '--')
plot(t,vibac4_data(:,4), 'Color', 'r', 'LineStyle', '--')
plot(t,vibac4_data(:,2), 'Color', 'b', 'LineStyle', '--')

set(gca, 'FontSize', 12)
xlabel('Time')
ylabel('acc')
% legend('filtered', 'true', 'Location', 'northwest')
title(['vib4']);
% xlim([0, 3])
ylim([-maxvalue , maxvalue ])

set(gcf, 'unit', 'centimeters', 'position', [35 10 7.8 6]);
set(gca,'unit', 'centimeters', 'position', [1.2, 1.5, 6, 4]);
xlh = get(gca,'xlabel');
xlv = get(gca,'ylabel');
xlh.Position(2) = xlh.Position(2) + 0.007;
xlv.Position(1) = xlv.Position(1) + 0.005;
ax = gca; % current axes
ax.XRuler.TickLabelGapOffset = -3; % negative numbers move the ticklabels down (positive -> up)
ax.YRuler.TickLabelGapOffset = -1.5; % negative numbers move the ticklabels right (negative -> left)
fileName = char(startDate)+"_acc4.png";  % 为文件指定名称
fullFileName = fullfile(folderName, fileName);  % 创建完整的文件路径
print(gcf, fullFileName,'-r300','-dpng');
% exportgraphics(gcf, 'test2.png','Resolution',300);
% copygraphics(gcf)

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
ylim([0, 15])

set(gcf, 'unit', 'centimeters', 'position', [35 10 7.8 6]);
set(gca,'unit', 'centimeters', 'position', [1.2, 1.5, 6, 4]);
xlh = get(gca,'xlabel');
xlv = get(gca,'ylabel');
xlh.Position(2) = xlh.Position(2) + 0.007;
xlv.Position(1) = xlv.Position(1) + 0.005;
ax = gca; % current axes
ax.XRuler.TickLabelGapOffset = -3; % negative numbers move the ticklabels down (positive -> up)
ax.YRuler.TickLabelGapOffset = -1.5; % negative numbers move the ticklabels right (negative -> left)
fileName = char(startDate)+"_wind.png";  % 为文件指定名称
fullFileName = fullfile(folderName, fileName);  % 创建完整的文件路径
print(gcf, fullFileName,'-r300','-dpng');
% exportgraphics(gcf, 'test2.png','Resolution',300);
% copygraphics(gcf)

[figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
plot(f, magnitude)
set(gca, 'FontSize', 12)
xlabel('Frequency (Hz)')
ylabel('amplitude')
title("Dominant frequency:"+num2str(f_dom(k1)))
xlim([0, 0.5])


set(gcf, 'unit', 'centimeters', 'position', [35 10 7.8 6]);
set(gca,'unit', 'centimeters', 'position', [1.2, 1.5, 6, 4]);
xlh = get(gca,'xlabel');
xlv = get(gca,'ylabel');
xlh.Position(2) = xlh.Position(2) + 0.007;
xlv.Position(1) = xlv.Position(1) + 0.005;
ax = gca; % current axes
ax.XRuler.TickLabelGapOffset = -3; % negative numbers move the ticklabels down (positive -> up)
ax.YRuler.TickLabelGapOffset = -1.5; % negative numbers move the ticklabels right (negative -> left)
fileName = char(startDate)+"_fft.png";  % 为文件指定名称
fullFileName = fullfile(folderName, fileName);  % 创建完整的文件路径
print(gcf, fullFileName,'-r300','-dpng');
% exportgraphics(gcf, 'test2.png','Resolution',300);
% copygraphics(gcf)

close all
figureIdx=0;
end

folderName = 'results';  % 为目标文件夹指定名称

% 确保文件夹存在
if ~exist(folderName, 'dir')
    mkdir(folderName);  % 如果不存在，则创建文件夹
end

% 将 f_dom 添加到 vivTable 表格中作为一个新列
vivTable.f_dom = f_dom';

% 显示更新后的表格
disp(vivTable);

% 假设您的表格名为 vivTable
filename = 'vivTable.xlsx';  % 您想要保存的文件名
writetable(vivTable,fullfile(folderName, filename));  % 将表格保存到Excel文件中








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