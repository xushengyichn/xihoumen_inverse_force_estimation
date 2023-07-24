%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%Date: 2023-06-05 10:56:30
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2023-06-05 16:33:01
%FilePath: \code\Main.m
%Description: Inverse calculation of vortex-induced force of Xihoumen Bridge
%
%Copyright (c) 2023 by ${git_name_email}, All Rights Reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
close all
%% General Parameters
addpath(genpath("F:\git\ssm_tools"))
addpath(genpath("D:\Users\xushe\Documents\GitHub\ssm_tools"))
addpath(genpath("D:\Users\xushe\Documents\GitHub\Function_shengyi_package"))
addpath(genpath("F:\git\Function_shengyi_package"))
subStreamNumberDefault = 2132;
run('InitScript.m');
%% 0 绘图参数
fig_bool = ON;
num_figs_in_row = 5; %每一行显示几个图
figPos = figPosMedium; %图的大小，参数基于InitScript.m中的设置
%设置图片间隔
gap_between_images = [0, 0];
figureIdx = 0;

%% finite element model information for visualization
% /OUTPUT, 'nodes', txt
% NLIST
% /OUTPUT
%
% /OUTPUT, 'elements', txt
% ELIST
% /OUTPUT
[nodes_info,X,Y,Z]=plot_FEM('nodes.txt','elements.txt');


%% plot
%plot finite element model
[figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
% Define color, point size and marker for scatter plot
color = 'b'; % Color blue
marker = 'o'; % Circle marker
size = 15; % Size of marker
scatter3(nodes_info.X, nodes_info.Y, nodes_info.Z, size, color, marker, 'filled'); % scatter plot for nodes
% Define line color for beams
beamColor = 'k'; % Black color
hold on; % hold current plot
% Plot each beam
% Plot all lines at once
line(X, Y, Z, 'Color', beamColor, 'LineWidth', lineWidthThin);
% Add grids
grid on;  
% Add labels to each axis
xlabel('X'); ylabel('Y'); zlabel('Z'); 
% Add title
title('My FEM Model'); 
% View in 3D
view(45, 35.264); 
% Use equal scaling
axis equal; 
% Add light
camlight('headlight');
% Use gouraud lighting
lighting gouraud;
hold off;
