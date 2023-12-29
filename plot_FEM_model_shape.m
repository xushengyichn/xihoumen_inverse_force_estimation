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
addpath(genpath("F:\git\xihoumen_inverse_force_estimation\FEM_model"))
addpath(genpath("C:\Users\xushe\OneDrive\NAS云同步\Drive\0博士研究生\3大论文\研究内容\研究内容 3：气动力模型及参数识别；\反算气动力"))

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
[X,Y,Z,nodes,elements]=plot_FEM('nodes.txt','elements.txt');

%% modal shape

KMatrix = 'KMatrix.matrix';
MMatrix = 'MMatrix.matrix';
NODEONDECK = 'nodeondeck.txt';
Mapping = 'KMatrix.mapping';
nodegap = 'nodegap.txt';

output = ImportMK(100, 'KMatrix.matrix', 'MMatrix.matrix', 'nodeondeck.txt', 'KMatrix.mapping', 'nodegap.txt');
node_modal_shape=[nodes.NODE, nodes.X, nodes.Y, nodes.Z];
% node_modal_shape_original = node_modal_shape;
Mapping_data = output.Mapping;
modal_shape_all = output.eig_vec;
mode_num  =23;
modal_shape_plot = modal_shape_all(:,mode_num);
modal_shape_plot = modal_shape_plot/max(abs(modal_shape_plot))*50;


[node_modal_shape,nodes_displacement,nodes_displacementMagnitude,X_modal_shape,Y_modal_shape,Z_modal_shape]= update_node(node_modal_shape,modal_shape_plot,Mapping_data,elements);




%% plot
%plot finite element model
[figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
% Define color, point size and marker for scatter plot
color = 'b'; % Color blue
marker = 'o'; % Circle marker
size = 15; % Size of marker
scatter3(nodes.X, nodes.Y, nodes.Z, size, color, marker, 'filled'); % scatter plot for nodes
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
% view(45, 35.264); 
view(1, -1);
% Use equal scaling
axis equal; 
% Add light
camlight('headlight');
% Use gouraud lighting
lighting gouraud;
hold off;
figureIdx=figureIdx-1

%plot finite element model
[figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
% Define color, point size and marker for scatter plot
color = 'b'; % Color blue
marker = 'o'; % Circle marker
size = 15; % Size of marker
scatter3(nodes.X, nodes.Y, nodes.Z, size, color, marker, 'filled'); % scatter plot for nodes
% Define line color for beams
beamColor = 'k'; % Black color
hold on; % hold current plot
% Plot each beam
% Plot all lines at once
line(X, Y, Z, 'Color', beamColor, 'LineWidth', lineWidthThin);
% Add grids
grid off; % 或者添加这行来明确关闭网格
% Add labels to each axis
% xlabel('X'); ylabel('Y'); zlabel('Z'); 
axis off;
% Add title
% title('My FEM Model'); 
% View in 3D
% view(45, 35.264); 
view(1, -1);
% Use equal scaling
axis equal; 
% Add light
camlight('headlight');
% Use gouraud lighting
lighting gouraud;
hold off;
% 设置背景透明
% set(gca, 'Color', 'none'); % 去除坐标轴背景
% set(gcf, 'Color', 'none'); % 去除整个图形的背景
% 去除刻度轴
set(gca, 'xtick', [], 'ytick', [], 'ztick', []);
print -clipboard -dbitmap
figureIdx=figureIdx-1

% plot modal shape
[figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
% Define color, point size and marker for scatter plot
color = 'b'; % Color blue
marker = 'o'; % Circle marker
size = 15; % Size of marker
% scatter3(nodes_info.X, nodes_info.Y, nodes_info.Z, size, color, marker, 'filled'); % scatter plot for nodes
hold on
% scatter3(node_modal_shape(:,2),node_modal_shape(:,3),node_modal_shape(:,4), size, color, marker, 'filled'); % scatter plot for nodes
% Create a scatter plot using displacementMagnitude as color data
scatterHandle = scatter3(node_modal_shape(:,2), node_modal_shape(:,3), node_modal_shape(:,4), size, nodes_displacementMagnitude, marker, 'filled');

% Apply a colormap
colormap(jet);  % or use another colormap if you prefer

% Add a colorbar
colorbar;

% If needed, you can set the limits of the color scale
% caxis([minDisplacement, maxDisplacement]);  % replace with actual min and max if needed

% Define line color for beams
beamColor = 'k'; % Black color
hold on; % hold current plot
% Plot each beam
% Plot all lines at once
line(X_modal_shape, Y_modal_shape, Z_modal_shape, 'Color', beamColor, 'LineWidth', lineWidthThin);
% Add grids
grid on;  
% Add labels to each axis
xlabel('X'); ylabel('Y'); zlabel('Z'); 
% Add title
title('My FEM Model'); 
% View in 3D
view(45, 35.264); 
% view(0,0); 
% Use equal scaling
axis equal; 
% Add light
camlight('headlight');
% Use gouraud lighting
lighting gouraud;
hold off;

% plot modal shape
[figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
% Define color, point size and marker for scatter plot
color = 'b'; % Color blue
marker = 'o'; % Circle marker
size = 15; % Size of marker
% scatter3(nodes_info.X, nodes_info.Y, nodes_info.Z, size, color, marker, 'filled'); % scatter plot for nodes
hold on
% scatter3(node_modal_shape(:,2),node_modal_shape(:,3),node_modal_shape(:,4), size, color, marker, 'filled'); % scatter plot for nodes
% Create a scatter plot using displacementMagnitude as color data
scatterHandle = scatter3(node_modal_shape(:,2), node_modal_shape(:,3), node_modal_shape(:,4), size, nodes_displacementMagnitude, marker, 'filled');

% Apply a colormap
colormap(jet);  % or use another colormap if you prefer

% Add a colorbar
colorbar;

% If needed, you can set the limits of the color scale
% caxis([minDisplacement, maxDisplacement]);  % replace with actual min and max if needed

% Define line color for beams
beamColor = 'k'; % Black color
hold on; % hold current plot
% Plot each beam
% Plot all lines at once
line(X_modal_shape, Y_modal_shape, Z_modal_shape, 'Color', beamColor, 'LineWidth', lineWidthThin);
% Add grids
grid on;  
% Add labels to each axis
xlabel('X'); ylabel('Y'); zlabel('Z'); 
% Add title
title('My FEM Model'); 
% View in 3D
% view(45, 35.264); 
view(0,0); 
% Use equal scaling
axis equal; 
% Add light
camlight('headlight');
% Use gouraud lighting
lighting gouraud;
hold off;