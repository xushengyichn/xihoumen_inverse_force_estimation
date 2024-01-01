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




%% plot

%plot finite element model
[figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
% Define color, point size and marker for scatter plot
color = 'b'; % Color blue
marker = 'o'; % Circle marker
size = 15; % Size of marker
% scatter3(nodes.X, nodes.Y, nodes.Z, size, color, marker, 'filled'); % scatter plot for nodes
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


%% plot FEM with different element
close all
[X,Y,Z,nodes,elements]=plot_FEM_with_element('nodes.txt','elements.txt');
width= 14;
height=3.51;
r = 1.6;
r2 = 0.8;
numSides=20;
width_t= 10;
height_t=10;
width_c= 1.5;
height_c=3.51;
% 获取所有字段名
fields = fieldnames(X);
h = figure
k1 = 1;
fieldName = fields{k1};
X_temp = X.(fieldName);
Y_temp = Y.(fieldName);
Z_temp = Z.(fieldName);
for k2 = 1:length(X_temp)
    [vertices,faces]=drawRectangularTube(X_temp(:,k2), Y_temp(:,k2), Z_temp(:,k2), width, height);
        % 绘制矩形管道
    patch('Vertices', vertices, 'Faces', faces, 'FaceColor', '#d8802f','EdgeColor','none');
end

k1 = 7;
fieldName = fields{k1};
X_temp = X.(fieldName);
Y_temp = Y.(fieldName);
Z_temp = Z.(fieldName);
for k2 = 1:length(X_temp)
    [vertices,faces]=drawRectangularTube(X_temp(:,k2), Y_temp(:,k2), Z_temp(:,k2), width_c, height_c);
        % 绘制矩形管道
    patch('Vertices', vertices, 'Faces', faces, 'FaceColor', '#d8802f','EdgeColor','none');
end

k1 = 2;
fieldName = fields{k1};
X_temp = X.(fieldName);
Y_temp = Y.(fieldName);
Z_temp = Z.(fieldName);
for k2 = 1:length(X_temp)
    [vertices,faces]=drawCylindricalTube(X_temp(:,k2), Y_temp(:,k2), Z_temp(:,k2), r, numSides);
        % 绘制矩形管道
    patch('Vertices', vertices, 'Faces', faces, 'FaceColor', '#9fa0a0','EdgeColor','none');
end



k1 = 3;
fieldName = fields{k1};
X_temp = X.(fieldName);
Y_temp = Y.(fieldName);
Z_temp = Z.(fieldName);
for k2 = 1:length(X_temp)
    [vertices,faces]=drawCylindricalTube(X_temp(:,k2), Y_temp(:,k2), Z_temp(:,k2), r, numSides);
        % 绘制矩形管道
    patch('Vertices', vertices, 'Faces', faces, 'FaceColor', '#9fa0a0','EdgeColor','none');
end


k1 = 4;
fieldName = fields{k1};
X_temp = X.(fieldName);
Y_temp = Y.(fieldName);
Z_temp = Z.(fieldName);
for k2 = 1:length(X_temp)
    [vertices,faces]=drawCylindricalTube(X_temp(:,k2), Y_temp(:,k2), Z_temp(:,k2), r, numSides);
        % 绘制矩形管道
    patch('Vertices', vertices, 'Faces', faces, 'FaceColor', '#9fa0a0','EdgeColor','none');
end


k1 = 5;
fieldName = fields{k1};
X_temp = X.(fieldName);
Y_temp = Y.(fieldName);
Z_temp = Z.(fieldName);
for k2 = 1:length(X_temp)
    [vertices,faces]=drawCylindricalTube(X_temp(:,k2), Y_temp(:,k2), Z_temp(:,k2), r2, numSides);
        % 绘制矩形管道
    patch('Vertices', vertices, 'Faces', faces, 'FaceColor', '#9fa0a0','EdgeColor','none');
end

k1 = 6;
fieldName = fields{k1};
X_temp = X.(fieldName);
Y_temp = Y.(fieldName);
Z_temp = Z.(fieldName);
for k2 = 1:length(X_temp)
    [vertices,faces]=drawRectangularTube(X_temp(:,k2), Y_temp(:,k2), Z_temp(:,k2), width_t, height_t);
        % 绘制矩形管道
    patch('Vertices', vertices, 'Faces', faces, 'FaceColor', '#9fa0a0','EdgeColor','none');
end



k1 = 10;
fieldName = fields{k1};
X_temp = X.(fieldName);
Y_temp = Y.(fieldName);
Z_temp = Z.(fieldName);
for k2 = 1:length(X_temp)
    [vertices,faces]=drawRectangularTube(X_temp(:,k2), Y_temp(:,k2), Z_temp(:,k2), width_t, height_t);
        % 绘制矩形管道
    patch('Vertices', vertices, 'Faces', faces, 'FaceColor', '#9fa0a0','EdgeColor','none');
end




grid off; % 或者添加这行来明确关闭网格
% Add labels to each axis
% xlabel('X'); ylabel('Y'); zlabel('Z'); 
axis off;
% Add title
% title('My FEM Model'); 
% View in 3D
% view(45, 35.264); 
view(25, 35.264); 
% view(1, -1);
% Use equal scaling
axis equal; 
% Add light
camlight('headlight');
camlight('left');
camlight('right');
% Use gouraud lighting
lighting phong;
hold off;
% 设置背景透明
% set(gca, 'Color', 'none'); % 去除坐标轴背景
% set(gcf, 'Color', 'none'); % 去除整个图形的背景
% 去除刻度轴
set(gca, 'xtick', [], 'ytick', [], 'ztick', []);
print -clipboard -dbitmap
figureIdx=figureIdx-1
set(h, 'Position', [100, 100, 1920, 1080]); % 设置为更高的分辨率
set(h, 'PaperPositionMode', 'manual', 'PaperPosition', [0 0 20 20]);

set(gcf, 'Renderer', 'painters'); % 使用 'painters' 渲染器
print(h, 'MyHighQualityOutput.pdf', '-dpdf', '-r2400');


function [vertices,faces]=drawRectangularTube(X, Y, Z, width, height)
    % 计算方向向量
    direction = [X(2) - X(1), Y(2) - Y(1), Z(2) - Z(1)];
    direction = direction / norm(direction);
    
    % 计算垂直于方向向量的两个向量
    if direction(1) == 0 && direction(2) == 0
        v1 = cross(direction, [1, 0, 0]);
    else
        v1 = cross(direction, [0, 0, 1]);
    end
    v1 = v1 / norm(v1);
    v2 = cross(direction, v1);
    v2 = v2 / norm(v2);

    % 计算矩形的四个角点
    corners = [v1 * width / 2 + v2 * height / 2; 
               -v1 * width / 2 + v2 * height / 2;
               -v1 * width / 2 - v2 * height / 2;
                v1 * width / 2 - v2 * height / 2];
    
    % 创建管道的八个顶点
    vertices = [X(1), Y(1), Z(1)] + corners;
    vertices = [vertices; [X(2), Y(2), Z(2)] + corners];

    % 定义矩形面
    faces = [1 2 6 5; 2 3 7 6; 3 4 8 7; 4 1 5 8; 1 2 3 4; 5 6 7 8];


end


function [vertices,faces]=drawCylindricalTube(X, Y, Z, radius, numSides)
    % 计算方向向量
    direction = [X(2) - X(1), Y(2) - Y(1), Z(2) - Z(1)];
    direction = direction / norm(direction);
    
    % 计算垂直于方向向量的两个向量
    if direction(1) == 0 && direction(2) == 0
        v1 = cross(direction, [1, 0, 0]);
    else
        v1 = cross(direction, [0, 0, 1]);
    end
    v1 = v1 / norm(v1);
    v2 = cross(direction, v1);
    v2 = v2 / norm(v2);

    % 生成圆柱的顶点
    theta = linspace(0, 2*pi, numSides);
    circlePoints = radius * [cos(theta); sin(theta)];
    verticesTop = circlePoints' * [v1; v2] + repmat([X(1), Y(1), Z(1)], numSides, 1);
    verticesBottom = circlePoints' * [v1; v2] + repmat([X(2), Y(2), Z(2)], numSides, 1);
    vertices = [verticesTop; verticesBottom];

    % 定义圆柱体的面
    faces = [];
    for i = 1:numSides-1
        faces = [faces; i i+1 i+1+numSides i+numSides];
    end
    faces = [faces; numSides 1 1+numSides 2*numSides];

    % 绘制圆柱形管道
    % patch('Vertices', vertices, 'Faces', faces, 'FaceColor', 'blue');
end