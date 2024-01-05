%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%Date: 2023-06-05 10:56:30
%LastEditors: ShengyiXu xushengyichn@outlook.com
%LastEditTime: 2024-01-05 11:26:06
%FilePath: \Exercises-for-Techniques-for-estimation-in-dynamics-systemsc:\Users\xushe\OneDrive\NAS云同步\Drive\0博士研究生\4其他科研工作\2023年12月19日反算涡激力\manuscript\matlab\HMS_3D\HMS_3D.m
%Description: Inverse calculation of vortex-induced force of Xihoumen Bridge
%
%Copyright (c) 2023 by ${git_name_email}, All Rights Reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
close all
%% General Parameters
addpath('FEM_model\')
run('CommonCommand.m');


%% finite element model information for visualization
% /OUTPUT, 'nodes', txt
% NLIST
% /OUTPUT
%
% /OUTPUT, 'elements', txt
% ELIST
% /OUTPUT
h = figure

%% plot FEM with different element
[X,Y,Z,nodes,elements]=plot_FEM_with_element('nodes.txt','elements.txt','if_modeshape',false);
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

k1 = 1;
fieldName = fields{k1};
X_temp = X.(fieldName);
Y_temp = Y.(fieldName);
Z_temp = Z.(fieldName);
for k2 = 1:length(X_temp)
    [vertices,faces]=drawRectangularTube(X_temp(:,k2), Y_temp(:,k2), Z_temp(:,k2), width, height);
        % 绘制矩形管道
    patch('Vertices', vertices, 'Faces', faces, 'FaceColor', 'black','EdgeColor','none', 'FaceAlpha', 0.5);
end

k1 = 7;
fieldName = fields{k1};
X_temp = X.(fieldName);
Y_temp = Y.(fieldName);
Z_temp = Z.(fieldName);
for k2 = 1:length(X_temp)
    [vertices,faces]=drawRectangularTube(X_temp(:,k2), Y_temp(:,k2), Z_temp(:,k2), width_c, height_c);
        % 绘制矩形管道
    patch('Vertices', vertices, 'Faces', faces, 'FaceColor', 'black','EdgeColor','none', 'FaceAlpha', 0.5);
end

%% plot FEM with different element with modal shape
output = ImportMK(100, 'KMatrix.matrix', 'MMatrix.matrix', 'nodeondeck.txt', 'KMatrix.mapping', 'nodegap.txt');
output.mode_num= 3;
output.amplify = 50;
disp(output.Freq(output.mode_num))
[X,Y,Z,nodes,elements]=plot_FEM_with_element('nodes.txt','elements.txt','if_modeshape',true,'FEM_data',output);

width= 14;
height=20;
r = 1.6;
r2 = 0.8;
numSides=20;
width_t= 10;
height_t=10;
width_c= 1.5;
height_c=3.51;
% 获取所有字段名
fields = fieldnames(X);

% for k1= 1:length(fields)
%     fieldName = fields{k1};
%     Y.(fieldName) = Y.(fieldName) ;
% end


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

%% set the ball
%生成球
% 塔1 21300 塔 22950
%% left
leftoffset = -50;
x_center = 21300+(22950-21300)/4;
y_center = -15.699;
z_center = 85+leftoffset;
radius = 10;
% 生成球的坐标
[X, Y, Z] = sphere;
% 调整球的坐标
X = X * radius + x_center;
Y = Y * radius + y_center;
Z = Z * radius + z_center;
% 绘制球
hold on
surf(X, Y, Z, 'FaceColor', 'red', 'EdgeColor', 'none','FaceAlpha',0.5);  % 设置表面颜色为红色，边缘颜色无

x_center = 21300+(22950-21300)/4;
y_center = -15.699;
z_center = 65+leftoffset ;
radius = 10;
% 生成球的坐标
[X, Y, Z] = sphere;
% 调整球的坐标
X = X * radius + x_center;
Y = Y * radius + y_center;
Z = Z * radius + z_center;
% 绘制球
hold on
surf(X, Y, Z, 'FaceColor', 'blue', 'EdgeColor', 'none','FaceAlpha',0.5);  % 设置表面颜色为红色，边缘颜色无

x_center = 21300+(22950-21300)/4;
y_center = 15.699;
z_center = 85+leftoffset;
radius = 10;
% 生成球的坐标
[X, Y, Z] = sphere;
% 调整球的坐标
X = X * radius + x_center;
Y = Y * radius + y_center;
Z = Z * radius + z_center;
% 绘制球
hold on
surf(X, Y, Z, 'FaceColor', 'red', 'EdgeColor', 'none','FaceAlpha',0.5);  % 设置表面颜色为红色，边缘颜色无

x_center = 21300+(22950-21300)/4;
y_center = 15.699;
z_center = 65+leftoffset ;
radius = 10;
% 生成球的坐标
[X, Y, Z] = sphere;
% 调整球的坐标
X = X * radius + x_center;
Y = Y * radius + y_center;
Z = Z * radius + z_center;
% 绘制球
hold on
surf(X, Y, Z, 'FaceColor', 'blue', 'EdgeColor', 'none','FaceAlpha',0.5);  % 设置表面颜色为红色，边缘颜色无

%% middle
middleoffset = 3;
x_center = 21300+(22950-21300)/4*2;
y_center = -15.699;
z_center = 85+middleoffset;
radius = 10;
% 生成球的坐标
[X, Y, Z] = sphere;
% 调整球的坐标
X = X * radius + x_center;
Y = Y * radius + y_center;
Z = Z * radius + z_center;
% 绘制球
hold on
surf(X, Y, Z, 'FaceColor', 'red', 'EdgeColor', 'none','FaceAlpha',0.5);  % 设置表面颜色为红色，边缘颜色无

x_center = 21300+(22950-21300)/4*2;
y_center = -15.699;
z_center = 65+middleoffset ;
radius = 10;
% 生成球的坐标
[X, Y, Z] = sphere;
% 调整球的坐标
X = X * radius + x_center;
Y = Y * radius + y_center;
Z = Z * radius + z_center;
% 绘制球
hold on
surf(X, Y, Z, 'FaceColor', 'blue', 'EdgeColor', 'none','FaceAlpha',0.5);  % 设置表面颜色为红色，边缘颜色无

x_center = 21300+(22950-21300)/4*2;
y_center = 15.699;
z_center = 85+middleoffset;
radius = 10;
% 生成球的坐标
[X, Y, Z] = sphere;
% 调整球的坐标
X = X * radius + x_center;
Y = Y * radius + y_center;
Z = Z * radius + z_center;
% 绘制球
hold on
surf(X, Y, Z, 'FaceColor', 'red', 'EdgeColor', 'none','FaceAlpha',0.5);  % 设置表面颜色为红色，边缘颜色无

x_center = 21300+(22950-21300)/4*2;
y_center = 15.699;
z_center = 65+middleoffset ;
radius = 10;
% 生成球的坐标
[X, Y, Z] = sphere;
% 调整球的坐标
X = X * radius + x_center;
Y = Y * radius + y_center;
Z = Z * radius + z_center;
% 绘制球
hold on
surf(X, Y, Z, 'FaceColor', 'blue', 'EdgeColor', 'none','FaceAlpha',0.5);  % 设置表面颜色为红色，边缘颜色无

%% right
rightoffset = 43;
x_center = 21300+(22950-21300)/4*3;
y_center = -15.699;
z_center = 85+rightoffset;
radius = 10;
% 生成球的坐标
[X, Y, Z] = sphere;
% 调整球的坐标
X = X * radius + x_center;
Y = Y * radius + y_center;
Z = Z * radius + z_center;
% 绘制球
hold on
surf(X, Y, Z, 'FaceColor', 'red', 'EdgeColor', 'none','FaceAlpha',0.5);  % 设置表面颜色为红色，边缘颜色无

x_center = 21300+(22950-21300)/4*3;
y_center = -15.699;
z_center = 65+rightoffset ;
radius = 10;
% 生成球的坐标
[X, Y, Z] = sphere;
% 调整球的坐标
X = X * radius + x_center;
Y = Y * radius + y_center;
Z = Z * radius + z_center;
% 绘制球
hold on
surf(X, Y, Z, 'FaceColor', 'blue', 'EdgeColor', 'none','FaceAlpha',0.5);  % 设置表面颜色为红色，边缘颜色无

x_center = 21300+(22950-21300)/4*3;
y_center = 15.699;
z_center = 85+rightoffset;
radius = 10;
% 生成球的坐标
[X, Y, Z] = sphere;
% 调整球的坐标
X = X * radius + x_center;
Y = Y * radius + y_center;
Z = Z * radius + z_center;
% 绘制球
hold on
surf(X, Y, Z, 'FaceColor', 'red', 'EdgeColor', 'none','FaceAlpha',0.5);  % 设置表面颜色为红色，边缘颜色无

x_center = 21300+(22950-21300)/4*3;
y_center = 15.699;
z_center = 65 +rightoffset;
radius = 10;
% 生成球的坐标
[X, Y, Z] = sphere;
% 调整球的坐标
X = X * radius + x_center;
Y = Y * radius + y_center;
Z = Z * radius + z_center;
% 绘制球
hold on
surf(X, Y, Z, 'FaceColor', 'blue', 'EdgeColor', 'none','FaceAlpha',0.5);  % 设置表面颜色为红色，边缘颜色无

%% set camera and light
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
light('Color',[1 1 1],'Position',[22000 -20000 75],'Style','local');

% Use gouraud lighting
lighting phong;
material('dull'); % 或者 'shiny', 'metal'
hold off;


% 设置摄像机的位置
campos([22000 -20000 75]);
% % 设置摄像机的观察目标
camtarget([22000 0 75]);
% % 设置摄像机的“上”方向
% camup([upX upY upZ]);
% % 设置摄像机的视角宽度
camva(5);
% 选择投影类型
camproj('orthographic');
% camproj('perspective');
% 更新图形窗口
drawnow;



% % 设置背景透明
set(gca, 'Color', 'none'); % 去除坐标轴背景
set(gcf, 'Color', 'none'); % 去除整个图形的背景
% 去除刻度轴
set(gca, 'xtick', [], 'ytick', [], 'ztick', []);
% print -clipboard -dbitmap
% figureIdx=figureIdx-1
set(h, 'Position', [100, 100, 1920, 1080]); % 设置为更高的分辨率
% set(h, 'PaperPositionMode', 'manual', 'PaperPosition', [0 0 20 20]);

set(gcf, 'Renderer', 'opengl');
set(gcf,'Renderermode','manual');
% opengl('hardware');

% print(h, 'MyHighQualityOutput.pdf', '-dpdf', '-r1200');
scriptName = mfilename;
filename = scriptName+".png";
% exportgraphics(h, filename, 'Resolution', 300,'BackgroundColor', 'none');



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