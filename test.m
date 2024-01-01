clc
clear
close all
% 测试数据
X = [0; 1];
Y = [0; 1];
Z = [0; 1];
width = 0.01; % 矩形宽度
height = 2; % 矩形高度

% figure;
% drawRectangularTube(X, Y, Z, width, height);
% axis equal;
% view(3);
drawCylindricalTube([0; 5], [0; 0], [0; 0], 1, 20);

% patch('Vertices', vertices, 'Faces', faces, 'FaceColor', 'blue');
function drawRectangularTube(X, Y, Z, width, height)
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

    % 绘制矩形管道
    patch('Vertices', vertices, 'Faces', faces, 'FaceColor', 'blue');
end


function drawCylindricalTube(X, Y, Z, radius, numSides)
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
    patch('Vertices', vertices, 'Faces', faces, 'FaceColor', 'blue');
end
