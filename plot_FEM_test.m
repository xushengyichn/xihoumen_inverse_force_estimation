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

%% Drawing Parameters
num_figs_in_row = 3; %number of figures in a row
figPos = figPosLarge; %The size of the graph, the parameters are based on the settings in InitScript.m
gap_between_images = [100 100]; %The gap between the images, the first is the horizontal gap, the second is the vertical gap
figureIdx = 0; %The serial number of the figure

%% 1 Finite Element Model

% finite element model information for visualization
% /OUTPUT, 'nodes', txt
% NLIST
% /OUTPUT
%
% /OUTPUT, 'elements', txt
% ELIST
% /OUTPUT
nodes = readtable('nodes.txt'); %This file is exported from ANSYS using the code above
nodes = rmmissing(nodes);
nodes_temp = nodes.Y;
nodes.Y = nodes.Z;
nodes.Z = nodes_temp;
elements = readtable('elements.txt');
elements = removevars(elements, 'Var9');
elements = rmmissing(elements);
% remove nodes that are not used by any element
% Extract all unique node numbers from the 7th and 8th columns
connectedNodeNumbers = unique([elements.Var7; elements.Var8]);
% Initialize a new empty table for connected nodes
connectedNodes = table();

% Loop through all node numbers in the 'nodes' table
for i = 1:height(nodes)
    % Check if the node number is in the list of connected node numbers
    if ismember(nodes.NODE(i), connectedNodeNumbers)
        % If it is, append the node to the 'connectedNodes' table
        connectedNodes = [connectedNodes; nodes(i, :)];
    end
end

% Now, 'connectedNodes' contains only nodes connected to an element
nodes = connectedNodes;


% Number of beams
numBeams = size(elements, 1);

for k1 = 1:numBeams
    % Get node indices for the i-th beam
    node1Index = elements{k1, 7};
    node2Index = elements{k1, 8};

    % Get the node coordinates
    node1(k1, :) = nodes{nodes{:, 1} == node1Index, 2:4};
    node2(k1, :) = nodes{nodes{:, 1} == node2Index, 2:4};
end
% Prepare matrix inputs for line function
X = [node1(:,1), node2(:,1)].'; % X coordinates of line ends
Y = [node1(:,2), node2(:,2)].'; % Y coordinates of line ends
Z = [node1(:,3), node2(:,3)].'; % Z coordinates of line ends

% modal shape information for visualization
% read the stiffness matrix and mass matrix from ANSYS
KMatrix = 'KMatrix.matrix';
MMatrix = 'MMatrix.matrix';
NODEONDECK = 'NODEONDECK2.txt';
Mapping = 'KMatrix.mapping';
NODEDISTANCE = 'NODEDISTANCE.txt';
NodeCol = 2;
if ~exist('modal_info.mat', 'file')
[output] = modal_info(KMatrix, MMatrix, Mapping, NODEONDECK, NODEDISTANCE, NodeCol);
save modal_info.mat output
else
load modal_info.mat
end
node_modal_shape=[nodes.NODE, nodes.X, nodes.Y, nodes.Z];
node_modal_shape_original = node_modal_shape;
Mapping_data = output.Mapping;
modal_shape_all = output.eig_vec;
mode_num  =3;
modal_shape_plot = modal_shape_all(:,mode_num);
modal_shape_plot = modal_shape_plot/max(abs(modal_shape_plot))*50;
for k1= 1:length(modal_shape_plot)
    Node_temp= Mapping_data{k1,2};
    Node_DOF= Mapping_data{k1,3};
    if Node_DOF == "UX"
        Seq = find(node_modal_shape(:,1)==Node_temp);
        node_modal_shape(Seq,2) = node_modal_shape(Seq,2)+modal_shape_plot(k1);
    end
    if Node_DOF == "UY"
        Seq = find(node_modal_shape(:,1)==Node_temp);
        node_modal_shape(Seq,4) = node_modal_shape(Seq,4)+modal_shape_plot(k1);
    end%为了在matlab中绘图方便，这里注意一下，matlab中的Y轴对应的是Z轴，Z轴对应的是Y轴
    if Node_DOF == "UZ"
        Seq = find(node_modal_shape(:,1)==Node_temp);
        node_modal_shape(Seq,3) = node_modal_shape(Seq,3)+modal_shape_plot(k1);
    end

end
nodes_displacement = node_modal_shape(:,2:4)-node_modal_shape_original(:,2:4);
nodes_displacementMagnitude = sqrt(sum(nodes_displacement.^2, 2));

% Number of beams
numBeams = size(elements, 1);

for k1 = 1:numBeams
    % Get node indices for the i-th beam
    node1Index = elements{k1, 7};
    node2Index = elements{k1, 8};

    % Get the node coordinates
    node1(k1, :) = node_modal_shape(node_modal_shape(:, 1) == node1Index, 2:4);
    node2(k1, :) = node_modal_shape(node_modal_shape(:, 1) == node2Index, 2:4);
end
% Prepare matrix inputs for line function
X_modal_shape = [node1(:,1), node2(:,1)].'; % X coordinates of line ends
Y_modal_shape = [node1(:,2), node2(:,2)].'; % Y coordinates of line ends
Z_modal_shape = [node1(:,3), node2(:,3)].'; % Z coordinates of line ends


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
view(45, 35.264); 
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
% scatter3(nodes.X, nodes.Y, nodes.Z, size, color, marker, 'filled'); % scatter plot for nodes
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