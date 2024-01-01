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
function [X_struct,Y_struct,Z_struct,nodes,elements]=plot_FEM_with_element(nodes_table,elements_table)


%% 1 Finite Element Model

% finite element model information for visualization
% /OUTPUT, 'nodes', txt
% NLIST
% /OUTPUT
%
% /OUTPUT, 'elements', txt
% ELIST
% /OUTPUT
nodes = readtable(nodes_table); %This file is exported from ANSYS using the code above
nodes = rmmissing(nodes);
nodes_temp = nodes.Y;
nodes.Y = nodes.Z;
nodes.Z = nodes_temp;
elements = readtable(elements_table);
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
% nodes_info = nodes;

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

    % 获取elements第二列的唯一值
    uniqueMaterials = unique(elements.Var2);
    
    % 初始化结构体
    X_struct = struct();
    Y_struct = struct();
    Z_struct = struct();

    % 遍历每种材料或类型
    for i = 1:length(uniqueMaterials)
        material = uniqueMaterials(i);
        
        % 提取对应材料的元素
        materialElements = elements(elements.Var2 == material, :);

        % 初始化节点坐标数组
        materialNode1 = [];
        materialNode2 = [];

        % 遍历每个元素
        for k1 = 1:size(materialElements, 1)
            % Get node indices for the k1-th element
            node1Index = materialElements{k1, 7};
            node2Index = materialElements{k1, 8};

            % Get the node coordinates
            materialNode1(k1, :) = nodes{nodes{:, 1} == node1Index, 2:4};
            materialNode2(k1, :) = nodes{nodes{:, 1} == node2Index, 2:4};
        end

        % 存储到结构体中
        X_struct.(['mat', num2str(material)]) = [materialNode1(:,1), materialNode2(:,1)].';
        Y_struct.(['mat', num2str(material)]) = [materialNode1(:,2), materialNode2(:,2)].';
        Z_struct.(['mat', num2str(material)]) = [materialNode1(:,3), materialNode2(:,3)].';
    end
end
