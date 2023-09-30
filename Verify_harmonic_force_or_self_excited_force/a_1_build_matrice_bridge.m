%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%Date: 2023-01-25 10:57:55
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2023-01-25 14:03:49
%FilePath: \20230124优化问题\a_1_build_matrice_bridge.m
%Description: 建立桥梁的质量、刚度、阻尼矩阵
%
%Copyright (c) 2023 by xushengyichn 54436848+xushengyichn@users.noreply.github.com, All Rights Reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [output] = a_1_build_matrice_bridge(input)
%a_1_build_matrice 建立质量、刚度、阻尼矩阵
%   input: 输入结构体
%   output: 输出结构体
%   读取input中的输入数据
% number_of_modes_to_control=input.number_of_modes_to_control;
number_of_modes_to_consider=input.number_of_modes_to_consider;
% number_of_tmds=input.number_of_tmds;
% TMDs_mass=TMDs_massinput.TMDs_mass;
% TMDs_frequencyinput.TMDs_frequency;
% TMDs_damping_ratio=input.TMDs_damping_ratio;
% TMDs_location=input.TMDs_location;

% numberofTMD=number_of_tmds;% 所需要计算的TMD的数量. The number of TMDs to calculate.
% nTMD=numberofTMD;% 所需要计算的TMD的数量. The number of TMDs to calculate.
nModes = number_of_modes_to_consider; % 所需要计算的模态的数量. The number of modes to calculate.

modeinfo = load('modeinfo_all.mat'); % 读取未安装TMD时的模态信息. Read the mode information.
number_of_modes_to_consider=number_of_modes_to_consider;% 所需要考虑的模态数量. The number of modes to consider.

% 生成计算用质量刚度矩阵以及模态信息
MM_eq = zeros(number_of_modes_to_consider, number_of_modes_to_consider);
KK_eq = zeros(number_of_modes_to_consider, number_of_modes_to_consider);
eig_val = zeros(number_of_modes_to_consider, number_of_modes_to_consider);
eig_vec = zeros(length(modeinfo.eig_vec(:, 1)), number_of_modes_to_consider);

for k1 = 1:number_of_modes_to_consider
    MM_eq(k1,k1)=modeinfo.MM_eq(k1,k1);
    KK_eq(k1,k1)=modeinfo.KK_eq(k1,k1);
    eig_val(k1,k1)=modeinfo.eig_val(k1,k1);
    eig_vec(:,k1)=modeinfo.eig_vec(:,k1);
end


% 导入矩阵序号对应节点自由度关系。导入KCM中 *.mapping 文件的任意一个即可。他们是一样的。
% The imported matrix number corresponds to the node degree of freedom relationship. Import any one of the *.mapping files in KCM. They are the same.
KMmapping = importmappingmatrix('KMatrix.mapping');
UYNode = sort(KMmapping.Node(KMmapping.DOF == 'UY'));

% 计算桥面节点的振型向量
nodeondeck = importdata('nodeondeck.txt');
mode = zeros(length(nodeondeck), nModes);

for k1 = 1:length(nodeondeck)
    position_index = KMmapping.MatrixEqn(find(and(KMmapping.Node == nodeondeck(k1), KMmapping.DOF == 'UY')));

    if isempty(position_index)
        mode(k1, :) = zeros(1, nModes);
    else
        mode(k1, :) = eig_vec(position_index, :);
    end

end

[phideckmax, location] = max(mode); %计算桥面每个模态的最大振型值
nodegap = importdata('nodegap.txt');

for k1 = 1:length(nodegap) - 1
    nodedis(k1, 1) = nodegap(k1 + 1) - nodegap(k1);
end

clear k1
% 模态振型比节点间距多一个点，所以将两节点的模态取平均值计算振型积分
for k1 = 1:length(nodegap) - 1
    modecal(k1, :) = (mode(k1 + 1, :) + mode(k1, :)) / 2;
end

omeg = sqrt(diag(eig_val));
Freq = omeg / (2 * pi);

output.MM_eq=MM_eq;
output.KK_eq=KK_eq;
output.eig_val=eig_val;
output.eig_vec=eig_vec;
output.nodeondeck=nodeondeck;
output.mode=mode;
output.nodegap=nodegap;
output.modecal=modecal;
output.nodedis=nodedis;
output.Freq = Freq;
output.phideckmax=phideckmax;
end


