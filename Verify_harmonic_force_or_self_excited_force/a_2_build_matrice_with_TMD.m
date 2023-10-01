%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%Date: 2023-01-25 10:57:55
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2023-01-25 14:04:43
%FilePath: \20230124优化问题\a_2_build_matrice_with_TMD.m
%Description: 建立质量、刚度、阻尼矩阵，考虑多项式非线性气动力
%
%Copyright (c) 2023 by xushengyichn 54436848+xushengyichn@users.noreply.github.com, All Rights Reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [output] = a_2_build_matrice_with_TMD(input)
%a_1_build_matrice 建立质量、刚度、阻尼矩阵
%   input: 输入结构体
%   output: 输出结构体
%   读取input中的输入数据
number_of_modes_to_control = input.number_of_modes_to_control;
number_of_modes_to_consider = input.number_of_modes_to_consider;
number_of_tmds = input.number_of_tmds;
TMDs_mass = input.TMDs_mass;
TMDs_frequency=input.TMDs_frequency;
TMDs_damping_ratio = input.TMDs_damping_ratio;
TMDs_location = input.TMDs_location;
modal_damping_ratios = input.modal_damping_ratios;
nodeondeck=input.nodeondeck;
mode=input.mode;
nodegap=input.nodegap;
modecal=input.modecal;

% 变量替换（目的是简化变量名）Variable substitution (to simplify variable names)
numberofTMD = number_of_tmds; % 所需要计算的TMD的数量. The number of TMDs to calculate.
nTMD = numberofTMD; % 所需要计算的TMD的数量. The number of TMDs to calculate.
nModes = number_of_modes_to_consider; % 所需要计算的模态的数量. The number of modes to calculate.
xTMD = TMDs_location; % TMD的位置. The location of TMD.
mTMD = TMDs_mass; % TMD的质量. The mass of TMD.
omegatmd = 2 * pi * TMDs_frequency;
kTMD = mTMD .* omegatmd .^ 2;  % TMD的刚度. The stiffness of TMD.
cTMD = 2 * TMDs_damping_ratio .* sqrt(mTMD .* kTMD); % TMD的阻尼. The damping of TMD.

MM_eq = input.MM_eq;
KK_eq = input.KK_eq;
eig_val = input.eig_val;
eig_vec = input.eig_vec;

omeg = sqrt(diag(eig_val));
Freq = omeg / (2 * pi);

%% 构建包含TMD的MCK矩阵
% % 导入矩阵序号对应节点自由度关系。导入KCM中 *.mapping 文件的任意一个即可。他们是一样的。
% % The imported matrix number corresponds to the node degree of freedom relationship. Import any one of the *.mapping files in KCM. They are the same.
% KMmapping = importmappingmatrix('KMatrix.mapping');
% UYNode = sort(KMmapping.Node(KMmapping.DOF == 'UY'));

% % 计算桥面节点的振型向量
% nodeondeck = importdata('nodeondeck.txt');
% mode = zeros(length(nodeondeck), nModes);

% for k1 = 1:length(nodeondeck)
%     position_index = KMmapping.MatrixEqn(find(and(KMmapping.Node == nodeondeck(k1), KMmapping.DOF == 'UY')));

%     if isempty(position_index)
%         mode(k1, :) = zeros(1, nModes);
%     else
%         mode(k1, :) = eig_vec(position_index, :);
%     end

% end

% [phideckmax, location] = max(mode); %计算桥面每个模态的最大振型值
% nodegap = importdata('nodegap.txt');

% for k1 = 1:length(nodegap) - 1
%     nodedis(k1, 1) = nodegap(k1 + 1) - nodegap(k1);
% end

% clear k1
% % 模态振型比节点间距多一个点，所以将两节点的模态取平均值计算振型积分
% for k1 = 1:length(nodegap) - 1
%     modecal(k1, :) = (mode(k1 + 1, :) + mode(k1, :)) / 2;
% end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 计算考虑TMDs的质量、刚度、阻尼矩阵
% Calculate mass, stiffness, damping matrix considering TMDs
%% 考虑TMD振动响应的模态叠加法
%% Modal Superposition Method Considering TMD Vibration Response
matrixsize = nTMD + nModes;
phiTMD = zeros(nTMD, nModes);
% phiTMD row:TMD for each loaction column:the mode shape at the each
% location of tmd
% phiTMD 行：TMD的位置 列：TMD位置的模态振型
for t1 = 1:nTMD

    for t2 = 1:nModes
        [~, index] = sort(abs(nodegap - xTMD(t1))); %查找与xTMD最接近的点的排序
        xResult = nodegap(index(1:2)); %获取最接近的两个点的x坐标
        mode2nodes = mode(index(1:2), 1:nModes); %获取两个点坐标的y值
        phi_result = interp1(xResult, mode2nodes, xTMD(t1), 'linear', 'extrap'); %插值以后任意点的振型
        %         disp(phi_result)
        phiTMD(t1, t2) = phi_result(t2);

    end

end

% 创建质量矩阵
% Create the Mass matrix
MM = zeros(matrixsize, matrixsize);
MM(1:nModes, 1:nModes) = MM_eq;

for k1 = nModes + 1:matrixsize
    MM(k1, k1) = mTMD(k1 - nModes);
end
clear k1
% 创建刚度矩阵
% Create the Stiffness Matrix
KK = zeros(size(MM, 1), size(MM, 2));
KK1 = zeros(size(MM, 1), size(MM, 2));
KK1(1:nModes, 1:nModes) = KK_eq;

for k1 = nModes + 1:matrixsize
    KK1(k1, k1) = kTMD(k1 - nModes);
end

KK2 = zeros(matrixsize, matrixsize);
for k1 = 1:nModes % 第k1行

    for k2 = 1:nModes % 第k2列

        for k3 = 1:nTMD % 第k3个TMD
            KK2(k1, k2) = KK2(k1, k2) + kTMD(k3) * phiTMD(k3, k1) * phiTMD(k3, k2);
        end

    end

end
clear k1 k2 k3
for k1 = 1:nModes

    for k2 = 1:nTMD
        KK2(k1, k2 + nModes) = -kTMD(k2) * phiTMD(k2, k1);
        KK2(k2 + nModes, k1) = -kTMD(k2) * phiTMD(k2, k1);
    end

end

KK = KK1 + KK2;
clear k1 k2

% 创建阻尼矩阵
% Create the Damping Matrix
CC = zeros(size(MM, 1), size(MM, 2));
CC1 = zeros(size(MM, 1), size(MM, 2));
CC2 = zeros(size(MM, 1), size(MM, 2));

for k1 = 1:matrixsize

    if k1 <= nModes
        CC1(k1, k1) = modal_damping_ratios(k1) * 4 * pi * MM(k1, k1) * Freq(k1);
    elseif k1 > nModes
        CC1(k1, k1) = cTMD(k1 - nModes);
    end

end

for k1 = 1:nModes % 第k1行

    for k2 = 1:nModes % 第k2列

        for k3 = 1:nTMD % 第k3个TMD
            CC2(k1, k2) = CC2(k1, k2) + cTMD(k3) * phiTMD(k3, k1) * phiTMD(k3, k2);
        end

    end

end

clear k1 k2

for k1 = 1:nModes

    for k2 = 1:nTMD
        CC2(k1, k2 + nModes) = -cTMD(k2) * phiTMD(k2, k1);
        CC2(k2 + nModes, k1) = -cTMD(k2) * phiTMD(k2, k1);
    end

end

CC = CC1 + CC2;
clear k1 k2

output.Mass = MM;
output.Stiffness = KK;
output.Damping = CC;
output.Freq = Freq;
output.mode = mode;
output.phiTMD = phiTMD;
output.modecal = modecal;
end