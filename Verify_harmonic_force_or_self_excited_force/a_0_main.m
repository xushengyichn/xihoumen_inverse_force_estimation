%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%Date: 2023-01-25 10:39:30
%LastEditors: ShengyiXu xushengyichn@outlook.com
%LastEditTime: 2023-09-30 21:06:45
%FilePath: \Exercises-for-Techniques-for-estimation-in-dynamics-systemsf:\git\xihoumen_inverse_force_estimation\Verify_harmonic_force_or_self_excited_force\a_0_main.m
%Description: 本函数为主函数，分别调用各个子函数完成桥梁在多阶级涡振下的响应总和计算
%
%Copyright (c) 2023 by xushengyichn 54436848+xushengyichn@users.noreply.github.com, All Rights Reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [result] = a_0_main(number_of_modes_to_control,number_of_modes_to_consider,number_of_tmds,modal_damping_ratios,t_length,TMDs_mass,TMDs_frequency,TMDs_damping_ratio,TMDs_location)

input.number_of_modes_to_control = number_of_modes_to_control;
input.number_of_modes_to_consider = number_of_modes_to_consider;
input.number_of_tmds = number_of_tmds;
input.TMDs_mass = TMDs_mass;
input.TMDs_frequency = TMDs_frequency;
input.TMDs_damping_ratio = TMDs_damping_ratio;
input.TMDs_location = TMDs_location;

% 建立质量、刚度矩阵
% Build mass and stiffness matrices

[output] = a_1_build_matrice_bridge(input);

MM_eq = output.MM_eq;
KK_eq = output.KK_eq;
eig_val = output.eig_val;
eig_vec = output.eig_vec;
nodeondeck = output.nodeondeck;
mode = output.mode;
nodegap = output.nodegap;
modecal = output.modecal;
nodedis = output.nodedis;
Freq = output.Freq;
phideckmax = output.phideckmax;

% 建立含有TMD的质量、刚度、阻尼矩阵
% Establish mass, stiffness, and damping matrices including TMD
input.MM_eq = MM_eq;
input.KK_eq = KK_eq;
input.eig_val = eig_val;
input.eig_vec = eig_vec;
input.modal_damping_ratios = modal_damping_ratios;
input.nodeondeck = nodeondeck;
input.mode = mode;
input.nodegap = nodegap;
input.modecal = modecal;
[output] = a_2_build_matrice_with_TMD(input);

MM = output.Mass;
KK = output.Stiffness;
CC = output.Damping;

% 获取气动力参数
% Apply aerodynamic force
girderindex = 1;
ExpNames = [
    'SZTD-110-case2-22.3-fasan-2401';
    'SZTD-110-case2-22.3-fasan-2501';
    'SZTD-110-case2-22.3-fasan-2601';
    'SZTD-110-case2-22.3-fasan-2701';
    'SZTD-110-case2-22.3-fasan-2801';
    'SZTD-110-case2-22.3-fasan-2901';
    'SZTD-110-case2-22.3-fasan-3101';
    ]; %记录文件名

for k1 = 1
    % 选择SZTD-110-case2-22.3-fasan-2401工况
    ExpName = ExpNames(k1, :);
end

input.ExpName = ExpName;
input.girderindex = girderindex;
[output] = a_3_VEF_parameter(input);
a = output.a;
a1 = a(1);
a2 = a(2);
a3 = a(3);
a4 = a(4);
a5 = a(5);
ReducedFrequency = output.ReducedFrequency;
rho = 1.225;

% 施加气动力并进行计算
% Apply aerodynamic force and calculate
matrixsize = number_of_tmds + number_of_modes_to_consider;
D = 20; %断面参考宽度 Cross-section Reference Width
h = 0.01; % Time step
gamma = 1/2; % Parameter in the Newmark algorithm
beta = 1/4; % Parameter in the Newmark algorithm
% Newmark-beta法中采用真实位移
% Use real displacement in Newmark-beta method

nModes = number_of_modes_to_consider;

for k1 = 1:nModes
    integral_1(k1) = sum(abs(modecal(:, k1)) .* nodedis);
    integral_2(k1) = sum(modecal(:, k1) .^ 2 .* nodedis);
    integral_3(k1) = sum(modecal(:, k1) .^ 2 .* abs(modecal(:, k1)) .* nodedis);
    integral_4(k1) = sum(modecal(:, k1) .^ 4 .* nodedis);
    integral_5(k1) = sum(modecal(:, k1) .^ 2 .* abs(modecal(:, k1)) .^ 3 .* nodedis);
    integral_6(k1) = sum(modecal(:, k1) .^ 6 .* nodedis);
end

modes_integral_1 = integral_1(number_of_modes_to_control);
modes_integral_2 = integral_2(number_of_modes_to_control);
modes_integral_3 = integral_3(number_of_modes_to_control);
modes_integral_4 = integral_4(number_of_modes_to_control);
modes_integral_5 = integral_5(number_of_modes_to_control);
modes_integral_6 = integral_6(number_of_modes_to_control);

dis_all_modes = zeros(number_of_modes_to_consider, 1);

for k1 = number_of_modes_to_control
    mode_number = k1; %气动力施加的模态
    % 预先计算的模态积分
    % Modal integral pre-calculated
    [~,k2]=find(number_of_modes_to_control==k1);
    mode_integral_1 = modes_integral_1(k2);
    mode_integral_2 = modes_integral_2(k2);
    mode_integral_3 = modes_integral_3(k2);
    mode_integral_4 = modes_integral_4(k2);
    mode_integral_5 = modes_integral_5(k2);
    mode_integral_6 = modes_integral_6(k2);

    % mode_number = k1; %气动力施加的模态
    % % 预先计算的模态积分
    % % Modal integral pre-calculated
    % mode_integral_1 = modes_integral_1(k1);
    % mode_integral_2 = modes_integral_2(k1);
    % mode_integral_3 = modes_integral_3(k1);
    % mode_integral_4 = modes_integral_4(k1);
    % mode_integral_5 = modes_integral_5(k1);
    % mode_integral_6 = modes_integral_6(k1);

    % 气动力参数

    for j = 1:nModes
        m_modal(j) = MM(j, j);
    end

    U = 2 * pi * Freq(k1) * D / ReducedFrequency; % 风速
    omega0 = 2 * pi * Freq(mode_number); %无风振动频率
    m = MM(mode_number, mode_number); %质量
    b1 = rho * U * D * a1 / m_modal(mode_number);
    b2 = rho * U * a2 / m_modal(mode_number);
    b3 = rho * U * a3 / D / m_modal(mode_number);
    b4 = rho * U * a4 / D ^ 2 / m_modal(mode_number);
    b5 = rho * U * a5 / D ^ 3 / m_modal(mode_number);

    gfun = @(u, udot) bridge_damper(u, udot, b1, b2, b3, b4, b5, MM, CC, KK, ...
        mode_integral_2, mode_integral_3, mode_integral_4, mode_integral_5, ...
        mode_integral_6, gamma, beta, h, matrixsize, mode_number); % Handle to the nonlinear function

    % TODO: 还未考虑气动刚度
    u0 = zeros(matrixsize, 1); % Initial displacement;
    udot0 = zeros(matrixsize, 1); % Initial velocity;
    u0max = 0.001; % Initial displacement of the first mode;
    u0(mode_number) = u0max / phideckmax(mode_number); % Initial displacement of the first mode;

    iter = 1;
    flag_convergence = 0;
    t_length_temp = t_length;
    iter_num_max = 3;

    while and(iter <= iter_num_max, flag_convergence == 0)
        t = 0:h: t_length_temp; % Time

        p = zeros(matrixsize, length(t)); %Initialize external load
        [u, udot, u2dot] = nonlinear_newmark_krenk(gfun, MM, p, u0, udot0, gamma, beta, h); % Solve the response by the Nonlinear Newmark algorithm
        %计算桥梁响应
        dis_temp = zeros(size(u, 2), nModes); % 分别表示时间长度和模态数量

        % for k01 = 1:nModes
        %     dis_temp(:, k01) = u(k01, :);
        % end
        dis_temp(:, 1:nModes) = u(1:nModes, :)';

        dis_TMD(:, 1:number_of_tmds) = u(nModes + 1:nModes + number_of_tmds, :)';

        npoints = size(eig_vec, 1);
        seq_length = length(u(1, :));
        seqs_allmodes = zeros(npoints, seq_length);

        for k01 = 1:nModes
            seqs_allmodes = seqs_allmodes + eig_vec(:, k01) .* u(k01, :);
        end
        seqs_allmodes_max = max(seqs_allmodes, [], 2);
        %寻找最大值和位置
        [max_value, max_index] = max(seqs_allmodes_max);
        seqs_allmodes_max_point = seqs_allmodes(max_index, :);
        seqs_allmodes = [];
        dis1 = max(seqs_allmodes_max_point(round(end / 8 * 6, 0):round(end / 8 * 7, 0)));
        dis2 = max(seqs_allmodes_max_point(round(end / 8 * 7, 0):round(end , 0)));


        if or(abs((dis1 - dis2) / dis2) < 0.02, dis2 < 1e-2)
            disp("计算收敛")
            flag_convergence = 1;

            %             seqs_all=[];
            dis_beam_max = std(seqs_allmodes_max_point(round(end / 8 * 6, 0):end)) * sqrt(2);

            %                 dis_TMDs_max = std(dis_TMD(round(end / 8 * 6, 0):end)) * sqrt(2);

            if iter >= iter_num_max
                flag_iter = 0;
            else
                flag_iter = 1;
            end

            dis_all_modes(k1) = dis_beam_max;
            clear dis_TMD dis_temp
        else
            disp("计算未收敛，增加计算时间")
            iter = iter + 1;
            t_length_temp = t_length_temp+450;
            clear dis_TMD dis_temp
        end
        dis_beam_max = std(seqs_allmodes_max_point(round(end / 8 * 6, 0):end)) * sqrt(2);
        dis_all_modes(k1) = dis_beam_max;
    end
    %          figure
    %          plot(t,seqs_allmodes_max_point)
    %          close all
end

dis_all_modes_sum = sum(dis_all_modes);
disp("多个模态振动位移之和" + dis_all_modes_sum + "m")
result.dis_all_modes_sum = dis_all_modes_sum;
result.dis_all_modes=dis_all_modes;
result.t = t;
result.u = u;
result.udot = udot;
result.u2dot = u2dot;

F_viv = rho*U^2*D*(a1*modes_integral_2+a2*modes_integral_3*abs(u/D)+a3*modes_integral_4*(u/D).^2+a4*modes_integral_5*abs(u/D).^3+a5*modes_integral_6*(u/D).^4).*udot/U;
result.F_viv = F_viv;


 [u_re udot_re u2dot_re] = NewmarkInt(t,MM,CC,KK,F_viv,1/2,1/4,u0,udot0);
result.u_re = u_re;
result.udot_re = udot_re;
result.u2dot_re = u2dot_re;
end
