%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com Date:
%LastEditors: xushengyichn xushengyichn@outlook.com
%LastEditTime: 2023-07-24 01:55:19
%Description: 计算简支梁施加荷载后的动力响应，并反算出荷载（分别按照集中力和模态力反算）
%
%Copyright (c) 2023 by ${git_name_email}, All Rights Reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
addpath(genpath("F:\git\ssm_tools"))
addpath(genpath("D:\Users\xushe\Documents\GitHub\ssm_tools"))
addpath(genpath("D:\Users\xushe\Documents\GitHub\Function_shengyi_package"))
addpath(genpath("F:\git\Function_shengyi_package"))
subStreamNumberDefault = 2132;
run('InitScript.m');

%% 0 绘图参数
fig_bool = ON;
num_figs_in_row = 5; %每一行显示几个图
figPos = figPosSmall; %图的大小，参数基于InitScript.m中的设置
%设置图片间隔
gap_between_images = [0, 0];
figureIdx = 0;

%% 1 有限元模型
% 读入ANSYS梁桥模型质量刚度矩阵  MCK矩阵 Import MCK matrix from ANSYS
% 将ANSYS中的稀疏矩阵处理为完全矩阵 Handling sparse matrices in ANSYS as full matrices
% 0 number of modes to consider

modesel= [2,3,5,6,7,9,15,21,23,29,33,39,44,45];
nmodes = length(modesel);
ns = nmodes * 2;
Result = ImportMK(nmodes, 'KMatrix.matrix', 'MMatrix.matrix', 'nodeondeck.txt', 'KMatrix.mapping', 'nodegap.txt','modesel',modesel);
mode_deck = Result.mode_deck;
mode_deck_re = Result.mode_deck_re;
node_loc = Result.node_loc;
Freq = Result.Freq;
MM_eq = Result.MM_eq;
KK_eq = Result.KK_eq;
mode_vec = Result.mode_vec;
nodeondeck = Result.nodeondeck;
Mapping_data = Result.Mapping;
% CC_eq = 0.1 * MM_eq + 0.005 * KK_eq; %人为指定瑞利阻尼
zeta = ones(size(modesel))*0.3/100;
omega = diag(2*pi*Freq);
CC_eq = 2.*MM_eq.*omega.*zeta;

%% 2 生成计算荷载
T = 300; dt = 0.01; fs = 1 / dt; t = 0:dt:T;

% Create force
lambda = 0.001; sigma_p_2 = 1; sigma_p = sqrt(sigma_p_2); f1 = Freq(1); A1 = 1000000;
[P] = quasiperiod_force(lambda, sigma_p, f1, dt, t);


% 更简化的荷载测试用代码
% A1=10;f1 = Freq(1);phi1=0;
% p = A1 * sin(2 * pi * f1 * t + phi1)*10;

%荷载施加位置
loc_p = [1500];



% 确定荷载位置和对应振型大小
p_node = FindNodewithLocation(loc_p, node_loc,nodeondeck);
p_node_list = p_node(:);
p_matrix_seq = node2matrixseq(p_node_list,Mapping_data);
% 计算荷载数量
np  = size(p_node(:) , 1);

p = A1 * P;
p = repmat(p,np,1);
%% 3 响应计算
% establish continuous time matrices
phi = mode_vec; %模态向量 每一列是一个模态
C = CC_eq; K = KK_eq; M = MM_eq;

Gamma = C; % 对应文献中表述的符号
omega2 = K;

S_p = zeros(size(phi, 1), np);


for k1 = 1:np
    S_p(p_node_list(k1), k1) = 1;
end

% accelerometer location
% loc_acc = [1000,1500,1700];
loc_acc = [1000,1500,1700];
loc_vel = [1000,1500,1700];
loc_dis = [1000,1500,1700];
[S_a, S_v, S_d, n_sensors] = sensor_selection(loc_acc, loc_vel, loc_dis, node_loc, phi,nodeondeck,Mapping_data);
% establish continuous time matrices
[A_c, B_c, G_c, J_c] = ssmod_c(nmodes, np, omega2, Gamma, phi, S_p, S_a, S_v, S_d);

% establish discrete time matrices
N = length(t);
[A_d, B_d, G_d, J_d, ~] = ssmod_c2d(A_c, B_c, G_c, J_c, dt);

% 生成带噪音的观测数据
Q =10 ^ (-8) * eye(ns);
R =10 ^ (-6) * eye(n_sensors);
S = zeros(ns, n_sensors);
x0 = zeros(ns, 1);
[xn, yn, xn_true] = CalResponse(A_d, B_d, G_d, J_d, p, Q, R, N, x0, ns, n_sensors);
x = modal2physical(xn, phi);
x_true = modal2physical(xn_true, phi);


% 绘制位移图
[X,Y,Z,nodes,elements]=plot_FEM('nodes.txt','elements.txt');

output = ImportMK(100, 'KMatrix.matrix', 'MMatrix.matrix', 'nodeondeck.txt', 'KMatrix.mapping', 'nodegap.txt');
node_modal_shape=[nodes.NODE, nodes.X, nodes.Y, nodes.Z];


% x_dis = std(x(1:end/2,:),0,2);
% x_dis_plot =x_dis/max(x_dis)*50;

x_dis = x(1:end/2,7001);
x_dis_plot=x_dis*5000;

[node_modal_shape_plot,nodes_displacement,nodes_displacementMagnitude,X_modal_shape,Y_modal_shape,Z_modal_shape]= update_node(node_modal_shape,x_dis_plot,Mapping_data,elements);



% 
% 

%% 4 知道荷载位置反算荷载大小
% Augmented Kalman Filter with latent force model
% y = yn;
% latent force model
lambdas = [0.01] * ones(1, np);
sigma_ps = [100] * ones(1, np);

omega_0 = 2 * pi * f1;
[F_c, L_c, H_c, sigma_w12] = ssmod_quasiperiod_coninue(lambdas, sigma_ps, omega_0, np);
Q_xd = Q;
[~, ~, ~, ~, Fad, ~, Had, ~, Qad] = ssmod_lfm_aug(A_c, B_c, G_c, J_c, F_c, H_c, L_c, Q_xd, sigma_w12, dt);
%TODO:Q_d calculation
A_a = Fad;
G_a = Had;
Q_a = Qad;
R_a = R;
yn_a = yn;


x_ak = zeros(ns + np * (2), 1);
P_ak = 10 ^ (1) * eye(ns + np * (2));
[x_k_k, x_k_kmin, P_k_k, P_k_kmin] = KalmanFilterNoInput(A_a, G_a, Q_a, R_a, yn_a, x_ak, P_ak);
% [x_k_k, P_k_k] = RTSFixedInterval(A_a, x_k_k, x_k_kmin, P_k_k, P_k_kmin);
xa_history = x_k_k;
pa_history = P_k_k;

x_filt_original = xa_history(1:ns, :);
x_filt = modal2physical(xa_history(1:ns, :), phi);
% Px_filt = abs([phi, zeros(size(phi)); zeros(size(phi)), phi]) * pa_history(1:ns, :);
H_d = H_c;
p_filt = H_d * xa_history(ns + 1:end, :);
Pp_filt = H_d * pa_history(ns + 1:end, :);

% 滤波
f_low = 0.95 * f1;
f_high = 1.05 * f1;
p_filt = bandpass(p_filt, [f_low f_high], fs);

% 分析预测荷载频率
fs = 1 / dt;

for k1 = 1:np
    [f_filt(k1, :), magnitude_filt(k1, :)] = fft_transform(fs, p_filt(k1, :));
end

for k1 = 1:np
    [f_true(k1, :), magnitude_true(k1, :)] = fft_transform(fs, p(k1, :));
end

%% 5 反算模态力
np_m = nmodes;
B_c_m = [zeros(nmodes, np_m); ...
             eye(np_m, np_m)];
J_c_m = [S_a * phi];
B_d_m = A_c \ (A_d - eye(size(A_d))) * B_c_m;
J_d_m = J_c_m;

lambdas_m = [0.01] * ones(1, np_m);

sigma_ps_m = [100] * ones(1, np_m);

[F_c_m, L_c_m, H_c_m, sigma_w_m12] = ssmod_quasiperiod_coninue(lambdas_m, sigma_ps_m, omega_0, np_m);

[~, ~, ~, ~, Fad_m, ~, Had_m, ~, Qad_m] = ssmod_lfm_aug(A_c, B_c_m, G_c, J_c_m, F_c_m, H_c_m, L_c_m, Q_xd, sigma_w_m12, dt);
A_a_m = Fad_m;
G_a_m = Had_m;
Q_a_m = Qad_m;
R_a_m = R;
yn_a = yn;
NN = N;
xa_history = zeros(ns + np_m * (2), NN);
pa_history = zeros(ns + np_m * (2), NN);

x_ak = zeros(ns + np_m * (2), 1);
P_ak = 10 ^ (1) * eye(ns + np_m * (2));

% G_a=G_a_m; A_a=A_a_m; Q_a=Q_a_m;
[x_k_k, x_k_kmin, P_k_k, P_k_kmin] = KalmanFilterNoInput(A_a_m, G_a_m, Q_a_m, R_a_m, yn_a, x_ak, P_ak);
% [x_k_k, P_k_k] = RTSFixedInterval(A_a_m, x_k_k, x_k_kmin, P_k_k, P_k_kmin);
xa_history = x_k_k;
pa_history = P_k_k;

H_d_m = H_c_m;
p_filt_m = H_d_m * xa_history(ns + 1:end, :);
Pp_filt_m = H_d_m * pa_history(ns + 1:end, :);

% 滤波
f_low = 0.95 * f1;
f_high = 1.05 * f1;
p_filt_m(1, :) = bandpass(p_filt_m(1, :), [f_low f_high], fs);
p_filt_m(2, :) = bandpass(p_filt_m(2, :), [f_low f_high], fs);
p_filt_m(3, :) = bandpass(p_filt_m(3, :), [f_low f_high], fs);

% 分析预测模态力频率
for k1 = 1:nmodes
    [f_p_filt_m(k1, :), magnitude_filt_m(k1, :)] = fft_transform(fs, p_filt_m(k1, :));
end

p_m_real = phi' * S_p * p;

for k1 = 1:nmodes
    [f_p_real_m(k1, :), magnitude_real_m(k1, :)] = fft_transform(fs, p_m_real(k1, :));
end

%% 6 virtual sensoring
loc_acc_v = [25; 50; 75];
loc_vel_v = [25; 50; 75];
loc_dis_v = [25; 50; 75];

[S_a_v, S_v_v, S_d_v, n_sensors_v] = sensor_selection(loc_acc, loc_vel, loc_dis, node_loc, phi,nodeondeck,Mapping_data);

G_c_v = [S_d_v * phi - S_a_v * phi * omega2, S_v_v * phi - S_a_v * phi * Gamma];
J_c_v = [S_a_v * phi * phi' * S_p];

h_hat = G_c_v * x_filt_original + J_c_v * p_filt;
h_hat_true = G_c_v * xn_true + J_c_v * p;

%% plot

if fig_bool == ON


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
    scatterHandle = scatter3(node_modal_shape_plot(:,2), node_modal_shape_plot(:,3), node_modal_shape_plot(:,4), size, nodes_displacementMagnitude, marker, 'filled');
    
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


%     for k1 = 1:nmodes
%         [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
%         plot(node_loc,mode_deck_re(:,k1),'Color', 'r', 'LineWidth', lineWidthThin)
%         xlabel('location(m)')
%         ylabel('modal_shape')
%     end
    


    [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
    subplot(1, 3, 1)
    plot(t, h_hat(1, :), 'Color', 'r', 'LineWidth', lineWidthThin);
    hold on
    plot(t, h_hat_true(1, :), 'Color', 'b', 'LineStyle', '--', 'LineWidth', lineWidthThin);
    plot(t, yn(1, :), 'Color', 'g', 'LineStyle', '-.', 'LineWidth', lineWidthThin);
    set(gca, 'FontSize', 12)
    xlabel('time (s)')
    ylabel('displacement (m)')
    legend('filtered', 'true', 'measure', 'Location', 'northwest')
    title("acc at "+num2str(loc_acc_v(1)));
    subplot(1, 3, 2)
    hLineObj = plot(t, h_hat(2, :), 'Color', 'r');
    set(gca, 'FontSize', 12)
    hold on
    hLineObj = plot(t, h_hat_true(2, :), 'Color', 'b', 'LineStyle', '--');
    set(gca, 'FontSize', 12)

    xlabel('time (s)')
    ylabel('displacement (m)')
    legend('filtered', 'true', 'Location', 'northwest')
    title("acc at "+num2str(loc_acc_v(2)));
    subplot(1, 3, 3)
    hLineObj = plot(t, h_hat(3, :), 'Color', 'r');
    hold on
    hLineObj = plot(t, h_hat_true(3, :), 'Color', 'b', 'LineStyle', '--');
    set(gca, 'FontSize', 12)
    xlabel('time (s)')
    ylabel('displacement (m)')
    legend('filtered', 'true', 'Location', 'northwest')
    title("acc at "+num2str(loc_acc_v(3)));
    set(hFigure, 'name', 'estimate acceleration', 'Numbertitle', 'off');

    [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
    subplot(1, 3, 1)
    hLineObj = plot(t, h_hat(7, :), 'Color', 'r');
    set(hLineObj, 'LineWidth', lineWidthThin);
    hold on
    hLineObj = plot(t, h_hat_true(7, :), 'Color', 'b', 'LineStyle', '--');
    set(hLineObj, 'LineWidth', lineWidthThin);
    xlabel('time (s)')
    ylabel('displacement (m)')
    legend('filtered', 'true', 'Location', 'northwest')
    title("dis at "+num2str(loc_dis_v(1)));
    subplot(1, 3, 2)
    hLineObj = plot(t, h_hat(8, :), 'Color', 'r');
    hold on
    hLineObj = plot(t, h_hat_true(8, :), 'Color', 'b', 'LineStyle', '--');
    set(gca, 'FontSize', 12)
    xlabel('time (s)')
    ylabel('displacement (m)')
    legend('filtered', 'true', 'Location', 'northwest')
    title("dis at "+num2str(loc_dis_v(2)));
    subplot(1, 3, 3)
    hLineObj = plot(t, h_hat(9, :), 'Color', 'r');
    hold on
    hLineObj = plot(t, h_hat_true(9, :), 'Color', 'b', 'LineStyle', '--');
    set(gca, 'FontSize', 12)
    xlabel('time (s)')
    ylabel('displacement (m)')
    legend('filtered', 'true', 'Location', 'northwest')
    title("dis at "+num2str(loc_dis_v(3)));
    set(hFigure, 'name', 'estimate displacement', 'Numbertitle', 'off');

    [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);

    subplot(1, 3, 1)
    plot(t, h_hat(4, :), 'Color', 'r')
    hold on
    plot(t, h_hat_true(4, :), 'Color', 'b', 'LineStyle', '--')
    set(gca, 'FontSize', 12)
    xlabel('time (s)')
    ylabel('velocity (m/s)')
    legend('filtered', 'true', 'Location', 'northwest')
    title("vel at "+num2str(loc_vel_v(1)));
    subplot(1, 3, 2)
    plot(t, h_hat(5, :), 'Color', 'r')
    hold on
    plot(t, h_hat_true(5, :), 'Color', 'b', 'LineStyle', '--')
    set(gca, 'FontSize', 12)
    xlabel('time (s)')
    ylabel('velocity (m/s)')
    legend('filtered', 'true', 'Location', 'northwest')
    title("vel at "+num2str(loc_vel_v(2)));
    subplot(1, 3, 3)
    plot(t, h_hat(6, :), 'Color', 'r')
    hold on
    plot(t, h_hat_true(6, :), 'Color', 'b', 'LineStyle', '--')
    set(gca, 'FontSize', 12)
    xlabel('time (s)')
    ylabel('velocity (m/s)')
    legend('filtered', 'true', 'Location', 'northwest')
    title("vel at "+num2str(loc_vel_v(3)));
    set(hFigure, 'name', 'estimate velocity', 'Numbertitle', 'off');

    [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);

    for k1 = 1:np
        subplot(1, np, k1)
        plot(t, p_filt(k1, :), 'Color', 'r')
        hold on
        plot(t, p(k1, :), 'Color', 'b', 'LineStyle', '--')
        set(gca, 'FontSize', 12)
        t_fill = [t(1:NN), fliplr(t(1:NN))];
        f_upper_bound = p_filt(k1, :) + sqrt(Pp_filt(k1, k1, :));
        f_lower_bound = p_filt(k1, :) - sqrt(Pp_filt(k1, k1, :));
        f_fill = [f_upper_bound(1:NN), fliplr(f_lower_bound(1:NN))];
        fill(t_fill, f_fill, 'r', 'FaceAlpha', 0.3, 'EdgeColor', 'none')
        xlabel('time (s)')
        ylabel('force (N)')
        legend('filtered', 'true', 'Location', 'northwest')
        title(['p', num2str(k1)]);
        ylim([-100000, 100000])
    end

    set(hFigure, 'name', 'estimate force', 'Numbertitle', 'off');

    [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);

    for k1 = 1:np
        subplot(1, np, k1)
        plot(f_filt(k1, :), magnitude_filt(k1, :), 'Color', 'r')
        hold on
        plot(f_true(k1, :), magnitude_true(k1, :), 'Color', 'b', 'LineStyle', '--')
        set(gca, 'FontSize', 12)
        xlabel('frequency (Hz)')
        ylabel('magnitude (N)')
        legend('filtered', 'true', 'Location', 'northwest')
        title(['p', num2str(k1)]);
        xlim([0, 3])
        % ylim([0, 50])
    end

    set(hFigure, 'name', 'force frequency', 'Numbertitle', 'off');

    [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);

    for k1 = 1:nmodes
        subplot(1, nmodes, k1)
        plot(t, p_filt_m(k1, :), 'Color', 'r')
        hold on
        plot(t, p_m_real(k1, :), 'Color', 'b', 'LineStyle', '--')
        xlabel('time (s)')
        ylabel('Modal force (N)')
        set(gca, 'FontSize', 12)
        legend('filtered modal force', 'Location', 'northwest')
        title(['mode ', num2str(k1)]);
    end

    set(hFigure, 'name', 'modal force estimation', 'Numbertitle', 'off');

    [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);

    for k1 = 1:nmodes
        subplot(1, nmodes, k1)
        plot(f_p_filt_m(k1, :), magnitude_filt_m(k1, :), 'Color', 'r')
        hold on
        plot(f_p_real_m(k1, :), magnitude_real_m(k1, :), 'Color', 'b', 'LineStyle', '--')
        xlabel('frequency (Hz)')
        ylabel('magnitude (N)')
        set(gca, 'FontSize', 12)
        legend('filtered modal force frequency', 'Location', 'northwest')
        title(['mode ', num2str(k1)]);
        xlim([0, 3])
        % ylim([0, 50])
    end

    set(hFigure, 'name', 'filtered modal force frequency', 'Numbertitle', 'off');
end

function node = FindNodewithLocation(loc, node_loc,nodeondeck)
    %myFun - Description
    %
    % Syntax: result = findnode(input)
    %
    % Long description
    node_seq=zeros(1,length(loc));
    for k1 = 1:length(loc)
        node_seq(k1) = find(node_loc >= loc(k1), 1);
    end
    node = nodeondeck(node_seq,:);
end

function [S_a, S_v, S_d, n_sensors] = sensor_selection(loc_acc, loc_vel, loc_dis, node_loc, phi,nodeondeck,Mapping_data)

    if ~isempty(loc_acc)
        acc_node = FindNodewithLocation(loc_acc, node_loc,nodeondeck);
        acc_node_list=acc_node(:);
    else
        acc_node_list=[];
    end

    if ~isempty(loc_vel)
        vel_node = FindNodewithLocation(loc_vel, node_loc,nodeondeck);
        vel_node_list=vel_node(:);
    else
        vel_node_list=[];
    end

    if ~isempty(loc_dis)
        dis_node = FindNodewithLocation(loc_dis, node_loc,nodeondeck);
        dis_node_list = dis_node(:);
    else
        dis_node_list=[];
    end

    n_acc = length(acc_node_list);
    n_vel = length(vel_node_list);
    n_dis = length(dis_node_list);
    acc_matrix_seq = node2matrixseq(acc_node_list,Mapping_data);
    vel_matrix_seq = node2matrixseq(vel_node_list,Mapping_data);
    dis_matrix_seq = node2matrixseq(dis_node_list,Mapping_data);

    n_sensors = n_acc + n_vel + n_dis;


    S_a = zeros(n_sensors, size(phi, 1));
    S_v = zeros(n_sensors, size(phi, 1));
    S_d = zeros(n_sensors, size(phi, 1));

    for k1 = 1:n_sensors

        if k1 <= n_acc
            S_a(k1, acc_matrix_seq(k1)) = 1;
        elseif k1 <= n_acc + n_vel
            S_v(k1, vel_matrix_seq(k1 - n_acc)) = 1;
        else
            S_d(k1, dis_matrix_seq(k1 - n_acc - n_vel)) = 1;
        end

    end

end

function [xn, yn, xn_true] = CalResponse(A_d, B_d, G_d, J_d, p, Q, R, N, x0, ns, n_sensors)
    w = sqrt(Q) * randn(ns, N);
    v = sqrt(R) * randn(n_sensors, N);
    xn = x0;
    xn_true = x0;

    for t1 = 1:N - 1
        xn(:, t1 + 1) = A_d * xn(:, t1) + B_d * p(:, t1) + w(:, t1);
        yn(:, t1) = G_d * xn(:, t1) + J_d * p(:, t1) + v(:, t1);
        xn_true(:, t1 + 1) = A_d * xn_true(:, t1) + B_d * p(:, t1);
    end

    yn(:, end + 1) = G_d * xn(:, end) + J_d * p(:, end) + v(:, end);
end

function x = modal2physical(xn, phi)
%     x = [phi, zeros(size(phi)); zeros(size(phi)), phi] * xn;
    xn_top = xn(1:end/2, :);
    xn_bottom = xn(end/2+1:end, :);
    
    x_part_top = phi * xn_top;
    x_part_bottom = phi * xn_bottom;
    
    x = [x_part_top; x_part_bottom];
end

function [F_c, L_c, H_c, sigma_w12] = ssmod_quasiperiod_coninue(lambdas, sigma_ps, omega_0, np)
    % create arrays for latent force model
    F_c_array = zeros(2, 2, np);
    L_c_array = zeros(2, 2, np);
    H_c_array = zeros(1, 2, np);
    sigma_w = zeros(np, 1);

    for k1 = 1:np
        [F_c_array(:, :, k1), L_c_array(:, :, k1), H_c_array(:, :, k1), sigma_w(k1)] = ssmod_quasiperiod(lambdas(k1), sigma_ps(k1), omega_0);
    end

    F_c = [];
    L_c = [];
    H_c = [];

    for k1 = 1:np
        F_c = blkdiag(F_c, F_c_array(:, :, k1));
        L_c = blkdiag(L_c, L_c_array(:, :, k1));
        H_c = blkdiag(H_c, H_c_array(:, :, k1));
    end

    sigma_w12 = [];

    for k1 = 1:np
        sigma_w12 = blkdiag(sigma_w12, sigma_w(k1) * eye(2));
    end

end

function [A_c, B_c, G_c, J_c] = ssmod_c(nmodes, np, omega2, Gamma, phi, S_p, S_a, S_v, S_d)
    A_c = [zeros(nmodes), eye(nmodes); ...
               -omega2, -Gamma];
    B_c = [zeros(nmodes, np); ...
               phi' * S_p];
    G_c = [S_d * phi - S_a * phi * omega2, ...
               S_v * phi - S_a * phi * Gamma];
    J_c = [S_a * phi * phi' * S_p];
end

function [P]=quasiperiod_force(lambda,sigma_p,f,dt,t)
    f1=f;
    omega_0=2*pi*f1;
    [F_c,L_c,H_c,sigma_w] = ssmod_quasiperiod(lambda, sigma_p, omega_0);
    sigma_w_2 = sigma_w ^ 2;
    F_d = expm(F_c * dt);
    H_d = H_c;
    Qwc = diag([sigma_w_2]);
    Qwd = dt * L_c * Qwc * L_c';
    N = length(t);
    s = zeros(2, N);
    s(:,1) = zeros(2, 1);
    w = randn(1, N) .* sqrt(diag(Qwd));
    for k1 = 1:N
        s(:,k1 + 1) = F_d * s(:,k1) +w(:,k1);
        P(:,k1) = H_d * s(:,k1);
    end
    P=rescale(P)-0.5;
end

function matrix_seq=node2matrixseq(node_list,KMmapping)
    UYNode=KMmapping.Node(KMmapping.DOF == 'UY');
    UYMatrix=KMmapping.MatrixEqn(KMmapping.DOF == 'UY');
    for k1 = 1:length(node_list)
        indices(k1) = find(UYNode==node_list(k1));
    end
    matrix_seq = UYMatrix(indices);
end
