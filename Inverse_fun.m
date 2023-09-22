%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com Date:
%LastEditors: ShengyiXu xushengyichn@outlook.com
%LastEditTime: 2023-07-29 22:42:26
%Description: 计算简支梁施加荷载后的动力响应，并反算出荷载（分别按照集中力和模态力反算）
%
%Copyright (c) 2023 by ${git_name_email}, All Rights Reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function result = Inverse_fun(lambdas_m,sigma_ps_m,modesel)
%% 1 有限元模型
% 读入ANSYS梁桥模型质量刚度矩阵  MCK矩阵 Import MCK matrix from ANSYS
% 将ANSYS中的稀疏矩阵处理为完全矩阵 Handling sparse matrices in ANSYS as full matrices
% 0 number of modes to consider

modesel= modesel;
% modesel= [23];
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




%% 3 传感器布置
% establish continuous time matrices
phi = mode_vec; %模态向量 每一列是一个模态
C = CC_eq; K = KK_eq; M = MM_eq;

Gamma = C; % 对应文献中表述的符号
omega2 = K;

% accelerometer location
loc_acc = [578+1650/4;578+1650/2;578+1650/4*3];
loc_vel = [];
loc_dis = [];
[S_a, S_v, S_d, n_sensors] = sensor_selection(loc_acc, loc_vel, loc_dis, node_loc, phi,nodeondeck,Mapping_data);
% establish continuous time matrices
[A_c, B_c, G_c, J_c] = ssmod_c_mode(nmodes, omega2, Gamma, phi, S_a, S_v, S_d);

% 生成带噪音的观测数据
Q =10 ^ (-8) * eye(ns);
R =10 ^ (-6) * eye(n_sensors);
S = zeros(ns, n_sensors);
x0 = zeros(ns, 1);


%% vibration
% vibac2_name=["2013-02-06 00-vibac2.txt";"2013-02-06 01-vibac2.txt";"2013-02-06 02-vibac2.txt"];
% vibac3_name=["2013-02-06 00-VIBac3.txt";"2013-02-06 01-VIBac3.txt";"2013-02-06 02-VIBac3.txt"];
% vibac4_name=["2013-02-06 00-VIBac4.txt";"2013-02-06 01-VIBac4.txt";"2013-02-06 02-VIBac4.txt"];
acc_names=["主跨1/4","主跨1/2","主跨3/4"];
% vibac2_data=read_vib_data(vibac2_name);
% vibac3_data=read_vib_data(vibac3_name);
% vibac4_data=read_vib_data(vibac4_name);
vibration_data = importdata("vibration_data.mat");
vibac2_data = vibration_data.vibac2_data;
vibac3_data = vibration_data.vibac3_data;
vibac4_data = vibration_data.vibac4_data;


dt = vibac2_data(2,1)-vibac2_data(1,1);

[f, magnitude] = fft_transform(1/0.02,vibac2_data(:,2));
[~,seq_temp]= max(magnitude);
f_max=f(seq_temp);
disp("Frequency regarding the maximum amplitude:"+num2str(f_max))
% figure
% plot(f, magnitude)
omega_0 = 2 * pi * Freq;
yn(1,:)=vibac2_data(:,2);
yn(2,:)=vibac2_data(:,4);
yn(3,:)=vibac3_data(:,2);
yn(4,:)=vibac3_data(:,4);
yn(5,:)=vibac4_data(:,2);
yn(6,:)=vibac4_data(:,4);

fs = 1/dt;
% establish discrete time matrices

startDate = datetime(2013,2,6,0,0,0);
t = startDate:seconds(dt):startDate+(length(vibac2_data(:,1))-1)*seconds(dt);
N = length(t);
[A_d, B_d, G_d, J_d, ~] = ssmod_c2d(A_c, B_c, G_c, J_c, dt);
%% wind

% winds1_name =["wind_property_result_s1.txt"];
% winds3_name =["wind_property_result_s3.txt"];
% winds4_name =["wind_property_result_s4.txt"];
% winds1_data = read_wind_data(winds1_name);
% winds3_data = read_wind_data(winds3_name);
% winds4_data = read_wind_data(winds4_name);
% 
% % startDate = datetime(2013,2,6,0,0,0);
% endDate = datetime(2013,2,6,3,0,0);
% winds1_trim = winds1_data(winds1_data.DateTime >= startDate & winds1_data.DateTime <= endDate, :);
% winds3_trim = winds3_data(winds3_data.DateTime >= startDate & winds3_data.DateTime <= endDate, :);
% winds4_trim = winds4_data(winds4_data.DateTime >= startDate & winds4_data.DateTime <= endDate, :);
% 




%% 5 反算模态力
Q_xd = Q;
np_m = nmodes;
B_c_m = [zeros(nmodes, np_m); ...
    eye(np_m, np_m)];
J_c_m = [S_a * phi];
B_d_m = A_c \ (A_d - eye(size(A_d))) * B_c_m;
J_d_m = J_c_m;

lambdas_m = lambdas_m;

sigma_ps_m = sigma_ps_m;

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
[x_k_k, x_k_kmin, P_k_k, P_k_kmin,result] = KalmanFilterNoInput(A_a_m, G_a_m, Q_a_m, R_a_m, yn_a, x_ak, P_ak, 'debugstate', true);
% % [x_k_k, P_k_k] = RTSFixedInterval(A_a_m, x_k_k, x_k_kmin, P_k_k, P_k_kmin);
% xa_history = x_k_k;
% pa_history = P_k_k;
% x_filt_original = xa_history(1:ns, :);
% H_d_m = H_c_m;
% p_filt_m = H_d_m * xa_history(ns + 1:end, :);
% Pp_filt_m = H_d_m * pa_history(ns + 1:end, :);
% 
% % % 滤波
% % f_low = 0.95 * 0.328241;
% % f_high = 1.05 * 0.328241;
% % p_filt_m(1, :) = bandpass(p_filt_m(1, :), [f_low f_high], fs);
% % p_filt_m(2, :) = bandpass(p_filt_m(2, :), [f_low f_high], fs);
% % p_filt_m(3, :) = bandpass(p_filt_m(3, :), [f_low f_high], fs);
% % p_filt_m(9, :) = bandpass(p_filt_m(9, :), [f_low f_high], fs);
% % 分析预测模态力频率
% for k1 = 1:nmodes
%     [f_p_filt_m(k1, :), magnitude_filt_m(k1, :)] = fft_transform(fs, p_filt_m(k1, :));
% end
% 
% 
% %% 6 virtual sensoring
% loc_acc_v = [578+1650/4;578+1650/2;578+1650/4*3];
% loc_vel_v = [];
% loc_dis_v = [];
% 
% [S_a_v, S_v_v, S_d_v, n_sensors_v] = sensor_selection(loc_acc, loc_vel, loc_dis, node_loc, phi,nodeondeck,Mapping_data);
% 
% G_c_v = [S_d_v * phi - S_a_v * phi * omega2, S_v_v * phi - S_a_v * phi * Gamma];
% J_c_v = [S_a_v * phi ];
% 
% h_hat = G_c_v * x_filt_original + J_c_v * p_filt_m;
% 
% %% 7 重构涡振响应
% p_reconstruct = p_filt_m;
% % p_reconstruct([1 2 3 4 5 6 7 8 10 11 12 13 14],:)=0;
% [~, yn_reconstruct, ~] = CalResponse(A_d, B_d, G_d, J_d, p_reconstruct, 0, 0, N, x0, ns, n_sensors);
% 
% %% 8 marginal likelihood
% 
% logL=result.logL;
% logSk=result.logSk;
% logek=result.logek;





%% plot

if 0
    for k1 = 1:length(acc_names)
    [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);

    % subplot(1, nmodes, k1)
    plot(t, yn(2*k1-1, :), 'Color', 'r')
    hold on
    plot(t, yn(2*k1, :), 'Color', 'b')

    xlabel('time (s)')
    ylabel('acc (1e-3*g)')
    set(gca, 'FontSize', 12)
    legend('filtered modal force', 'Location', 'northwest')
    title([acc_names(k1)]);
    ylim([-75, 75])
    end

    for k1 = 1:length(acc_names)
        [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
    
        % subplot(1, nmodes, k1)
        plot(t, yn_reconstruct(2*k1-1, :), 'Color', 'r')
        hold on
        plot(t, yn_reconstruct(2*k1, :), 'Color', 'b')
    
        xlabel('time (s)')
        ylabel('acc (1e-3*g)')
        set(gca, 'FontSize', 12)
        legend('filtered modal force', 'Location', 'northwest')
        title([acc_names(k1)]+"重构");
        ylim([-75, 75])
    end

    for k1 = 1:nmodes
    [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);

        % subplot(1, nmodes, k1)
        plot(t, p_filt_m(k1, :), 'Color', 'r')
        hold on
        %         plot(t, p_m_real(k1, :), 'Color', 'b', 'LineStyle', '--')
        xlabel('time (s)')
        ylabel('Modal force (N)')
        set(gca, 'FontSize', 12)
        legend('filtered modal force', 'Location', 'northwest')
        title(['mode ', num2str(modesel(k1))]);
        ylim([-1e5, 1e5])
    end

    set(hFigure, 'name', 'modal force estimation', 'Numbertitle', 'off');

    for k1 = 1:nmodes
    [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
        % subplot(1, nmodes, k1)
        plot(f_p_filt_m(k1, :), magnitude_filt_m(k1, :), 'Color', 'r')
        hold on
        %         plot(f_p_real_m(k1, :), magnitude_real_m(k1, :), 'Color', 'b', 'LineStyle', '--')
        xlabel('frequency (Hz)')
        ylabel('magnitude (N)')
        set(gca, 'FontSize', 12)
        legend('filtered modal force frequency', 'Location', 'northwest')
        title(['mode ', num2str(modesel(k1))]);
        xlim([0, 3])
        % ylim([0, 50])
    end

    set(hFigure, 'name', 'filtered modal force frequency', 'Numbertitle', 'off');
end



% 
% [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
% 
% for k1 = 1
%     plot(t, p_filt_m(k1, :), 'Color', 'r')
%     hold on
%     %         plot(t, p_m_real(k1, :), 'Color', 'b', 'LineStyle', '--')
%     xlabel('time (s)')
%     ylabel('Modal force (N)')
%     set(gca, 'FontSize', 12)
%     legend('filtered modal force', 'Location', 'northwest')
%     title(['mode ', num2str(modesel(k1))]);
% end
% 
% set(hFigure, 'name', 'filtered modal force frequency', 'Numbertitle', 'off');
% 
% 
% 
% [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
% 
% for k1 = 1
% 
%     plot(f_p_filt_m(k1, :), magnitude_filt_m(k1, :), 'Color', 'r')
%     hold on
%     %         plot(f_p_real_m(k1, :), magnitude_real_m(k1, :), 'Color', 'b', 'LineStyle', '--')
%     xlabel('frequency (Hz)')
%     ylabel('magnitude (N)')
%     set(gca, 'FontSize', 12)
%     legend('filtered modal force frequency', 'Location', 'northwest')
%     title(['mode ', num2str(modesel(k1))]);
%     xlim([0, 3])
%     % ylim([0, 50])
% end
% 
% set(hFigure, 'name', 'filtered modal force frequency', 'Numbertitle', 'off');
% 
% 




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
    acc_node_list=reshape(permute(acc_node , [2 1]), [], 1);  % 交错重塑
    % acc_node_list=acc_node(:);
else
    acc_node_list=[];
end

if ~isempty(loc_vel)
    vel_node = FindNodewithLocation(loc_vel, node_loc,nodeondeck);
    vel_node_list=reshape(permute(vel_node , [2 1]), [], 1);
else
    vel_node_list=[];
end

if ~isempty(loc_dis)
    dis_node = FindNodewithLocation(loc_dis, node_loc,nodeondeck);
    dis_node_list = reshape(permute(dis_node , [2 1]), [], 1);
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

function x = modal2physical_node(xn, phi, node_list,KMmapping)
% 节点在matrix中的编号
matrixseq=node2matrixseq(node_list,KMmapping);
phi_trim = phi(matrixseq,:);
x = [phi_trim, zeros(size(phi_trim)); zeros(size(phi_trim)), phi_trim] * xn;
end

function [F_c, L_c, H_c, sigma_w12] = ssmod_quasiperiod_coninue(lambdas, sigma_ps, omega_0, np)
% create arrays for latent force model
F_c_array = zeros(2, 2, np);
L_c_array = zeros(2, 2, np);
H_c_array = zeros(1, 2, np);
sigma_w = zeros(np, 1);

for k1 = 1:np
    [F_c_array(:, :, k1), L_c_array(:, :, k1), H_c_array(:, :, k1), sigma_w(k1)] = ssmod_quasiperiod(lambdas(k1), sigma_ps(k1), omega_0(k1));
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

function [A_c, B_c, G_c, J_c] = ssmod_c_mode(nmodes, omega2, Gamma, phi, S_a, S_v, S_d)
A_c = [zeros(nmodes), eye(nmodes); ...
    -omega2, -Gamma];
B_c = [zeros(nmodes, nmodes); ...
    eye(nmodes,nmodes)];
G_c = [S_d * phi - S_a * phi * omega2, ...
    S_v * phi - S_a * phi * Gamma];
J_c = [S_a * phi ];
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
indices=[];
for k1 = 1:length(node_list)
    indices(k1) = find(UYNode==node_list(k1));
end
if isempty(indices)
    matrix_seq=[];
else
    matrix_seq = UYMatrix(indices);
end
end

function [vibac_data]=read_vib_data(vibac_name)
vibac_data=[];
for i=1:length(vibac_name)
    vibac_data_temp=importdata(vibac_name(i));
    vibac_data=[vibac_data;vibac_data_temp];
end
end

function data = read_wind_data(wind_name)
data = readtable(wind_name);
data.DateTime = datetime(data.day, 'InputFormat', 'yyyy-MM-dd') + hours(data.hour) + minutes(data.x10min*10);
data.DateTime.Format = 'yyyy-MM-dd HH:mm';
data.day = [];
data.hour = [];
data.x10min = [];
data = movevars(data, 'DateTime', 'Before', data.Properties.VariableNames{1});
end
