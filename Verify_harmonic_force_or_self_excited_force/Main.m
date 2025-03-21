%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%Date: 2023-01-25 15:24:39
%LastEditors: xushengyichn xushengyichn@outlook.com
%LastEditTime: 2023-12-21 21:16:30
%FilePath: /ssm_tools_sy/Users/xushengyi/Documents/GitHub/xihoumen_inverse_force_estimation/Verify_harmonic_force_or_self_excited_force/Main.m
%Description: 本脚本功能为输入一组TMD参数以及气动力控制参数，获得多阶模态涡振的响应
%
%Copyright (c) 2023 by xushengyichn 54436848+xushengyichn@users.noreply.github.com, All Rights Reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
close all
run("CommonCommand.m")
% addpath(genpath("F:\git\Function_shengyi_package"))
% addpath(genpath("F:\git\ssm_tools"))
% addpath(genpath("F:\git\HHT-Tutorial\HuangEMD"))
% addpath(genpath("D:\matlab\Mingjie杰哥参数识别程序\HHT-Tutorial-master"))
% subStreamNumberDefault = 2132;
% run('InitScript.m');
% %% 0 绘图参数
% fig_bool = ON;
% num_figs_in_row = 10; %每一行显示几个图
% figPos = figPosSmall; %图的大小，参数基于InitScript.m中的设置
% %设置图片间隔
% gap_between_images = [0, 0];
% figureIdx = 0;

%% 输入预先设定的参数
% Preset parameters

number_of_modes_to_control = [1]; % 需要控制的模态数 The number of modes to be controlled
number_of_modes_to_consider = 1; % 考虑的总模态数 The total number of modes considered
number_of_tmds = 0; % 考虑的总TMD数 The total number of TMDs
total_tmd_mass_ratio = 0.02; % 总质量比 The total mass ratio
mass_six_span = 10007779.7; % 深中通道非通航桥六跨连续梁质量 The mass of 6-span continuous beam of the non-navigational bridge of the Zhenzhong-Link
total_tmd_mass = total_tmd_mass_ratio * mass_six_span; % 总质量 The total mass
modal_damping_ratios = ones(1, number_of_modes_to_consider) * 0.003; % 模态阻尼 The damping of modes
t_length = 100; % 设定计算时间长度

% 输入初始参数
% Input initial parameters
TMDs_mass = 0; % TMD质量 The mass of TMDs
TMDs_frequency = 0; % TMD频率 The frequency of TMDs
TMDs_damping_ratio = 0; % TMD阻尼 The damping of TMDs
TMDs_location = 0; % TMD位置 The location of TMDs

[result] = a_0_main(number_of_modes_to_control, number_of_modes_to_consider, number_of_tmds, modal_damping_ratios, t_length, TMDs_mass, TMDs_frequency, TMDs_damping_ratio, TMDs_location);
t=result.t;
dt = t(2)-t(1);
u=result.u;
udot=result.udot;
u2dot=result.u2dot;
F_viv=result.F_viv;

u_re=result.u_re;
udot_re=result.udot_re;
u2dot_re=result.u2dot_re;


u02=u(end);
udot02=udot(end);
F_viv2 = zeros(size(F_viv));
MM = result.MM;
CC = result.CC;
KK = result.KK;
[u_2,udot_2,u2dot_2] = NewmarkInt(t,MM,CC,KK,F_viv2,1/2,1/4,u02,udot02);

u = [u,u_2];
udot=[udot,udot_2];
u2dot_2=[u2dot,u2dot_2];

u_re=[u_re,u_2];
udot_re=[udot_re,udot_2];
u2dot_re=[u2dot_re,u2dot_2];
F_viv = [F_viv,F_viv2];
t = 0:0.01: 200+0.01; % Time
[f, magnitude] = fft_transform(100,u);
[f_re, magnitude_re] = fft_transform(100,u_re);



%% Kalman filter
modesel= [1];
nmodes = length(modesel);
ns = nmodes * 2;
Result = ImportMK(nmodes, 'KMatrix.matrix', 'MMatrix.matrix', 'nodeondeck.txt', 'KMatrix.mapping', 'nodegap.txt','modesel',modesel);
mode_deck = Result.mode_deck;
mode_deck_re = Result.mode_deck_re;
node_loc = Result.node_loc;
Freq = Result.Freq;
MM_eq = Result.MM_eq;
KK_eq = Result.KK_eq;
nodegap = Result.nodegap;

mode_vec = Result.mode_vec;
nodeondeck = Result.nodeondeck;
Mapping_data = Result.Mapping;
% CC_eq = 0.1 * MM_eq + 0.005 * KK_eq; %人为指定瑞利阻尼
zeta = modal_damping_ratios;
omega = diag(2*pi*Freq);
CC_eq = 2.*MM_eq.*omega.*zeta;

phi = mode_vec; %模态向量 每一列是一个模态
C = CC_eq; K = KK_eq; M = MM_eq;


Gamma = C; % 对应文献中表述的符号
omega2 = K;

loc_acc = [660/4*1;660/4*2;660/4*3];
loc_vel = [];
loc_dis = [];
[S_a, S_v, S_d, n_sensors] = sensor_selection(loc_acc, loc_vel, loc_dis, node_loc, phi,nodeondeck,Mapping_data);
% establish continuous time matrices
[A_c, B_c, G_c, J_c] = ssmod_c_mode(nmodes, omega2, Gamma, phi, S_a, S_v, S_d);

[A_d, B_d, G_d, J_d, ~] = ssmod_c2d(A_c, B_c, G_c, J_c, dt);
% 生成带噪音的观测数据
Q =10 ^ (-8) * eye(ns);
R =10 ^ (-5) * eye(n_sensors);
S = zeros(ns, n_sensors);
x0 = zeros(ns, 1);
N = length(t);
%% calculate response

[xn_true, yn_true,~] = CalResponse(A_d, B_d, G_d, J_d, F_viv, 0, 0, N, x0, ns, n_sensors);
[xn, yn,~] = CalResponse(A_d, B_d, G_d, J_d, F_viv, Q, R, N, x0, ns, n_sensors);


%% perform kalman filter
Q_xd = Q;
np_m = nmodes;
B_c_m = [zeros(nmodes, np_m); ...
    eye(np_m, np_m)];
J_c_m = [S_a * phi];
B_d_m = A_c \ (A_d - eye(size(A_d))) * B_c_m;
J_d_m = J_c_m;


lambda = 1e-03;
sigma_p = 10000;
lambdas_m = [lambda] * ones(1, np_m);
sigma_ps_m = [sigma_p] * ones(1, np_m);

for k1 = 1:nmodes
    omega_0(k1)=Freq(k1)*2*pi;
end
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

% % % % % % % G_a=G_a_m; A_a=A_a_m; Q_a=Q_a_m;
[x_k_k, x_k_kmin, P_k_k, P_k_kmin,result] = KalmanFilterNoInput(A_a_m, G_a_m, Q_a_m, R_a_m, yn_a, x_ak, P_ak, 'debugstate', true,'showtext',false);
xa_history = x_k_k;
pa_history = P_k_k;
x_filt_original = xa_history(1:ns, :);
H_d_m = H_c_m;
p_filt_m = H_d_m * xa_history(ns + 1:end, :);
Pp_filt_m = H_d_m * pa_history(ns + 1:end, :);

%% recalculate the response with estimated force
[xn_re_kalman, yn_re_kalman,~] = CalResponse(A_d, B_d, G_d, J_d, p_filt_m, Q, R, N, x0, ns, n_sensors);

[f_re_kalman, magnitude_re_kalman] = fft_transform(100,xn_re_kalman(1,:));

%% Calculate aerodynamic force work
fs = 1/0.01;
% t = 0:1/fs:100;

% 
[p,fd,td] = pspectrum(x_k_k(1,:),t,'spectrogram','FrequencyResolution',0.05);

[ifq1,t1] = instfreq(p,fd,td);
% [ifq2,t2] = instfreq(u_re,fs);



% 首先，我们对于t1范围内的t值执行插值
inside_indices = t >= min(t1) & t <= max(t1);
t_inside = t(inside_indices);
ifq1_interpolated_inside = interp1(t1, ifq1, t_inside, 'linear');

% 然后，我们对于t1范围外的t值保留最近的ifq1值
t_outside_left = t(t < min(t1));
t_outside_right = t(t > max(t1));
ifq1_extrapolated_left = repmat(ifq1(1), size(t_outside_left));
ifq1_extrapolated_right = repmat(ifq1(end), size(t_outside_right));

% 最后，我们将这些值合并到一个数组中
ifq1_interpolated = [ifq1_extrapolated_left  ifq1_interpolated_inside  ifq1_extrapolated_right];

% 现在，ifq1_interpolated 是一个与t具有相同长度的数组，其中包含通过插值和外插获得的ifq1值。



% x_analytic = hilbert(u_re);
% 
% inst_amplitude = abs(x_analytic);
% figure;
% plot(t, inst_amplitude);
% hold on
% plot(t,u_re)
% xlabel('Time (s)');
% ylabel('Instantaneous Amplitude');
% title('Instantaneous Amplitude vs Time');
% 


Fa = p_filt_m;
vel = x_k_k(2,:);
dis = x_k_k(1,:);
t = t;


% 找到峰值
d = 100;
p = 0;
[peaks, locs] = findpeaks(Fa,'MinPeakDistance', d,'MinPeakProminence', p);



total_plots = 16; % 或任何你需要的子图数量
current_plot = 1;
num_figs_in_row = [];
figWidthFactor = 1.5;
%     figPosition = [1080*2.5,100];
figPosition = [100,100];
newfigure = false;
holdon = true;

% 画图来验证峰值检测的准确性
create_subplot(@plot, total_plots, current_plot, {t, Fa}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor,'figPosition',figPosition,'newfigure',true,'firstfigure',true);
hold on;
create_subplot(@plot, total_plots, current_plot, {t(locs), peaks, 'ro'}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor,'figPosition',figPosition,'newfigure',newfigure,'firstfigure',false);
title('Peak detection');
xlabel('Time');
ylabel('Force');
current_plot = current_plot + 1;

for k1 = 1:length(locs)-1
    dt = t(locs(k1+1))-t(locs(k1));
    work(k1) = trapz(t(locs(k1):locs(k1+1)),Fa(locs(k1):locs(k1+1)).*vel(locs(k1):locs(k1+1)));
    dis_k1 = dis(locs(k1):locs(k1+1));
    amp(k1) = (max(dis_k1)-min(dis_k1))/2;
    f(k1) = (ifq1_interpolated(locs(k1))+ifq1_interpolated(locs(k1+1)))/2;
    omega(k1) = 2*pi*f(k1);
    c(k1) = work(k1)/pi/amp(k1)^2/omega(k1);
    zeta(k1)=c(k1)/2/omega(k1);
end
amp_filt_kalman=amp;
zeta_filt_kalman=zeta;
create_subplot(@plot, total_plots, current_plot, {amp_filt_kalman,zeta_filt_kalman}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor,'figPosition',figPosition,'newfigure',newfigure,'firstfigure',false);
current_plot = current_plot + 1;



Fa = F_viv;
vel = udot;
dis = u;
t = t;


% 找到峰值
d = 100;
p = 0;
[peaks, locs] = findpeaks(Fa,'MinPeakDistance', d,'MinPeakProminence', p);

% 画图来验证峰值检测的准确性
create_subplot(@plot, total_plots, current_plot, {t, Fa}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor,'figPosition',figPosition,'newfigure',newfigure,'firstfigure',false);
hold on;
create_subplot(@plot, total_plots, current_plot, {t(locs), peaks, 'ro'}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor,'figPosition',figPosition,'newfigure',newfigure,'firstfigure',false);
title('Peak detection');
xlabel('Time');
ylabel('Force');
current_plot = current_plot + 1;


for k1 = 1:length(locs)-1
    dt = t(locs(k1+1))-t(locs(k1));
    work(k1) = trapz(t(locs(k1):locs(k1+1)),Fa(locs(k1):locs(k1+1)).*vel(locs(k1):locs(k1+1)));
    dis_k1 = dis(locs(k1):locs(k1+1));
    amp(k1) = (max(dis_k1)-min(dis_k1))/2;
    f(k1) = 1/dt;
    omega(k1) = 2*pi*f(k1);
    c(k1) = work(k1)/pi/amp(k1)^2/omega(k1);
    zeta(k1)=c(k1)/2/omega(k1);
end
amp_real=amp;
zeta_real=zeta;


create_subplot(@plot, total_plots, current_plot, {amp_real,zeta_real}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor,'figPosition',figPosition,'newfigure',newfigure,'firstfigure',false);
current_plot = current_plot + 1;


%% Calculate aerodynamic damping ratio
% zeta_polynomial = F_viv./udot./2./MM_eq./(2*pi*Freq).^2;
% zeta_filt_kalman = p_filt_m./udot./2./MM_eq./(2*pi*Freq).^2;

%% Calculate the damping ratio based on displacement 




[ex, frex] = ee(x_k_k(1,:), 1 / fs); %经验包络法求瞬时频率和瞬时振幅 % empirical envelope method to find instantaneous frequency and instantaneous amplitude
% [ex, frex] = ee(u_re, 1 / fs); %经验包络法求瞬时频率和瞬时振幅 % empirical envelope method to find instantaneous frequency and instantaneous amplitude
% figure
% plot(t, ex);
% hold on
% plot(t,x_k_k(1,:))
% xlabel('Time (s)');
% ylabel('Instantaneous Amplitude');
% title('Instantaneous Amplitude vs Time');



create_subplot(@plot, total_plots, current_plot, {t,frex}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor,'figPosition',figPosition,'newfigure',newfigure,'firstfigure',false);
current_plot = current_plot + 1;

omgx = frex*2*pi;
for k1 = 1:length(ex)-1
    epsx(k1)=log(ex(k1)/ex(k1+1))/omgx(k1)*fs;
end

epsx = [epsx epsx(end)];


create_subplot(@plot, total_plots, current_plot, {t,epsx}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor,'figPosition',figPosition,'newfigure',newfigure,'firstfigure',false);
grid
current_plot = current_plot + 1;


create_subplot(@plot, total_plots, current_plot, {x_k_k(1,:),x_k_k(2,:)}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor,'figPosition',figPosition,'newfigure',newfigure,'firstfigure',false);
grid
current_plot = current_plot + 1;




%% plot




create_subplot(@plot, total_plots, current_plot, {t, u, 'LineWidth', 2}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor,'figPosition',figPosition,'newfigure',newfigure,'firstfigure',false);
hold on
create_subplot(@plot, total_plots, current_plot, {t, u_re, 'LineWidth', 2}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor,'figPosition',figPosition,'newfigure',newfigure,'firstfigure',false);
xlabel('Time (s)')
ylabel('Modal Displacement')
legend('True','Recalculated using VIV force')
current_plot = current_plot + 1;

create_subplot(@plot, total_plots, current_plot, {t, F_viv, 'LineWidth', 2}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor,'figPosition',figPosition,'newfigure',newfigure,'firstfigure',false);
hold on
create_subplot(@plot, total_plots, current_plot, {t, p_filt_m, 'LineWidth', 2}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor,'figPosition',figPosition,'newfigure',newfigure,'firstfigure',false);
xlabel('Time (s)')
ylabel('Modal VIV force')
legend('True','Estimated')
current_plot = current_plot + 1;

create_subplot(@plot, total_plots, current_plot, {f, log(magnitude), 'LineWidth', 2}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor,'figPosition',figPosition,'newfigure',newfigure,'firstfigure',false);
hold on 
create_subplot(@plot, total_plots, current_plot, {f_re, log(magnitude_re), 'LineWidth', 2}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor,'figPosition',figPosition,'newfigure',newfigure,'firstfigure',false);
create_subplot(@plot, total_plots, current_plot, {f_re_kalman, log(magnitude_re_kalman), 'LineWidth', 2}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor,'figPosition',figPosition,'newfigure',newfigure,'firstfigure',false);
xlabel('Frequency (Hz)')
ylabel('Log Magnitude')
legend('True','Recalculated using VIV force','Recalculated using VIV force estimated by Kalman Filter')
current_plot = current_plot + 1;

create_subplot(@plot, total_plots, current_plot, {nodegap,mode_deck, 'LineWidth', 2}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor,'figPosition',figPosition,'newfigure',newfigure,'firstfigure',false);
xlabel('Location (m)')
ylabel('Mode shape')
current_plot = current_plot + 1;


create_subplot(@plot, total_plots, current_plot, {t,xn_true(1,:), 'LineWidth', 2}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor,'figPosition',figPosition,'newfigure',newfigure,'firstfigure',false);
hold on
create_subplot(@plot, total_plots, current_plot, {t,xn(1,:), 'LineWidth', 2}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor,'figPosition',figPosition,'newfigure',newfigure,'firstfigure',false);
create_subplot(@plot, total_plots, current_plot, {t,xn_re_kalman(1,:),'LineWidth',2}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor,'figPosition',figPosition,'newfigure',newfigure,'firstfigure',false);
create_subplot(@plot, total_plots, current_plot, {t,xa_history(1,:),'LineWidth',2}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor,'figPosition',figPosition,'newfigure',newfigure,'firstfigure',false);
xlabel('Time (s)')
ylabel('Modal Displacement')
legend('True','Calculate using VIV force','Recalculated using VIV force estimated by Kalman filter','Kalman filter estimated')
current_plot = current_plot + 1;


create_subplot(@scatter, total_plots, current_plot, {ex,epsx,'green'}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor,'figPosition',figPosition,'newfigure',newfigure,'firstfigure',false);
hold on
create_subplot(@scatter, total_plots, current_plot, {amp_filt_kalman,-zeta_filt_kalman,'o'}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor,'figPosition',figPosition,'newfigure',newfigure,'firstfigure',false);
create_subplot(@scatter, total_plots, current_plot, {amp_real,-zeta_real,'k'}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor,'figPosition',figPosition,'newfigure',newfigure,'firstfigure',false);
xlabel('Displacement')
ylabel('Aerodynamic damping ratio')
ylim([-0.2 0])
legend("calculated by displacement","calculated by force","real")
current_plot = current_plot + 1;


%% functions

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

function [A_c, B_c, G_c, J_c] = ssmod_c_mode(nmodes, omega2, Gamma, phi, S_a, S_v, S_d)
A_c = [zeros(nmodes), eye(nmodes); ...
    -omega2, -Gamma];
B_c = [zeros(nmodes, nmodes); ...
    eye(nmodes,nmodes)];
G_c = [S_d * phi - S_a * phi * omega2, ...
    S_v * phi - S_a * phi * Gamma];
J_c = [S_a * phi ];
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


function [ amp,fre ] = ee( data,dt )
%EE Summary of this function goes here
%   Detailed explanation goes here
[normalizeddata,amp]=splinenormalizeep(data);
frecarrier=gradient(normalizeddata,dt);
[normalizedfrecarrier,ampfrecarrier]=splinenormalizeep(frecarrier);
fre=ampfrecarrier/2/pi;
end