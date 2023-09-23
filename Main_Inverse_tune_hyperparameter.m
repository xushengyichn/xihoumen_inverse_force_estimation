%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com Date:
%LastEditors: ShengyiXu xushengyichn@outlook.com
%LastEditTime: 2023-09-22 21:44:02
%Description: 计算简支梁施加荷载后的动力响应，并反算出荷载（分别按照集中力和模态力反算）
%
%Copyright (c) 2023 by ${git_name_email}, All Rights Reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc; clear; close all;
addpath(genpath("F:\git\ssm_tools"))
addpath(genpath("D:\Users\xushe\Documents\GitHub\ssm_tools"))
addpath(genpath("D:\Users\xushe\Documents\GitHub\Function_shengyi_package"))
addpath(genpath("F:\git\Function_shengyi_package"))
addpath(genpath("C:\Users\xushe\OneDrive\NAS云同步\Drive\0博士研究生\3大论文\研究内容\研究内容 3：气动力模型及参数识别；\反算气动力"))

subStreamNumberDefault = 2132;
run('InitScript.m');

%% 0 绘图参数
fig_bool = ON;
num_figs_in_row = 10; %每一行显示几个图
figPos = figPosSmall; %图的大小，参数基于InitScript.m中的设置
%设置图片间隔
gap_between_images = [0, 0];
figureIdx = 0;

%% 1 set hyperparameters

modesel= [2,3,5,6,7,9,15,21,23,29,33,39,44,45];
nmodes = length(modesel);
np_m = nmodes;
n1 = 100;
n2 = 100;

lambdas_m_list = logspace(-2, 2, n1);
sigma_ps_m_list = linspace(10,1000,n2);
[X, Y] = meshgrid(lambdas_m_list, sigma_ps_m_list);
combinations = [reshape(X, [], 1), reshape(Y, [], 1)];

numIterations = size(combinations,1);

if isempty(gcp('nocreate'))
    parpool();
end

b = ProgressBar(numIterations, ...
    'IsParallel', true, ...
    'WorkerDirectory', pwd(), ...
    'Title', 'Parallel 2' ...
    );
b.setup([], [], []);
parfor k1 = 1:numIterations
    lambdas_m = [combinations(k1,1)] * ones(1, np_m);
    sigma_ps_m = [combinations(k1,2)] * ones(1, np_m);

    result= Inverse_fun(lambdas_m,sigma_ps_m,modesel);
    logL(k1)=result.logL;
    logSk(k1) = result.logSk;
    logek(k1)=result.logek;

    % USE THIS FUNCTION AND NOT THE STEP() METHOD OF THE OBJECT!!!
    updateParallel([], pwd);
end

b.release();

Z_1 = reshape(logL, n2, n1);
% Z_2 = reshape(logSk, 50, 50);
% Z_3 = reshape(logek, 50, 50);

save tune_hyperparameter_result

%% plot

if fig_bool == ON

    [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);

    
    contourf(X, Y, Z_1);  % 绘制等高线图
    set(gca, 'XScale', 'log');
    xlabel('lambdas');
    ylabel('sigma_ps');
    colorbar;  % 添加颜色栏
    title('Contour Plot');

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
