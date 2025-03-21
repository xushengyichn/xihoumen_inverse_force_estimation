%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: ShengyiXu xushengyichn@outlook.com
%Date: 2023-09-27 12:52:24
%LastEditors: xushengyichn xushengyichn@outlook.com
%LastEditTime: 2023-10-03 00:08:39
%FilePath: \Exercises-for-Techniques-for-estimation-in-dynamics-systemsd:\git\xihoumen_inverse_force_estimation\simulation signal\Main.m
%Description: 使用模拟荷载信号，测试反算气动力的准确性，以及参数选择的合理性
%
%Copyright (c) 2023 by ${git_name_email}, All Rights Reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clc;clear; close all
addpath(genpath("F:\git\ssm_tools"))
addpath(genpath("D:\Users\xushe\Documents\GitHub\ssm_tools"))
addpath(genpath("D:\Users\xushe\Documents\GitHub\Function_shengyi_package"))
addpath(genpath("F:\git\Function_shengyi_package"))
addpath(genpath("C:\Users\xushe\OneDrive\NAS云同步\Drive\0博士研究生\3大论文\研究内容\研究内容 3：气动力模型及参数识别；\反算气动力"))
addpath(genpath("D:\OneDrive\NAS云同步\Drive\0博士研究生\3大论文\研究内容\研究内容 3：气动力模型及参数识别；\反算气动力"))
addpath(genpath("F:\git\xihoumen_inverse_force_estimation"))
addpath(genpath("D:\git\xihoumen_inverse_force_estimation"))
addpath(genpath("/Users/xushengyi/Documents/GitHub/ssm_tools_sy"))
addpath(genpath("/Users/xushengyi/Documents/GitHub/Function_shengyi_package/"))
addpath(genpath("/Users/xushengyi/Documents/GitHub/xihoumen_inverse_force_estimation/"))
addpath(genpath("D:\XSY\20230930"))

subStreamNumberDefault = 2132;
run('InitScript.m');

%% 0 绘图参数
fig_bool = ON;
num_figs_in_row = 5; %每一行显示几个图
figPos = figPosSmall; %图的大小，参数基于InitScript.m中的设置
%设置图片间隔
gap_between_images = [0, 0];
figureIdx = 0;

parfor_loop_two_variables = OFF;
parfor_loop_three_variables = OFF;
parfor_loop_Q_R = OFF;

%% Generate force
modesel= [2,3,5,6,7,9,15,21,23,29,33,39,44,45];
nmodes = length(modesel);
% modesel= [23];
mode_resonant=23;
f1=0.321813289457002;
mode_resonant_seq = find(modesel==mode_resonant);
T = 10000; dt = 0.01; fs = 1 / dt; t = 0:dt:T;N = length(t);
lambda = 0.00002;  sigma_p =20;


Modal_Force_all = zeros(length(modesel),N);
[Modal_Force] = quasiperiod_force(lambda, sigma_p, f1, dt, t);
Modal_Force_all(mode_resonant_seq ,:)=Modal_Force;

% white_noise = randn(nmodes,N)*1;
% Modal_Force_all = Modal_Force_all+white_noise;
% rms(Modal_Force);


%% Finite element model



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



%% place the accelerometer

phi = mode_vec; %模态向量 每一列是一个模态
C = CC_eq; K = KK_eq; M = MM_eq;

Gamma = C; % 对应文献中表述的符号
omega2 = K;

% accelerometer location
% loc_acc= [578+1650/4*3;578+1650/2;578+1650/4];
loc_acc = [578+1650/4;578+1650/2;578+1650/4*3];
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

%% read real vibration data(just for compare)
acc_names=["Main span 1/4","Main span 1/2","Main span 3/4"];

result = viv2013(4);


% vibac2_name=result.vibac2_name;
vibac3_name=result.vibac3_name;
% vibac4_name=result.vibac4_name;
startDate = result.startDate;
endDate = result.endDate;


% vibac2_data=read_vib_data(vibac2_name);
vibac3_data=read_vib_data(vibac3_name);
% vibac4_data=read_vib_data(vibac4_name);


t_data = result.startDate:seconds(1/50):result.endDate-seconds(1/50);



%% calculate response

[xn_true, yn_true,~] = CalResponse(A_d, B_d, G_d, J_d, Modal_Force_all, 0, 0, N, x0, ns, n_sensors);
[xn, yn,~] = CalResponse(A_d, B_d, G_d, J_d, Modal_Force_all, Q, R, N, x0, ns, n_sensors);


% disp("RMS of accelration in the middle span is: "+num2str(rms(vibac3_data(:,2))))
% disp("MAX of accelration in the middle span is: "+num2str(max(abs(vibac3_data(:,2)))))
% disp("Reconstruct RMS of accelration in the middle span is: "+num2str(rms(yn_true(3,:))))
% disp("Reconstruct MAX of accelration in the middle span is: "+num2str(max(abs(yn_true(3,:)))))


%% perform kalman filter
Q_xd = Q;
np_m = nmodes;
B_c_m = [zeros(nmodes, np_m); ...
    eye(np_m, np_m)];
J_c_m = [S_a * phi];
B_d_m = A_c \ (A_d - eye(size(A_d))) * B_c_m;
J_d_m = J_c_m;


% lambdas_m = [lambda] * ones(1, np_m);
% % 
% sigma_ps_m = [sigma_p] * ones(1, np_m);


% lambdas_m = [3.59381e-07] * ones(1, np_m);
% 
% sigma_ps_m = [44.5556] * ones(1, np_m);

lambdas_m = [3.59381e-07] * ones(1, np_m);

sigma_ps_m = [50] * ones(1, np_m);

% omega_0= f1*2*pi;

% omega_variation = 0.9889;
omega_variation = 1;

for k1 = 1:nmodes
    omega_0(k1)=Freq(k1)*2*pi*omega_variation;
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



%% Virtual sensoring
loc_acc_v = [578+1650/2];
loc_vel_v = [578+1650/2];
loc_dis_v = [578+1650/2];

[S_a_v, S_v_v, S_d_v, n_sensors_v] = sensor_selection(loc_acc_v, loc_vel_v, loc_dis_v, node_loc, phi,nodeondeck,Mapping_data);

G_c_v = [S_d_v * phi - S_a_v * phi * omega2, S_v_v * phi - S_a_v * phi * Gamma];
J_c_v = [S_a_v * phi ];

h_hat = G_c_v * x_filt_original + J_c_v * p_filt_m;


%% 相位分析（假定都是单频信号）

% Acceleration
data1 = h_hat(1,:);
data2 = Modal_Force;

[f, magnitude] = fft_transform(fs,data1);
[~,loc]=max(magnitude);
f_max_data1=f(loc);
[f, magnitude] = fft_transform(fs,data2);
[~,loc]=max(magnitude);
f_max_data2=f(loc);
if f_max_data1-f_max_data1>1e-2
    warning("Frequency of the two signal are not coincide");
end

phase_lag_acc = fft_phase_lag(fs,data1,data2,f_max_data1)/pi*180;


% velocity
data1 = h_hat(3,:);
data2 = Modal_Force;

[f, magnitude] = fft_transform(fs,data1);
[~,loc]=max(magnitude);
f_max_data1=f(loc);
[f, magnitude] = fft_transform(fs,data2);
[~,loc]=max(magnitude);
f_max_data2=f(loc);
if f_max_data1-f_max_data1>1e-2
    warning("Frequency of the two signal are not coincide");
end

phase_lag_vel = fft_phase_lag(fs,data1,data2,f_max_data1)/pi*180;


% displacement
data1 = h_hat(5,:);
data2 = Modal_Force;

[f, magnitude] = fft_transform(fs,data1);
[~,loc]=max(magnitude);
f_max_data1=f(loc);
[f, magnitude] = fft_transform(fs,data2);
[~,loc]=max(magnitude);
f_max_data2=f(loc);
if f_max_data1-f_max_data1>1e-2
    warning("Frequency of the two signal are not coincide");
end

phase_lag_dis = fft_phase_lag(fs,data1,data2,f_max_data1)/pi*180;


%% recalculate the response with estimated force
[xn_re, yn_re,~] = CalResponse(A_d, B_d, G_d, J_d, p_filt_m, Q, R, N, x0, ns, n_sensors);


%% parfor loop two variables
if parfor_loop_two_variables == ON
n1 = 10;
n2 = 10;
lambdas_m_list = logspace(-8,-1, n1);
sigma_ps_m_list = linspace(1,50,n2);
[X, Y] = meshgrid(lambdas_m_list, sigma_ps_m_list);
combinations = [reshape(X, [], 1), reshape(Y, [], 1)];
numIterations = size(combinations,1);

if isempty(gcp('nocreate'))
    parpool();
end

b = ProgressBar(numIterations, ...
    'IsParallel', true, ...
    'WorkerDirectory', pwd(), ...
    'Title', 'Parallel 1' ...
    );
b.setup([], [], []);



parfor k1 = 1:numIterations
% for k1 = 1:numIterations
    parameters = combinations(k1,:);
    lambdas_m = parameters(1) * ones(1, np_m);
    sigma_ps_m = parameters(2)  * ones(1, np_m);
    % omega_0= f1*2*pi;
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
    [x_k_k, x_k_kmin, P_k_k, P_k_kmin,result] = KalmanFilterNoInput(A_a_m, G_a_m, Q_a_m, R_a_m, yn_a, x_ak, P_ak, 'debugstate', true,'showtext',false);
    logL(k1)=result.logL;
    logSk(k1) = result.logSk;
    logek(k1)=result.logek;

    % USE THIS FUNCTION AND NOT THE STEP() METHOD OF THE OBJECT!!!
    updateParallel([], pwd);
end

[~,MaxIdx]= max(logL);
lambda_Max = combinations(MaxIdx,1);
sigma_p_Max = combinations(MaxIdx,2);

b.release();

Z_1 = reshape(logL, n2, n1);
Z_2 = reshape(logSk, n2, n1);
Z_3 = reshape(logek, n2, n1);
end

%% parfor loop two variables Q R
if parfor_loop_Q_R == ON
n1 = 10;
n2 = 10;
Q_list = logspace(-10,-1, n1);
R_list = logspace(-10,-1, n2);
[X, Y] = meshgrid(Q_list, R_list);
combinations = [reshape(X, [], 1), reshape(Y, [], 1)];
numIterations = size(combinations,1);

if isempty(gcp('nocreate'))
    parpool();
end

b = ProgressBar(numIterations, ...
    'IsParallel', true, ...
    'WorkerDirectory', pwd(), ...
    'Title', 'Parallel 1' ...
    );
b.setup([], [], []);



parfor k1 = 1:numIterations
% for k1 = 1:numIterations
    parameters = combinations(k1,:);
    lambdas_m = lambda * ones(1, np_m);
    sigma_ps_m =sigma_p  * ones(1, np_m);
    % omega_0= f1*2*pi;
    [F_c_m, L_c_m, H_c_m, sigma_w_m12] = ssmod_quasiperiod_coninue(lambdas_m, sigma_ps_m, omega_0, np_m);
    
    Q_xd_temp = parameters(1) * eye(ns);
    R_temp = parameters(2) * eye(n_sensors);
    [~, ~, ~, ~, Fad_m, ~, Had_m, ~, Qad_m] = ssmod_lfm_aug(A_c, B_c_m, G_c, J_c_m, F_c_m, H_c_m, L_c_m,Q_xd_temp, sigma_w_m12, dt);
    A_a_m = Fad_m;
    G_a_m = Had_m;
    Q_a_m = Qad_m;
    R_a_m = R_temp;
    yn_a = yn;
    NN = N;
    xa_history = zeros(ns + np_m * (2), NN);
    pa_history = zeros(ns + np_m * (2), NN);
    
    x_ak = zeros(ns + np_m * (2), 1);
    P_ak = 10 ^ (1) * eye(ns + np_m * (2));
    
    % G_a=G_a_m; A_a=A_a_m; Q_a=Q_a_m;
    [x_k_k, x_k_kmin, P_k_k, P_k_kmin,result] = KalmanFilterNoInput(A_a_m, G_a_m, Q_a_m, R_a_m, yn_a, x_ak, P_ak, 'debugstate', true,'showtext',false);
    logL(k1)=result.logL;
    logSk(k1) = result.logSk;
    logek(k1)=result.logek;

    % USE THIS FUNCTION AND NOT THE STEP() METHOD OF THE OBJECT!!!
    updateParallel([], pwd);
end

[~,MaxIdx]= max(logL);
Q_Max = combinations(MaxIdx,1);
R_Max = combinations(MaxIdx,2);

b.release();

Z_1 = reshape(logL, n2, n1);
Z_2 = reshape(logSk, n2, n1);
Z_3 = reshape(logek, n2, n1);
end



%% parfor loop three variables
if parfor_loop_three_variables == ON
    n1 = 100;
    n2 = 100;
    n3 = 10;
    lambdas_m_list = logspace(-10,-1, n1);
    sigma_ps_m_list = linspace(5,100,n2);
    omega_0_list = linspace(0.9,1.1,n3);
    [X, Y,Z] = meshgrid(lambdas_m_list, sigma_ps_m_list,omega_0_list);
    combinations = [reshape(X, [], 1), reshape(Y, [], 1),reshape(Z, [], 1)];
    numIterations = size(combinations,1);
    
    if isempty(gcp('nocreate'))
        parpool();
    end
    
    b = ProgressBar(numIterations, ...
        'IsParallel', true, ...
        'WorkerDirectory', pwd(), ...
        'Title', 'Parallel 1' ...
        );
    b.setup([], [], []);
    
    
    
    parfor k1 = 1:numIterations
    % for k1 = 1:numIterations
        parameters = combinations(k1,:);
        lambdas_m = parameters(1) * ones(1, np_m);
        sigma_ps_m = parameters(2)  * ones(1, np_m);
        omega_0_m= parameters(3)  * omega_0;
        [F_c_m, L_c_m, H_c_m, sigma_w_m12] = ssmod_quasiperiod_coninue(lambdas_m, sigma_ps_m, omega_0_m, np_m);
        
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
        [x_k_k, x_k_kmin, P_k_k, P_k_kmin,result] = KalmanFilterNoInput(A_a_m, G_a_m, Q_a_m, R_a_m, yn_a, x_ak, P_ak, 'debugstate', true,'showtext',false);
        logL(k1)=result.logL;
        logSk(k1) = result.logSk;
        logek(k1)=result.logek;
    
        % USE THIS FUNCTION AND NOT THE STEP() METHOD OF THE OBJECT!!!
        updateParallel([], pwd);
    end
    
    [~,MaxIdx]= max(logL);
    lambda_Max = combinations(MaxIdx,1);
    sigma_p_Max = combinations(MaxIdx,2);
    
    b.release();
    
    % Z_1 = reshape(logL, n2, n1);
    % Z_2 = reshape(logSk, n2, n1);
    % Z_3 = reshape(logek, n2, n1);
    end

%% plot the figures
if fig_bool == ON
    for k1 = 1:nmodes
    [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
    hold on
    plot(t,p_filt_m(k1,:))
    plot(t,Modal_Force_all(k1,:))
    legend("filter","true")
    xlabel('t (s)')
    ylabel('Modal force')
    xlim([4000,4100])
    end

    [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
    plot(t,yn(3,:))
    hold on
    plot(t,yn_true(3,:))
    plot(t,yn_re(3,:))
    xlabel('t (s)')
    ylabel('Acceleration')
    legend("measurement","true","reconstruct")
    ylim([-0.5,0.5])
    title("RMS: "+num2str(rms(yn_true(3,:)))+"MAX: "+num2str(max(abs(yn_true(3,:)))))
    copygraphics(gcf)
        
    [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
    plot(t_data,vibac3_data(:,2))
    xlabel('t (s)')
    ylabel('Acceleration')
    ylim([-0.5,0.5])
    title("RMS: "+num2str(rms(vibac3_data(:,2)))+"MAX: "+num2str(max(abs(vibac3_data(:,2)))))
    copygraphics(gcf)

    [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
    plot(t,Modal_Force)
    xlabel('t (s)')
    title("Modal Force");
    copygraphics(gcf)

    
    [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
    plot(t,h_hat(1,:))
    hold on 
    plot(t,Modal_Force*0.1)
    xlabel('t (s)')
    title("Acceleration, phase lag:"+num2str(phase_lag_acc));
    xlim([4000,4100])

    [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
    plot(t,h_hat(3,:))
        hold on 
    plot(t,Modal_Force*0.1)
    title("Velocity, phase lag:"+num2str(phase_lag_vel));
    xlim([4000,4100])

    [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
    plot(t,h_hat(5,:))
        hold on 
    plot(t,Modal_Force*0.1)
    title("Displacement, phase lag:"+num2str(phase_lag_dis));
    xlim([4000,4100])


end

if parfor_loop_two_variables == ON



    [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
    % n = 256;
    % x = linspace(0, 1, n)';
    % baseColormap = parula(n);
    % customMap = [baseColormap(:,1), baseColormap(:,2), baseColormap(:,3)].*x.^2;  % 四次多项式映射
    % colormap(customMap);
    contourf(X, Y, Z_1,1000,'LineColor','none');   % 绘制等高线图

    % clim([3.0e+07 3.1512e+07]);
    hold on
    plot(lambda,sigma_p,'ro')
    plot(lambda_Max, sigma_p_Max, 'bo')  % 绘制最大值的位置
    set(gca, 'XScale', 'log');
    xlabel('Q');
    ylabel('R');
    colorbar;  % 添加颜色栏

    title('logL');
    

    [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
    test = sort(Z_1(:));
    test_min = min(test);
    test_max = max(test);
    test = test(2:5:end,:);
    clev = [test_min;test;test_max];
    hc = contourfcmap(log(X), Y, Z_1,clev,jet(length(clev)-1), 'lo', [.8 .8 .8], 'hi', [.2 .2 .2], 'cbarloc', 'eastoutside', 'method', 'calccontour',   'evencb', true);
    hold on
    plot(log(lambda),sigma_p,'ro')
    plot(log(lambda_Max), sigma_p_Max, 'bo')  % 绘制最大值的位置
    % logXTicks = -8:-1;  % 对数值
    % xticks(logXTicks);
    % newXTickLabels = "10^" + string(logXTicks);
    % xticklabels(newXTickLabels);
    title('logL');
    xlabel('log(lambdas)');
    ylabel('sigma_ps');


        [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
    contourf(X, Y, Z_2);  % 绘制等高线图
    set(gca, 'XScale', 'log');
    xlabel('lambdas');
    ylabel('sigma_ps');
    colorbar;  % 添加颜色栏
    title('logSk');


    [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
    contourf(X, Y, Z_3);  % 绘制等高线图
    set(gca, 'XScale', 'log');
    xlabel('lambdas');
    ylabel('sigma_ps');
    colorbar;  % 添加颜色栏
    title('logek');


    % save('result_two_variables.mat','X','Y','Z_1','Z_2','Z_3','logL','logSk','logek','lambda','sigma_p','lambda_Max','sigma_p_Max')
end

if parfor_loop_Q_R == ON



    [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
    % n = 256;
    % x = linspace(0, 1, n)';
    % baseColormap = parula(n);
    % customMap = [baseColormap(:,1), baseColormap(:,2), baseColormap(:,3)].*x.^2;  % 四次多项式映射
    % colormap(customMap);
    contourf(X, Y, Z_1,1000,'LineColor','none');   % 绘制等高线图

    % clim([3.0e+07 3.1512e+07]);
    hold on
    plot(10^(-8),10^(-5),'ro')
    plot(Q_Max, R_Max, 'bo')  % 绘制最大值的位置
    set(gca, 'XScale', 'log');
    set(gca, 'YScale', 'log');
    xlabel('Q');
    ylabel('R');
    colorbar;  % 添加颜色栏

    title('logL');
    

    [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
    test = sort(Z_1(:));
    test_min = min(test);
    test_max = max(test);
    test = test(2:5:end,:);
    clev = [test_min;test;test_max];
    hc = contourfcmap(log(X), log(Y), Z_1,clev,jet(length(clev)-1), 'lo', [.8 .8 .8], 'hi', [.2 .2 .2], 'cbarloc', 'eastoutside', 'method', 'calccontour',   'evencb', true);
    hold on
    plot(log(10^(-8)),log(10^(-5)),'ro')
    plot(log(Q_Max), log(R_Max), 'bo')  % 绘制最大值的位置
    % logXTicks = -8:-1;  % 对数值
    % xticks(logXTicks);
    % newXTickLabels = "10^" + string(logXTicks);
    % xticklabels(newXTickLabels);
    title('logL');
    xlabel('log(Q)');
    ylabel('log(R)');


    [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
    contourf(X, Y, Z_2);  % 绘制等高线图
    set(gca, 'XScale', 'log');
    xlabel('Q');
    ylabel('R');
    colorbar;  % 添加颜色栏
    title('logSk');


    [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
    contourf(X, Y, Z_3);  % 绘制等高线图
    set(gca, 'XScale', 'log');
    xlabel('Q');
    ylabel('R');
    colorbar;  % 添加颜色栏
    title('logek');


    % save('result_two_variables.mat','X','Y','Z_1','Z_2','Z_3','logL','logSk','logek','lambda','sigma_p','lambda_Max','sigma_p_Max')
end

if parfor_loop_three_variables == ON

    [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
    scatter3(X(:),Y(:),Z(:),10,logL(:))
    set(gca, 'XScale', 'log');
    colorbar;  % 添加颜色栏
    title('logL');
    hold on 
    scatter3(lambda,sigma_p,1,30,"red",'filled')
    [~,maxLogL_idx] = max(logL(:));
    x_list = X(:);
    y_list = Y(:);
    z_list = Z(:);
    x_max = x_list(maxLogL_idx);
    y_max = y_list(maxLogL_idx);
    z_max = z_list(maxLogL_idx);
    scatter3(x_max,y_max,z_max,30,"cyan",'filled')
    xlabel('lambda')
    ylabel('sigma_p')
    zlabel('omega variation')


    x = [-8*10^(9)  0 3e7 3.05e7 3.06e7 3.07e7 3.08e7 3.09e7 3.10e7 3.11e7 3.12e7 3.13e7 3.14e7 3.1489e7];
    Nx = length(x);
    c_lim = [min(x) max(x)];
    dx = min(diff(x));
    y = c_lim(1):dx:c_lim(2);
    for k=1:Nx-1, y(y>x(k) & y<=x(k+1)) = x(k+1); end % NEW
    cmap = colormap(parula(Nx));
    cmap2 = [...
    interp1(x(:),cmap(:,1),y(:)) ...
    interp1(x(:),cmap(:,2),y(:)) ...
    interp1(x(:),cmap(:,3),y(:)) ...
    ];

    [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
    x_list = X(:);
    y_list = Y(:);
    z_list = Z(:);
    logL_list = logL(:);

    scatter3(x_list,y_list,z_list,10,logL_list)
    colormap(cmap2)
    % clim(c_lim)
    set(gca, 'XScale', 'log');
    colorbar;  % 添加颜色栏
    title('logL');
    hold on 
    scatter3(lambda,sigma_p,1,30,"red",'filled')
    [~,maxLogL_idx] = max(logL(:));

    x_max = x_list(maxLogL_idx);
    y_max = y_list(maxLogL_idx);
    z_max = z_list(maxLogL_idx);
    scatter3(x_max,y_max,z_max,30,"cyan",'filled')
    xlabel('lambda')
    ylabel('sigma_p')
    zlabel('omega variation')

    z_list_unique = unique(z_list);
    for k2 = 1:length(z_list_unique)
        [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
        x_list = X(:);
        y_list = Y(:);
        z_list = Z(:);
        logL_list = logL(:);
        index = Z(:) == z_list_unique(k2);
    
        scatter3(x_list(index),y_list(index),z_list(index),10,logL_list(index))
        colormap(cmap2)
        % clim(c_lim)
        set(gca, 'XScale', 'log');
        colorbar;  % 添加颜色栏
        title('logL');
        hold on 
        scatter3(lambda,sigma_p,1,30,"red",'filled')
        [~,maxLogL_idx] = max(logL(:));
    
        x_max = x_list(maxLogL_idx);
        y_max = y_list(maxLogL_idx);
        z_max = z_list(maxLogL_idx);
        scatter3(x_max,y_max,z_max,30,"cyan",'filled')
        xlabel('lambda')
        ylabel('sigma_p')
        zlabel('omega variation')
        zlim([0.9,1.1])
    end
    % save('result_three_variables.mat','X','Y','Z','logL','logSk','logek','lambda','sigma_p','x_max','y_max','z_max')
end


% save('result.mat','X','Y','Z','logL','logSk','logek')
%% necessary functions
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

function [A_c, B_c, G_c, J_c] = ssmod_c_mode(nmodes, omega2, Gamma, phi, S_a, S_v, S_d)
A_c = [zeros(nmodes), eye(nmodes); ...
    -omega2, -Gamma];
B_c = [zeros(nmodes, nmodes); ...
    eye(nmodes,nmodes)];
G_c = [S_d * phi - S_a * phi * omega2, ...
    S_v * phi - S_a * phi * Gamma];
J_c = [S_a * phi ];
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

function [vibac_data]=read_vib_data(vibac_name)
vibac_data=[];
for i=1:length(vibac_name)
    vibac_data_temp=importdata(vibac_name(i));
    vibac_data_temp(:,2:4)=vibac_data_temp(:,2:4)/1000*9.8;%第一列是是时间，不要修改
    vibac_data=[vibac_data;vibac_data_temp];
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