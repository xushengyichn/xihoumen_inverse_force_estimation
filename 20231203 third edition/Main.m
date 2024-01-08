clc; clear; close all;
run('CommonCommand.m');
tic
%% 输入参数
% input.num_figs_in_row = 12; %每一行显示几个图
% input.figPos = figPos; %图的大小，参数基于InitScript.m中的设置
%设置图片间隔
% input.ON =ON;
% input.OFF =OFF;
% input.gap_between_images = [0, 0];
% input.figureIdx = 0;
% n = 13;
% n = 4;
n = 4;
[result] = viv2013(n, OFF);
startDate_global = result.startDate;
endDate_global = result.endDate;
input_data.start_time = startDate_global;
input_data.end_time = endDate_global;

% input.lambda_VIV = 10 ^ (-8.725260766057426);
% input.sigma_p_VIV = 4.594524428349437e+04;
% input.omega_0_variation_VIV = 1.001880394311541;
% input.Q_value = 10 ^ (-5.056871357602549);
% input.sigma_noise = 10 ^ (-1.301006929986168);
% input.sigma_buff = 10 ^ (0.070362659875298);

input_data.lambda_VIV = 10 ^ (-3.059366447657016);
input_data.sigma_p_VIV = 5.682743296566430e+03;
input_data.Q_value = 10 ^ (-5.004258502467064);
input_data.sigma_buff = 10 ^ (1.237973827837277);

input_data.omega_0_variation_VIV = 1;
input_data.sigma_noise = 10 ^ (-1.301006929986168);

%% logL优化参数
modesel = [2, 3, 5, 6, 7, 9, 15, 21, 23, 29, 33, 39, 44, 45];
% modesel = [2,3,23];

VIV_mode_seq = find(modesel == 23);
% modesel = 23;
input_data.modesel = modesel;
input_data.VIV_mode_seq = VIV_mode_seq;
nVIV = length(VIV_mode_seq);
input_data.nVIV = nVIV;

showtext = true;
showplot = false;


%% plot logL with changing parameters
if 1
    lambda_VIVs = logspace(-4, -0, 10);
    sigma_p_VIVs = linspace(1e1, 1e4, 10);

    [lambda_VIVs, sigma_p_VIVs] = meshgrid(lambda_VIVs, sigma_p_VIVs);
    variables = [lambda_VIVs(:), sigma_p_VIVs(:)];

    logLs = zeros(size(variables, 1), 1);
    logSks = zeros(size(variables, 1), 1);
    logeks = zeros(size(variables, 1), 1);

    parfor k1 = 1:size(variables, 1)
        input_data_temp = input_data;
        input_data_temp.lambda_VIV = variables(k1, 1);
        input_data_temp.sigma_p_VIV = variables(k1, 2);
        input_data_temp.omega_0_variation_VIV = 1;
        input_data_temp.Q_value = 10 ^ (-5.004258502467064);
        input_data_temp.sigma_noise = 10 ^ (-1.301006929986168);
        input_data_temp.sigma_buff = 10 ^ (1.237973827837277);
        [result_Main] = KalmanMain(input_data_temp, 'showtext', false, 'showplot', false, 'filterstyle', 'nofilter');
        logLs(k1) = result_Main.logL;
        logSks(k1) = result_Main.logSk;
        logeks(k1) = result_Main.logek;
    end

    logLs = reshape(logLs, size(lambda_VIVs));
    logSks = reshape(logSks, size(lambda_VIVs));
    logeks = reshape(logeks, size(lambda_VIVs));

    figure
    surf(log(lambda_VIVs), sigma_p_VIVs, logLs)
    xlabel('lambda_VIV')
    ylabel('sigma_p_VIV')
    zlabel('logL')

    figure
    surf(log(lambda_VIVs), sigma_p_VIVs, logSks)
    xlabel('lambda_VIV')
    ylabel('sigma_p_VIV')
    zlabel('logSk')

    figure
    surf(log(lambda_VIVs), sigma_p_VIVs, logeks)
    xlabel('lambda_VIV')
    ylabel('sigma_p_VIV')
    zlabel('logek')







end


%% Apply Kalman Filter
% [result_Main] = KalmanMain(input, 'showtext', true, 'showplot', false, 'filterstyle', 'fft', 'f_keep', 0.33 * [0.9, 1.1]);
[result_Main] = KalmanMain(input_data, 'showtext', showtext, 'showplot', showplot, 'filterstyle', 'nofilter');
logL = result_Main.logL;
logSk = result_Main.logSk;
logek = result_Main.logek;
% display logL logSk logek
dispstr = sprintf("logL = %f, logSk = %f, logek = %f", logL, logSk, logek);

if showtext
    disp(dispstr)
end

mode_deck = result_Main.mode_deck;

if 0 %优化所有参数
    %% 导入直接积分获得的涡激力
    ft_directint = importdata("DirectIntegration.mat");
    %% optimization logL to get the maximum with changing lambda sigma_p omega_0_variation Q_value R_value
    % 在调用 ga 函数之前，您可以这样设置 external_params：
    external_params.modesel = [2, 3, 5, 6, 7, 9, 15, 21, 23, 29, 33, 39, 44, 45];
    external_params.acc_dir = input_data.acc_dir;
    external_params.VIV_mode_seq = VIV_mode_seq;
    external_params.nVIV = nVIV;
    external_params.ft_directint = ft_directint;
    % external_params.modesel = [23];
    % 定义参数的范围
    lb = [-12, 1e2, 0.98, -10, -8, 0]; % 这里的值是假设的，请根据您的情况进行修改
    ub = [-4, 1e5, 1.02, -5, -1, 2]; % 这里的值也是假设的
    
    % 定义整数和连续变量
    IntCon = []; % 如果没有整数变量，否则提供整数变量的索引
    
    options = optimoptions('ga', 'MaxGenerations', 100, 'Display', 'iter', 'UseParallel', true);
    [x, fval] = ga(@(params) fitnessFunction(params, external_params), 6, [], [], [], [], lb, ub, [], IntCon, options);
    % 保存结果
    save('optimization_results.mat', 'x', 'fval');
    
end

if 0 % 选择参数进行优化
    %% 导入直接积分获得的涡激力
    ft_directint = importdata("DirectIntegration.mat");
    %% optimization logL to get the maximum with changing lambda sigma_p omega_0_variation Q_value R_value
    % 在调用 ga 函数之前，您可以这样设置 external_params：
    external_params.modesel = [2, 3, 5, 6, 7, 9, 15, 21, 23, 29, 33, 39, 44, 45];
    external_params.acc_dir = input_data.acc_dir;
    external_params.VIV_mode_seq = VIV_mode_seq;
    external_params.nVIV = nVIV;
    external_params.ft_directint = ft_directint;
    external_params.omega_0_variation_VIV = input_data.omega_0_variation_VIV;
    external_params.sigma_noise = input_data.sigma_noise;
    % external_params.modesel = [23];
    % 定义参数的范围
    lb = [-4, 1e1, -10, 0]; % 这里的值是假设的，请根据您的情况进行修改
    ub = [0, 1e4,  -5,  2]; % 这里的值也是假设的
    
    % 定义整数和连续变量
    IntCon = []; % 如果没有整数变量，否则提供整数变量的索引
    
    options = optimoptions('ga', 'MaxGenerations', 100, 'Display', 'iter', 'UseParallel', true);
    [x, fval] = ga(@(params) fitnessFunction_sel(params, external_params), 4, [], [], [], [], lb, ub, [], IntCon, options);
    % 保存结果
    save('optimization_results.mat', 'x', 'fval');
    
    input_data.lambda_VIV = 10 ^ (x(1));
    input_data.sigma_p_VIV = x(2);
    input_data.Q_value = 10 ^ (x(3));
    input_data.sigma_buff = 10 ^ (x(4));
    
    
end

%% Recalcualte the modal displacement
nmodes = result_Main.nmodes;
Result = ImportMK(nmodes, 'KMatrix.matrix', 'MMatrix.matrix', 'nodeondeck.txt', 'KMatrix.mapping', 'nodegap.txt', 'modesel', modesel);
mode_deck = Result.mode_deck; mode_deck_re = Result.mode_deck_re; node_loc = Result.node_loc; nodeondeck = Result.nodeondeck;
KMmapping = Result.Mapping;
nodegap = Result.nodegap;
mode_vec = Result.mode_vec;
loc_acc = [1403];
% loc_acc = [990.5; 1403; 1815.5];
node_shape = FindModeShapewithLocation(loc_acc, node_loc, nodeondeck, KMmapping, nodegap, mode_vec);

h_hat = result_Main.h_hat;
MM = result_Main.MM;
CC = result_Main.CC;
KK = result_Main.KK;
t_temp = result_Main.t;
t_temp = (datenum(t_temp) - datenum('1970-01-01 00:00:00')) * 86400; % Convert to seconds since epocht = result_Main.t;
t_temp = t_temp - t_temp(1); % Subtract the first element of t from all elements of t
p_filt_m = result_Main.p_filt_m;
u0 = zeros(nVIV, 1);
udot0 = zeros(nVIV, 1);
[u udot u2dot] = NewmarkInt(t_temp, MM(VIV_mode_seq, VIV_mode_seq), CC(VIV_mode_seq, VIV_mode_seq), KK(VIV_mode_seq, VIV_mode_seq), p_filt_m, 1/2, 1/4, u0, udot0);

input_data.u = u;
input_data.udot = udot;
input_data.u2dot = u2dot;

%% replace the force
yn = result_Main.yn;
yn_reconstruct = result_Main.yn_reconstruct;
% disp_dir = acc2dsip(yn(1, :), 50);
% F_direct = MM.*yn(1,:)/node_shape+CC.*disp_dir.vel/node_shape + KK.*disp_dir.disp/node_shape;
%
% input.p_filt_m=F_direct;

%% Caldamping ratio
fields = fieldnames(result_Main);

for i = 1:numel(fields)
    input_data.(fields{i}) = result_Main.(fields{i});
end

% input = result_Main;
input_data.ncycle = 10;

[result_Damping] = Cal_aero_damping_ratio(input_data, 'showplot', false, 'filterstyle', 'fft');
% [result_Damping] = Cal_aero_damping_ratio(input, 'showplot', false, 'filterstyle', 'nofilter');

amp_cell = result_Damping.amp_cell;
zeta_all_cell = result_Damping.zeta_all_cell;
top_freqs = result_Damping.top_freqs;
t_cycle_mean_cell = result_Damping.t_cycle_mean_cell;

%% read wind data

% [result] = viv2013(n, OFF);
start_time = input_data.start_time;
end_time = input_data.end_time;
wind_dir = input_data.wind_dir;
[result_wind] = read_wind_data(start_time, end_time, wind_dir);

t = result_Main.t;
F_filter = p_filt_m;

amp_temp = amp_cell{1}{1};
t_cycle_mean_temp = t_cycle_mean_cell{1}{1};
duration = minutes(10);


start_time.Format = 'dd_MMM_yyyy_HH_mm_ss';
formatted_start_time = string(start_time);

end_time.Format = 'dd_MMM_yyyy_HH_mm_ss';
formatted_end_time = string(end_time);


formatted_start_time = strrep(strrep(strrep(formatted_start_time, ':', '-'), ' ', '_'), '/', '-');
formatted_end_time = strrep(strrep(strrep(formatted_end_time, ':', '-'), ' ', '_'), '/', '-');

filename = "windspeed_result_"+formatted_start_time+"_"+formatted_end_time+".mat";

if exist(filename, 'file') == 2
    load(filename);
else
    windspeed_result = read_wind_data_onetime(t_cycle_mean_temp, duration, input_data.wind_dir_all, [1 1 0 0 0 0]);
    save(filename,'windspeed_result',"start_time","end_time");
end



%% plot
if fig_bool
    % 定义总子图数量
    total_plots = 20; % 或任何你需要的子图数量
    current_plot = 1;
    num_figs_in_row = [];
    figWidthFactor = 1.5;
    %     figPosition = [1080*2.5,100];
    figPosition = [100, 100];
    newfigure = true;
    holdon = false;
    
    %% modal force
    for k1 = 1:nVIV
        
        if k1 == 1
            create_subplot(@plot, total_plots, current_plot, {t, F_filter(k1, :)}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', newfigure, 'firstfigure', true, 'holdon', holdon);
            hold on
            NN = length(t);
            tnew=t';
            t_fill = [tnew(1:NN) fliplr(tnew(1:NN))];
            Pp_filt_m_cov = result_Main.Pp_filt_m;
            Pp_filt_m_cov_upper_bound =  F_filter(k1, :)+sqrt(Pp_filt_m_cov );
            Pp_filt_m_cov_lower_bound =  F_filter(k1, :)-sqrt(Pp_filt_m_cov );
            h_hat_fill = [Pp_filt_m_cov_upper_bound fliplr(Pp_filt_m_cov_lower_bound)];
            [t_fill, h_hat_fill] = reduceDataPoints(t_fill, h_hat_fill, 10);
            hold on
            create_subplot(@fill, total_plots, current_plot, { t_fill, h_hat_fill,'red','FaceAlpha',0.5,'EdgeColor','none'}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'holdon', holdon);


        else
            create_subplot(@plot, total_plots, current_plot, {t, F_filter(k1, :)}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', newfigure, 'holdon', holdon);
        end
        
        legend("Force from Kalman Filter")
        title("modal force for " + "mode"+modesel(VIV_mode_seq(k1)));
        current_plot = current_plot + 1;
        ylim([-100, 100])
    end
    [f, magnitude] = fft_transform(50,F_filter(1, :));
    create_subplot(@plot, total_plots, current_plot, {f, magnitude}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', newfigure, 'holdon', holdon);
    title("Force in the frequency domain")
    xlabel("Frequency")
    ylabel("Amplitude")
    current_plot = current_plot + 1;
    save temp.mat t F_filter
    %% installation of the sensors and the rms of the sensors
    % 三个传感器位置的振型大小
    loc_acc = result_Main.loc_acc;
    loc_acc_shape = FindModeShapewithLocation(loc_acc, node_loc, nodeondeck, KMmapping, nodegap, mode_vec);
    max_loc_acc_shape = max(abs(loc_acc_shape(:, VIV_mode_seq)));
    % 三个传感器位置的振动rms大小，并基于振型大小进行缩放
    yn_rms = rms(yn, 2);
    max_acc = max(yn_rms);
    scale_factor = max_acc / max_loc_acc_shape;
    yn_rms_scaled = yn_rms / scale_factor;
    create_subplot(@plot, total_plots, current_plot, {node_loc, mode_deck(:, VIV_mode_seq)}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'holdon', holdon);
    hold on
    create_subplot(@scatter, total_plots, current_plot, {loc_acc, loc_acc_shape(:, VIV_mode_seq)}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'holdon', holdon);
    create_subplot(@scatter, total_plots, current_plot, {loc_acc, yn_rms_scaled([1, 3, 5], :)}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'holdon', holdon);
    legend("modal shape", "modal shape at the sensors location", "scaled rms value of the sensors", 'Location', 'southeast')
    title("Mode shape and the locaiton of the sensors")
    current_plot = current_plot + 1;
    xlabel('Time (s)')
    ylabel('Dimensionless displacement')
    
    %% reconstructed data comparison (reconstructed vs kalman filter vitural sensoring)
    create_subplot(@plot, total_plots, current_plot, { t, h_hat(15, :),t, node_shape(VIV_mode_seq) * u}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'holdon', holdon);
    legend( 'filtered from kalman filter',"Kalman filter reconstruct")
    % title("reconstructed displacement vs kalman filter vitural sensoring")
    title("Dis (1/2span)")
    current_plot = current_plot + 1;
    xlabel('Time (s)')
    ylabel('Displacement (m)')
    
    create_subplot(@plot, total_plots, current_plot, { t, h_hat(9, :),t, node_shape(VIV_mode_seq) * udot}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'holdon', holdon);
    legend( 'filtered from kalman filter',"Kalman filter reconstruct")
    % title("reconstructed velocity vs kalman filter vitural sensoring")
    title("Vel (1/2span)")
    current_plot = current_plot + 1;
    xlabel('Time (s)')
    ylabel('Velocity (m/s)')
    
    create_subplot(@plot, total_plots, current_plot, { t, h_hat(3, :),t, node_shape(VIV_mode_seq) * u2dot}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'holdon', holdon);
    legend( 'filtered from kalman filter',"Kalman filter reconstruct")
    % title("reconstructed acceleration vs kalman filter vitural sensoring")
    title("Acc (1/2span)")
    current_plot = current_plot + 1;
    xlabel('Time (s)')
    ylabel('Acceleration (m/s^2)')
    

    %% kalman filter vitural sensoring and the covariance
    create_subplot(@plot, total_plots, current_plot, { t, h_hat(15, :)}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'holdon', holdon);
    NN = length(t);
    tnew=t';
    t_fill = [tnew(1:NN) fliplr(tnew(1:NN))];
    h_hat_cov = result_Main.h_hat_covariance(15,15);
    h_hat_upper_bound = h_hat(15, :)+sqrt(h_hat_cov);
    h_hat_lower_bound = h_hat(15, :)-sqrt(h_hat_cov);
    h_hat_fill = [h_hat_upper_bound fliplr(h_hat_lower_bound)];
    hold on
    [t_fill, h_hat_fill] = reduceDataPoints(t_fill, h_hat_fill, 10);
    create_subplot(@fill, total_plots, current_plot, { t_fill, h_hat_fill,'red','FaceAlpha',0.5,'EdgeColor','none'}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'holdon', holdon);
   
    legend( 'filtered from kalman filter')
    % title("reconstructed displacement vs kalman filter vitural sensoring")
    title("Dis (1/2span)")
    current_plot = current_plot + 1;
    xlabel('Time (s)')
    ylabel('Displacement (m)')
    
    create_subplot(@plot, total_plots, current_plot, { t, h_hat(9, :)}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'holdon', holdon);
    tnew=t';
    t_fill = [tnew(1:NN) fliplr(tnew(1:NN))];
    h_hat_cov = result_Main.h_hat_covariance(9,9);
    h_hat_upper_bound = h_hat(9, :)+sqrt(h_hat_cov);
    h_hat_lower_bound = h_hat(9, :)-sqrt(h_hat_cov);
    h_hat_fill = [h_hat_upper_bound fliplr(h_hat_lower_bound)];
    [t_fill, h_hat_fill] = reduceDataPoints(t_fill, h_hat_fill, 10);
    hold on
    create_subplot(@fill, total_plots, current_plot, { t_fill, h_hat_fill,'red','FaceAlpha',0.5,'EdgeColor','none'}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'holdon', holdon);

    legend( 'filtered from kalman filter')
    % title("reconstructed velocity vs kalman filter vitural sensoring")
    title("Vel (1/2span)")
    current_plot = current_plot + 1;
    xlabel('Time (s)')
    ylabel('Velocity (m/s)')
    
    create_subplot(@plot, total_plots, current_plot, { t, h_hat(3, :)}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'holdon', holdon);
    tnew=t';
    t_fill = [tnew(1:NN) fliplr(tnew(1:NN))];
    h_hat_cov = result_Main.h_hat_covariance(3,3);
    h_hat_upper_bound = h_hat(3, :)+sqrt(h_hat_cov);
    h_hat_lower_bound = h_hat(3, :)-sqrt(h_hat_cov);
    h_hat_fill = [h_hat_upper_bound fliplr(h_hat_lower_bound)];
    [t_fill, h_hat_fill] = reduceDataPoints(t_fill, h_hat_fill, 10);
    hold on
    create_subplot(@fill, total_plots, current_plot, { t_fill, h_hat_fill,'red','FaceAlpha',0.5,'EdgeColor','none'}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'holdon', holdon);

    legend( 'filtered from kalman filter')
    % title("reconstructed acceleration vs kalman filter vitural sensoring")
    title("Acc (1/2span)")
    current_plot = current_plot + 1;
    xlabel('Time (s)')
    ylabel('Acceleration (m/s^2)')
    %% wind speed during the period
    % beta_deg_mean_UA5 = result_wind.resultsTable_UA5.beta_deg_mean;
    % beta_deg_mean_UA6 = result_wind.resultsTable_UA6.beta_deg_mean;
    % % judge if the wind direction of both UA5 and UA6 is from 45°-225°
    % for k1 = 1:length(beta_deg_mean_UA5)
    %
    %     if and(beta_deg_mean_UA5(k1) > 45, beta_deg_mean_UA5(k1) < 225) && and(beta_deg_mean_UA6(k1) > 45, beta_deg_mean_UA6(k1) < 225)
    %         U_sel(k1) = result_wind.resultsTable_UA5.U(k1);
    %         % judge if the wind direction of both UA5 and UA6 is from 225°-360° or 0°-45°
    %     elseif or(and(beta_deg_mean_UA5(k1) > 225, beta_deg_mean_UA5(k1) < 360), and(beta_deg_mean_UA5(k1) > 0, beta_deg_mean_UA5(k1) < 45)) && or(and(beta_deg_mean_UA6(k1) > 225, beta_deg_mean_UA6(k1) < 360), and(beta_deg_mean_UA6(k1) > 0, beta_deg_mean_UA6(k1) < 45))
    %         U_sel(k1) = result_wind.resultsTable_UA6.U(k1);
    %     else
    %
    %         if abs(result_wind.resultsTable_UA5.alpha_deg_mean(k1)) < abs(result_wind.resultsTable_UA6.alpha_deg_mean)
    %             U_sel(k1) = result_wind.resultsTable_UA5.U(k1);
    %         else
    %             U_sel(k1) = result_wind.resultsTable_UA6.U(k1);
    %         end
    %
    %         disp("该时间点两个风速仪风向不一致，取风攻角较小的值")
    %     end
    %
    % end
    %
    % create_subplot(@scatter, total_plots, current_plot, {result_wind.resultsTable_UA6.Time_Start, U_sel}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure,'holdon',holdon);
    % xlabel('Time (s)')
    % ylabel('Wind speed (m/s)')
    % title("Wind speed vs. Time")
    % current_plot = current_plot + 1;
    
    %% measured, filtered and recalcuated acceleration
    create_subplot(@plot, total_plots, current_plot, {t, yn(1, :), t, h_hat(1, :), t, yn_reconstruct(1, :)}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'holdon', holdon);
    xlabel('Time (s)')
    ylabel('Acceleration (m/s^2)')
    title("Acceleration vs. Time 1/4 span")
    legend("measure", "kalman filter", "reconstruct using force from kalman filter")
    current_plot = current_plot + 1;
    
    create_subplot(@plot, total_plots, current_plot, {t, yn(3, :), t, h_hat(3, :), t, yn_reconstruct(3, :)}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'holdon', holdon);
    xlabel('Time (s)')
    ylabel('Acceleration (m/s^2)')
    title("Acceleration vs. Time 1/2 span")
    legend("measure", "kalman filter", "reconstruct using force from kalman filter")
    current_plot = current_plot + 1;
    
    create_subplot(@plot, total_plots, current_plot, {t, yn(5, :), t, h_hat(5, :), t, yn_reconstruct(5, :)}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'holdon', holdon);
    xlabel('Time (s)')
    ylabel('Acceleration (m/s^2)')
    title("Acceleration vs. Time 3/4 span")
    legend("measure", "kalman filter", "reconstruct using force from kalman filter")
    current_plot = current_plot + 1;
    
    %% Damping ratio calculation
    amp_temp = amp_cell{1}{1};
    % t_cycle_mean_temp = t_cycle_mean_cell{1}{1};
    m_cycle = input_data.ncycle; %cycles to be averaged
    zetam = zeros(1, length(t_cycle_mean_temp)); % Pre-allocate zetam with zeros
    
    for k1 = 1:length(t_cycle_mean_temp)
        
        if k1 <= m_cycle / 2 % Beginning boundary
            start_idx = 1;
            end_idx = start_idx + m_cycle;
        elseif k1 > length(t_cycle_mean_temp) - m_cycle / 2 % Ending boundary
            end_idx = length(t_cycle_mean_temp);
            start_idx = end_idx - m_cycle;
        else % Middle
            start_idx = k1 - floor(m_cycle / 2);
            end_idx = k1 + floor(m_cycle / 2);
        end
        
        deltam = log(amp_temp(start_idx) / amp_temp(end_idx));
        
        zetam(k1) = sqrt(deltam ^ 2 / (4 * m_cycle ^ 2 * pi ^ 2 + deltam ^ 2));
        
        if deltam > 0
            zetam(k1) = abs(zetam(k1));
        else
            zetam(k1) = -abs(zetam(k1));
        end
        
    end
    
    zeta_structure = result_Main.zeta;
    zeta1 = zeta_all_cell{1}{1};
    zeta2 = zetam - zeta_structure(VIV_mode_seq);
    
    target = norm(zeta1 - zeta2);
    disp(target)
    
    for k1 = 1:nVIV
        
        for k2 = 1:length(top_freqs{k1})
            % 假设 t_cycle_mean_cell{k1}{k2} 是一个包含 datetime 对象的数组
            datetimeArray = t_cycle_mean_cell{k1}{k2};
            
            % 提取第一个 datetime 对象作为参考点
            referenceDatetime = datetimeArray(1);
            
            % 计算每个 datetime 对象相对于参考点的秒数
            secondsFromReference = seconds(datetimeArray - referenceDatetime);
            
            % 现在，secondsFromReference 包含相对于第一个时间戳的秒数
            
            create_subplot(@scatter, total_plots, current_plot, {amp_temp * max(mode_deck(:, VIV_mode_seq(k1))), zetam - zeta_structure(VIV_mode_seq)}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', false, 'holdon', holdon);
            hold on
            create_subplot(@scatter, total_plots, current_plot, {amp_cell{k1}{k2} * max(mode_deck(:, VIV_mode_seq(k1))), zeta_all_cell{k1}{k2}, [], secondsFromReference, 'filled'}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', false, 'holdon', holdon);
            
            current_plot = current_plot + 1;
            
            % 设置 colormap
            colormap('jet')
            colorbar
            
            hold on
            plot([0, 0.15], [-0.003, -0.003])
            % scatter(ex,epsx,'green')
            str = "Mode : %d, Frequency : %.2f Hz";
            title(sprintf(str, modesel(k1), top_freqs{k1}(k2)));
            xlim([0.05, 0.12])
            ylim([-0.5, 0.5] / 100)
            xlabel("Amplitude(m)")
            ylabel("Damping ratio")
        end
        
    end
    
    %% wind speed with amplitude and damping ratio
    
    UA5 = windspeed_result.UA1;
    UA6 = windspeed_result.UA2;
    
    beta_deg_mean_UA5 = UA5.beta_deg_mean;
    beta_deg_mean_UA6 = UA6.beta_deg_mean;
    
    % judge if the wind direction of both UA5 and UA6 is from 45°-225°
    for k1 = 1:length(beta_deg_mean_UA5)
        
        if and(beta_deg_mean_UA5(k1) > 45, beta_deg_mean_UA5(k1) < 225) && and(beta_deg_mean_UA6(k1) > 45, beta_deg_mean_UA6(k1) < 225)
            U_sel(k1) = UA5.U(k1);
            AoA_sel(k1) = UA5.alpha_deg_mean(k1);
            
            % judge if the wind direction of both UA5 and UA6 is from 225°-360° or 0°-45°
        elseif or(and(beta_deg_mean_UA5(k1) > 225, beta_deg_mean_UA5(k1) < 360), and(beta_deg_mean_UA5(k1) > 0, beta_deg_mean_UA5(k1) < 45)) && or(and(beta_deg_mean_UA6(k1) > 225, beta_deg_mean_UA6(k1) < 360), and(beta_deg_mean_UA6(k1) > 0, beta_deg_mean_UA6(k1) < 45))
            U_sel(k1) = UA6.U(k1);
            AoA_sel(k1) = UA6.alpha_deg_mean(k1);
        else
            
            if abs(UA5.alpha_deg_mean(k1)) < abs(UA6.alpha_deg_mean)
                U_sel(k1) = UA5.U(k1);
                AoA_sel(K1) = UA5.alpha_deg_mean(k1);
            else
                U_sel(k1) = UA6.U(k1);
                AoA_sel(K1) = UA6.alpha_deg_mean(k1);
            end
            
            disp("该时间点两个风速仪风向不一致，取风攻角较小的值")
        end
        
    end
    
    create_subplot(@scatter, total_plots, current_plot, {t_cycle_mean_temp, U_sel}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'holdon', holdon);
    xlabel('Time (s)')
    ylabel('Wind speed (m/s)')
    title("Wind speed vs. Time")
    current_plot = current_plot + 1;
    
    create_subplot(@scatter, total_plots, current_plot, {t_cycle_mean_temp, AoA_sel}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'holdon', holdon);
    xlabel('Time (s)')
    ylabel('AoA (deg)')
    title("AoA vs. Time")
    current_plot = current_plot + 1;
    
    % figure
    amp_filter = amp_cell{1}{1} * max(mode_deck(:, VIV_mode_seq(1)));
    zeta_filter = zeta_all_cell{1}{1};
    % C=rescale(AoA_sel,10,100);
    S = rescale(secondsFromReference, 10, 100);
    C = AoA_sel;
    % scatter3(amp_filter,U_sel,zeta_filter,S,C)
    
    create_subplot(@scatter3, total_plots, current_plot, {amp_filter, U_sel, zeta_filter, S, C}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'holdon', holdon);
    xlabel('Amp. (m)')
    ylabel('Wind speed (m/s)')
    title("Damping ratio")
    current_plot = current_plot + 1;
    
    % 设置 colormap
    
    colorbar
    colormap('jet'); % 选择一个颜色映射，比如 'parula'
    % clim([-2 0]); % 设置颜色映射的数据范围为 -3 到 +3
    
    zlim([-0.01, 0])
    xlim([0.01, 0.09])
    ylim([5, 12])
    % 定义平面的四个角的 x, y, 和 z 坐标
    x = [min(amp_filter), max(amp_filter), max(amp_filter), min(amp_filter)];
    y = [min(U_sel), min(U_sel), max(U_sel), max(U_sel)];
    z = [-0.003, -0.003, -0.003, -0.003]; % 所有点都在 z = 0.003 高度
    
    % 创建透明平面
    patch(x, y, z, 'blue', 'FaceAlpha', 0.3); % 设置颜色和透明度
    
    amp_filterdis = amp_temp * max(mode_deck(:, VIV_mode_seq(1)));
    zeta_filterdis = zetam - zeta_structure(VIV_mode_seq);
    create_subplot(@scatter3, total_plots, current_plot, {amp_filterdis, U_sel, zeta_filterdis, S, C}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'holdon', holdon);
    xlabel('Amp. (m)')
    ylabel('Wind speed (m/s)')
    title("Damping ratio using filtered dis")
    current_plot = current_plot + 1;
    
    % 设置 colormap
    
    colorbar
    colormap('jet'); % 选择一个颜色映射，比如 'parula'
    % clim([-2 0]); % 设置颜色映射的数据范围为 -3 到 +3
    
    zlim([-0.01, 0])
    xlim([0.01, 0.09])
    ylim([5, 12])
    % 定义平面的四个角的 x, y, 和 z 坐标
    x = [min(amp_filter), max(amp_filter), max(amp_filter), min(amp_filter)];
    y = [min(U_sel), min(U_sel), max(U_sel), max(U_sel)];
    z = [-0.003, -0.003, -0.003, -0.003]; % 所有点都在 z = 0.003 高度
    
    % 创建透明平面
    patch(x, y, z, 'blue', 'FaceAlpha', 0.3); % 设置颜色和透明度
    
end
holdon = true
toc
%% functions
function target = fitnessFunction(params, external_params)
input.lambda_VIV = 10 ^ params(1);
input.sigma_p_VIV = params(2);
input.omega_0_variation_VIV = params(3);
input.Q_value = 10 ^ params(4);
input.sigma_noise = 10 ^ params(5);
input.sigma_buff = 10 ^ params(6);

input.modesel = external_params.modesel;

n = 4;
[result] = viv2013(n, false);
input.start_time = result.startDate;
input.end_time = result.endDate;
input.acc_dir = external_params.acc_dir;
input.VIV_mode_seq = external_params.VIV_mode_seq;
input.nVIV = external_params.nVIV;

% result_Main = KalmanMain(input, 'showtext', false, 'showplot', false, 'filterstyle', 'fft', 'f_keep', 0.33 * [0.9, 1.1]);
[result_Main] = KalmanMain(input, 'showtext', false, 'showplot', false, 'filterstyle', 'nofilter');
fields = fieldnames(result_Main);

logL = result_Main.logL;
logek = result_Main.logek;
%     target = -logL;% 因为 ga 试图最小化函数，所以取负数
target = abs(logek); % 因为 ga 试图最小化函数，所以取负数

end


function target = fitnessFunction_sel(params, external_params)
input.lambda_VIV = 10 ^ params(1);
input.sigma_p_VIV = params(2);
input.Q_value = 10 ^ params(3);
input.sigma_buff = 10 ^ params(4);

input.modesel = external_params.modesel;

n = 4;
[result] = viv2013(n, false);
input.start_time = result.startDate;
input.end_time = result.endDate;
input.acc_dir = external_params.acc_dir;
input.VIV_mode_seq = external_params.VIV_mode_seq;
input.nVIV = external_params.nVIV;
input.omega_0_variation_VIV = external_params.omega_0_variation_VIV;
input.sigma_noise = external_params.sigma_noise;

% result_Main = KalmanMain(input, 'showtext', false, 'showplot', false, 'filterstyle', 'fft', 'f_keep', 0.33 * [0.9, 1.1]);
[result_Main] = KalmanMain(input, 'showtext', false, 'showplot', false, 'filterstyle', 'nofilter');
fields = fieldnames(result_Main);

logL = result_Main.logL;
logek = result_Main.logek;
target = -logL;% 因为 ga 试图最小化函数，所以取负数
% target = abs(logek); % 因为 ga 试图最小化函数，所以取负数

end

function [x_reduced, y_reduced] = reduceDataPoints(x, y, factor)
    % 确保降采样因子是正整数
    factor = max(1, round(factor));
    
    % 进行降采样
    indices = 1:factor:length(x); % 每隔 'factor' 个点取一个点
    x_reduced = x(indices);
    y_reduced = y(indices);
end