% 在被调用的脚本中
if ~exist('skipClear', 'var')
    clc; clear; close all;
end

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
% n = 4;
% [result] = viv2013(n, OFF);
% startDate_global = result.startDate;
% endDate_global = result.endDate;
if ~exist("VIV_sel",'var')
    VIV_sel = 3;
end
opts = detectImportOptions('viv_in_the_paper.csv');
% opts = detectImportOptions('vivData.csv');

% 设置日期时间格式
% 假设日期时间格式为 'MM/dd/yyyy HH:mm'，请根据您的实际情况进行调整
opts = setvartype(opts, 'startDate', 'datetime'); % 确保变量类型为 datetime
opts = setvartype(opts, 'endDate', 'datetime'); % 确保变量类型为 datetime
opts = setvartype(opts, 'startDate_update', 'datetime'); % 确保变量类型为 datetime
opts = setvartype(opts, 'endDate_update', 'datetime'); % 确保变量类型为 datetime

opts = setvaropts(opts, 'startDate', 'InputFormat', 'MM/dd/yyyy HH:mm');
opts = setvaropts(opts, 'endDate', 'InputFormat', 'MM/dd/yyyy HH:mm');
opts = setvaropts(opts, 'startDate_update', 'InputFormat', 'MM/dd/yyyy HH:mm');
opts = setvaropts(opts, 'endDate_update', 'InputFormat', 'MM/dd/yyyy HH:mm');

vivTable = readtable('viv_in_the_paper.csv',opts);
% vivTable = readtable('vivData.csv',opts);

start_time = vivTable.startDate(VIV_sel);
end_time = vivTable.endDate(VIV_sel);

startDate_global = start_time;
endDate_global = end_time;
input_data.start_time = startDate_global;
input_data.end_time = endDate_global;
input_data.VIV_sel = VIV_sel;

% input.lambda_VIV = 10 ^ (-8.725260766057426);
% input.sigma_p_VIV = 4.594524428349437e+04;
% input.omega_0_variation_VIV = 1.001880394311541;
% input.Q_value = 10 ^ (-5.056871357602549);
% input.sigma_noise = 10 ^ (-1.301006929986168);
% input.sigma_buff = 10 ^ (0.070362659875298);


% input_data.lambda_VIV = 10 ^ (-3.942762294993564);
% input_data.sigma_p_VIV = 2.045644004352163e+02;
% input_data.sigma_buff = 10 ^ (0.057229481634271);
% 
% 
% % input_data.Q_value = 10 ^ (-7.743718510318171);
% input_data.Q_value = 10 ^ (-9.023298502893022);
% input_data.omega_0_variation_VIV = 1;
% input_data.sigma_noise = 10 ^ (-1.3);
% % input_data.sigma_noise = 10 ^ (0);


input_data.lambda_VIV = 10 ^ (-3.972197664915089);
input_data.sigma_p_VIV = 19.811759673791514;
input_data.sigma_buff = 10 ^ (0.5411);


% input_data.Q_value = 10 ^ (-7.743718510318171);
input_data.Q_value = 10 ^ (-6);
input_data.omega_0_variation_VIV = 1;
input_data.sigma_noise = 10 ^ (-1.3);
% input_data.sigma_noise = 10 ^ (0);


%% logL优化参数
% modeall = [2, 3, 5, 6, 7, 13, 20, 22, 27, 33, 39, 43, 44, 46, 49, 52, 59, 61, 65, 69, 71, 76, 78, 83, 87, 94, 98];
modeall = [2, 3, 5, 6, 7, 13, 20, 22, 27, 33];
moderemove = [4,9,10];
modeall(moderemove)=[];


modesel = modeall;
if ~exist("modelupdate",'var')
    modelupdate = true;
end

input_data.modelupdate = modelupdate;
% modesel = [2,3,5,6,13,20,22,33];%去掉了两个广义坐标为0的点模态


VIV_mode_seq = find(modesel == 22);

input_data.modesel = modesel;
input_data.VIV_mode_seq = VIV_mode_seq;
nVIV = length(VIV_mode_seq);
input_data.nVIV = nVIV;

showtext = true;
showplot = false;


%% plot logL with changing parameters
if 0

    lambda_VIVs = logspace(-10,-1,20);
    % sigma_p_VIVs = linspace(1e1,1e5,10);
    sigma_p_VIVs = logspace(-1,6,20);

    [lambda_VIVs, sigma_p_VIVs] = meshgrid(lambda_VIVs, sigma_p_VIVs);
    variables = [lambda_VIVs(:), sigma_p_VIVs(:)];

    logLs = zeros(size(variables, 1), 1);
    logSks = zeros(size(variables, 1), 1);
    logeks = zeros(size(variables, 1), 1);


    numIterations = length(variables);

% Instantiate the object with the 'IsParallel' switch set to true
b = ProgressBar(numIterations, ...
    'IsParallel', true, ...
    'WorkerDirectory', pwd(), ...
    'Title', 'Parallel 2' ...
    );

% ALWAYS CALL THE SETUP() METHOD FIRST!!!
b.setup([], [], []);

    parfor k1 = 1:size(variables, 1)
        input_data_temp = input_data;
        input_data_temp.lambda_VIV = variables(k1, 1);
        input_data_temp.sigma_p_VIV = variables(k1, 2);
        input_data_temp.omega_0_variation_VIV = 1;

        input_data_temp.Q_value = 10 ^ (-6);
        input_data_temp.sigma_noise = 10 ^ (-1.3);
        input_data_temp.sigma_buff = 10 ^ (0.541113482743460);
        [result_Main] = KalmanMain(input_data_temp,'shouldFilterYn', true, 'showtext', false, 'showplot', false, 'filterstyle', 'nofilter');
        logLs(k1) = result_Main.logL;
        logSks(k1) = result_Main.logSk;
        logeks(k1) = result_Main.logek;

            % USE THIS FUNCTION AND NOT THE STEP() METHOD OF THE OBJECT!!!
    updateParallel([], pwd);
    end

    b.release();
    
    logLs = reshape(logLs, size(lambda_VIVs));
    logSks = reshape(logSks, size(lambda_VIVs));
    logeks = reshape(logeks, size(lambda_VIVs));

    figure
    surf(log10(lambda_VIVs), log10(sigma_p_VIVs), logLs)
    xlabel('lambda_VIV')
    ylabel('sigma_p_VIV')
    zlabel('logL')
    hold on
    scatter3(log10(10 ^ (-3.942762294993564)),log10(2.045644004352163e+02),9256033.339209)

    figure
    surf(log10(lambda_VIVs), log10(sigma_p_VIVs), logSks)
    xlabel('lambda_VIV')
    ylabel('sigma_p_VIV')
    zlabel('logSk')

    figure
    surf(log10(lambda_VIVs), log10(sigma_p_VIVs), logeks)
    xlabel('lambda_VIV')
    ylabel('sigma_p_VIV')
    zlabel('logek')


    figure
    contourf(log10(lambda_VIVs), log10(sigma_p_VIVs), logLs)
    hold on
    scatter(log10(10^(-3.942762294993564)), log10(2.045644004352163e+02), 'filled');

    return

end



%% Apply Kalman Filter


if 0 %优化所有参数
    %% 导入直接积分获得的涡激力
    ft_directint = importdata("DirectIntegration.mat");
    %% optimization logL to get the maximum with changing lambda sigma_p omega_0_variation Q_value R_value
    % 在调用 ga 函数之前，您可以这样设置 external_params：
    external_params.modesel = modesel;
    external_params.acc_dir = input_data.acc_dir;
    external_params.VIV_mode_seq = VIV_mode_seq;
    external_params.nVIV = nVIV;
    external_params.ft_directint = ft_directint;
    external_params.modelupdate = modelupdate;
    external_params.VIV_sel = VIV_sel;
    external_params.start_time=startDate_global;
    external_params.end_time=endDate_global;

    % 定义参数的范围
    lb = [-4, 1e0, 0.98, -6, -3, 0]; % 这里的值是假设的，请根据您的情况进行修改
    ub = [0, 1e4, 1.02, -4, -1, 2]; % 这里的值也是假设的

    % 定义整数和连续变量
    IntCon = []; % 如果没有整数变量，否则提供整数变量的索引

    options = optimoptions('ga', 'MaxGenerations', 50, 'Display', 'iter', 'UseParallel', true);
    [x, fval] = ga(@(params) fitnessFunction(params, external_params), 6, [], [], [], [], lb, ub, [], IntCon, options);
    % 保存结果
    save('optimization_results.mat', 'x', 'fval');

    input_data.lambda_VIV = 10 ^ (x(1));
    input_data.sigma_p_VIV = x(2);
    input_data.omega_0_variation_VIV = x(3);
    input_data.sigma_buff = 10 ^ (x(6));
    input_data.Q_value = 10 ^ (x(4));

    input_data.sigma_noise = 10 ^ (x(5));



end

if 0 % 选择参数进行优化
    %% 导入直接积分获得的涡激力
    ft_directint = importdata("DirectIntegration.mat");
    %% optimization logL to get the maximum with changing lambda sigma_p omega_0_variation Q_value R_value
    % 在调用 ga 函数之前，您可以这样设置 external_params：
    external_params.modesel = modesel;
    external_params.acc_dir = input_data.acc_dir;
    external_params.VIV_mode_seq = VIV_mode_seq;
    external_params.nVIV = nVIV;
    external_params.ft_directint = ft_directint;
    external_params.omega_0_variation_VIV = input_data.omega_0_variation_VIV;
    external_params.sigma_noise = input_data.sigma_noise;
    external_params.modelupdate = modelupdate;
    external_params.VIV_sel = VIV_sel;
    external_params.start_time=startDate_global;
    external_params.end_time=endDate_global;
    % 定义参数的范围
    lb = [-5, 1e0, -10, 0]; % 这里的值是假设的，请根据您的情况进行修改
    ub = [0, 1e4,  -1,  2]; % 这里的值也是假设的

    % 定义整数和连续变量
    IntCon = []; % 如果没有整数变量，否则提供整数变量的索引

    options = optimoptions('ga', 'MaxGenerations', 30, 'Display', 'iter', 'UseParallel', true);
    [x, fval] = ga(@(params) fitnessFunction_sel(params, external_params), 4, [], [], [], [], lb, ub, [], IntCon, options);
    % 保存结果
    save('optimization_results.mat', 'x', 'fval');

    input_data.lambda_VIV = 10 ^ (x(1));
    input_data.sigma_p_VIV = x(2);
    input_data.Q_value = 10 ^ (x(3));
    input_data.sigma_buff = 10 ^ (x(4));


end


if 0 % 选择参数进行优化
    %% 导入直接积分获得的涡激力
    ft_directint = importdata("DirectIntegration.mat");
    %% optimization logL to get the maximum with changing lambda sigma_p omega_0_variation Q_value R_value
    % 在调用 ga 函数之前，您可以这样设置 external_params：
    external_params.modesel = modesel;
    external_params.acc_dir = input_data.acc_dir;
    external_params.VIV_mode_seq = VIV_mode_seq;
    external_params.nVIV = nVIV;
    external_params.ft_directint = ft_directint;
    external_params.omega_0_variation_VIV = input_data.omega_0_variation_VIV;
    external_params.sigma_noise = input_data.sigma_noise;
    external_params.Q_value = input_data.Q_value;
    external_params.modelupdate = modelupdate;
    external_params.VIV_sel = VIV_sel;
    external_params.start_time=startDate_global;
    external_params.end_time=endDate_global;
    % 定义参数的范围
    lb = [-4, 1e0,  0]; % 这里的值是假设的，请根据您的情况进行修改
    ub = [0, 1e2,  2]; % 这里的值也是假设的

    % 定义整数和连续变量
    IntCon = []; % 如果没有整数变量，否则提供整数变量的索引

    options = optimoptions('ga', 'MaxGenerations', 30, 'Display', 'iter', 'UseParallel', true);
    [x, fval] = ga(@(params) fitnessFunction_sel2(params, external_params), 3, [], [], [], [], lb, ub, [], IntCon, options);
    % 保存结果
    save('optimization_results.mat', 'x', 'fval');

    input_data.lambda_VIV = 10 ^ (x(1));
    input_data.sigma_p_VIV = x(2);
    input_data.sigma_buff = 10 ^ (x(3));


end

if 0 % 选择参数进行优化
    %% 导入直接积分获得的涡激力
    %% optimization logL to get the maximum with changing lambda sigma_p omega_0_variation Q_value R_value
    % 在调用 ga 函数之前，您可以这样设置 external_params：
    input_data.sigma_p_VIV = 55;
    external_params.modesel = modesel;
    external_params.acc_dir = input_data.acc_dir;
    external_params.VIV_mode_seq = VIV_mode_seq;
    external_params.nVIV = nVIV;
    % external_params.ft_directint = ft_directint;
    external_params.omega_0_variation_VIV = input_data.omega_0_variation_VIV;
    external_params.sigma_noise = input_data.sigma_noise;
    external_params.Q_value = input_data.Q_value;
    external_params.sigma_buff = input_data.sigma_buff;
    external_params.sigma_p_VIV = input_data.sigma_p_VIV;
    external_params.modelupdate = modelupdate;
    external_params.VIV_sel = VIV_sel;
    external_params.start_time=startDate_global;
    external_params.end_time=endDate_global;
    % 定义参数的范围
    lb = [-5]; % 这里的值是假设的，请根据您的情况进行修改
    ub = [0]; % 这里的值也是假设的

    % 定义整数和连续变量
    IntCon = []; % 如果没有整数变量，否则提供整数变量的索引

    options = optimoptions('ga', 'MaxGenerations', 30, 'Display', 'iter', 'UseParallel', true);
    [x, fval] = ga(@(params) fitnessFunction_sel3(params, external_params), 1, [], [], [], [], lb, ub, [], IntCon, options);
    % 保存结果
    save('optimization_results_only_lambda.mat', 'x', 'fval');

    input_data.lambda_VIV = 10 ^ (x(1));


end


if 0 % 选择参数进行优化
    %% 导入直接积分获得的涡激力
    %% optimization logL to get the maximum with changing lambda sigma_p omega_0_variation Q_value R_value
    % 在调用 ga 函数之前，您可以这样设置 external_params：
    input_data.lambda_VIV = 10 ^ (-3);
    external_params.modesel = modesel;
    external_params.acc_dir = input_data.acc_dir;
    external_params.VIV_mode_seq = VIV_mode_seq;
    external_params.nVIV = nVIV;
    % external_params.ft_directint = ft_directint;
    external_params.omega_0_variation_VIV = input_data.omega_0_variation_VIV;
    external_params.sigma_noise = input_data.sigma_noise;
    external_params.Q_value = input_data.Q_value;
    external_params.sigma_buff = input_data.sigma_buff;
    external_params.lambda_VIV = input_data.lambda_VIV;
    external_params.VIV_sel = VIV_sel;
    external_params.start_time=startDate_global;
    external_params.end_time=endDate_global;
    % 定义参数的范围
    lb = [10^0]; % 这里的值是假设的，请根据您的情况进行修改
    ub = [10^3]; % 这里的值也是假设的

    % 定义整数和连续变量
    IntCon = []; % 如果没有整数变量，否则提供整数变量的索引

    options = optimoptions('ga', 'MaxGenerations', 30, 'Display', 'iter', 'UseParallel', true);
    [x, fval] = ga(@(params) fitnessFunction_sel4(params, external_params), 1, [], [], [], [], lb, ub, [], IntCon, options);
    % 保存结果
    save('optimization_results_only_sigma.mat', 'x', 'fval');

    input_data.sigma_p_VIV = x(1);


end

%% Apply kalman filter
% [result_Main] = KalmanMain(input, 'showtext', true, 'showplot', false, 'filterstyle', 'fft', 'f_keep', 0.33 * [0.9, 1.1]);

[result_Main] = KalmanMain(input_data, 'shouldFilterYn', true,'showtext', showtext, 'showplot', showplot, 'filterstyle', 'nofilter');
logL = result_Main.logL;
logSk = result_Main.logSk;
logek = result_Main.logek;
% display logL logSk logek
dispstr = sprintf("logL = %f, logSk = %f, logek = %f", logL, logSk, logek);

if showtext
    disp(dispstr)
end

mode_deck = result_Main.mode_deck;

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


start_time.Format = 'yyyy_MM_dd_HH_mm_ss';
formatted_start_time = string(start_time);

end_time.Format = 'yyyy_MM_dd_HH_mm_ss';
formatted_end_time = string(end_time);


formatted_start_time = strrep(strrep(strrep(formatted_start_time, ':', '-'), ' ', '_'), '/', '-');
formatted_end_time = strrep(strrep(strrep(formatted_end_time, ':', '-'), ' ', '_'), '/', '-');

filename = "windspeed_result_"+formatted_start_time+"_"+formatted_end_time+".mat";

if exist(filename, 'file') == 2
    load(filename);
else
    windspeed_result = read_wind_data_onetime(t_cycle_mean_temp, duration, input_data.wind_dir_all, [1 1 1 1 1 1]);
    save(filename,'windspeed_result',"start_time","end_time");
end

%% Damping Ratio calculation
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


    end

end

%% wind speed with amplitude and damping ratio

UAa = windspeed_result.UA1;
UAb = windspeed_result.UA2;

beta_deg_mean_UAa = UAa.beta_deg_mean;
beta_deg_mean_UAb = UAb.beta_deg_mean;

% judge if the wind direction of both UA5 and UA6 is from 45°-225°
for k1 = 1:length(beta_deg_mean_UAa)

    if and(beta_deg_mean_UAa(k1) > 45, beta_deg_mean_UAa(k1) < 225) && and(beta_deg_mean_UAb(k1) > 45, beta_deg_mean_UAb(k1) < 225)
        U_sel(k1) = UAa.U(k1);
        AoA_sel(k1) = UAa.alpha_deg_mean(k1);
        TI_u_sel(k1) = UAa.TI_u(k1);
        TI_v_sel(k1) = UAa.TI_v(k1);
        TI_w_sel(k1) = UAa.TI_w(k1);


        % judge if the wind direction of both UA5 and UA6 is from 225°-360° or 0°-45°
    elseif or(and(beta_deg_mean_UAa(k1) > 225, beta_deg_mean_UAa(k1) < 360), and(beta_deg_mean_UAa(k1) > 0, beta_deg_mean_UAa(k1) < 45)) && or(and(beta_deg_mean_UAb(k1) > 225, beta_deg_mean_UAb(k1) < 360), and(beta_deg_mean_UAb(k1) > 0, beta_deg_mean_UAb(k1) < 45))
        U_sel(k1) = UAb.U(k1);
        AoA_sel(k1) = UAb.alpha_deg_mean(k1);
        TI_u_sel(k1) = UAb.TI_u(k1);
        TI_v_sel(k1) = UAb.TI_v(k1);
        TI_w_sel(k1) = UAb.TI_w(k1);
    else

        if abs(UAa.alpha_deg_mean(k1)) < abs(UAb.alpha_deg_mean(k1))
            U_sel(k1) = UAa.U(k1);
            AoA_sel(k1) = UAa.alpha_deg_mean(k1);
            TI_u_sel(k1) = UAa.TI_u(k1);
            TI_v_sel(k1) = UAa.TI_v(k1);
            TI_w_sel(k1) = UAa.TI_w(k1);
        else
            U_sel(k1) = UAb.U(k1);
            AoA_sel(k1) = UAb.alpha_deg_mean(k1);
            TI_u_sel(k1) = UAb.TI_u(k1);
            TI_v_sel(k1) = UAb.TI_v(k1);
            TI_w_sel(k1) = UAb.TI_w(k1);
        end

        disp("该时间点两个风速仪风向不一致，取风攻角较小的值")
    end

end

U_sel_loc_1 = U_sel;
AoA_sel_1 = AoA_sel;
beta_deg_mean_UA1a=beta_deg_mean_UAa;
beta_deg_mean_UA1b=beta_deg_mean_UAb;

TI_u_sel_1 = TI_u_sel;
TI_v_sel_1 = TI_v_sel;
TI_w_sel_1 = TI_w_sel;


UAa = windspeed_result.UA3;
UAb = windspeed_result.UA4;

beta_deg_mean_UAa = UAa.beta_deg_mean;
beta_deg_mean_UAb = UAb.beta_deg_mean;

% judge if the wind direction of both UA5 and UA6 is from 45°-225°
for k1 = 1:length(beta_deg_mean_UAa)

    if and(beta_deg_mean_UAa(k1) > 45, beta_deg_mean_UAa(k1) < 225) && and(beta_deg_mean_UAb(k1) > 45, beta_deg_mean_UAb(k1) < 225)
        U_sel(k1) = UAa.U(k1);
        AoA_sel(k1) = UAa.alpha_deg_mean(k1);
        TI_u_sel(k1) = UAa.TI_u(k1);
        TI_v_sel(k1) = UAa.TI_v(k1);
        TI_w_sel(k1) = UAa.TI_w(k1);


        % judge if the wind direction of both UA5 and UA6 is from 225°-360° or 0°-45°
    elseif or(and(beta_deg_mean_UAa(k1) > 225, beta_deg_mean_UAa(k1) < 360), and(beta_deg_mean_UAa(k1) > 0, beta_deg_mean_UAa(k1) < 45)) && or(and(beta_deg_mean_UAb(k1) > 225, beta_deg_mean_UAb(k1) < 360), and(beta_deg_mean_UAb(k1) > 0, beta_deg_mean_UAb(k1) < 45))
        U_sel(k1) = UAb.U(k1);
        AoA_sel(k1) = UAb.alpha_deg_mean(k1);
        TI_u_sel(k1) = UAb.TI_u(k1);
        TI_v_sel(k1) = UAb.TI_v(k1);
        TI_w_sel(k1) = UAb.TI_w(k1);
    else

        if abs(UAa.alpha_deg_mean(k1)) < abs(UAb.alpha_deg_mean)
            U_sel(k1) = UAa.U(k1);
            AoA_sel(K1) = UAa.alpha_deg_mean(k1);
            TI_u_sel(k1) = UAa.TI_u(k1);
            TI_v_sel(k1) = UAa.TI_v(k1);
            TI_w_sel(k1) = UAa.TI_w(k1);
        else
            U_sel(k1) = UAb.U(k1);
            AoA_sel(K1) = UAb.alpha_deg_mean(k1);
            TI_u_sel(k1) = UAb.TI_u(k1);
            TI_v_sel(k1) = UAb.TI_v(k1);
            TI_w_sel(k1) = UAb.TI_w(k1);
        end

        disp("该时间点两个风速仪风向不一致，取风攻角较小的值")
    end

end


U_sel_loc_2 = U_sel;
AoA_sel_2 = AoA_sel;
beta_deg_mean_UA2a=beta_deg_mean_UAa;
beta_deg_mean_UA2b=beta_deg_mean_UAb;

TI_u_sel_2 = TI_u_sel;
TI_v_sel_2 = TI_v_sel;
TI_w_sel_2 = TI_w_sel;

UAa = windspeed_result.UA5;
UAb = windspeed_result.UA6;

beta_deg_mean_UAa = UAa.beta_deg_mean;
beta_deg_mean_UAb = UAb.beta_deg_mean;

% judge if the wind direction of both UA5 and UA6 is from 45°-225°
for k1 = 1:length(beta_deg_mean_UAa)

    if and(beta_deg_mean_UAa(k1) > 45, beta_deg_mean_UAa(k1) < 225) && and(beta_deg_mean_UAb(k1) > 45, beta_deg_mean_UAb(k1) < 225)
        U_sel(k1) = UAa.U(k1);
        AoA_sel(k1) = UAa.alpha_deg_mean(k1);
        TI_u_sel(k1) = UAa.TI_u(k1);
        TI_v_sel(k1) = UAa.TI_v(k1);
        TI_w_sel(k1) = UAa.TI_w(k1);


        % judge if the wind direction of both UA5 and UA6 is from 225°-360° or 0°-45°
    elseif or(and(beta_deg_mean_UAa(k1) > 225, beta_deg_mean_UAa(k1) < 360), and(beta_deg_mean_UAa(k1) > 0, beta_deg_mean_UAa(k1) < 45)) && or(and(beta_deg_mean_UAb(k1) > 225, beta_deg_mean_UAb(k1) < 360), and(beta_deg_mean_UAb(k1) > 0, beta_deg_mean_UAb(k1) < 45))
        U_sel(k1) = UAb.U(k1);
        AoA_sel(k1) = UAb.alpha_deg_mean(k1);
        TI_u_sel(k1) = UAb.TI_u(k1);
        TI_v_sel(k1) = UAb.TI_v(k1);
        TI_w_sel(k1) = UAb.TI_w(k1);
    else

        if abs(UAa.alpha_deg_mean(k1)) < abs(UAb.alpha_deg_mean)
            U_sel(k1) = UAa.U(k1);
            AoA_sel(K1) = UAa.alpha_deg_mean(k1);
            TI_u_sel(k1) = UAa.TI_u(k1);
            TI_v_sel(k1) = UAa.TI_v(k1);
            TI_w_sel(k1) = UAa.TI_w(k1);
        else
            U_sel(k1) = UAb.U(k1);
            AoA_sel(K1) = UAb.alpha_deg_mean(k1);
            TI_u_sel(k1) = UAb.TI_u(k1);
            TI_v_sel(k1) = UAb.TI_v(k1);
            TI_w_sel(k1) = UAb.TI_w(k1);
        end

        disp("该时间点两个风速仪风向不一致，取风攻角较小的值")
    end

end


U_sel_loc_3 = U_sel;
AoA_sel_3 = AoA_sel;
beta_deg_mean_UA3a=beta_deg_mean_UAa;
beta_deg_mean_UA3b=beta_deg_mean_UAb;

TI_u_sel_3 = TI_u_sel;
TI_v_sel_3 = TI_v_sel;
TI_w_sel_3 = TI_w_sel;

amp_filter = amp_cell{1}{1} * max(mode_deck(:, VIV_mode_seq(1)));
zeta_filter = zeta_all_cell{1}{1};
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
        title("modal force for " + "mode"+modesel(VIV_mode_seq));
        current_plot = current_plot + 1;
        ylim([-100, 100])
    end
    [f, magnitude] = fft_transform(50,F_filter(1, :));
    create_subplot(@plot, total_plots, current_plot, {f, magnitude}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', newfigure, 'holdon', holdon);
    title("Force in the frequency domain")
    xlabel("Frequency")
    ylabel("Amplitude")
    current_plot = current_plot + 1;
    save temp2.mat t F_filter
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


    create_subplot(@scatter, total_plots, current_plot, {amp_temp * max(mode_deck(:, VIV_mode_seq)), zetam - zeta_structure(VIV_mode_seq)}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', false, 'holdon', holdon);
    hold on
    create_subplot(@scatter, total_plots, current_plot, {amp_cell{1}{1} * max(mode_deck(:, VIV_mode_seq)), zeta_all_cell{1}{1}, [], secondsFromReference, 'filled'}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', false, 'holdon', holdon);

    current_plot = current_plot + 1;

    % 设置 colormap
    colormap('jet')
    colorbar

    hold on
    plot([0, 0.15], [-0.003, -0.003])
    % scatter(ex,epsx,'green')
    str = "Mode : %d, Frequency : %.2f Hz";
    title(sprintf(str, modesel(k1), top_freqs{k1}(k2)));
    xlim([0, 0.12])
    ylim([-0.5, 0.5] / 100)
    xlabel("Amplitude(m)")
    ylabel("Damping ratio")
    

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

    %% Damping ratio calculation with wind speed
    create_subplot(@scatter, total_plots, current_plot, {amp_cell{1}{1} * max(mode_deck(:, VIV_mode_seq)), zeta_all_cell{1}{1}, [], U_sel, 'filled'}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', false, 'holdon', holdon);

    current_plot = current_plot + 1;

    % 设置 colormap
    colormap('jet')
    colorbar

    hold on
    plot([0, 0.15], [-0.003, -0.003])
    % scatter(ex,epsx,'green')
    str = "Mode : %d, Frequency : %.2f Hz with wind speed";
    % title(sprintf(str, modesel(k1), top_freqs{k1}(k2)));
    xlim([0.00, 0.12])
    ylim([-0.5, 0.5] / 100)
    xlabel("Amplitude(m)")
    ylabel("Damping ratio")
    %% Damping ratio calculation with wind speed
    create_subplot(@scatter, total_plots, current_plot, {amp_cell{1}{1} * max(mode_deck(:, VIV_mode_seq)), zeta_all_cell{1}{1}, [], AoA_sel, 'filled'}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'firstfigure', false, 'holdon', holdon);

    current_plot = current_plot + 1;

    % 设置 colormap
    colormap('jet')
    colorbar

    hold on
    plot([0, 0.15], [-0.003, -0.003])
    % scatter(ex,epsx,'green')
    str = "Mode : %d, Frequency : %.2f Hz with attack of angle";
    title(sprintf(str, modesel(k1), top_freqs{k1}(k2)));
    xlim([0.00, 0.12])
    ylim([-0.5, 0.5] / 100)
    xlabel("Amplitude(m)")
    ylabel("Damping ratio")


    
    % C=rescale(AoA_sel,10,100);
    S = rescale(zeta_filter, 10, 100);
    C = AoA_sel;
    % C = secondsFromReference;
    % scatter3(amp_filter,U_sel,zeta_filter,S,C)

    create_subplot(@scatter3, total_plots, current_plot, {amp_filter, U_sel,  zeta_filter,[],C}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'newfigure', newfigure, 'holdon', holdon);
    xlabel('Amp. (m)')
    ylabel('Wind speed (m/s)')
    title("Damping ratio （Amp,U,damping ratio,,AOA）")
    current_plot = current_plot + 1;

    % 设置 colormap

    colorbar
    colormap('jet'); % 选择一个颜色映射，比如 'parula'
    % clim([-2 0]); % 设置颜色映射的数据范围为 -3 到 +3

    % zlim([-0.01, 0])
    % xlim([0.01, 0.09])
    % ylim([5, 12])
    % 定义平面的四个角的 x, y, 和 z 坐标
    x = [min(amp_filter), max(amp_filter), max(amp_filter), min(amp_filter)];
    y = [min(U_sel), min(U_sel), max(U_sel), max(U_sel)];
    z = [-0.003, -0.003, -0.003, -0.003]; % 所有点都在 z = 0.003 高度

    % 创建透明平面
    patch(x, y, z, 'blue', 'FaceAlpha', 0.3); % 设置颜色和透明度
    holdon = true
    toc


end
%% save the result
% eg: result_starttime_VIV_sel.mat VIV_sel is a variable
startDate_global.Format = 'yyyy_MM_dd_HH_mm';
zeta = result_Main.zeta;

if input_data.modelupdate == true
    str_modelupdate = "updatedmodel";
else
    str_modelupdate = "nomodelupdate";
end
filename = sprintf('result_%s_%s_%s.mat', startDate_global, num2str(VIV_sel),str_modelupdate);
save(filename, 'VIV_sel',"amp_filter","zeta_filter","datetimeArray","VIV_mode_seq","zeta","U_sel_loc_1","AoA_sel_1","U_sel_loc_2","AoA_sel_2","U_sel_loc_3","AoA_sel_3","beta_deg_mean_UA1a","beta_deg_mean_UA1b","beta_deg_mean_UA2a","beta_deg_mean_UA2b","beta_deg_mean_UA3a","beta_deg_mean_UA3b","t_cycle_mean_temp"...
    ,"TI_u_sel_1","TI_u_sel_2","TI_u_sel_3","TI_v_sel_1","TI_v_sel_2","TI_v_sel_3","TI_w_sel_1","TI_w_sel_2","TI_w_sel_3")

return
%% plot state estimate
close all
modal_state = result_Main.x_k_k(1:(end-2)/2,:);



for k1 = 1:size(modal_state,1)
    figure
    plot(t,modal_state(k1,:)*max(abs(result_Main.mode_deck(:,k1))))

    [f, magnitude] = fft_transform(50,modal_state(k1,:)*max(abs(result_Main.mode_deck(:,k1))));
    figure
    semilogy(f, magnitude)
    xlim([0,1])
    ylim([1e-10,0.5])
end





tilefigs([4 14],'',[],[100 0],[100 150],[100 100])
    %% functions
    function target = fitnessFunction(params, external_params)
    input.lambda_VIV = 10 ^ params(1);
    input.sigma_p_VIV = params(2);
    input.omega_0_variation_VIV = params(3);
    input.Q_value = 10 ^ params(4);
    input.sigma_noise = 10 ^ params(5);
    input.sigma_buff = 10 ^ params(6);
    input.modelupdate =external_params.modelupdate;
    input.modesel = external_params.modesel;
    input.VIV_sel=external_params.VIV_sel;
    % n = 4;
    % [result] = viv2013(n, false);
    input.start_time =    external_params.start_time;
    input.end_time = external_params.end_time;
    input.acc_dir = external_params.acc_dir;
    input.VIV_mode_seq = external_params.VIV_mode_seq;
    input.nVIV = external_params.nVIV;

    % result_Main = KalmanMain(input, 'showtext', false, 'showplot', false, 'filterstyle', 'fft', 'f_keep', 0.33 * [0.9, 1.1]);
    [result_Main] = KalmanMain(input, 'showtext', false, 'showplot', false, 'filterstyle', 'nofilter');
    fields = fieldnames(result_Main);

    logL = result_Main.logL;
    logek = result_Main.logek;
    target = -logL;% 因为 ga 试图最小化函数，所以取负数
    % target = abs(logek); % 因为 ga 试图最小化函数，所以取负数

    end


    function target = fitnessFunction_sel(params, external_params)
    input.lambda_VIV = 10 ^ params(1);
    input.sigma_p_VIV = params(2);
    input.Q_value = 10 ^ params(3);
    input.sigma_buff = 10 ^ params(4);
    input.modelupdate =external_params.modelupdate;
    input.modesel = external_params.modesel;
    input.VIV_sel=external_params.VIV_sel;
    % n = 4;
    % [result] = viv2013(n, false);
    input.start_time =    external_params.start_time;
    input.end_time = external_params.end_time;
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

    function target = fitnessFunction_sel2(params, external_params)
    input.lambda_VIV = 10 ^ params(1);
    input.sigma_p_VIV = params(2);
    input.sigma_buff = 10 ^ params(3);
    input.modelupdate =external_params.modelupdate;
    input.modesel = external_params.modesel;
    input.VIV_sel=external_params.VIV_sel;
    % n = 4;
    % [result] = viv2013(n, false);
    input.start_time =    external_params.start_time;
    input.end_time = external_params.end_time;
    input.acc_dir = external_params.acc_dir;
    input.VIV_mode_seq = external_params.VIV_mode_seq;
    input.nVIV = external_params.nVIV;
    input.omega_0_variation_VIV = external_params.omega_0_variation_VIV;
    input.sigma_noise = external_params.sigma_noise;
    input.Q_value =external_params.Q_value ;
    % result_Main = KalmanMain(input, 'showtext', false, 'showplot', false, 'filterstyle', 'fft', 'f_keep', 0.33 * [0.9, 1.1]);
    [result_Main] = KalmanMain(input, 'showtext', false, 'showplot', false, 'filterstyle', 'nofilter');
    fields = fieldnames(result_Main);

    logL = result_Main.logL;
    logek = result_Main.logek;
    target = -logL;% 因为 ga 试图最小化函数，所以取负数
    % target = abs(logek); % 因为 ga 试图最小化函数，所以取负数

    end


    function target = fitnessFunction_sel3(params, external_params)
    input.lambda_VIV = 10 ^ params(1);
    input.modelupdate =external_params.modelupdate;
    input.modesel = external_params.modesel;
    input.VIV_sel=external_params.VIV_sel;
    % n = 4;
    % [result] = viv2013(n, false);
    input.start_time =    external_params.start_time;
    input.end_time = external_params.end_time;
    input.acc_dir = external_params.acc_dir;
    input.VIV_mode_seq = external_params.VIV_mode_seq;
    input.nVIV = external_params.nVIV;
    input.omega_0_variation_VIV = external_params.omega_0_variation_VIV;
    input.sigma_noise = external_params.sigma_noise;
    input.Q_value =external_params.Q_value ;
    input.sigma_p_VIV =external_params.sigma_p_VIV ;
    input.sigma_buff =external_params.sigma_buff ;
    % result_Main = KalmanMain(input, 'showtext', false, 'showplot', false, 'filterstyle', 'fft', 'f_keep', 0.33 * [0.9, 1.1]);
    [result_Main] = KalmanMain(input, 'showtext', false, 'showplot', false, 'filterstyle', 'nofilter');
    fields = fieldnames(result_Main);

    logL = result_Main.logL;
    logek = result_Main.logek;
    target = -logL;% 因为 ga 试图最小化函数，所以取负数
    % target = abs(logek); % 因为 ga 试图最小化函数，所以取负数

    end

    function target = fitnessFunction_sel4(params, external_params)
    input.sigma_p_VIV = params(1);
    input.modelupdate =external_params.modelupdate;
    input.modesel = external_params.modesel;
    input.VIV_sel=external_params.VIV_sel;
    % n = 4;
    % % [result] = viv2013(n, false);
    input.start_time =    external_params.start_time;
    input.end_time = external_params.end_time;
    input.acc_dir = external_params.acc_dir;
    input.VIV_mode_seq = external_params.VIV_mode_seq;
    input.nVIV = external_params.nVIV;
    input.omega_0_variation_VIV = external_params.omega_0_variation_VIV;
    input.sigma_noise = external_params.sigma_noise;
    input.Q_value =external_params.Q_value ;
    input.lambda_VIV =external_params.lambda_VIV ;
    input.sigma_buff =external_params.sigma_buff ;
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