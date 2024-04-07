%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: ShengyiXu xushengyichn@outlook.com
%Date: 2023-10-09 22:23:15
%LastEditors: xushengyichn xushengyichn@outlook.com
%LastEditTime: 2024-03-25 22:11:02
%FilePath: \Exercises-for-Techniques-for-estimation-in-dynamics-systemsd:\git\xihoumen_inverse_force_estimation\20240320 fourth edition\KalmanMain.m
%Description: 加上更多模态，不要只留下单一模态，看看能不能起到滤波的作用
%
%Copyright (c) 2023 by ${git_name_email}, All Rights Reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [result_Main] = KalmanMain(input,varargin)

if nargin == 0
    % clc; clear; close all;
    % addpath(genpath("F:\git\ssm_tools\"))
    % addpath(genpath("F:\git\Function_shengyi_package\"))
    % addpath(genpath("F:\git\xihoumen_inverse_force_estimation\FEM_model\"))
    % addpath(genpath("C:\Users\xushengyi\Documents\Github\"))
    %
    %
    % subStreamNumberDefault = 2132;
    %
    
    % input.ON = params.ON;
    % input.OFF = params.OFF;
    %% 0 绘图参数
    
    % input.num_figs_in_row = 12; %每一行显示几个图
    % input.figPos = params.figPosSmall; %图的大小，参数基于InitScript.m中的设置
    %设置图片间隔
    % input.gap_between_images = [0, 0];
    % input.figureIdx = 0;
    % input.fig_bool = params.ON;
    % input.displayText = params.ON;
    run('CommonCommand.m');
    params = Init_fun();
    n = 4;
    [result] = viv2013(n, params.OFF);
    input.start_time = result.startDate;
    input.end_time = result.endDate;
    
    % input.start_time = datetime('2013-02-06 01:30:00', 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
    % input.end_time = datetime('2013-02-06 01:45:59', 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
    
    
    % input.acc_dir = "F:\test\result";
    
    % input.acc_dir = "Z:\Drive\Backup\SHENGYI_HP\F\test\result";
    
    
    input.lambda = 10 ^ (-4.987547778158018);
    input.sigma_p = 1.411528858719115e+04;
    % input.omega_0_variation =1;
    % input.Q_value = 10 ^ (-8);
    % input.R_value = 10 ^ (-6);
    input.Q_value = 10 ^ (-1);
    % input.R_value = 10 ^ (-8);
    
    input.sigma_buff = 10;
    input.sigma_noise = 10e-4;
    
    input.lambda_VIV = 10 ^ (-4.987547778158018);
    input.sigma_p_VIV = 1.411528858719115e+04;
    input.omega_0_variation_VIV =1;
    
    % input.lambda_matern = 10 ^ (-4.987547778158018);
    % input.sigma_p_matern = 1.411528858719115e+04;
    
    % input.modesel= [2,3,5,6,7,9,15,21,23,29,33,39,44,45];
    input.modesel= [23,45];
    input.VIV_mode_seq = find(input.modesel ==23);
    % KalmanMain(input,'showtext', false,'showplot',true,'shouldFilterYn', true,'shouldFilterp_filt_m', true);
    KalmanMain(input,'showtext', false,'showplot',false)
    return;
end
p = inputParser;
addParameter(p, 'showtext', true, @islogical)
addParameter(p, 'shouldFilterYn', false, @islogical)
addParameter(p, 'shouldFilterp_filt_m', false, @islogical)
addParameter(p, 'showplot', true, @islogical)
addParameter(p, 'num_figs_in_row', 12, @isnumeric)
addParameter(p, 'figPos', [100,100,400,300], @isnumeric)
addParameter(p, 'gap_between_images', [0,0], @isnumeric)
addParameter(p, 'figureIdx', 0, @isnumeric)
addParameter(p, 'f_keep', 0.33*[0.9,1.1], @isnumeric)
addParameter(p,'filterstyle','nofilter',@ischar);% fft use the function by myself, bandpass use the function in matlab. fft is much faster than bandpass but may not be accurate
parse(p, varargin{:});
filterstyle = p.Results.filterstyle;
showtext = p.Results.showtext;
shouldFilterYn = p.Results.shouldFilterYn;
shouldFilterp_filt_m = p.Results.shouldFilterp_filt_m;
showplot = p.Results.showplot;
f_keep = p.Results.f_keep;
num_figs_in_row = p.Results.num_figs_in_row;
figPos = p.Results.figPos;
gap_between_images = p.Results.gap_between_images;
figureIdx = p.Results.figureIdx;

if showtext
    showtext_char = 'yes';
else
    showtext_char = 'no';
end

%% 1 读取数据
start_time = input.start_time;
end_time = input.end_time;

acc_dir = input.acc_dir;

lambda_VIV = input.lambda_VIV;
sigma_p_VIV = input.sigma_p_VIV;
omega_0_variation_VIV = input.omega_0_variation_VIV;

% lambda_matern_noVIV = input.lambda_matern;
% sigma_p_matern_noVIV = input.sigma_p_matern;

modesel= input.modesel;
VIV_mode_seq = input.VIV_mode_seq;

sigma_buff = input.sigma_buff;
sigma_noise = input.sigma_noise;

% lambda = input.lambda;
% sigma_p = input.sigma_p;
% omega_0_variation = input.omega_0_variation;
Q_value = input.Q_value;
% R_value = input.R_value;
R_value = sigma_noise^2;




fig_bool = showplot;
ON = true;
OFF = false;

[Acc_Data] = read_acceleration_data(start_time, end_time, acc_dir);

switch filterstyle
    case 'fft'
        timeDifferences = diff(Acc_Data.mergedData.Time); % Returns a duration array
        dt = seconds(timeDifferences(1)); % Converts the first duration value to seconds and assigns to dt
        fs = 1 / dt;
        Acc_Data.mergedData.AC2_1 = fft_filter(fs, Acc_Data.mergedData.AC2_1, f_keep)';
        Acc_Data.mergedData.AC2_2 = fft_filter(fs, Acc_Data.mergedData.AC2_2, f_keep)';
        Acc_Data.mergedData.AC2_3 = fft_filter(fs, Acc_Data.mergedData.AC2_3, f_keep)';
        Acc_Data.mergedData.AC3_1 = fft_filter(fs, Acc_Data.mergedData.AC3_1, f_keep)';
        Acc_Data.mergedData.AC3_2 = fft_filter(fs, Acc_Data.mergedData.AC3_2, f_keep)';
        Acc_Data.mergedData.AC3_3 = fft_filter(fs, Acc_Data.mergedData.AC3_3, f_keep)';
        Acc_Data.mergedData.AC4_1 = fft_filter(fs, Acc_Data.mergedData.AC4_1, f_keep)';
        Acc_Data.mergedData.AC4_2 = fft_filter(fs, Acc_Data.mergedData.AC4_2, f_keep)';
        Acc_Data.mergedData.AC4_3 = fft_filter(fs, Acc_Data.mergedData.AC4_3, f_keep)';
        
        
        
    case 'bandpass'
        timeDifferences = diff(Acc_Data.mergedData.Time); % Returns a duration array
        dt = seconds(timeDifferences(1)); % Converts the first duration value to seconds and assigns to dt
        fs = 1 / dt;
        Acc_Data.mergedData.AC2_1 = bandpass(Acc_Data.mergedData.AC2_1', f_keep, fs)';
        Acc_Data.mergedData.AC2_2 = bandpass(Acc_Data.mergedData.AC2_2', f_keep, fs)';
        Acc_Data.mergedData.AC2_3 = bandpass(Acc_Data.mergedData.AC2_3', f_keep, fs)';
        Acc_Data.mergedData.AC3_1 = bandpass(Acc_Data.mergedData.AC3_1', f_keep, fs)';
        Acc_Data.mergedData.AC3_2 = bandpass(Acc_Data.mergedData.AC3_2', f_keep, fs)';
        Acc_Data.mergedData.AC3_3 = bandpass(Acc_Data.mergedData.AC3_3', f_keep, fs)';
        Acc_Data.mergedData.AC4_1 = bandpass(Acc_Data.mergedData.AC4_1', f_keep, fs)';
        Acc_Data.mergedData.AC4_2 = bandpass(Acc_Data.mergedData.AC4_2', f_keep, fs)';
        Acc_Data.mergedData.AC4_3 = bandpass(Acc_Data.mergedData.AC4_3', f_keep, fs)';

        
    case 'nofilter'
        
        
    otherwise
        error('Invalid filter style. Choose either "fft" or "bandpass".')
end

[uniqueTimestamps, ia, ~] = unique(Acc_Data.mergedData.Time, 'stable');

% Find duplicates by checking the difference in lengths
if showtext
    if length(uniqueTimestamps) < height(Acc_Data.mergedData)
        disp('There are duplicate timestamps:');
        
        % Get duplicated indices
        duplicatedIndices = setdiff(1:height(Acc_Data.mergedData), ia);
        
        % Display duplicates
        for i = 1:length(duplicatedIndices)
            disp(Acc_Data.mergedData.Time(duplicatedIndices(i)));
        end
        
    else
        disp('No duplicates found.');
    end
end


% select which wind profile to use

%% 判断涡振模态
% acc_1 = Acc_Data.mergedData.AC3_1;
% t = Acc_Data.mergedData.Time;
% % [p, fd, td] = pspectrum(acc_1, t, 'spectrogram', 'FrequencyResolution', 0.005);
% % instfreq(p, fd, td);
%
% % fs = 50;
% % [f, magnitude] = fft_transform(fs, acc_1);
% % figure
% % plot(f, magnitude)


%% 2 有限元模型
% 读入ANSYS梁桥模型质量刚度矩阵  MCK矩阵 Import MCK matrix from ANSYS
% 将ANSYS中的稀疏矩阵处理为完全矩阵 Handling sparse matrices in ANSYS as full matrices

% modesel= [2,3,5,6,7,9,15,21,23,29,33,39,44,45];
% modesel = [23];
nmodes = length(modesel); ns = nmodes * 2;
Result = ImportMK(nmodes, 'KMatrix.matrix', 'MMatrix.matrix', 'nodeondeck.txt', 'KMatrix.mapping', 'nodegap.txt', 'modesel', modesel,'showtext',showtext);
mode_deck = Result.mode_deck; mode_deck_re = Result.mode_deck_re; node_loc = Result.node_loc;
Freq = Result.Freq;
MM_eq = Result.MM_eq; KK_eq = Result.KK_eq;

% zeta = ones(size(modesel)) * 0.3/100;
% zeta = ones(size(modesel)) * 0.0/100;
% zeta = importdata("zeta_update.mat");
zeta = importdata("zeta_update0406.mat");

% opts = detectImportOptions('viv_in_the_paper.csv');
% % opts = detectImportOptions('vivData.csv');
% % 设置日期时间格式
% % 假设日期时间格式为 'MM/dd/yyyy HH:mm'，请根据您的实际情况进行调整
% opts = setvartype(opts, 'startDate', 'datetime'); % 确保变量类型为 datetime
% opts = setvartype(opts, 'endDate', 'datetime'); % 确保变量类型为 datetime
% opts = setvartype(opts, 'startDate_update', 'datetime'); % 确保变量类型为 datetime
% opts = setvartype(opts, 'endDate_update', 'datetime'); % 确保变量类型为 datetime
% opts = setvaropts(opts, 'startDate', 'InputFormat', 'MM/dd/yyyy HH:mm');
% opts = setvaropts(opts, 'endDate', 'InputFormat', 'MM/dd/yyyy HH:mm');
% opts = setvaropts(opts, 'startDate_update', 'InputFormat', 'MM/dd/yyyy HH:mm');
% opts = setvaropts(opts, 'endDate_update', 'InputFormat', 'MM/dd/yyyy HH:mm');
% 
% vivTable = readtable('viv_in_the_paper.csv',opts);

opts = delimitedTextImportOptions("NumVariables", 25);

% Specify range and delimiter
opts.DataLines = [2, Inf];
opts.Delimiter = ",";

% Specify column names and types
opts.VariableNames = ["caseNumber", "startDate", "endDate", "Fre", "startDate_update", "endDate_update", "bad_data", "k_value", "sensor_selection", "delay2", "delay3", "mode1_SSI", "mode1_seq", "mode2_SSI", "mode2_seq", "mode3_SSI", "mode3_seq", "mode4_SSI", "mode4_seq", "mode5_SSI", "mode5_seq", "mode6_SSI", "mode6_seq", "mode7_SSI", "mode7_seq"];
opts.VariableTypes = ["double", "datetime", "datetime", "double", "datetime", "datetime", "double", "double", "categorical", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double", "double"];

% Specify file level properties
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule = "read";

% Specify variable properties
opts = setvaropts(opts, "sensor_selection", "EmptyFieldRule", "auto");
opts = setvaropts(opts, "startDate", "InputFormat", "MM/dd/yyyy HH:mm", "DatetimeFormat", "preserveinput");
opts = setvaropts(opts, "endDate", "InputFormat", "MM/dd/yyyy HH:mm", "DatetimeFormat", "preserveinput");
opts = setvaropts(opts, "startDate_update", "InputFormat", "MM/dd/yyyy HH:mm", "DatetimeFormat", "preserveinput");
opts = setvaropts(opts, "endDate_update", "InputFormat", "MM/dd/yyyy HH:mm", "DatetimeFormat", "preserveinput");

% Import the data
vivTable = readtable('viv_in_the_paper.csv', opts);

%% modal updating
if input.modelupdate

    % vivTable = readtable('vivData.csv',opts);
    
    start_time = vivTable.startDate_update(input.VIV_sel);
    end_time = vivTable.endDate_update(input.VIV_sel);
    % start_time = datetime('2013-02-04 23:15:00', 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
    % end_time = datetime('2013-02-04 23:45:00', 'InputFormat', 'yyyy-MM-dd HH:mm:ss'); % Example time range
    
    % 替换日期时间字符串中的冒号和其他特殊字符
    start_time.Format = 'MM_dd_yyyy_HH_mm_ss';
    formatted_start_time = string(start_time);
    
    end_time.Format = 'MM_dd_yyyy_HH_mm_ss';
    formatted_end_time = string(end_time);
    
    
    formatted_start_time = strrep(strrep(strrep(formatted_start_time, ':', '-'), ' ', '_'), '/', '-');
    formatted_end_time = strrep(strrep(strrep(formatted_end_time, ':', '-'), ' ', '_'), '/', '-');
    
    filename = "v2_Modal_updating_"+formatted_start_time+"_"+formatted_end_time+".mat";
    
    if exist(filename, 'file') == 2
        Modal_updating_file = load(filename);
    else
        error("The file "+filename+" does not exist. Please run the modal updating script first. The script is named as Operational_Modal_Analysis.m")
    end
    
    Modal_updating_data = Modal_updating_file.table_fre;
    
    for k1 = 1:size(Modal_updating_data,1)
        idx_ModalFEM_temp = Modal_updating_data.idx_ModalFEM(k1);
        replece_idx = find(abs(Freq - idx_ModalFEM_temp)<0.001);
        if ~isempty(replece_idx)
            Freq(replece_idx) = Modal_updating_data.frequency(k1);
            % zeta(replece_idx) = Modal_updating_data.damping_ratio(k1);
            KK_eq(replece_idx,replece_idx) = MM_eq(replece_idx,replece_idx) * (2 * pi * Freq(replece_idx))^2;
            
            
            % Freq(VIV_mode_seq) = 0.328240;
            % % Freq = 0.326;
            % KK_eq(VIV_mode_seq) = MM_eq(VIV_mode_seq) * (2 * pi * Freq(VIV_mode_seq))^2;
            if showtext
                disp("The "+num2str(replece_idx)+"th mode (original frequency: "+num2str(Result.Freq(replece_idx))+"Hz) is replaced by the "+num2str(k1)+"th mode in the modal updating process. The frequency is set as "+num2str(Modal_updating_data.frequency(k1))+"Hz and the damping ratio is set as "+num2str(Modal_updating_data.damping_ratio(k1)*100)+"%");
            end
        end
    end
end

%% temporary revise for the higher mode
% Freq(16) = 0.656;
% KK_eq(16,16) = MM_eq(16,16) * (2 * pi * Freq(16))^2;
% 
% Freq(23) = 0.985;
% KK_eq(23,23) = MM_eq(23,23) * (2 * pi * Freq(23))^2;

%%




mode_vec = Result.mode_vec;
nodeondeck = Result.nodeondeck;
Mapping_data = Result.Mapping;


if showtext
    disp("Damping ratio of the structure is set as "+num2str(zeta));
end
% TODO: test begin
% Freq(4)=0.15
% TODO: test end
omega = diag(2 * pi * Freq);
CC_eq = 2 .* MM_eq .* omega .* zeta;

phi = mode_vec; %模态向量 每一列是一个模态
% TODO: test begin
% phi(:,4)=phi(:,4)*0.00001
% TODO: test end
C = CC_eq; K = KK_eq; M = MM_eq;

Gamma = C; % 对应文献中表述的符号
omega2 = K;

%% 3 传感器布置
% accelerometer location
% loc_acc= [578+1650/4*3;578+1650/2;578+1650/4];
% loc_acc = [1403];

% loc_acc = [1403];
loc_acc = [990.5; 1403; 1815.5];
loc_vel = [];
loc_dis = [];

timeDifferences = diff(Acc_Data.mergedData.Time); % Returns a duration array
dt = seconds(timeDifferences(1)); % Converts the first duration value to seconds and assigns to dt
fs = 1 / dt;
acc_names = ["Main span 1/4", "Main span 1/2", "Main span 3/4"];
% yn(1, :) = Acc_Data.mergedData.AC2_1 / 1000 * 9.8;
% yn(2, :) = Acc_Data.mergedData.AC2_3 / 1000 * 9.8;
% yn(3, :) = Acc_Data.mergedData.AC3_1 / 1000 * 9.8;
% yn(4, :) = Acc_Data.mergedData.AC3_3 / 1000 * 9.8;
% yn(5, :) = Acc_Data.mergedData.AC4_1 / 1000 * 9.8;
% yn(6, :) = Acc_Data.mergedData.AC4_3 / 1000 * 9.8;
% 

acc_result=Acc_Data.mergedData;
AC2_1 = acc_result.AC2_1';
AC2_3 = acc_result.AC2_3';
AC3_1 = acc_result.AC3_1';
AC3_3 = acc_result.AC3_3';
AC4_1 = acc_result.AC4_1';
AC4_3 = acc_result.AC4_3';

%% selection of the sensors' data
sensor_sel= string(vivTable.sensor_selection(input.VIV_sel));
sensor_sel = strsplit(sensor_sel, ';'); % 以分号为分隔符分割字符串
sensor_sel = str2double(sensor_sel); % 将字符串数组转换为double数组

% 检查1和2是否在数组中
contains1 = ismember(1, sensor_sel);
contains2 = ismember(2, sensor_sel);
AC2 = sel_sensor(AC2_1,AC2_3,contains1,contains2,showtext);

contains1 = ismember(3, sensor_sel);
contains2 = ismember(4, sensor_sel);
AC3 = sel_sensor(AC3_1,AC3_3,contains1,contains2,showtext);

contains1 = ismember(5, sensor_sel);
contains2 = ismember(6, sensor_sel);
AC4 = sel_sensor(AC4_1,AC4_3,contains1,contains2,showtext);

yn = [AC2;AC2;AC3;AC3;AC4;AC4]/ 1000 * 9.8;

delay2= int32(vivTable.delay2(input.VIV_sel));
delay3= int32(vivTable.delay3(input.VIV_sel));

yn(3,:)=circshift(yn(3,:),-delay2/dt);
yn(4,:)=circshift(yn(4,:),-delay2/dt);


yn(5,:)=circshift(yn(5,:),-delay3/dt);
yn(6,:)=circshift(yn(6,:),-delay3/dt);


% yn(1, :) = (Acc_Data.mergedData.AC2_1 / 1000 * 9.8+Acc_Data.mergedData.AC2_3 / 1000 * 9.8)/2;
% yn(2, :) = (Acc_Data.mergedData.AC2_1 / 1000 * 9.8+Acc_Data.mergedData.AC2_3 / 1000 * 9.8)/2;
% yn(3, :) = (Acc_Data.mergedData.AC3_1 / 1000 * 9.8+Acc_Data.mergedData.AC3_3 / 1000 * 9.8)/2;
% yn(4, :) = (Acc_Data.mergedData.AC3_1 / 1000 * 9.8+Acc_Data.mergedData.AC3_3 / 1000 * 9.8)/2;
% % yn(5, :) = (Acc_Data.mergedData.AC4_1 / 1000 * 9.8+Acc_Data.mergedData.AC4_3 / 1000 * 9.8)/2;
% % yn(6, :) = (Acc_Data.mergedData.AC4_1 / 1000 * 9.8+Acc_Data.mergedData.AC4_3 / 1000 * 9.8)/2;
% yn(5, :) = Acc_Data.mergedData.AC4_3 / 1000 * 9.8;
% yn(6, :) = Acc_Data.mergedData.AC4_3 / 1000 * 9.8;
% figure
% plot(yn(1, :))
% hold on 
% plot(yn(2, :))
% 
% plot(yn(3, :))
% plot(yn(4, :))
% plot(yn(5, :))
% plot(yn(6, :))
% legend("1","2","3","4","5","6")
% yn(1, :) = Acc_Data.mergedData.AC3_1 / 1000 * 9.8;
% yn(2, :) = Acc_Data.mergedData.AC3_3 / 1000 * 9.8;

% yn(1, :) = (Acc_Data.mergedData.AC3_1 / 1000 * 9.8+Acc_Data.mergedData.AC3_3 / 1000 * 9.8)/2;
% yn(2, :) = (Acc_Data.mergedData.AC3_1 / 1000 * 9.8+Acc_Data.mergedData.AC3_3 / 1000 * 9.8)/2;

if shouldFilterYn == true
    yn_temp = yn;
    clear yn
    lowpassfreq = 0.36;
    yn(1,:) = lowpass(yn_temp(1,:),lowpassfreq,fs,Steepness=0.99);
    yn(2,:) = lowpass(yn_temp(2,:),lowpassfreq,fs,Steepness=0.99);
    yn(3,:) = lowpass(yn_temp(3,:),lowpassfreq,fs,Steepness=0.99);
    yn(4,:) = lowpass(yn_temp(4,:),lowpassfreq,fs,Steepness=0.99);
    yn(5,:) = lowpass(yn_temp(5,:),lowpassfreq,fs,Steepness=0.99);
    yn(6,:) = lowpass(yn_temp(6,:),lowpassfreq,fs,Steepness=0.99);
    [f, magnitude] =fft_transform(fs,yn_temp(1,:));
    figure
    semilogy(f,magnitude)
    xlim([0,1])
    [f, magnitude] =fft_transform(fs,yn(1,:));
    figure
    semilogy(f,magnitude)
    xlim([0,1])
    figure
    plot(uniqueTimestamps,yn(1,:))
    hold on
    plot(uniqueTimestamps,yn(2,:))
end


[S_a, S_v, S_d, n_sensors] = sensor_selection(loc_acc, loc_vel, loc_dis, node_loc, phi, nodeondeck, Mapping_data);
% establish continuous time matrices
[A_c, B_c_buff, G_c, J_c_buff] = ssmod_c_mode(nmodes, omega2, Gamma, phi, S_a, S_v, S_d);

[B_c_VIV, J_c_VIV] = ssmod_c_mode_VIV(nmodes,  phi, S_a, VIV_mode_seq);

maxvalue = max(max(abs(yn)));

% establish discrete time matrices

[A_d, B_d_buff, G_d, J_d_buff, ~] = ssmod_c2d(A_c, B_c_buff, G_c, J_c_buff, dt);
[~, B_d_VIV, ~, J_d_VIV, ~] = ssmod_c2d(A_c, B_c_VIV, G_c, J_c_VIV, dt);

%% 4 反算模态力
C_buff = sigma_buff^2*eye(nmodes);



Q = Q_value * eye(ns);
R = R_value * eye(n_sensors);
% TODO: test begin
% R(1,1)=10;
% R(2,2)=10;
% R(5,5)=10;
% R(6,6)=10;
% TODO: test end 
x0 = zeros(ns, 1);

Q_xd = Q;
np_m = nmodes;
% B_c_m = [zeros(nmodes, np_m); ...
%              eye(np_m, np_m)];
B_c_m = B_c_VIV;
J_c_m = J_c_VIV;
B_d_m = B_d_VIV;
J_d_m = J_d_VIV;
% TODO: Here the B_c_m, J_c_m, B_d_m, J_d_m are equal to B_c, J_c, B_d, J_d, which is can be revised in the future

% set the kernal parameters for the latent force model
% quasiperiodic kernel
for k1 = 1:np_m
    if k1 == VIV_mode_seq
        lambda_quasi_periodic(k1) = lambda_VIV;
        sigma_ps_quasi_periodic(k1) = sigma_p_VIV;
        omega_0_quasi_periodic(k1) = 2 * pi * Freq(VIV_mode_seq)*omega_0_variation_VIV;
    else
        lambda_quasi_periodic(k1) = 0;
        sigma_ps_quasi_periodic(k1) = 0;
        omega_0_quasi_periodic(k1) = 0;
    end
end

% % matern kernel
% for k1 = 1:np_m
%     if k1 == VIV_mode_seq
%         lambda_matern(k1) = 0;
%         sigma_ps_matern(k1) = 0;
%     else
%         lambda_matern(k1) = lambda_matern_noVIV;
%         sigma_ps_matern(k1) = sigma_p_matern_noVIV;
%     end
% end


% lambdas_m = [lambda] * ones(1, np_m);
% sigma_ps_m = [sigma_p] * ones(1, np_m);
% omega_0 = 2 * pi * Freq*omega_0_variation* ones(1, np_m);

nVIV = length(VIV_mode_seq);
[F_c_m, L_c_m, H_c_m, sigma_w_m12] = ssmod_quasiperiod_coninue(lambda_quasi_periodic, sigma_ps_quasi_periodic, omega_0_quasi_periodic, np_m,VIV_mode_seq);

% [E_c_m,K_c_m,T_c_m,sigma_z_m12] = ssmod_matern_coninue(lambda_matern, sigma_ps_matern, np_m,VIV_mode_seq);

% [~, ~, ~, ~, Fad_m, ~, Gad_m, ~, Qad_m]=ssmod_lfm_aug_matern_and_quasiperiod(A_c, B_c_m, G_c, J_c_m, F_c_m, H_c_m, L_c_m,E_c_m,T_c_m,K_c_m, Q_xd, sigma_w_m12,sigma_z_m12, dt);
[~, ~, ~, ~, Fad_m, ~, Gad_m, ~, Qad_m] = ssmod_lfm_aug(A_c, B_c_m, G_c, J_c_m, F_c_m, H_c_m, L_c_m, Q_xd, sigma_w_m12, dt);
A_a_m = Fad_m;
G_a_m = Gad_m;
Q_a_m = Qad_m;
% Q_buff_c= B_c_buff * C_buff * B_c_buff';
% Q_buff_d = Q_buff_c*dt;
Q_buff_d = B_d_buff * C_buff * B_d_buff'/dt;
Q_a_m(1:ns,1:ns)=Q_a_m(1:ns,1:ns)+Q_buff_d;
R_a_m = J_d_buff*C_buff*J_d_buff'/dt+ R;
S_a_m = [B_d_buff*C_buff*J_d_buff'/dt;zeros(size(L_c_m,1),size(B_d_buff*C_buff*J_d_buff',2))];
% TODO: here might be wrong
yn_a = yn;
N = length(Acc_Data.mergedData.Time);
NN = N;
xa_history = zeros(ns + np_m * (2), NN);
pa_history = zeros(ns + np_m * (2), NN);

x_ak = zeros(ns + nVIV * (2), 1);
P_ak = 10 ^ (1) * eye(ns + nVIV * (2));
p_det= zeros(size(A_a_m,1),size(yn_a,2));
B_a_m = zeros(size(A_a_m));
J_a_m = zeros(size(G_a_m));

% [x_k_k, x_k_kmin, P_k_k, P_k_kmin, result] = KalmanFilterNoInput(A_a_m, G_a_m, Q_a_m, R_a_m, yn_a, x_ak, P_ak, 'debugstate', true,'showtext',showtext);
% [x_k_k,x_k_kmin,P_k_k,P_k_kmin,K_k_ss]=KalmanFilter(A_a_m, G_a_m, Q_a_m, R_a_m,S_a_m, yn_a, x_ak, P_ak,'steadystate','no');
[x_k_k,x_k_kmin,P_k_k,P_k_kmin,K_k_ss,M_k_ss,Omega_k_ss,result]=KF(A_a_m,B_a_m, G_a_m,J_a_m, Q_a_m, R_a_m,S_a_m, yn_a,p_det, x_ak, P_ak,'steadystate',true,'showtext',showtext,'debugstate',true);
% [x_k_k, P_k_k] = RTSFixedInterval(A_a_m, x_k_k, x_k_kmin, P_k_k, P_k_kmin);
% [x_k_k,P_k_k]=RTSSmoother(A_a_m,x_k_k,x_k_kmin,P_k_k,P_k_kmin,'showtext',false);
[x_k_k,P_k_k]=RTSS(A_a_m-S_a_m/R_a_m*G_a_m,x_k_k,x_k_kmin,P_k_k,P_k_kmin,'showtext',showtext);
xa_history = x_k_k;
pa_history = P_k_k;
x_filt_original = xa_history(1:ns, :);
H_d_m = H_c_m;
p_filt_m = H_d_m * xa_history(ns + 1:end, :);


Pp_filt_m = H_d_m * pa_history(ns + 1:end,ns + 1:end)*H_d_m';
% [f, magnitude] = fft_transform(fs, x_k_k(1, :));


%% test if discrete the equation before augment it has the same results.
if 0
    Ac = A_c;
    Bc = B_c_m;
    Hc = H_c_m;
    Fc = F_c_m;

    Fa_c = [Ac,Bc*Hc;zeros(size(Fc,1),size(Ac,2)),Fc];

    [Fa_d,~,~,~]=ssmod_c2d(Fa_c,[],[],[],dt);
    
    [Ad,Bd,~,~,~]=ssmod_c2d(Ac,Bc,[],[],dt);
    Hd=Hc;
    [Fd,~,~,~]=ssmod_c2d(Fc,[],[],[],dt);

    Fa_d_2=[Ad,Bd*Hd;zeros(size(Fd,1),size(Ad,2)),Fd];
    test=Fa_d_2-Fa_d;

end

%% 5 fft and bandpass filter for the estimated modal force
for k1 = 1:length(VIV_mode_seq)
    if shouldFilterp_filt_m == true
        [f, magnitude] = fft_transform(fs, x_k_k(k1,:));
        [~, idx] = max(magnitude);
        f_keep_temp = [f(idx) * 0.8, f(idx) * 1.2];
        % p_filt_m = fft_filter(fs, p_filt_m, f_keep_temp);
        p_filt_m = bandpass(p_filt_m', f_keep_temp, fs)';
    end
    [f_p_filt_m(k1, :), magnitude_filt_m(k1, :)] = fft_transform(fs, p_filt_m(k1, :));
end

%% 6 virtual sensoring
% loc_acc_v = [990.5; 1403; 1815.5];
loc_acc_v = loc_acc;
loc_vel_v = loc_acc;
loc_dis_v = loc_acc;
% loc_acc_v = [1403];
% loc_acc_v = [578+1650/4*3;578+1650/2;578+1650/4];
% loc_vel_v = [1403];
% loc_dis_v = [1403];
% loc_vel_v = [];
% loc_dis_v = [];
[S_a_v, S_v_v, S_d_v, n_sensors_v] = sensor_selection(loc_acc_v, loc_vel_v, loc_dis_v, node_loc, phi, nodeondeck, Mapping_data);

[B_c_v_VIV, J_c_v_VIV] = ssmod_c_mode_VIV(nmodes,  phi, S_a_v, VIV_mode_seq);

G_c_v = [S_d_v * phi - S_a_v * phi * omega2, S_v_v * phi - S_a_v * phi * Gamma];
% J_c_v = [S_a_v * phi];

% h_hat = G_c_v * x_filt_original + J_c_v * p_filt_m;
h_hat = G_c_v * x_filt_original + J_c_v_VIV * p_filt_m;

%calculate covariance
G_a_v=[G_c_v,J_c_v_VIV*H_c_m ];
h_hat_covariance = G_a_v*P_k_k*G_a_v';

%% 7 重构涡振响应
p_reconstruct = p_filt_m;
% p_reconstruct([1 2 3 ],:)=0;
[~, yn_reconstruct, ~] = CalResponse(A_d, B_d_VIV, G_d, J_d_VIV, p_reconstruct, 0, 0, N, x0, ns, n_sensors);

% recalculate logek (this part of code will never be used due to it is against the methodology)
% nt = size(yn_reconstruct,2);
% logek=0;
% for k1 = 1:nt
%     ek = yn(:,k1)-yn_reconstruct(:,k1);
%     logek_i = -0.5*ek.'*result.invSk*ek;
%     % logek_i = -0.5*ek.'*ek;
%     logek = logek + logek_i;
% end

% fft
% [f_origin, magnitude_origin] = fft_transform(1 / dt, yn(3, :));
% [f_re, magnitude_re] = fft_transform(1 / dt, yn_reconstruct(3, :));

%% 8 marginal likelihood
%     logL = result.logL;
logSk = result.logSk;
logek = result.logek;
invSk = result.invSk;
logL = logSk +logek;

t = Acc_Data.mergedData.Time;

%% output
result_Main.logL = logL;
result_Main.logSk = logSk;
result_Main.logek = logek;
result_Main.invSk = invSk;
result_Main.t = t;
result_Main.p_filt_m = p_filt_m;
result_Main.x_k_k = x_k_k;
result_Main.nmodes = nmodes;
result_Main.Freq = Freq;
result_Main.fs = fs;
result_Main.mode_deck = mode_deck;
result_Main.yn=yn;
result_Main.h_hat = h_hat;
result_Main.yn_reconstruct = yn_reconstruct;
result_Main.x_filt_original=x_filt_original;
result_Main.MM = MM_eq;
result_Main.CC = CC_eq;
result_Main.KK = KK_eq;
result_Main.A_d = A_d;
% result_Main.B_d = B_d;
result_Main.G_d = G_d;
% result_Main.J_d = J_d;
result_Main.N = N;
result_Main.x0 = x0;
result_Main.ns = ns;
result_Main.n_sensors = n_sensors;
result_Main.S_a = S_a;
result_Main.phi = phi;
result_Main.loc_acc = loc_acc;
result_Main.zeta = zeta;
result_Main.P_k_k = P_k_k;
result_Main.Pp_filt_m = Pp_filt_m;
result_Main.h_hat_covariance=h_hat_covariance;

if fig_bool == ON
    
    for k1 = 1:length(loc_acc)
        [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
        
        % subplot(1, nmodes, k1)
        plot(t, yn(2 * k1 - 1, :), 'Color', 'r')
        hold on
        plot(t, yn(2 * k1, :), 'Color', 'b')
        
        xlabel('time (s)')
        ylabel('acc (m/s^2)')
        set(gca, 'FontSize', 12)
        legend('left','right', 'Location', 'northwest')
        title([num2str(loc_acc(k1))]);
        ylim([-maxvalue, maxvalue])
    end
    
    for k1 = 1:length(loc_acc)
        [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
        
        % subplot(1, nmodes, k1)
        plot(t, h_hat(2 * k1 - 1, :), 'Color', 'r')
        hold on
        plot(t, h_hat(2 * k1, :), 'Color', 'b')
        
        xlabel('time (s)')
        ylabel('acc (m/s^2)')
        set(gca, 'FontSize', 12)
        legend('left','right', 'Location', 'northwest')
        title([num2str(loc_acc(k1))] + "filter");
        ylim([-maxvalue, maxvalue])
    end
    
    for k1 = 1:length(loc_acc)
        [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
        
        % subplot(1, nmodes, k1)
        plot(t, yn_reconstruct(2 * k1 - 1, :), 'Color', 'r')
        hold on
        plot(t, yn_reconstruct(2 * k1, :), 'Color', 'b')
        
        xlabel('time (s)')
        ylabel('acc (m/s^2)')
        set(gca, 'FontSize', 12)
        legend('left','right', 'Location', 'northwest')
        title([num2str(loc_acc(k1))] + "recalculate");
        ylim([-maxvalue, maxvalue])
    end
    
    for k1 = 1:nmodes
        [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
        
        % subplot(1, nmodes, k1)
        plot(t, p_filt_m(k1, :), 'Color', 'r')
        hold on
        %         plot(t, p_m_real(k1, :), 'Color', 'b', 'LineStyle', '--')
        xlabel('time (s)')
        ylabel('Modal force ')
        set(gca, 'FontSize', 12)
        legend('filtered modal force', 'Location', 'northwest')
        title(['mode ', num2str(modesel(k1))]);
        % ylim([-1e3, 1e3])
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
        xlim([0, 0.5])
        % ylim([0, 50])
    end
    
    % set(hFigure, 'name', 'filtered modal force frequency', 'Numbertitle', 'off');
    % % figureIdx=0;
    % [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
    % plot(f_re, magnitude_re)
    % hold on
    % plot(f_origin, magnitude_origin)
    % legend("cal", "measure")
    % xlim([0, 0.5])
    %
    % [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
    % plot(f_re, log(magnitude_re))
    % hold on
    % plot(f_origin, log(magnitude_origin))
    % legend("cal", "measure")
    % xlim([0, 0.5])
    
    % [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
    % plot(t2, ifq1)
    % title("inst frequency")
    
    % [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
    % plot(f, magnitude)
    % title("frequency of estimated vibration")
    
    % instfreq(p, fd, td);
    
    
    
end




end

function data = sel_sensor(data1,data2,contains1,contains2,showtext)
% 基于条件执行不同的操作
if contains1 && ~contains2
    % 数组仅包含1的操作
    if showtext
        disp('Array contains only contains 1.');
    end
    data=data1;
elseif ~contains1 && contains2
    % 数组仅包含2的操作
    if showtext
    disp('Array contains only contains 2.');
    end
    data=data2;
elseif contains1 && contains2
    % 数组同时包含1和2的操作
    if showtext
        disp('Array contains both contains 1 and contains 2.');
    end
    data = (data1+data2)/2;
else
    % 数组既不包含1也不包含2的操作
    error('Array contains neither contains 1 nor contains 2.');
end
end

