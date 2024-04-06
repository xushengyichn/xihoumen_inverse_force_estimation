clc;clear;close all
run('CommonCommand.m');
VIV_sel=3;
input.VIV_sel =VIV_sel;
modeall = [2, 3, 5, 6, 7, 13, 20, 22, 27, 33];
moderemove = [4,9,10];
modeall(moderemove)=[];
modesel = modeall;
nmodes = length(modesel);
showtext = true;
%% Load data
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
start_time = vivTable.startDate(VIV_sel);
end_time = vivTable.endDate(VIV_sel);
startDate_global = start_time;
endDate_global = end_time;
input_data.start_time = startDate_global;
input_data.end_time = endDate_global;
input_data.VIV_sel = VIV_sel;

Result = ImportMK(nmodes, 'KMatrix.matrix', 'MMatrix.matrix', 'nodeondeck.txt', 'KMatrix.mapping', 'nodegap.txt', 'modesel', modesel,'showtext',showtext);
mode_deck = Result.mode_deck; mode_deck_re = Result.mode_deck_re; node_loc = Result.node_loc;
Freq = Result.Freq;
MM_eq = Result.MM_eq; KK_eq = Result.KK_eq;
zeta = importdata("zeta_update0406.mat");

input.modelupdate=true;
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

if showtext
    disp("Damping ratio of the structure is set as "+num2str(zeta));
end

mode_vec = Result.mode_vec;
nodeondeck = Result.nodeondeck;
Mapping_data = Result.Mapping;

omega = diag(2 * pi * Freq);
CC_eq = 2 .* MM_eq .* omega .* zeta;

phi = mode_vec; %模态向量 每一列是一个模态
C = CC_eq; K = KK_eq; M = MM_eq;
Gamma = C; % 对应文献中表述的符号
omega2 = K;

%% 3 传感器布置
loc_acc = [];
loc_vel = [];
loc_dis = [990.5; 1403; 1815.5];

dt = 0.02;
fs = 1 / dt;

[S_a, S_v, S_d, n_sensors] = sensor_selection(loc_acc, loc_vel, loc_dis, node_loc, phi, nodeondeck, Mapping_data);
% establish continuous time matrices
[A_c, B_c_buff, G_c, J_c_buff] = ssmod_c_mode(nmodes, omega2, Gamma, phi, S_a, S_v, S_d);

[A_d, B_d_buff, G_d, J_d_buff, ~] = ssmod_c2d(A_c, B_c_buff, G_c, J_c_buff, dt);


%% simulation
sigma_buff = 10 ^ (0.5411);
t = 0:dt:1000;
p_buff = randn(nmodes, length(t))*sigma_buff;

% verify the correctness of the simulation
std1 = std(p_buff(1, :));
std2 = std(p_buff(2, :));
std3 = std(p_buff(3, :));
std4 = std(p_buff(4, :));
std5 = std(p_buff(5, :));
std6 = std(p_buff(6, :));
std7 = std(p_buff(7, :));

disp("The standard deviation of the random excitation is as follows:")
disp(std1)
disp(std2)
disp(std3)
disp(std4)
disp(std5)
disp(std6)
disp(std7)

x = zeros(nmodes*2, 1);

for k1 = 1:length(t)
    p = p_buff(:, k1);
    x_new = A_d * x + B_d_buff * p;
    y = G_d * x + J_d_buff * p;
    x=x_new;
    y_buff(:, k1) = y;
end



[f1, magnitude1] = fft_transform(fs,y_buff(1,:));
[f2, magnitude2] = fft_transform(fs,y_buff(3,:));
[f3, magnitude3] = fft_transform(fs,y_buff(5,:));



realdatapath="D:\Users\xushe\Documents\Paper\VIV_force_estimation\658209c08d0dc70f30a78d4e\python\plotdata\Modaldisplacement.mat";
% 从论文绘图文件中找到该绘图数据，基于不同电脑路径有所不同
real_data = importdata(realdatapath);
h_hat = real_data.h_hat;
dis1=h_hat(13,:);
dis2=h_hat(15,:);
dis3=h_hat(17,:);

[f1_measure, magnitude1_measure] = fft_transform(fs,dis1);
[f2_measure, magnitude2_measure] = fft_transform(fs,dis2);
[f3_measure, magnitude3_measure] = fft_transform(fs,dis3);

figure
semilogy(f1, magnitude1)
hold on
semilogy(f2, magnitude2)
semilogy(f3, magnitude3)
xlim([0,0.5])
legend("simulation1","simulation2","simulation3")

figure
semilogy(f1_measure, magnitude1_measure)
hold on
semilogy(f2_measure, magnitude2_measure)
semilogy(f3_measure, magnitude3_measure)
xlim([0,0.5])
legend("m1","m2","m3")

figure
semilogy(f1, magnitude1)
hold on
semilogy(f1_measure, magnitude1_measure)
xlim([0,0.5])
legend("s1","m1")


figure
semilogy(f2, magnitude2)
hold on
semilogy(f2_measure, magnitude2_measure)
xlim([0,0.5])
legend("s1","m1")


figure
semilogy(f3, magnitude3)
hold on
semilogy(f3_measure, magnitude3_measure)
xlim([0,0.5])
legend("s1","m1")


figure
semilogy(f1_measure, magnitude1_measure)
hold on
semilogy(f1, magnitude1,LineWidth=2)
xlim([0,0.5])
legend("m1","s1")
title("1/4 span")


figure
semilogy(f2_measure, magnitude2_measure)
hold on
semilogy(f2, magnitude2,LineWidth=2)
xlim([0,0.5])
legend("m2","s2")
title("1/2 span")

figure
semilogy(f3_measure, magnitude3_measure)
hold on
semilogy(f3, magnitude3,LineWidth=2)
xlim([0,0.5])
legend("m3","s3")
title("3/4 span")

tilefigs