clc; clear; close all;
run('CommonCommand.m');

n = 4;
[result] = viv2013(n, OFF);
startDate_global = result.startDate;
endDate_global = result.endDate + hours(1);
input.start_time = startDate_global;
input.end_time = endDate_global;

start_time = input.start_time;
end_time = input.end_time;

acc_dir = input.acc_dir;

[Acc_Data] = read_acceleration_data(start_time, end_time, acc_dir);

f_keep = 0.33 * [0.9, 1.1];
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


yn(1, :) = (Acc_Data.mergedData.AC2_1 / 1000 * 9.8+Acc_Data.mergedData.AC2_3 / 1000 * 9.8)/2;
yn(2, :) = (Acc_Data.mergedData.AC2_1 / 1000 * 9.8+Acc_Data.mergedData.AC2_3 / 1000 * 9.8)/2;
yn(3, :) = (Acc_Data.mergedData.AC3_1 / 1000 * 9.8+Acc_Data.mergedData.AC3_3 / 1000 * 9.8)/2;
yn(4, :) = (Acc_Data.mergedData.AC3_1 / 1000 * 9.8+Acc_Data.mergedData.AC3_3 / 1000 * 9.8)/2;
yn(5, :) = (Acc_Data.mergedData.AC4_1 / 1000 * 9.8+Acc_Data.mergedData.AC4_3 / 1000 * 9.8)/2;
yn(6, :) = (Acc_Data.mergedData.AC4_1 / 1000 * 9.8+Acc_Data.mergedData.AC4_3 / 1000 * 9.8)/2;


modesel = 23;
nmodes = length(modesel); ns = nmodes * 2;
Result = ImportMK(nmodes, 'KMatrix.matrix', 'MMatrix.matrix', 'nodeondeck.txt', 'KMatrix.mapping', 'nodegap.txt', 'modesel', modesel, 'showtext', false);
node_loc = Result.node_loc;
mode_vec = Result.mode_vec;
nodeondeck = Result.nodeondeck;
Mapping_data = Result.Mapping;
phi = mode_vec; %模态向量 每一列是一个模态
nodegap = Result.nodegap;

xi = 0.3/100;
Omega = 2 * pi * Result.Freq;
Fs = 50;
MM = 1;
KK = MM * Omega ^ 2;
CC = 2 * MM * Omega * xi;


N = length(yn(1, :));
t = (0:N - 1) / Fs;

loc_acc = [989; 1403; 1783];

loc_vel = [];
loc_dis = [];
[S_a, S_v, S_d, n_sensors] = sensor_selection(loc_acc, loc_vel, loc_dis, node_loc, phi, nodeondeck, Mapping_data);

nodeshape = FindModeShapewithLocation(loc_acc, node_loc, nodeondeck, Mapping_data, nodegap, mode_vec);


%% fft
N = length(yn(1, :)); % 信号长度

yw = fftshift(fft(yn, [], 2), 2);

freq = linspace(-0.5, 0.5, length(t)) * Fs;
omega = 2 * pi * freq;

Hw = 1 ./ (-omega .^ 2 + 2 * 1i * xi * Omega * omega + Omega ^ 2);

fw = pinv(S_a * phi) * yw .* 1 ./ (-omega .^ 2 .* Hw);

fw(abs(freq) < 0.3) = 0;

ft = ifft(ifftshift(fw));
ft_real = real(ft);
ft_imag = imag(ft);



total_plots = 10; % 或任何你需要的子图数量
current_plot = 1;
num_figs_in_row = [];


create_subplot(@plot, total_plots, current_plot, {t, ft}, 'num_figs_in_row', num_figs_in_row);
title("ft");
current_plot = current_plot + 1;



