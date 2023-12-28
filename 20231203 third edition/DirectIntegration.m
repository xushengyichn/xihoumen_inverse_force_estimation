clc; clear; close all;
run('CommonCommand.m');

n = 4;
[result] = viv2013(n, OFF);
startDate_global = result.startDate;
endDate_global = result.endDate;
input_data.start_time = startDate_global;
input_data.end_time = endDate_global;

start_time = input_data.start_time;
end_time = input_data.end_time;

acc_dir = input_data.acc_dir;

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

% loc_acc = [989; 1403; 1783];
loc_acc = [990.5; 1403; 1815.5];

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
% fw(abs(freq) < 0.03) = 0;

ft = ifft(ifftshift(fw));
ft_real = real(ft);
ft_imag = imag(ft);



total_plots = 2; % 或任何你需要的子图数量
current_plot = 1;
num_figs_in_row = [];
newfigure = false;
firstfigure=true;

data = importdata("temp.mat");
t= data.t;
F_filter = data.F_filter;


create_subplot(@plot, total_plots, current_plot, {t,F_filter}, 'num_figs_in_row', num_figs_in_row,'newfigure',newfigure);
hold on
create_subplot(@plot, total_plots, current_plot, {t, ft}, 'num_figs_in_row', num_figs_in_row,'newfigure',newfigure,'firstfigure',firstfigure);
legend("filter","direct ingetration")
title("ft");
current_plot = current_plot + 1;


save DirectIntegration.mat ft


%% frequency domian analysis
fs = 50;
[f1, magnitude_ft] =fft_transform(fs,ft);
[f2, magnitude_F_filter] =fft_transform(fs,F_filter);

create_subplot(@plot, total_plots, current_plot, {f1, magnitude_ft}, 'num_figs_in_row', num_figs_in_row,'newfigure',newfigure);
title("magnitude_ft");
hold on
create_subplot(@plot, total_plots, current_plot, {f2, magnitude_F_filter}, 'num_figs_in_row', num_figs_in_row,'newfigure',newfigure);
xlim([0,0.5])

figure
current_plot =  1;
[f, magnitude] =fft_transform(fs,yn(3, :));
create_subplot(@plot, total_plots, current_plot, {t, yn(3, :)}, 'num_figs_in_row', num_figs_in_row,'newfigure',true);
title("Time domain");
xlabel("Time")
ylabel("Acceleration (m/s^2)")
print -clipboard -dbitmap

current_plot = current_plot + 1;
create_subplot(@plot, total_plots, current_plot, {f, magnitude}, 'num_figs_in_row', num_figs_in_row,'newfigure',true);
xlim([0,0.5])
set(gca,'YScale','log')
title("Frequency domain");
xlabel("Frequency (Hz)")
ylabel("Amplitude")
current_plot = current_plot + 1;
