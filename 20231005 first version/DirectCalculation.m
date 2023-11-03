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

% yn(1, :) = (Acc_Data.mergedData.AC3_1 / 1000 * 9.8+Acc_Data.mergedData.AC3_3 / 1000 * 9.8)/2;
% yn(2, :) = (Acc_Data.mergedData.AC3_1 / 1000 * 9.8+Acc_Data.mergedData.AC3_3 / 1000 * 9.8)/2;
%

yn(1, :) = Acc_Data.mergedData.AC2_1 / 1000 * 9.8;
yn(2, :) = Acc_Data.mergedData.AC2_3 / 1000 * 9.8;
yn(3, :) = Acc_Data.mergedData.AC3_1 / 1000 * 9.8;
yn(4, :) = Acc_Data.mergedData.AC3_3 / 1000 * 9.8;
yn(5, :) = Acc_Data.mergedData.AC4_1 / 1000 * 9.8;
yn(6, :) = Acc_Data.mergedData.AC4_3 / 1000 * 9.8;

%%

modesel = 23;
nmodes = length(modesel); ns = nmodes * 2;
Result = ImportMK(nmodes, 'KMatrix.matrix', 'MMatrix.matrix', 'nodeondeck.txt', 'KMatrix.mapping', 'nodegap.txt', 'modesel', modesel, 'showtext', false);
node_loc = Result.node_loc;
mode_vec = Result.mode_vec;
nodeondeck = Result.nodeondeck;
Mapping_data = Result.Mapping;
phi = mode_vec; %模态向量 每一列是一个模态
nodegap = Result.nodegap;

%%
xi = 0.3/100;
Omega = 2 * pi * 0.328;
Fs = 50;
MM = 1;
KK = MM * Omega ^ 2;
CC = 2 * MM * Omega * xi;

N = length(yn(1, :));
t = (0:N - 1) / Fs;

% loc_acc = [1403];
loc_acc = [990.5; 1403; 1815.5];
loc_vel = [];
loc_dis = [];
[S_a, S_v, S_d, n_sensors] = sensor_selection(loc_acc, loc_vel, loc_dis, node_loc, phi, nodeondeck, Mapping_data);

nodeshape = FindModeShapewithLocation(loc_acc, node_loc, nodeondeck, Mapping_data, nodegap, mode_vec);

%% Simulated signal
% omega_sim = Omega;
% f_sim =50* sin(omega_sim.*t);
% [u_sim udot_sim u2dot_sim] = NewmarkInt(t, MM, CC, KK, f_sim, 1/2, 1/4, 0, 0);
%
% yn(1, :) = u2dot_sim*nodeshape;
% yn(2, :) = u2dot_sim*nodeshape;


%% fft
N = length(yn(1, :)); % 信号长度
yw(1, :) = fftshift(fft(yn(1, :)));
yw(2, :) = fftshift(fft(yn(2, :)));
yw(3, :) = fftshift(fft(yn(3, :)));
yw(4, :) = fftshift(fft(yn(4, :)));
yw(5, :) = fftshift(fft(yn(5, :)));
yw(6, :) = fftshift(fft(yn(6, :)));

freq = linspace(-0.5, 0.5, length(t)) * Fs;
omega = 2 * pi * freq;

% yw = [yw1;yw2];

%%

Hw = 1 ./ (-omega .^ 2 + 2 * 1i * xi * Omega * omega + Omega ^ 2);

fw = pinv(S_a * phi) * yw .* 1 ./ (-omega .^ 2 .* Hw);

% for k1 = 1:length(omega)
%     fw(k1)= pinv(S_a*phi)*yw(:,k1).*1./(-omega(k1).^2.*Hw(k1));
% end

fw(abs(freq) < 0.3) = 0;

% figure
% plot(freq, fw)
% 
% figure
% plot(freq, pinv(S_a * phi) * yw)

% figure
% plot(freq, Hw)


ft = ifft(ifftshift(fw));
ft_real = real(ft);
ft_imag = imag(ft);

%% recalculate

[u_re, udot_re, u2dot_re] = NewmarkInt(t, MM, CC, KK, ft_real, 1/2, 1/4, 0, 0);

%% direct integral
if 0
    disp_dir = acc2dsip(yn(1, :), 50);
    F_direct = MM .* yn(1, :) / nodeshape + CC .* disp_dir.vel / nodeshape + KK .* disp_dir.disp / nodeshape;

    [u_direct udot_direct u2dot_direct] = NewmarkInt(t, MM, CC, KK, F_direct, 1/2, 1/4, 0, 0);

end

%% plot

% 定义总子图数量
total_plots = 7; % 或任何你需要的子图数量
current_plot = 1;
num_figs_in_row=5;
% 画第一个子图
create_subplot(@plot, total_plots, current_plot, {t, u2dot_re},'num_figs_in_row',num_figs_in_row);

legend("use ft_real to reclaculate the acceleration");
title("reconstructed displacement");
current_plot = current_plot + 1;

create_subplot(@plot, total_plots, current_plot, {freq, yw(1, :)},'num_figs_in_row',num_figs_in_row);
legend("yw(1, :)");
title("frequency domian of y");
current_plot = current_plot + 1;

create_subplot(@plot, total_plots, current_plot, {freq, fw},'num_figs_in_row',num_figs_in_row);
legend("fw");
title("frequency domian of f");
current_plot = current_plot + 1;

create_subplot(@plot, total_plots, current_plot, {freq, Hw},'num_figs_in_row',num_figs_in_row);
legend("Hw");
title("frequency response function");
current_plot = current_plot + 1;

create_subplot(@plot, total_plots, current_plot,  {t, ft},'num_figs_in_row',num_figs_in_row);
legend("ft");
title("time domian of f");
current_plot = current_plot + 1;


create_subplot(@plot, total_plots, current_plot,  {t, ft_real},'num_figs_in_row',num_figs_in_row);
legend("ft_real");
title("time domian of ft real part");
current_plot = current_plot + 1;

create_subplot(@plot, total_plots, current_plot,{ t, ft_imag},'num_figs_in_row',num_figs_in_row);
legend("ft_imag");
title("time domian of ft imagin part");
current_plot = current_plot + 1;


create_subplot(@plot, total_plots, current_plot,{t, ft_real, t, ft_imag},'num_figs_in_row',num_figs_in_row);
legend("ft_real","ft_imag");
title("time domian of ft imagin part");
ylim([-100, 100])
current_plot = current_plot + 1;


% [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
% % plot(t,yn(1,:)/nodeshape)
% hold on
% % plot(t,u2dot_direct)
% plot(t, u2dot_re)

% [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
% plot(freq, yw(1, :))

% [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
% plot(freq, fw)

% [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
% plot(t, ft)

% [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
% plot(t, ft_real)

% [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
% plot(t, ft_imag)

% [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
% plot(t, ft_real)
% hold on
% plot(t, ft_imag)

% [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
% % plot(t,f_sim)
% hold on
% % plot(t,F_direct)
% plot(t, ft_real)
% ylim([-100, 100])

%% 3 Reconstruct displacement from acceleration Ole's method frequency domain integration
% X_u2dot = fftshift(fft(u2dot_measure));
% freq = linspace(-0.5, 0.5, length(t)) * fs;
% X_u = X_u2dot ./ (- (freq * 2 * pi) .^ 2);
% X_u(abs(freq) > 5) = 0;
%
% X_u(abs(freq) < 0.5) = 0;
% u_int3 = ifft(ifftshift(X_u));

