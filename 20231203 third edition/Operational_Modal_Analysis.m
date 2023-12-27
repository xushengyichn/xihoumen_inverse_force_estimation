clc; clear; close all;
run('CommonCommand.m');

%% load VIV data
% n=4;
% displayDates=true;
% [result]=viv2013(n,displayDates);


%% load acceleration and wind data
start_time = datetime('2013-02-04 22:00:00', 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
end_time = datetime('2013-02-05 00:00:00', 'InputFormat', 'yyyy-MM-dd HH:mm:ss'); % Example time range
% start_time = result.startDate;
% end_time = result.endDate;
[acc_result] = read_acceleration_data(start_time,end_time,input.acc_dir);
[wind_result] = read_wind_data(start_time,end_time,input.wind_dir);
acc_result=acc_result.mergedData;

UA1 = wind_result.resultsTable_UA1;
UA2 = wind_result.resultsTable_UA2;
UA3 = wind_result.resultsTable_UA3;
UA4 = wind_result.resultsTable_UA4;
UA5 = wind_result.resultsTable_UA5;
UA6 = wind_result.resultsTable_UA6;


%% resample
% Øyvind's instruction: 2Hz lowpass filter + resample to 5 or 10 Hz
fs_original = 1/(seconds(acc_result.Time(2)-acc_result.Time(1)));
T_original = acc_result.Time;
AC2_1 = acc_result.AC2_1';
AC2_3 = acc_result.AC2_3';
AC3_1 = acc_result.AC3_1';
AC3_3 = acc_result.AC3_3';
AC4_1 = acc_result.AC4_1';
AC4_3 = acc_result.AC4_3';

lowpassfreq = 2;
AC2_1 = lowpass(AC2_1,lowpassfreq,fs_original);
AC2_3 = lowpass(AC2_3,lowpassfreq,fs_original);
AC3_1 = lowpass(AC3_1,lowpassfreq,fs_original);
AC3_3 = lowpass(AC3_3,lowpassfreq,fs_original);
AC4_1 = lowpass(AC4_1,lowpassfreq,fs_original);
AC4_3 = lowpass(AC4_3,lowpassfreq,fs_original);

fs = 5;
[P,Q]=rat(fs/fs_original);
AC2_1 = resample(AC2_1,P,Q);
AC2_3 = resample(AC2_3,P,Q);
AC3_1 = resample(AC3_1,P,Q);
AC3_3 = resample(AC3_3,P,Q);
AC4_1 = resample(AC4_1,P,Q);
AC4_3 = resample(AC4_3,P,Q);
T_new = T_original(1):seconds(1/fs):T_original(end);


y = [AC2_1;AC2_3;AC3_1;AC3_3;AC4_1;AC4_3];

%% SSI

y_ref =y([2 3],:);
fs = fs;
blockrows = 100
nb =20

order = 1:50;

[lambda,Phi_id,cov_omega,cov_xi,cov_phi]=covssi_REF_unc(y,y_ref,fs,blockrows,nb,'order',order);

%% test
% 
% y1 = y(3,:);
% [f, magnitude] = fft_transform(fs,y1);
% figure
% plot(f,log(magnitude))
% 
% figure
% plot(f,magnitude)

%%

close all

[omega_id,xi_id]=eig2modal(lambda);

stabplot(omega_id,xi_id,order,Phi_id,'cov_w',cov_omega,'cov_xi',cov_xi,'std_w_tol',0.1,'std_xi_tol',1);

%%

close all
% 定义总子图数量
total_plots = 50; % 或任何你需要的子图数量
current_plot = 1;
num_figs_in_row = [];
figWidthFactor = 1.5;
%     figPosition = [1080*2.5,100];
figPosition = [100, 100];
newfigure = true;
holdon = false;
figure
for k1 = 1:length(order)
    create_subplot(@scatter, total_plots, current_plot, {1:length(omega_id{k1}),omega_id{k1}}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', newfigure, 'firstfigure', false, 'holdon', holdon);
    current_plot = current_plot+1;
end

figure
for k1 = 1:length(order)
    create_subplot(@scatter, total_plots, current_plot, {1:length(omega_id{k1}),omega_id{k1}}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', newfigure, 'firstfigure', false, 'holdon', true);
    hold on
end

figure
current_plot = 1;
for k1 = 1:length(order)
    create_subplot(@scatter, total_plots, current_plot, {1:length(xi_id{k1}),xi_id{k1}}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', newfigure, 'firstfigure', false, 'holdon', holdon);
    ylim([0,1/100]);
    current_plot = current_plot+1;
end

figure
for k1 = 1:length(order)
    create_subplot(@scatter, total_plots, current_plot, {1:length(xi_id{k1}),xi_id{k1}}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', newfigure, 'firstfigure', false, 'holdon', true);
    hold on
end
ylim([0,1/100]);



%% plot
% 定义总子图数量
total_plots = 12; % 或任何你需要的子图数量
current_plot = 1;
num_figs_in_row = [];
figWidthFactor = 1.5;
%     figPosition = [1080*2.5,100];
figPosition = [100, 100];
newfigure = true;
holdon = false;

create_subplot(@plot, total_plots, current_plot, {acc_result.Time, acc_result.AC2_1}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', newfigure, 'firstfigure', true, 'holdon', holdon);
hold on
create_subplot(@plot, total_plots, current_plot, {T_new, AC2_1}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'holdon', holdon);
ylim([-60,60]);xlabel("Time");ylabel("Acceleration");title("AC2_1")
current_plot = current_plot+1;

create_subplot(@plot, total_plots, current_plot, {acc_result.Time, acc_result.AC2_3}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', newfigure, 'firstfigure', false, 'holdon', holdon);
hold on
create_subplot(@plot, total_plots, current_plot, {T_new, AC2_3}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'holdon', holdon);
ylim([-60,60]);xlabel("Time");ylabel("Acceleration");title("AC2_3")
current_plot = current_plot+1;

create_subplot(@plot, total_plots, current_plot, {acc_result.Time, acc_result.AC3_1}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', newfigure, 'firstfigure', false, 'holdon', holdon);
hold on
create_subplot(@plot, total_plots, current_plot, {T_new, AC3_1}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'holdon', holdon);
ylim([-60,60]);xlabel("Time");ylabel("Acceleration");title("AC3_1")
current_plot = current_plot+1;

create_subplot(@plot, total_plots, current_plot, {acc_result.Time, acc_result.AC3_3}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', newfigure, 'firstfigure', false, 'holdon', holdon);
hold on
create_subplot(@plot, total_plots, current_plot, {T_new, AC3_3}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'holdon', holdon);
ylim([-60,60]);xlabel("Time");ylabel("Acceleration");title("AC3_3")
current_plot = current_plot+1;

create_subplot(@plot, total_plots, current_plot, {acc_result.Time, acc_result.AC4_1}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', newfigure, 'firstfigure', false, 'holdon', holdon);
hold on
create_subplot(@plot, total_plots, current_plot, {T_new, AC4_1}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'holdon', holdon);
ylim([-60,60]);xlabel("Time");ylabel("Acceleration");title("AC4_1")
current_plot = current_plot+1;

create_subplot(@plot, total_plots, current_plot, {acc_result.Time, acc_result.AC4_3}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', newfigure, 'firstfigure', false, 'holdon', holdon);
hold on
create_subplot(@plot, total_plots, current_plot, {T_new, AC4_3}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'holdon', holdon);
ylim([-60,60]);xlabel("Time");ylabel("Acceleration");title("AC4_3")
current_plot = current_plot+1;

create_subplot(@plot, total_plots, current_plot, {UA1.Time_Start,UA1.U}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', false, 'firstfigure', false, 'holdon', holdon);
ylim([0,15]);xlabel("Time");ylabel("Wind Speed");title("UA1")
current_plot = current_plot+1;

create_subplot(@plot, total_plots, current_plot, {UA2.Time_Start,UA2.U}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', false, 'firstfigure', false, 'holdon', holdon);
ylim([0,15]);xlabel("Time");ylabel("Wind Speed");title("UA2")
current_plot = current_plot+1;

create_subplot(@plot, total_plots, current_plot, {UA3.Time_Start,UA3.U}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', false, 'firstfigure', false, 'holdon', holdon);
ylim([0,15]);xlabel("Time");ylabel("Wind Speed");title("UA3")
current_plot = current_plot+1;

create_subplot(@plot, total_plots, current_plot, {UA4.Time_Start,UA4.U}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', false, 'firstfigure', false, 'holdon', holdon);
ylim([0,15]);xlabel("Time");ylabel("Wind Speed");title("UA4")
current_plot = current_plot+1;

create_subplot(@plot, total_plots, current_plot, {UA5.Time_Start,UA5.U}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', false, 'firstfigure', false, 'holdon', holdon);
ylim([0,15]);xlabel("Time");ylabel("Wind Speed");title("UA5")
current_plot = current_plot+1;

create_subplot(@plot, total_plots, current_plot, {UA6.Time_Start,UA6.U}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', false, 'firstfigure', false, 'holdon', holdon);
ylim([0,15]);xlabel("Time");ylabel("Wind Speed");title("UA6")
current_plot = current_plot+1;

holdon = true;
