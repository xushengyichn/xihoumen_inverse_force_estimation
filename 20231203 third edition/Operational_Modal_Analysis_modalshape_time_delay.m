clc; clear; close all;
run('CommonCommand.m');

%% Read viv Data
% 创建导入选项对象
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

%% load acceleration and wind data
VIV_sel = 11;
% start_time = vivTable.startDate_update(VIV_sel);
% end_time = vivTable.endDate_update(VIV_sel);
start_time = vivTable.startDate(VIV_sel);
end_time = vivTable.endDate(VIV_sel);
disp(vivTable.startDate(VIV_sel))

[acc_result] = read_acceleration_data(start_time,end_time,input_data.acc_dir);
[wind_result] = read_wind_data(start_time,end_time,input_data.wind_dir);
acc_result=acc_result.mergedData;

UA1 = wind_result.resultsTable_UA1;
UA2 = wind_result.resultsTable_UA2;
UA3 = wind_result.resultsTable_UA3;
UA4 = wind_result.resultsTable_UA4;
UA5 = wind_result.resultsTable_UA5;
UA6 = wind_result.resultsTable_UA6;

% dt = seconds(acc_result.Time(2)-acc_result.Time(1));
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
T_new = 1:1:length(AC2_1);


figure
scatter(UA1.Time_Start,UA1.U)
hold on
scatter(UA2.Time_Start,UA2.U)
scatter(UA3.Time_Start,UA3.U)
scatter(UA4.Time_Start,UA4.U)
scatter(UA5.Time_Start,UA5.U)
scatter(UA6.Time_Start,UA6.U)

figure
plot(T_new,AC2_1)
hold on
plot(T_new,AC2_3)
plot(T_new,AC3_1)
plot(T_new,AC3_3)
plot(T_new,AC4_1)
plot(T_new,AC4_3)
legend('AC2_1','AC2_3','AC3_1','AC3_3','AC4_1','AC4_3')
dt = 1/fs;
%% selection of the sensors' data
sensor_selection= string(vivTable.sensor_selection(VIV_sel));
sensor_selection = strsplit(sensor_selection, ';'); % 以分号为分隔符分割字符串
sensor_selection = str2double(sensor_selection); % 将字符串数组转换为double数组

% 检查1和2是否在数组中
contains1 = ismember(1, sensor_selection);
contains2 = ismember(2, sensor_selection);
AC2 = sel_sensor(AC2_1,AC2_3,contains1,contains2);

contains1 = ismember(3, sensor_selection);
contains2 = ismember(4, sensor_selection);
AC3 = sel_sensor(AC3_1,AC3_3,contains1,contains2);

contains1 = ismember(5, sensor_selection);
contains2 = ismember(6, sensor_selection);
AC4 = sel_sensor(AC4_1,AC4_3,contains1,contains2);

y = [AC2;AC3;AC4];



figure
plot(T_new,AC2)
hold on
plot(T_new,AC3)
plot(T_new,AC4)
legend('AC2','AC3','AC4')

%% SSI
if 1
    y_ref =y([1,2,3],:);
    fs = fs;
    blockrows = 100
    nb = 20
    
    order = 1:50;
    
    [lambda,Phi_id,cov_omega,cov_xi,cov_phi]=covssi_REF_unc(y,y_ref,fs,blockrows,nb,'order',order);
    [omega_id,xi_id]=eig2modal(lambda);
    stabplot(omega_id,xi_id,order,Phi_id,'cov_omega',cov_omega,'cov_xi',cov_xi,'std_omega_tol',0.1,'std_xi_tol',1);

end
%% sync the sensors
if 1
% Select the VIV mode and lambda
phi_prime=[];
lambda1=[];
% for k1 = 1:7
%     columnName = sprintf('mode%d_SSI', k1); % 动态生成列名
%     columnName_seq = sprintf('mode%d_seq', k1); % 动态生成列名
%     if or(vivTable.(columnName)(VIV_sel)==0,isnan(vivTable.(columnName)(VIV_sel)))
%         continue
%     else
%         phi_prime=[phi_prime,Phi_id{vivTable.(columnName)(VIV_sel)}(:,vivTable.(columnName_seq)(VIV_sel))];
%         lambda1 = [lambda1;lambda{vivTable.(columnName)(VIV_sel)}(vivTable.(columnName_seq)(VIV_sel))];
%     end
% end
switch VIV_sel
    case 1 
        phi_prime=Phi_id{42};
        phi_prime = phi_prime(:,2:7);
        lambda1=lambda{42};
        lambda1 = lambda1(2:7);
    case 2 
        phi_prime=Phi_id{40};
        phi_prime = phi_prime(:,1:7);
        lambda1=lambda{40};
        lambda1 = lambda1(1:7);
    case 3
        phi_prime=Phi_id{40};
        phi_prime = phi_prime(:,2:8);
        lambda1=lambda{40};
        lambda1 = lambda1(2:8);
    case 4
        phi_prime=Phi_id{40};
        phi_prime = phi_prime(:,2:7);
        lambda1=lambda{40};
        lambda1 = lambda1(2:7);
    case 5
        phi_prime=Phi_id{41};
        phi_prime = phi_prime(:,2:6);
        lambda1=lambda{41};
        lambda1 = lambda1(2:6);
    case 6
        phi_prime=Phi_id{43};
        phi_prime = phi_prime(:,3:7);
        lambda1=lambda{43};
        lambda1 = lambda1(3:7);
    case 7
        phi_prime=Phi_id{38};
        phi_prime = phi_prime(:,[1:3,6]);
        lambda1=lambda{38};
        lambda1 = lambda1([1:3,6]);
    case 8
        phi_prime=Phi_id{40};
        phi_prime = phi_prime(:,[3:5,9]);
        lambda1=lambda{40};
        lambda1 = lambda1([3:5,9]);
    case 9
        phi_prime=Phi_id{38};
        phi_prime = phi_prime(:,[1:5]);
        lambda1=lambda{38};
        lambda1 = lambda1([1:5]);
    case 10
        phi_prime=Phi_id{42};
        phi_prime = phi_prime(:,[1:6]);
        lambda1=lambda{42};
        lambda1 = lambda1([1:6]);
    case 11
        phi_prime=Phi_id{41};
        phi_prime = phi_prime(:,[1:3,6:7]);
        lambda1=lambda{41};
        lambda1 = lambda1([1:3,6:7]);
end

[dt_opt,psi_opt,mpcw_opt]=sync_mpcw(phi_prime,lambda1,1,{[2] [3]},10/dt,dt,'plot',true,'plotobj',true)

str= num2str(dt_opt(1))+","+num2str(dt_opt(2));
disp(str)

disp('maximize mpcw:')
fprintf('Signal 2 delay relative to Signal 1: %f seconds\n', dt_opt(1));
fprintf('Signal 3 delay relative to Signal 1: %f seconds\n', dt_opt(2));

end
%% set period
period = 2.5
%% sync the sensors using VIV data
if 1 
disp('no filter:')
[delay2_seconds, delay3_seconds,results] = sync_cov(AC2, AC3, AC4, fs,period);


% y(2,:)=circshift(y(2,:),-delay2_seconds/dt);
% y(3,:)=circshift(y(3,:),-delay3_seconds/dt);
% 
% 
% figure
% plot(results.lags2*dt,results.corr2)
% figure
% plot(results.lags3*dt,results.corr3)
% 
% str= num2str(delay2_seconds)+","+num2str(delay3_seconds);
% disp(str)
% 
% [delay2_seconds_verify, delay3_seconds_verify,results] = sync_cov(y(1,:), y(2,:), y(3,:), fs,period);

end

%% sync the sensors using VIV data
if 1 
    disp('after applying lowpass filter 0.1Hz')
    lowpassfreq = 0.1;
    y_lowpass(1,:) = lowpass(y(1,:),lowpassfreq,fs,'Steepness',0.9);
    y_lowpass(2,:) = lowpass(y(2,:),lowpassfreq,fs,'Steepness',0.9);
    y_lowpass(3,:) = lowpass(y(3,:),lowpassfreq,fs,'Steepness',0.9);

    [f, magnitude] = fft_transform(fs,y(1,:));
    figure
    plot(f,magnitude)
    
    [f, magnitude] = fft_transform(fs,y_lowpass(1,:));
    figure
    plot(f,magnitude)

    [delay2_seconds, delay3_seconds,results] = sync_cov(y_lowpass(1,:), y_lowpass(2,:), y_lowpass(3,:), fs,period);
    

    figure
    hold on
    plot(y_lowpass(1,:))
    plot(y_lowpass(2,:))
    plot(y_lowpass(3,:))

    figure
    plot(results.lags2*dt,results.corr2)
    figure
    plot(results.lags3*dt,results.corr3)

end


%% sync the sensors using VIV data
if 1 
    disp('after applying highpass filter 0.35Hz')
    highpassfreq = 0.35;
    y_highpass(1,:) = highpass(y(1,:),highpassfreq,fs,'Steepness',0.9);
    y_highpass(2,:) = highpass(y(2,:),highpassfreq,fs,'Steepness',0.9);
    y_highpass(3,:) = highpass(y(3,:),highpassfreq,fs,'Steepness',0.9);

    [f, magnitude] = fft_transform(fs,y(1,:));
    figure
    plot(f,magnitude)
    
    [f, magnitude] = fft_transform(fs,y_highpass(1,:));
    figure
    plot(f,magnitude)

    [delay2_seconds, delay3_seconds,results] = sync_cov(y_highpass(1,:), y_highpass(2,:), y_highpass(3,:), fs,period);
    
    figure
    hold on
    plot(y_highpass(1,:))
    plot(y_highpass(2,:))
    plot(y_highpass(3,:))

    figure
    plot(results.lags2*dt,results.corr2)
    figure
    plot(results.lags3*dt,results.corr3)

end

%% sync the sensors using VIV data
if 1 
    disp('after applying bandpass filter [0.3 0.35]Hz')
    bandpassfreq = [0.3,0.35];

    y_bandpass(1,:) = bandpass(y(1,:),bandpassfreq,fs,'Steepness',0.9);
    y_bandpass(2,:) = bandpass(y(2,:),bandpassfreq,fs,'Steepness',0.9);
    y_bandpass(3,:) = bandpass(y(3,:),bandpassfreq,fs,'Steepness',0.9);

    [f, magnitude] = fft_transform(fs,y(1,:));
    figure
    plot(f,magnitude)
    
    [f, magnitude] = fft_transform(fs,y_bandpass(1,:));
    figure
    plot(f,magnitude)

    [delay2_seconds, delay3_seconds,results] = sync_cov(y_bandpass(1,:), y_bandpass(2,:), y_bandpass(3,:), fs,period);
    
    figure
    hold on
    plot(y_bandpass(1,:))
    plot(y_bandpass(2,:))
    plot(y_bandpass(3,:))

    figure
    plot(results.lags2*dt,results.corr2)
    figure
    plot(results.lags3*dt,results.corr3)

end
    
tilefigs

%% functions
function data = sel_sensor(data1,data2,contains1,contains2)
% 基于条件执行不同的操作
if contains1 && ~contains2
    % 数组仅包含1的操作
    disp('Array contains only contains 1.');
    data=data1;
elseif ~contains1 && contains2
    % 数组仅包含2的操作
    disp('Array contains only contains 2.');
    data=data2;
elseif contains1 && contains2
    % 数组同时包含1和2的操作
    disp('Array contains both contains 1 and contains 2.');
    data = (data1+data2)/2;
else
    % 数组既不包含1也不包含2的操作
    error('Array contains neither contains 1 nor contains 2.');
end
end

function MAC_value = calculate_MAC(phi1, phi2)
% 确保振型是列向量
phi1 = phi1(:);
phi2 = phi2(:);

% 计算MAC值
numerator = abs(phi1'*phi2)^2;
denominator = (phi1'*phi1) * (phi2'*phi2);
MAC_value = numerator / denominator;
end
