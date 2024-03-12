clc; clear; close all;
run('CommonCommand.m');

%% Read viv Data
% 创建导入选项对象
opts = detectImportOptions('viv_in_the_paper.csv');

% 设置日期时间格式
% 假设日期时间格式为 'MM/dd/yyyy HH:mm'，请根据您的实际情况进行调整
opts = setvaropts(opts, 'startDate', 'InputFormat', 'MM/dd/yyyy HH:mm');
opts = setvaropts(opts, 'endDate', 'InputFormat', 'MM/dd/yyyy HH:mm');
opts = setvaropts(opts, 'startDate_update', 'InputFormat', 'MM/dd/yyyy HH:mm');
opts = setvaropts(opts, 'endDate_update', 'InputFormat', 'MM/dd/yyyy HH:mm');

vivTable = readtable('viv_in_the_paper.csv',opts);

%% load acceleration and wind data
VIV_sel = 3;
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
phi_prime=Phi_id{40};
phi_prime = phi_prime(:,2:8);
lambda1=lambda{40};
lambda1 = lambda1(2:8);

[dt_opt,psi_opt,mpcw_opt]=sync_mpcw(phi_prime,lambda1,1,{[2] [3]},10/dt,dt,'plot',true)

str= num2str(dt_opt(1))+","+num2str(dt_opt(2));
disp(str)


end
%% sync the sensors using VIV data
if 1 

[delay2_seconds, delay3_seconds,results] = sync_cov(AC2, AC3, AC4, fs,10);


y(2,:)=circshift(y(2,:),-delay2_seconds/dt);
y(3,:)=circshift(y(3,:),-delay3_seconds/dt);


figure
plot(results.lags2/dt,results.corr2)
figure
plot(results.lags3/dt,results.corr3)

str= num2str(delay2_seconds)+","+num2str(delay3_seconds);
disp(str)

[delay2_seconds_verify, delay3_seconds_verify,results] = sync_cov(y(1,:), y(2,:), y(3,:), fs,600);

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
