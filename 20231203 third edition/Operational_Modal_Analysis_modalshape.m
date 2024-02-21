clc; clear; close all;
run('CommonCommand.m');

%% Read viv Data
% 创建导入选项对象
opts = detectImportOptions('viv_in _the_paper.csv');

% 设置日期时间格式
% 假设日期时间格式为 'MM/dd/yyyy HH:mm'，请根据您的实际情况进行调整
opts = setvaropts(opts, 'startDate', 'InputFormat', 'MM/dd/yyyy HH:mm');
opts = setvaropts(opts, 'endDate', 'InputFormat', 'MM/dd/yyyy HH:mm');
opts = setvaropts(opts, 'startDate_update', 'InputFormat', 'MM/dd/yyyy HH:mm');
opts = setvaropts(opts, 'endDate_update', 'InputFormat', 'MM/dd/yyyy HH:mm');

vivTable = readtable('viv_in _the_paper.csv',opts);

%% load acceleration and wind data
VIV_sel = 11;
start_time = vivTable.startDate_update(VIV_sel);
end_time = vivTable.endDate_update(VIV_sel);
% start_time = vivTable.startDate(VIV_sel);
% end_time = vivTable.endDate(VIV_sel);
disp(vivTable.startDate(VIV_sel))
% start_time = datetime('2013-02-04 23:15:00', 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
% end_time = datetime('2013-02-04 23:45:00', 'InputFormat', 'yyyy-MM-dd HH:mm:ss'); % Example time range
% start_time = result.startDate;
% end_time = result.endDate;
[acc_result] = read_acceleration_data(start_time,end_time,input_data.acc_dir);
[wind_result] = read_wind_data(start_time,end_time,input_data.wind_dir);
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

y_ref =y([1,2,3],:);
fs = fs;
blockrows = 100
nb = 20

order = 1:50;

[lambda,Phi_id,cov_omega,cov_xi,cov_phi]=covssi_REF_unc(y,y_ref,fs,blockrows,nb,'order',order);
[omega_id,xi_id]=eig2modal(lambda);

%%



orders = order;
w_array = omega_id;
xi_array = xi_id;
w_tol = 0.01;
xi_tol = 0.05;
cov_w = cov_omega;
if isempty(cov_w)
    for k=1:length(orders)
        std_w_array{k}=zeros(size(w_array{k}));
    end
    exist_std_w=false;
else
    for k=1:length(orders)
        std_w_array{k}=cov_w{k}.^0.5;
    end
    exist_std_w=true;
end
cov_xi = cov_xi;
if isempty(cov_xi)
    for k=1:length(orders)
        std_xi_array{k}=zeros(size(xi_array{k}));
    end
    exist_std_xi=false;
else
    for k=1:length(orders)
        std_xi_array{k}=cov_xi{k}.^0.5;
    end
    exist_std_xi=true;
end

std_w_tol = 0.1;
std_xi_tol = 1;
xi_bounds = [0 0.2];
phi_array = Phi_id;

stabLevel = 1;
mpcw_tol = 20;


[w_array_stable,w_array_unstable,...
    xi_array_stable,xi_array_unstable,...
    std_w_array_stable,std_w_array_unstable,...
    std_xi_array_stable,std_xi_array_unstable,...
    phi_array_stable,phi_array_unstable,...
    order_array_stable,order_array_unstable...
    ]=filter_stable_modes(w_array,xi_array,w_tol,xi_tol,std_w_array,std_xi_array,std_w_tol,std_xi_tol,xi_bounds,phi_array,orders,stabLevel,'mpcw_tol',mpcw_tol);


%% deal with the stable plot data

all_w = []; % 存储所有 w_temp 值
all_orders = []; % 存储每个 w_temp 对应的 order
all_xi = []; % 存储每个 w_temp 对应的 xi_temp
all_std_w = [];
all_std_xi = [];

% 汇集所有 w_temp 和对应的 order
for k1 = 1:length(orders)
    w_temp = w_array_stable{k1};
    xi_temp = xi_array_stable{k1};
    std_w_temp = std_w_array_stable{k1};
    std_xi_temp = std_xi_array_stable{k1};
    all_w = [all_w; w_temp];
    all_xi = [all_xi; xi_temp];
    all_std_w = [all_std_w; std_w_temp];
    all_std_xi = [all_std_xi; std_xi_temp];
    all_orders = [all_orders; k1*ones(length(w_temp), 1)];
end
all_std_fre = all_std_w/(2*pi);
all_w_orders = [all_w,all_orders];
rng(1); % For reproducibility
% k_value = 14;
k_value = vivTable.k_value(VIV_sel);
[idx,w_kmean] = kmeans(all_w,k_value);

% 计算每个组的均值
unique_groups = unique(idx); % 获取所有唯一的组
group_means = arrayfun(@(g) mean(all_w(idx == g)), unique_groups);

% 对组均值进行排序，并获取排序后的索引
[~, sort_order] = sort(group_means);

% 生成新的 idx，映射排序后的组
new_idx = zeros(size(idx));
for i = 1:length(sort_order)
    new_idx(idx == unique_groups(sort_order(i))) = i;
end



list_fre = [];
list_xi = [];
list_std_w = [];
list_std_xi = [];
all_fre = all_w/(2*pi);
fre_kmean = w_kmean/(2*pi);
fre_kmean_sort = sort(fre_kmean);
for k1 = 1:k_value
    idx_temp = find(fre_kmean == fre_kmean_sort(k1));
    list_fre{k1} = all_fre(idx == idx_temp);
    list_xi{k1} = all_xi(idx == idx_temp);
    list_std_w{k1} = all_std_w(idx == idx_temp);
    list_std_xi{k1} = all_std_xi(idx == idx_temp);
    list_std_fre{k1} = list_std_w{k1}/(2*pi);
    mean_fre(k1,1) = mean(list_fre{k1});
    mean_xi(k1,1) = mean(list_xi{k1});
end

% % create a table | sequence | frequency | damping ratio |
% table_fre = table();
% table_fre.sequence = (1:k_value)';
% table_fre.frequency = mean_fre;
% table_fre.damping_ratio = mean_xi;
% table_fre.omega = mean_fre*2*pi;

%% calculate the error of the modes
fre_error=[];
xi_error=[];
for k1 = 1:length(unique_groups)
    list_fre_temp = list_fre{k1};
    list_std_fre_temp = list_std_fre{k1};
    fre_min = min(list_fre_temp-list_std_fre_temp);
    fre_max = max(list_fre_temp+list_std_fre_temp);
    fre_error(:,k1)=[fre_min;fre_max];
end

for k1 = 1:length(unique_groups)
    list_xi_temp = list_xi{k1};
    list_std_xi_temp = list_std_xi{k1};
    xi_min = min(list_xi_temp-list_std_xi_temp);
    xi_max = max(list_xi_temp+list_std_xi_temp);
    xi_error(:,k1)=[xi_min;xi_max];
end

%% modal updating

modesel = [2, 3, 5, 7, 13, 20, 22];
nmodes = length(modesel);
Result = ImportMK(nmodes, 'KMatrix.matrix', 'MMatrix.matrix', 'nodeondeck.txt', 'KMatrix.mapping', 'nodegap.txt', 'modesel', modesel,'showtext',true);
Fre_FEM = Result.Freq;

% create a table | sequence | frequency
table_FEM = table();
table_FEM.sequence = (1:length(Fre_FEM))';
table_FEM.frequency = Fre_FEM;

%% collect data
frequency=[];
frequency_sigma =[];
damping_ratio= [];
damping_ratio_sigma= [];
omega=[];
idx_ModalFEM=[];
MPC = [];
MPCW = [];
for k1 = 1:length(modesel)
    columnName = sprintf('mode%d_SSI', k1); % 动态生成列名
    columnName_seq = sprintf('mode%d_seq', k1); % 动态生成列名
    if or(vivTable.(columnName)(VIV_sel)==0,isnan(vivTable.(columnName)(VIV_sel)))
        continue
    else
        order_sel = vivTable.(columnName)(VIV_sel);
        mode_sel = vivTable.(columnName_seq)(VIV_sel);
        idx_ModalFEM = [idx_ModalFEM;Fre_FEM(k1)];
        % frequency
        temp_list = w_array{order_sel};
        temp = temp_list(mode_sel);
        frequency_temp = temp/2/pi;
        frequency = [frequency;frequency_temp];
        omega=[omega;temp]

        % frequency covariance
        temp_list = cov_omega{order_sel};
        temp = temp_list(mode_sel);
        temp = sqrt(temp);
        frequency_sigma_temp = temp/2/pi;
        frequency_sigma = [frequency_sigma;frequency_sigma_temp];

        % damping ratio
        temp_list = xi_array{order_sel};
        temp = temp_list(mode_sel);
        xi_temp = temp;
        damping_ratio = [damping_ratio;xi_temp];

        % damping ratio covariance
        temp_list = cov_xi{order_sel};
        temp = temp_list(mode_sel);
        temp = sqrt(temp);
        xi_sigma_temp = temp;
        damping_ratio_sigma = [damping_ratio_sigma;xi_sigma_temp];    

        % MPCW,MPC
         [MPCW_temp,~,~,MPC_temp]=rotate_modes(phi_array{order_sel}(:,mode_sel));
         MPC= [MPC;MPC_temp];
         MPCW =[MPCW;MPCW_temp];
    end
end

% create a table | sequence | frequency | damping ratio |
table_fre = table();
table_fre.sequence = (1:length(omega))';
table_fre.frequency = frequency;
table_fre.damping_ratio = damping_ratio;
table_fre.omega = omega;
table_fre.frequency_cov = frequency_sigma;
table_fre.damping_ratio_cov = damping_ratio_sigma;
table_fre.idx_Fre_FEM=idx_ModalFEM;
table_fre.MPC=MPC;
table_fre.MPCW = MPCW;

%% test
% close all
% FEM_mode_sel = 1;
% mode_deck_plot = (Result.mode_deck(:,FEM_mode_sel));
% FEM_freq = Result.Freq(FEM_mode_sel);
% figure('Position', [[100, 100], 960, 540]);
% plot(Result.node_loc,mode_deck_plot)
% 
% sensor_loc = [578+1650/4,578+1650*2/4,578+1650*3/4];
% % find the modal shape of the sensor location
% [~,idx_sensor_loc] = min(abs(Result.node_loc-sensor_loc));
% phi_sensor_loc = mode_deck_plot(idx_sensor_loc);
% hold on
% scatter(sensor_loc,phi_sensor_loc,'filled')

%% plot
if 1



    stabplot(omega_id,xi_id,order,Phi_id,'cov_w',cov_omega,'cov_xi',cov_xi,'std_w_tol',0.1,'std_xi_tol',1);
    
    for k1 = 1:length(modesel)
        columnName = sprintf('mode%d_SSI', k1); % 动态生成列名
        columnName_seq = sprintf('mode%d_seq', k1); % 动态生成列名
        if or(vivTable.(columnName)(VIV_sel)==0,isnan(vivTable.(columnName)(VIV_sel)))
            continue
        else
            % animation
            Modal_analysis_mode_sel = vivTable.(columnName_seq)(VIV_sel);
            order_sel = vivTable.(columnName)(VIV_sel);
            w_sel = w_array{order_sel}(Modal_analysis_mode_sel);
            freq_sel = w_sel/(2*pi);
            xi_sel = xi_array{order_sel}(Modal_analysis_mode_sel);
            phi_sel = phi_array{order_sel}(:,Modal_analysis_mode_sel);
            
            loc = [1,2,3]';
            plot_modal_shape_animation(phi_sel,'loc',loc)
            
            FEM_mode_sel = k1;
            FEM_freq = Result.Freq(FEM_mode_sel);

            disp("FEM mode "+FEM_mode_sel+" : Frequency = "+FEM_freq+" Hz")
            disp("Mode shape of order "+order_sel+" and mode "+Modal_analysis_mode_sel+" : Frequency = "+freq_sel+" Hz, Damping ratio = "+xi_sel)
            
            % complex plane
            plot_mode_complex_plane(phi_sel,'rotate','yes');
            
            % bridge deck mode comparison
            mode_deck_plot = (Result.mode_deck(:,k1));
            sensor_loc = [578+1650/4,578+1650*2/4,578+1650*3/4];
            mode_deck = Result.mode_deck; mode_deck_re = Result.mode_deck_re; node_loc = Result.node_loc; nodeondeck = Result.nodeondeck;
            KMmapping = Result.Mapping;
            nodegap = Result.nodegap;
            mode_vec = Result.mode_vec;
            mode_vec = mode_vec(:,k1)
            phi_sensor_loc = FindModeShapewithLocation(sensor_loc,node_loc,nodeondeck,KMmapping,nodegap,mode_vec);
            % phi_sensor_loc = mode_deck_plot(idx_sensor_loc);
            figure('Position', [[100, 100], 960, 540]);
            plot(Result.node_loc,mode_deck_plot)
            hold on
            scatter(sensor_loc,phi_sensor_loc,'filled')

            
            

            phi_id = real(phi_sel);
            [~,loc_temp]= max(phi_id);
            phi_id_scale = phi_id*(phi_sensor_loc(loc_temp)/phi_id(loc_temp));
            scatter(sensor_loc,phi_id_scale,'filled')
            
            seq_in_table =find(table_fre.idx_Fre_FEM==FEM_freq);
            
            MAC_value(seq_in_table,1) = calculate_MAC(phi_sensor_loc, phi_sel);

            disp("Mode shape of FEM mode "+FEM_mode_sel+" : " + ...
                "Frequency = "+table_fre.frequency(seq_in_table)+" Hz, " + ...
                "Damping ratio = "+table_fre.damping_ratio(seq_in_table) + ", " + ...
                "MPCW = " + table_fre.MPCW(seq_in_table) +", " + ...
                "MPC = " +table_fre.MPC(seq_in_table)+", " + ...
                "MAC = " +MAC_value(seq_in_table));

            % test = 1;
         
            close all
        end

    end

    table_fre.MAC_value = MAC_value;
    % TODO：check
    %
    %% plot
    % 定义总子图数量
    total_plots = 16; % 或任何你需要的子图数量
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

    create_subplot(@gscatter, total_plots, current_plot, {all_w,all_orders,new_idx}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', false, 'firstfigure', false, 'holdon', holdon);
    legend('off')
    xlabel("Omega(rad/s)")
    current_plot = current_plot+1;

    for k1 = 1:k_value
        create_subplot(@scatter, total_plots, current_plot, {list_fre{k1},k1*ones(length(list_fre{k1}),1)}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', false, 'firstfigure', false, 'holdon', holdon);
        hold on
    end
    xlabel("Frequency(Hz)");ylabel("Cluster");title("Frequency")
    current_plot = current_plot+1;

    for k1 = 1:k_value
        create_subplot(@scatter, total_plots, current_plot, {list_xi{k1},k1*ones(length(list_xi{k1}),1)}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', false, 'firstfigure', false, 'holdon', holdon);
        hold on
    end
    xlabel("Damping ratio");ylabel("Cluster");title("Damping ratio")
    current_plot = current_plot+1;

    for k1 = 1:k_value
        create_subplot(@scatter, total_plots, current_plot, {k1*ones(length(list_xi{k1}),1),list_xi{k1}}, 'num_figs_in_row', num_figs_in_row, 'figWidthFactor', figWidthFactor, 'figPosition', figPosition, 'newfigure', false, 'firstfigure', false, 'holdon', holdon);
        hold on
    end
    ylabel("Damping ratio");xlabel("Cluster");title("Damping ratio")
    current_plot = current_plot+1;



    figure
    h = gscatter(all_fre,all_xi,new_idx);
    legend("off")
    xlim([0,0.5])

    hold on
    % 获取每个组的颜色，以便误差线与散点颜色匹配
    colors = get(h, 'Color');
    if iscell(colors)
        colors = cell2mat(colors);
    end

    % 绘制每个组的矩形
    for k1 = 1:size(fre_error, 2)
        fre_min = fre_error(1, k1);
        fre_max = fre_error(2, k1);
        xi_min = xi_error(1, k1);
        xi_max = xi_error(2, k1);

        % 矩形的位置和尺寸
        pos = [fre_min, xi_min, fre_max - fre_min, xi_max - xi_min];

        % 绘制矩形
        rectangle('Position', pos, 'EdgeColor', colors(k1, :), 'LineWidth', 1.5);
    end

    scatter(mean_fre,mean_xi)
    scatter(Result.Freq,0.3/100*ones(size(Result.Freq)))
    % 关闭 hold 状态
    hold off;

    holdon = true;
end

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
