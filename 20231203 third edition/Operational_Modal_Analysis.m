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


y = [AC2_1;AC2_3;AC3_1;AC3_3;AC4_1;AC4_3];

%% SSI

y_ref =y([2 3],:);
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

% figure
% for k1 = 1:length(orders)
%     w_temp = w_array_stable{k1};
%     if isempty(w_temp)
%         continue;
%     end
%     y_temp = k1*ones(length(w_temp),1);
%     scatter(w_temp,y_temp);
%     hold on
% end

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
k_value = 20;
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

% create a table | sequence | frequency | damping ratio |
table_fre = table();
table_fre.sequence = (1:k_value)';
table_fre.frequency = mean_fre;
table_fre.damping_ratio = mean_xi;
table_fre.omega = mean_fre*2*pi;

%% modal updating
modesel = [2, 3, 5, 6, 7, 9, 15, 21, 23, 29, 33, 39, 44, 45];
nmodes = length(modesel);
Result = ImportMK(nmodes, 'KMatrix.matrix', 'MMatrix.matrix', 'nodeondeck.txt', 'KMatrix.mapping', 'nodegap.txt', 'modesel', modesel,'showtext',true);
Fre_FEM = Result.Freq;

% create a table | sequence | frequency
table_FEM = table();
table_FEM.sequence = (1:length(Fre_FEM))';
table_FEM.frequency = Fre_FEM;

% for k1 = 1:k_value
%     table_sel = table_fre(k1,:);
%     disp(table_FEM);
%     dispstr = sprintf('Cluster %d : Frequency = %f, please select the closest frequency from the above list',k1,table_sel.frequency);
%     disp(dispstr)
%     prompt = 'Please input the index of the closest frequency: ';
%     input_temp = input(prompt);
%     if input_temp == 0
%         idx_ModalFEM(k1,1) = 0; % if the frequency is out of range, set the index to 0
%     else
%         idx_ModalFEM(k1,1) = Fre_FEM(input_temp);
%     end
% end

idx_ModalFEM = zeros(size(table_fre, 1), 1); % 初始化 idx_ModalFEM

for k1 = 1:k_value
    table_sel = table_fre(k1,:);
    frequency = table_sel.frequency; % 当前频率
    
    % 寻找最接近的频率
    if frequency > max(Fre_FEM)
        % 如果频率超过了table_FEM的最大值
        idx_ModalFEM(k1) = 0;
    else
        % 计算与Fre_FEM中每个频率的差值，并找到最小的差值
        [~, closestIdx] = min(abs(Fre_FEM - frequency));
        % 记录最接近的频率值
        idx_ModalFEM(k1) = Fre_FEM(closestIdx);
        Fre_FEM(closestIdx)=[]; % remove duplicated value
    end
end

% 将结果添加到table_fre中
table_fre.idx_ModalFEM = idx_ModalFEM;


% 替换日期时间字符串中的冒号和其他特殊字符
start_time.Format = 'dd_MMM_yyyy_HH_mm_ss';
formatted_start_time = string(start_time);

end_time.Format = 'dd_MMM_yyyy_HH_mm_ss';
formatted_end_time = string(end_time);


formatted_start_time = strrep(strrep(strrep(formatted_start_time, ':', '-'), ' ', '_'), '/', '-');
formatted_end_time = strrep(strrep(strrep(formatted_end_time, ':', '-'), ' ', '_'), '/', '-');

filename = "Modal_updating_"+formatted_start_time+"_"+formatted_end_time+".mat";
save(filename,'table_fre',"start_time","end_time");
disp(table_fre)
%% test
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
%% plot
if 0
    
    close all
    
    
    
    stabplot(omega_id,xi_id,order,Phi_id,'cov_w',cov_omega,'cov_xi',cov_xi,'std_w_tol',0.1,'std_xi_tol',1);
    
    
    
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
    % close all
    
    
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
    
    holdon = true;
end