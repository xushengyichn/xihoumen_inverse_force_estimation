clc; clear; close all;
run('CommonCommand.m');

opts = detectImportOptions('vivData.csv');

% 设置日期时间格式
% 假设日期时间格式为 'MM/dd/yyyy HH:mm'，请根据您的实际情况进行调整
opts = setvartype(opts, 'startDate', 'datetime'); % 确保变量类型为 datetime
opts = setvartype(opts, 'endDate', 'datetime'); % 确保变量类型为 datetime
opts = setvartype(opts, 'startDate_update', 'datetime'); % 确保变量类型为 datetime
opts = setvartype(opts, 'endDate_update', 'datetime'); % 确保变量类型为 datetime

opts = setvaropts(opts, 'startDate', 'InputFormat', 'MM/dd/yyyy HH:mm');
opts = setvaropts(opts, 'endDate', 'InputFormat', 'MM/dd/yyyy HH:mm');
opts = setvaropts(opts, 'startDate_update', 'InputFormat', 'MM/dd/yyyy HH:mm');
opts = setvaropts(opts, 'endDate_update', 'InputFormat', 'MM/dd/yyyy HH:mm');

vivTable = readtable('vivData.csv',opts);

FEM_fre = [0.0938;0.1013;0.1319;0.1779;0.2273;0.2718;0.3213];
VIV_sels = [2;3;4;5;6;7;8;9;10;12;16;17;18;19];

%% mode 1
updated_frequency_mode1=[];
updated_damping_ratio_mode1=[];
for k1 = 1:length(VIV_sels)
    VIV_sel = VIV_sels(k1);
    start_time = vivTable.startDate_update(VIV_sel);
    end_time = vivTable.endDate_update(VIV_sel);
    % start_time = datetime('2013-02-04 23:15:00', 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
    % end_time = datetime('2013-02-04 23:45:00', 'InputFormat', 'yyyy-MM-dd HH:mm:ss'); % Example time range
    
    % 替换日期时间字符串中的冒号和其他特殊字符
    start_time.Format = 'MM_dd_yyyy_HH_mm_ss';
    formatted_start_time = string(start_time);
    
    end_time.Format = 'MM_dd_yyyy_HH_mm_ss';
    formatted_end_time = string(end_time);
    
    
    formatted_start_time = strrep(strrep(strrep(formatted_start_time, ':', '-'), ' ', '_'), '/', '-');
    formatted_end_time = strrep(strrep(strrep(formatted_end_time, ':', '-'), ' ', '_'), '/', '-');
    
    filename = "Modal_updating_"+formatted_start_time+"_"+formatted_end_time+".mat";
    
    if exist(filename, 'file') == 2
        Modal_updating_file = load(filename);
    else
        error("The file "+filename+" does not exist. Please run the modal updating script first. The script is named as Operational_Modal_Analysis.m")
    end

    % find the row index of the VIV mode
    VIV_index = find(abs(Modal_updating_file.table_fre.idx_ModalFEM-FEM_fre(1))<1e-3);
    if ~isempty(VIV_index)
        updated_frequency_mode1= [updated_frequency_mode1;Modal_updating_file.table_fre.frequency(VIV_index)];
        updated_damping_ratio_mode1 = [updated_damping_ratio_mode1;Modal_updating_file.table_fre.damping_ratio(VIV_index)];
    end
end

%% mode 2
updated_frequency_mode2=[];
updated_damping_ratio_mode2=[];
for k1 = 1:length(VIV_sels)
    VIV_sel = VIV_sels(k1);
    start_time = vivTable.startDate_update(VIV_sel);
    end_time = vivTable.endDate_update(VIV_sel);
    % start_time = datetime('2013-02-04 23:15:00', 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
    % end_time = datetime('2013-02-04 23:45:00', 'InputFormat', 'yyyy-MM-dd HH:mm:ss'); % Example time range
    
    % 替换日期时间字符串中的冒号和其他特殊字符
    start_time.Format = 'MM_dd_yyyy_HH_mm_ss';
    formatted_start_time = string(start_time);
    
    end_time.Format = 'MM_dd_yyyy_HH_mm_ss';
    formatted_end_time = string(end_time);
    
    
    formatted_start_time = strrep(strrep(strrep(formatted_start_time, ':', '-'), ' ', '_'), '/', '-');
    formatted_end_time = strrep(strrep(strrep(formatted_end_time, ':', '-'), ' ', '_'), '/', '-');
    
    filename = "Modal_updating_"+formatted_start_time+"_"+formatted_end_time+".mat";
    
    if exist(filename, 'file') == 2
        Modal_updating_file = load(filename);
    else
        error("The file "+filename+" does not exist. Please run the modal updating script first. The script is named as Operational_Modal_Analysis.m")
    end

    % find the row index of the VIV mode
    VIV_index = find(abs(Modal_updating_file.table_fre.idx_ModalFEM-FEM_fre(2))<1e-3);
    if ~isempty(VIV_index)
        updated_frequency_mode2= [updated_frequency_mode2;Modal_updating_file.table_fre.frequency(VIV_index)];
        updated_damping_ratio_mode2 = [updated_damping_ratio_mode2;Modal_updating_file.table_fre.damping_ratio(VIV_index)];
    end
end

%% mode 3
updated_frequency_mode3=[];
updated_damping_ratio_mode3=[];
for k1 = 1:length(VIV_sels)
    VIV_sel = VIV_sels(k1);
    start_time = vivTable.startDate_update(VIV_sel);
    end_time = vivTable.endDate_update(VIV_sel);
    % start_time = datetime('2013-02-04 23:15:00', 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
    % end_time = datetime('2013-02-04 23:45:00', 'InputFormat', 'yyyy-MM-dd HH:mm:ss'); % Example time range
    
    % 替换日期时间字符串中的冒号和其他特殊字符
    start_time.Format = 'MM_dd_yyyy_HH_mm_ss';
    formatted_start_time = string(start_time);
    
    end_time.Format = 'MM_dd_yyyy_HH_mm_ss';
    formatted_end_time = string(end_time);
    
    
    formatted_start_time = strrep(strrep(strrep(formatted_start_time, ':', '-'), ' ', '_'), '/', '-');
    formatted_end_time = strrep(strrep(strrep(formatted_end_time, ':', '-'), ' ', '_'), '/', '-');
    
    filename = "Modal_updating_"+formatted_start_time+"_"+formatted_end_time+".mat";
    
    if exist(filename, 'file') == 2
        Modal_updating_file = load(filename);
    else
        error("The file "+filename+" does not exist. Please run the modal updating script first. The script is named as Operational_Modal_Analysis.m")
    end

    % find the row index of the VIV mode
    VIV_index = find(abs(Modal_updating_file.table_fre.idx_ModalFEM-FEM_fre(3))<1e-3);
    if ~isempty(VIV_index)
        updated_frequency_mode3= [updated_frequency_mode3;Modal_updating_file.table_fre.frequency(VIV_index)];
        updated_damping_ratio_mode3 = [updated_damping_ratio_mode3;Modal_updating_file.table_fre.damping_ratio(VIV_index)];
    end
end

%% mode 4
updated_frequency_mode4=[];
updated_damping_ratio_mode4=[];
for k1 = 1:length(VIV_sels)
    VIV_sel = VIV_sels(k1);
    start_time = vivTable.startDate_update(VIV_sel);
    end_time = vivTable.endDate_update(VIV_sel);
    % start_time = datetime('2013-02-04 23:15:00', 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
    % end_time = datetime('2013-02-04 23:45:00', 'InputFormat', 'yyyy-MM-dd HH:mm:ss'); % Example time range
    
    % 替换日期时间字符串中的冒号和其他特殊字符
    start_time.Format = 'MM_dd_yyyy_HH_mm_ss';
    formatted_start_time = string(start_time);
    
    end_time.Format = 'MM_dd_yyyy_HH_mm_ss';
    formatted_end_time = string(end_time);
    
    
    formatted_start_time = strrep(strrep(strrep(formatted_start_time, ':', '-'), ' ', '_'), '/', '-');
    formatted_end_time = strrep(strrep(strrep(formatted_end_time, ':', '-'), ' ', '_'), '/', '-');
    
    filename = "Modal_updating_"+formatted_start_time+"_"+formatted_end_time+".mat";
    
    if exist(filename, 'file') == 2
        Modal_updating_file = load(filename);
    else
        error("The file "+filename+" does not exist. Please run the modal updating script first. The script is named as Operational_Modal_Analysis.m")
    end

    % find the row index of the VIV mode
    VIV_index = find(abs(Modal_updating_file.table_fre.idx_ModalFEM-FEM_fre(4))<1e-3);
    if ~isempty(VIV_index)
        updated_frequency_mode4= [updated_frequency_mode4;Modal_updating_file.table_fre.frequency(VIV_index)];
        updated_damping_ratio_mode4 = [updated_damping_ratio_mode4;Modal_updating_file.table_fre.damping_ratio(VIV_index)];
    end
end


%% mode 5
updated_frequency_mode5=[];
updated_damping_ratio_mode5=[];
for k1 = 1:length(VIV_sels)
    VIV_sel = VIV_sels(k1);
    start_time = vivTable.startDate_update(VIV_sel);
    end_time = vivTable.endDate_update(VIV_sel);
    % start_time = datetime('2013-02-04 23:15:00', 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
    % end_time = datetime('2013-02-04 23:45:00', 'InputFormat', 'yyyy-MM-dd HH:mm:ss'); % Example time range
    
    % 替换日期时间字符串中的冒号和其他特殊字符
    start_time.Format = 'MM_dd_yyyy_HH_mm_ss';
    formatted_start_time = string(start_time);
    
    end_time.Format = 'MM_dd_yyyy_HH_mm_ss';
    formatted_end_time = string(end_time);
    
    
    formatted_start_time = strrep(strrep(strrep(formatted_start_time, ':', '-'), ' ', '_'), '/', '-');
    formatted_end_time = strrep(strrep(strrep(formatted_end_time, ':', '-'), ' ', '_'), '/', '-');
    
    filename = "Modal_updating_"+formatted_start_time+"_"+formatted_end_time+".mat";
    
    if exist(filename, 'file') == 2
        Modal_updating_file = load(filename);
    else
        error("The file "+filename+" does not exist. Please run the modal updating script first. The script is named as Operational_Modal_Analysis.m")
    end

    % find the row index of the VIV mode
    VIV_index = find(abs(Modal_updating_file.table_fre.idx_ModalFEM-FEM_fre(5))<1e-3);
    if ~isempty(VIV_index)
        updated_frequency_mode5= [updated_frequency_mode5;Modal_updating_file.table_fre.frequency(VIV_index)];
        updated_damping_ratio_mode5 = [updated_damping_ratio_mode5;Modal_updating_file.table_fre.damping_ratio(VIV_index)];
    end
end

%% mode 6
updated_frequency_mode6=[];
updated_damping_ratio_mode6=[];
for k1 = 1:length(VIV_sels)
    VIV_sel = VIV_sels(k1);
    start_time = vivTable.startDate_update(VIV_sel);
    end_time = vivTable.endDate_update(VIV_sel);
    % start_time = datetime('2013-02-04 23:15:00', 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
    % end_time = datetime('2013-02-04 23:45:00', 'InputFormat', 'yyyy-MM-dd HH:mm:ss'); % Example time range
    
    % 替换日期时间字符串中的冒号和其他特殊字符
    start_time.Format = 'MM_dd_yyyy_HH_mm_ss';
    formatted_start_time = string(start_time);
    
    end_time.Format = 'MM_dd_yyyy_HH_mm_ss';
    formatted_end_time = string(end_time);
    
    
    formatted_start_time = strrep(strrep(strrep(formatted_start_time, ':', '-'), ' ', '_'), '/', '-');
    formatted_end_time = strrep(strrep(strrep(formatted_end_time, ':', '-'), ' ', '_'), '/', '-');
    
    filename = "Modal_updating_"+formatted_start_time+"_"+formatted_end_time+".mat";
    
    if exist(filename, 'file') == 2
        Modal_updating_file = load(filename);
    else
        error("The file "+filename+" does not exist. Please run the modal updating script first. The script is named as Operational_Modal_Analysis.m")
    end

    % find the row index of the VIV mode
    VIV_index = find(abs(Modal_updating_file.table_fre.idx_ModalFEM-FEM_fre(6))<1e-3);
    if ~isempty(VIV_index)
        updated_frequency_mode6= [updated_frequency_mode6;Modal_updating_file.table_fre.frequency(VIV_index)];
        updated_damping_ratio_mode6 = [updated_damping_ratio_mode6;Modal_updating_file.table_fre.damping_ratio(VIV_index)];
    end
end


%% VIV mode
for k1 = 1:length(VIV_sels)
    VIV_sel = VIV_sels(k1);
    start_time = vivTable.startDate_update(VIV_sel);
    end_time = vivTable.endDate_update(VIV_sel);
    % start_time = datetime('2013-02-04 23:15:00', 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
    % end_time = datetime('2013-02-04 23:45:00', 'InputFormat', 'yyyy-MM-dd HH:mm:ss'); % Example time range
    
    % 替换日期时间字符串中的冒号和其他特殊字符
    start_time.Format = 'MM_dd_yyyy_HH_mm_ss';
    formatted_start_time = string(start_time);
    
    end_time.Format = 'MM_dd_yyyy_HH_mm_ss';
    formatted_end_time = string(end_time);
    
    
    formatted_start_time = strrep(strrep(strrep(formatted_start_time, ':', '-'), ' ', '_'), '/', '-');
    formatted_end_time = strrep(strrep(strrep(formatted_end_time, ':', '-'), ' ', '_'), '/', '-');
    
    filename = "Modal_updating_"+formatted_start_time+"_"+formatted_end_time+".mat";
    
    if exist(filename, 'file') == 2
        Modal_updating_file = load(filename);
    else
        error("The file "+filename+" does not exist. Please run the modal updating script first. The script is named as Operational_Modal_Analysis.m")
    end

    % find the row index of the VIV mode
    VIV_index = find(abs(Modal_updating_file.table_fre.idx_ModalFEM-FEM_fre(7))<1e-3);
    updated_frequency(k1,1) = Modal_updating_file.table_fre.frequency(VIV_index);
    updated_damping_ratio(k1,1) = Modal_updating_file.table_fre.damping_ratio(VIV_index);
end

result_table = table(updated_frequency,updated_damping_ratio);
mean_xi=mean(result_table.updated_damping_ratio);
disp(result_table)
result.updated_frequency=updated_frequency;
result.updated_damping_ratio = updated_damping_ratio;

zeta(1,1) = mean(updated_damping_ratio_mode1);
zeta(2,1) = mean(updated_damping_ratio_mode2);
zeta(3,1) = mean(updated_damping_ratio_mode3);
zeta(4,1) = mean(updated_damping_ratio_mode4);
zeta(5,1) = mean(updated_damping_ratio_mode5);
zeta(6,1) = mean(updated_damping_ratio_mode6);
zeta(7,1) = mean(updated_damping_ratio);

frequency(1,1) = mean(updated_frequency_mode1);
frequency(2,1) = mean(updated_frequency_mode2);
frequency(3,1) = mean(updated_frequency_mode3);
frequency(4,1) = mean(updated_frequency_mode4);
frequency(5,1) = mean(updated_frequency_mode5);
frequency(6,1) = mean(updated_frequency_mode6);
frequency(7,1) = mean(updated_frequency);

result_table_mode = table(frequency,zeta);
disp(result_table_mode)

save("zeta_update.mat","zeta")

