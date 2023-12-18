function [result] = read_wind_data_onetime(timestamps,duration,dirName)
     % 如果没有输入参数，执行以下测试或调试代码
    if nargin == 0
        clc;clear;close all;
        disp('Running tests...');
        run("CommonCommand.m")
        timestamp1 = datetime('2013-03-01 00:01:00', 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
        timestamp2 = datetime('2013-03-01 00:02:00', 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
        timestamps = [timestamp1 ;timestamp2];
        duration = minutes(10);
        dirName = 'F:\test\result_wind\';
        [result] = read_wind_data_onetime(timestamps,duration,dirName);
        disp('Tests completed.');
        return;
    end
    
    % 生成要读取的文件名列表
    fileNames = {};
    
    % 计算时间序列的整体开始和结束时间
    start_time_whole = min(timestamps) - duration / 2;
    end_time_whole = max(timestamps) + duration / 2;
    
    % 计算开始和结束时间之间的所有小时
    hours = floor(datenum(start_time_whole)*24):floor(datenum(end_time_whole)*24);
    % 初始化文件列表
    fileNames = cell(length(hours), 1);
    
    % 生成文件名
    for i = 1:length(hours)
        % 将小时转换回日期时间
        date_time = datestr(hours(i)/24, 'yyyy-mm-dd HH');
        
        % 生成文件名
        fileNames{i} = [date_time '-UANua.mat'];
    end

    
    %     % 循环读取每个文件
      whole_table = table();

    for f = 1:length(fileNames)
        matFile = fullfile(dirName, fileNames{f});
        if ~isfile(matFile)
            % continue; % 如果文件不存在，则跳过
            error("文件不存在");
        end
        data = load(matFile);

    % 假设每个文件中的变量名为 'windData'
        if isfield(data, 'mergedData')
            whole_table = [whole_table; data.mergedData]; % 合并表格
        else
            error("变量 'mergedData' 在文件中不存在: " + matFile);
        end
    end
    result_table=[];
    UA1 =[];
    UA2 =[];
    UA3 =[];
    UA4 =[];
    UA5 =[];
    UA6 =[];

    for k1  = 1:length(timestamps)
        timestamp = timestamps(k1);
        start_time = timestamp - duration / 2;
        end_time = timestamp + duration / 2;
        % 筛选出位于 start_time 和 end_time 之间的数据
        % 假设表格中有一个名为 'Timestamp' 的时间戳列
        segment_temp = whole_table(whole_table.Time >= start_time & whole_table.Time <= end_time, :);
        time_series = segment_temp.Time;
        U_x_1 = segment_temp.UA1_x;%north
        U_y_1 = segment_temp.UA1_y;%west
        U_z_1 = segment_temp.UA1_z;%up

        U_x_2 = segment_temp.UA2_x;%north
        U_y_2 = segment_temp.UA2_y;%west
        U_z_2 = segment_temp.UA2_z;%up

        U_x_3 = segment_temp.UA3_x;%north
        U_y_3 = segment_temp.UA3_y;%west
        U_z_3 = segment_temp.UA3_z;%up

        U_x_4 = segment_temp.UA4_x;%north
        U_y_4 = segment_temp.UA4_y;%west
        U_z_4 = segment_temp.UA4_z;%up

        U_x_5 = segment_temp.UA5_x;%north
        U_y_5 = segment_temp.UA5_y;%west
        U_z_5 = segment_temp.UA5_z;%up

        U_x_6 = segment_temp.UA6_x;%north
        U_y_6 = segment_temp.UA6_y;%west
        U_z_6 = segment_temp.UA6_z;%up

        % UA3是靠南安装的（假定西堠门桥45°走向，测量风向为45°-225°），UA4是靠北安装的（假定西堠门桥45°走向，测量风向为225°-360°，0°-45°）

        result_1=cal_wind_property(U_x_1,U_y_1,U_z_1,45);
        resultsTable_UA1 = struct2table(result_1); % Convert struct result to table


        result_2=cal_wind_property(U_x_2,U_y_2,U_z_2,45);
        resultsTable_UA2 = struct2table(result_2); % Convert struct result to table


        result_3=cal_wind_property(U_x_3,U_y_3,U_z_3,45);
        resultsTable_UA3 = struct2table(result_3); % Convert struct result to table


        result_4=cal_wind_property(U_x_4,U_y_4,U_z_4,45);
        resultsTable_UA4 = struct2table(result_4); % Convert struct result to table

        result_5=cal_wind_property(U_x_5,U_y_5,U_z_5,45);
        resultsTable_UA5 = struct2table(result_5); % Convert struct result to table


        result_6=cal_wind_property(U_x_6,U_y_6,U_z_6,45);
        resultsTable_UA6 = struct2table(result_6); % Convert struct result to table

        resultsTable_UA1.timestamp = timestamp;
        resultsTable_UA2.timestamp = timestamp;
        resultsTable_UA3.timestamp = timestamp;
        resultsTable_UA4.timestamp = timestamp;
        resultsTable_UA5.timestamp = timestamp;
        resultsTable_UA6.timestamp = timestamp;

        resultsTable_UA1 = movevars(resultsTable_UA1, 'timestamp', 'Before', resultsTable_UA1.Properties.VariableNames{1});
        resultsTable_UA2 = movevars(resultsTable_UA2, 'timestamp', 'Before', resultsTable_UA2.Properties.VariableNames{1});
        resultsTable_UA3 = movevars(resultsTable_UA3, 'timestamp', 'Before', resultsTable_UA3.Properties.VariableNames{1});
        resultsTable_UA4 = movevars(resultsTable_UA4, 'timestamp', 'Before', resultsTable_UA4.Properties.VariableNames{1});
        resultsTable_UA5 = movevars(resultsTable_UA5, 'timestamp', 'Before', resultsTable_UA5.Properties.VariableNames{1});
        resultsTable_UA6 = movevars(resultsTable_UA6, 'timestamp', 'Before', resultsTable_UA6.Properties.VariableNames{1});

        UA1 = [UA1;resultsTable_UA1];
        UA2 = [UA2;resultsTable_UA2];
        UA3 = [UA3;resultsTable_UA3];
        UA4 = [UA4;resultsTable_UA4];
        UA5 = [UA5;resultsTable_UA5];
        UA6 = [UA6;resultsTable_UA6];


        
        % result_table_temp=[resultsTable_UA1;resultsTable_UA2;resultsTable_UA3;resultsTable_UA4;resultsTable_UA5;resultsTable_UA6;];
        
        % % 假设 strArray 是一个单元数组，长度与 T 的行数相同
        % strArray = {'UA1'; 'UA2'; 'UA3'; 'UA4'; 'UA5';'UA6'}; % 根据需要填充
        
        % % 添加到表格
        % result_table_temp.property = strArray;
        % result_table_temp.timestamp =  repmat(timestamp, height(result_table_temp), 1);
        % % 将 NewColumn 移动到第一列
        % result_table_temp = movevars(result_table_temp, 'property', 'Before', result_table_temp.Properties.VariableNames{1});
        % result_table_temp = movevars(result_table_temp, 'timestamp', 'Before', result_table_temp.Properties.VariableNames{1});

        % result_table=[result_table;result_table_temp];
    end
        
    result.UA1 = UA1;
    result.UA2 = UA2;
    result.UA3 = UA3;
    result.UA4 = UA4;
    result.UA5 = UA5;
    result.UA6 = UA6;

end
