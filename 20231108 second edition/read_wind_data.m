function [result] = read_wind_data(start_time,end_time,dirName)
     % 如果没有输入参数，执行以下测试或调试代码
    if nargin == 0
        clc;clear;close all;
        disp('Running tests...');
        start_time = datetime('2013-03-01 00:00:00', 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
        end_time = datetime('2013-03-05 00:00:00', 'InputFormat', 'yyyy-MM-dd HH:mm:ss'); % 改为跨年时间
        dirName = 'F:\test\result_wind_10min\';
        [result] = read_wind_data(start_time,end_time,dirName);
        disp('Tests completed.');
        return;
    end
    
    % 生成要读取的文件名列表
    fileNames = {};
    
    years = year(start_time):year(end_time);
    for y = years
        if y == year(start_time) % 起始年份
            start_month = month(start_time);
        else
            start_month = 1;
        end
        
        if y == year(end_time) % 结束年份
            end_month = month(end_time);
        else
            end_month = 12;
        end
    
        for m = start_month:end_month
            fileName = sprintf('windData_results_%d_%02d_10min.mat', y, m);
            fileNames{end+1} = fileName;
        end
    end


        % 循环读取每个文件
    tableNames = {'resultsTable_UA1', 'resultsTable_UA2', 'resultsTable_UA3', 'resultsTable_UA4', 'resultsTable_UA5', 'resultsTable_UA6'};
    result = struct();

    for f = 1:length(fileNames)
        matFile = fullfile(dirName, fileNames{f});
        if ~isfile(matFile)
            continue; % 如果文件不存在，则跳过
        end
        data = load(matFile);

        for i = 1:length(tableNames)
            tblName = tableNames{i};
            if ~isfield(data, tblName)
                continue; % 如果表格不存在，则跳过
            end
            
            tableData = data.(tblName);


            timestamps = tableData.Time_Start;

            % 查找在给定时间范围内的数据
            mask = (timestamps >= start_time) & (timestamps <= end_time);
            filteredData = tableData(mask, :);

            % 将筛选后的数据合并到结果结构体中
            if ~isfield(result, tblName)
                result.(tblName) = filteredData;
            else
                result.(tblName) = [result.(tblName); filteredData];
            end
        end
    end
end
