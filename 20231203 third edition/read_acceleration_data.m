function [result] = read_acceleration_data(start_time,end_time,dirName)
     % 如果没有输入参数，执行以下测试或调试代码
        if nargin == 0
            clc;clear;close all;
            disp('Running tests...');
            start_time = datetime('2013-02-15 13:00:00', 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
            end_time = datetime('2013-02-17 05:00:00', 'InputFormat', 'yyyy-MM-dd HH:mm:ss'); % Example time range
            dirName = 'F:\test\result\';
            [result] = read_acceleration_data(start_time,end_time,dirName);
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
                if y == year(start_time) && m == month(start_time)
                    start_day = day(start_time);
                    start_hour = hour(start_time);
                else
                    start_day = 1;
                    start_hour = 0;
                end
                
                if y == year(end_time) && m == month(end_time)
                    end_day = day(end_time);
                    end_hour = hour(end_time);
                else
                    end_day = eomday(y, m);  % Get the last day of the month
                    end_hour = 23;
                end
                
                for d = start_day:end_day
                    if d == end_day
                        last_hour = end_hour;
                    else
                        last_hour = 23;
                    end
                    
                    for h = start_hour:last_hour
                        fileName = sprintf('%d-%02d-%02d %02d.mat', y, m, d, h);
                        fileNames{end+1} = fileName;
                    end
                    start_hour = 0;  % Reset for the next day
                end
            end
        end




        % 循环读取每个文件
    tableNames = {'mergedData'};
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


            timestamps = tableData.Time;

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
