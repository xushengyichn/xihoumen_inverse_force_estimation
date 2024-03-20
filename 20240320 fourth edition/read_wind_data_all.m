function [result] = read_wind_data_all(start_time,end_time,dirName)
     % 如果没有输入参数，执行以下测试或调试代码
    if nargin == 0
        clc;clear;close all;
        disp('Running tests...');
        start_time = datetime('2013-03-01 00:00:00', 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
        end_time = datetime('2013-03-01 03:00:00', 'InputFormat', 'yyyy-MM-dd HH:mm:ss'); % 改为跨年时间
        dirName = 'F:\test\result_wind';
        [result] = read_wind_data_all(start_time,end_time,dirName);
        disp('Tests completed.');
        return;
    end
    
      % 生成要读取的文件名列表
    fileNames = {};
    
    % 计算时间序列的整体开始和结束时间
    start_time_whole = start_time;
    end_time_whole = end_time;
    
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
    result = whole_table;

end
