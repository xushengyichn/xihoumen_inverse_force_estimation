function filtered_signals = filterSignals(input_signal, frequencies, fs, bandwidth)
    % filterSignals - 对输入信号进行带通滤波
    % 
    % Inputs:
    % input_signal - 待滤波的信号
    % frequencies - 一个包含要保留的频率的数组
    % fs - 采样频率
    % bandwidth - 带宽
    %
    % Outputs:
    % filtered_signals - 一个cell数组，包含滤波后的信号

    if nargin < 4
        bandwidth = 0.01; % 如果未指定带宽，则使用默认值
    end

    filtered_signals = cell(1, length(frequencies));
    
    for k = 1:length(frequencies)
        freq = frequencies(k); % 获取当前频率
        
        % 定义带通滤波器的频率范围
        low_freq = max(0, freq - bandwidth/2);
        high_freq = min(fs/2, freq + bandwidth/2); 
        
        % 进行带通滤波
        filtered_signal = bandpass(input_signal, [low_freq, high_freq], fs);
        
        % 保存滤波后的信号
        filtered_signals{k} = filtered_signal;
    end
end
