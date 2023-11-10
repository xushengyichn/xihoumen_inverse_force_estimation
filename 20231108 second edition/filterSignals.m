function filtered_signals = filterSignals(input_signal, frequencies, fs, bandwidth, varargin)
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
    
    p = inputParser;
    addParameter(p,'showtext',false,@islogical);
    addParameter(p,'showplot',false,@islogical);
    addParameter(p,'filterstyle','fft',@ischar);% fft use the function by myself, bandpass use the function in matlab. fft is much faster than bandpass but may not be accurate
    parse(p,varargin{:});
    showtext = p.Results.showtext;
    showplot = p.Results.showplot;
    filterstyle = p.Results.filterstyle;


    filtered_signals = cell(1, length(frequencies));
    
    for k = 1:length(frequencies)
        freq = frequencies(k); % 获取当前频率
        
        % 定义带通滤波器的频率范围
        low_freq = max(0, freq - bandwidth/2);
        high_freq = min(fs/2, freq + bandwidth/2); 
        
        % 进行带通滤波
        switch filterstyle
            case 'fft'
                filtered_signal = fft_filter(fs, input_signal, [low_freq, high_freq]);
            case 'bandpass'
                filtered_signal = bandpass(input_signal, [low_freq, high_freq], fs);
            case 'nofilter'
            filtered_signal = input_signal;
            otherwise
                error('Invalid filter style. Choose either "fft" or "bandpass".')
        end
        
        % 保存滤波后的信号
        filtered_signals{k} = filtered_signal;
    end
end
