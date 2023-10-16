function [ifq_interpolated] = instfreq_samelength(fs, x_k_k, f_keep, t, varargin)
    % YOURFUNCTIONNAME - 提供函数的简短描述
    % 这个函数根据提供的输入计算ifq_interpolated,通过外插保证时间长度不变
    %
    % Inputs:
    % fs - 输入描述
    % x_k_k - 输入描述，一个向量
    % f_keep - 输入描述
    % t - 输入描述
    % varargin - 可选的频率分辨率参数
    %
    % Outputs:
    % ifq_interpolated - 输出描述

    % % 默认频率分辨率
    % freqRes = 0.005;
    
    % % 如果提供了可选参数，则使用提供的频率分辨率
    % if nargin > 4
    %     freqRes = varargin{1};
    % end

    p = inputParser;
    addParameter(p,'showtext',true,@islogical);
    addParameter(p,'showplot',true,@islogical);
    addParameter(p,'freqRes',0.005,@isnumeric);% default freqRes is 0.005Hz
    addParameter(p,'filterstyle','fft',@ischar);% fft use the function by myself, bandpass use the function in matlab. fft is much faster than bandpass but may not be accurate
    parse(p,varargin{:});
    showtext = p.Results.showtext;
    showplot = p.Results.showplot;
    filterstyle = p.Results.filterstyle;
    freqRes = p.Results.freqRes;
    
    % 您的代码
        % 进行带通滤波
    switch filterstyle
        case 'fft'
            x_k_k_filtered = fft_filter(fs, x_k_k', f_keep);
        case 'bandpass'
            x_k_k_filtered = bandpass(x_k_k', f_keep, fs)';
        otherwise
            error('Invalid filter style. Choose either "fft" or "bandpass".')
    end
    
    [p, fd, td] = pspectrum(x_k_k_filtered, t, 'spectrogram', 'FrequencyResolution', freqRes);
    [ifq, t1] = instfreq(p, fd, td);
    
    inside_indices = t >= min(t1) & t <= max(t1);
    t_inside = t(inside_indices);
    ifq_interpolated_inside = interp1(t1, ifq, t_inside, 'linear');
    
    t_outside_left = t(t < min(t1));
    t_outside_right = t(t > max(t1));
    ifq_extrapolated_left = repmat(ifq(1), size(t_outside_left));
    ifq_extrapolated_right = repmat(ifq(end), size(t_outside_right));
    
    ifq_interpolated = [ifq_extrapolated_left(:)', ifq_interpolated_inside(:)', ifq_extrapolated_right(:)'];
end
