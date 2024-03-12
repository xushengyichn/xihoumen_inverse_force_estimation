function [delay2_seconds, delay3_seconds,results] = sync_cov(signal1, signal2, signal3, sampleRate,maxDelaySeconds)
    % 定义最大延时时间（秒）

    % 将最大延时时间转换为样本数
    maxDelaySamples = round(maxDelaySeconds * sampleRate);
    
    % 计算信号1与信号2的互相关，限制lags范围
    [corr2, lags2] = xcorr(signal2, signal1, maxDelaySamples, 'coeff');
    % 找到互相关的最大值位置
    [~, I2] = max(corr2);
    % 延时是最大值位置对应的lag
    delay2_samples = lags2(I2);

    % 计算信号1与信号3的互相关，限制lags范围
    [corr3, lags3] = xcorr(signal3, signal1, maxDelaySamples, 'coeff');
    % 找到互相关的最大值位置
    [~, I3] = max(corr3);
    % 延时是最大值位置对应的lag
    delay3_samples = lags3(I3);

    % 将延时从样本数转换为秒
    delay2_seconds = int32(delay2_samples / sampleRate);
    delay3_seconds = int32(delay3_samples / sampleRate);
    results.corr2=corr2;
    results.corr3=corr3;
    results.lags2=lags2;
    results.lags3=lags3;

% 打印延时结果（以秒为单位）
    fprintf('Signal 2 delay relative to Signal 1: %f seconds\n', delay2_seconds);
    fprintf('Signal 3 delay relative to Signal 1: %f seconds\n', delay3_seconds);
end