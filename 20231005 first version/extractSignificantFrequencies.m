function [significantFreqs, significantPeaks, significantLocs] = extractSignificantFrequencies(fs, Fa,varargin)
     p = inputParser;
    addParameter(p,'showtext',true,@islogical);
    addParameter(p,'showplot',true,@islogical);
    parse(p,varargin{:});
    showtext = p.Results.showtext;
    showplot = p.Results.showplot;
    
    [f, magnitude] = fft_transform(fs, Fa);
    df = f(2) - f(1);
    d_temp = 0.01 / df;
    [peaks, locs] = findpeaks(magnitude, 'MinPeakDistance', d_temp);
    
    significantFreqs = [];
    significantPeaks = [];
    significantLocs = [];
    
    if ~isempty(peaks)
        [sortedPeaks, sortIndex] = sort(peaks, 'descend');
        sortedLocs = locs(sortIndex);
        
        % Retain peaks that are above 5% of the highest peak
        threshold = 0.05 * sortedPeaks(1);
        aboveThresholdIdx = find(sortedPeaks >= threshold);
        
        % Keep the peaks and locations above the threshold
        significantPeaks = sortedPeaks(aboveThresholdIdx);
        significantLocs = sortedLocs(aboveThresholdIdx);
        
        % Save frequencies of the peaks that meet the conditions
        significantFreqs = f(significantLocs);
    end
    
    % 如果 shouldPlot 为真，则生成图形
    if showplot
        figure;
        plot(f, magnitude);
        hold on;
        scatter(f(significantLocs), significantPeaks);
        hold off;
        xlim([0, 0.5]);
        title('Significant Frequencies');
    end
end
