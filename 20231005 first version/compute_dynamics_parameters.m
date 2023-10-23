%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: ShengyiXu xushengyichn@outlook.com
%Date: 2023-10-13 01:48:38
%LastEditors: ShengyiXu xushengyichn@outlook.com
%LastEditTime: 2023-10-13 20:46:01
%FilePath: \Exercises-for-Techniques-for-estimation-in-dynamics-systemsf:\git\xihoumen_inverse_force_estimation\20231005 first version\compute_dynamics_parameters.m
%Description: 
%
%Copyright (c) 2023 by ${git_name_email}, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [result] = compute_dynamics_parameters(ncycle, t, Fa_temp, vel_temp, freq_temp, dis_temp , locs,varargin)
    p = inputParser;
    addParameter(p,'showtext',false,@islogical);
    addParameter(p,'showplot',false,@islogical);
    parse(p,varargin{:});
    showtext = p.Results.showtext;
    showplot = p.Results.showplot;


    
    k2 = 1;
    
    overlap = 0.9; % 设置重叠的比例，例如0.5代表50%重叠
    step = round(ncycle * (1 - overlap)); % 计算每步应该移动的距离

    for k1 = 1:step:(length(locs) - ncycle + 1) % 更新k1的迭代方式
    % for k1 = 1:ncycle:(length(locs) - ncycle + 1)
        if k1 + ncycle > length(locs)
            break;
        end

        dt_duration = t(locs(k1 + ncycle)) - t(locs(k1));
        dt = seconds(dt_duration);
        
        t_cycle_mean(k2) = mean(t(locs(k1):locs(k1 + ncycle))); 

        timeIntervals = seconds(t(locs(k1):locs(k1 + ncycle)) - t(locs(k1)));

        work(k2) = trapz(timeIntervals, Fa_temp(locs(k1):locs(k1 + ncycle)) .* vel_temp(locs(k1):locs(k1 + ncycle)));
        
        sum_amplitudes = 0;
        for cycle_offset = 0:(ncycle - 1)
            current_dis = dis_temp(locs(k1 + cycle_offset):locs(k1 + cycle_offset + 1));
            current_amp = (max(current_dis) - min(current_dis)) / 2;
            % current_amp = rms(current_dis)*sqrt(2);
            sum_amplitudes = sum_amplitudes + current_amp;
        end

        amp(k2) = sum_amplitudes / ncycle;
        
        f(k2) = (freq_temp(locs(k1)) + freq_temp(locs(k1 + ncycle))) / 2;
        omega(k2) = 2 * pi * f(k2);
        % omega(k2) = 2 * pi * mean(freq_temp);
        c(k2) = work(k2) / pi / amp(k2) ^ 2 / omega(k2) / ncycle;
        zeta_all(k2) = -c(k2) / 2 / omega(k2);
        
        timestamp_cycle(k2) = mean(t(locs(k1):locs(k1 + ncycle)));
        % wind_color(k2) = interp1(Wind_Data.resultsTable_UA4.Time_Start, Wind_Data.resultsTable_UA4.U, timestamp_cycle(k2));
        
        k2 = k2 + 1;

        if showplot
            plot(t(locs(k1):locs(k1 + ncycle)),Fa_temp(locs(k1):locs(k1 + ncycle)) ...
                /max(abs(Fa_temp(locs(k1):locs(k1 + ncycle)))))
            hold on
            plot(t(locs(k1):locs(k1 + ncycle)),vel_temp(locs(k1):locs(k1 + ncycle)) ...
                /max(abs(vel_temp(locs(k1):locs(k1 + ncycle)))))
            title("Fa vs. vel, zeta = "+num2str(zeta_all(k2-1)))
            legend("Fa","vel")
            hold off
            % pause
        end
    end
    
    result.amp = amp;
    result.zeta_all = zeta_all;
    result.work = work;
    result.t_cycle_mean = t_cycle_mean;
end
