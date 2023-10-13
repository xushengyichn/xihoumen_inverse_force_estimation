function [result] = compute_dynamics_parameters(ncycle, t, Fa_temp, vel_temp, freq_temp, dis_temp , locs, Wind_Data)


    
    k2 = 1;
    
    for k1 = 1:ncycle:(length(locs) - ncycle + 1)
        if k1 + ncycle > length(locs)
            break;
        end

        dt_duration = t(locs(k1 + ncycle)) - t(locs(k1));
        dt = seconds(dt_duration);

        timeIntervals = seconds(t(locs(k1):locs(k1 + ncycle)) - t(locs(k1)));

        work(k2) = trapz(timeIntervals, Fa_temp(locs(k1):locs(k1 + ncycle)) .* vel_temp(locs(k1):locs(k1 + ncycle)));
        
        sum_amplitudes = 0;
        for cycle_offset = 0:(ncycle - 1)
            current_dis = dis_temp(locs(k1 + cycle_offset):locs(k1 + cycle_offset + 1));
            current_amp = (max(current_dis) - min(current_dis)) / 2;
            sum_amplitudes = sum_amplitudes + current_amp;
        end

        amp(k2) = sum_amplitudes / ncycle;
        
        f(k2) = (freq_temp(locs(k1)) + freq_temp(locs(k1 + ncycle))) / 2;
        omega(k2) = 2 * pi * f(k2);
        c(k2) = work(k2) / pi / amp(k2) ^ 2 / omega(k2) / ncycle;
        zeta_all(k2) = -c(k2) / 2 / omega(k2);
        
        timestamp_cycle(k2) = mean(t(locs(k1):locs(k1 + ncycle)));
        wind_color(k2) = interp1(Wind_Data.resultsTable_UA4.Time_Start, Wind_Data.resultsTable_UA4.U, timestamp_cycle(k2));
        
        k2 = k2 + 1;
    end
    
    result.amp = amp;
    result.zeta_all = zeta_all;

end
