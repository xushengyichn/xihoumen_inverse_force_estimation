function [result_Main]=Cal_aero_damping_ratio(input,varargin)
    
% input.ncycle = 1;%计算气动阻尼时n个周期算一次阻尼比

% input.wind_dir = "F:\test\result_wind_10min";
% input.wind_dir = "Z:\Drive\Backup\SHENGYI_HP\F\test\result_wind_10min";
% [Wind_Data] = read_wind_data(start_time, end_time, wind_dir);
 %% 9 calculate aerodynamic damping ratio

    ncycle = input.ncycle;
    t = input.t;
    p_filt_m = input.p_filt_m;
    x_k_k = input.x_k_k;
    nmodes = input.nmodes;
    Freq = input.Freq;
    fs = input.fs;
    nVIV = input.nVIV;
    VIV_mode_seq = input.VIV_mode_seq;

    p = inputParser;
    addParameter(p,'showtext',true,@islogical);
    addParameter(p,'showplot',true,@islogical);
    addParameter(p,'filterstyle','fft',@ischar);% fft use the function by myself, bandpass use the function in matlab. fft is much faster than bandpass but may not be accurate
    parse(p,varargin{:});
    showtext = p.Results.showtext;
    showplot = p.Results.showplot;
    filterstyle = p.Results.filterstyle;

    f_keep = [Freq(VIV_mode_seq) * 0.9, Freq(VIV_mode_seq) * 1.1];



    Fa = p_filt_m;
    dis = x_k_k(VIV_mode_seq, :);
    vel = x_k_k(nmodes + VIV_mode_seq, :);

    % dis = input.u;
    % vel = input.udot;


    

    switch filterstyle
        case 'fft'
        
        for k1 = 1:size(Fa, 1)
            Fa_filtered(k1,:) = fft_filter(fs, Fa(k1,:), f_keep);
        end
        for k1 = 1:size(Fa, 1)
            vel_filtered(k1,:) = fft_filter(fs, vel(k1,:), f_keep);
        end
        for k1 = 1:size(Fa, 1)
            dis_filtered(k1,:) = fft_filter(fs, dis(k1,:), f_keep);
        end
        
        case 'bandpass'
        Fa_filtered = bandpass(Fa', [min(Freq) * 0.9, max(Freq) * 1.1], fs)';
        vel_filtered = bandpass(vel', [min(Freq) * 0.9, max(Freq) * 1.1], fs)';
        dis_filtered = bandpass(dis', [min(Freq) * 0.9, max(Freq) * 1.1], fs)';


        case 'nofilter'
        Fa_filtered = Fa;
        vel_filtered = vel;
        dis_filtered = dis;

        otherwise
        error('Invalid filter style. Choose either "fft" or "bandpass".')
    end

    if showplot
        figure
        plot(t,Fa(1,:))
        hold on
        plot(t,Fa_filtered(1,:))
        xlabel('Time(s)')
        ylabel("Modal force")
        legend("Kalman Filter","Bandpass")
        [f2, magnitude2] = fft_transform(fs,Fa_filtered(1,:));
        % [f2, magnitude2] = fft_transform(fs,Fa_filtered(1,:));
        [f1, magnitude1] = fft_transform(fs,Fa(1,:));
        
        figure
        plot(f1, magnitude1)
        hold on
        plot(f2, magnitude2)
        % plot(f3, magnitude3)
        xlim([0, max(Freq) * 1.1])
        xlabel('Frequency(Hz)')
        ylabel("Value")
        legend("Kalman Filter","Bandpass")
    end

    if showplot
        figure
        plot(t,vel(1,:))
        hold on
        plot(t,vel_filtered(1,:))
        xlabel('Time(s)')
        ylabel("Velocity")
        legend("Kalman Filter","Bandpass")
        [f2, magnitude2] = fft_transform(fs,vel(1,:));
        % [f2, magnitude2] = fft_transform(fs,Fa_filtered(1,:));
        [f1, magnitude1] = fft_transform(fs,vel_filtered(1,:));
        
        figure
        plot(f1, magnitude1)
        hold on
        plot(f2, magnitude2)
        % plot(f3, magnitude3)
        xlim([0, max(Freq) * 1.1])
        xlabel('Frequency(Hz)')
        ylabel("Value")
        legend("Kalman Filter","Bandpass")
    end

    % 找到峰值设定保存频率成分的变量

    top_freqs = cell(1, nVIV);
    for k1 = 1:nVIV
        [top_freqs{k1}, ~, ~] = extractSignificantFrequencies(fs, vel_filtered(k1, :), 'showplot', showplot);
        % [top_freqs{k1}, ~, ~] = extractSignificantFrequencies(fs, Fa_filtered(k1, :), 'showplot', showplot);
        [top_freqs_vel{k1}, ~, ~] = extractSignificantFrequencies(fs, vel_filtered(k1, :),'showplot', showplot);
    end

    ifq_interpolated_allmodes = cell(1, nVIV);
    for k1=1:nVIV
        for k2 = 1:length(top_freqs{k1})
            f_keep_temp = [top_freqs{k1}(k2) * 0.9, top_freqs{k1}(k2) * 1.1];
            [ifq_interpolated_allmodes{k1}{k2}] = instfreq_samelength(fs, Fa(k1,:), f_keep_temp, t, 'showplot', showplot, 'showtext', showtext);
        end
    end

    bandwidth = 0.01; %根据需要调整带宽
    filtered_Fa = cell(1,nVIV);
    for k1 = 1:nVIV
        frequencies = top_freqs{k1};
        Fa_current = Fa_filtered(k1, :);
        filtered_Fa{k1} =filterSignals(Fa_current, frequencies, fs, bandwidth,'filterstyle','fft','showplot', showplot);
    end

    filtered_vel = cell(1,nVIV);
    for k1 = 1:nVIV
        frequencies = top_freqs{k1};
        vel_current = vel_filtered(k1, :);
        filtered_vel{k1} =filterSignals(vel_current, frequencies, fs, bandwidth,'filterstyle','fft','showplot', showplot);
    end

    filtered_dis = cell(1,nVIV);
    for k1 = 1:nVIV
        frequencies = top_freqs{k1};
        dis_current = dis_filtered(k1, :);
        filtered_dis{k1} =filterSignals(dis_current, frequencies, fs, bandwidth,'filterstyle','fft','showplot', showplot);
    end


    peaks_locs_cell = cell(nVIV, 1); % 初始化一个cell数组以保存peaks和locs的结构


    % 寻找不同模态不同频率力信号的周期
    for k1 = 1:nVIV
        % filtered_Fa_current = filtered_Fa{k1}; % 获取当前模态的滤波数据
        % filtered_dis_current = filtered_dis{k1}; % 获取当前模态的滤波数据
        filtered_vel_current = filtered_vel{k1}; % 获取当前模态的滤波数据
        
        peaks_locs_struct = struct(); % 初始化一个结构体以保存peaks和locs       
        for k2 = 1:length(top_freqs{k1})
            f_temp = top_freqs{k1}(k2);
            T_temp = 1/f_temp;
            d = T_temp * 50*0.9;
            pp = 0;
            [peaks, locs] = findpeaks(filtered_vel_current{k2}, ...
                                      'MinPeakDistance', d, ...
                                      'MinPeakProminence', pp);
            % [peaks, locs] = findpeaks(filtered_Fa_current{k2}, ...
            %                           'MinPeakDistance', d, ...
            %                           'MinPeakProminence', pp);
            % [peaks, locs] = findpeaks(filtered_dis_current{k2}, ...
            %               'MinPeakDistance', d, ...
            %               'MinPeakProminence', pp);
            % 将peaks和locs保存为结构体的字段
            peaks_locs_struct(k2).peaks = peaks;
            peaks_locs_struct(k2).locs = locs;
            
            if showplot
            figure
            plot(t,filtered_Fa_current{k2})
            % plot(t,filtered_dis_current{k2})
            hold on 
            scatter(t(locs),peaks)
            hold off
            str = sprintf('mode freq %.2f Hz',Freq(VIV_mode_seq(k1)));
            title(str)
            end
        end 
        peaks_locs_cell{k1} = peaks_locs_struct; % 将结构体保存在cell数组中
    end


    amp_cell = cell(size(filtered_Fa)); % 初始化用于保存amp数组的cell
    zeta_all_cell = cell(size(filtered_Fa)); % 初始化用于保存zeta_all数组的cell

    for i = 1:length(filtered_Fa)
        % i = 9
        force_mode_signals = filtered_Fa{i}; % 获取当前模态下的力信号
        vel_filtered_mode_signal = filtered_vel{i}; % 获取当前模式下的速度信号
        dis_filtered_mode_signal = filtered_dis{i}; % 获取当前模式下的位移信号
        ifq_interpolated_mode = ifq_interpolated_allmodes{i}; % 获取当前模式下的ifq_interpolated
        peaks_locs_cell_mode =peaks_locs_cell{i};

        amp_mode = {}; % 初始化用于保存当前模式amp的cell
        zeta_all_mode = {}; % 初始化用于保存当前模式zeta_all的cell
        
        for j = 1:length(force_mode_signals) % 遍历当前模式下的每个信号
            freq_temp = ifq_interpolated_mode{j};
            % 
            % Fa_temp = force_mode_signals{j};
            % vel_temp = vel_filtered_mode_signal{j};
            % dis_temp = dis_filtered_mode_signal{j};
            % 
            Fa_temp = Fa(i,:);
            vel_temp = vel(i,:);
            dis_temp = dis(i,:);

            locs = peaks_locs_cell_mode(j).locs;
            [result] = compute_dynamics_parameters(ncycle, t, Fa_temp, vel_temp, freq_temp, dis_temp , locs,'showplot', showplot);
            
          

            amp_mode{j} = result.amp; % 将此信号的amp数组保存到当前模式的cell中
            zeta_all_mode{j} = result.zeta_all; % 将此信号的zeta_all数组保存到当前模式的cell中
            t_cycle_mean_mode{j} = result.t_cycle_mean;
            work_mode{j} = result.work;
            % [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
            % scatter(result.amp*max(mode_deck(:,i)), result.zeta_all)
            % titlestr= sprintf('mode freq %.2f Hz, mode sel %d, freq %.2f Hz',Freq(i),modesel(i),top_freqs{i}(j));
            % title(titlestr)
            

        end
        
        amp_cell{i} = amp_mode; % 将当前模式的amp cell保存到总cell中
        zeta_all_cell{i} = zeta_all_mode; % 将当前模式的zeta_all cell保存到总cell中
        t_cycle_mean_cell{i}=t_cycle_mean_mode;
        work_cell{i}=work_mode;
    end



    result_Main.amp_cell = amp_cell;
    result_Main.zeta_all_cell = zeta_all_cell;
    result_Main.top_freqs = top_freqs;
    result_Main.top_freqs_vel = top_freqs_vel;
    result_Main.peaks_locs_cell = peaks_locs_cell;
    result_Main.filtered_Fa = filtered_Fa;
    result_Main.filtered_vel = filtered_vel;
    result_Main.filtered_dis = filtered_dis;
    result_Main.ifq_interpolated_allmodes = ifq_interpolated_allmodes;
    result_Main.t_cycle_mean_cell=t_cycle_mean_cell;
    result_Main.work_cell=work_cell;
    
end