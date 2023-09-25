function [uu, vv, ww, average_wind] = cal_windspeed(data, samp_freq, tag)
% tag: UA1 3 5 --> 0 UA2 4 6 -->1
average_wind = [];
uu = [];
vv = [];
ww = [];

slice_len = samp_freq*60*10;
data(:, 2) = -data(:, 2);  %Modify the y direction

len = size(data, 1)/slice_len;
for j = 1:len
    x = data((j-1)*slice_len+1:j*slice_len, 1);  
    y = data((j-1)*slice_len+1:j*slice_len, 2);
    z = data((j-1)*slice_len+1:j*slice_len, 3); 
    x(x > 40 | x < -40) = 0;
    y(y > 40 | y < -40) = 0;
    z(z > 40 | z < -40) = 0;
    
    x_mean = mean(x);
    y_mean = mean(y);
    U = sqrt(x_mean^2 + y_mean^2);
    if U == 0
        continue
    end
    if y_mean > 0
        beta_mean = acosd(x_mean/U);
    else
        beta_mean = 360 - acosd(x_mean/U);
    end
    u = x*cos(beta_mean) + y*sin(beta_mean);
    v = -x*sin(beta_mean) + y*cos(beta_mean);
    w = z;
    alpha_mean = atand(mean(w)/U);
    

    
    sigmau=std(u);
    a=0;
    
  for i=1:19200
    if u(i,1)-U>5*sigmau
       a=a+1;

    else
        a=a+0;
    end 
  end
  
  if a>0
            continue
        end
  
%% data quality control
if ~tag
        if  U < 0
    
            continue
        end
    else
        if  U < 0
            continue
        end
end
    
%     if ~tag
%         if U < 8
%             continue
%         end
%     else
%         if U < 8
%             continue
%         end
%     end
    [Su, Fu] = psd_pwelch(u, 32);
    Su = Fu.*Su/(U^2);
    lg_su = log10(Su);
    lg_fu = log10(Fu);
    
 


    % peak values delete
    smooth_su = smooth(lg_su, 0.1);
    [~, locs] = findpeaks(smooth_su, lg_fu, 'NPeaks', 2, 'SortStr', 'descend');
    
   
    if locs(2) >= log10(2) || locs(2) >= log10(2)
     continue
  end
    
  
    
   
%%
    clear Su Fu  x y z   lg_su lg_fu
    uu = [uu ,u];
    vv = [vv, v];
    ww = [ww, w];
    average_wind = [average_wind, [U; beta_mean; alpha_mean]];
    
    
end



