function [uu, vv, ww, average_wind] = cal_windspeed(data, samp_freq, tag)
% tag: UA1 3 5 --> 0 UA2 4 6 -->1
average_wind = [];
uu = [];
vv = [];
ww = [];

slice_len = samp_freq*60*10;
%data(:, 2) = -data(:, 2); % Modify the y direction
data(:, 1) = -data(:, 1); % Modify the y direction
len = size(data, 1)/slice_len;
for j = 1:len
    x = data((j-1)*slice_len+1:j*slice_len, 1);  
    y = data((j-1)*slice_len+1:j*slice_len, 2);
    z = data((j-1)*slice_len+1:j*slice_len, 3); 
     x(x > 60 | x < -60) = 0;
     y(y > 60 | y < -60) = 0;
     z(z > 60 | z < -60) = 0;
    
    x_mean = mean(x);
    y_mean = mean(y);


    U = sqrt(x_mean^2 + y_mean^2);
    
    if y_mean > 0
        beta_mean = acosd(x_mean/U);
    else
        beta_mean = 360 - acosd(x_mean/U);
    end
     u = x*cos(beta_mean) + y*sin(beta_mean);
     v = -x*sin(beta_mean) + y*cos(beta_mean);
   % u=x*cos(pi/4)-y*cos(pi/4);
   % v=x*cos(pi/4)+y*cos(pi/4);
    w = z;
   % if beta_mean < 225 && beta_mean > 45
   %    U=U*cos((135-beta_mean)/180*pi);
   %  else
   %     U=U*cos((beta_mean-315)/180*pi);
   %  end

  

    alpha_mean = atand(mean(w)/U);
    

    
    
%% data quality control
%if ~tag
 %       if beta_mean < 60 || beta_mean > 220 || U < 0
  %          continue
  %      end
  %  else
  %     if beta_mean > 30 && beta_mean < 240 || U < 0
   %        continue
   %     end
%end
    
%     if ~tag
%         if U < 8
%             continue
%         end
%     else
%         if U < 8
%             continue
%         end
%     end
    
    
  %  if locs(2) >= log10(2) || locs(2) >= log10(2)
   %     continue
  %  end
    
  
    
   
%%
    clear Su Fu  x y z   lg_su lg_fu
    uu = [uu ,u];
    vv = [vv, v];
    ww = [ww, w];
    average_wind = [average_wind, [U; beta_mean; alpha_mean]];
    
    
end



