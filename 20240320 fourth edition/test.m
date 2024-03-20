%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: ShengyiXu xushengyichn@outlook.com
%Date: 2023-10-18 12:22:32
%LastEditors: ShengyiXu xushengyichn@outlook.com
%LastEditTime: 2023-10-23 11:39:53
%FilePath: \Exercises-for-Techniques-for-estimation-in-dynamics-systemsf:\git\xihoumen_inverse_force_estimation\20231005 first version\test.m
%Description: 
%
%Copyright (c) 2023 by ${git_name_email}, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
clear
close all


MM = 1;
f =1;
omega = 2*pi*f;
KK = omega^2*MM;
zeta = 0.3/100;
CC = 2*zeta*omega*MM;
fs = 1/0.01;
t = 0:1/fs:100;
% no force
p = 0*t;
[u udot u2dot] = NewmarkInt(t, MM, CC, KK, p, 1/2, 1/4, 1, 0);

[peaks, locs] = findpeaks(u);


amp_temp = peaks;
t_cycle_mean_temp = t(locs);
m_cycle=2;

for k1=1:length(t_cycle_mean_temp)-2
   
    deltam = log(amp_temp(k1+2)/amp_temp(k1));
    
    zetam(k1) = sqrt(deltam^2 / (4*m_cycle^2*pi^2 + deltam^2));
    
    if deltam > 0
        zetam(k1) = abs(zetam(k1));
    else
        zetam(k1) = -abs(zetam(k1));
    end
end
zetam =[zetam,zetam(end-1),zetam(end)]
% for k1=1:length(t_cycle_mean_temp)
%     if k1 <= m_cycle/2 % Beginning boundary
%         start_idx = 1;
%         end_idx = start_idx + m_cycle;
%     elseif k1 > length(t_cycle_mean_temp) - m_cycle/2 % Ending boundary
%         end_idx = length(t_cycle_mean_temp);
%         start_idx = end_idx - m_cycle;
%     else % Middle
%         start_idx = k1 - floor(m_cycle/2);
%         end_idx = k1 + floor(m_cycle/2);
%     end
% 
%     deltam = log(amp_temp(start_idx)/amp_temp(end_idx));
% 
%     zetam(k1) = sqrt(deltam^2 / (4*m_cycle*pi^2 + deltam^2));
% 
%     if deltam > 0
%         zetam(k1) = abs(zetam(k1));
%     else
%         zetam(k1) = -abs(zetam(k1));
%     end
% end

figure 
plot(t, u)
hold on
scatter(t(locs),peaks)

figure 
plot(t_cycle_mean_temp, zetam)
ylim([-0.5,0]/100)


%%
clc;clear;close all
test1 = importdata("temp.mat");
test2 = importdata("temp2.mat");

t = test1.t;
f1 = test1.F_filter;
f2 = test2.F_filter;
figure
plot(t,f1)
hold on 
plot(t,f2)
legend("temp1","temp2")

