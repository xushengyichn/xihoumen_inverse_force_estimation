clc
clear
close all

% n1 = 100;
% n2 = 100;
% n3 = 50;
% 
% 
% lambdas_m_list = linspace(100,200, n1);
% sigma_ps_m_list = linspace(200,300,n2);
% omega_0_list = linspace(400,500,n3);
% 
% [X, Y, W] = meshgrid(lambdas_m_list, sigma_ps_m_list, omega_0_list);
% combinations = [reshape(X, [], 1), reshape(Y, [], 1), reshape(W, [], 1)];
% numIterations = size(combinations,1);
% 
% for k1 = 1:numIterations
%     lambda_m = combinations(k1,1);
%     sigma_ps_m = combinations(k1,2);
%     omega_0 = combinations(k1,3);
% 
%     logL(k1) = lambda_m+sigma_ps_m+omega_0;
% end
% 
% 
% Z_1 = reshape(logL, n1, n2,n3);
% 
% figure
% scatter3(combinations(:,1), combinations(:,2), combinations(:,3), 10, logL, 'filled')


% fs = 50;
% t = 0:1/fs:100;
% f =0.5;
% phaselag = 2.8-pi
% x1 = cos(2*pi*f*t - phaselag);
% x2 = cos(2*pi*f*t )*100;
% 
% 
% 
% [phaseDiff] = fft_phase_lag(fs,x1,x2,f)

figure
% x = [0 0.001 0.005 0.01 0.05 1 2 3 4 5];
% Nx = length(x);
% clim = [min(x) max(x)];
% dx = min(diff(x));
% y = clim(1):dx:clim(2);
% 
% cmap = colormap(jet(Nx));
% cmap2 = [...
% interp1(x(:),cmap(:,1),y(:)) ...
% interp1(x(:),cmap(:,2),y(:)) ...
% interp1(x(:),cmap(:,3),y(:)) ...
% ];
S = peaks(30);
surf(abs(S))


colorbar



figure
x = [0 0.00005 0.0001 0.0002 0.0003 0.005 10];
Nx = length(x);
c_lim = [min(x) max(x)];
dx = min(diff(x));
y = c_lim(1):dx:c_lim(2);
for k=1:Nx-1, y(y>x(k) & y<=x(k+1)) = x(k+1); end % NEW
cmap = colormap(parula(Nx));
cmap2 = [...
interp1(x(:),cmap(:,1),y(:)) ...
interp1(x(:),cmap(:,2),y(:)) ...
interp1(x(:),cmap(:,3),y(:)) ...
];
% S = peaks(30);
surf(abs(S))
colormap(cmap2)
clim(c_lim)
colorbar
