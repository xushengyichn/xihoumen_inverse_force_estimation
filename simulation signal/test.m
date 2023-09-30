clc
clear
close all

n1 = 100;
n2 = 100;
n3 = 50;


lambdas_m_list = linspace(100,200, n1);
sigma_ps_m_list = linspace(200,300,n2);
omega_0_list = linspace(400,500,n3);

[X, Y, W] = meshgrid(lambdas_m_list, sigma_ps_m_list, omega_0_list);
combinations = [reshape(X, [], 1), reshape(Y, [], 1), reshape(W, [], 1)];
numIterations = size(combinations,1);

for k1 = 1:numIterations
    lambda_m = combinations(k1,1);
    sigma_ps_m = combinations(k1,2);
    omega_0 = combinations(k1,3);
    
    logL(k1) = lambda_m+sigma_ps_m+omega_0;
end


Z_1 = reshape(logL, n1, n2,n3);

figure
scatter3(combinations(:,1), combinations(:,2), combinations(:,3), 10, logL, 'filled')