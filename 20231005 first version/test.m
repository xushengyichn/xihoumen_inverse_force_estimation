clc
clear
close all

lambda=1e-1;

t=0:0.01:100;

y = exp(-1*lambda*abs(t));

figure;
plot(t,y)

hold on

lambda=1e-2;

t=0:0.01:100;

y = exp(-1*lambda*abs(t));

plot(t,y)

