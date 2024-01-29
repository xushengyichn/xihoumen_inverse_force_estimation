clc;clear;close all
skipClear =true;


VIV_sels = [2,3,4,5,6,7,8,9,10,11,12,16,17,18,19,22];
% VIV_sels = [7,8,9,10,11,12,16,17,18,19,22];

for k1 = 1:length(VIV_sels)
    VIV_sel = VIV_sels(k1);
    fig_bool = false;
    modelupdate = false;
    run("Main.m");
    clearvars('-except', 'VIV_sels','skipClear','fig_bool');
end


for k1 = 1:length(VIV_sels)
    VIV_sel = VIV_sels(k1);
    fig_bool = false;
    modelupdate = true;
    run("Main.m");
    clearvars('-except', 'VIV_sels','skipClear','fig_bool');
end