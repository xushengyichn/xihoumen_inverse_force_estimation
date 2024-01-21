skipClear =true;


VIV_sels = [3,4];

for k1 = 1:length(VIV_sels)
    VIV_sel = VIV_sels(k1);
    fig_bool = false;
    run("Main.m");
    clearvars('-except', 'VIV_sels','skipClear','fig_bool');
end
