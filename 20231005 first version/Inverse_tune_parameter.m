clc; clear; close all;
addpath(genpath("F:\git\ssm_tools\"))
addpath(genpath("F:\git\Function_shengyi_package\"))
addpath(genpath("F:\git\xihoumen_inverse_force_estimation\FEM_model\"))
subStreamNumberDefault = 2132;

params = Init_fun();
input.ON = params.ON;
input.OFF = params.OFF;
%% 0 绘图参数
fig_bool = params.ON;
input.fig_bool = params.OFF;
input.num_figs_in_row = 12; %每一行显示几个图
input.figPos = params.figPosSmall; %图的大小，参数基于InitScript.m中的设置
%设置图片间隔
input.gap_between_images = [0, 0];
input.figureIdx = 0;


input.omega_0_variation =1;
n = 4;
[result] = viv2013(n, params.OFF);
input.start_time = result.startDate;
input.end_time = result.endDate;
% start_time = datetime('2013-04-03 16:24:00', 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
% end_time = datetime('2013-04-08 15:30:59', 'InputFormat', 'yyyy-MM-dd HH:mm:ss');
input.wind_dir = "F:\test\result_wind_10min";
input.acc_dir = "F:\test\result";


figureIdx = 0;
num_figs_in_row = 12; %每一行显示几个图
figPos = params.figPosSmall; %图的大小，参数基于InitScript.m中的设置
gap_between_images = [0, 0];
%% 1 set hyperparameters
n1 = 10;
n2 = 10;

lambdas_m_list = logspace(-15,-1, n1);
sigma_ps_m_list = linspace(1e2,1e5,n2);
[X, Y] = meshgrid(lambdas_m_list, sigma_ps_m_list);
combinations = [reshape(X, [], 1), reshape(Y, [], 1)];

numIterations = size(combinations,1);

if isempty(gcp('nocreate'))
    parpool();
end

b = ProgressBar(numIterations, ...
    'IsParallel', true, ...
    'WorkerDirectory', pwd(), ...
    'Title', 'Parallel 2' ...
    );
b.setup([], [], []);

parfor k1 = 1:numIterations
    localInput = input;  % Make a local copy
    combinations_temp = combinations(k1,:);
    localInput.lambada = combinations_temp(1);
    localInput.sigma_p = combinations_temp(2);

    [result_Main] = Main(localInput,'showtext',false);
    logL(k1)=result_Main.logL;
    logSk(k1) = result_Main.logSk;
    logek(k1)=result_Main.logek;
    updateParallel([], pwd);
end


b.release();

Z_1 = reshape(logL, n2, n1);
Z_2 = reshape(logSk, n2, n1);
Z_3 = reshape(logek, n2, n1);

save temp

%% plot
if fig_bool == params.ON


    % x = [-8*10^(10),  -3*10^(10), -2*10^(10), -1.5*10^(10), -1.46*10^(10),-1.449*10^(10),-1.448*10^(10),-1.447*10^(10),-1.446*10^(10),-1.445*10^(10),-1.444*10^(10),-1.443*10^(10),-1.442*10^(10),-1.441*10^(10),-1.44*10^(10),-1.43*10^(10) ];
    % Identify max value from data
    max_val = max(logL);
    
    % Set up the intervals
    % I've taken the maximum value from your data and set up intervals that focus granularity around it
    x = [min(logL), linspace(-1.5*10^(10), max_val, 100)];
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

    [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);

    colormap(cmap2)
    contourf(X, Y, Z_1);  % 绘制等高线图
    set(gca, 'XScale', 'log');
    xlabel('lambdas');
    ylabel('sigma_ps');
    colorbar;  % 添加颜色栏
    title('logL');

[figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
    contourf(X, Y, Z_2);  % 绘制等高线图
    set(gca, 'XScale', 'log');
    xlabel('lambdas');
    ylabel('sigma_ps');
    colorbar;  % 添加颜色栏
    title('logSk');


    [figureIdx, figPos_temp, hFigure] = create_figure(figureIdx, num_figs_in_row, figPos, gap_between_images);
    contourf(X, Y, Z_3);  % 绘制等高线图
    set(gca, 'XScale', 'log');
    xlabel('lambdas');
    ylabel('sigma_ps');
    colorbar;  % 添加颜色栏
    title('logek');
end