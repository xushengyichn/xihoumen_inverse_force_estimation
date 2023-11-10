computer_name = getenv('COMPUTERNAME');
if strcmp(computer_name,'SHENGYI_HP')
    addpath(genpath("F:\git\ssm_tools\"))
    addpath(genpath("F:\git\Function_shengyi_package\"))
    addpath(genpath("F:\git\xihoumen_inverse_force_estimation\FEM_model\"))
    addpath(genpath("F:\git\HHT-Tutorial\"))
elseif strcmp(computer_name,'mac')
    addpath(genpath("/Users/xushengyi/Documents/GitHub/Function_shengyi_package"))
    addpath(genpath("/Users/xushengyi/Documents/GitHub/ssm_tools_sy"))
    addpath(genpath("/Users/xushengyi/Documents/GitHub/xihoumen_inverse_force_estimation/FEM_model"))
    addpath(genpath("/Users/xushengyi/Documents/GitHub/HHT-Tutorial"))
elseif strcmp(computer_name,'ROG-SHENGYI')
    addpath(genpath("D:\Users\xushe\Documents\GitHub\Function_shengyi_package"))
    addpath(genpath("D:\Users\xushe\Documents\GitHub\ssm_tools"))
    addpath(genpath("D:\git\xihoumen_inverse_force_estimation\FEM_model"))
    addpath(genpath("D:\Users\xushe\Documents\GitHub\xihoumen_data_extract"))
    addpath(genpath("D:\Users\xushe\Documents\GitHub\HHT-Tutorial\"))
elseif strcmp(computer_name,'NTNU08916')
    addpath(genpath("C:\Users\shengyix\Documents\GitHub\Function_shengyi_package"))
    addpath(genpath("C:\Users\shengyix\Documents\GitHub\ssm_tools_sy"))
    addpath(genpath("C:\Users\shengyix\Documents\GitHub\xihoumen_inverse_force_estimation\FEM_model"))
    addpath(genpath("C:\Users\shengyix\Documents\GitHub\xihoumen_data_extract"))
else
    error("Please add path first.")
end


subStreamNumberDefault = 2132;

run("InitScript.m")

%% 0 绘图参数
fig_bool = ON;
num_figs_in_row = 6; %每一行显示几个图
figPos = figPosSmall; %图的大小，参数基于InitScript.m中的设置
%设置图片间隔
gap_between_images = [0, 0];
figureIdx = 0;


computer_name = getenv('COMPUTERNAME');
if strcmp(computer_name,'SHENGYI_HP')
    input.acc_dir = "F:\test\result";
    input.wind_dir = "F:\test\result_wind_10min";
elseif strcmp(computer_name,'mac')
    input.acc_dir = "/Users/xushengyi/Documents/xihoumendata/acc";
    input.wind_dir = "/Users/xushengyi/Documents/xihoumendata/wind";
elseif strcmp(computer_name,'ROG-SHENGYI')
    input.acc_dir = "D:\xihoumendata\acc";
    input.wind_dir = "D:\xihoumendata\wind";
elseif strcmp(computer_name,'ketizu')
    input.acc_dir = "Z:\Drive\Backup\SHENGYI_HP\F\test\result";
    input.wind_dir = "Z:\Drive\Backup\SHENGYI_HP\F\test\result_wind_10min";
elseif strcmp(computer_name,'NTNU08916')
     input.acc_dir = "C:\Users\shengyix\Documents\xihoumendata\acc";
    input.wind_dir = "C:\Users\shengyix\Documents\xihoumendata\wind";
else
    error("Please add data folder first.")
end