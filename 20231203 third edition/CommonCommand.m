% 获取操作系统类型
[osType, ~] = computer;
switch osType
    case 'PCWIN' % 对于 Windows 32位系统
        % 执行 Windows 专用命令
        computer_name = getenv('COMPUTERNAME');
    case 'PCWIN64' % 对于 Windows 64位系统
        % 执行 Windows 专用命令
        computer_name = getenv('COMPUTERNAME');
    case 'GLNXA64' % 对于 Linux 系统
        % 执行 Linux 专用命令
        [status, result] = system('hostname');
    case 'MACA64' % 对于 Mac 系统
        % 执行 Mac 专用命令
        computer_name = getenv('LOGNAME');% mac获取登录名
    otherwise
        % 对于其他或未知系统
        disp('Unknown operating system');
        result = '';
end


if strcmp(computer_name,'SHENGYI_HP')
    addpath(genpath("F:\git\ssm_tools\"))
    addpath(genpath("F:\git\Function_shengyi_package\"))
    addpath(genpath("F:\git\xihoumen_inverse_force_estimation\FEM_model\"))
    addpath(genpath("F:\git\xihoumen_data_extract\"))
    addpath(genpath("F:\git\HHT-Tutorial\"))
elseif strcmp(computer_name,'xushengyi')
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
elseif strcmp(computer_name,'WIN-JFOFTCAS8GU')
    addpath(genpath("D:\XSY\ForceEstimation"))
elseif strcmp(computer_name,'WIN-IUMUERP66UT')
    addpath(genpath("C:\Users\xushengyi\Documents\Github"))
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


[osType, ~] = computer;
switch osType
    case 'PCWIN' % 对于 Windows 32位系统
        % 执行 Windows 专用命令
        computer_name = getenv('COMPUTERNAME');
    case 'PCWIN64' % 对于 Windows 64位系统
        % 执行 Windows 专用命令
        computer_name = getenv('COMPUTERNAME');
    case 'GLNXA64' % 对于 Linux 系统
        % 执行 Linux 专用命令
        [status, result] = system('hostname');
    case 'MACA64' % 对于 Mac 系统
        % 执行 Mac 专用命令
        computer_name = getenv('LOGNAME');% mac获取登录名
    otherwise
        % 对于其他或未知系统
        disp('Unknown operating system');
        result = '';
end
if strcmp(computer_name,'SHENGYI_HP')
    input.acc_dir = "F:\test\result";
    input.wind_dir = "F:\test\result_wind_10min";
    input.wind_dir_all = "F:\test\result_wind";
elseif strcmp(computer_name,'xushengyi')
    input.acc_dir = "/Users/xushengyi/Documents/xihoumendata/acc";
    input.wind_dir = "/Users/xushengyi/Documents/xihoumendata/wind";
elseif strcmp(computer_name,'ROG-SHENGYI')
    input.acc_dir = "D:\xihoumendata\acc";
    input.wind_dir = "D:\xihoumendata\wind";
elseif strcmp(computer_name,'WIN-IUMUERP66UT')
    input.acc_dir = "Z:\Drive\Backup\SHENGYI_HP\F\test\result";
    input.wind_dir = "Z:\Drive\Backup\SHENGYI_HP\F\test\result_wind_10min";
    input.wind_dir_all = "Z:\Drive\Backup\SHENGYI_HP\F\test\result_wind";
elseif strcmp(computer_name,'NTNU08916')
    input.acc_dir = "C:\Users\shengyix\Documents\xihoumendata\acc";
    input.wind_dir = "C:\Users\shengyix\Documents\xihoumendata\wind";
elseif strcmp(computer_name,'WIN-JFOFTCAS8GU')
    input.acc_dir = "Z:\Drive\Backup\SHENGYI_HP\F\test\result";
    input.wind_dir = "Z:\Drive\Backup\SHENGYI_HP\F\test\result_wind_10min";
else
    error("Please add data folder first.")
end