%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%Date: 2023-01-25 10:57:55
%LastEditors: xushengyichn 54436848+xushengyichn@users.noreply.github.com
%LastEditTime: 2023-01-25 13:47:58
%FilePath: \20230124优化问题\a_3_VEF_parameter.m
%Description: 获取气动力参数（改代码仅用于双幅桥案例）以及实验数据
%
%Copyright (c) 2023 by xushengyichn 54436848+xushengyichn@users.noreply.github.com, All Rights Reserved.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [output] = a_3_VEF_parameter(input)
    ExpName = input.ExpName; % 试验名称
    girderindex = input.girderindex; % 设置读取某一幅桥的数据

    % 试验数据
      
    my_table = load('SZTD110_testfile.mat');
    my_table = my_table.my_table;
    fname = ExpName;
    chanum = 8; %记录通道数
    filename = strsplit(fname, '-');
    casenumber = cell2mat(filename(3));
    spacing = str2double(cell2mat(filename(4)));
    type = cell2mat(filename(5));
    Rspeed = cell2mat(filename(6));

    casenumber = str2double(casenumber(5:end));
Rspeed = Rspeed(1:end - 1);
Rspeed = str2double(Rspeed);
    if strcmp(type, 'fasan')
        type = 1;
    else

        if strcmp(type, 'shuaijian')
            type = 2;
        else
            error("实验数据文件名可能有误，请确认信号状态为衰减或发散，程序终止")
        end

    end

      % 判断工况风攻角
      if casenumber <= 17 || and(casenumber >= 32, casenumber <= 35) || and(casenumber >= 40, casenumber <= 16) || casenumber == 55
        AOA = 3;
    else

        if and(casenumber >= 18, casenumber <= 25) || and(casenumber >= 36, casenumber <= 37) || and(casenumber >= 47, casenumber <= 50) || casenumber == 54
            AOA = 0;
        else

            if and(casenumber >= 26, casenumber <= 31) || and(casenumber >= 38, casenumber <= 39) || and(casenumber >= 51, casenumber <= 53) || and(casenumber >= 56, casenumber <= 57)
                AOA = -3;
            end

        end

    end

    if girderindex == 1
        girder = 'up';
        %导入试验数据
        Name = my_table.casename;
        isexist = find(Name == fname);
        a1 = my_table.up_parameter_a1(isexist);
        a2 = my_table.up_parameter_a2(isexist);
        a3 = my_table.up_parameter_a3(isexist);
        a4 = my_table.up_parameter_a4(isexist);
        a5 = my_table.up_parameter_a5(isexist);
        a = [a1 a2 a3 a4 a5];
        m_exp = my_table.up_m(isexist); %试验中单位长度模型的质量
        H4 = my_table.up_parameter_H4(isexist); % 气动刚度
        Fren_vibration_withwind = my_table.up_Fren_vibration_withwind(isexist);
        F0 = my_table.up_Fre_vibration(isexist); % Frequency without wind
        Zeta0 = my_table.up_dltx_zeta0(isexist); % damping ratio without wind
        ReducedFrequency = my_table.up_ReducedFre(isexist); % 折减频率
    else
        girder = 'down';
        %导入试验数据
        Name = my_table.casename;
        isexist = find(Name == fname);
        a1 = my_table.down_parameter_a1(isexist);
        a2 = my_table.down_parameter_a2(isexist);
        a3 = my_table.down_parameter_a3(isexist);
        a4 = my_table.down_parameter_a4(isexist);
        a5 = my_table.down_parameter_a5(isexist);
        a = [a1 a2 a3 a4 a5];
        m_exp = my_table.down_m(isexist); %试验中单位长度模型的质量
        H4 = my_table.down_parameter_H4(isexist); % 气动刚度
        Fren_vibration_withwind = my_table.down_Fren_vibration_withwind(isexist);
        F0 = my_table.down_Fre_vibration(isexist); % Frequency without wind
        Zeta0 = my_table.down_dltx_zeta0(isexist); % damping ratio without wind
        ReducedFrequency = my_table.down_ReducedFre(isexist); % 折减频率
    end

    output.a = a;
    output.m_exp = m_exp;
    output.H4 = H4;
    output.Fren_vibration_withwind = Fren_vibration_withwind;
    output.F0 = F0;
    output.Zeta0 = Zeta0;
    output.ReducedFrequency = ReducedFrequency;

end