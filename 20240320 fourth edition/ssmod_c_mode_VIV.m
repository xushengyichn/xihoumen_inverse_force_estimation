%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: ShengyiXu xushengyichn@outlook.com
%Date: 2023-12-07 16:48:37
%LastEditors: ShengyiXu xushengyichn@outlook.com
%LastEditTime: 2023-12-07 16:48:54
%FilePath: \Exercises-for-Techniques-for-estimation-in-dynamics-systemsf:\git\xihoumen_inverse_force_estimation\20231203 third edition\ssmod_c_mode_VIV.m
%Description: 
%
%Copyright (c) 2023 by ${git_name_email}, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [B_c_VIV, J_c_VIV] = ssmod_c_mode_VIV(nmodes,  phi, S_a, VIV_mode_seq)
        % A_c = [zeros(nmodes), eye(nmodes); ...
                   % -omega2, -Gamma];
        B_c_VIV = zeros(nmodes*2,length(VIV_mode_seq));
        for k1 = 1:length(VIV_mode_seq)
            B_c_VIV(nmodes+VIV_mode_seq(k1),k1)=1;
        end
        % G_c = [S_d * phi - S_a * phi * omega2, ...
                   % S_v * phi - S_a * phi * Gamma];
        J_c_VIV = [S_a * phi(:,VIV_mode_seq)];
    end