%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Author: xushengyichn xushengyichn@outlook.com
%Date: 2023-12-03 22:49:01
%LastEditors: xushengyichn xushengyichn@outlook.com
%LastEditTime: 2023-12-03 23:03:32
%FilePath: /ssm_tools_sy/Users/xushengyi/Documents/GitHub/xihoumen_inverse_force_estimation/20231203 third edition/ssmod_matern_coninue.m
%Description: 
%
%Copyright (c) 2023 by ${git_name_email}, All Rights Reserved. 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function [F_c, L_c, H_c, sigma_w12] = ssmod_matern_coninue(lambdas, sigma_ps, np,VIV_mode_seq)
        % create arrays for latent force model
        F_c_array = zeros(1, 1, np);
        L_c_array = zeros(1, 1, np);
        H_c_array = zeros(1, 1, np);
        sigma_w = zeros(np, 1);
        
        for k1 = 1:np
            if k1==VIV_mode_seq
                
                F_c_array(:, :, k1)=zeros(1, 1);
                L_c_array(:, :, k1)=zeros(1, 1);
                H_c_array(:, :, k1)=[0];
                sigma_w(k1) = 0;
            else
                [F_c_array(:, :, k1), L_c_array(:, :, k1), H_c_array(:, :, k1), sigma_w(k1)] = ssmod_matern(lambdas(k1),0,sigma_ps(k1));
            end
        end

        F_c = [];
        L_c = [];
        H_c = [];

        for k1 = 1:np
            F_c = blkdiag(F_c, F_c_array(:, :, k1));
            L_c = blkdiag(L_c, L_c_array(:, :, k1));
            H_c = blkdiag(H_c, H_c_array(:, :, k1));
        end

        sigma_w12 = [];

        for k1 = 1:np
            sigma_w12 = blkdiag(sigma_w12, sigma_w(k1) * eye(2));
        end

    end