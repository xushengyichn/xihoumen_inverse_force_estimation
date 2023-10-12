    function [F_c, L_c, H_c, sigma_w12] = ssmod_quasiperiod_coninue(lambdas, sigma_ps, omega_0, np)
        % create arrays for latent force model
        F_c_array = zeros(2, 2, np);
        L_c_array = zeros(2, 2, np);
        H_c_array = zeros(1, 2, np);
        sigma_w = zeros(np, 1);

        for k1 = 1:np
            [F_c_array(:, :, k1), L_c_array(:, :, k1), H_c_array(:, :, k1), sigma_w(k1)] = ssmod_quasiperiod(lambdas(k1), sigma_ps(k1), omega_0(k1));
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