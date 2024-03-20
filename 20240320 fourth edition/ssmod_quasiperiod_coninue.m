    function [F_c, L_c, H_c, sigma_w12] = ssmod_quasiperiod_coninue(lambdas, sigma_ps, omega_0, np,VIV_mode_seq)
        % create arrays for latent force model
        F_c_array = zeros(2, 2, np);
        L_c_array = zeros(2, 2, np);
        H_c_array = zeros(1, 2, np);
        sigma_w = zeros(np, 1);
        
        for k1 = 1:np
            if k1==VIV_mode_seq
                [F_c_array(:, :, k1), L_c_array(:, :, k1), H_c_array(:, :, k1), sigma_w(k1)] = ssmod_quasiperiod(lambdas(k1), sigma_ps(k1), omega_0(k1));
            else
                F_c_array(:, :, k1)=zeros(2, 2);
                L_c_array(:, :, k1)=zeros(2, 2);
                H_c_array(:, :, k1)=[0, 0];
                sigma_w(k1) = 0;
            end
        end

        F_c = [];
        L_c = [];
        H_c = [];
        
        nVIV = length(VIV_mode_seq);
        H_c = zeros(nVIV, nVIV*2); % 2 is the kernel size for quasi periodic
        F_c = zeros(nVIV*2, nVIV*2);
        L_c = zeros(nVIV*2, nVIV*2);

        for k1 = 1:np
            idx = find(VIV_mode_seq == k1);
            if ~isempty(idx)
                F_c(idx*2-1:idx*2,idx*2-1:idx*2) = F_c_array(:, :, k1);
                L_c(idx*2-1:idx*2,idx*2-1:idx*2) = L_c_array(:, :, k1);
                H_c(idx,idx*2-1:idx*2) = H_c_array(:, :, k1);
            end
        end

        sigma_w12 = zeros(nVIV*2, nVIV*2);

        for k1 = 1:np
            idx = find(VIV_mode_seq == k1);
            if ~isempty(idx)
                sigma_w12(idx*2-1:idx*2,idx*2-1:idx*2) = sigma_w(k1) * eye(2);
            end
            
        end

    end