function [P]=quasiperiod_force(lambda,sigma_p,f,dt,t)
    f1=f;
    omega_0=2*pi*f1;
    [F_c,L_c,H_c,sigma_w] = ssmod_quasiperiod(lambda, sigma_p, omega_0);
    sigma_w_2 = sigma_w ^ 2;
    F_d = expm(F_c * dt);
    H_d = H_c;
    Qwc = diag([sigma_w_2]);
    Qwd = dt * L_c * Qwc * L_c';
    N = length(t);
    s = zeros(2, N);
    s(:,1) = zeros(2, 1);
    w = randn(1, N) .* sqrt(diag(Qwd));
    for k1 = 1:N
        s(:,k1 + 1) = F_d * s(:,k1) +w(:,k1);
        P(:,k1) = H_d * s(:,k1);
    end
    % P=rescale(P)-0.5;
end