function [xn, yn, xn_true] = CalResponse(A_d, B_d, G_d, J_d, p, Q, R, N, x0, ns, n_sensors)
    w = sqrt(Q) * randn(ns, N);
    v = sqrt(R) * randn(n_sensors, N);
    xn = x0;
    xn_true = x0;

    for t1 = 1:N - 1
        xn(:, t1 + 1) = A_d * xn(:, t1) + B_d * p(:, t1) + w(:, t1);
        yn(:, t1) = G_d * xn(:, t1) + J_d * p(:, t1) + v(:, t1);
        xn_true(:, t1 + 1) = A_d * xn_true(:, t1) + B_d * p(:, t1);
    end

    yn(:, end + 1) = G_d * xn(:, end) + J_d * p(:, end) + v(:, end);
end