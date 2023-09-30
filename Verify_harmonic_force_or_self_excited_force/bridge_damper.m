function [g, ks] = bridge_damper(u, udot, b1, b2, b3, b4, b5, MM, CC, KK, mode_integral_2, mode_integral_3, mode_integral_4, mode_integral_5, mode_integral_6, gamma, beta, h, matrixsize, mode_number)

    for k01 = 1:matrixsize
        g(k01) = 0;

        if k01 == mode_number

            for k02 = 1:matrixsize
                g(k01) = g(k01) + CC(k01, k02) * udot(k02) + KK(k01, k02) * u(k02);
            end

            g(k01) = g(k01) - (b1 * mode_integral_2 + b2 * abs(u(k01)) * mode_integral_3 + b3 * u(k01) ^ 2 * mode_integral_4 + b4 * abs(u(k01)) ^ 3 * mode_integral_5 + b5 * u(k01) ^ 4 * mode_integral_6) * udot(k01);
        else

            for k02 = 1:matrixsize
                g(k01) = g(k01) + CC(k01, k02) * udot(k02) + KK(k01, k02) * u(k02);
            end

        end

    end

    g = g';
    Kc = CC;
    Kc(mode_number, mode_number) = Kc(mode_number, mode_number) - (b1 * mode_integral_2 + b2 * abs(u(mode_number)) * mode_integral_3 + b3 * u(mode_number) ^ 2 * mode_integral_4 + b4 * abs(u(mode_number)) ^ 3 * mode_integral_5 + b5 * u(mode_number) ^ 4 * mode_integral_6);

    Kk = KK;
    Kk(mode_number, mode_number) = Kk(mode_number, mode_number) - udot(mode_number) * b2 * sign(u(mode_number)) * mode_integral_3 - 2 * b3 * udot(mode_number) * u(mode_number) * mode_integral_4 - 3 * udot(mode_number) * b4 * abs(u(mode_number)) ^ 2 * mode_integral_5 * sign(u(mode_number)) - 4 * b5 * udot(mode_number) * u(mode_number) ^ 3 * mode_integral_6;

    ks = Kk + gamma * h / beta / h ^ 2 * Kc + 1 / beta / h ^ 2 * MM;
end
