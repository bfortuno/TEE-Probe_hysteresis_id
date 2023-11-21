function result = hysteresis_model(theta_i, theta_i_minus_1, phi)
    if (theta_i >= theta_i_minus_1) && (theta_i > -56 && theta_i <= -40)
        a = -4.336e-5;
        b = -6.107e-3;
        c = -2.866e-1;
        d = -4.282;
        e = 0;
        f = -4.8e-4;
        g = -5.2e-4;
    elseif (theta_i >= theta_i_minus_1) && (theta_i > -40 && theta_i <= 0)
        a = -1.258e-6;
        b = 2.377e-5;
        c = -6.791e-4;
        d = 4.914e-2;
        e = 0;
        f = -4.8e-4;
        g = -5.2e-4;
    elseif (theta_i >= theta_i_minus_1) && (theta_i > 0 && theta_i <= 34)
        a = -1.258e-6;
        b = 2.377e-5;
        c = -6.791e-4;
        d = 4.914e-2;
        e = 0;
        f = 0;
        g = 0;
    elseif (theta_i >= theta_i_minus_1) && (theta_i > 34 && theta_i <= 100)
        a = 2.662e-8;
        b = -5.81e-6;
        c = -4.703e-4;
        d = 5.316e-2;
        e = -100;
        f = -1.6e-4;
        g = 5.6e-4;
    elseif (theta_i <= theta_i_minus_1) && (theta_i <= 13 && theta_i > -50)
        a = -8.707e-7;
        b = -4.805e-5;
        c = -4.784e-4;
        d = -1.456e-2;
        e = 0;
        f = 0;
        g = 0;
    elseif (theta_i <= theta_i_minus_1) && (theta_i > 13 && theta_i <= 100)
        a = -1.517e-7;
        b = 2.034e-5;
        c = -2.903e-3;
        d = 1.35e-2;
        e = -13;
        f = -6e-4;
        g = -3e-4;
    else
        result = 0;
        return;
    end
    % theta_i = theta_i * pi / 180;
    result = a * theta_i^3 + b * theta_i^2 + c * theta_i + d + (f * phi^2 + g * phi) * (theta_i + e);
end
