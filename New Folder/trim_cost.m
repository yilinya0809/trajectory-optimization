function trim_cost = trim_cost(q)
    LC62();
    h = 10;
    V = 45;

    alp = q(1);
    F_r = q(2);
    F_p = q(3);

    vel = V * [cos(alp); sin(alp)];
    [Fx_w, Fz_w] = get_Fw(h, vel);
    Vx_dot = (F_p + Fx_w) / m - g * sin(alp);
    Vz_dot = (-F_r + Fz_w) / m + g* cos(alp);
    vel_dot = [Vx_dot; Vz_dot];

    trim_cost = vel_dot.' * vel_dot;
end