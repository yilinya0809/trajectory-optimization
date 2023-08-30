function trim_cost = trim_cost(q)
    LC62();

    alp = q(1); 
    F_r = q(2);
    F_p = q(3);
    theta = q(4);

    global h V;
    vel = V * [cos(alp); sin(alp)];
    [Fx_w, Fz_w] = get_Fw(h, vel);
    Vx_dot = (F_p + Fx_w)/m - g * sin(theta);
    Vz_dot = (-F_r + Fz_w)/m + g* cos(theta);
    vel_dot = [Vx_dot; Vz_dot];

    trim_cost = vel_dot.' * vel_dot;

end