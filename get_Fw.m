function [Fx_w, Fz_w] = get_Fw(z, vel)
    LC62();
    Vx = vel(1); Vz = vel(2);
    V = norm(vel);

    rho = get_rho(-z);
    qbar = 0.5 * rho * V^2;
    alp = atan2(Vz, Vx);
    [CL, CD, ~] = aero_coeff(alp);
    
    Fx_w = -qbar * S * CD;
    Fz_w = -qbar * S * CL;

end