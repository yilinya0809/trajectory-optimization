function [Fx_w, Fz_w] = get_Fw(z, vel)
    LC62();
    Vx = vel(1); Vz = vel(2);
    V = norm(vel);

%     rho = get_rho(-z);
    rho = 1.2241;
    qbar = 0.5 * rho * V^2;
    alp = atan2(Vz, Vx);
    [CL, CD, ~] = aero_coeff(alp);
    CL = 0; CD = 0;
    
    Fx_w = -qbar * S * CD;
    Fz_w = -qbar * S * CL;

end