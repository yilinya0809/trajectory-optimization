function [Fx_w, Fz_w] = get_Fw(z, vel)
    LC62();
    
    Vx = vel(1); Vz = vel(2);
    V = norm(vel);

    qbar = 0.5 * rho * V^2;
    alpha = atan2(Vz, Vx);

    alpha = 0;
    % CL = interp1(alp, C_L, alpha, 'spline');
    % CD = interp1(alp, C_D, alpha, 'spline');
    CL = 0.1931;
    CD = 0.0617;
    % CL = 0;
    % CD = 0;
    
    Fx_w = -qbar * S * CD;
    Fz_w = -qbar * S * CL;

end