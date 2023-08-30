function [CL, CD, Cm] = aero_coeff(alpha)
    LC62();
    
    CL_interp = griddedInterpolant(alp, C_L);
    CD_interp = griddedInterpolant(alp, C_D);
    Cm_interp = griddedInterpolant(alp, C_m);

    CL = CL_interp(alpha);
    CD = CD_interp(alpha);
    Cm = Cm_interp(alpha);
end