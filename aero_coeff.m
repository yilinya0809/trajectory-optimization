function [CL, CD, Cm] = aero_coeff(alpha)
    import casadi.*
    LC62();
    
    CL_interp = casadi.interpolant('CL_interp','bspline',{alp},C_L);
    CD_interp = casadi.interpolant('CD_interp','bspline',{alp},C_D);
    Cm_interp = casadi.interpolant('Cm_interp','bspline',{alp},C_m);
    CL = full(CL_interp(alpha));
    CD = full(CD_interp(alpha));
    Cm = full(Cm_interp(alpha));
end