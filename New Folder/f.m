function dX = f(X, U)
LC62

z = X(1);
vel = X(2:3);

Fr = U(1);
Fp = U(2);
theta = U(3);

Vx = vel(1); Vz = vel(2);
VT = norm(vel);

qbar = 0.5 * rho * VT^2;
alpha = atan2(Vz, Vx);
clgrid = casadi.interpolant('CLGRID', 'bspline', {alp}, C_L);
cdgrid = casadi.interpolant('CDGRID', 'bspline', {alp}, C_D);
CL = full(clgrid(alpha));
CD = full(cdgrid(alpha));

Fx_w = -qbar * S * CD;
Fz_w = -qbar * S * CL;

R = [cos(theta), sin(theta);
    -sin(theta), cos(theta)];

Vx_dot = (Fp + Fx_w) / m - g * sin(theta);
Vz_dot = (-Fr + Fz_w) / m + g * cos(theta);

dpos = R * vel;
dvel = [Vx_dot; Vz_dot];

dX = [dpos(2); dvel];
end