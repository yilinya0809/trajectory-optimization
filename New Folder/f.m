function dX = f(Xt, Ut)

LC62();

z = Xt(1);
vel = Xt(2:3);

[Fx_w, Fz_w] = get_Fw(z, vel);
F_r = Ut(1);
F_p = Ut(2);
theta = Ut(3);

R = [cos(theta), -sin(theta);
     sin(theta),  cos(theta)];

Vx_dot = (F_p + Fx_w) / m - g * sin(theta);
Vz_dot = (-F_r + Fz_w) / m + g * cos(theta);

pos_dot = R * vel;
vel_dot = [Vx_dot; Vz_dot];

dX = [pos_dot(2); vel_dot];
end