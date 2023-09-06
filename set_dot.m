function X_dot = set_dot(Xt, Ut)

LC62();


pos = Xt(1:2);
vel = Xt(3:4);

[Fx_w, Fz_w] = get_Fw(pos(2), vel);
F_r = Ut(1);
F_p = Ut(2);
theta = Ut(3);

R = [cos(theta), sin(theta);
    -sin(theta), cos(theta)];

Vx_dot = (F_p + Fx_w)/m - g*sin(theta);
Vz_dot = (-F_r + Fz_w)/m + g*cos(theta);

pos_dot = R * vel;
vel_dot = [Vx_dot; Vz_dot];

X_dot = [pos_dot; vel_dot];
end