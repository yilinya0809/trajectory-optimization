function [X_trim, U_trim] = get_trim(h, V)
global h V;
LC62();

Fr_max = polyval(th_r, 1) * g/1000; % 159.2089
Fp_interp = griddedInterpolant(cmd, th_p);
Fp_max = Fp_interp(1); % 91.5991

q0 = [0.0 0 0 0];
lb = [0.0 0 0 0];
ub = [deg2rad(20) Fr_max Fp_max pi/2];

q = fmincon(@trim_cost, q0, [], [], [], [], lb, ub);

alp_trim = q(1);
Fr_trim = q(2);
Fp_trim = q(3);
theta_trim = q(4);

Vx_trim = V * cos(alp_trim);
Vz_trim = V * sin(alp_trim);
vel_trim = [Vx_trim; Vz_trim];

X_trim = [h; vel_trim];
U_trim = [Fr_trim; Fp_trim; rad2deg(theta_trim)];

end