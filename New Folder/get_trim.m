function [X_trim, U_trim] = get_trim(h, V)
LC62();

% Fr_max = polyval(th_r, 1) * g / 1000; % 159.2089
Fr_max = 0; % 159.2089, 0.0064
Fp_max = interp1(cmd, th_p, 1, 'spline');
% Fp_max = 0; % 82.8813, 0.0318

q0 = [0.0 0 0];
lb = [0.0 0 0];
ub = [deg2rad(20) Fr_max Fp_max];

q = fmincon(@trim_cost, q0, [], [], [], [], lb, ub);

alp_trim = q(1);
Fr_trim = q(2);
Fp_trim = q(3);

Vx_trim = V * cos(alp_trim);
Vz_trim = V * sin(alp_trim);
vel_trim = [Vx_trim; Vz_trim];

X_trim = [-h; vel_trim];
U_trim = [Fr_trim; Fp_trim; alp_trim];

end