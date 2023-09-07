function [X_trim, U_trim] = get_trim(h_trim, VT_trim)

Fr_max = 0;
Fp_max = 1000;

q0 = [0 50 0];
lb = [0 0 0];
ub = [deg2rad(20) Fr_max Fp_max];

q = fmincon(@(q)trim_cost(q, h_trim, VT_trim), q0, [], [], [], [], lb, ub);

alp_trim = q(1);
vel_trim = VT_trim * [cos(alp_trim); sin(alp_trim)];

Fr_trim = q(2);
Fp_trim = q(3);
theta_trim = alp_trim;

X_trim = [-h_trim; vel_trim];
U_trim = [Fr_trim; Fp_trim; theta_trim];

end