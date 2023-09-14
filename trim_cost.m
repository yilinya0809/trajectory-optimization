function trim_cost = trim_cost(q, h_trim, VT_trim)

alp_trim = q(1);
vel_trim = VT_trim * [cos(alp_trim); sin(alp_trim)];

Fr_trim = q(2);
Fp_trim = q(3);
theta_trim = alp_trim;

X_trim = [-h_trim; vel_trim];
U_trim = [Fr_trim; Fp_trim; theta_trim];

dX = f(X_trim, U_trim);
dvel = dX(2:end);

W = diag([1 1]);

trim_cost = dvel.' * W * dvel;
end