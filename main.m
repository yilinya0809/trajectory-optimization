% main
LC62();

h = 10;
V = 45;

[X_trim, U_trim] = get_trim(h, V);
z_trim = X_trim(1);
vel_trim = X_trim(2:3);

Fr_trim = U_trim(1);
Fp_trim = U_trim(2);
theta_trim = U_trim(3);

global dt
dt = 0.01;
tf = 10;

U0 = [0; 0; 0];