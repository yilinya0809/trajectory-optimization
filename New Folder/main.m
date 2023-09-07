% main
%% Get trim
clc; clear; close all;

import casadi.*
LC62;

h_trim = 10;
VT_trim = 45;

[X_trim, U_trim] = get_trim(h_trim, VT_trim);

test_X = X_trim;
test_U = U_trim;

test_dX = f(test_X, test_U);


%% Optimization

dt = 0.05;
tf = 20;
N = tf / dt;

opti = casadi.Opti();

X = opti.variable(3, N+1);
z = X(1, :);
Vx = X(2, :);
Vz = X(3, :);

U = opti.variable(3, N);
Fr = U(1, :);
Fp = U(2, :);
theta = U(3, :);

opti.minimize(norm(Fr) + norm(Fp) + 500 * norm(theta));

for k = 1:N
   k1 = f(X(:, k),               U(:, k));
   k2 = f(X(:, k) + dt / 2 * k1, U(:, k));
   k3 = f(X(:, k) + dt / 2 * k2, U(:, k));
   k4 = f(X(:, k) + dt * k3,     U(:, k));
   x_next = X(:, k) + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
   opti.subject_to(X(:, k+1) == x_next); % close the gaps
end    

opti.subject_to(0 <= Fr <= 500);
opti.subject_to(0 <= Fp <= 150);
opti.subject_to(-deg2rad(50) <= theta <= deg2rad(50));
opti.subject_to(X_trim(1) - 1 <= z <= X_trim(1) + 1);

opti.subject_to(z(1) == X_trim(1));
opti.subject_to(Vx(1) == 20);
opti.subject_to(Vz(1) == 0);
opti.subject_to(Fr(1) == m * g);
opti.subject_to(Fp(1) == 0);
opti.subject_to(theta(1) == deg2rad(-10));

opti.subject_to(z(end) == X_trim(1));
opti.subject_to(Vx(end) == X_trim(2));
opti.subject_to(Vz(end) == X_trim(3));
opti.subject_to(Fr(end) == U_trim(1));
opti.subject_to(Fp(end) == U_trim(2));
opti.subject_to(theta(end) == U_trim(3));

opti.set_initial(z, X_trim(1));
opti.set_initial(Vx, X_trim(2) / 2);
opti.set_initial(Vz, X_trim(3) / 2);
opti.set_initial(Fr, m * g / 2);
opti.set_initial(Fp, U_trim(2) / 2);
opti.set_initial(theta, U_trim(3) / 2);

opti.solver('ipopt');
sol = opti.solve();


%% Plot

figure(1)
plot(sol.value(z));


figure(2)
subplot(2,1,1)
plot(sol.value(Vx));

subplot(2,1,2)
plot(sol.value(Vz));


figure(3)
subplot(3,1,1)
plot(sol.value(Fr))

subplot(3,1,2)
plot(sol.value(Fp))

subplot(3,1,3)
plot(sol.value(theta))