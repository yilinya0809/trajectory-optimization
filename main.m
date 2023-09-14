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

% dt = 0.02; tf = 20; 
% N = tf / dt;
N = 100;

opti = casadi.Opti();
   
X = opti.variable(3, N+1);
z = X(1,:);
Vx = X(2,:);
Vz = X(3,:);

U = opti.variable(3, N);
Fr = U(1,:);
Fp = U(2,:);
theta = U(3,:);
T = opti.variable(); % final time

% objective
W = diag([1 1 4000]);  % input weight
opti.minimize(dot(U, W * U) + 100 * T^2);

% dynamic constraints
dt = T / N;
for k = 1:N
   k1 = f(X(:, k),               U(:, k));
   k2 = f(X(:, k) + dt / 2 * k1, U(:, k));
   k3 = f(X(:, k) + dt / 2 * k2, U(:, k));
   k4 = f(X(:, k) + dt * k3,     U(:, k));
   x_next = X(:, k) + dt / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
   opti.subject_to(X(:, k+1) == x_next); % close the gaps
end    

% input constraints
Fr_max = 6 * th_r_max; % 955.2534
Fp_max = 2 * th_p_max; % 183.1982
opti.subject_to(0 <= Fr <= Fr_max);
opti.subject_to(0 <= Fp <= Fp_max);
opti.subject_to(-deg2rad(50) <= theta <= deg2rad(50));
z_eps = 2;
opti.subject_to(X_trim(1) - z_eps <= z <= X_trim(1) + z_eps);
opti.subject_to(0 <= T <= 20);


% initial boundary conditions
opti.subject_to(z(1) == X_trim(1));
opti.subject_to(Vx(1) == 0);
opti.subject_to(Vz(1) == 0);
opti.subject_to(Fr(1) == m * g);
opti.subject_to(Fp(1) == 0);
opti.subject_to(theta(1) == deg2rad(0));

% terminal boundary conditions
opti.subject_to(z(end) == X_trim(1));
opti.subject_to(Vx(end) == X_trim(2));
opti.subject_to(Vz(end) == X_trim(3));
opti.subject_to(Fr(end) == U_trim(1));
opti.subject_to(Fp(end) == U_trim(2));
opti.subject_to(theta(end) == U_trim(3));



% initial values for solver
load('init_traj.mat','X_init_traj','U_init_traj');
opti.set_initial(z, X_init_traj(1,:));
opti.set_initial(Vx, X_init_traj(2,:));
opti.set_initial(Vz, X_init_traj(3,:));

opti.set_initial(Fr, U_init_traj(1,:));
opti.set_initial(Fp, U_init_traj(2,:));
opti.set_initial(theta, U_init_traj(3,:));

% opti.set_initial(z, X_trim(1));
% opti.set_initial(Vx, X_trim(2) / 2);
% opti.set_initial(Vz, X_trim(3) / 2);
% opti.set_initial(Fr, m * g / 2);
% opti.set_initial(Fp, U_trim(2) / 2);
% opti.set_initial(theta, U_trim(3) / 2);
opti.set_initial(T, 20);

% Vx_init_traj = zeros(1,N+1);
% for i = 1:N+1
%     Vx_init_traj(i) = X_trim(2) / (N+1) * i;
% end
% 
% Fr_init_traj = zeros(1, N);
% theta_init_traj = zeros(1, N);
% 
% for i = 1:N
%     Fr_init_traj(i) = (U_trim(1) - m*g) / N * i + m*g;
%     theta_init_traj(i) = U_trim(3) / N * i;
% end
% 
% opti.set_initial(z, X_trim(1));
% opti.set_initial(Vx, Vx_init_traj);
% opti.set_initial(Vz, X_trim(3) / 2);
% opti.set_initial(Fr, Fr_init_traj);
% opti.set_initial(Fp, U_trim(2) / 2);
% opti.set_initial(theta, theta_init_traj);



opti.solver('ipopt');
sol = opti.solve();


%% Plot
tf = sol.value(T);
tspan = linspace(0, tf, N+1);

figure(1)
plot(tspan, -sol.value(z), 'k'), grid on
xlabel('Time, s')
ylabel('Altitude, m')


figure(2)
subplot(2,1,1)
plot(tspan, sol.value(Vx), 'k'), grid on;
xlabel('Time, s')
ylabel('V_x, m/s')

subplot(2,1,2)
plot(tspan, sol.value(Vz), 'k'), grid on;
xlabel('Time, s')
ylabel('V_z, m/s')


figure(3)
subplot(3,1,1)
plot(tspan(1:end-1), sol.value(Fr), 'k'), grid on
xlabel('Time, s')
ylabel('Rotor thrust, N')

subplot(3,1,2)
plot(tspan(1:end-1), sol.value(Fp), 'k'), grid on
xlabel('Time, s')
ylabel('Pusher thrust, N')

subplot(3,1,3)
plot(tspan(1:end-1), rad2deg(sol.value(theta)), 'k'), grid on
xlabel('Time, s')
ylabel('\theta, deg')