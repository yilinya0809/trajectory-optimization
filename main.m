% main
clc; clear; close all;

import casadi.*

LC62();
global h V
h = 10; 
V = 45;

[X_trim, U_trim] = get_trim(h, V);
test_X_dot = set_dot([0; X_trim], U_trim);

%% Optimization
global dt
dt = 0.02; tf = 10; 
N = tf / dt; t = 0:dt:tf;

opti = casadi.Opti();

X = opti.variable(4, N+1);
x = X(1,:);
z = X(2,:);
Vx = X(3,:);
Vz = X(4,:);

U = opti.variable(3, N);
Fr = U(1,:);
Fp = U(2,:);
theta = U(3,:);


opti.subject_to(Fr >= 0);
opti.subject_to(Fp >= 0);
opti.subject_to(-deg2rad(20) < theta < deg2rad(20));

opti.subject_to(X(:,1) == [0;-h;0;0]);
opti.subject_to(U(:,1) == [0;0;0]);


for k=1:N
   k1 = set_dot(X(:,k),         U(:,k));
   k2 = set_dot(X(:,k)+dt/2*k1, U(:,k));
   k3 = set_dot(X(:,k)+dt/2*k2, U(:,k));
   k4 = set_dot(X(:,k)+dt*k3,   U(:,k));
   x_next = X(:,k) + dt/6*(k1+2*k2+2*k3+k4);
   opti.subject_to(X(:,k+1)==x_next); % close the gaps
end

opti.subject_to(X(2:4,end) == X_trim);
opti.subject_to(U(:,end) == U_trim);    

J = performance_index(N, X_trim, X, U);
opti.minimize(J);

% Linear traj
z_init_traj = zeros(1, N+1);
Vx_init_traj = zeros(1,N+1);
Vz_init_traj = zeros(1,N+1);
for i = 1:N+1
    z_init_traj(i) = -h;
    Vx_init_traj(i) = X_trim(2) / (N+1) * i;
    Vz_init_traj(i) = X_trim(3) / (N+1) * i;
end

Fr_init_traj = zeros(1, N);
Fp_init_traj = zeros(1, N);
theta_init_traj = zeros(1, N);

for i = 1:N
    Fr_init_traj(i) = U_trim(1) / N * i;
    Fp_init_traj(i) = U_trim(2) / N * i;
    theta_init_traj(i) = U_trim(3) / N * i;
end

% LQR traj
% z_init_traj = X_lqr(:,2).';
% Vx_init_traj = X_lqr(:,3).';
% Vz_init_traj = X_lqr(:,4).';
% 
% Fr_init_traj = U_lqr(:,1).';
% Fp_init_traj = U_lqr(:,2).';
% theta_init_traj = U_lqr(:,3).';
% 
% 
% opti.set_initial(z, z_init_traj);
% opti.set_initial(Vx, Vx_init_traj);
% opti.set_initial(Vz, Vz_init_traj);
% opti.set_initial(Fr, Fr_init_traj);
% opti.set_initial(Fp, Fp_init_traj);
% opti.set_initial(theta, theta_init_traj);

opti.solver('ipopt');
sol = opti.solve();


%% Plot

figure(1)
subplot(2,1,1)
plot(sol.value(x));

subplot(2,1,2)
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







