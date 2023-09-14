clc; clear; close all;

X_ndi = readmatrix('X_traj.csv');
U_ndi = readmatrix('U_traj.csv');

ts = X_ndi(1, :);
z_ndi = X_ndi(2, :);
Vx_ndi = X_ndi(3, :);
Vz_ndi = X_ndi(4, :);

Fr_ndi = U_ndi(2, :);
Fp_ndi = U_ndi(3, :);
theta_ndi = U_ndi(4, :);

N = 100; Tf = 20 / 0.01;
dt = Tf / N;

z_init_traj = zeros(1, N+1);
Vx_init_traj = zeros(1, N+1);
Vz_init_traj = zeros(1, N+1);

for i = 0:N
    z_init_traj(i+1) = z_ndi(dt * i + 1);
    Vx_init_traj(i+1) = Vx_ndi(dt * i + 1);
    Vx_init_traj(i+1) = Vz_ndi(dt * i + 1);
end

Fr_init_traj = zeros(1, N);
Fp_init_traj = zeros(1, N);
theta_init_traj = zeros(1, N);

for i = 0:N-1
    Fr_init_traj(i+1) = Fr_ndi(dt * i + 1);
    Fp_init_traj(i+1) = Fp_ndi(dt * i + 1);
    theta_init_traj(i+1) = theta_ndi(dt * i + 1);
end

X_init_traj = vertcat(z_init_traj, Vx_init_traj, Vz_init_traj);
U_init_traj = vertcat(Fr_init_traj, Fp_init_traj, theta_init_traj);

save("init_traj.mat",'X_init_traj','U_init_traj');