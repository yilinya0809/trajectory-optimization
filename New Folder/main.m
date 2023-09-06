% main
%% Get trim
clc; clear; close all;

import casadi.*
LC62();

h = 10;
V = 45;

[X_trim, U_trim] = get_trim(h, V);

test_X = X_trim;
test_U = U_trim;

test_dX = f(test_X, test_U);


%% Optimization

dt = 0.02;
tf = 10;
N = tf / dt;

opti = casadi.Opti();

X = opti.variable(3, N+1);
z = X(1,:);
Vx = X(2,:);
Vz = X(3,:);

U = opti.variable(3, N);
Fr = U(1,:);
Fp = U(2,:);
theta = U(3,:);

opti.minimize(sum(Fp + Fr) + 100*sum(theta));

for k=1:N
   k1 = f(X(:,k),         U(:,k));
   k2 = f(X(:,k)+dt/2*k1, U(:,k));
   k3 = f(X(:,k)+dt/2*k2, U(:,k));
   k4 = f(X(:,k)+dt*k3,   U(:,k));
   x_next = X(:,k) + dt/6*(k1+2*k2+2*k3+k4);
   opti.subject_to(X(:,k+1)==x_next); % close the gaps
end    

% opti.subject_to(Fr >= 0);
% opti.subject_to(Fp >= 0);
opti.subject_to(-deg2rad(80) <= theta <= deg2rad(80));

% opti.subject_to(z(1) == -h);
% opti.subject_to(Fr(1) == 160);
% opti.subject_to(Fp(1) == 0);
% opti.subject_to(theta(1) == 0);
% 
% opti.subject_to(z(end) == X_trim(1));
% opti.subject_to(Vx(end) == X_trim(2));
% opti.subject_to(Vz(end) == X_trim(3));
% opti.subject_to(Fr(end) == U_trim(1));
% opti.subject_to(Fp(end) == U_trim(2));
% opti.subject_to(theta(end) == U_trim(3));
% 
% opti.set_initial(z(N+1), X_trim(1)/2);
% opti.set_initial(Vx(N+1), X_trim(2)/2);
% opti.set_initial(Vz(N+1), X_trim(3)/2);
% opti.set_initial(Fr(N), 200);
% opti.set_initial(Fp(N), 80);
% opti.set_initial(theta(N), deg2rad(20));

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