%% Get trim & initial trajectory by LQR
import casadi.*
LC62();

h_trim = 10; 
VT_trim = 45;
[X_trim_, U_trim] = get_trim(h, V);

% X_trim = [0; X_trim_];

test_X_dot = f(X_trim, U_trim);

global dt
dt = 0.02; tf = 10; 
N = tf / dt; t = 0:dt:tf;

ptrb = 1e-9;
[A, B] = linearization(X_trim, U_trim, ptrb);
Q = diag([0.05,4,4,1]);
R = diag([1,1,1]);

[K, ~, ~] = lqr(A, B, Q, R);

X0 = [0; 10; 0; 0];
X_lqr = zeros(4, length(t)); X_lqr(:, 1) = X0;
U_lqr = zeros(3, length(t)); 

for i = 2:length(t)
    U_lqr(:, i-1) = U_trim - K * (X_lqr(:, i-1) - X_trim);
    X_dot = A*X_lqr(:, i-1) + B*U_lqr(:,i-1);
    X_lqr(:, i) = X_lqr(:, i-1) + dt * X_dot;
end

% plot LQR
figure(1)
subplot(2,2,1)
plot(t, X_lqr(1, :)); hold on
plot(t, X_trim(1) * ones(length(t))); hold off
ylabel('x');

subplot(2,2,2)
plot(t, X_lqr(2, :)); hold on
plot(t, X_trim(2) * ones(length(t))); hold off 
ylabel('z')

subplot(2,2,3)
plot(t, X_lqr(3, :)); hold on
plot(t, X_trim(3) * ones(length(t))); hold off 
ylabel('Vx')

subplot(2,2,4)
plot(t, X_lqr(4, :)); hold on
plot(t, X_trim(4) * ones(length(t))); hold off 
ylabel('Vz')


figure(2)
subplot(3,1,1)
plot(t, U_lqr(1, :)); hold on
plot(t, U_trim(1) * ones(length(t))); hold off
ylabel('F_r')

subplot(3,1,2)
plot(t, U_lqr(2, :)); hold on
plot(t, U_trim(2) * ones(length(t))); hold off
ylabel('F_p')

subplot(3,1,3)
plot(t, U_lqr(3, :)); hold on
plot(t, U_trim(3) * ones(length(t))); hold off
ylabel('theta')