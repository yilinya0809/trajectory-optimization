% main
import casadi.*
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
dt = 0.02;
tf = 10;
N = tf / dt;

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

opti.subject_to(X(:,1) == [0;0;0;0]);
opti.subject_to(U(:,1) == [0;0;0]);


for k=1:N
   k1 = set_dot(X(:,k),         U(:,k));
   k2 = set_dot(X(:,k)+dt/2*k1, U(:,k));
   k3 = set_dot(X(:,k)+dt/2*k2, U(:,k));
   k4 = set_dot(X(:,k)+dt*k3,   U(:,k));
   x_next = X(:,k) + dt/6*(k1+2*k2+2*k3+k4);
   opti.subject_to(X(:,k+1)==x_next); % close the gaps
end

opti.subject_to(X(:,end) == X_trim);
opti.subject_to(U(:,end) == U_trim);    

opti.minimize(performance_index(N, X_trim, X, U));

opti.set_initial(X, [0;0;0;0]);
opti.set_initial(U, [0;0;0]);

opti.solver('ipopt');
sol = opti.solve();


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







