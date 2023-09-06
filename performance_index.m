function J = performance_index(N, X_trim, X, U)
    import casadi.*

    z_trim = X_trim(1);
    vel_trim = X_trim(2:3);
    z_f = full(X(2, end));
    vel_f = full(X(3:4, end));
    Phi_f = norm(z_f - z_trim) + norm(vel_f - vel_trim);
    
    J = Phi_f;
    
    % Running cost
    for i = 1:N
        Vx = X(3, i); 
        Vz = X(4, i);
        F_r = U(1, i);
        F_p = U(2, i);
    
        energy = F_r * (-Vz) + F_p * Vx;
        J = J + energy;
    end


end