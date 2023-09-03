function J = performance_index(N, X_trim, X, U)
    
    z_trim = X_trim(1);
    vel_trim = X_trim(2:3);
    z_f = X(end, 2);
    vel_f = X(end, 3:4);
    Phi_f = norm(z_f - z_trim) + norm(vel_f - vel_trim);
    
    J = Phi_f;
    
    % Running cost
    for i = 1:N
        Vx = X(i, 3); 
        Vz = X(i, 4);
        F_r = U(i, 1);
        F_p = U(i, 2);
    
        energy = F_r * Vz + F_p * Vx;
        J = J + energy;
    end


end