function [A, B] = linearization(X, U, ptrb)

    n = length(X);
    m = length(U);
    A = zeros(n, n);
    B = zeros(n, m);

    for i = 1:n
        ptrbvec_X = zeros(n,1);
        ptrbvec_X(i) = ptrb;
        X_ptrb = X + ptrbvec_X;
        dfdX = (set_dot(X_ptrb, U) - set_dot(X, U)) / ptrb;
        for j = 1:n
            A(j, i) = dfdX(j);
        end
    end

    for i = 1:m
        ptrbvec_U = zeros(m,1);
        ptrbvec_U(i) = ptrb;
        U_ptrb = U + ptrbvec_U;
        dfdU = (set_dot(X, U_ptrb) - set_dot(X, U)) / ptrb;
        for j = 1:n
            B(j, i) = dfdU(j);
        end
    end

end