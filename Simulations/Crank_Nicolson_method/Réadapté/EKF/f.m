function X_k_1 = f(T_k, theta_k, U, deltaX, deltaT, N)
    %returns X_k+1 the augmented matrice for the next time step, size(N+3xN+3 matrice)
    %A_ext matrice
    A_ext = [A_mat(theta_k, deltaX, deltaT, N)  zeros(N,3);
             zeros(3,N) eye(3)];
    %U_ext matrice
    U_ext = [U; zeros(3,1)];
    %
    X_k = [T_k; theta_k];
    X_k_1 = A_ext*X_k + U_ext;
end