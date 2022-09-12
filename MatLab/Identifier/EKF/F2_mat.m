function F2 = F2_mat(X_k, N_layers, V,deltaX,deltaT)
    %X_k = [T_k; Theta_k] where theta_k = [alpha_k, Dc_k, UL_k]
    %T_k = [T1_k, T2_k, ..., TN_k, T_amb ]
    [~, Z2] = Matrix(N_layers, V, deltaX, deltaT,X_k(N_layers+2), X_k(N_layers+3), X_k(N_layers+4));
    
    %Jacobian of the function f1 given theta_k's elements
    a = 0.5/(deltaX^2);
    Jac_f2_theta_k = zeros(N_layers+1, 3);
    for i = 2:N_layers-1
        Jac_f2_theta_k(i,1) = X_k(N_layers+3)*a*(X_k(i-1)- 2*X_k(i) + X_k(i+1));
        Jac_f2_theta_k(i,2) = X_k(N_layers+2)*a*(X_k(i-1)+ 2*X_k(i) + X_k(i+1));
        Jac_f2_theta_k(i,3) = -X_k(i)/2;
    end
    %F1 is a compond matrix
    F2 = [Z2    Jac_f2_theta_k;
        zeros(3,N_layers+1) eye(3)];
end