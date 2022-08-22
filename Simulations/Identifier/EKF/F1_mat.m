function F1 = F1_mat(X_k, N_layers, V,deltaX,deltaT,eps, Tank, heatState, HE)
    %X_k = [T_k; Theta_k] where theta_k = [alpha_k, Dc_k, UL_k]
    %T_k = [T1_k, T2_k, ..., TN_k, T_amb ]
    [Z1, ~, ~] = Matrix__(N_layers, V,deltaX,deltaT,eps, Tank, heatState, HE);
    
    %Jacobian of the function f1 given theta_k's elements
    a = 0.5/(deltaX^2);
    Jac_f1_theta_k = zeros(N_layers+1, 3);
    for i = 2:N_layers-1
        Jac_f1_theta_k(i,1) = X_k(N_layers+3)*(-X_k(i-1)+ 2*X_k(i) - X_k(i+1));
        Jac_f1_theta_k(i,2) = X_k(N_layers+2)*(-X_k(i-1)+ 2*X_k(i) - X_k(i+1));
        Jac_f1_theta_k(i,3) = X_k(i)/2 - X_k(N_layers+1);
    end
    %F1 is a compond matrix
    F1 = [Z1    Jac_f1_theta_k;
        zeros(3,N_layers+1) eye(3)];
end