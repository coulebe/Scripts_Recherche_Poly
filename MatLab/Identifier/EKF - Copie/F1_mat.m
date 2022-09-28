function F1 = F1_mat(X_k, N_layers,deltaX,deltaT, eps, Tank)
    %X_k = [T_k; V_k] 
    %T_k = [T1_k, T2_k, ..., TN_k, T_amb ]
    [Z1, ~] = Matrix(N_layers, X_k(N_layers+2), deltaX, deltaT, eps, Tank);
    
    %Jacobian of the function f1 given theta_k's elements

    Jac_f1_V_k = zeros(N_layers+1, 1);
    for i = 2:N_layers-1
        Jac_f1_V_k(i,1) = (X_k(i+1) - X_k(i-1))/(4*deltaX);
    end
    %F1 is a compond matrix
    F1 = [Z1    Jac_f1_V_k;
        zeros(1,N_layers+1) 1];
end