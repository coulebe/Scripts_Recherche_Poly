function F1 = F1_mat(X_k, N_layers, V, deltaX,deltaT, eps, Tank)
    %V_k = [T1_k, T2_k, ..., TN_k, T_amb ]
    [Z1, ~] = Matrix(N_layers, V, deltaX, deltaT, eps, Tank);
    
    %Jacobian of the function f1 given theta_k's elements


    %F1 is a compond matrix
    F1 = Z1;
end