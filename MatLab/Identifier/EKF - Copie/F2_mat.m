function F2 = F2_mat(X_k, N_layers, V, deltaX,deltaT, eps, Tank)
    %X_k = [T1_k, T2_k, ..., TN_k, T_amb ]
    [~, Z2] = Matrix(N_layers, V, deltaX, deltaT, eps, Tank);

    F2 = Z2;
end