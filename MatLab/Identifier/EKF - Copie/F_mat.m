function F = F_mat(X_k, N_layers, deltaX,deltaT,Q_vec, eps, Tank)
    %X_k = [T_k; V_k]
    %T_k = [T1_k, T2_k, ..., TN_k, T_amb ]
    
    %calculate f2(X_k,Q)
    [~, Z2] = Matrix(N_layers, X_k(N_layers+2), deltaX, deltaT, eps, Tank);
    B = [Z2, zeros(N_layers+1, 1);
         zeros(1,N_layers+1), 1];
    C = [Q_vec; zeros(2,1)];%Dans notre cas Q_vec ne contient que les N couches alors que notre état augmenté contient les N couches+ la température ambiante + V_k
    f2 = B*X_k + C;
    
    %Now we can calculate F1 applied on f2 
    F1_of_f2 = F1_mat(f2, N_layers, deltaX,deltaT, eps, Tank);
    
    %We calculate F2 applied on X_k and Q
    F2_of_X_k = F2_mat(X_k, N_layers, deltaX,deltaT, eps, Tank);
    
    F  = F1_of_f2 \F2_of_X_k;
    
end