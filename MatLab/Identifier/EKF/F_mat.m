function F = F_mat(X_k, N_layers, V,deltaX,deltaT,Q_vec)
    %X_k = [T_k; Theta_k] where theta_k = [alpha_k, Dc_k, UL_k]
    %T_k = [T1_k, T2_k, ..., TN_k, T_amb ]
    
    %calculate f2(X_k,Q)
    [~, Z2] = Matrix(N_layers, V, deltaX, deltaT,X_k(N_layers+2), X_k(N_layers+3), X_k(N_layers+4));
    B = [Z2, zeros(N_layers+1, 3);
         zeros(3,N_layers+1), eye(3)];
    C = [Q_vec; zeros(4,1)];%Dans notre cas Q_vec ne contient que les N couches alors que notre état augmenté contient les N couches+ la température ambiante + les 3 paramètres
    f2 = B*X_k + C;
    
    %Now we can calculate F1 applied on f2 
    F1_of_f2 = F1_mat(f2, N_layers, V,deltaX,deltaT);
    
    %We calculate F2 applied on X_k and Q
    F2_of_X_k = F2_mat(X_k, N_layers, V,deltaX,deltaT);
    
    F  = F1_of_f2 \F2_of_X_k;
    
end