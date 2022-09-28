function X_k_plus_1 = f(X_k, N_layers,  deltaX, deltaT, eps, Tank, Q_vec, T_in)
    %Our System's dynamic based on its PDE
    %returns X_k+1 the augmented matrice for the next time step, size(N+4xN+4 matrice)
    %X_k = [T_k; V_k]
    %T_k = [T1_k; T2_k; ..., TN_k; T_amb ]
    [Z1, Z2] = Matrix(N_layers, X_k(N_layers+2), deltaX, deltaT, eps, Tank);
    A = [Z1, zeros(N_layers+1, 1);
         zeros(1,N_layers+1), 1];
    B = [Z2, zeros(N_layers+1, 1);
         zeros(1,N_layers+1), 1];
    C = [Q_vec; zeros(2,1)]; %Dans notre cas Q_vec ne contient que les N couches alors que notre état augmenté contient les N couches+ la température ambiante + V_k
    
    X_k_plus_1 = A\(B*X_k + C);
    if(X_k(N_layers+2) > 0)
        X_k_plus_1(1) = T_in; 
    else
       X_k_plus_1(1) = (4*X_k_plus_1(2) - X_k_plus_1(3))/3;
    end
end