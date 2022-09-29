function X_k_plus_1 = f(X_k, N_layers, V, deltaX, deltaT, eps, Tank, Q_vec, T_in)
    %Our System's dynamic based on its PDE
    %returns X_k+1 the augmented matrice for the next time step, size(N+4xN+4 matrice)

    %X_k = [T1_k; T2_k; ...; TN_k; T_amb ]
    [Z1, Z2] = Matrix(N_layers, V, deltaX, deltaT, eps, Tank);

    C = [Q_vec; 0]; %Dans notre cas Q_vec ne contient que les N couches alors que notre état augmenté contient les N couches+ la température ambiante
    
    X_k_plus_1 = Z1\(Z2*X_k + C);
    if(V > 0)
        X_k_plus_1(1) = T_in; 
    else
       X_k_plus_1(1) = (4*X_k_plus_1(2) - X_k_plus_1(3))/3;
    end
end