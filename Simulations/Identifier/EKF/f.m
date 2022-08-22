function X_k_plus_1 = f(X_k, N_layers, V,deltaX,deltaT,eps, Tank, heatState, HE, T_in)
    %Our System's dynamic based on its PDE
    %returns X_k+1 the augmented matrice for the next time step, size(N+4xN+4 matrice)
    %X_k = [T_k; Theta_k] where theta_k = [alpha_k, Dc_k, UL_k]
    %T_k = [T1_k, T2_k, ..., TN_k, T_amb ]
    [Z1, Z2, Z3] = Matrix__(N_layers, V,deltaX,deltaT,eps, Tank, heatState, HE);
    A = [Z1, zeros(N_layers+1, 3);
         zeros(3,N_layers+1), eye(3)];
    B = [Z2, zeros(N_layers+1, 3);
         zeros(3,N_layers+1), eye(3)];
    C = [Z3; zeros(3,1)];
    
    X_k_plus_1 = A\(B*X_k + C);
    if(V > 0)
        X_k_plus_1(1) = T_in; 
    else
       X_k_plus_1(1) = (4*X_k_plus_1(2) - X_k_plus_1(3))/3;
    end
end