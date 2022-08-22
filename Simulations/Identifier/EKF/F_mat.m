function F = F_mat(X_k, N_layers, V,deltaX,deltaT,eps, Tank, heatState, HE)
    %X_k = [T_k; Theta_k] where theta_k = [alpha_k, Dc_k, UL_k]
    %T_k = [T1_k, T2_k, ..., TN_k, T_amb ]
    
    %calculate f2(X_k,Q)
    [~, Z2, Z3] = Matrix__(N_layers, V,deltaX,deltaT,eps, Tank, heatState, HE);
    B = [Z2, zeros(N_layers+1, 3);
         zeros(3,N_layers+1), eye(3)];
    C = [Z3; zeros(3,1)];
    f2 = B*X_k + C;
    
    %Now we can calculate F1 applied on f2 
    F1_of_f2 = F1_mat(f2, N_layers, V,deltaX,deltaT,eps, Tank, heatState, HE);
    
    %We calculate F2 applied on X_k and Q
    F2_of_X_k = F2_mat(X_k, N_layers, V,deltaX,deltaT,eps, Tank, heatState, HE);
    
    F  = F1_of_f2 \F2_of_X_k;
    
end