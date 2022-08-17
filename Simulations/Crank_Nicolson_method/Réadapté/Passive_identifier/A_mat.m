function A = A_mat(deltaT, deltaX, Dc_hat_k_plus_1, UL_hat_k_plus_1, N, a)
    A = zeros(N,N);
    A(1,:) = [-3, 4, -1, zeros(1,N-3)];
    A(N,:) = [zeros(1,N-3), 1, -4, 3];
    alpha = deltaT*Dc_hat_k_plus_1/(deltaX^2);
    for i = 2:N-1
        A(i,i-1:i+1) = [alpha, 1+2*alpha+UL_hat_k_plus_1+a, alpha];
    end
    
end