function B = B_mat(theta,  N)
    %theta = [alpha; DC; UL]
    B = zeros(N, N+1);
    B(1:N, 1:N)= eye(N);
    B(:,end) = theta(3)*[0; ones(N-2,1); 0];
    
    
end