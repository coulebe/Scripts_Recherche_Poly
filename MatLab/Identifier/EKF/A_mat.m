function A = A_mat(theta, deltaX, deltaT, N)
    %theta = [alpha, DC, UL]
    A = zeros(N, N);
    temp = [-3, 4, -1, zeros(1, N-3)];
    A(1,:) = temp;
    temp = [zeros(1, N-3),1, -4, 3 ];
    A(N,:) = temp;
    a = theta(1)*theta(2)/(deltaX^2);
    for i = 0:N-3
       temp = [zeros(1,i), a, -2*a-theta(3)+1/deltaT, a, zeros(1, N-i-3)];
       A(i+2,:) =temp;
    end
    A = deltaT * A;
end