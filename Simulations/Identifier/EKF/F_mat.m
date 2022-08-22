function F = F_mat(theta, deltaX, deltaT,  N, Temp)
    %theta = [alpha; DC; UL]
    %Temp(N);
    A = A_mat(theta, deltaX, deltaX, N);
    %
    M = zeros(N, 3);
    for i =2:N-1
        a = Temp(i-1) - 2*Temp(i) + Temp(i+1)/(deltaX^2); 
        M(i,:) = [theta(1)*a, theta(2)*a, -Temp(i)*theta(3)];
    end
    M = M*deltaT;
    %
    F = [A, M; zeros(3,N), eye(3)];
end