clear ;
load('Donnee_sim.mat');
data = data - T_amb;
%%
%Variance matrices
Q = eye(N+3,N+3); %Process noise variance
R = eye(N,N); %measurements noise variance
%Observation matrix/H
H_ = [eye(N), zeros(N,3)];
% H_ = zeros(N, N+3);
% H_(1, 1) = 1;
% H_(N, N) = 1;
%initialisation
Temp_hat = data(:,1);
theta_hat = 1*[1e2;1e3;1e-3];
P_k = [eye(N), zeros(N,3); zeros(3,N), 1e3*eye(3)];
%Loop
[~, T] = size(data);
Theta = zeros(3, T);
Theta(:,1) = theta_hat;
for i = 2:T
    %Observtaion
    Y_k = data(:,i);
    %Coefficients
    U = U_mat(Tank, HeatElem, N, H(:,i), deltaX, deltaT);
    
    %Prediction step
    X_ = f(Temp_hat, theta_hat, U, deltaX, deltaT, N);
    
    %Error covariance propagation
    F = F_mat(theta_hat, deltaX,deltaT, N, Temp_hat);
    P_k_ = F* P_k* (F.') + Q;
    %Kalman Gain matrix
    K_k = P_k_ * H_' * (H_ * P_k_ * H_' + R)^-1 ;
    
    %Update step
    Y_ = Y_k - H_*X_;
    X_k_hat = X_ + K_k*Y_;
    P_k = (eye(N+3) - K_k*H_)*P_k_;
    
    Temp_hat = X_(1:N);
    theta_hat = X_(N+1:end);
    Theta(:,i) = theta_hat;
end
