clear ;
load('Donnee_sim_Null_Flow.mat');
addpath("D:\Alfred\Cours\Projet_Recherche\Poly\Scripts_Recherche_Poly\Simulations\Crank_Nicolson_method\Réadapté");
data = [data; T_amb*ones(1,Max_Simulation_Count)];
% data = data + normrnd(0, 0.1, size(data)); % Add some noise
%%
tic
%Variance matrices
Q = [eye(N+2)   zeros(N+2,2);
    zeros(2,N+2)    1e-8*eye(2)]; %Process noise variance
R = eye(N+1,N+1); %measurements noise variance
%Observation matrix/H
H_ = [eye(N+1), zeros(N+1,3)];

%initialisation
[~, T] = size(data);

X_k_hat = zeros(N+4,T);
theta_0 = 1*[100;0.5e-7;1e-7];
X_k_hat(:,1) = [data(:,1); theta_0];

P_k = [eye(N+2), zeros(N+2,2); zeros(2,N+2), 1e3*eye(2)];
%%
%Loop

for i = 2:T

    %Observtaion
    Y_k = data(:,i);
    
    %Prediction step
    X_ = f(X_k_hat(:,i-1), N, V_vec(i), deltaX, deltaT, eps, Tank, Power(:,i), HeatElem, T_in);
    
    %Error covariance propagation
    F = F_mat(X_k_hat(:,i-1), N, V_vec(i), deltaX, deltaT, eps, Tank, Power(:,i), HeatElem);
    P_k_ = F* P_k* (F.') + Q;
    %Kalman Gain matrix
    K_k = P_k_ * H_' / (H_ * P_k_ * H_' + R);
    
    %Update step
    Y_ = Y_k - H_*X_;
    X_k_hat(:,i) = X_ + K_k*Y_;
    P_k = (eye(N+4) - K_k*H_)*P_k_;

end
toc
%%
Theta = X_k_hat(N+2:N+4,:);