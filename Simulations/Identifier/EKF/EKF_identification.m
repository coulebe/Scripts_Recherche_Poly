clear ;
load('Donnee_sim_Null_Flow.mat');
% addpath("D:\Alfred\Cours\Projet_Recherche\Poly\Scripts_Recherche_Poly\Simulations\Crank_Nicolson_method\Réadapté");
data = [data; T_amb*ones(1,Max_Simulation_Count)];
% data = data + normrnd(0, 0.1, size(data)); % Add some noise
%%
tic
%Variance matrices
Q = [eye(N+1)   zeros(N+1,3);
    zeros(3,N+1)    diag([1, 1e-7, 1e-7])]; %Process noise variance
R = eye(N+1,N+1); %measurements noise variance
%Observation matrix/H
H_ = [eye(N+1), zeros(N+1,3)];
% H_ = zeros(N+1, N+4);
% H_(1,1) = 1;
% H_(5,5) = 1;
% H_(N,N) = 1;
% H_(N+1,N+1) = 1;

%initialisation
[~, T] = size(data);

X_k_hat = zeros(N+4,T);
theta_0 = [100;1e-7;1e-7];
X_k_hat(:,1) = [mean(data(:,1))*ones(size(data(:,1))); theta_0];

P_k = [eye(N+2), zeros(N+2,2); zeros(2,N+2), 1e4*eye(2)];
%%
%Loop

for i = 2:1000

    %Observtaion
    Y_k = data(:,i);
    
    %Prediction step
    X_ = f(X_k_hat(:,i-1), N, V_vec(i), deltaX, deltaT, Q_mat(:,i), T_in);
    
    %Error covariance propagation
    F = F_mat(X_k_hat(:,i-1), N, V_vec(i), deltaX, deltaT, Q_mat(:,i));
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
Theta = X_k_hat(N+2:end,:);

%%
%Result plot
figure()
plot(Time, Theta(1,:));
hold on;
grid on;
yline(eps);
xlim([Time(1)  Time(1000)]);
legend('\alpha estimated', '\alpha')
%
figure()
plot(Time, Theta(2,:));
hold on;
grid on;
yline(Tank.Dc);
xlim([Time(1)  Time(1000)]);
legend('Dc estimated', 'Dc')
%
figure()
plot(Time, Theta(3,:));
hold on;
grid on;
yline(Tank.UL);
xlim([Time(1)  Time(1000)]);
legend('UL estimated', 'UL')