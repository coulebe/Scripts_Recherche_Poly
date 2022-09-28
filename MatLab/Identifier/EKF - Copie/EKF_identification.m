clear ;
data = load('Donnee_sim_Null_Power_VarFlow'); %Pour l'instant, j'importe manuellement les données de simulation
% data = data + normrnd(0, 0.1, size(data)); % Add some noise
%%
tic
[N,T] = size(data.data);
deltaX = data.Layers(2) - data.Layers(1);
deltaT = data.Time(2) - data.Time(1);

%Variance matrices
W = [eye(N+1)   zeros(N+1,1);
    zeros(1,N+1)    0*1e-7]; %Process noise variance
R = 0.1*eye(N+1,N+1); %measurements noise variance 
%Observation matrix/H
H = [eye(N+1), zeros(N+1,1)];
% H_ = zeros(N+1, N+4);
% for i = 1:HeatElem.N
%     pos = HELayer(N, deltaX, HeatElem.Thermos(i));
%     H_(pos, pos) = 1;
% end
%initialisation

X_k_hat = zeros(N+2,T);
V_0 = 0;
X_k_hat(:,1) = [data.data(:,1); data.T_amb; V_0];

P_k = 1*[eye(N+1), zeros(N+1,1); zeros(1,N+1), 1e4];
%%
%Loop

for i = 2:T

    %Observtaion
    Y_k = [data.data(:,i); data.T_amb];
    
    %Prediction step
    X_ = f(X_k_hat(:,i-1), N, deltaX, deltaT, data.eps, data.Tank, data.Q_mat(:,i), data.T_in);
    
    %Error covariance propagation
    F = F_mat(X_k_hat(:,i-1), N, deltaX, deltaT, data.Q_mat(:,i), data.eps, data.Tank);
    P_k_ = F* P_k* (F.') + W;
    %Kalman Gain matrix
    K_k = P_k_ * H' / (H * P_k_ * H' + R);
    
    %Update step
    Y_ = Y_k - H*X_;
    X_k_hat(:,i) = X_ + K_k*Y_;
    P_k = (eye(N+2) - K_k*H)*P_k_;

end
toc
%%
V = X_k_hat(N+2,:);

%%
%Result plot
figure()
plot(data.Time, V);
hold on;
grid on;
plot(data.Time, data.V_vec);
legend('V estimated', 'V')
xlabel('Time $(h)$','Interpreter','Latex','FontSize',12)
ylabel('Draw Rate (m/s)')
%%
for i=1:10
    figure()
    plot(data.Time, data.data(i,:));
    hold on;
    grid on;
    plot(data.Time, X_k_hat(i,:));
    legend(strcat('T_',num2str(i) ,' estimated'), strcat('T_',num2str(i)))
    xlabel('Time $(h)$','Interpreter','Latex','FontSize',12')
    ylabel('°C')
end