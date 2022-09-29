clear ;
data = load('Donnee_sim_Q_and_V'); %Pour l'instant, j'importe manuellement les données de simulation
% data = data + normrnd(0, 0.1, size(data)); % Add some noise
HeatElem = struct('n_eff', {}, 'Power', {}, 'Positions',{}, 'Thermos', {}, 'N', {});
HeatElem(1).n_eff = 0.95; %Efficiency of the heating elements
HeatElem(1).Power = 6e3; %Watt % Electrical power delivered by each heating element
HeatElem(1).Positions = [0.2975; 0.7735];
HeatElem(1).Thermos = [0.2975; 0.7735];
HeatElem(1).N = 2;


%%
tic
[N,T] = size(data.data);
deltaX = data.Layers(2) - data.Layers(1);
deltaT = data.Time(2) - data.Time(1);
data.T_in = 25;% Je triche un peu
%Variance matrices
W = eye(N+1)   ; %Process noise variance
R = zeros(HeatElem.N+1); %measurements noise variance 
%Observation matrix/H
% H = [eye(N+1), zeros(N+1,1)];
H = zeros(HeatElem.N+1,N+1);
for i = 1:HeatElem.N
    pos = HELayer(N, deltaX, HeatElem.Thermos(i));
    H(i, pos) = 1;
    R(i,i)  = 0.1;
end
H(end, N+1) = 1;
R(end,end) = 0.1;
%initialisation

X_k_hat = zeros(N+1,T);

X_k_hat(:,1) = [data.data(:,1); data.T_amb];

P_k = eye(N+1);
%%
%Loop

for i = 2:T

    %Observtaion
    Y_k = H*[data.data(:,i); data.T_amb];
    
    %Prediction step
    X_ = f(X_k_hat(:,i-1), N, data.V_vec(i,1), deltaX, deltaT, data.eps, data.Tank, data.Q_mat(:,i), data.T_in);
    
    %Error covariance propagation
    F = F_mat(X_k_hat(:,i-1), N, data.V_vec(i), deltaX, deltaT, data.Q_mat(:,i), data.eps, data.Tank);
    P_k_ = F* P_k* (F.') + W;
    %Kalman Gain matrix
    K_k = P_k_ * H' / (H * P_k_ * H' + R);
    
    %Update step
    Y_ = Y_k - H*X_;
    X_k_hat(:,i) = X_ + K_k*Y_;
    P_k = (eye(N+1) - K_k*H)*P_k_;

end
toc

%%
for i=1:10
    figure()
    plot(data.Time, X_k_hat(i,:), '*');
    hold on;
    grid on;
    plot(data.Time,data.data(i,:));
    legend(strcat('T_',num2str(i) ,' estimated'), strcat('T_',num2str(i)))
    xlabel('Time $(h)$','Interpreter','Latex','FontSize',12')
    ylabel('°C')
end