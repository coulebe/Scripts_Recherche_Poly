clear
load('Data_Sim_passive.mat');
%%
%Updating coefficients
gamma = 0.1; % Temperature upadte
Phi_1 = 1e-7; %Dc update
Phi_2 = 1e-7; %Ul update
%
% data = data;%-T_amb;
%
T_hat = zeros(size(data));
%
Dc_hat = zeros(size(Time));
%
UL_hat = zeros(size(Time));
%
T_xx = spaceDiff2(data, deltaX);
%%
%Initialisation
T_hat(:,1) = data(:,1);

Dc_hat(1) = 1e-6;

UL_hat(1) = 1e-6;
%%
%Update

for k = 2:count
    alpha = deltaT*(gamma * numFuncNorm(data(:,k-1), Layers)).^2;
    Err = data(:,k-1) - T_hat(:,k-1);
    
    %Dc
    x = trapz(Layers, eps * T_xx(:,k-1).*Err);
    Dc_hat(k) = Dc_hat(k-1) + deltaT * Phi_1 * proj(x, Dc_hat(k-1), 0, (deltaX^2)/(2*eps*deltaT));
    
    %UL
    x = trapz(Layers, (data(:,k-1)-T_amb).*Err);
    UL_hat(k) = UL_hat(k-1) + deltaT * Phi_2 * proj(x, UL_hat(k-1), 0, 1);
    
    %T
    A = A_mat(deltaT, deltaX, Dc_hat(k), UL_hat(k), N, alpha);
    
    T_hat(:,k) = A\T_hat(:,k-1) + alpha * A\data(:,k) + deltaT*A\Q_mat(:,k);
    
end
%%
yline(Tank.Dc);
hold on;
grid on;
plot(Time, Dc_hat, '*');
legend('Dc', 'Dc_hat')
