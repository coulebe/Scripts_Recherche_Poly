%% Définition des variables
clear
tic
%%
%System Parameters
%Tank
Tank = struct('Vol', {}, 'H', {}, 'Cv', {},'Rho', {}, 'Dc', {}, 'UL', {}, 'UL_', {});
Tank(1).Vol = 112e-3; %m^3 %Tank Volume
Tank(1).H = 1.19; %m %Tank Height 
Tank(1).Cv = 4185.5; %J/(kg T)%Heat capacity of water 
Tank(1).Rho = 1e3; %kg/m^3 %Density of water
Tank(1).Dc = 0.14e-6; %m²/s % Thermal diffusity coefficient
Tank(1).UL = 6.3588e-7; %s^-1 %Thermal losses coefficient 
Tank(1).UL_ = 1.2382e-6; %s^-1 %Thermal losses coefficient on boundaries 

%Heating Element
HeatElem = struct('n_eff', {}, 'Power', {}, 'Positions',{}, 'Thermos', {}, 'N', {});
HeatElem(1).n_eff = 0.95; %Efficiency of the heating elements
HeatElem(1).Power = 6e3; %Watt % Electrical power delivered by each heating element
HeatElem(1).Positions = [0.2975; 0.7735];
HeatElem(1).Thermos = [0.2975; 0.7735];
HeatElem(1).N = 2;
%%
%Schedule of water drawing
Draw_Tab = [ %% Draw_start(h) Draw_Duration(min) Draw_Debit(l/min)
%     2.5   40  3 ;
%     5 15  6;
    ];
%%
%Simulation parameters
%Time
T_Simulation_hours = 10;% Hours
deltaT = 1; %s
Simulation_Time = 3600*T_Simulation_hours; 
Max_Simulation_Count = Simulation_Time/deltaT + 1; % seconds = 30 minutes
%Space
N = 10; % Number of layers
deltaX = Tank.H/N; %m
%Conditions
T_tank = 25; %Tank's initial temperature
T_Target = 60; % Target temperature



%Simulation parameters
%For Later,  Think about convergence and stability conditions of the Crank-Nicolson's 
%Scheme if needed


%
eps = 150;
HeatElem.N = 2; %Number of heating elements
% HE_Positions_index = HELayer(N, deltaX, HeatElem.Positions); 
% HE_Thermos_index = HELayer(N, deltaX, HeatElem.Thermos);
%% Simulation
count = 0;
Draw = Draw_Tab;
T_amb = 25;
T_in = 25;
Tm = [ones(N,1)*T_tank;T_amb];
V = 0;

Power = zeros(HeatElem.N,Max_Simulation_Count);
PowerTotal = zeros(1,Max_Simulation_Count);
Time = zeros(1,Max_Simulation_Count);
data = zeros(N,Max_Simulation_Count);
E = 0;
Q_mat = zeros(N,Max_Simulation_Count);
while (count < Max_Simulation_Count)

    if(~isempty(Draw))
        if(count > (Draw(1,1) * 60*60/deltaT))
            V = Draw(1,3)*1e-3*Tank.H/(Tank.Vol*60); %m/s
            if(count > (Draw(1,1) * 60*60/deltaT + Draw(1,2)*60/deltaT))
                V = 0;
                Draw(1,:) = [];
            end
        else
            V = 0;
        end
    else
        V = 0;
    end
%Check if a heatint element need to be activated 
    heatState = 1*PowerState(Tm,HeatElem,T_Target,deltaX, N);


    [Z1, Z2, Z3] = Matrix_(N, V,deltaX,deltaT,eps, Tank, heatState, HeatElem);


    Am = Z1\Z2;
    Bm = Z1\Z3;
    Tm =  Am *Tm + Bm ;
    %Boundary conditions
    if V == 0
        Tm(1) = (4*Tm(2)-Tm(3))/3;
    else
        Tm(1) = T_in;
%         Tm(N) = (4*Tm(N-1) - Tm(N-2))/3;
    end
    Tm(N) = (4*Tm(N-1) - Tm(N-2))/3;
    count = count+1;
    Power(:,count) = heatState;
    Q_mat(:,count) = Z3(1:end-1);
    PowerTotal(:,count) = sum(heatState);
    Time(:,count) = count*deltaT/3600;
    data(:,count) = Tm(1:end-1);
    E = E + HeatElem.Power*(sum(heatState))*deltaT/3600/1000; %KJ
end
%% Figure de la température dans le réservoir
% figure ()
% for i=1:(N-2)
%     subplot(N-2,1,i);
%     if (i == 3)
%         plot(Time,data(N-i,:));%,'r', Time, Power(2,:)*20,'-k');
%     elseif (i == 7)
%         plot(Time,data(N-i,:));%,'r', Time, Power(1,:)*20,'-k');
%     else
%         plot(Time,data(N-i,:));
%     end
%     legend(strcat('Couche n = ',num2str(N-i')))
%     ylabel('Temp $^{\circ}$ C','Interpreter','Latex','FontSize',9');
%     xlim([0 T_Simulation_hours])
%     ylim([15 65])
% end
% xlabel('Time $(h)$','Interpreter','Latex','FontSize',12');
%% Figure combinée;
% figure (10)
% for i=1:(n-6)
%     subplot(n-6,1,i);
%     plot(Time,data(n-i,:));
%     hold on
%     plot(pdeTime,pdepeU(:,n-i));
%     legend(strcat('Couche n = ',num2str(n-i')))
%     ylabel('Temp $^{\circ}$ C','Interpreter','Latex','FontSize',9');
%     xlim([0 T_Simulation_hours])
%     ylim([15 65])
% end
% xlabel('Time $(h)$','Interpreter','Latex','FontSize',12');
%% Figure en 3D de l'évolution de la température
figure()
Layers = linspace(0,Tank.H,N);
surf(Time,Layers,data)
xlim([0 T_Simulation_hours])
% ylim([2 9])
shading interp;
colormap(jet(30000));
caxis([0, 100]);
rotate3d on;
xlabel('Time $(h)$','Interpreter','Latex')
ylabel('Couche $n$','Interpreter','Latex')
zlabel('Temp $^{\circ}$ C','Interpreter','Latex')
hc=colorbar();
title(hc,'$^{\circ}$ C','Interpreter','Latex');
%%
toc