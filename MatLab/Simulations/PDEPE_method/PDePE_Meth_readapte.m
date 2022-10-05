clear 
tic
%%
%Export
% global Draw_Tab;
% global iterator;
%%
%System parameters
%Tank
Tank = struct('Vol', {}, 'H', {}, 'Cv', {},'Rho', {}, 'Dc', {}, 'UL', {}, 'UL_', {});
Tank(1).Vol = 112e-3; %m^3 %Tank Volume
Tank(1).H = 1.19; %m %Tank Height 
Tank(1).Cv = 4185.5; %%Heat capacity of water 
Tank(1).Rho = 1e3; %kg/m^3 %Density of water
Tank(1).Dc = 0.14e-6; %m²/s % Thermal diffusity coefficient
Tank(1).UL = 6.3588e-7; %s^-1 %Thermal losses coefficient 
Tank(1).UL_ = 1.2382e-6; %s^-1 %Thermal losses coefficient on boundaries 
%Heating Elements
HeatElem = struct('n_eff', {}, 'Power', {}, 'Positions', {}, 'Thermos', {},'N', {});
HeatElem(1).n_eff = 0.95; %Efficiency of the heating elements
HeatElem(1).Power = 6e3; %Watt % Electrical power delivered by each heating element
HeatElem(1).Positions = [0.2975; 0.7735];
HeatElem(1).Thermos = [0.2975; 0.7735];
HeatElem(1).N = 2;
% iterator = 0;
%%
%Schedule of water drawing
Draw_Tab = [ %% Draw_start(h) Draw_Duration(min) Draw_Debit(l/min)
%     2.5   40  3 ;
%     5 15  6;
    ];
%%
%Simulation parameters
%Time
deltaT = 5; %s
sim_time = 5;%h
t_max = 3600*sim_time;
nb_point = t_max/deltaT + 1;
%Space
n = 100;
xVector = linspace(0,Tank.H,n); %Discretization of space
%Conditions
T_tank = 25; %Initial temperature in the tank 
T_amb = 25; %Ambiant temparature
T_in = 25; %Inlet water temperature
T_target = 60; 
eps = 150;
heatState = zeros(HeatElem.N,1);
%%
%Simulation
initial = [T_tank;ones(n-2,1)*T_tank;T_tank]';
currentTemp = initial;
t0 = 0; step = deltaT; tf = step; tn = 3;
t = linspace(t0,tf,tn);
% u = zeros(nb_point, n);
% u(1,:)= initial;
u = initial;
option = odeset('RelTol',1e-2,'AbsTol',1e-5);
% HS = [t0; heatState];

%
for z=1:nb_point-1
    pdefunc = @(x,t,T,dTdx) pdefun_(x,t,T,dTdx, Tank, HeatElem, T_amb, heatState, eps, n, Draw_Tab);
    bcfunc = @(xl,ul,xr,ur,t) bcfun_(xl,ul,xr,ur,t, Tank, T_in, currentTemp, n, Draw_Tab);
    icfunc = @(x) icfun_(x,xVector, initial);
    sol=pdepe(0,pdefunc,icfunc,bcfunc,xVector,t,option);
    currentTemp = sol(end,:);
    u=[u; sol(end,:)];
%     u(z+1, :) = sol(end, :);
    initial = sol(end,:);
    heatState = PowerState_(u(end, :),HeatElem,T_target,xVector);

%     H = [tf;heatState];
%     HS = [HS, H];

    t0=t0+step;
    tf=tf+step;
    t=linspace(t0,tf,tn);
end

%% Figure de la temp rature spatiale en fonction du temps
figure();
tg = linspace(0,tf,nb_point);
surf(tg/3600,xVector,u.');
xlim([0 sim_time]);
ylim([0 max(xVector)]);
% surf(x,t,u);
shading interp;
colormap(jet(300));
rotate3d on;
caxis([0, 100]);
hc=colorbar();
title(hc,'$^{\circ}$ C','Interpreter','Latex');
%view(45,30);
%title('Reaction-Advection-Diffusion Equation','FontSize',12');
ylabel('Distance $x$ (m)','Interpreter','Latex','FontSize',12');
xlabel('Time $t$ (h)','Interpreter','Latex','FontSize',12');
zlabel('Solution $u(x,t)$','Interpreter','Latex','FontSize',12');
%% Figure de la temp rature spatiale en fonction du temps
figure ()
for i=1:n
    subplot(n,1,i);
    if (i == 3)
        plot(tg/3600,u(:,i),'r');
    elseif (i == 7)
        plot(tg/3600,u(:,i),'r');
    else
        plot(tg/3600,u(:,i));
    end
    hold on;
    grid on;
    legend(strcat('Couche n = ',num2str(n-i+1)))
    ylabel('Temp $^{\circ}$ C','Interpreter','Latex','FontSize',9');
    xlim([0 sim_time])
%     ylim([15 65])
end
xlabel('Time $(h)$','Interpreter','Latex','FontSize',12');
%%
% %Figure montrant la qu'il n'y a pas d'activation simultanée des éléments
% %chauffants
% figure()
% for i =HeatElem.N:-1:1  
%     subplot(HeatElem.N,1,HeatElem.N-i+1);
%     stem(HS(1,:)/3600, HS(i+1,:));
%     hold on;
%     grid on;
%     legend(strcat('Element chauffant', num2str(HeatElem.N-i+1)));
%     ylim([0 1]);
% %     xlim([0 500]);
%     ylabel('Etat')
%     
% end
% 
% figure()
% stem(HS(1,:)/3600, HS(2,:) & HS(3,:))
% hold on;
%     grid on;
%     title("Vérification de simultanéité d'activation");
%     ylim([0 1]);
% %     xlim([0 10]);
%     ylabel('Etat')
%%
toc


