% echo on;
clear 
tic
sim_time = 10;%h
global T_tank;
global T_amb;
global T_in;
global iterator;
global T_target;
global xVector;
global eps;
global n;
global heatState;
global Draw_Tab;
%
Draw_Tab = [ %% Draw_start(h) Draw_Duration(min) Draw_Debit(l/min)
    2.5   40  3 ;
    5 15  6;
    ];
%
global Tank;
Tank = struct('Vol', {}, 'H', {}, 'Cv', {},'Rho', {}, 'Dc', {}, 'UL', {}, 'UL_', {});
Tank(1).Vol = 112e-3; %m^3 %Tank Volume
Tank(1).H = 1.19; %m %Tank Height 
Tank(1).Cv = 4185.5; %J/(kg T)%Heat capacity of water 
Tank(1).Rho = 1e3; %kg/m^3 %Density of water
Tank(1).Dc = 0.14e-6; %m²/s % Thermal diffusity coefficient
Tank(1).UL = 6.3588e-7; %s^-1 %Thermal losses coefficient 
Tank(1).UL_ = 1.2382e-6; %s^-1 %Thermal losses coefficient on boundaries 
%
global HeatElem;
HeatElem = struct('n_eff', {}, 'Power', {}, 'Positions', {}, 'N', {}, 'Positions_index', {});
HeatElem(1).n_eff = 0.95; %Efficiency of the heating elements
HeatElem(1).Power = 6e3; %Watt % Electrical power delivered by each heating element
HeatElem(1).Positions = [0.2975; 0.7735];
HeatElem(1).N = 2;
HeatElem.Positions_index = [3,7];
iterator = 0;
deltaT = 15;

n = 10;
t_max = 3600*sim_time;
nb_point = t_max/deltaT + 1;
t_0 = 0;
x = linspace(0,1.19,n);
xVector = x;
T_tank = 25;
T_amb = 25;
T_in = 25;
T_target = 60;
eps = 150;

global currentTemp;
global initial;
heatState = zeros(HeatElem.N,1);

initial = [T_tank;ones(n-2,1)*T_tank;T_tank]';
currentTemp = initial;
t0 = 0; step = deltaT; tf = step; tn = 3;
t = linspace(t0,tf,tn);
zf = t_max/step;
u=[];
option = odeset('RelTol',1e-2,'AbsTol',1e-5);

%
HS = [t0; heatState];
%
for z=1:nb_point
    sol=pdepe(0,@pdefun2,@icfun,@bcfun,x,t,option);
    currentTemp = smooth(sol(end,:));
%     currentTemp = sol(end,:);
    u=[u; sol(1:end-1,:)];
    initial = smooth(sol(end,:));
%     initial = sol(end,:);
    heatState = PowerState_(u(end, :),HeatElem,T_target);
%     H = [tf; heatState];
%     HS = [HS, H];


    t0=t0+step;
    tf=tf+step;
    t=linspace(t0,tf,tn);
end

%% Figure de la temp rature spatiale en fonction du temps
figure();
tg = linspace(0,tf,length(u));
surf(tg/3600,x,u.');
xlim([0 sim_time]);
ylim([0 max(x)]);
% surf(x,t,u);
shading interp;
colormap(jet(300));
rotate3d on;
caxis([0, 100]);
hc=colorbar();
title(hc,'$^{\circ}$ C','Interpreter','Latex');
%view(45,30);
%title('Reaction-Advection-Diffusion Equation','FontSize',12');
ylabel('Distance $x$','Interpreter','Latex','FontSize',12');
xlabel('Time $t$','Interpreter','Latex','FontSize',12');
zlabel('Solution $u(x,t)$','Interpreter','Latex','FontSize',12');
%% Figure de la temp rature spatiale en fonction du temps
pdeTime = tg/3600;
pdepeU = u;
figure ()
for i=1:(n-2)
    subplot(n-2,1,i);
    if (i == 3)
        plot(pdeTime,pdepeU(:,n-i),'r');
    elseif (i == 7)
        plot(pdeTime,pdepeU(:,n-i),'r');
    else
        plot(pdeTime,pdepeU(:,n-i));
    end
    hold on;
    grid on;
    legend(strcat('Couche n = ',num2str(n-i')))
    ylabel('Temp $^{\circ}$ C','Interpreter','Latex','FontSize',9');
    xlim([0 sim_time])
    ylim([15 65])
end
xlabel('Time $(h)$','Interpreter','Latex','FontSize',12');
%%
%Figure montrant la qu'il n'y a pas d'activation simultanée des éléments
%chauffants
% figure()
% for i =1:HeatElem.N  
%     subplot(HeatElem.N,1,i);
%     stem(HS(1,:), HS(i+1,:));
%     hold on;
%     grid on;
%     legend(strcat('Element chauffant', num2str(i)));
%     ylim([0 1]);
%     xlim([0 500]);
%     ylabel('Etat')
%     
% end

% figure()
% stem(HS(1,:), HS(2,:) & HS(3,:))
% hold on;
%     grid on;
%     title("Vérification de simultanéité d'activation");
%     ylim([0 1]);
% %     xlim([0 10]);
%     ylabel('Etat')
toc


