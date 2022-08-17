clear 
tic
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
%%
%Schedule of water drawing
Draw_Tab = [ %% Draw_start(h) Draw_Duration(min) Draw_Debit(l/min)
    2.5   40  3 ;
    5 15  6;
    ];
%Simulation parameters
%Time
deltaT = 5; %s
sim_time = 10;%h
%Space
n = 10;
%Conditions
T_tank = 60; %Initial temperature in the tank 
T_amb = 25; %Ambiant temparature
T_in = 25; %Inlet water temperature
T_target = 60; 
eps = 150;
%%
[tsol, xVector, sol] = PDEPE_Meth(Tank, HeatElem, Draw_Tab, deltaT, sim_time,n, T_tank, T_amb, T_in, T_target, eps  );
%%
%% Figure de la temp rature spatiale en fonction du temps
figure();
surf(tsol/3600,xVector,sol.');
xlim([0 sim_time]);
ylim([0 max(xVector)]);
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
        plot(tsol/3600,sol(:,i),'r');
    elseif (i == 7)
        plot(tsol/3600,sol(:,i),'r');
    else
        plot(tsol/3600,sol(:,i));
    end
    hold on;
    grid on;
    legend(strcat('Couche n = ',num2str(n-i+1)))
    ylabel('Temp $^{\circ}$ C','Interpreter','Latex','FontSize',9');
    xlim([0 sim_time])
%     ylim([15 65])
end

%%
toc