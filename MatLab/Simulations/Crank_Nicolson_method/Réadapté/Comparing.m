clear 
tic
%%
%System parameters
%Tank
Tank = struct('Vol', {}, 'H', {}, 'Cv', {},'Rho', {}, 'Dc', {}, 'UL', {}, 'UL_', {}, 'A', {});
Tank(1).H = 1.12; %m %Tank Height 
Tank(1).A = 0.1; %m²%Tank cross-sectionnal area
Tank(1).Vol = Tank.H * Tank.A; %m^3 %Tank Volume
Tank(1).Cv = 4181.3; %%Heat capacity of water 
Tank(1).Rho = 1e3; %kg/m^3 %Density of water
Tank(1).Dc = 1.8/(Tank.Rho*Tank.Cv); %m²/s % Thermal diffusity coefficient
% Tank(1).UL = 4*1.5*sqrt(pi/Tank.A)/(Tank.Rho*Tank.Cv*n); %s^-1 %Thermal losses coefficient 
% Tank(1).UL_ = 1.2382e-6; %s^-1 %Thermal losses coefficient on boundaries 
%Heating Elements
HeatElem = struct('n_eff', {}, 'Power', {}, 'Positions', {}, 'Thermos', {},'N', {});
HeatElem(1).n_eff = 0.95; %Efficiency of the heating elements
HeatElem(1).Power = 0*6e3; %Watt % Electrical power delivered by each heating element
HeatElem(1).Positions = [0.2; 0.72];
HeatElem(1).Thermos = [0.4; 0.92];
HeatElem(1).N = 2;
%%
%Schedule of water drawing
% Draw_Tab = [ %% Draw_start(h) Draw_Duration(min) Draw_Debit(l/min)
%     2.5   40  3 ;
%     5 15  6;
%     ];
Draw_Tab = readmatrix('Draw_tab_var.csv');
Draw_Tab = Draw_Tab(:,2:end);
%Simulation parameters
%Time
deltaT = 5; %s
sim_time = 5;%h
%Space
n = 10;
%Conditions
T_tank = 70; %Initial temperature in the tank 
T_amb = 25; %Ambiant temparature
T_in = 25; %Inlet water temperature
T_target = 70; 
eps = 150;
%%
Tank(1).UL = 4*1.5*sqrt(pi/Tank.A)/(Tank.Rho*Tank.Cv*n); %s^-1 %Thermal losses coefficient 
Tank(1).UL_ = Tank.UL*(1 + (n/(4*Tank.H)) * sqrt(Tank.A/pi) ); %s^-1 %Thermal losses coefficient on boundaries 
%%
[tsol_CN, xVector_CN, sol_CN] = CN_Meth(Tank, HeatElem, Draw_Tab, deltaT, sim_time,n, T_tank, T_amb, T_in, T_target, eps  );
addpath("D:\Alfred\Cours\Projet_Recherche\Poly\Scripts_Recherche_Poly\MatLab\Simulations\PDEPE_method\Readapte");
[tsol_PDEPE, xVector_PDEPE, sol_PDEPE] = PDEPE_Meth(Tank, HeatElem, Draw_Tab, deltaT, sim_time,n, T_tank, T_amb, T_in, T_target, eps  );


%% Figure de la temp rature spatiale en fonction du temps_CN
figure();
surf(tsol_CN/3600,xVector_CN,sol_CN);
xlim([0 sim_time]);
ylim([0 max(xVector_CN)]);
shading interp;
colormap(jet(300));
rotate3d on;
% caxis([0, 100]);
hc=colorbar();
title(hc,'$^{\circ}$ C','Interpreter','Latex');
%view(45,30);
title("Résultat de l'équation par Crank-Nicholson",'FontSize',12');
ylabel('Distance $x$ (m)','Interpreter','Latex','FontSize',12');
xlabel('Time $t$ (h)','Interpreter','Latex','FontSize',12');
zlabel('Solution $u(x,t)$','Interpreter','Latex','FontSize',12');
%%
figure();
surf(tsol_PDEPE/3600,xVector_PDEPE,sol_PDEPE);
xlim([0 sim_time]);
ylim([0 max(xVector_PDEPE)]);
shading interp;
colormap(jet(300));
rotate3d on;
% caxis([0, 100]);
hc=colorbar();
title(hc,'$^{\circ}$ C','Interpreter','Latex');
%view(45,30);
title("Résultat de l'équation par PDEPE",'FontSize',12');
ylabel('Distance $x$ (m)','Interpreter','Latex','FontSize',12');
xlabel('Time $t$ (h)','Interpreter','Latex','FontSize',12');
zlabel('Solution $u(x,t)$','Interpreter','Latex','FontSize',12');
%% Figure de la temp rature spatiale en fonction du temps
figure ()
for i=1:n
    subplot(n,1,n-i+1);
%     figure ();

    if (i == 3) || (i == 7)
        plot(tsol_CN/3600,sol_CN(i,:),'r');
        hold on;
        plot(tsol_PDEPE/3600,sol_PDEPE(i,:),'g');
    else
        plot(tsol_CN/3600,sol_CN(i,:));
        hold on;
        plot(tsol_PDEPE/3600,sol_PDEPE(i,:));

    end
    grid on;
    title(strcat('Couche n = ',num2str(i)))
    legend('Crank-Nicolson', 'PDEPE');
    ylabel('Temp $^{\circ}$ C','Interpreter','Latex','FontSize',9');
    xlim([0 sim_time])
%     ylim([15 65])
end

%%
toc