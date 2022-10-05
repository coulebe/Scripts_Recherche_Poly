clc
clear 
close all
addpath Crank_Nicolson_method
addpath PDEPE_method

%%
%System parameters
load('Tank Without name.mat')
%Tank
% Tank = struct('Vol', {}, 'H', {}, 'Cv', {},'Rho', {}, 'Dc', {}, 'UL', {}, 'UL_', {}, 'A', {});
% Tank(1).H = 1.12; %m %Tank Height 
% Tank(1).A = 0.1; %m²%Tank cross-sectionnal area
% Tank(1).Vol = Tank.H * Tank.A; %m^3 %Tank Volume
% Tank(1).Cv = 4181.3; %%Heat capacity of water 
% Tank(1).Rho = 1e3; %kg/m^3 %Density of water
% Tank(1).Dc = 1.8/(Tank.Rho*Tank.Cv); %m²/s % Thermal diffusity coefficient

%Heating Elements
% HeatElem = struct('n_eff', {}, 'Power', {}, 'Positions', {}, 'Thermos', {},'N', {});
% HeatElem(1).n_eff = 0.95; %Efficiency of the heating elements
% HeatElem(1).Power = 6e3; %Watt % Electrical power delivered by each heating element
% HeatElem(1).Positions = [0.2975; 0.7735];
% HeatElem(1).Thermos = [0.2975; 0.7735];
% HeatElem(1).N = min([size(HeatElem.Positions, 1) size(HeatElem.Thermos, 1)]);
%%
%Schedule of water drawing
% Draw_Tab = [ %% Draw_start(h) Draw_Duration(min) Draw_Debit(l/min)
%     2.5   40  3 ;
%     5 15  6;
%     ];

[file, path] = uigetfile(...
                            {'*.csv', 'Draw Tab(*.csv)'}, ...
                            'Select a file containing a draw tab for the simulation');
file = fullfile(path, file);                   
Draw_Tab = readmatrix(file);

Draw_Tab = Draw_Tab(:,2:end);
%Simulation parameters
%Time
deltaT = 5; %s
sim_time = 5;%h
%Space
n = 100;
%Conditions
T_tank = 25; %Initial temperature in the tank 
T_amb = 25; %Ambiant temparature
T_in = 25; %Inlet water temperature
T_target = 60; 
eps = 150;
%%
% Tank(1).UL = 4*1.5*sqrt(pi/Tank.A)/(Tank.Rho*Tank.Cv); %s^-1 %Thermal losses coefficient 
Tank(1).UL_ = Tank.UL*(1 + (n/(4*Tank.H)) * sqrt(Tank.A/pi) ); %s^-1 %Thermal losses coefficient on boundaries 
%%
fprintf ("Crank-Nicholson method: ")
tic
[tsol_CN, xVector_CN, sol_CN] = CN_Meth(Tank, HeatElem, Draw_Tab, deltaT, sim_time,n, T_tank, T_amb, T_in, T_target, eps );
toc

% addpath("D:\Alfred\Cours\Projet_Recherche\Poly\Scripts_Recherche_Poly\MatLab\Simulations\PDEPE_method");
fprintf ("PDEPE method: ")
tic
[tsol_PDEPE, xVector_PDEPE, sol_PDEPE] = PDEPE_Meth(Tank, HeatElem, Draw_Tab, deltaT, sim_time,n, T_tank, T_amb, T_in, T_target, eps  );
toc

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
% figure ()
% for i=1:n
%     subplot(n,1,n-i+1);
% %     figure ();
%     plot(tsol_CN/3600,abs(sol_CN(i,:)-sol_PDEPE(i,:))*100./sol_PDEPE(i,:));
%     hold on;
% %     plot(tsol_PDEPE/3600,);
% 
% 
%     grid on;
% 
%     legend(strcat('Couche n = ',num2str(i)));
%     xlim([0 sim_time])
%     xlabel('Time $t$ (h)','Interpreter','Latex','FontSize',9');
%     
% end
% sgtitle('Crank-Nicolson-PDEPE relative error(%)','Interpreter','Latex','FontSize',20')
%%
figure();
surf(tsol_PDEPE/3600,xVector_PDEPE,abs(sol_CN - sol_PDEPE));
xlim([0 sim_time]);
ylim([0 max(xVector_PDEPE)]);
shading interp;
colormap(jet(300));
rotate3d on;
% caxis([0, 100]);
hc=colorbar();
title(hc,'$^{\circ}$ C','Interpreter','Latex');
%view(45,30);
title("Crank-Nicolson-PDEPE absolute error",'FontSize',12');
ylabel('Distance $x$ (m)','Interpreter','Latex','FontSize',12');
xlabel('Time $t$ (h)','Interpreter','Latex','FontSize',12');
zlabel('Solution $u(x,t)$','Interpreter','Latex','FontSize',12');