%% Définition des variables
deltaT = 1; %s
eps = 150;
n = 10; % Nombre de couches
H = 1.19; % Hauteur du tank (m)
debit = 3; %6; % L/min
drawDuration = 30; %3.5; % minutes
V_tank = 112; %L
debit_SI = debit*H/(V_tank*60); %m/s
deltaX = H/n; %m
a = 1*0.14e-6; %m2/s
k = 1*6.3588e-7; %l/s
k_ = 1*1.2382e-6; %l/s
pos1 = 3;%4; %Position de l' lment chauffant au fond du chauffe-eau
pos2 = 7;%10; %Position de l' lment chauffant au sommet du chauffe-eau
on = 1;
off = 0;
m_tank = V_tank; %kg
m_i = m_tank/n;
n_eff = 0.95; % Efficiency = 95%
P_ele = 1*6000; %watt
c = 4185.5; %Water specific heat capacity
T_tank = 25; %Temperature du tank
T_Target = 60; % Target temperature
T_Simulation_hours = 5;% Hours
Simulation_Time = 3600*T_Simulation_hours; % seconds = 30 minutes
Max_Simulation_Count = Simulation_Time/deltaT + 1; % seconds = 30 minutes
%% Calcul de la temp rature du apr s le d but de l'op ration
count = 0;
count2 = 1;
T_amb = 25;
T_out = 25;
Tm = [T_amb;ones(n-2,1)*T_tank;T_out];
V = 0;
Power = zeros(2,Max_Simulation_Count);
PowerTotal = zeros(1,Max_Simulation_Count);
Time = zeros(1,Max_Simulation_Count);
data = zeros(n,Max_Simulation_Count);
E = 0;
offset = 0;
while (count < Max_Simulation_Count)
    if (count > (0*60*60/deltaT) && count <(0*60*60/deltaT+drawDuration*60/deltaT))
        V = debit_SI*0;
    elseif (count > (1*60*60/deltaT) && count <(1*60*60/deltaT+1*drawDuration*60/deltaT))
        V = debit_SI*0;
    elseif (count > (2.5*60*60/deltaT) && count <(2.5*60*60/deltaT+0.5*drawDuration*60/deltaT))
        V = debit_SI*1;
    elseif (count > (3*60*60/deltaT) && count <(3*60*60/deltaT+1*drawDuration*60/deltaT))
        V = debit_SI*0;
    elseif (count > (4*60*60/deltaT) && count <(4*60*60/deltaT+1*drawDuration*60/deltaT))
        V = debit_SI*0;
    elseif (count > (5*60*60/deltaT) && count <(5*60*60/deltaT+1*drawDuration*60/deltaT))
        V = debit_SI*0;
    else
        V = 0;
    end
    [heatState1,heatState2] = getPowerState(Tm,pos1,pos2,T_Target,offset);

    [e,d,q1,q2,q3,q4,q5,q6] = gen_q(deltaX,V,deltaT,eps);
    [Z1, Z2, Z3] = genMatrix(n, e, d, deltaT, q1, q2, q3, q4, q5, q6, heatState1, heatState2, pos1,pos2);
    Am = Z1\Z2;
    Bm = Z1\Z3;
    Tm1 = Am*Tm+Bm;
    Tm = Tm1;
    if V == 0
        Tm(1) = Tm(2);
        Tm(n) = Tm(n-1);
    else
        Tm(1) = 25;
        Tm(n) = Tm(n-1);
    end
    count = count+1;
    Power(:,count) = [heatState1;heatState2];
    PowerTotal(:,count) = heatState1 + heatState2;
    Time(:,count) = count*deltaT/3600;
    data(:,count) = Tm;
    E = E + P_ele*(heatState1+heatState2)*deltaT/3600/1000;
end
%% Figure de la temp rature dans le r servoir
figure (3)
for i=1:(n-2)
    subplot(n-2,1,i);
    if (i == 3)
        plot(Time,data(n-i,:));%,'r', Time, Power(2,:)*20,'-k');
    elseif (i == 7)
        plot(Time,data(n-i,:));%,'r', Time, Power(1,:)*20,'-k');
    else
        plot(Time,data(n-i,:));
    end
    legend(strcat('Couche n = ',num2str(n-i')))
    ylabel('Temp $^{\circ}$ C','Interpreter','Latex','FontSize',9');
    xlim([0 T_Simulation_hours])
    ylim([15 65])
end
xlabel('Time $(h)$','Interpreter','Latex','FontSize',12');
%% Figure combin e;
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
%% Figure en 3D de l' volution de la temp rature
figure(5)
Layers = linspace(1,10,10);
surf(Time,Layers,data)
xlim([0 T_Simulation_hours])
ylim([2 9])
shading interp;
colormap(jet(30000));
caxis([0, 100]);
rotate3d on;
xlabel('Time $(h)$','Interpreter','Latex')
ylabel('Couche $n$','Interpreter','Latex')
zlabel('Temp $^{\circ}$ C','Interpreter','Latex')
hc=colorbar();
title(hc,'$^{\circ}$ C','Interpreter','Latex');