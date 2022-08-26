function [tsol, xVector, sol] = CN_Meth(Tank, HeatElem, Draw_Tab, deltaT, sim_time,N_layer, T_tank, T_amb, T_in, T_target, eps  )
    %This function resolve the pdepe of the temperature in a water tank by
    %using Crank-Nicolson method
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Parameters%%%%%%%%%%%%%%%%%%%%%
    %
    %Tank: A structure containing the tank caracteristics(
    %       Vol: Volume(m^3), H: Heigth(m), Cv: Heat capacity of the fluid(J/(kg °K)), 
    %       Rho: Fluid density(kg/m^3), Dc: Thermal diffusity coefficient(m²/s), 
    %       UL: Thermal losses coefficient(s^-1), UL_ = Thermal losses coefficient the boudaries)
    %
    %HeatElem: A strauture containing the Heating Elements caracteristics(
    %       n_eff: Efficiency of the heating elements
    %       Power: Electrical power delivred by each element(Watt),
    %       Positions: Tab containing the position of the heating element(Ascending)
    %       Thermos: Tab containing the position of the the thermostat used for heating elements control(Ascending)
    %       N: Number of heating elements
    %       Positions and thermos must have same length = N)
    %
    %Draw_tab : Tab containing the draw planned during the simulation. Each
    %line contain the hour when  the draw start(h), its duration(min), and
    % its debit(l/min). The tab is ascending given the hour of start.
    %
    %deltaT : Time step(s)
    %sim_time: Duration of the simulation(h)
    %N_layers: number of layers that discretize the space
    %T_tank: initial temperature in the tank(°C)
    %T_amb: Ambient temperaure(°C)
    %T_in: inlet fluid temperaure(°C)
    %T_target: Target temperature for the heating element control(°C)
    %eps coefficient that replace the thermal expansion coefficient
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Output%%%%%%%%%%%%%%%%%%%%%%%%%
    %tsol: solution time vector(1xM)(sec)
    %xVector: Space vector (1xN_layers)(m)
    %sol: Solution matrice(MxN)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Time variable
    Max_Simulation_Count = sim_time*3600/deltaT + 1;
    %Space
    xVector = linspace(0,Tank.H,N_layer); %Discretization of space
    %Simulation parameters
    %Space
    deltaX = Tank.H/N_layer; %m

    %% Simulation
    count = 0;
    Tm = [T_amb;ones(N_layer-1,1)*T_tank;T_tank];

    tsol = zeros(1,Max_Simulation_Count);
    sol = zeros(N_layer,Max_Simulation_Count);
    while (count < Max_Simulation_Count)

        if(~isempty(Draw_Tab))
            if(count > (Draw_Tab(1,1) * 60*60/deltaT))
                V = Draw_Tab(1,3)*1e-3*Tank.H/(Tank.Vol*60); %m/s

                if(count > (Draw_Tab(1,1) * 60*60/deltaT + Draw_Tab(1,2)*60/deltaT))
                    V = 0;
                    Draw_Tab(1,:) = [];
                end
            else
                V = 0;
            end
        else
            V = 0;
        end

        heatState = PowerState(Tm,HeatElem,T_target,deltaX, N_layer);


        [Z1, Z2, Z3] = Matrix_(N_layer, V,deltaX,deltaT,eps, Tank, heatState, HeatElem);


        Am = Z1\Z2;
        Bm = Z1\Z3;
        Tm = Am*Tm+Bm;
        if V == 0
            Tm(1) = Tm(2);
            Tm(N_layer) = Tm(N_layer-1);
        else
            Tm(1) = T_in;
            Tm(N_layer) = Tm(N_layer-1);
        end
        count = count+1;
        tsol(:,count) = count*deltaT;
        sol(:,count) = Tm(1:end-1);
    end

end
