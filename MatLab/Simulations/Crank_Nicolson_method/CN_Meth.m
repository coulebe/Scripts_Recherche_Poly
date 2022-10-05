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


    tsol = zeros(1,Max_Simulation_Count);
    sol = zeros(N_layer+1,Max_Simulation_Count);
    sol(:, 1) = [ones(N_layer,1)*T_tank;T_amb];
    for count = 2:Max_Simulation_Count

        V = drawRate(count*deltaT, Tank, Draw_Tab);
        heatState = PowerState(sol(:, count-1),HeatElem,T_target,deltaX, N_layer);


        [Z1, Z2, Z3] = Matrix(N_layer, V,deltaX,deltaT,eps, Tank, heatState, HeatElem);


        Am = Z1\Z2;
        Bm = Z1\Z3;
        sol(:, count) = Am*sol(:, count-1)+Bm;
        if V == 0
            sol(1, count) = (4*sol(2, count)-sol(3, count))/3;
        else
            sol(1, count) = T_in;
        end
        tsol(:,count) = count*deltaT;
    end
    sol = sol(1:end-1, :);
end
