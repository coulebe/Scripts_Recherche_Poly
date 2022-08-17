function [tsol, xVector, sol] = PDEPE_Meth(Tank, HeatElem, Draw_Tab, deltaT, sim_time,N_layer, T_tank, T_amb, T_in, T_target, eps  )
    %This function resolve the pdepe of the temperature in a water tank by
    %using pdepe method
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
    %sol: Solution matrice(NxM)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Time variable
    t_max = 3600*sim_time;
    nb_point = t_max/deltaT + 1;
    %Space
    xVector = linspace(0,Tank.H,N_layer); %Discretization of space
    heatState = zeros(HeatElem.N,1);
    %%
    %Simulation
    initial = [T_tank;ones(N_layer-2,1)*T_tank;T_tank]';
    currentTemp = initial;
    t0 = 0; step = deltaT; tf = step; tn = 3;
    t = linspace(t0,tf,tn);
    sol= initial;
    option = odeset('RelTol',1e-2,'AbsTol',1e-5);

    for z=1:nb_point-1
        pdefunc = @(x,t,T,dTdx) pdefun_(x,t,T,dTdx, Tank, HeatElem, T_amb, heatState, eps, N_layer, Draw_Tab);
        bcfunc = @(xl,ul,xr,ur,t) bcfun_(xl,ul,xr,ur,t, Tank, T_in, currentTemp, N_layer, Draw_Tab);
        icfunc = @(x) icfun_(x,xVector, initial);
        u=pdepe(0,pdefunc,icfunc,bcfunc,xVector,t,option);
        currentTemp = u(end,:);
        sol=[sol; u(end,:)];
        initial = u(end,:);
        heatState = PowerState_(sol(end, :),HeatElem,T_target,xVector);
        t0=t0+step;
        tf=tf+step;
        t=linspace(t0,tf,tn);
    end
    tsol = linspace(0, tf, nb_point);
    sol = sol.';
    

end
