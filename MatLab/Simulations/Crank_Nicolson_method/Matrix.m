function [Z1, Z2, Z3] = Matrix(N_Layers, V,deltaX,deltaT,eps, Tank, heatState, HE)
    %We suppose the heatState has the same size than HE.Positions_index
    m_i = Tank.Vol*Tank.Rho/N_Layers;

    
    e = eps * Tank.Dc/(2 * (deltaX^2));
    d = V/(4 * deltaX);
    q1 = e + d;
    q2 = 1/deltaT + 2*e +Tank.UL/(2*N_Layers);
    q3 = d - e ;
    q4 = 1/deltaT - 2*e - Tank.UL/(2*N_Layers);

 %%   
    % Generate Z1
    Z1 = zeros(N_Layers+1,N_Layers+1);

    Z1(1,1) = 1;

    for i=2:( N_Layers-1)
        Z1(i,i-1:i+1) = [-q1,q2,q3];
        Z1(i,N_Layers+1) = -Tank.UL/N_Layers;
    end

    Z1(N_Layers,N_Layers-2:N_Layers) = [1, -4, 3];

    Z1(N_Layers+1,N_Layers+1) = 1;
%%    
    % Generate Z2
    Z2 = zeros(N_Layers+1,N_Layers+1);

    Z2(1,1) = 1;

    for i=2:(N_Layers-1)
        Z2(i,i-1:i+1) = [q1,q4,-q3];
    end

    Z2(N_Layers+1,N_Layers+1) = 1;
  %%  
    % Generate Z3
    Z3 = zeros(N_Layers+1,1);
    Positions_index = HELayer(N_Layers, deltaX, HE.Positions); 
    for i = 1:HE.N
        pos = Positions_index(i);
        Z3(pos,:) = HE.n_eff*HE.Power/(m_i*Tank.Cv)*heatState(i); %K/s
    end
end