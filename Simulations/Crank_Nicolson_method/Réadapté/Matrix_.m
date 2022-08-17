function [Z1, Z2, Z3] = Matrix_(N_layers, V,deltaX,deltaT,eps, Tank, heatState, HE)
    %We suppose the heatState has the same size than HE.Positions_index
    m_i = Tank.Vol*Tank.Rho/N_layers;

    
    [e,d,q1,q2,q3,q4,q5,q6] = coeffs(deltaX,V,deltaT,eps, Tank);

    
    % Generate Z1
    Z1 = zeros(N_layers+1,N_layers+1);
    temp = [1, zeros(1,N_layers)];
 
    Z1(1,:) = temp;
    temp = [-4*d, 1/deltaT, zeros(1,N_layers-2), -Tank.UL_];
    Z1(2,:) = temp;
    for i=1:( N_layers-3)
        temp = [zeros(1,i),-q1,q2,q3,zeros(1,N_layers-3-i),-Tank.UL];
        Z1(i+2,:) = temp;
    end
%     temp = [zeros(1,N_layers-2), q6, q5, 0];
    temp = [zeros(1,N_layers-1), 1/deltaT, -Tank.UL_];
%     temp = [zeros(1,N_layers-1),-1, 1];
    Z1(N_layers,:) = temp;
    temp = [zeros(1,N_layers), 1];
    Z1(N_layers+1,:) = temp;
    
    % Generate Z2
    Z2 = zeros(N_layers+1,N_layers+1);
    temp = [1, zeros(1,N_layers)];
    Z2(1,:) = temp;
    temp = [0, q5, e, zeros(1,N_layers-2)];
    Z2(2,:) = temp;
    for i=1:(N_layers-3)
        temp = [zeros(1,i),q1,q4,-q3,zeros(1,N_layers-2-i)];
        Z2(i+2,:) = temp;
    end
    temp = [zeros(1,N_layers-2), q6, q5, 0];
    Z2(N_layers,:) = temp;
    temp = [zeros(1,N_layers), 1];
    Z2(N_layers+1,:) = temp;
    
    % Generate Z3
    Z3 = zeros(N_layers+1,1);
    Positions_index = HELayer(N_layers, deltaX, HE.Positions); 
    for i = 1:HE.N
        pos = Positions_index(i);
        Z3(pos,:) = HE.n_eff*HE.Power/(m_i*Tank.Cv)*heatState(i); %K/s
    end
end