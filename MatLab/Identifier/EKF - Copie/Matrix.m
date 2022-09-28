function [Z1, Z2] = Matrix(N_layers, V,deltaX,deltaT,eps, Tank)
    %We suppose the heatState has the same size than HE.Positions_index
    
    e = eps * Tank.Dc/(2 * (deltaX^2));
    d = V/(4 * deltaX);
    q1 = e + d;
    q2 = 1/deltaT + 2*e +Tank.UL/2;
    q3 = d - e ;
    q4 = 1/deltaT - 2*e - Tank.UL/2;

 %%   
    % Generate Z1
    Z1 = zeros(N_layers+1,N_layers+1);

    Z1(1,1) = 1;

    for i=2:( N_layers-1)
        Z1(i,i-1:i+1) = [-q1,q2,q3];
        Z1(i,N_layers+1) = -Tank.UL;
    end

    Z1(N_layers,N_layers-2:N_layers) = [1, -4, 3];

    Z1(N_layers+1,N_layers+1) = 1;
%%    
    % Generate Z2
    Z2 = zeros(N_layers+1,N_layers+1);

    Z2(1,1) = 1;

    for i=2:(N_layers-1)
        Z2(i,i-1:i+1) = [q1,q4,-q3];
    end

    Z2(N_layers+1,N_layers+1) = 1;
  
end