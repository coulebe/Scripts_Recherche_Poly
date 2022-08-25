function [Z1, Z2] = Matrix(N_layers, V,deltaX,deltaT,eps, Dc, UL)



    
    [~,~,q1,q2,q3,q4,] = coeffs(deltaX,V,deltaT,eps, Dc, UL);

 %%   
    % Generate Z1
    Z1 = zeros(N_layers+1,N_layers+1);

 
    Z1(1,1) = 1;

    for i=0:( N_layers-3)
        temp = [zeros(1,i),-q1,q2,q3,zeros(1,N_layers-3-i),-UL];
        Z1(i+2,:) = temp;
    end

%     temp = [zeros(1,N_layers-1), 1/deltaT, -Tank.UL_];

    Z1(N_layers,:) = [zeros(1,N_layers - 3), 1, -4, 3, 0];
%     temp = [zeros(1,N_layers), 1];
    Z1(N_layers+1,N_layers+1) = 1;
%%    
    % Generate Z2
    Z2 = zeros(N_layers+1,N_layers+1);
%     temp = [1, zeros(1,N_layers)];
    Z2(1,1) = 1;
%     temp = [0, q5, e, zeros(1,N_layers-2)];
%     Z2(2,:) = temp;
    for i=0:(N_layers-3)
        temp = [zeros(1,i),q1,q4,-q3,zeros(1,N_layers-2-i)];
        Z2(i+2,:) = temp;
    end
%     temp = [zeros(1,N_layers-2), q6, q5, 0];
%     Z2(N_layers,N_layers) = temp;
%     temp = [zeros(1,N_layers), 1];
    Z2(N_layers+1,N_layers+1) = 1;
  
end