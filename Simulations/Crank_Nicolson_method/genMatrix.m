function [Z1, Z2, Z3] = genMatrix(n, e, d, deltaT, q1, q2, q3, q4, q5, q6, heatState1, heatState2,pos1, pos2)
    k = 1*6.3588e-7; %l/s
    k_ = 1*1.2382e-6; %l/s
    m_tank = 112;
    m_i = m_tank/n;
    n_eff = 0.95; % Efficiency = 95%
    P_ele = 1*6000; %watt
    c = 4185.5; %Water specific heat capacity
    % Generate Z1
    Z1 = zeros(n,n);
    temp = [1, zeros(1,n-1)];
    Z1(1,:) = temp;
    temp = [-4*d, 1/deltaT, zeros(1,n-3), -k_];
    Z1(2,:) = temp;
    for i=1:(n-4)
        temp = [zeros(1,i),-q1,q2,q3,zeros(1,n-4-i),-k];
        Z1(i+2,:) = temp;
    end
    temp = [zeros(1,n-2), 1/deltaT, -k_];
    Z1(n-1,:) = temp;
    temp = [zeros(1,n-1), 1];
    Z1(n,:) = temp;
    % Generate Z2
    Z2 = zeros(n,n);
    temp = [1, zeros(1,n-1)];
    Z2(1,:) = temp;
    temp = [0, q5, e, zeros(1,n-3)];
    Z2(2,:) = temp;
    for i=1:(n-4)
        temp = [zeros(1,i),q1,q4,-q3,zeros(1,n-3-i)];
        Z2(i+2,:) = temp;
    end
    temp = [zeros(1,n-3), q6, q5, 0];
    Z2(n-1,:) = temp;
    temp = [zeros(1,n-1), 1];
    Z2(n,:) = temp;
    % Generate Z3
    Z3 = zeros(n,1);
    Z3(pos1,:) = n_eff*P_ele/(m_i*c)*heatState1;
    Z3(pos2,:) = n_eff*P_ele/(m_i*c)*heatState2;
end