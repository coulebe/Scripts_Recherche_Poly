function Pow = PowerState(Tm,HE,T_Target, deltaX, N)
    %PowerState function is to indicate what heat element need to be
    %activated. The priority order is the higher is the heat element, the
    %most it is prior
    %Tm is the current temperature of each layer
    %
    Pow = zeros(HE.N,1); %1  on, 0 = off
    Thermos_index = HELayer(N, deltaX, HE.Thermos); 
    for i = HE.N:-1:1 %If the highest layer is too cold, we don't need to check if we have to heat the others
        if(Tm(Thermos_index(i)) < T_Target)
            Pow(i) = 1;
            break
        end
    end
end