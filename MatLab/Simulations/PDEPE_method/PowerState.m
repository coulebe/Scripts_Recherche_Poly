function Pow = PowerState(Tm,HE,T_Target)
    %PowerState function is to indicate what heat element need to be
    %activated. The priority order is the higher is the heat element, the
    %most it is prior
    %Tm is the current temperature of each layer
    %positions is the position index of the heating elements
    %
    Pow = zeros(HE.N,1); %1  on, 0 = off
    for i = HE.N:-1:1 %If the highest layer is too cold, we don't need to check if we have to heat the others
        if(Tm(HE.Positions_index(i)) < (T_Target))
            Pow(i) = 1;
            break
        else
           Pow(i) = 0; 
        end
    end
end