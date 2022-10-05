function Pow = PowerState_(Tm,HE,T_Target,x )
    %PowerState function is to indicate what heat element need to be
    %activated. The priority order is the higher is the heat element, the
    %most it is prior
    %Tm is the current temperature of each layer
    %positions is the position index of the heating elements
    %
%     a = HE.N;
    Pow = zeros(HE.N,1); %1  on, 0 = off
    for i = HE.N:-1:1 %If the highest layer is too cold, we don't need to check if we have to heat the others
            [~,I] = min(abs(x - HE.Thermos(i)));
        if(Tm(I) < T_Target)
            Pow(i) = 1;
            break
        else
           Pow(i) = 0; 
        end
    end
end