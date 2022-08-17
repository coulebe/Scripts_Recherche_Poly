function [heatState1,heatState2] = getPowerState(Tm,pos1,pos2,T_Target,offset)
    %getPowerState Summary of this function goes here
    % Detailed explanation goes here
    on = 1;
    off = 0;
    if (Tm(pos1+offset) > (T_Target-0))
        heatState1 = off;
    else
        heatState1 = on;
    end
    if (Tm(pos2+offset) > (T_Target-0))
        heatState2 = off;
    else
        heatState2 = on;
    end
    if (heatState1 == heatState2)
        if (heatState2 == on)
            heatState1 = off;
        end
    end
end