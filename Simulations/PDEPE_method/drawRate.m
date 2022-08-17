% Water Draws function
function V = drawRate(t)
    global Draw_Tab;
    global Tank;
    if(~isempty(Draw_Tab))
        if(t > (Draw_Tab(1,1) * 60*60)) 
            V = Draw_Tab(1,3)*1e-3*Tank.H/(Tank.Vol*60); %m/s
        
            if(t > (Draw_Tab(1,1) * 60*60 + Draw_Tab(1,2)*60))
                V = 0;
                Draw_Tab(1,:) = [];
            end
        else
            V = 0;
        end 
    else
        V = 0;
    end
end