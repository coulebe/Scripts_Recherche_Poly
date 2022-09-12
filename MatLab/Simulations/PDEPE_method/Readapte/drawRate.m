% Water Draws function
function V = drawRate(t, Tank, Draw_Tab)
%     global Draw_Tab;
%     if(~isempty(Draw_Tab))
%         if(t > (Draw_Tab(1,1) * 60*60)) 
%             V = Draw_Tab(1,3)*1e-3*Tank.H/(Tank.Vol*60); %m/s
%         
%             if(t > (Draw_Tab(1,1) * 60*60 + Draw_Tab(1,2)*60))
%                 V = 0;
%                 Draw_Tab(1,:) = [];
%             end
%         else
%             V = 0;
%         end 
%     else
%         V = 0;
%     end
    if(~isempty(Draw_Tab))
        [I, ~] = size(Draw_Tab);
        for i = I:-1:1
            if(t > (Draw_Tab(i,1) * 60*60))
                break;
            end
        end
        if(t > (Draw_Tab(i,1) * 60*60)) && (t < (Draw_Tab(i,1) * 60*60 + Draw_Tab(i,2)*60))
            V = Draw_Tab(1,3)*1e-3*Tank.H/(Tank.Vol*60); %m/s   
        else
            V = 0;
        end
    else
        V = 0;
    end
end