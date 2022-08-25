function tab_pos = HELayer(N_layer, DeltaX, pos)
    %this function wil determine at which layer belong each heating
    %element, 
    %Pos is a (nx1) array have the positions of each heating element n <= Nlayer
    %N_layer is the number of layer and DeltaX the space step
    %We'll suppose that we can just have one heating element per layer for
    %the moment
    n = size(pos);
    n = n(1);
    tab_pos = -1*ones(n,1);
    count = 0;
    for i = 1:n
        for j = 1:N_layer
            if (DeltaX*(j-1) <= pos(i)) && (pos(i) < DeltaX*j)
                tab_pos(i) = j;
                count = count +1;
                break
            end
        end
    end
    if(count ~= n)
        error('Some HE positions are out of boundaries')
    end
    tab_pos = unique(tab_pos);
    
end