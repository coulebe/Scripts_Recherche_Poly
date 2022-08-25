% Initial condition
function u0=icfun_(x, xVector, initial)
    %We calculate the solution by each step of time, so initial condition
    %will change at each loop and it is stored on the initial tab

    [~,I] = min(abs(xVector - x));
    u0=initial(I);
end
