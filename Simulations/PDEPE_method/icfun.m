% Initial condition
function u0=icfun(x)
    %We calculate the solution by each step of time, so initial condition
    %will change at each loop and it is stored on the initial tab
    global xVector;
    global initial;
    [M,I] = min(abs(xVector - x));
    u0=initial(I);
end
