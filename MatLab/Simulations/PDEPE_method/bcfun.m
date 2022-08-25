function [pl,ql,pr,qr] = bcfun(xl,ul,xr,ur,t)
    global currentTemp;
    global T_in;  global n;
    V = drawRate(t);
    if (V == 0)
        pl = ul - currentTemp(2);
        ql = 0;
        pr = ur - currentTemp(n-1);
        qr = 0;
    else
        pl = ul - T_in;
        ql = 0;
        pr = ur - currentTemp(n-1);
        qr = 0;
    end
end