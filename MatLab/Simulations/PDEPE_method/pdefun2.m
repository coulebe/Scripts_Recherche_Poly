% Source element chauf, F fonction
function[c,f,s]=pdefun2(x,t,T,dTdx)
    global n;
    global Tank;
    global HeatElem;
    global T_amb;
    global heatState;
    global eps;
    
    deltaX = Tank.H/(n-1);
    m_i = Tank.Vol*Tank.Rho/n;
    Q = 0;

%     theta = 0.46;
%     d = 2*sqrt(Tank.Vol/(Tank.H*pi));
%     g = 9.8;
%     beta = 207e-6;
    %c
    c=1;
    %f
    f=eps*(Tank.Dc)*dTdx;
%     f = Tank.Dc*dTdx + (2/3) * ((theta*d)^2) * sqrt(- g * beta * dTdx^3);
    H = zeros(HeatElem.N,1);
    for i = HeatElem.N:-1:1
        x1 = HeatElem.Positions(i)-deltaX;
        x2 = HeatElem.Positions(i)+deltaX;
        if(x > x1 && x < x2)
            Q = HeatElem.Power * HeatElem.n_eff*heatState(i)/(m_i*Tank.Cv);
            break;
            
        end
    end
    

    %s
    V = drawRate(t);
    s= -dTdx*(V) - Tank.UL*T + Tank.UL*T_amb + Q;

    global iterator;
    iterator = iterator + 1;
    

end