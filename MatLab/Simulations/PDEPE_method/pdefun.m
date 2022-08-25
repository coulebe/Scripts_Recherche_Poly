% Source element chauf, F fonction
function[c,f,s]=pdefun(x,t,T,dTdx)
    global n;
    global Tank;
    global HeatElem;
    global T_amb;
    global T_target;
    global HS;
%     global eps;
    
    deltaX = Tank.H/(n-1);
    m_i = Tank.Vol*Tank.Rho/n;
    Q = 0;%

%     theta = 0.46;
%     d = 2*sqrt(Tank.Vol/(Tank.H*pi));
%     g = 9.8;
%     beta = 207e-6;
    %c
    c=1;
    %f
    f=1*(Tank.Dc)*dTdx;
%     f = Tank.Dc*dTdx + (2/3) * ((theta*d)^2) * sqrt(- g * beta * dTdx^3);
    H = zeros(HeatElem.N,1);
    for i = HeatElem.N:-1:1
        x1 = HeatElem.Positions(i)-deltaX;
        x2 = HeatElem.Positions(i)+deltaX;
        if(x > x1 && x < x2)
            if(T < T_target)
                Q = HeatElem.Power * HeatElem.n_eff/(m_i*Tank.Cv);
                H(i) = 1;
                break;
            end
        end
    end
    H = [t;H];
    HS = [HS, H];
    %s
    V = drawRate(t);
    s= -dTdx*(V) - Tank.UL*T + Tank.UL*T_amb + Q;

    global iterator;
    iterator = iterator + 1;
    

end