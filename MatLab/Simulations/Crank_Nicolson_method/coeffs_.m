function [e,d,q1,q2,q3,q4] = coeffs_(deltaX,V,deltaT,eps, Tank)
    e = eps*Tank.Dc/(2*deltaX*deltaX);
    d = V/(4*deltaX);
    q1 = e + d;
    q2 = 1/deltaT + 2*e + Tank.UL/2;
    q3 = d - e;
    q4 = 1/deltaT - 2*e - Tank.UL/2;
    
end