function [e,d,q1,q2,q3,q4,q5,q6,q7,q8,q9] = coeffs_(deltaX,V,deltaT,eps, Tank)
    e = eps*Tank.Dc/(2*deltaX*deltaX);
    d = V/(4*deltaX);
    q1 = e + d;
    q2 = 1/deltaT + 2*e + Tank.UL/2;
    q3 = d - e;
    q4 = 1/deltaT - 2*e - Tank.UL/2;
    q5 = 1/deltaT - e + Tank.UL_;
    q6 = 1/deltaT + e - Tank.UL_;
    q7 = q5 - 3*d;
    q8 = q6 + 3*d;
    q9 = 2*e + 4*d;
end