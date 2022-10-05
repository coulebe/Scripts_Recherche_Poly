function [e,d,q1,q2,q3,q4,q5,q6] = coeffs(deltaX,V,deltaT,eps, Tank, N_Layers)
    e = eps*Tank.Dc/(deltaX*deltaX);
    d = V/(4*deltaX);
    q1 = e/2 + d;
    q2 = 1/deltaT + e + Tank.UL/(2*N_Layers);
    q3 = d - e/2;
    q4 = 1/deltaT - e - Tank.UL/(2*N_Layers);
    q5 = 1/deltaT - e - 4*d - Tank.UL_/N_Layers;
    q6 = e + 4*d;
end