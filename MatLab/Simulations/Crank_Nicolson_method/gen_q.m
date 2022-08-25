function [e,d,q1,q2,q3,q4,q5,q6] = coeffs(deltaX,V,deltaT,eps, )
    a = 1*0.14e-6; %m2/s
    k = 1*6.3588e-7; %l/s
    k_ = 1*1.2382e-6; %l/s
    e = eps*a/(deltaX*deltaX);
    d = V/(4*deltaX);
    q1 = e/2 + d;
    q2 = 1/deltaT + e + k/2;
    q3 = d - e/2;
    q4 = 1/deltaT - e - k/2;
    q5 = 1/deltaT - e - 4*d - k_;
    q6 = e + 4*d;
end