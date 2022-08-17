clear;
load('Data_UL_Identification');
%%
d2Tdx2 = spaceDiff2WoBoundary(data, deltaX);
dTdt = diff(data(2:9,:), 1, 2)/deltaT;
[I,J] = size(dTdt);
mat = zeros(I,J);
for i = 1:I
    for j = 1:J
        mat(i,j) = (-dTdt(i,j) +(data(i+1,j) - T_amb)*Tank.UL)/(d2Tdx2(i,j)*eps);
    end
end

