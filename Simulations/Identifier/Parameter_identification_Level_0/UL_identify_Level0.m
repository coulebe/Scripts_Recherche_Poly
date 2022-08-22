clear;
load('Data_UL_Identification');
%%
d2Tdx2 = spaceDiff2WoBoundary(data, deltaX);
dTdt = diff(data(2:9,:), 1, 2)/deltaT;
[I,J] = size(dTdt);
mat = zeros(I,J);
for i = 1:I
    for j = 1:J
        mat(i,j) = (-dTdt(i,j) + d2Tdx2(i,j)*eps*Tank.Dc)/(data(i+1,j) - T_amb);
    end
end

UL_mat = mean(mat/10, 2);
