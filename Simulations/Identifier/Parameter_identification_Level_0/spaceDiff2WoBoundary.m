function M = spaceDiff2WoBoundary(Mat, deltaX)
%Mat is is data of temperature
%Rows representant space discretization and columns time discretization
[I,J] = size(Mat);
M = zeros(I-2, J);
for i = 1:I-2
    for j = 1:J
        M(i,j) = Mat(i,j) - 2*Mat(i+1,j) + Mat(i+2,j);
    end
end
M = M/(deltaX^2);
end