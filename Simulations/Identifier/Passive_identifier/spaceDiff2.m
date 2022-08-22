function M = spaceDiff2(Mat, deltaX)
%Mat is is data of temperature
%Rows representant space discretization and columns time discretization
[I,J] = size(Mat);
M = zeros(I, J);
for i = 1:I
    for j = 1:J
        if i == 1
            M(i,j) = Mat(i,j) - 2*Mat(i+1,j) + Mat(i+2,j);     
        elseif i == I
            M(i,j) = Mat(i,j) - 2*Mat(i-1,j) + Mat(i-2,j);
        else
            M(i,j) = Mat(i-1,j) - 2*Mat(i,j) + Mat(i+1,j);
        end
    end
end
M = M/(deltaX^2);
end