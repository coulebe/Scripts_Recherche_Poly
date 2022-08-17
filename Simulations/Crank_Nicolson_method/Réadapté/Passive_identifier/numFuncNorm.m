function norme = numFuncNorm(data, x)
    %Calculate the norm of a represented by its data and x, the vector
    %represanting the space where the function is defined(Don't work if x
    %is R). data and x must have the same length
    data_sq = data.^2;
    norme = sqrt( trapz(x, data_sq));
end