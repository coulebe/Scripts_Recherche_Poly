function U = U_mat(Tank, HE, N, heatState, deltaX, deltaT)
    %size(Nx1 matrice)
    Positions_index = HELayer(N, deltaX, HE.Positions);
    m_i = Tank.Vol*Tank.Rho/N;
    %%
    U = zeros(N, 1);
    for i = 1:HE.N
        pos = Positions_index(i);
        U(pos,:) = HE.n_eff*HE.Power/(m_i*Tank.Cv)*heatState(i);
    end
    U = deltaT * U;
end