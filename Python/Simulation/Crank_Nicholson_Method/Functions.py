import numpy as np

#%%
def HELayer(N_layer, deltaX, pos_vec):
    # this function wil determine at which layer belong each heating
    # element, 
    # Pos is a (nx1) array have the positions of each heating element n <= Nlayer
    # N_layer is the number of layer and DeltaX the space step
    # We'll suppose that we can just have one heating element per layer for
    # the moment
    N = np.size(pos_vec)
    count = 0
    tab_pos = -1*np.ones(n)
    for i in range(N):
        for j in range(N_layer):
            if j*deltaX <= pos_vec[i] < (j+1)*deltaX:
              tab_pos[i] = j
              count += 1
              break
          
    if(count != n):
        return np.unique(tab_pos)
    else:
        print("Some HE positions are out of boundaries")
        return None
    
    
#%%
def matrix(N_layer, V, deltaX, deltaT, eps, heatState, Tank, HE):
    
    ##Coefficients
    m_i = Tank.Vol*Tank.Rho/N_layer
    e = eps + Tank.Dc/(2 * (deltaX**2))
    d = V/(4 * deltaX)
    q1 = e+ d
    q2 = 1/deltaT + 2*e +Tank.UL/2
    q3 = d -e 
    q4 = 1/deltaT - 2*e - Tank.UL/2
    
    
    #G#enerate Z1
    Z1 = np.zeros((N_layer+1, N_layer+1))
    Z1[0,0] = 1
    
    for i in range(1, N_layer - 1):
        Z1[i, i-1:i+2] = [-q1, q2, q3]
        Z1[i, N_layer+1] = -Tank.UL
    Z1[N_layer, N_layer - 2:N_layer+1] = [1, -4, 3]
     
    Z1[N_layer+1,N_layer+1] = 1
    
    ##Generate Z2
    Z2 = np.zeros((N_layer+1, N_layer+1))
    Z1[0,0] = 1
    
    for i in range(1, N_layer - 1):
        Z1[i, i-1:i+2] = [q1, q4, -q3]
   
    Z1[N_layer+1,N_layer+1] = 1
    
    ##Generate Z3
    Z3 = np.zeros((N_layer+1, 1))
    Positions index = HELayer(N_layer, deltaX, HE.Positions)
    for i in range():
        pos = Positions_index[i]
        Z3[pos, :] = HE.n_eff * HE.Power/(m_i * Tank.Cv) * heatState[i] 
    return Z1, Z2, Z3