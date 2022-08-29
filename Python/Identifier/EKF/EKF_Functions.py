#%%Import
import numpy as np


#%%
#%%Matrix
def matrix(N_layer, V, deltaX, deltaT, eps, Dc, UL):
    '''
    Used to construct the matrices for the EWH's temperature dynamics, for our identifier
    Parameter: 
        N_layers: number of layers (space discretization)
        V: Water draw 'debit' (m/s)
        deltaX: Space step(m)
        deltaT: Time step(s)
        eps: ɑ(cf our equation)
        Dc: Thermal diffusity coefficient (m²/s)
        UL: Thermal losses coefficient (s⁻¹) 
    Output:
        Z1, Z2, Z3 such as Z1 T(k+1) = Z2 T(k) + Z3
        '''
    ##Coefficients
    e = eps * Dc/(2 * (deltaX**2))
    d = V/(4 * deltaX)
    q1 = e+ d
    q2 = 1/deltaT + 2*e +UL/(2*N_layer)
    q3 = d -e 
    q4 = 1/deltaT - 2*e - UL/(2*N_layer)
    
    
    #G#enerate Z1
    Z1 = np.zeros((N_layer+1, N_layer+1))
    Z1[0,0] = 1
    
    for i in range(1, N_layer - 1):
        Z1[i, i-1:i+2] = [-q1, q2, q3]
        Z1[i, N_layer] = -UL/N_layer
    Z1[N_layer-1, N_layer - 3:N_layer] = [1, -4, 3]
     
    Z1[N_layer,N_layer] = 1
    
    ##Generate Z2
    Z2 = np.zeros((N_layer+1, N_layer+1))
    Z2[0,0] = 1
    
    for i in range(1, N_layer - 1):
        Z2[i, i-1:i+2] = [q1, q4, -q3]
   
    Z2[N_layer,N_layer] = 1
    
    
    return Z1, Z2

#%%
def f(X_k, N_layers, V, deltaX, deltaT, Q_vec, T_in):
    '''
    Our System dynamics based on its PDE
    Parameters:
        X_k: (N+4)Vector. Augmented state at step k. X_k = [T_k; Theta_k] where theta_k = [alpha_k, Dc_k, UL_k] 
            and T_k = [T1_k, T2_k, ..., TN_k, T_amb ]. 
        N_layers: number of layers (space discretization)
        V: Water draw 'debit' (m/s)
        deltaX: Space step(m)
        deltaT: Time step(s)
        Q_vec: (N)Vector containing the values of the energy (°C/s)
        T_in: inlet fluid temperature(°C)
    Output:
        X_k_plus_1: (N+4)Vector. Augmented state at step k+1.'''
    Z1, Z2 = matrix( N_layers, V, deltaX, deltaT,X_k[N_layers+1], X_k[N_layers+2], X_k[N_layers+3] )
    #Coumpound matrices for augmented state
    A = np.concatenate( (\
                    np.concatenate( ( Z1, np.zeros((N_layers+1, 3)) ), axis = 1), \
                    np.concatenate( ( np.zeros((3, N_layers+1)),  np.eye(3)), axis = 1)    \
                     ), axis = 0)
    
    B = np.concatenate( (\
                    np.concatenate( ( Z2, np.zeros((N_layers+1, 3)) ), axis = 1), \
                    np.concatenate( ( np.zeros((3, N_layers+1)),  np.eye(3)), axis = 1)    \
                     ), axis = 0) 
    C = np.concatenate( (Q_vec, np.zeros(4)), axis = 0 )
    
    #
    X_k_plus_1 = np.linalg.inv(A)@(B@X_k + C)
    #Conditions for first element
    if V > 0:
        X_k_plus_1[0] = T_in
    else:
        X_k_plus_1[0] = (4*X_k_plus_1[1] - X_k_plus_1[2])/3

    return X_k_plus_1
#%%
def F1_mat(X_k, N_layers, V,deltaX,deltaT):
    '''
    
    Parameters:
        X_k: (N+4)Vector. Augmented state at step k. X_k = [T_k; Theta_k] where theta_k = [alpha_k, Dc_k, UL_k] 
            and T_k = [T1_k, T2_k, ..., TN_k, T_amb ]. 
        N_layers: number of layers (space discretization)
        V: Water draw 'debit' (m/s)
        deltaX: Space step(m)
        deltaT: Time step(s)

    Output:
        F1: Jacobian of f1 function applied on X_k.'''
    
    Z1 = matrix( N_layers, V, deltaX, deltaT,X_k[N_layers+1], X_k[N_layers+2], X_k[N_layers+3] )[0]
    a = 0.5/(deltaX**2)
    Jac_f1_theta_k = np.zeros((N_layers+1, 3))
    
    for i in range(1,N_layers-1):
        Jac_f1_theta_k[i,0] = X_k[N_layers+2]*a*(-X_k[i-1]+ 2*X_k[i] - X_k[i+1]);
        Jac_f1_theta_k[i,1] = X_k[N_layers+1]*a*(-X_k[i-1]+ 2*X_k[i] - X_k[i+1]);
        Jac_f1_theta_k[i,2] = (X_k[i]/2 - X_k[N_layers+1])/10;
    F1 = np.concatenate( (\
                    np.concatenate( ( Z1, Jac_f1_theta_k ), axis = 1), \
                    np.concatenate( ( np.zeros((3, N_layers+1)),  np.eye(3)), axis = 1)    \
                     ), axis = 0)
    return F1

#%%
def F2_mat(X_k, N_layers, V,deltaX,deltaT):
    '''
    Parameters:
        X_k: (N+4)Vector. Augmented state at step k. X_k = [T_k; Theta_k] where theta_k = [alpha_k, Dc_k, UL_k] 
            and T_k = [T1_k, T2_k, ..., TN_k, T_amb ]. 
        N_layers: number of layers (space discretization)
        V: Water draw 'debit' (m/s)
        deltaX: Space step(m)
        deltaT: Time step(s)

    Output:
        F2: Jacobian of f2 function applied on X_k.'''
    
    Z2 = matrix( N_layers, V, deltaX, deltaT,X_k[N_layers+1], X_k[N_layers+2], X_k[N_layers+3] )[1]
    a = 0.5/(deltaX**2)
    Jac_f2_theta_k = np.zeros((N_layers+1, 3))
    
    for i in range(1,N_layers-1):
        Jac_f2_theta_k[i,0] = X_k[N_layers+2]*a*(-X_k[i-1]+ 2*X_k[i] - X_k[i+1]);
        Jac_f2_theta_k[i,1] = X_k[N_layers+1]*a*(-X_k[i-1]+ 2*X_k[i] - X_k[i+1]);
        Jac_f2_theta_k[i,2] = - X_k[N_layers+1]/10;
    F2 = np.concatenate( (\
                    np.concatenate( ( Z2, Jac_f2_theta_k ), axis = 1), \
                    np.concatenate( ( np.zeros((3, N_layers+1)),  np.eye(3)), axis = 1)    \
                     ), axis = 0)
    return F2

#%%
def F_mat(X_k, N_layers, V,deltaX,deltaT,Q_vec):
    '''
    Parameters:
        X_k: (N+4)Vector. Augmented state at step k. X_k = [T_k; Theta_k] where theta_k = [alpha_k, Dc_k, UL_k] 
            and T_k = [T1_k, T2_k, ..., TN_k, T_amb ]. 
        N_layers: number of layers (space discretization)
        V: Water draw 'debit' (m/s)
        deltaX: Space step(m)
        deltaT: Time step(s)
        Q_vec: (N)Vector containing the values of the energy (°C/s)

    Output:
        F: Jacobian of f function applied on X_k.'''
    Z2 = matrix( N_layers, V, deltaX, deltaT,X_k[N_layers+1], X_k[N_layers+2], X_k[N_layers+3] )[1]
    #calculate f2(X_k,Q_vec)
    B = np.concatenate( (\
                    np.concatenate( ( Z2, np.zeros((N_layers+1, 3)) ), axis = 1), \
                    np.concatenate( ( np.zeros((3, N_layers+1)),  np.eye(3)), axis = 1)    \
                     ), axis = 0) 
    C = np.concatenate( (Q_vec, np.zeros(4)), axis = 0 )    
    
    f2 = B@X_k + C
    
    #Now we can calculate F1 applied on f2(X_k,Q_vec)
    F1_of_f2 = F1_mat(f2,N_layers, V, deltaX, deltaT)
    #Now we can calculate F2 applied on X_k
    F2_of_X_k = F2_mat(X_k,N_layers, V, deltaX, deltaT)
    
    F = np.linalg.inv(F1_of_f2)@F2_of_X_k
    
    return F
    
        
        