import numpy as np
#%%
#%%Objects
# Tank = struct('Vol', {}, 'H', {}, 'Cv', {},'Rho', {}, 'Dc', {}, 'UL', {}, 'UL_', {});
#HeatElem = struct('n_eff', {}, 'Power', {}, 'Positions',{}, 'Thermos', {}, 'N', {});
class Tank:
    '''Class containning the Water tank property
    Attributs:
        Vol: Volume (m³)
        H: Height (m)
        Cv: Heat capacity coefficient ( J/(kg T) )
        Rho: Density (kg/m³)
        Dc: Thermal diffusity coefficient (m²/s)
        UL: Thermal losses coefficient (s⁻¹)'''
    def  __init__(self, Vol_, H_, UL_, Dc_ =0.14558e-6, Rho_ = 1e3, Cv_ =  4185.5):
        self.vol = Vol_
        self.H = H_
        self.Cv = Cv_
        self.Rho = Rho_
        self.Dc = Dc_
        self.UL = UL_

class HeatElem:
    '''Class containning the Heating Elements (HE) informations
    Attributs:
        n_eff: HE efficiency coefficient
        Power: Electrical power delivered by each HE(Watt)
        Positions: Position of each HE
        Thermos: Positions of Thermostat for eah HE
        N: Number of HE
    Assumptions:
        N = size(Position) = size(Thermos)
        Thermos and Positions elements are entered already sorted from the lowest to the highest HE. For 0 <= i< N, Position[i]'s and Thermos[i]'s values are near'''
    def __init__(self, n_eff_, Power_, Positions_, Thermos_, N_):
        self.n_eff = n_eff_
        self.Power = Power_
        self.Positions = Positions_
        self.Thermos = Thermos_
        self.N = N_
        



#%%Functions
#%%
def HELayer(N_layer, deltaX, pos_vec):
    '''this function wil determine at which layer belong each heating
    element
    We'll suppose that we can just have one heating element per layer for
    the moment
    Parametes: 
          N_layers : number of layers (space discretization)
          deltaX: Space step(m)
          Pos_vec (nx1) array having the positions of each heating element n <= Nlayer
    Output: tab_pos a vector with the index of each heating element'''
    
    N = np.size(pos_vec)
    count = 0
    tab_pos = -1*np.ones(N)
    for i in range(N):
        for j in range(N_layer):
            if j*deltaX <= pos_vec[i] < (j+1)*deltaX:
              tab_pos[i] = j
              count += 1
              break
          
    if(count != N):
        return np.unique(tab_pos)
    else:
        print("Some HE positions are out of boundaries")
        return None
    
    
#%%
def matrix(N_layer, V, deltaX, deltaT, eps, heatState, Tank, HE):
    '''
    Used to construct the matrices for the EWH's temperature dynamics
    Parameter: 
        N_layers: number of layers (space discretization)
        V: Water draw 'debit' (m/s)
        deltaX: Space step(m)
        deltaT: Time step(s)
        eps: ɑ(cf our equation)
        heatState: Describe the state(ON/OFF) of each heating element
        Tank: Object having the tank Caracteristics(Cf higher for more details)
        HE: Object having the Heating elements caracteristics(Cf higher for more details)
    Output:
        Z1, Z2, Z3 such as Z1 T(k+1) = Z2 T(k) + Z3
        '''
    ##Coefficients
    m_i = Tank.Vol*Tank.Rho/N_layer
    e = eps + Tank.Dc/(2 * (deltaX**2))
    d = V/(4 * deltaX)
    q1 = e+ d
    q2 = 1/deltaT + 2*e +Tank.UL/(2*N_layer)
    q3 = d -e 
    q4 = 1/deltaT - 2*e - Tank.UL/(2*N_layer)
    
    
    #G#enerate Z1
    Z1 = np.zeros((N_layer+1, N_layer+1))
    Z1[0,0] = 1
    
    for i in range(1, N_layer - 1):
        Z1[i, i-1:i+2] = [-q1, q2, q3]
        Z1[i, N_layer+1] = -Tank.UL/N_layer
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
    Positions_index = HELayer(N_layer, deltaX, HE.Positions) #Consider that HE position were already sorted
    for i in range(HE.N):
        pos = Positions_index[i]
        Z3[pos, :] = HE.n_eff * HE.Power/(m_i * Tank.Cv) * heatState[i]  #K/s
    return Z1, Z2, Z3

#%%
def PowerState(T_vec, HE, Target, deltaX, N_layers):
    '''
    Mean to indicate what heating elements need to be activated. The higher is the heating element's position
    the most prior it is. 
    Parameter: 
        T_vec: Vector cotainning temperature at each layer
        HE: Object having the Heating elements caracteristics(Cf higher for more details)
        Target: The reference that for the heating elemt activation
        deltaX: Space step(m)
        N_layers: number of layers (space discretization)
    output: 
        heatState: Vector describing the state(ON/OFF) of each heating element
        '''
    heatState = np.zeros(HE.N) #1= ON/ 0 = OFF
    Thermos_index = HELayer(N_layers, deltaX, HE.Thermos) #Consider that thermos positions were already sorted 
    for i in reversed(range(HE.N)):
        if T_vec[Thermos_index] < Target:
            heatState[i] = 1
            break
    return heatState

#%%

def CN_meth(Tank, HeatElem, Draw_Tab, deltaT, sim_time,N_layer, T_init, T_amb, T_in, T_target, eps  ):
    '''
    This function resolve the pdepe of the temperature in a water tank by using Crank-Nicolson discretization scheme
    Parameters:
        Tank: Object having the tank Caracteristics(Cf higher for more details)
        HeatElem: Object having the Heating elements caracteristics(Cf higher for more details)
        Draw_tab: Tab containing the draw planned during the simulation. Each 
                line contain the hour when  the draw start(h), its duration(min), and
                its debit(l/min). The tab is ascending given the hour of start.
        deltaT: Time step(s)
        sim_time: Duration of the simulation(h)
        N_layers: number of layers (space discretization)
        T_init: initial temperature in the tank(°C)
        T_amb: Ambient temperaure(°C)
        T_in: inlet fluid temperature(°C)
        T_target: Target temperature for the heating element control(°C)
        eps coefficient that replace the thermal expansion coefficient
    Outputs:
        tsol: solution time vector(1xM)(sec)
        xVector: Space vector (1xN_layers)(m)
        sol: Solution matrice(MxN)
        
        
        '''
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    