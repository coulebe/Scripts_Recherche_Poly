import numpy as np
import matplotlib.pyplot as plt
import sys
import scipy.io
import EKF_Functions as EKF_F
#%%
sys.path.append("D:/Alfred/Cours/Projet_Recherche/Poly/Scripts_Recherche_Poly/Python/Simulation/Crank_Nicholson_Method")
from CN_Functions import Tank, HeatElem
#%%Extracton des données 
#%%
mat = scipy.io.loadmat('Donnee_sim_Null_Flow.mat')
#%%
data = mat.get('data')
deltaT = mat.get('deltaT')[0,0]
deltaX = mat.get('deltaX')[0,0]
eps = mat.get('eps')[0,0]
# HeatElem = mat.get('HeatElem')[0,0]
Xvector = mat.get('Layers')[0,:]
N_layers = mat.get('N')[0,0]
Q_mat = mat.get('Q_mat')
T_amb = mat.get('T_amb')[0,0]
T_in = mat.get('T_in')[0,0]
# Tank = mat.get('Tank')[0,0]
tsol = mat.get('Time')[0,:]
V_vec = mat.get('V_vec')[:,0]

Tank_ = Tank(112e-3, 1.19, 6.3588e-6)
HE = HeatElem(.95, 6e3, np.array([[0.2975], [.7735]]), \
              np.array([[0.2975], [.7735]]), 2)
    
    
#%%Compound matrices
# = np.concatenate( (\
#                 np.concatenate( ( ,  ), axis = 1), \
#                 np.concatenate( ( ,  ), axis = 1)    \
#                  ), axis = 0)
#%%

T = np.shape(data)[1]
#Rajout de T_amb pour la matrice data
data = np.concatenate((data, T_amb*np.ones((1, T)) ), axis = 0)

#Variance matrices
W = np.concatenate( (\
                np.concatenate( (np.eye(N_layers+1) , np.zeros((N_layers+1,3)) ), axis = 1), \
                np.concatenate( (np.zeros((3,N_layers+1)) ,  np.diag([1, 1e-7, 1e-7])), axis = 1)    \
                 ), axis = 0)  #Process noise covariance matrice
R = np.eye(N_layers+1)
#
H = np.concatenate( (np.eye(N_layers+1) , np.zeros((N_layers+1,3)) ), axis = 1)


#%%EKF
#%%Initialisation
X_hat = np.zeros((N_layers+4, T))
theta_0 = np.array([eps*0.9, Tank_.Dc*0.9, Tank_.UL*0.9])

X_hat[:,0] = np.concatenate( (data[:,0] , theta_0 ), axis = 0)

P_k = np.concatenate( (\
                np.concatenate( ( np.eye(N_layers+1), np.zeros((N_layers+1,3)) ), axis = 1), \
                np.concatenate( ( np.zeros((3,N_layers+1)),  1e4*np.eye(3)), axis = 1)    \
                  ), axis = 0)
    

#%%Loop
for i in range(1,T):
    #Observation
    
    
    #Prediction
    X_ = EKF_F.f(X_hat[:,i-1], N_layers, V_vec[i], deltaX, deltaT, Q_mat[:,i], T_in)
    
    #Error Covariance propagation
    F = EKF_F.F_mat(X_hat[:,i-1], N_layers, V_vec[i], deltaX, deltaT, Q_mat[:,i])
    P_k_ = F @ P_k @ np.transpose(F) + W
    
    #Kalman Filter Gain
    K_k = P_k_  @ np.transpose(H) @ np.linalg.inv( H @ P_k_@ np.transpose(H) + R)
    
    #Update Step
    Y_ = data[:,i] - H @ X_
    
    X_hat[:,i] = X_ + K_k @ Y_
    P_k = (np.eye(N_layers+4) - K_k @ H) @ P_k_
    
Theta = X_hat[-3:,:]


#%%
#Result Plot  
plt.figure()
plt.plot(tsol/3600, Theta[0,:], label = '\epsilon _ hat')
plt.grid()
plt.hlines(eps, tsol[0], tsol[-1], label = '\epsilon')
plt.legend(loc = 'best')


plt.figure()
plt.plot(tsol/3600, Theta[1,:], label = 'Dc_ hat')
plt.grid()
plt.hlines(Tank_.Dc, tsol[0], tsol[-1], label = 'Dc')
plt.legend(loc = 'best')

    
plt.figure()
plt.plot(tsol/3600, Theta[2,:], label = 'UL_ hat')
plt.grid()
plt.hlines(Tank_.Dc, tsol[0], tsol[-1], label = 'UL')
plt.legend(loc = 'best')   
    
    
    
    
plt.show()
    
    
    
    
    
    
    
    
    

