#%%Imports
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
# from matplotlib.ticker import LinearLocator
import CN_Functions as CN


#%%Initialisation
#%%
Tank_ = CN.Tank(112e-3, 1.19, 6.3588e-6)
HE = CN.HeatElem(.95, 6e3, np.array([[0.2975], [.7735]]), \
              np.array([[0.2975], [.7735]]), 2)

DrawTab = np.array([
                    [2.5, 40, 3], \
                    [5, 15, 6]\
                        ])
#%%
deltaT = 1 #s
sim_time = 5 #h
N = 100
#%%
T_init = 25
T_amb = 25
T_in = 25
T_target = 60
eps = 150
#%%Simulation
#%%
tsol, xVector, sol = CN.CN_meth(Tank_, HE, DrawTab, deltaT, sim_time, N, T_init, T_amb, T_in, T_target, eps)


#%%Plotting
#%%

fig, ax = plt.subplots(subplot_kw={"projection": "3d"})
tsol_mesh, xVector_mesh = np.meshgrid(tsol/3600, xVector)

surf = ax.plot_surface(tsol_mesh, xVector_mesh, sol[0:N, :], cmap=cm.coolwarm, linewidth=0, antialiased=False)
fig.colorbar(surf, shrink=10, aspect=5)

plt.show()