#%%Imports
import numpy as np
import numpy.matlib as mb
import matplotlib.pyplot as plt
# from matplotlib.ticker import LinearLocator
import CN_Functions as CN
import pandas as pd 
import os
import tkinter.filedialog




#%%Initialisation
#%%
Tank_ = CN.Tank(112e-3, 1.19, 6.3588e-6)
HE = CN.HeatElem(.95, 6e3, np.array([[0.2975], [.7735]]), \
              np.array([[0.2975], [.7735]]), 2)
#%%
#Find the Draw file
fTyp = [("fichier de données (*.csv)", "*.csv")]
iDir = os.path.abspath(os.path.dirname(__file__))
file_name = tkinter.filedialog.askopenfilename(filetypes=fTyp, initialdir=iDir)

DrawTab_DF = pd.read_csv(file_name)
#%%
deltaT = 1 #s
sim_time = 5 #h
N = 10
#%%
T_init = 70
T_amb = 25
T_in = 25
T_target = 70
eps = 150
#%%Simulation
#%%
tsol, xVector, sol, Q_mat, vVec = CN.CN_meth(Tank_, HE, DrawTab_DF, deltaT, sim_time, N, T_init, T_amb, T_in, T_target, eps, He_Activ = False )


#%%Plotting
#%%

fig, ax = plt.subplots()
tsol_mesh, xVector_mesh = np.meshgrid(tsol/3600, xVector)
im = ax.contourf(tsol_mesh,xVector_mesh,  sol)
ax.set_ylabel('x(m)')
ax.set_xlabel('t(h)')
fig.colorbar(mappable=im)

plt.show()
#%%Save
vVecmesh = mb.repmat(vVec, N, 1)

dico = {'x' : xVector_mesh,  't' : tsol_mesh*3600, 'u' : sol, 'Q' : Q_mat, \
            'V': vVecmesh}
np.save('EWH_sim_NullPower_VarFlow_CN_data', dico)
#%%
