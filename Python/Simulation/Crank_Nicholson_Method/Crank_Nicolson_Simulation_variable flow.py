#%%Imports
import numpy as np
import numpy.matlib as mb
import matplotlib.pyplot as plt
# from matplotlib.ticker import LinearLocator
import CN_Functions as CN
import pandas as pd 
import os
import tkinter as tk
from tkinter import filedialog
# from matplotlib import cm




#%%Initialisation
#%%
Tank_ = CN.Tank(112e-3, 1.19, 6.3588e-6)
HE = CN.HeatElem(.95, 6e3, np.array([[0.2975], [.7735]]), \
              np.array([[0.2975], [.7735]]), 2)
#%%
#Find the Draw file
try:
    application_window = tk.Tk() 
    fTyp = [("fichier de données (*.csv)", "*.csv")]
    file_name = filedialog.askopenfilename(parent=application_window,
                                        initialdir=os.getcwd(),
                                        title="Please select your csv file containing the Draw tab:",
                                        filetypes=fTyp)


    application_window.destroy()
    
    DrawTab = pd.read_csv(file_name)
except:
    DrawTab = pd.DataFrame(columns = ['Start(h)', 'Duration(min)', 'Draw Rate(L/min)'])# If we want a nullFlow sim
#%%
deltaT = 1 #s
sim_time = DrawTab['Start(h)'][max(DrawTab.index)] +  DrawTab['Duration(min)'][max(DrawTab.index)]/60 + 0.5 #We'll end the simulation 30 min after the last draw
N = 10
#%%
T_init = 70
T_amb = 25
T_in = 25
T_target = 70
eps = 150
#%%Simulation
#%%
tsol, xVector, sol, Q_mat, vVec = CN.CN_meth(Tank_, HE, DrawTab, deltaT, sim_time, N, T_init, T_amb, T_in, T_target, eps, He_Activ = False )


#%%Plotting
#%%

fig, ax = plt.subplots()
tsol_mesh, xVector_mesh = np.meshgrid(tsol/3600, xVector)
im = ax.contourf(tsol_mesh,xVector_mesh,  sol)

# im = ax.plot_surface(tsol_mesh,xVector_mesh,  sol, cmap=cm.coolwarm,\
#                        linewidth=0, antialiased=False)
ax.set_ylabel('x(m)')
ax.set_xlabel('t(h)')
fig.colorbar(mappable=im)
# fig.colorbar(im, shrink=0.5, aspect=5)    
plt.show()
#%%Save
vVecmesh = mb.repmat(vVec, N, 1)

dico = {'x' : xVector_mesh,  't' : tsol_mesh*3600, 'u' : sol, 'Q' : Q_mat, \
            'V': vVecmesh}
np.save('EWH_sim_NullPower_CF_CN_data_bis', dico)
#%%
