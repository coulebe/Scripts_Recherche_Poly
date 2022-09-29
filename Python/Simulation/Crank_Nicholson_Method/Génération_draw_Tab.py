import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


#%%
start = 0 #Heure
stop = 4.5 # Heure
step_time = 0.5 # Heure
step_flow = 0.05 # L/min
#%%

Time_tab = np.arange(start, stop, step_time, dtype = float)


V = 0

flow_tab = np.zeros_like(Time_tab)
flow_tab[0] = V

middle = int(len(Time_tab)/2)

for i in range(middle+1):
    V  += step_flow
    flow_tab[i] = V

for i in range(middle+1, len(flow_tab)):
    V  -= step_flow
    flow_tab[i] = V
    
#%%
Duration = step_time*60 *np.ones_like(Time_tab) #Minutes

Draw_tab = np.transpose(np.array([Time_tab, Duration, flow_tab]))

def DrawRate(t, DrawTab):
    '''
    Mean to determine the draw rate at each moment given the draw tab
    Parameters:
        t: time (s)
        Tank: Object having the tank Caracteristics(Cf higher for more details)
        DrawTab: Tab having the schedule of all draw with ['Start_time(h), duration(min), Debit(L/min)']
    Outputs: Value of Draw Rate(m/s)
        '''
    I = np.size(DrawTab, 0)
    if I > 0:
        for i in reversed(range(I)):
            if t >= DrawTab[i,0] *3600:
                break
        if(t >= (DrawTab[i,0] *3600)) and (t <= (DrawTab[i,0] *3600 + DrawTab[i,1] *60)):
            V = DrawTab[i,2] 
        else:
            V = 0
    else:
         V = 0
                 
    return V

#%%
T = np.linspace(0, 5*3600, 1000)
Y = np.zeros_like(T)
for i in range(len(T)):
    Y[i] = DrawRate(T[i], Draw_tab)

plt.plot(T/3600, Y)
plt.grid()
plt.xlabel('t(h)')
plt.ylabel('V(L/min)')
plt.title('V(t) (L/min)')
plt.show()
    
#%%
df = pd.DataFrame(data = Draw_tab, columns = ['Start(h)', 'Duration(min)', 'Draw Rate(L/min)'])

df.to_csv('Draw_tab_var_bis.csv')
