import numpy as np

data =  np.load('EWH_sim_NullFlow_CN_data.npy', allow_pickle=True).item()
print(np.shape(data["Q"]))
print(np.shape(np.argwhere(data["Q"] == 0)))
# %%
