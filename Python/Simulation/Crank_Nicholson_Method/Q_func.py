import numpy as np
import torch


array = np.load('EWH_sim_NullFlow_CN_data.npy', allow_pickle=True).item()

coords = torch.from_numpy(np.stack((array["t"],array["x"]), axis=-1)).float()
Q = torch.from_numpy(array["Q"]).float()
Q_numpy = array["Q"]
#%%
sub_coord = coords[:, 490: 500, :]
coord = sub_coord.reshape([-1, sub_coord.shape[-1]])
t, x = coords[0,:,0], coords[:,0,1]

def Q_(data: torch.Tensor, Q_mat: torch.Tensor):
    samples = data.shape[0]
    Q_value = torch.ones_like(data[:,0:1])
    for n in np.arange(samples):
        i = torch.argmin(torch.abs( x - data[n,1]))
        j = torch.argmin(torch.abs( t - data[n,0]))
        
        Q_value[n] = Q_mat[i,j]
        
    
    return Q_value

Q_samp = Q_(coord, Q)

