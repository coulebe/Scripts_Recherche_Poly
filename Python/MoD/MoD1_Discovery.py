
# coding: utf-8

# # Model discovery MoD1

# In this notebook, we'll discover model of EWH1 (OSO300L) form the TB datat in MoD1, i.e., we have Temperature, draw debit(m³/s), and power consumption(W)

# In[2]:


# General imports
import numpy as np
import torch
import matplotlib.pylab as plt
import tkinter as tk
from tkinter import filedialog
import os
import pandas as pd

#Function file import
import DeepMod_Functions as DD_func
# DeepMoD functions


from deepymod import DeepMoD
from deepymod.data import Dataset, get_train_test_loader
from deepymod.data.samples import Subsample_random
from deepymod.model.func_approx import NN
from deepymod.model.constraint import LeastSquares, GradParams, Ridge
from deepymod.model.sparse_estimators import Threshold, PDEFIND
from deepymod.training import train
from deepymod.training.sparsity_scheduler import TrainTestPeriodic
from scipy.io import loadmat


import torch
from torch.autograd import grad
from itertools import combinations
from functools import reduce
from typing import Tuple
from deepymod.utils.types import TensorList
from deepymod import Library

from scipy.integrate import odeint


# Settings for reproducibility
np.random.seed(40)
torch.manual_seed(0)

# Configuring GPU or CPU
if torch.cuda.is_available():
    device = "cuda"
    torch.cuda.empty_cache()
else:
    device = "cpu"
print(device)


# Before we send the data into the Dataset format, create a plot to get an idea of the data:
# 

# In[2]:


#Asking the filename to the user(It will be a .npy file for the moment we'll see after if we need to import other fomat as .mat or ohter, think of a preprocess script for the essay datas at Canmet)
application_window = tk.Tk() 
fTyp = [("fichier de données (*.npy)", "*.npy")]
files = filedialog.askopenfilenames(parent=application_window,
                                    initialdir=os.getcwd(),
                                    title="Please select your npy file containing the datas and its coords:",
                                    filetypes=fTyp)


application_window.destroy()


for file_name in files:
    test_name = os.path.splitext(os.path.basename(file_name))[0]

    # We can now pass the function that loads the data to the DeePyMoD Dataset module, which loads the data, preprocesses it, subsamples it and then sends it to the right device.

    # 


    dataset = Dataset(
        load_function= DD_func.load_MoD1,
        subsampler=Subsample_random,
        subsampler_kwargs={"number_of_samples": 500000, 
        },
        load_kwargs= {"file": file_name},
        preprocess_kwargs={
            "normalize_coords": True,
            "normalize_data":  True
        }, 
        device=device,
        apply_normalize= DD_func.custom_normalize,
    )
    
    #Now we need to split our data into a train and test dataloaders for PyTorch
    train_dataloader, test_dataloader = get_train_test_loader(dataset, train_test_split=0.8)


    # ## Configuring DeepMoD
    # 
    # Configuration of the function approximator: Here the first argument is the number of input and the last argument the number of output layers.

    # 


    network = NN(5, [50, 50, 50, 50], 1)


    # Configuration of the library function: We select the custom build library we created earlier
    # * [$1, (T - T_{a}), \frac{\partial^{2} T}{\partial x^{2}}, V \frac{\partial T}{\partial x}, P$] 

    # 


    Library = DD_func.Library_MoD1()


    # Configuration of the sparsity estimator and sparsity scheduler used. In this case we use the most basic threshold-based Lasso estimator and a scheduler that asseses the validation loss after a given patience. If that value is smaller than 1e-8, the algorithm is converged.

    # 


    #Sparsity scheduler
    estimator = Threshold(1e-4)
    sparsity_scheduler = TrainTestPeriodic(periodicity=50, patience=200, delta=1e-5)

    #Configuration of the sparsity estimator
    constraint = LeastSquares()

    # Now, we instantiate the model and select the optimizer
    model = DeepMoD(network, Library, estimator, constraint).to(device)

    #Defining optimizer
    optimizer = torch.optim.Adam(
        model.parameters(), betas=(0.99, 0.99), amsgrad=True, lr=1e-3
    )


    # ## Run DeepMoD
    # 
    # We can now run DeepMoD using all the options we have set and the training data. We need to slightly preprocess the input data for the derivatives:

    # 


    train(
        model,
        train_dataloader,
        test_dataloader,
        optimizer,
        sparsity_scheduler,
        log_dir= os.getcwd() + "\\Logdir\\MoD1\\" ,
        max_iterations=100000,
        delta=1e-3,
        patience=200
    )


    # Now that DeepMoD has converged, it has found the following coefficients to not be zero: 

    # 


    model.sparsity_masks


    # And it found the following coefficients.




    model.estimator_coeffs()


    # Then the absolute value of the coefficients
    # * [$1, (T - T_{a}), \frac{\partial^{2} T}{\partial x^{2}}, V \frac{\partial T}{\partial x}, P$]
    # 



    print(model.library.norms[0].cpu().numpy() * model.estimator_coeffs()[0][:,0])
    coeff_norm = model.library.norms[0].cpu().numpy() * model.estimator_coeffs()[0][:,0]


    # We have now to denormalize it a second time because, we just have the relation for the normalized(data and coords) terms of the equation




    #Max and max of coords and datas
    coords, data = DD_func.load_MoD1(file_name)
    T_min, T_max = data.view(-1, data.shape[-1]).min(dim=0).values.numpy()[0], data.view(-1, data.shape[-1]).max(dim=0).values.numpy()[0]
    t_min, t_max = coords[0, :, 0].min(dim = 0).values.numpy(), coords[0, :, 0].max(dim = 0).values.numpy()
    x_min, x_max = coords[:, 0, 1].min(dim = 0).values.numpy(), coords[:, 0, 1].max(dim = 0).values.numpy()
    V_min, V_max = coords[:, :, 2].view(-1, 1).min(dim = 0).values.numpy()[0], coords[:, :, 2].view(-1, 1).max(dim = 0).values.numpy()[0]
    P_min, P_max = coords[:, :, 3].view(-1, 1).min(dim = 0).values.numpy()[0], coords[:, :, 3].view(-1, 1).max(dim = 0).values.numpy()[0]
    Ta_min, Ta_max = coords[:, :, 4].view(-1, 1).min(dim = 0).values.numpy()[0], coords[:, :, 4].view(-1, 1).max(dim = 0).values.numpy()[0]

    #Result dataframe with all the denormalized coefficients
    result = pd.DataFrame(columns= ["1", "T", "Ta", "T_xx", "V*T_x", "T_x","P"] )
    result.index.name = "Test_name"

    result.loc[test_name, "1"] = (T_max-T_min) * (coeff_norm[0] - coeff_norm[1]*((T_min/(T_max-T_min)) - (Ta_min/(Ta_max-Ta_min))) - (coeff_norm[4]*P_min/(P_max-P_min)) ) / (t_max-t_min)
    result.loc[test_name,"T"] = coeff_norm[1] / (t_max-t_min)
    result.loc[test_name,"Ta"] = -coeff_norm[1] * (T_max-T_min) / ((t_max-t_min)*(Ta_max-Ta_min))# Minus added bcause we have smthg(T-Ta)
    result.loc[test_name, "T_xx"] = coeff_norm[2] * ((x_max - x_min)**2) / (t_max-t_min)
    result.loc[test_name, "V*T_x"] = coeff_norm[3]*(x_max - x_min) / ((t_max-t_min) * (V_max - V_min))
    result.loc[test_name, "T_x"] = -coeff_norm[3]*(x_max - x_min)*V_min / ((t_max-t_min) * (V_max - V_min))
    result.loc[test_name, "P"] = coeff_norm[4]*(T_max-T_min) / ((t_max-t_min)*(P_max-P_min))




    prev = pd.read_csv("MoD1_results.csv", index_col=0)
    result = pd.concat([prev, result])
    result.to_csv("MoD1_results.csv")

