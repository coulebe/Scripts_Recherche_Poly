#!/usr/bin/env python
# coding: utf-8

# ## Model discovery for PDE
# 
# In this notebook we use the DeepMoD algorithm applied on the the PDE of an EWH. This version is just a test on the simulation datas. A final version will be developped when we have access to test bench again.
# 
# We start by importing the required DeepMoD functions:

# In[29]:


# General imports
import numpy as np
import torch
import matplotlib.pylab as plt

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
if False: #torch.cuda.is_available():
    device = "cuda"
else:
    device = "cpu"
print(device)


# Before we send the data into the Dataset format, create a plot to get an idea of the data:
# 

# In[30]:


# data = np.load('EWH_sim_NullFlow_CN_data.npy', allow_pickle=True).item()
# fig, ax = plt.subplots()
# im = ax.contourf(data['x'], data['t']/3600, np.real(data['u']))
# ax.set_xlabel('x(m)')
# ax.set_ylabel('t(h)')
# fig.colorbar(mappable=im)

# plt.show()


# We can now define a function that loads the data, makes torch tensors out of it and then returns it in a coords, data format. The shape of the samples will be (t,x) for the input and (u) for the dataset. Ensure that any array is not 1D, so an array with a single feature can be the shape (N,1) using reshape(-1,1) in numpy or unsqueeze(-1) in torch.

# In[31]:


def load_data():
    array = np.load('EWH_sim_NullFlow_CN_data.npy', allow_pickle=True).item()
    coords = torch.from_numpy(np.stack((array["t"],array["x"]), axis=-1)).float()
    data = torch.from_numpy(np.stack((np.real(array["u"]),array["Q"]), axis=-1)).float()
    return coords, data


# We can now pass the function that loads the data to the DeePyMoD Dataset module, which loads the data, preprocesses it, subsamples it and then sends it to the right device.

# In[32]:



# V = torch.from_numpy(np.real(data["V"])).unsqueeze(-1).float()

dataset = Dataset(
    load_data,
    subsampler=Subsample_random,
    subsampler_kwargs={"number_of_samples": 200},
     preprocess_kwargs={
        "normalize_coords": False,
        "normalize_data": False,
    },
    device=device,
)


# Now we need to split our data into a train and test dataloaders for PyTorch

# In[33]:


train_dataloader, test_dataloader = get_train_test_loader(dataset, train_test_split=0.8)


# # Setup a custom library
# 
# In this cell we'll create a custom build library.

# In[34]:


from torch.autograd import grad
from itertools import combinations, product
from functools import reduce
from deepymod.model.library import library_poly, library_deriv


# In[35]:




class Library_EWH(Library):
    """[summary]

    Args:
        Library ([type]): [description]
    """

    def __init__(self) -> None:
        """Calculates the temporal derivative a library/feature matrix consisting of
        1) polynomials up to order poly_order, i.e. u, u^2...
        2) derivatives up to order diff_order, i.e. u_x, u_xx
        3) cross terms of 1) and 2), i.e. $uu_x$, $u^2u_xx$

        Order of terms is derivative first, i.e. [$1, u_x, u, uu_x, u^2, ...$]

        Only works for 1D+1 data. Also works for multiple outputs but in that case doesn't calculate
        polynomial and derivative cross terms.

        Args:
            poly_order (int): maximum order of the polynomial in the library
            diff_order (int): maximum order of the differentials in the library
        """

        super().__init__()


    def library(
        self, input: Tuple[torch.Tensor, torch.Tensor]
    ) -> Tuple[TensorList, TensorList]:
        """Compute a 2D library up to given polynomial order with second order derivatives
         i.e. for poly_order=1: [$1, u_x, u_Q, u_{xx}, u_{QQ}, u_{xQ}, Q$] because is considered as an input of our neural network

        Args:
            input (Tuple[torch.Tensor, torch.Tensor]): A prediction u (n_samples, n_outputs) and spatiotemporal locations (n_samples, 3).(Q is the thrid input)

        Returns:
            Tuple[TensorList, TensorList]: The time derivatives and the thetas
            computed from the library and data.
        """
        prediction, data = input


     
        
        u = prediction[:,0:1]
        # # Gradients
        du = grad(
            u,
            data,
            grad_outputs=torch.ones_like(u),
            create_graph=True,
        )[0]
        u_t = du[:, 0:1]
        u_x = du[:, 1:2]

        du2 = grad(
            u_x, data, grad_outputs=torch.ones_like(u), create_graph=True
        )[0]
        u_xx = du2[:, 1:2]
        # time_deriv, theta_1D = library_deriv(data= data, prediction= u,max_order= 2)
        
        Q = prediction[:,1:2]


        theta = torch.cat((torch.ones_like(u), u,u_x, u_xx, Q), dim=1)

        # dT = grad(
        #     prediction,
        #     data,
        #     grad_outputs=torch.ones_like(prediction),
        #     create_graph=True,
        # )[0]
        # time_deriv = dT[:, 0:1]

        return [u_t, Q], [theta, theta]


# ## Configuring DeepMoD
# 
# Configuration of the function approximator: Here the first argument is the number of input and the last argument the number of output layers.

# In[36]:


network = NN(2, [50, 50, 50, 50], 2)


# Configuration of the library function: We select the custom build library we created earlier
# * [$1, u, u_x, u_{xx}, Q$] 

# In[37]:


library = Library_EWH()


# Configuration of the sparsity estimator and sparsity scheduler used. In this case we use the most basic threshold-based Lasso estimator and a scheduler that asseses the validation loss after a given patience. If that value is smaller than 1e-8, the algorithm is converged.

# In[38]:


estimator = Threshold(0.5)
sparsity_scheduler = TrainTestPeriodic(periodicity=50, patience=200, delta=1e-8)


# Configuration of the sparsity estimator

# In[39]:


constraint = LeastSquares()


# Now, we instantiate the model and select the optimizer

# In[40]:


model = DeepMoD(network, library, estimator, constraint).to(device)


# Defining optimizer

# In[41]:


optimizer = torch.optim.Adam(
    model.parameters(), betas=(0.99, 0.99), amsgrad=True, lr=1e-3
)


# ## Run DeepMoD
# 
# We can now run DeepMoD using all the options we have set and the training data. We need to slightly preprocess the input data for the derivatives:

# In[42]:


train(
    model,
    train_dataloader,
    test_dataloader,
    optimizer,
    sparsity_scheduler,
    log_dir="/data/deepymod/coupled_new/",
    max_iterations=100000,
    delta=1e-8,
    patience=100
)


# Now that DeepMoD has converged, it has found the following coefficients to not be zero: 

# In[ ]:


model.sparsity_masks


# And it found the following coefficients.

# In[ ]:


model.estimator_coeffs()


# Then the absolute value of the coefficients
# * [$1, u, u_{xx}, Q$] (with poly order = 1)
# 

# In[ ]:


print(model.library.norms)

