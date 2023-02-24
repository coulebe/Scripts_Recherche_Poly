#Python script containing the classes and function needed for the parametert identification

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
from deepymod.model.library import library_poly, library_deriv
from scipy.io import loadmat


import torch
from torch.autograd import grad
from itertools import combinations, product
from functools import reduce
from typing import Tuple
from deepymod.utils.types import TensorList
from deepymod import Library

from scipy.integrate import odeint



#%%
#General functions

def custom_normalize(X):
        """minmax Normalize the data. Per feature
        Args:
            X (torch.tensor): data to be minmax normalized
        Returns:
            (torch.tensor): minmaxed data"""
        X_norm = (X - X.view(-1, X.shape[-1]).min(dim=0).values) / (
            X.view(-1, X.shape[-1]).max(dim=0).values
            - X.view(-1, X.shape[-1]).min(dim=0).values
        ) 
        return X_norm
#%%
# MoD1
#Functions
def load_MoD1(file):
    array = np.load(file, allow_pickle=True).item()
    coords = torch.from_numpy(np.stack((array["t"],array["x"], array["V"], array["P"], array["Ta"]), axis=-1)).float()
    data = torch.from_numpy(np.real(array["T"])).unsqueeze(-1).float()
    return coords, data


#Set custom library
class Library_MoD1(Library):
    """[summary]

    Args:
        Library ([type]): [description]
    """

    def __init__(self) -> None:
        # """Calculates the temporal derivative a library/feature matrix consisting of
        # 1) polynomials up to order poly_order, i.e. u, u^2...
        # 2) derivatives up to order diff_order, i.e. u_x, u_xx
        # 3) cross terms of 1) and 2), i.e. $uu_x$, $u^2u_xx$

        # Order of terms is derivative first, i.e. [$1, u_x, u, uu_x, u^2, ...$]

        # Only works for 1D+1 data. Also works for multiple outputs but in that case doesn't calculate
        # polynomial and derivative cross terms.

        # Args:
        #     poly_order (int): maximum order of the polynomial in the library
        #     diff_order (int): maximum order of the differentials in the library
        # """

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

        # Gradients
        du = grad(
            prediction,
            data,
            grad_outputs=torch.ones_like(prediction),
            create_graph=True,
        )[0]
        u_t = du[:, 0:1]
        u_x = du[:, 1:2]

        du2 = grad(
            u_x, data, grad_outputs=torch.ones_like(prediction), create_graph=True
        )[0]
        u_xx = du2[:, 1:2]

        #We add V and multiply it time x, add Q and Ta
        V = data[:, 2:3]
        P = data[:, 3:4]
        Ta = data[:, 4:5]

        theta = torch.cat((torch.ones_like(prediction),  prediction-Ta ,u_xx, u_x * V, P), dim=1)

        return [u_t], [theta]


