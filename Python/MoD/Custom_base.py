""" Contains a custom version of the base class for the Dataset  (1 and 2 dimensional) and a function
     that takes a Pytorch tensor and converts it to a numpy array"""

import torch
import numpy as np
from numpy import ndarray
from numpy.random import default_rng
from deepymod.data.samples import Subsampler

from abc import ABC, ABCMeta, abstractmethod

from deepymod.data import Dataset, get_train_test_loader


class Dataset_EWH(Dataset):
    def __init__(
        self,
        load_function,
        Q, 
        V,
        apply_normalize=None,
        apply_noise=None,
        apply_shuffle=None,
        shuffle=True,
        subsampler: Subsampler = None,
        load_kwargs: dict = {},
        preprocess_kwargs: dict = {
            "random_state": 42,
            "noise_level": 0.0,
            "normalize_coords": False,
            "normalize_data": False,
        },
        subsampler_kwargs: dict = {},
        device: str = None,
    ):
        '''
        Same class as Dataset BUT, we'll introduce Q and V matrice for the prupose of our model
        Args:
            load_function (func):Must return torch tensors in the format coordinates, data
            Q: Matrix containing the  value of Q(x,t) in a tensor form(same shape as u)
            V: Matrix containing the  value of V(x,t)(even if V is the same along x) in a tensor form(same shape as u)
            shuffle (bool, optional): Shuffle the data. Defaults to True.
            apply_normalize (func)
            subsampler (Subsampler, optional): Add some subsampling function. Defaults to None.
            load_kwargs (dict, optional): kwargs to pass to the load_function. Defaults to {}.
            preprocess_kwargs (dict, optional): (optional) arguments to pass to the preprocess method. Defaults to { "random_state": 42, "noise_level": 0.0, "normalize_coords": False, "normalize_data": False, }.
            subsampler_kwargs (dict, optional): (optional) arguments to pass to the subsampler method. Defaults to {}.
            device (str, optional): which device to send the data to. Defaults to None.
        '''
        super()._init_(
        load_function,
        apply_normalize,
        apply_noise,
        apply_shuffle,
        shuffle,
        subsampler,
        load_kwargs,
        preprocess_kwargs,
        subsampler_kwargs,
        device,
        )

        self.Q = Q
        self.V = V
        