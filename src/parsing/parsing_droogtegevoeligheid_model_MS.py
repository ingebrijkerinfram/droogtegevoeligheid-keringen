# -*- coding: utf-8 -*-
"""
Created on Wed Aug 13 15:24:29 2025

@author: inge.brijker
"""

# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import warnings
import xarray as xr 

warnings.filterwarnings("ignore")

# Handmatige input
# path = r'C:\Users\marloes.slokker\Documents\Droogtemodel'
# grondsoort = "zand"  # 'zand' of 'klei'
# dijkvak_id = "Europoortkering"
# em_scen = 'ref'
# year = '2005'

def read_parameters(path, em_scen):
    ds_hurs = xr.open_dataset(f'{path}/hurs_{em_scen}_interp.nc', engine="netcdf4")
    ds_tas = xr.open_dataset(f'{path}/tas_{em_scen}_interp.nc', engine="netcdf4")
    ds_tasmax = xr.open_dataset(f'{path}/tasmax_{em_scen}_interp.nc', engine="netcdf4")
    ds_tasmin = xr.open_dataset(f'{path}/tasmin_{em_scen}_interp.nc', engine="netcdf4")
    ds_pr = xr.open_dataset(f'{path}/pr_{em_scen}_interp.nc', engine="netcdf4")
    ds_sfcwind = xr.open_dataset(f'{path}/sfcwind_{em_scen}_interp.nc', engine="netcdf4")
    ds_pet = xr.open_dataset(f'{path}/pet_{em_scen}_interp.nc', engine="netcdf4")
    ds_rsds = xr.open_dataset(f'{path}/rsds_{em_scen}_interp.nc', engine="netcdf4")
    
    return ds_hurs, ds_tas, ds_tasmax, ds_tasmin, ds_pr, ds_sfcwind, ds_pet, ds_rsds




