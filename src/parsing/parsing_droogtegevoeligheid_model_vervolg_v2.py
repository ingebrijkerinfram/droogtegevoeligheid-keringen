# -*- coding: utf-8 -*-
"""
Created on Wed Aug 13 15:24:29 2025

@author: inge.brijker
"""

# -*- coding: utf-8 -*-

import pandas as pd
import warnings
import xarray as xr 
from functools import lru_cache

warnings.filterwarnings("ignore")

@lru_cache(maxsize=None)
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

def read_coordinates(grid_path):
    file = pd.read_excel(f'{grid_path}/LatLong_grids.xlsx', sheet_name='NL_grid') # Nieuwe file geconnect van alle gridscellen in NL
    latitude = file['lat']
    longitude = file['long']
    dijk_id = file['ID'] # Was ''Dijk_ID''
    ID = file['ID']
    return latitude, longitude, dijk_id, ID
    
    




