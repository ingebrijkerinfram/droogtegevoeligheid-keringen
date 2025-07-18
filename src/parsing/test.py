# -*- coding: utf-8 -*-
"""
Created on Fri Jul 18 14:28:32 2025

@author: inge.brijker
"""

import xarray as xr

# Open alle bestanden met xarray
ds_hurs = xr.open_dataset('hurs_Hd_2050_interp (1).nc')
# ds_tas = xr.open_dataset('tas_Hd_2050_interp.nc')
# ds_tasmax = xr.open_dataset('tasmax_Hd_2050_interp.nc')
# ds_tasmin = xr.open_dataset('tasmin_Hd_2050_interp.nc')
# ds_pr = xr.open_dataset('pr_Hd_2050_interp.nc')
# ds_sfcwind = xr.open_dataset('sfcwind_Hd_2050_interp.nc')
# ds_pet = xr.open_dataset('pet_Hd_2050_interp.nc')
# ds_rsds = xr.open_dataset('rsds_Hd_2050_interp.nc')

# Selecteer alleen ens = 1 (de tweede realisatie)
hurs = ds_hurs['hurs'].isel(ens=1)
# tas = ds_tas['tas'].isel(ens=1)
# tasmax = ds_tasmax['tasmax'].isel(ens=1)
# tasmin = ds_tasmin['tasmin'].isel(ens=1)
# pr = ds_pr['pr'].isel(ens=1)
# sfcwind = ds_sfcwind['sfcwind'].isel(ens=1)
# pet = ds_pet['pet'].isel(ens=1)
# rsds = ds_rsds['rsds'].isel(ens=1)

# Zet arrays om naar numpy (optioneel naar list)
hurs_list = hurs.values.flatten().tolist()
# tas_list = tas.values.flatten().tolist()
# tasmax_list = tasmax.values.flatten().tolist()
# tasmin_list = tasmin.values.flatten().tolist()
# pr_list = pr.values.flatten().tolist()
# sfcwind_list = sfcwind.values.flatten().tolist()
# pet_list = pet.values.flatten().tolist()
# rsds_list = rsds.values.flatten().tolist()

# Latitude-array ophalen als lijst
lat_list = ds_hurs['lat'].values.flatten().tolist()

# Voorbeeldoutput
# print("tas_list[:10]:", tas_list[:10])
print("hurs_list[:10]:", hurs_list[:10])
# print("tasmax_list[:10]:", tasmax_list[:10])
# print("tasmin_list[:10]:", tasmin_list[:10])
# print("pr_list[:10]:", pr_list[:10])
# print("sfcwind_list[:10]:", sfcwind_list[:10])
# print("pet_list[:10]:", pet_list[:10])
# print("rsds_list[:10]:", rsds_list[:10])
# print("lat_list[:5]:", lat_list[:5])