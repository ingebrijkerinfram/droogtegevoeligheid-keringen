# -*- coding: utf-8 -*-
"""
Created on Fri Jul 18 14:28:32 2025

@author: inge.brijker
"""

from netCDF4 import Dataset
#from twdm import tqdm 

ds1 = Dataset('hurs_Hd_2050_interp (1).nc', mode='r')
ds2 = Dataset('tas_Hd_2050_interp.nc', mode='r')
ds3 = Dataset('tasmax_Hd_2050_interp.nc', mode='r')
ds4 = Dataset('tasmin_Hd_2050_interp.nc', mode='r')
ds5 = Dataset('pr_Hd_2050_interp.nc', mode='r')
ds6 = Dataset('sfcwind_Hd_2050_interp.nc', mode='r')
ds7 = Dataset('pet_Hd_2050_interp.nc', mode='r')
ds8 = Dataset('rsds_Hd_2050_interp.nc', mode='r')


# Temperatuur (tas)
# Maximum Temperatuur (tasmax)
# Minimum Temperatuur (tasmin)
# Neerslag (pr)
# Wind (sfcwind)
# Relatieve Luchtvochtigheid (hurs)
# Potentiële Evapotranspiratie (pet)
# Zonnestraling (rsds)


# Correct way with netCDF4
ens_array = ds1.variables['ens'][:]
ens_list = ens_array.flatten().tolist()

# hurs_array = ds1.variables['hurs'][:]         # this is a numpy array
# hurs_list = hurs_array.flatten().tolist()

# tas_array = ds2.variables['tas'][:]         # this is a numpy array
# tas_list = tas_array.flatten().tolist()

# tasmax_array = ds3.variables['tasmax'][:]         # this is a numpy array
# tasmax_list = tasmax_array.flatten().tolist()

# tasmin_array = ds4.variables['tasmin'][:]         # this is a numpy array
# tasmin_list = tasmin_array.flatten().tolist()

# pr_array = ds5.variables['pr'][:]         # this is a numpy array
# pr_list = pr_array.flatten().tolist()

# sfcwind_array = ds6.variables['sfcwind'][:]         # this is a numpy array
# sfcwind_list = sfcwind_array.flatten().tolist()

# pet_array = ds7.variables['pet'][:]         # this is a numpy array
# pet_list = pet_array.flatten().tolist()

# rsds_array = ds8.variables['rsds'][:]         # this is a numpy array
# rsds_list = rsds_array.flatten().tolist()


# print("taslist is:", tas_list[:10])
# print("hurslist is:", hurs_list[:10])
# print("tasmax_listt is:", tasmax_list[:10])
# print("tasmin_lis is:", tasmin_list[:10])
# print("pr_listt is:", pr_list[:10])
# print("sfcwind_lis is:", sfcwind_list[:10])
# print("pet_list is:", pet_list[:10])
# print("rsds_listt is:", rsds_list[:10])
print("ens_list is:", ens_list)
