# -*- coding: utf-8 -*-
"""
Created on Wed Aug 13 15:24:29 2025

@author: inge.brijker
"""

# -*- coding: utf-8 -*-

# import matplotlib.pyplot as plt
import numpy as np
# import pandas as pd
import warnings
# import xarray as xr 
from src.parsing.parsing_droogtegevoeligheid_model_MS import read_parameters, read_coordinates
from src.processing.processing_droogtegevoeligheid_model_MS1_LPC import load_weather_from_nc, watercalculations, verwelkingspunt, plot, analyse_per_jaar, select_extreme_cases, main_per_ensembles, main_all_ensembles
import time
warnings.filterwarnings("ignore")

# Starttijd meten
start_time = time.time()


# Handmatige input
path = r'C:\Users\ingeb\Documents\Infram\Droogtegevoeligheid\Scenario Data'
path_grid = r''C:\Users\ingeb\Documents\Infram\Droogtegevoeligheid\Grid'
input_path = r'C:\Users\ingeb\Documents\Infram\Droogtegevoeligheid\Scenario Data'
output_path = r'C:\Users\ingeb\Documents\Infram\Droogtegevoeligheid\Resultaten'


grondsoort = "zand"     # "zand" of "klei"
helling = 1/2           #  1/2, 1/3, 1/4, 1/5
em_scen = "Hn_2050"     # "ref", "Hn_2050", "Hd_2050", "Ln_2100", "Ld_2100", "Hn_2100", "Hd_2100"

#Stel hier de gewenste locatie in (coördinaten of index)
ens_idx = 1  #  realisatie
#TODO
#lat, lon en dijkvak_id controleren/kloppend maken

# lat = np.array([read_coordinates(grid_path=path_grid)[0][2], read_coordinates(grid_path=path_grid)[0][3], read_coordinates(grid_path=path_grid)[0][4]])  # kies gewenste lat index
# lon = np.array([read_coordinates(grid_path=path_grid)[1][2], read_coordinates(grid_path=path_grid)[1][3], read_coordinates(grid_path=path_grid)[1][4]])  # kies gewenste lon index
# dijkvak_id = np.array([str(read_coordinates(grid_path=path_grid)[2][2]), str(read_coordinates(grid_path=path_grid)[2][3]), str(read_coordinates(grid_path=path_grid)[2][4])])

lat = read_coordinates(grid_path=path_grid)[0][:]
lon = read_coordinates(grid_path=path_grid)[1][:]
dijkvak_id = read_coordinates(grid_path=path_grid)[2][:]

if em_scen == 'ref':
    jaar = "2005"
    jarenrange = range(1991, 2021)
else:
    jaar = em_scen.split("_")[1]
    if jaar == "2050":
        jarenrange = range(2036, 2066)
    elif jaar == "2100":
        jarenrange = range(2086, 2116)
    
if helling == 1/2:
    cot = 'cot2'
elif helling == 1/3:
    cot = 'cot3'
elif helling == 1/4:
    cot = 'cot4'
elif helling == 1/5:
    cot = 'cot5'

# main_per_ensembles(ens_idx, lat, lon, jaar)
main_all_ensembles(
    path=path,
    em_scen=em_scen,
    lat=lat,
    lon=lon,
    year=jaar,
    jaren=jarenrange,
    grondsoort=grondsoort,
    dijkvak_id=dijkvak_id,
    helling=helling, # <-- toegevoegd
    cot=cot,
    output_path=output_path
)

# Eindtijd meten
end_time = time.time()
elapsed_time = end_time - start_time
print(f"Script duurt {elapsed_time:.4f} seconden om te runnen.")