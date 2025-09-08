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
from src.parsing.parsing_droogtegevoeligheid_model_MS import read_parameters
from src.processing.processing_droogtegevoeligheid_model_MS1 import load_weather_from_nc, watercalculations, verwelkingspunt, living_plant_cover, plot, analyse_per_jaar, select_extreme_cases, main_per_ensembles, main_all_ensembles


warnings.filterwarnings("ignore")

# Handmatige input
path = r'C:\Users\marloes.slokker\Documents\Droogtemodel'
input_path = r'C:\Users\marloes.slokker\Documents\Droogtemodel'
# output_path = r'C:\Users\marloes.slokker\Infram BV\Infram Projecten - 000370 Droogtegevoeligheid keringen RWS\Uitvoering\Resultaten'
grondsoort = "zand"  # 'zand' of 'klei'
dijkvak_id = "Europoortkering"
em_scen = 'ref'
year = '2005'


#Stel hier de gewenste locatie in (coördinaten of index)
ens_idx = 1  #  realisatie
lon = 4.28  # kies gewenste lon index
lat = 51.9  # kies gewenste lat index
jaar = 2005
jarenrange = range(1991, 2021) 



# main_per_ensembles(ens_idx, lat, lon, jaar)
main_all_ensembles(
    path=path,
    em_scen=em_scen,
    lat=lat,
    lon=lon,
    jaren=jarenrange,
    dijkvak_id=dijkvak_id,   # <-- toegevoegd
    output_path=r"C:\Users\marloes.slokker\Infram BV\Infram Projecten - 000370 Droogtegevoeligheid keringen RWS\Uitvoering\Resultaten"
)