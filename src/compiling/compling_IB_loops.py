import numpy as np
import warnings
import time
from src.parsing.parsing_droogtegevoeligheid_model_IB import read_parameters, read_coordinates
from src.processing.processing_droogtegevoeligheid_model_IB_LPC import (
    load_weather_from_nc,
    watercalculations,
    verwelkingspunt,
    plot,
    analyse_per_jaar,
    select_extreme_cases,
    main_per_ensembles,
    main_all_ensembles,
)

warnings.filterwarnings("ignore")

#Starttijd van run
start_time = time.time()

# Handmatige input
path = r'C:\Users\inge.brijker\Documents\Droogtemodel'
path_grid = r'C:\Users\inge.brijker\Infram BV\Infram Projecten - 000370 Droogtegevoeligheid keringen RWS\Uitvoering\Methode'
input_path = r'C:\Users\inge.brijker\Documents\Droogtemodel'
output_path = r'C:\Users\inge.brijker\Infram BV\Infram Projecten - 000370 Droogtegevoeligheid keringen RWS\Uitvoering\Resultaten\Zuid'

#Coordinaten en dijkvakken
lat = read_coordinates(grid_path=path_grid)[0][:]
lon = read_coordinates(grid_path=path_grid)[1][:]
dijkvak_id = read_coordinates(grid_path=path_grid)[2][:]

# Lijsten met scenario's, grondsoorten en hellingen
em_scen_lijst = ["ref", "Hn_2050", "Hd_2050", "Ln_2100", "Ld_2100", "Hn_2100", "Hd_2100"]
grondsoorten = ["zand", "klei"]
hellingen = [1/2, 1/3, 1/4, 1/5]   

for em_scen in em_scen_lijst:
    if em_scen == 'ref':
        jaar = "2005"
        jarenrange = range(1991, 2021)
    else:
        jaar = em_scen.split("_")[1]
        if jaar == "2050":
            jarenrange = range(2036, 2066)
        elif jaar == "2100":
            jarenrange = range(2086, 2116)
        else:
            raise ValueError(f"Onbekend jaartal voor emissiescenario: {em_scen}")

    for grondsoort in grondsoorten:
        for helling in hellingen:
            if helling == 1/2:
                cot = 'cot2'
            elif helling == 1/3:
                cot = 'cot3'
            elif helling == 1/4:
                cot = 'cot4'
            elif helling == 1/5:
                cot = 'cot5'
            else:
                raise ValueError(f"Onbekende helling: {helling}")

            print("\n" + "=" * 80)
            print(f"Start run: scenario = {em_scen}, grondsoort = {grondsoort}, helling = {helling} ({cot})")
            print("=" * 80)

            main_all_ensembles(
                path=path,
                em_scen=em_scen,
                lat=lat,
                lon=lon,
                year=jaar,          
                jaren=jarenrange,   
                grondsoort=grondsoort,
                dijkvak_id=dijkvak_id,
                helling=helling,
                cot=cot,
                output_path=output_path
            )

end_time = time.time()
elapsed_time = end_time - start_time
print(f"\nScript duurt {elapsed_time:.2f} seconden om te runnen.")
print("Alle scenario's, grondsoorten en hellingen zijn doorgerekend.")
