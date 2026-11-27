# -*- coding: utf-8 -*-
"""
Created on Wed Aug 13 15:24:29 2025

@author: inge.brijker
"""

# -*- coding: utf-8 -*-

import numpy as np
import warnings
from src.padenv2noord import bepaal_paden, controleer_klimaatdata
from src.parsing.parsing_droogtegevoeligheid_model_vervolg_v2compare_noord  import read_parameters, read_coordinates
from src.processing.processing_droogtegevoeligheid_model_vervolg_v2compare_noord  import load_weather_from_nc, watercalculations, verwelkingspunt, plot, analyse_per_jaar, select_extreme_cases, main_per_ensembles, main_all_ensembles
from src.processing.processing_droogtegevoeligheid_model_vervolg_v2compare_noord import ModelOptions
from itertools import product
from dataclasses import asdict
import json
import time
import sys
warnings.filterwarnings("ignore")



# ================================ INSTELLINGEN ================================
# Alle combinaties van deze lijsten worden na elkaar doorgerekend. Je kan varieren in welke data je allemaal wilt rekenen
SCENARIOS    = ["Hd_2050"]   # "ref", "Hn_2050", "Hd_2050", "Ln_2100", "Ld_2100", "Hn_2100", "Hd_2100"
if __name__ == "__main__" and len(sys.argv) > 1:
    SCENARIOS = sys.argv[1].split(",")
HELLINGEN    = [1/2, 1/3, 1/4, 1/5]  # 1/2, 1/3, 1/4, 1/5
GRONDSOORTEN = ["klei", "zand"]      # "zand", "klei"

# Onafhankelijke schakelaars: meerdere waarden = alle combinaties.
# True/False voor bodem: True is de reeds aangeleverde parameterisatie.
BODEM_NIEUW = [False]
MAAIEN = [False]
NIEUWE_MORTALITEIT = [False]
LUCHTDRUK_KPA = 101.3
MAAIDATA_ZAND = ((5, 31), (8, 31))   # maand, dag (jaarlijks)
MAAIDATA_KLEI = ((5, 31), (8, 31))
MAAI_KC_FACTOR = 0.4473684210526316  # (0.275 + 0.15) / 0.95, ARVIOR-proxy
MAAI_HERSTEL = 0.005625      # fit op ARVIOR Kcb_target + nominale Ke; voorlopig
S50 = 68.16
CURVEBREEDTE = 1.5
NABIJ_NUL_PCT = 1.0         # operationele categorie, geen waterveiligheidsnorm

# Optie voor testrun
TESTRUN = False   # True = alleen 5 gridcellen en 2 jaar (output in map 'Testrun'), False = alles
# ==============================================================================



def run(em_scen="Hd_2050", helling=1/2, grondsoort="klei", testrun=True, options=None):
    """
    Deze functie draait model voor 1 combinatie. De volgende parameters kunnen gevarieerd worden:
 
    em_scen    : "ref", "Hn_2050", "Hd_2050", "Ln_2100", "Ld_2100", "Hn_2100", "Hd_2100"
    helling    : 1/2, 1/3, 1/4, 1/5
    grondsoort : "zand" of "klei"
    testrun    : True  = alleen de eerste 5 gridcellen en 2 jaar, output in de map 'Testrun'
                 False = alles, output direct in de map Output
    """
    options = options or ModelOptions()
    # Starttijd meten
    start_time = time.time()
    
    
    # Automatisch bepalen van paden. Dit gebeurt via src/paden.py waarin projectmap 000687 wordt aangehaald.
    paden = bepaal_paden() # Zie import van functie bepaal_paden. Deze is nieuw.
    path = paden["klimaatdata"]
    path_grid = paden["locaties"].parent
    input_path = paden["klimaatdata"]
    if testrun:
        output_path = paden["output"] / "Testrun" # Zodat testen en definitieve runs gescheiden blijven
    else:
        output_path = paden["output"]
    
    output_path = output_path / options.label()
    output_path.mkdir(parents=True, exist_ok=True)
    config_path = output_path / "modelinstellingen.json"
    config_json = json.dumps(asdict(options), indent=2, ensure_ascii=False)
    if config_path.exists() and config_path.read_text(encoding="utf-8") != config_json:
        raise ValueError(f"Andere parameters bestaan al in {output_path}; kies een aparte uitvoermap.")
    config_path.write_text(config_json, encoding="utf-8")

    ens_idx = 1  # realisatie, stel hier de gewenste locatie in (coördinaten of index)
    
    
    # Controle van input files
    controleer_klimaatdata(paden, em_scen)   # Met deze functie stopt het script met een  melding als er een bestand ontbreekt
    
    
    # Inlezen van coordinaten
    lat = read_coordinates(grid_path=path_grid)[0][:] #testrun, anders [:] voor alles rekenen, voor alle drie de regels.
    lon = read_coordinates(grid_path=path_grid)[1][:]
    dijkvak_id = read_coordinates(grid_path=path_grid)[2][:]
    
    
    # Jaren voor VOLLEDIGE run
    if em_scen == 'ref':
        jaar = "2005"
        jarenrange = range(1991, 2021)
    else:
        jaar = em_scen.split("_")[1]
        if jaar == "2050":
            jarenrange = range(2036, 2066)
        elif jaar == "2100":
            jarenrange = range(2086, 2116)
    
    # Jaren voor TESTRUN: alleen de eerste 5 gridcellen en de eerste 2 jaar
    if testrun:
        lat = lat[:5]
        lon = lon[:5]
        dijkvak_id = dijkvak_id[:5]
        jarenrange = jarenrange[:2]
    
    if helling == 1/2:
        cot = 'cot2'
    elif helling == 1/3:
        cot = 'cot3'
    elif helling == 1/4:
        cot = 'cot4'
    elif helling == 1/5:
        cot = 'cot5'
    
    
    # Runnen van proces berekeningen
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
        output_path=output_path,
        options=options
    )
    
    
    # Eindtijd meten
    end_time = time.time()
    elapsed_time = end_time - start_time
    print(f"Script duurt {elapsed_time:.4f} seconden om te runnen.")



# Alle combinaties; afzonderlijke uitvoermappen voorkomen overschrijven.
if __name__ == "__main__":
    for bodem, maaien, mortaliteit in product(
            BODEM_NIEUW, MAAIEN, NIEUWE_MORTALITEIT):
        options = ModelOptions(
            use_revised_soil=bodem, enable_mowing=maaien,
            enable_ks_mortality=mortaliteit,
            air_pressure_kpa=LUCHTDRUK_KPA,
            mow_dates_sand=MAAIDATA_ZAND, mow_dates_clay=MAAIDATA_KLEI,
            mow_kc_fraction=MAAI_KC_FACTOR, mow_recovery_rate=MAAI_HERSTEL,
            s50=S50, width=CURVEBREEDTE, near_zero_pct=NABIJ_NUL_PCT,
        )
        for em_scen in SCENARIOS:
            for grondsoort in GRONDSOORTEN:
                for helling in HELLINGEN:
                    run(em_scen, helling, grondsoort, TESTRUN, options)
