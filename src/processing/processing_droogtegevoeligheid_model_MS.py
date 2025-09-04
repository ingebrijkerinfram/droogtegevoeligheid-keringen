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

warnings.filterwarnings("ignore")

# # Handmatige input
path = r'C:\Users\marloes.slokker\Documents\Droogtemodel'
grondsoort = "zand"  # 'zand' of 'klei'
dijkvak_id = "Europoortkering"
em_scen = 'ref'
# year = '2005'

# #Stel hier de gewenste locatie in (coördinaten of index)
ens_idx = 1  #  realisatie
lat = 4.28  # kies gewenste lat index
lon = 51.9  # kies gewenste lon index
jaar = 2005
jarenrange = range(1991, 2021) 

#DikeGrass crop parameters getest door Thomas (voor nu)
#hier kan een if statement: if "zand", then:, else (clay)

SOIL_PARAMS = {
    "zand": {
        "Kcini": 0.5, "Kcmid": 0.8, "Kclate": 0.6,
        "h": 0.4, "Zr": 0.3,
        "Ofc": 0.2, "Owp": 0.05,
        "ini": 30, "dev": 0, "mid": 365, "late": 0,
        "pcrop": 0.45
    },
    "klei": {
        "Kcini": 0.5, "Kcmid": 0.9, "Kclate": 0.6,
        "h": 0.3, "Zr": 0.3,
        "Ofc": 0.35, "Owp": 0.17,
        "ini": 30, "dev": 0, "mid": 365, "late": 0,
        "pcrop": 0.4
    }
}

ds_hurs, ds_tas, ds_tasmax, ds_tasmin, ds_pr, ds_sfcwind, ds_pet, ds_rsds = read_parameters(path, em_scen)
dates = ds_hurs['time'].values

def load_weather_from_nc(ens=ens_idx, lat=lat, lon=lon, jaar=jaar):
    # Tijd filteren
    time_filter = ds_hurs['time'].dt.year == jaar

    # Tijdstappen ophalen (gefilterd)
    dates = ds_hurs['time'].where(time_filter, drop=True).values

    # Waarden ophalen met .sel(..., method='nearest')
    T = ds_tas['tas'].sel(ens=ens, lat=lat, lon=lon, method="nearest").where(time_filter, drop=True).values
    Tmin = ds_tasmin['tasmin'].sel(ens=ens, lat=lat, lon=lon, method="nearest").where(time_filter, drop=True).values
    Tmax = ds_tasmax['tasmax'].sel(ens=ens, lat=lat, lon=lon, method="nearest").where(time_filter, drop=True).values
    RH = ds_hurs['hurs'].sel(ens=ens, lat=lat, lon=lon, method="nearest").where(time_filter, drop=True).values
    Wind = ds_sfcwind['sfcwind'].sel(ens=ens, lat=lat, lon=lon, method="nearest").where(time_filter, drop=True).values
    P = ds_pet['pet'].sel(ens=ens, lat=lat, lon=lon, method="nearest").where(time_filter, drop=True).values
    p = ds_pr['pr'].sel(ens=ens, lat=lat, lon=lon, method="nearest").where(time_filter, drop=True).values

    # Berekende velden
    RHmin = (0.6108 * np.exp((17.27 * Tmin) / (Tmin + 237.3))) / (0.6108 * np.exp((17.27 * Tmax) / (Tmax + 237.3))) * 100
    dRH = (RH - RHmin) / 2
    RHmax = np.where(RH + dRH > 100, 100, np.where(dRH < 0, RHmin, RH + dRH))
    Tsoil = T * 0.949 + 0.6136
    

    # DataFrame bouwen
    df = pd.DataFrame({
        'Date': pd.to_datetime(dates),
        'Wind': Wind,
        'T': T,
        'Tmin': Tmin,
        'Tmax': Tmax, 
        'RH': RH,
        'RHmin': RHmin,
        'RHmax': RHmax,
        'dRH': dRH,
        'P': P,
        'p': p,
        'Tsoil': Tsoil
    })
    return df.dropna()

def watercalculations(df, soil):
    # soil is een dict met keys zoals in SOIL_PARAMS['zand']
    Kcini  = soil["Kcini"]; Kcmid = soil["Kcmid"]; Kclate = soil["Kclate"]
    Zr     = soil["Zr"];    Ofc   = soil["Ofc"];   Owp    = soil["Owp"]
    ini    = soil["ini"];   dev   = soil["dev"];   mid    = soil["mid"]; late = soil["late"]
    pcrop  = soil["pcrop"]

    df = df.copy()
    kolommen = ['Tmin', 'Tmax', 'T', 'RH', 'RHmin', 'RHmax', 'P', 'p', 'Wind']
    df[kolommen] = df[kolommen].apply(pd.to_numeric, errors='coerce')
    df = df.dropna(subset=kolommen).reset_index(drop=True)
    df['day_number'] = range(1, len(df) + 1)

    latitude = 51.99806
    phi = np.radians(latitude)
    helling = 1/3  #1/2, 1/3, 1/4, 1/5
    beta = np.arctan(helling)
    J = df['Date'].dt.dayofyear

    delta = 0.409 * np.sin((2 * np.pi / 365) * J - 1.39)
    omega_s = (np.pi / 2) - np.arctan((-np.tan(phi) * np.tan(delta)) / np.sqrt(np.maximum(1 - (np.tan(phi)**2) * (np.tan(delta)**2), 1e-5)))
   
    G_sc = 0.0820
    dr = 1 + 0.033 * np.cos(2 * np.pi * J / 365)
    Ra = (24 * 60 / np.pi) * G_sc * dr * (omega_s * np.sin(phi) * np.sin(delta) + np.cos(phi) * np.cos(delta) * np.sin(omega_s))

    r_zuid = ( np.cos(phi - beta) * np.cos(delta) * np.sin(omega_s)
    + omega_s * np.sin(phi - beta) * np.sin(delta)) / (
    np.cos(phi) * np.cos(delta) * np.sin(omega_s)
    + omega_s * np.sin(phi) * np.sin(delta))
        
    r_noord = 1
        
    # Ra_slope = Ra * r_zuid
    Ra_slope = Ra * r_noord
    Rso = 0.75 * Ra_slope
    Rs = 0.16 * np.sqrt(df['Tmax'] - df['Tmin']) * Ra_slope

    e0_Tmax = 0.6108 * np.exp((17.27 * df['Tmax']) / (df['Tmax'] + 237.3))
    e0_Tmin = 0.6108 * np.exp((17.27 * df['Tmin']) / (df['Tmin'] + 237.3))
    es = (e0_Tmax + e0_Tmin) / 2

    df['delta'] = delta
    df['omega_s'] = omega_s
    df['es'] = es
    df['Ra'] = Ra
    df['Rso'] = Rso

    ea = (df['RHmin'] * e0_Tmax + df['RHmax'] * e0_Tmin) / 200

    Tmax_K = df['Tmax'] + 273.16
    Tmin_K = df['Tmin'] + 273.16
    Rs_Rso = np.minimum(Rs / Rso, 1.0)
    sigma = 4.903e-9
    Rnl = sigma * ((Tmax_K**4 + Tmin_K**4) / 2) * (0.34 - 0.14 * np.sqrt(ea)) * (1.35 * Rs_Rso - 0.35)
    Rn = 0.77 * Rs - Rnl

    gamma = 0.665e-3 * df['P']
    delta_slope = 4098 * (0.6108 * np.exp((17.27 * df['T']) / (df['T'] + 237.3))) / ((df['T'] + 237.3)**2)
    ET0 = ((0.408 * delta_slope * Rn) + gamma * (900 / (df['T'] + 273)) * (0.75 * df['Wind']) * (es - ea)) / (delta_slope + gamma * (1 + 0.34 * 0.75 * df['Wind']))
    df['ET0'] = ET0

    df['Kc'] = np.nan; df['ETc'] = np.nan; df['Ks'] = np.nan
    df['Diend'] = np.nan; df['Distart'] = np.nan; df['RAW'] = np.nan

    TAW = 1000 * (Ofc - Owp) * Zr
    Distart = 0.0
    
    # TODO (aangepaste formule voor soil_water_start)
    soil_water_start = Ofc * np.sqrt(300)   # startvoorraad (mm)
    
    for i in range(len(df)):
        day = df.at[i, 'day_number']
        T   = df.at[i, 'T']
        RH  = df.at[i, 'RH']
        ET0 = df.at[i, 'ET0']
        p   = float(df.at[i, 'p'])
    
        # Kc per groeifase
        if day <= ini:
            Kc = Kcini
        elif day <= ini + dev and dev > 0:
            Kc = Kcini + ((day - ini) / dev) * (Kcmid - Kcini)
        elif day <= ini + dev + mid:
            Kc = Kcmid + ((0.04 * (T - 2)) - (0.004 * (RH - 45) * ((Zr / 3) ** 0.3)))
        elif late > 0:
            Kc = Kcmid + ((day - ini - dev - mid) / late) * (Kclate - Kcmid)
        else:
            Kc = Kclate
    
        ETc = Kc * ET0
        RAW = (pcrop + 0.04 * (5 - ETc)) * TAW
        denom = (TAW - RAW) if (TAW - RAW) != 0 else 1e-6
        Ks = max(0.0, min(1.0, (TAW - Distart) / denom))
    
        # Runoff
        soil_water_after_p = soil_water_start + p
        runoff = max(
            0.0,
            p - max(
                3.0,
                40.0 * np.exp(-1.5/3.0) * np.sqrt(max(0.0, 1.0 - min(soil_water_after_p, 120.0)/120.0))
            )
        )
    
        p_eff = p - runoff
    
        final_soil_water = soil_water_start + p_eff - (ETc * Ks)
        DP = max(0.0, final_soil_water - soil_water_start)    
        soil_water_next_day = final_soil_water - DP 
        soil_water_next_day = max(0.0, soil_water_next_day)
    
        # Deficit-update met p_eff
        Diend = Distart + (ETc * Ks) - p_eff
    
        df.at[i, 'Kc'] = Kc
        df.at[i, 'ETc'] = ETc
        df.at[i, 'Ks'] = Ks
        df.at[i, 'Diend'] = Diend
        df.at[i, 'Distart'] = Distart
        df.at[i, 'RAW'] = RAW
        df.at[i, 'soil_water_start'] = soil_water_start
        df.at[i, 'soil_water_after_p'] = soil_water_after_p
        df.at[i, 'runoff'] = runoff
        df.at[i, 'p_eff'] = p_eff
        df.at[i, 'final_soil_water'] = final_soil_water
        df.at[i, 'DP'] = DP
        df.at[i, 'soil_water_next_day'] = soil_water_next_day
    
        soil_water_start = soil_water_next_day
        Distart = 0.0 if (pd.isna(Diend) or Diend < 0) else Diend
    
    return df, TAW
    
    for i in range(len(df)):
        day = df.at[i, 'day_number']
        T = df.at[i, 'T']
        RH = df.at[i, 'RH']
        ET0 = df.at[i, 'ET0']
        p = df.at[i, 'p']
    
        # Kc bepalen per groeifase
        if day <= ini:
            Kc = Kcini
        elif day <= ini + dev:
            Kc = Kcini + ((day - ini) / dev) * (Kcmid - Kcini)
        elif day <= ini + dev + mid:
            Kc = Kcmid + ((0.04 * (T - 2)) - (0.004 * (RH - 45) * ((Zr / 3) ** 0.3)))
        else:
            Kc = Kcmid + ((day - ini - dev - mid) / late) * (Kclate - Kcmid)
    
        # ETc berekenen
        ETc = Kc * ET0
    
        # RAW - Readily Available Water
        RAW = (pcrop + 0.04 * (5 - ETc)) * TAW
    
        # Ks berekenen
        Ks = 1 if (TAW - Distart) / (TAW - RAW) > 1 else (TAW - Distart) / (TAW - RAW)

        # Diend berekenen
        Diend = Distart + (ETc * Ks) - p
        
        # Waarden opslaan in de DataFrame
        df.at[i, 'Kc'] = Kc
        df.at[i, 'ETc'] = ETc
        df.at[i, 'Ks'] = Ks
        df.at[i, 'Diend'] = Diend
        df.at[i, 'Distart'] = Distart
        df.at[i, 'RAW'] = RAW
        
        if pd.isna(Diend) or Diend < 0:
            Distart = 0
        else:
            Distart = Diend 
     # Runoff = 

    return df, TAW


def verwelkingspunt(df, TAW):
    grens = 0.8 * TAW
    max_dagen_onder_grens = 0
    huidig_lengte = 0

    onder_grens_dagen = []

    for i, waarde in enumerate(df['Distart']):
        if waarde >= grens:
            huidig_lengte += 1
            onder_grens_dagen.append(df['Date'].iloc[i])
            max_dagen_onder_grens = max(max_dagen_onder_grens, huidig_lengte)
        else:
            huidig_lengte = 0  

    totaal_onder_grens = len(onder_grens_dagen)
    eerste_dag = onder_grens_dagen[0] if onder_grens_dagen else None
    laatste_dag = onder_grens_dagen[-1] if onder_grens_dagen else None

    vegetatie_uitval = max_dagen_onder_grens > 45

    return max_dagen_onder_grens, vegetatie_uitval, totaal_onder_grens, eerste_dag, laatste_dag
    #check voor een reset value, aka als het 44 dagen is geweest en na twee dagen regen weer, reset het?

def living_plant_cover(ds):
    """
    Bereken de levende vegetatiebedekking (in %) op basis van aantal droge dagen.
    
    Parameters:
        ds (int or float): Aantal opeenvolgende droge dagen onder de drempel
        D100 (int): Aantal dagen waarbij 100% verwelkt is (default = 50)
        B (float): Vormparameter van de kromme (default = 2)
    """
    D100 = 45
    B = 2
    if ds < D100:
        return 100 * (1 - (ds / D100) ** B)
    else:
        return 0.0
    
    # cover = living_plant_cover(ds)
    # return cover
    

def plot(df, TAW, grondsoort, dijkvak_id, jaar=None, ens=None):
    verwelkingsgrens = 0.8 * TAW
    plt.figure(figsize=(12, 6))
    plt.plot(df['Date'], df['Distart'], label='Watertekort (Distart)')
    plt.plot(df['Date'], df['RAW'], label='Beschikbaar water (RAW)', linestyle=':')
    plt.axhline(verwelkingsgrens, linestyle='--', label='Verwelkingsgrens (80% TAW)')
    plt.title(f"Watertekort per dag – {dijkvak_id} {em_scen} ({grondsoort})")
    plt.xlabel("Datum"); plt.ylabel("Watertekort (mm)")
    plt.gca().invert_yaxis()
    plt.legend(); plt.grid(True); plt.tight_layout()
    
    #TODO
    #path controleren
    
    path = r'C:\Users\marloes.slokker\Infram BV\Infram Projecten - 000370 Droogtegevoeligheid keringen RWS\Uitvoering\Resultaten'
    fname = (f"{path}\{dijkvak_id}_{grondsoort}_jaar{jaar}_ens{ens}_{em_scen}_extreem.png"
             if jaar is not None and ens is not None
             else f"{path}\{dijkvak_id}_{grondsoort}.png")
    plt.savefig(fname, dpi=200, bbox_inches="tight")
    print(f"💾 Plot opgeslagen: {fname}")
    plt.show()
    plt.close()

def analyse_per_jaar(ens, lat, lon, jaar, grondsoort):
    soil = SOIL_PARAMS[grondsoort]
    df_weather = load_weather_from_nc(jaar=jaar, lat=lat, lon=lon, ens=ens)
    df_water, TAW = watercalculations(df_weather, soil)
    max_dagen, uitval, totaal_dagen, eerste_dag, laatste_dag = verwelkingspunt(df_water, TAW)
    cover = living_plant_cover(max_dagen)

    resultaten = [{
        'Locatie': dijkvak_id,
        'Latitude': lat,
        'Longitude': lon,
        'Jaar': jaar,
        'Ensemble': ens,
        'Grondsoort': grondsoort,
        'Langste uitvalperiode (dagen)': max_dagen,
        'Totaal dagen onder TAW': totaal_dagen,
        'Eerste dag onder TAW': eerste_dag.strftime('%Y-%m-%d') if eerste_dag else None,
        'Laatste dag onder TAW': laatste_dag.strftime('%Y-%m-%d') if laatste_dag else None,
        'Vegetatie-uitval': "Ja" if uitval else "Nee",
        'LivingPlantCover (%)': round(cover, 2)
    }]
    return pd.DataFrame(resultaten)

def select_extreme_cases(resultaten_df, lat, lon):
    """
    Selecteert de heftigste jaren + ensembles o.b.v. laagste vegetatiedekking
    en toont de kernresultaten.
    """
    min_cover = resultaten_df["LivingPlantCover (%)"].min()
    kandidaten = resultaten_df[resultaten_df["LivingPlantCover (%)"] == min_cover]

    max_dagen = kandidaten["Totaal dagen onder TAW"].max()
    zwaarste = kandidaten[kandidaten["Totaal dagen onder TAW"] == max_dagen]

    kolommen = [
        'Ensemble',
        'Jaar',
        'Langste uitvalperiode (dagen)',
        'Totaal dagen onder TAW',
        'Eerste dag onder TAW',
        'Laatste dag onder TAW',
        'Vegetatie-uitval',
        'LivingPlantCover (%)'
    ]
    print("Extreemste gevallen op basis van vegetatie-uitval:\n")
    with pd.option_context('display.max_rows', None, 'display.max_columns', None):
        print(zwaarste[kolommen])

    return zwaarste[kolommen]

def main_per_ensembles(ens=ens_idx, lat=lat, lon=lon, jaar=jaar):
    df_weather = load_weather_from_nc(jaar=jaar)
    df_water, TAW = watercalculations(df_weather)
    max_dagen, uitval, totaal_dagen, eerste_dag, laatste_dag = verwelkingspunt(df_water, TAW)
    cover = living_plant_cover(max_dagen)
   
    print(f"Langste aaneengesloten periode onder TAW: {max_dagen} dagen")
    print(f"Totaal aantal dagen onder TAW: {totaal_dagen}")
    
    if eerste_dag and laatste_dag:
        print(f"Eerste dag onder grens: {eerste_dag.strftime('%Y-%m-%d')}")
        print(f"Laatste dag onder grens: {laatste_dag.strftime('%Y-%m-%d')}")
    else:
        print("🟢 Geen enkele dag onder de TAW-grens")

    if uitval:
        print("Vegetatie-uitval gedetecteerd (>45 dagen aaneengesloten)")
    else:
        print("✅ Geen vegetatie-uitval")
    print(f"Droogte Periode = {max_dagen} dagen → LivingPlantCover = {cover:.2f}%")

    plot(df_water, TAW)
    analyse_per_jaar(ens_idx, lat, lon, jaar)
    return max_dagen, uitval, totaal_dagen, eerste_dag, laatste_dag

def main_all_ensembles(lat, lon, jaren):
    alle_resultaten = []
    resultaten_per_grondsoort = {"zand": [], "klei": []}

    # 1) Runs draaien en samenvatten
    for jaar in jaren:
        print(f"\n🗓️ Start analyse voor jaar {jaar}\n")
        for grond in ["zand", "klei"]:
            resultaten_per_jaar = []
            for ens in range(8):
                resultaten_df = analyse_per_jaar(ens, lat, lon, jaar, grond)
                resultaten_per_jaar.append(resultaten_df)
            jaar_df = pd.concat(resultaten_per_jaar, ignore_index=True)
            resultaten_per_grondsoort[grond].append(jaar_df)
            alle_resultaten.append(jaar_df)

    alle_resultaten_df = pd.concat(alle_resultaten, ignore_index=True)

    # 2) CSV per grondsoort met ALLE jaren
    for grond, frames in resultaten_per_grondsoort.items():
        if frames:
            df_out = pd.concat(frames, ignore_index=True)
            
            # TODO
            # Path controleren
            path = r'C:\Users\marloes.slokker\Infram BV\Infram Projecten - 000370 Droogtegevoeligheid keringen RWS\Uitvoering\Resultaten'
            outfile = f"{path}\{dijkvak_id}_{grond}_{em_scen}_alle_jaren.csv"
            df_out.to_csv(outfile, index=False, encoding="utf-8")
            # print(f"💾 Opgeslagen: {outfile}")

            # ---- Nieuw blok: kansberekening LivingPlantCover == 0 ----
            totaal = len(df_out)
            last_col = df_out.iloc[:, -1]  # laatste kolom
            n_zero = (last_col == 0).sum()
            n_high = ((last_col > 80) & (last_col <= 100)).sum()
            n_medium = ((last_col > 50) & (last_col <= 80)).sum()
            n_low = (last_col <= 50).sum()
            kans_zero = n_zero / totaal if totaal > 0 else 0
            kans_low = n_low / totaal if totaal > 0 else 0
            kans_medium = n_medium / totaal if totaal > 0 else 0
            kans_high = n_high / totaal if totaal > 0 else 0
            print(f"📉 Kans op 0% LivingPlantCover ({grond}): {n_zero}/{totaal} = {kans_zero:.2%}")
            print(f"📉 Kans op 0-50% LivingPlantCover ({grond}): {n_low}/{totaal} = {kans_low:.2%}")
            print(f"📉 Kans op 50-80% LivingPlantCover ({grond}): {n_medium}/{totaal} = {kans_medium:.2%}")
            print(f"📉 Kans op 80-100% LivingPlantCover ({grond}): {n_high}/{totaal} = {kans_high:.2%}")
            # ----------------------------------------------------------

    # 3) Extreemste gevallen per grondsoort vinden en PLOTTEN
    for grond in ["zand", "klei"]:
        print(f"\n🔥 Extreemste gevallen ({grond}):\n")
        df_g = alle_resultaten_df[alle_resultaten_df['Grondsoort'] == grond]
        if df_g.empty:
            print("Geen data voor deze grondsoort.")
            continue

        extreem_df = select_extreme_cases(df_g, lat=lat, lon=lon)  # print + tabel terug
        if extreem_df.empty:
            print("Geen extreme gevallen gevonden.")
            continue

        for _, row in extreem_df.iterrows():
            jaar_ext = int(row["Jaar"])
            ens_ext  = int(row["Ensemble"])

            # Data opnieuw laden & waterbalans herberekenen voor EXACT dit jaar/ens/grond
            df_weather = load_weather_from_nc(ens=ens_ext, lat=lat, lon=lon, jaar=jaar_ext)
            df_water, TAW = watercalculations(df_weather, SOIL_PARAMS[grond])

            # Automatisch opslaan + tonen
            plot(df_water, TAW, grond, dijkvak_id, jaar=jaar_ext, ens=ens_ext)
