# -*- coding: utf-8 -*-
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import warnings
import xarray as xr 

warnings.filterwarnings("ignore")

# Handmatige input
grondsoort = "zand"  # 'zand' of 'klei'
dijkvak_id = "Haven Breskens"

ds_hurs = xr.open_dataset('hurs_Hd_2050_interp (1).nc')
ds_tas = xr.open_dataset('tas_Hd_2050_interp.nc')
ds_tasmax = xr.open_dataset('tasmax_Hd_2050_interp.nc')
ds_tasmin = xr.open_dataset('tasmin_Hd_2050_interp.nc')
ds_pr = xr.open_dataset('pr_Hd_2050_interp.nc')
ds_sfcwind = xr.open_dataset('sfcwind_Hd_2050_interp.nc')
ds_pet = xr.open_dataset('pet_Hd_2050_interp.nc')
ds_rsds = xr.open_dataset('rsds_Hd_2050_interp.nc')

#DikeGrass crop parameters getest door Thomas (voor nu)
#hier kan een if statement: if "zand", then:, else (clay)
Kcini, Kcmid, Kclate = 0.5, 0.8, 0.6
h, Zr = 0.4, 0.3
Ofc, Owp = 0.2, 0.05
ini, dev, mid, late = 30, 0, 365, 0
pcrop = 0.45

#Stel hier de gewenste locatie in (coördinaten of index)
ens_idx = 7  #  realisatie
lat = 3.54933  # kies gewenste lat index
lon = 51.40097  # kies gewenste lon index
jaar = 2052

#Haal tijdreeks op voor 1 punt en 1 realisatie
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
    P = ds_pr['pr'].sel(ens=ens, lat=lat, lon=lon, method="nearest").where(time_filter, drop=True).values

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
        'p': P,
        'Tsoil': Tsoil
    })

    return df.dropna()

    return df.dropna()


def watercalculations(df):
    df = df.copy()
    kolommen = ['Tmin', 'Tmax', 'T', 'RH', 'RHmin', 'RHmax', 'P', 'p', 'Wind']
    df[kolommen] = df[kolommen].apply(pd.to_numeric, errors='coerce')
    df = df.dropna(subset=kolommen)
    df = df.reset_index(drop=True)  # << voeg deze regel toe
    df['day_number'] = range(1, len(df) + 1)

    latitude = 51.99806
    phi = np.radians(latitude)
    # TODO
    # check of radians goed overgezet is
    J = df['Date'].dt.dayofyear
    #Solar declination
    delta = 0.409 * np.sin((2 * np.pi / 365) * J - 1.39)
    tan_phi = np.tan(phi)
    tan_delta = np.tan(delta)
    #Sunset hour angle ws
    omega_s = (np.pi / 2) - np.arctan((-tan_phi * tan_delta) / np.sqrt(np.maximum(1 - (tan_phi**2) * (tan_delta**2), 0.00001)))
    #Extraterrestrial Radiation
    G_sc = 0.0820
    dr = 1 + 0.033 * np.cos(2 * np.pi * J / 365)
    Ra = (24 * 60 / np.pi) * G_sc * dr * (omega_s * np.sin(phi) * np.sin(delta) + np.cos(phi) * np.cos(delta) * np.sin(omega_s))

    #Hier moet nog correctie voor Ra in, maar dit liep niet lekker: 

    # # Correctie voor helling & orientatie
    # orientatie_dijk = 180  # zelf invoeren
    # beta = np.radians(Ra)  # hellingshoek in radialen (Ra hier als placeholder, normaal is beta een vaste hellinghoek)
    # gamma = np.radians(orientatie_dijk)
 
    # cos_theta_s = (
    #     np.sin(delta) * np.sin(phi) * np.cos(beta)
    #     - np.sin(delta) * np.cos(phi) * np.sin(beta) * np.cos(gamma)
    #     + np.cos(delta) * np.cos(phi) * np.cos(beta) * np.cos(omega_s)
    #     + np.cos(delta) * np.sin(phi) * np.sin(beta) * np.cos(gamma) * np.cos(omega_s)
    #     + np.cos(delta) * np.sin(beta) * np.sin(gamma) * np.sin(omega_s)
    # )
 
    # cos_theta_h = np.sin(phi) * np.sin(delta) + np.cos(phi) * np.cos(delta) * np.cos(omega_s)
    # ratio = np.where(cos_theta_h != 0, cos_theta_s / cos_theta_h, 0)
    # Ra_corrected = Ra * np.maximum(0, ratio)
 
    # Rso = 0.75 * Ra_corrected
    # Rs = 0.16 * np.sqrt(df['Tmax'] - df['Tmin']) * Ra_corrected

    # # Rso - Clear Sky Radiation
    Rso = 0.75 * Ra
    # # Rs- Hargreaves Estimate
    Rs = 0.16 * np.sqrt(df['Tmax'] - df['Tmin']) * Ra
    e0_Tmax = 0.6108 * np.exp((17.27 * df['Tmax']) / (df['Tmax'] + 237.3))
    e0_Tmin = 0.6108 * np.exp((17.27 * df['Tmin']) / (df['Tmin'] + 237.3))
    # Es - Mean saturated vapour pressure
    es = (e0_Tmax + e0_Tmin) / 2
    
    df['delta'] = delta
    df['omega_s'] = omega_s
    df['es'] = es
    df['Ra'] = Ra
    df['Rso'] = Rso
    # Actuele dampdruk gebaseerd op RHmin en RHmax
    ea = (df['RHmin'] * e0_Tmax + df['RHmax'] * e0_Tmin) / 200

    # Temperaturen omzetten naar Kelvin
    Tmax_K = df['Tmax'] + 273.16
    Tmin_K = df['Tmin'] + 273.16
    # Begrens Rs/Rso tot maximaal 1.0
    Rs_Rso = np.minimum(Rs / Rso, 1.0)
    sigma = 4.903e-9
   # Netto langgolvige straling
    Rnl = sigma * ((Tmax_K**4 + Tmin_K**4) / 2) * (0.34 - 0.14 * np.sqrt(ea)) * (1.35 * Rs_Rso - 0.35)
    #Rn - Net Radiation
    Rn = 0.77 * Rs - Rnl

    # Psychrometrische constante
    gamma = 0.665* 10**-3 * df['P']
    delta_slope = 4098 * (0.6108 * np.exp((17.27 * df['T']) / (df['T'] + 237.3))) / ((df['T'] + 237.3)**2)
    # Penman-Monteith ET0
    ET0 = ((0.408 * delta_slope * Rn) + gamma * (900 / (df['T'] + 273)) * (0.75 * df['Wind']) * (es - ea)) / (delta_slope + gamma * (1 + 0.34 * 0.75 * df['Wind']))
    df['ET0'] = ET0

    #Kc, weet niet of dit compleet goed loopt, dubbele check nodig
    df['Kc'] = np.nan
    df['ETc'] = np.nan
    df['Ks'] = np.nan
    df['Diend'] = np.nan
    df['Distart'] = np.nan
    df['RAW'] = np.nan

    # TAW - Total Available Water
    TAW = 1000 * (Ofc - Owp) * Zr
    Distart = 0  # Startwaarde
    
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
    #kijk of de verwelkingslijn hetzelfde is als in het excel model
    #voeg overige formules toe (runoff, correctie helling)


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
            huidig_lengte = 0  # reset bij boven grens

    # Totaal aantal dagen onder de grens (ongeacht aaneengeslotenheid)
    totaal_onder_grens = len(onder_grens_dagen)

    # Eerste en laatste dag onder grens
    eerste_dag = onder_grens_dagen[0] if onder_grens_dagen else None
    laatste_dag = onder_grens_dagen[-1] if onder_grens_dagen else None

    # Vegetatie-uitval als de langste aaneengesloten periode > 45 dagen
    vegetatie_uitval = max_dagen_onder_grens > 45

    return max_dagen_onder_grens, vegetatie_uitval, totaal_onder_grens, eerste_dag, laatste_dag
    #check voor een reset value, aka als het 44 dagen is geweest en na twee dagen regen weer, reset het?


def plot(df, TAW):
    #TAW = df['TAW'].iloc[0]
    verwelkingsgrens = 0.8 * TAW

    plt.figure(figsize=(12, 6))
    plt.plot(df['Date'], df['Distart'], label='Watertekort (Distart)', color='orange')
    # Beschikbaar water (RAW)
    plt.plot(df['Date'], df['RAW'], label='Beschikbaar water (RAW)', color='blue', linestyle=':')
    plt.axhline(verwelkingsgrens, linestyle='--', color='red', label='Verwelkingsgrens (80% TAW)')
    plt.title(f"Watertekort per dag – {dijkvak_id} ({grondsoort})")
    plt.xlabel("Datum")
    plt.ylabel("Watertekort (mm)")
    plt.gca().invert_yaxis()  # omgekeerde y-as
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

# def print_result_table(df):
#     """
#     Print een tabel met de eerste 10 dagen van de berekening.
#     Alles is nu correct gecheckt
#     """
#     kolommen = ['Date', 'ETc', 'ET0', 'Diend', 'RAW']
#     # Print de eerste 10 rijen als tabel
#     print("\n📊 Eerste 10 dagen van de waterbalansberekening:\n")
#     print(df[kolommen].head(10).to_string(index=False, justify='center'))

def analyse_per_jaar(ens, lat, lon, jaren):
    resultaten = []

    # for jaar in jaren:
    df_weather = load_weather_from_nc(jaar=jaar, lat=lat, lon=lon, ens=ens)
    df_water, TAW = watercalculations(df_weather)
    max_dagen, uitval, totaal_dagen, eerste_dag, laatste_dag = verwelkingspunt(df_water, TAW)

    resultaten.append({
        'Ensemble': ens,
        'Jaar': jaar,
        'Langste uitvalperiode (dagen)': max_dagen,
        'Totaal dagen onder TAW': totaal_dagen,
        'Eerste dag onder TAW': eerste_dag.strftime('%Y-%m-%d') if eerste_dag else None,
        'Laatste dag onder TAW': laatste_dag.strftime('%Y-%m-%d') if laatste_dag else None,
        'Vegetatie-uitval': "Ja" if uitval else "Nee"
    })

    return pd.DataFrame(resultaten)

def main(ens=ens_idx, lat=lat, lon=lon, jaar=jaar):
    df_weather = load_weather_from_nc(jaar=jaar)
    df_water, TAW = watercalculations(df_weather)

    # Haal alle relevante info uit verwelkingsanalyse
    max_dagen, uitval, totaal_dagen, eerste_dag, laatste_dag = verwelkingspunt(df_water, TAW)

    # Print resultaten
    print(f"Langste aaneengesloten periode onder TAW: {max_dagen} dagen")
    print(f"Totaal aantal dagen onder TAW: {totaal_dagen}")
    
    if eerste_dag and laatste_dag:
        print(f"Eerste dag onder grens: {eerste_dag.strftime('%Y-%m-%d')}")
        print(f"Laatste dag onder grens: {laatste_dag.strftime('%Y-%m-%d')}")
    else:
        print("🟢 Geen enkele dag onder de TAW-grens")

    if uitval:
        print("\U0001F33F Vegetatie-uitval gedetecteerd (>45 dagen aaneengesloten)")
    else:
        print("✅ Geen vegetatie-uitval")

    # Plot grafiek
    plot(df_water, TAW)
    analyse_per_jaar(ens_idx, lat, lon, jaar)
    return max_dagen, uitval, totaal_dagen, eerste_dag, laatste_dag

main(ens=ens_idx, lat=lat, lon=lon, jaar=jaar)
