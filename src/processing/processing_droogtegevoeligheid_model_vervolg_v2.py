
# -*- coding: utf-8 -*-

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import warnings
from functools import lru_cache
from src.parsing.parsing_droogtegevoeligheid_model_vervolg_v2  import read_parameters 


warnings.filterwarnings("ignore")

# # Handmatige input
ens_idx = 1  #  realisatie

#DikeGrass crop parameters getest door Thomas (voor nu)
#hier kan een if statement: if "zand", then:, else (clay)

# Gemeenschappelijke grasparameters (ARVIOR-dijkgras). De bodem beïnvloedt
# actuele verdamping via TAW/RAW/Ks, niet via een andere potentiële Kc.
GRASS_PARAMS = {
    "Kcini": 0.40, "Kcmid": 0.95, "Kclate": 0.85,
    "h": 0.25, "ini": 10, "dev": 20, "mid": 300, "late": 35,
    "pcrop": 0.60,
}
SOIL_PARAMS = {
    "zand": dict(GRASS_PARAMS, Zr=0.30, Ofc=0.12, Owp=0.05),
    "klei": dict(GRASS_PARAMS, Zr=0.40, Ofc=0.31, Owp=0.11),
}


# Onafhankelijke onderzoeksopties; hydrologische correcties zijn altijd actief.
from dataclasses import dataclass
from datetime import date

@dataclass(frozen=True)
class ModelOptions:
    use_revised_soil: bool = True
    enable_mowing: bool = False
    enable_ks_mortality: bool = False
    air_pressure_kpa: float = 101.3
    mow_dates_sand: tuple = ((5, 31), (8, 31))
    mow_dates_clay: tuple = ((5, 31), (8, 31))
    mow_kc_fraction: float = 0.4473684210526316
    # Fit op dagelijkse ARVIOR Kcb_target (Vecht/Heumen 2026) + nominale Ke
    # 0.15, genormaliseerd door Kc_mid=0.95, conditioneel op ARVIOR Ks.
    # Dit is een proxyfit, geen onafhankelijke validatie van totale Kc.
    mow_recovery_rate: float = 0.005625
    s50: float = 68.16
    width: float = 1.5
    near_zero_pct: float = 1.0

    def __post_init__(self):
        if not np.isfinite(self.air_pressure_kpa) or not 80 <= self.air_pressure_kpa <= 110:
            raise ValueError("Luchtdruk moet in kPa liggen (80–110)")
        if not (0 < self.mow_kc_fraction <= 1 and 0 <= self.mow_recovery_rate <= 1):
            raise ValueError("Ongeldige maai- of herstelparameter")
        if not (np.isfinite(self.s50) and self.s50 > 0 and np.isfinite(self.width)
                and self.width > 0 and 0 < self.near_zero_pct < 50):
            raise ValueError("Ongeldige mortaliteitsparameters")
        for dates in (self.mow_dates_sand, self.mow_dates_clay):
            for month, day in dates:
                date(2001, month, day)

    def label(self):
        return (f"bodem_{'nieuw' if self.use_revised_soil else 'historisch'}"
                f"_maaien_{int(self.enable_mowing)}_ksmortaliteit_{int(self.enable_ks_mortality)}"
                "_waterbalans_v4fix")

    def mowing_dates(self, grondsoort):
        return self.mow_dates_sand if grondsoort == "zand" else self.mow_dates_clay

DEFAULT_OPTIONS = ModelOptions()

def get_soil_params(grondsoort, options=DEFAULT_OPTIONS):
    soil = dict(SOIL_PARAMS[grondsoort])
    # Historische variant afgeleid uit de opmerkingen in de ontvangen code;
    # niet onafhankelijk geverifieerd als eerdere productieversie.
    if not options.use_revised_soil:
        if grondsoort == "zand":
            soil.update(Ofc=0.20, Owp=0.05)
        elif grondsoort == "klei":
            soil.update(Ofc=0.35, Owp=0.17, Zr=0.30)
    return soil


# KNMI23 daggemiddelden. Bij ontbrekende units geldt de expliciete
# standaard hieronder; onbekende eenheden worden niet stil geaccepteerd.
# https://www.klimaatscenarios-data.knmi.nl/downloads
# Voor rsds: W/m2 daggemiddelde -> MJ/m2/dag (*0.0864).
import re
import unicodedata

def _unit_key(unit):
    unit = unicodedata.normalize("NFKC", str(unit)).lower().strip()
    return re.sub(r"[\s^*()·]", "", unit.replace("−", "-").replace("°", "deg"))


def _convert_knmi(values, name, units=None):
    default = {"tas":"degC", "tasmin":"degC", "tasmax":"degC",
               "hurs":"%", "sfcwind":"m/s", "pr":"mm/day",
               "pet":"mm/day", "rsds":"W/m2"}
    u = _unit_key(units or default[name])
    x = np.asarray(values, dtype=float)
    if name in ("tas", "tasmin", "tasmax"):
        if u in {"k", "kelvin"}: return x - 273.15
        if u in {"degc", "c", "celsius", "degreecelsius", "degreescelsius", "degree_celsius", "degrees_celsius", "deg_c", "degrees_c", "degree_c"}: return x
    elif name == "hurs":
        if u in {"%", "percent", "percentage", "pct", "0-100"}: return x
        if u in {"1", "fraction", "0-1", "dimensionless"}: return x * 100
    elif name == "sfcwind":
        if u in {"m/s", "ms-1", "m/sec", "msecond-1", "meterpersecond"}: return x
        if u in {"km/h", "kmh-1"}: return x / 3.6
    elif name in ("pr", "pet"):
        if u in {"mm", "mm/day", "mmd-1", "mmday-1", "mm/d", "mm/24h", "kgm-2day-1", "kgm-2d-1", "kg/m2/day", "kgm-2"}: return x
        if u in {"kgm-2s-1", "kg/m2/s", "mms-1", "mm/s"}: return x * 86400
        if u in {"m/day", "md-1", "mday-1", "m/d", "m"}: return x * 1000
    elif name == "rsds":
        if u in {"w/m2", "wm-2", "watt/m2", "wattm-2", "watts/m2"}: return x * 0.0864
        if u in {"mj/m2/day", "mjm-2day-1", "mjm-2d-1", "mj/m2/d", "mj/m2", "mjm-2"}: return x
        if u in {"j/m2/day", "jm-2day-1", "jm-2d-1", "j/m2/d", "j/m2", "jm-2"}: return x / 1e6
        if u in {"j/cm2", "jcm-2", "j/cm2/day", "jcm-2day-1", "jcm-2d-1"}: return x * 0.01
    raise ValueError(f"Onbekende KNMI-eenheid voor {name}: {units!r}")


def _saturation_vapour_pressure(t):
    return 0.6108 * np.exp(17.27 * t / (t + 237.3))


def load_weather_from_nc(path, em_scen, ens, lat, lon, jaar):
    datasets = read_parameters(path, em_scen)
    names = ("hurs", "tas", "tasmax", "tasmin", "pr", "sfcwind", "pet", "rsds")
    data = {}
    dates = None
    for name, ds in zip(names, datasets):
        da = ds[name].sel(ens=ens, lat=lat, lon=lon, method="nearest")
        da = da.where(da['time'].dt.year == int(jaar), drop=True)
        d = pd.to_datetime(da['time'].values)
        if dates is None:
            dates = d
        elif not dates.equals(d):
            raise ValueError(f"Afwijkende tijdas in KNMI-bestand {name}")
        data[name] = _convert_knmi(da.values, name, da.attrs.get("units"))
    T, Tmin, Tmax = data['tas'], data['tasmin'], data['tasmax']
    RH = data['hurs']
    # FAO-56: ea bij uitsluitend gemiddelde RH. RHmin en RHmax zijn
    # benaderingen bij gelijkblijvende dampdruk, alleen voor Kc-correctie.
    es = (_saturation_vapour_pressure(Tmin) + _saturation_vapour_pressure(Tmax)) / 2
    ea = es * RH / 100
    RHmin = np.clip(100 * ea / _saturation_vapour_pressure(Tmax), 0, 100)
    RHmax = np.clip(100 * ea / _saturation_vapour_pressure(Tmin), 0, 100)
    df = pd.DataFrame({
        'Date': dates, 'Wind': data['sfcwind'], 'T': T, 'Tmin': Tmin, 'Tmax': Tmax,
        'RH': RH, 'RHmin': RHmin, 'RHmax': RHmax, 'dRH': (RHmax-RHmin)/2,
        'P': data['pet'], 'p': data['pr'], 'Tsoil': T * 0.949 + 0.6136,
        'Rs_horizontal': data['rsds'],  # MJ/m2/dag; horizontaal oppervlak
    })
    return df

def watercalculations(df, soil, lat, helling, options=DEFAULT_OPTIONS, mowing_dates=()):
    # soil is een dict met keys zoals in SOIL_PARAMS['zand']
    Kcini  = soil["Kcini"]; Kcmid = soil["Kcmid"]; Kclate = soil["Kclate"]
    Zr     = soil["Zr"];    Ofc   = soil["Ofc"];   Owp    = soil["Owp"]
    ini    = soil["ini"];   dev   = soil["dev"];   mid    = soil["mid"]; late = soil["late"]
    pcrop  = soil["pcrop"]

    df = df.copy()

    kolommen = ['Tmin', 'Tmax', 'T', 'RH', 'p', 'Wind', 'Rs_horizontal']
    df[kolommen] = df[kolommen].apply(pd.to_numeric, errors='coerce')
    
    
    ####################################################################
    # NIEUW STUKJE
    ####################################################################
    # Controle op ontbrekende / niet-numerieke waarden
    if not np.isfinite(df[kolommen]).all().all():
        raise ValueError('Ontbrekende of ongeldige KNMI-dagwaarden')
    
    
    # ---------------------------------------------------------
    # Correctie voor dagen waarop Tmin > Tmax in de KNMI-data
    # ---------------------------------------------------------
    mask_temp = df['Tmin'] > df['Tmax']
    
    if mask_temp.any():
        print(
            f"Waarschuwing: {mask_temp.sum()} dagen met Tmin > Tmax. "
            "Tmin en Tmax worden voor deze dagen omgewisseld."
        )
    
        # Bewaar oorspronkelijke Tmin tijdelijk
        tmin_origineel = df.loc[mask_temp, 'Tmin'].copy()
    
        # Verwissel Tmin en Tmax
        df.loc[mask_temp, 'Tmin'] = df.loc[mask_temp, 'Tmax']
        df.loc[mask_temp, 'Tmax'] = tmin_origineel
    
    
    # Controle overige klimaatvariabelen
    if (
       # not df['RH'].between(0, 100).all()
        not (df[['Wind', 'p', 'Rs_horizontal']] < 0).any().any()
    ):
        raise ValueError(
            'Ongeldige klimaatinvoer; controleer RH, Wind, neerslag en straling'
        )
    
    
    df = df.reset_index(drop=True)
    df['day_number'] = df['Date'].dt.dayofyear

    ####################################################################
    # EINDE NIEUW STUKJE (ter vervanging van onderstaande)
    ####################################################################

    ############ LET OP ONDERSTAANDE IS VERVANGEN DOOR BOVENSTAANDE IVM EEN AANTAL WAARDEN DIE MAX WAREN AANGEGEVEN WAAR MINWAARDE WAS.
    # df = df.copy()
    # kolommen = ['Tmin', 'Tmax', 'T', 'RH', 'p', 'Wind', 'Rs_horizontal']
    # df[kolommen] = df[kolommen].apply(pd.to_numeric, errors='coerce')
    # if not np.isfinite(df[kolommen]).all().all():
    #     raise ValueError('Ontbrekende of ongeldige KNMI-dagwaarden')
    
    # if (df['Tmin'] > df['Tmax']).any() or not df['RH'].between(0, 100).all() or (df[['Wind','p','Rs_horizontal']] < 0).any().any():
    #     raise ValueError('Ongeldige klimaatinvoer; controleer eenheden')
    # df = df.reset_index(drop=True)
    # df['day_number'] = df['Date'].dt.dayofyear
    ############ einde oude stukje
    
    #TODO
    # Latitude moet een range worden.
    
    latitude = lat
    phi = np.radians(latitude)
    helling = helling  #1/2, 1/3, 1/4, 1/5
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
    
    # Orientatie
    ORIENTATIE = "noord"   # "zuid" of "noord"
    Ra_slope = Ra * (r_zuid if ORIENTATIE == "zuid" else r_noord)
    
    Rso = 0.75 * Ra_slope
    # Gebruik gemeten/gemodelleerde horizontale KNMI-straling; behoud de
    # bestaande zuid/noord-taludfactor (geen nieuwe oriëntatiemodule).
    Rs = df['Rs_horizontal'] * (r_zuid if ORIENTATIE == "zuid" else r_noord)
    df['Rs_slope'] = Rs

    e0_Tmax = 0.6108 * np.exp((17.27 * df['Tmax']) / (df['Tmax'] + 237.3))
    e0_Tmin = 0.6108 * np.exp((17.27 * df['Tmin']) / (df['Tmin'] + 237.3))
    es = (e0_Tmax + e0_Tmin) / 2

    df['delta'] = delta
    df['omega_s'] = omega_s
    df['es'] = es
    df['Ra'] = Ra
    df['Rso'] = Rso

    ea = es * df['RH'] / 100.0
    df['RHmin'] = (100 * ea / e0_Tmax).clip(0, 100)
    df['RHmax'] = (100 * ea / e0_Tmin).clip(0, 100)

    Tmax_K = df['Tmax'] + 273.16
    Tmin_K = df['Tmin'] + 273.16
    Rs_Rso = (Rs / Rso).clip(0.3, 1.0)
    sigma = 4.903e-9
    Rnl = sigma * ((Tmax_K**4 + Tmin_K**4) / 2) * (0.34 - 0.14 * np.sqrt(ea)) * (1.35 * Rs_Rso - 0.35)
    Rn = 0.77 * Rs - Rnl

    # PET is geen luchtdruk; gebruik de ingestelde atmosferische druk in kPa.
    gamma = 0.665e-3 * options.air_pressure_kpa
    delta_slope = 4098 * (0.6108 * np.exp((17.27 * df['T']) / (df['T'] + 237.3))) / ((df['T'] + 237.3)**2)
    
    # TODO
    # factor windsnelheid (wind, van 10m mean naar 2m mean)
    
    factor_10m_to_2m = 4.87 / np.log(67.8 * 10 - 5.42)
    u2 = factor_10m_to_2m * df['Wind']
    ET0 = ((0.408 * delta_slope * Rn) + gamma * (900 / (df['T'] + 273)) * u2 * (es - ea)) / (delta_slope + gamma * (1 + 0.34 * u2))
    ET0 = ET0.clip(lower=0.0)
    df['ET0'] = ET0
    
    
    ####################### nieuw ############################
    TAW = 1000 * (Ofc - Owp) * Zr
    Distart = 0.0

    # TODO (aangepaste formule voor soil_water_start)
    soil_water_start = TAW
    # Correctie: voorraad is plantbeschikbaar water boven WP (mm); start op FC.
    if not np.isfinite(TAW) or TAW <= 0:
        raise ValueError("TAW moet positief zijn")

    # Sneller: de kolommen één keer omzetten naar numpy-arrays (in plaats van df.at in de loop)
    runoff_factor = 40.0 * np.exp(-1.5 * helling)   # constant, dus buiten de loop
    day_arr = df['day_number'].to_numpy()
    u2_arr = u2.to_numpy()
    RHmin_arr = df['RHmin'].to_numpy()
    ET0_arr = df['ET0'].to_numpy()
    p_arr   = df['p'].to_numpy(dtype=float)

    namen = ['Kc', 'ETc', 'Ks', 'Diend', 'Distart', 'RAW', 'soil_water_start',
             'soil_water_after_p', 'runoff', 'p_eff', 'final_soil_water', 'DP',
             'soil_water_next_day']
    uit = {naam: np.full(len(df), np.nan) for naam in namen}
    for naam in ('ETa', 'VWC_rootzone_pct', 'water_balance_residual_mm'):
        uit[naam] = np.full(len(df), np.nan)
    if df.empty or not df['Date'].diff().dropna().eq(pd.Timedelta(days=1)).all():
        raise ValueError('Gecorrigeerde waterbalans vereist aaneengesloten dagelijkse weerdata')
    if not np.isfinite(df[['ET0','p']]).all().all() or (df['p'] < 0).any():
        raise ValueError('Ongeldige ET0 of neerslag')

    # Kc is hier een enkele gewascoefficient (geen gescheiden Kcb/Ke).
    # De factor is een voorlopige reductie van totale Kc, geen biomassamodel.
    if options.enable_mowing:
        dates_arr = list(zip(df['Date'].dt.month, df['Date'].dt.day))
        mow_dates = set(mowing_dates)
        canopy_factor = 1.0
        uit['Kc_unmown'] = np.full(len(df), np.nan)
        uit['mowing_factor'] = np.full(len(df), np.nan)
        uit['mowing_today'] = np.zeros(len(df), dtype=bool)

    for i in range(len(df)):
        day = day_arr[i]
        ET0 = ET0_arr[i]
        p   = p_arr[i]

        # Maai-event vóór gewasverdamping; Kc herstelt geleidelijk.
        if options.enable_mowing:
            cut = dates_arr[i] in mow_dates
            if cut:
                canopy_factor = min(canopy_factor, options.mow_kc_fraction)
        # Kc per groeifase
        if day <= ini:
            Kc = Kcini
        elif day <= ini + dev and dev > 0:
            Kc = Kcini + ((day - ini) / dev) * (Kcmid - Kcini)
        elif day <= ini + dev + mid:
            Kc = Kcmid
        elif late > 0:
            Kc = Kcmid + min(1.0, (day - ini - dev - mid) / late) * (Kclate - Kcmid)
        else:
            Kc = Kclate

        # FAO-56 klimaatcorrectie op midden- en eindfase; lineair ingroeien
        # tijdens de ontwikkelingsfase. Dezelfde parameters voor beide bodems.
        climate = (0.04 * (np.clip(u2_arr[i], 1.0, 6.0) - 2)
                   - 0.004 * (np.clip(RHmin_arr[i], 20.0, 80.0) - 45))
        phase_weight = np.clip((day - ini) / max(1, dev), 0, 1)
        full_height = soil['h']
        kc_base = Kc
        Kc = max(0.0, kc_base + phase_weight * climate * (full_height / 3) ** 0.3)
        if options.enable_mowing:
            uit['Kc_unmown'][i] = Kc
            uit['mowing_factor'][i] = canopy_factor
            uit['mowing_today'][i] = cut
            # Benadering voor enkelvoudige Kc, geen biomassamodel. De factor
            # volgt ARVIOR Kcb_target plus nominale Ke=0.15, geen gemeten Ke.
            # 0.275 + 0.15 bij 10 cm, 0.80 + 0.15 bij volgroeid gras.
            Kc *= canopy_factor
        ETc = max(0.0, Kc * ET0)
        RAW = np.clip(pcrop + 0.04 * (5 - ETc), 0.1, 0.8) * TAW
        # Runoff
        soil_water_after_p = soil_water_start + p
        runoff = max(
            0.0,
            p - max(
                3.0,
                runoff_factor * np.sqrt(max(0.0, 1.0 - min(soil_water_after_p, 120.0)/120.0))
            )
        )

        p_eff = p - runoff

        # Neerslag vult eerst de wortelzone aan; Ks na infiltratie,
        # vóór ET. Dagoverschot boven FC draineert na ET.
        available = soil_water_start + p_eff
        d_before_et = max(0.0, TAW - available)
        Ks = float(np.clip((TAW - d_before_et) / (TAW - RAW), 0.0, 1.0))
        eta = min(available, ETc * Ks)
        final_soil_water = available - eta
        DP = max(0.0, final_soil_water - TAW)
        soil_water_next_day = final_soil_water - DP
        Diend = TAW - soil_water_next_day
        uit['ETa'][i] = eta
        uit['VWC_rootzone_pct'][i] = 100 * (Owp + soil_water_next_day / (1000 * Zr))
        uit['water_balance_residual_mm'][i] = soil_water_start + p - runoff - eta - DP - soil_water_next_day

        # Resultaten in de arrays zetten (in plaats van df.at)
        uit['Kc'][i] = Kc
        uit['ETc'][i] = ETc
        uit['Ks'][i] = Ks
        uit['Diend'][i] = Diend
        uit['Distart'][i] = Distart
        uit['RAW'][i] = RAW
        uit['soil_water_start'][i] = soil_water_start
        uit['soil_water_after_p'][i] = soil_water_after_p
        uit['runoff'][i] = runoff
        uit['p_eff'][i] = p_eff
        uit['final_soil_water'][i] = final_soil_water
        uit['DP'][i] = DP
        uit['soil_water_next_day'][i] = soil_water_next_day

        soil_water_start = soil_water_next_day
        Distart = 0.0 if (pd.isna(Diend) or Diend < 0) else Diend
        if options.enable_mowing:
            # herstel van verdampingscapaciteit, uitsluitend voor volgende dag
            canopy_factor = min(1.0, canopy_factor + options.mow_recovery_rate * Ks * (1.0 - canopy_factor))

    # Na de loop alles in één keer terugzetten in de DataFrame
    for naam in uit:
        df[naam] = uit[naam]

    return df, TAW


def verwelkingspunt(df, TAW, D100=45, B=2):
    """
    Bereken kenmerken van verwelkingspunt en levende vegetatiebedekking.
    
    Parameters:
        df (pd.DataFrame): DataFrame met kolommen 'Date' en 'Distart'
        TAW (float): Total Available Water (gebruikt voor de grens)
        D100 (int): Aantal dagen waarbij 100% verwelkt is (default = 45)
        B (float): Vormparameter van de kromme (default = 2)

    Returns:
        df (pd.DataFrame): DataFrame met extra kolommen:
                           - 'dagen_boven_grens'
                           - 'living_plant_cover'
        dict: Dictionary met kernresultaten:
              - max_dagen_boven_grens
              - vegetatie_uitval (bool)
              - totaal_boven_grens
              - eerste_dag_boven_grens
              - laatste_dag_boven_grens
              - maatgevend_cover
              - dag_maatgevend_cover
              - cover_1_april
              - cover_1_oktober
    """
    grens = 0.8 * TAW 
    max_dagen_boven_grens = 0
    huidig_lengte = 0
    
    boven_grens_dagen = []
    reeks_per_dag = []
    cover_per_dag = []
    cover_0_bereikt = False
    cover_below_50_reached = False
    last_cover = 100.0  # startwaarde (kan nooit hoger worden later)
    
    for i, waarde in enumerate(df['Distart']):
        # Bepaal dagen boven grens
        if waarde >= grens:
            huidig_lengte += 1
            boven_grens_dagen.append(df['Date'].iloc[i])
            max_dagen_boven_grens = max(max_dagen_boven_grens, huidig_lengte)
        else:
            huidig_lengte = 0
        reeks_per_dag.append(huidig_lengte)
    
        # Bereken living plant cover
        ds = huidig_lengte
        if cover_0_bereikt:
            cover = 0.0
        else:
            if ds < D100:
                cover = 100 * (1 - (ds / D100) ** B)
            else:
                cover = 0.0
                cover_0_bereikt = True

        # Als cover < 50 bereikt is, nooit meer omhoog laten gaan
        if cover_below_50_reached:
            cover = min(cover, last_cover)  # max blijft het vorige getal
        if cover < 50:
            cover_below_50_reached = True
    
        last_cover = cover
        cover_per_dag.append(cover)


    # Zet kolommen in dataframe
    df['dagen_boven_grens'] = reeks_per_dag
    df['living_plant_cover'] = cover_per_dag

    # Extra resultaten
    totaal_boven_grens = len(boven_grens_dagen)
    eerste_dag = boven_grens_dagen[0] if boven_grens_dagen else None
    laatste_dag = boven_grens_dagen[-1] if boven_grens_dagen else None
    vegetatie_uitval = max_dagen_boven_grens > 45

    # Maatgevende cover en dag
    maatgevend_cover = min(cover_per_dag)
    dag_maatgevend_cover = df['Date'].iloc[cover_per_dag.index(maatgevend_cover)]

    # Cover op 1 oktober (None als die datum niet voorkomt)
    cover_1_oktober = None
    if not df[df['Date'].dt.month.eq(10) & df['Date'].dt.day.eq(1)].empty:
        cover_1_oktober = df.loc[(df['Date'].dt.month == 10) & (df['Date'].dt.day == 1), 'living_plant_cover'].iloc[0]

    resultaten = {
        "max_dagen_boven_grens": max_dagen_boven_grens,
        "vegetatie_uitval": vegetatie_uitval,
        "totaal_boven_grens": totaal_boven_grens,
        "eerste_dag_boven_grens": eerste_dag,
        "laatste_dag_boven_grens": laatste_dag,
        "maatgevend_cover": maatgevend_cover,
        "dag_maatgevend_cover": dag_maatgevend_cover,
        "cover_1_oktober": cover_1_oktober,
        "cover_per_dag": cover_per_dag,
    }

    return df, resultaten
    


def vegetatierespons(df, TAW, options=DEFAULT_OPTIONS):
    """Selecteer oude TAW-functie of cumulatieve Ks-functie.

    TAW-diagnostiek blijft beschikbaar onder beide varianten; deze stuurt
    in de nieuwe variant de mortaliteit niet. Jaarreset blijft gehandhaafd.
    """
    df, resultaten = verwelkingspunt(df, TAW, D100=45, B=2)
    if not options.enable_ks_mortality:
        return df, resultaten
    if df.empty:
        raise ValueError("Geen dagwaarden voor mortaliteit")
    ks = df["Ks"].to_numpy(dtype=float)
    dates = pd.to_datetime(df["Date"])
    if not np.isfinite(ks).all():
        raise ValueError("Niet-eindige Ks in mortaliteitsreeks")
    if not dates.diff().dropna().eq(pd.Timedelta(days=1)).all():
        raise ValueError("Ontbrekende of dubbele dagen in Ks-reeks")
    # De jaarlijkse runner start op 1 januari; geen herstel van stressgeheugen.
    s = 1.0 - np.clip(ks, 0.0, 1.0)
    S = np.cumsum(s)
    exponent = np.clip((S - options.s50) / options.width, -700, 700)
    cover = 100.0 / (1.0 + np.exp(exponent))
    df["stressdag_equivalent"] = s
    df["stress_equivalent_cumulatief"] = S
    df["mortaliteit_pct"] = 100.0 - cover
    df["living_plant_cover"] = cover
    j = int(np.argmin(cover))
    okt = (dates.dt.month == 10) & (dates.dt.day == 1)
    resultaten.update(
        maatgevend_cover=float(cover[j]),
        dag_maatgevend_cover=df["Date"].iloc[j],
        vegetatie_uitval=bool(round(float(cover[j]), 2) <= options.near_zero_pct),
        cover_1_oktober=float(cover[np.flatnonzero(okt)[0]]) if okt.any() else None,
        cover_per_dag=cover.tolist(),
        cumulatieve_stress=float(S[-1]),
        max_mortaliteit=float(100.0-cover[j]),
    )
    return df, resultaten

def run_dagmodel(df_weather, grondsoort, lat, helling, options=DEFAULT_OPTIONS):
    """Eén centraal pad voor alle jaarlijkse runs en extreme gevallen."""
    soil = get_soil_params(grondsoort, options)
    df, TAW = watercalculations(df_weather, soil, lat, helling,
                                options=options, mowing_dates=options.mowing_dates(grondsoort))
    df, resultaten = vegetatierespons(df, TAW, options)
    return df, TAW, resultaten


def plot(df, TAW, grondsoort, dijkvak_id, em_scen, cot, year, jaar=None, ens=None):
    verwelkingsgrens = 0.8 * TAW
    plt.figure(figsize=(12, 6))
    plt.plot(df['Date'], df['Distart'], label='Watertekort (Distart)')
    plt.plot(df['Date'], df['RAW'], label='Beschikbaar water (RAW)', linestyle=':')
    plt.axhline(verwelkingsgrens, linestyle='--', label='Verwelkingsgrens (80% TAW)')
    plt.title(f"Watertekort per dag – {em_scen} ({grondsoort}), helling = {cot}, ens = {ens}, dijk_id = {dijkvak_id}")
    plt.xlabel("Datum"); plt.ylabel("Watertekort (mm)")
    plt.gca().invert_yaxis()
    plt.legend(); plt.grid(True); plt.tight_layout()
    
    #TODO
    #path controleren
    
    output_path = r'C:\Users\ingeb\Documents\Infram\Droogtegevoeligheid\Scenario Data'
    fname = (f"{output_path}\{year}\{em_scen}_{grondsoort}_{cot}_{jaar}_ens{ens}_dijkid{dijkvak_id}_extreem.png"
             if jaar is not None and ens is not None
             else f"{output_path}\{dijkvak_id}_{grondsoort}.png")
    plt.savefig(fname, dpi=200, bbox_inches="tight")
    print(f"💾 Plot opgeslagen: {em_scen}_{grondsoort}_{cot}_{jaar}_ens{ens}_dijkid{dijkvak_id}_extreem.png")
    plt.show()
    plt.close()

def analyse_per_jaar(path, em_scen, ens, lat, lon, jaar, grondsoort, dijkvak_id, helling, options=DEFAULT_OPTIONS):
    df_weather = load_weather_from_nc(path=path, em_scen=em_scen, ens=ens, lat=lat, lon=lon, jaar=jaar)
    df_water, TAW, resultaten_verwelkingspunt = run_dagmodel(df_weather, grondsoort, lat, helling, options)
    max_dagen = resultaten_verwelkingspunt['max_dagen_boven_grens']
    totaal_dagen = resultaten_verwelkingspunt['totaal_boven_grens']
    eerste_dag = resultaten_verwelkingspunt['eerste_dag_boven_grens']
    laatste_dag = resultaten_verwelkingspunt['laatste_dag_boven_grens']
    uitval = resultaten_verwelkingspunt['vegetatie_uitval']
    cover = resultaten_verwelkingspunt['maatgevend_cover']
    dag_cover = resultaten_verwelkingspunt['dag_maatgevend_cover']
    cover_1oktober = resultaten_verwelkingspunt['cover_1_oktober']
    
    p = df_water['p'].sum()
    ETC = df_water['ETc'].sum()
    runoff = df_water['runoff'].sum()

    resultaten = [{
        'Locatie': dijkvak_id,
        'Latitude': lat,
        'Longitude': lon,
        'Jaar': jaar,
        'Ensemble': ens,
        'Grondsoort': grondsoort,
        'Totale neerslag (mm)': p,
        'Totale verdamping (mm)': ETC,
        'Totale runoff (mm)': runoff,
        'Langste uitvalperiode (dagen)': max_dagen,
        'Totaal dagen onder TAW': totaal_dagen,
        'Eerste dag onder TAW': eerste_dag.strftime('%Y-%m-%d') if eerste_dag else None,
        'Laatste dag onder TAW': laatste_dag.strftime('%Y-%m-%d') if laatste_dag else None,
        'Vegetatie-uitval': "Ja" if uitval else "Nee",
        'LivingPlantCover (%)': round(cover, 2),
        'Bijbehorende datum': dag_cover,
        'LPC [01-10] (%)': cover_1oktober,
    }]
    resultaten[0].update({
        'Totale verdamping (mm)': df_water['ETa'].sum(),
        'Totale potentiële gewasverdamping ETc (mm)': ETC,
        'Totale drainage (mm)': df_water['DP'].sum(),
        'Max absolute waterbalansresidu (mm)': df_water['water_balance_residual_mm'].abs().max(),
    })
    if options.enable_ks_mortality:
        resultaten[0].update({
            'Cumulatieve Ks-stress (dag-equivalenten)': resultaten_verwelkingspunt['cumulatieve_stress'],
            'Maximale mortaliteit (%)': resultaten_verwelkingspunt['max_mortaliteit'],
            'Mortaliteitsmethode': 'cumulatieve Ks',
        })
    return pd.DataFrame(resultaten)

def select_extreme_cases(resultaten_df, lat, lon, dijkvak_id):
    """
    Selecteert de heftigste jaren + ensembles o.b.v. laagste vegetatiedekking
    en toont de kernresultaten.
    """
    
    resultaten_df = resultaten_df[resultaten_df['Locatie'] == dijkvak_id]
    min_cover = resultaten_df["LivingPlantCover (%)"].min()
    kandidaten = resultaten_df[resultaten_df["LivingPlantCover (%)"] == min_cover]

    tie_col = ("Cumulatieve Ks-stress (dag-equivalenten)"
               if "Cumulatieve Ks-stress (dag-equivalenten)" in kandidaten
               else "Totaal dagen onder TAW")
    max_dagen = kandidaten[tie_col].max()
    zwaarste = kandidaten[kandidaten[tie_col] == max_dagen]

    kolommen = [
        'Locatie', 
        'Latitude',
        'Longitude',
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

import os

def main_per_ensembles(path, em_scen, ens_idx, lat, lon, jaar, grondsoort, dijkvak_id, helling, cot, output_path=None, options=DEFAULT_OPTIONS):

    # Weerdata laden
    df_weather = load_weather_from_nc(path=path, em_scen=em_scen, ens=ens_idx, lat=lat, lon=lon, jaar=jaar)

    # Waterbalans berekenen
    df_water, TAW, resultaten_verwelkingspunt = run_dagmodel(df_weather, grondsoort, lat, helling, options)
    max_dagen = resultaten_verwelkingspunt['max_dagen_boven_grens']
    totaal_dagen = resultaten_verwelkingspunt['totaal_boven_grens']
    eerste_dag = resultaten_verwelkingspunt['eerste_dag_boven_grens']
    laatste_dag = resultaten_verwelkingspunt['laatste_dag_boven_grens']
    uitval = resultaten_verwelkingspunt['vegetatie_uitval']
    cover = resultaten_verwelkingspunt['maatgevend_cover']
    # dag_cover = resultaten_verwelkingspunt['dag_maatgevend_cover']
    cover_1oktober = resultaten_verwelkingspunt['cover_1_oktober']

    print(f"Het percentage LPC op 1 oktober is {cover_1oktober:.2f} %")
    print(f"Langste aaneengesloten periode onder TAW: {max_dagen} dagen")
    print(f"Totaal aantal dagen onder TAW: {totaal_dagen}")

    if eerste_dag and laatste_dag:
        print(f"Eerste dag onder grens: {eerste_dag.strftime('%Y-%m-%d')}")
        print(f"Laatste dag onder grens: {laatste_dag.strftime('%Y-%m-%d')}")
    else:
        print("🟢 Geen enkele dag onder de TAW-grens")

    if uitval:
        print("Vegetatie-uitval gedetecteerd (LPC nabij nul)" if options.enable_ks_mortality else "Vegetatie-uitval gedetecteerd (>45 dagen aaneengesloten)")
    else:
        print("✅ Geen vegetatie-uitval")

    print(f"Droogte Periode = {max_dagen} dagen → LivingPlantCover = {cover:.2f}%")

    # Plot opslaan
    if output_path:
        plot(df_water, TAW, grondsoort, dijkvak_id, em_scen, cot, year=jaar, jaar=jaar, ens=ens_idx)

    # Resultaten teruggeven
    resultaten_df = analyse_per_jaar(path, em_scen, ens_idx, lat, lon, jaar, grondsoort, dijkvak_id, helling, options)
    return resultaten_df


def main_all_ensembles(path, em_scen, lat, lon, year, jaren, grondsoort, dijkvak_id, helling, cot, output_path=None, options=DEFAULT_OPTIONS):
    alle_resultaten = []
    kans_resultaten_dict = {}
    
    print(f"🔄 Start berekening voor emissiescenario {em_scen}, voor {grondsoort} met helling {cot}.")

    # 1) Runs draaien en samenvatten
    for la, lo, dv_id in zip(lat, lon, dijkvak_id):
        print("📍 ")
        print(f"\n ⏯️ Start analyse voor locatie met latitude: {la}, longitude: {lo} en dijk_id: {dv_id}.")
        resultaten_locatie = []
        
        for jaar in jaren:
            print(f"\r⏳ ... Analyseren voor jaar {jaar}", end="", flush=True)
            resultaten_per_jaar = []
            for ens in range(8):
                resultaten_df = analyse_per_jaar(
                    path, em_scen, ens, la, lo, jaar, grondsoort, dv_id, helling, options
                )
                resultaten_per_jaar.append(resultaten_df)
    
            jaar_df = pd.concat(resultaten_per_jaar, ignore_index=True)
            resultaten_locatie.append(jaar_df)    
            
        print()
        
        alle_resultaten_perlocatie = pd.concat(resultaten_locatie, ignore_index=True)
        alle_resultaten.append(alle_resultaten_perlocatie)
            
        # Kansberekeningen
        totaal = len(alle_resultaten_perlocatie)
        if totaal > 0:
            last_col = alle_resultaten_perlocatie["LivingPlantCover (%)"]

            # De logistieke functie bereikt wiskundig nooit exact nul.
            zero_label = f"≤{options.near_zero_pct:g}%" if options.enable_ks_mortality else "0%"
            n_zero = (last_col <= options.near_zero_pct).sum() if options.enable_ks_mortality else (last_col == 0).sum()
            n_high = ((last_col >= 70) & (last_col <= 100)).sum()
            n_medium = ((last_col >= 50) & (last_col < 70)).sum()
            n_low = (last_col < 50).sum()

            print(f"📉 Kans op {zero_label} LivingPlantCover ({grondsoort}): {n_zero}/{totaal} = {n_zero/totaal:.2%}")
            print(f"📉 Kans op 0-50% LivingPlantCover ({grondsoort}): {n_low}/{totaal} = {n_low/totaal:.2%}")
            print(f"📉 Kans op 50-70% LivingPlantCover ({grondsoort}): {n_medium}/{totaal} = {n_medium/totaal:.2%}")
            print(f"📉 Kans op 70-100% LivingPlantCover ({grondsoort}): {n_high}/{totaal} = {n_high/totaal:.2%}") 
            
            # 📊 opslaan voor Excel
            df_kans_maatgevend = pd.DataFrame({
                "Type LPC": "Maatgevend",
                "Categorie (%)": [zero_label, "0-50%", "50-70%", "70-100%"],
                "Aantal": [n_zero, n_low, n_medium, n_high],
                "Kans (%)": [n_zero/totaal*100, n_low/totaal*100, n_medium/totaal*100, n_high/totaal*100]
            })
            
            # --- LivingPlantCover op 1 oktober ---
            okt_col = alle_resultaten_perlocatie['LPC [01-10] (%)']

            n_zero_okt = (okt_col.round(2) <= options.near_zero_pct).sum() if options.enable_ks_mortality else (okt_col == 0).sum()
            n_high_okt = ((okt_col >= 70) & (okt_col <= 100)).sum()
            n_medium_okt = ((okt_col >= 50) & (okt_col < 70)).sum()
            n_low_okt = (okt_col < 50).sum()

            print(f"📉 Kans op {zero_label} LivingPlantCover op 1 okt ({grondsoort}): {n_zero_okt}/{totaal} = {n_zero_okt/totaal:.2%}")
            print(f"📉 Kans op 0-50% LivingPlantCover op 1 okt ({grondsoort}): {n_low_okt}/{totaal} = {n_low_okt/totaal:.2%}")
            print(f"📉 Kans op 50-70% LivingPlantCover op 1 okt ({grondsoort}): {n_medium_okt}/{totaal} = {n_medium_okt/totaal:.2%}")
            print(f"📉 Kans op 70-100% LivingPlantCover op 1 okt ({grondsoort}): {n_high_okt}/{totaal} = {n_high_okt/totaal:.2%}") 

            df_kans_okt = pd.DataFrame({
                "Type LPC": "1 oktober",
                "Categorie (%)": [zero_label, "0-50%", "50-70%", "70-100%"],
                "Aantal": [n_zero_okt, n_low_okt, n_medium_okt, n_high_okt],
                "Kans (%)": [n_zero_okt/totaal*100, n_low_okt/totaal*100, n_medium_okt/totaal*100, n_high_okt/totaal*100]
            })
            
            # Opslaan in dictionary (beide datasets onder elkaar)
            kans_resultaten_dict[dv_id] = pd.concat([df_kans_maatgevend, df_kans_okt], ignore_index=True)
            
        print(f"\n ⏹️ Einde analyse voor locatie met latitude: {la}, longitude: {lo} en dijk_id: {dv_id}.")

    print(f"\n ✅ Einde berekening voor emissiescenario {em_scen}, voor {grondsoort} met helling {cot}.")
    
    # Alles samenvoegen
    alle_resultaten_df = pd.concat(alle_resultaten, ignore_index=True)
    # print(alle_resultaten_df)

    # 2) CSV met ALLE jaren (per opgegeven grondsoort)
    if output_path:
        os.makedirs(output_path, exist_ok=True)

        outfile = os.path.join(
            output_path, f"{em_scen}_{grondsoort}_{cot}_N_alle_jaren.csv"
        )
        alle_resultaten_df.to_csv(outfile, index=False, encoding="utf-8")
        print(f"💾 Opgeslagen: {em_scen}_{grondsoort}_{cot}_N_alle_jaren.csv")
        
        # Excel met meerdere tabbladen per dijk_id
        
        # outfile_excel = os.path.join(output_path, f"{em_scen}_{grondsoort}_{cot}_kans_LCP_per_locatie.xlsx")
        # with pd.ExcelWriter(outfile_excel, engine='openpyxl') as writer:
        #     for dv_id, df_kans in kans_resultaten_dict.items():
        #         # elk tabblad krijgt als naam de dijk_id (max 31 tekens voor Excel)
        #         sheet_name = str(dv_id)[:31]
        #         df_kans.to_excel(writer, sheet_name=sheet_name, index=False)
        # print(f"💾 Excel met LivingPlantCover per locatie opgeslagen: {em_scen}_{grondsoort}_{cot}_kans_LCP_per_locatie.xlsx")

        kans_lijst = []
        for la, lo, dv_id in zip(lat, lon, dijkvak_id):
           df_kans = kans_resultaten_dict[dv_id].copy()
           # kolommen vooraan toevoegen
           df_kans.insert(0, "Longitude", lo)
           df_kans.insert(0, "Latitude", la)
           df_kans.insert(0, "Locatie_ID", dv_id)
           kans_lijst.append(df_kans)

        if kans_lijst:
           kans_samengevoegd = pd.concat(kans_lijst, ignore_index=True)

           outfile_csv = os.path.join(
               output_path, f"{em_scen}_{grondsoort}_{cot}_kans_LCP_per_locatie.csv"
           )
           kans_samengevoegd.to_csv(outfile_csv, index=False, encoding="utf-8")
           print(f"💾 CSV met LivingPlantCover per locatie opgeslagen: {em_scen}_{grondsoort}_{cot}_kans_LCP_per_locatie.csv")
    
    # 3) Extreemste gevallen per grondsoort vinden en PLOTTEN
    
    for la, lo, dv_id in zip(lat, lon, dijkvak_id):
        print(f"\n🔥 Extreemste gevallen voor dijk ID {dv_id} ({grondsoort}):\n")
        df_g = alle_resultaten_df[alle_resultaten_df['Grondsoort'] == grondsoort]
        if df_g.empty:
            print("Geen data voor deze grondsoort.")
        
        extreem_df = select_extreme_cases(df_g, lat=la, lon=lo, dijkvak_id=dv_id)
        if extreem_df.empty:
            print("Geen extreme gevallen gevonden.")
        
        for _, row in extreem_df.iterrows():
            jaar_ext = int(row["Jaar"])
            ens_ext = int(row["Ensemble"])
        
            df_weather = load_weather_from_nc(path, em_scen, ens_ext, lat=la, lon=lo, jaar=jaar_ext)
            df_water, TAW, resultaten_verwelkingspunt = run_dagmodel(df_weather, grondsoort, la, helling, options)
            # max_dagen = resultaten_verwelkingspunt['max_dagen_boven_grens']
            # totaal_dagen = resultaten_verwelkingspunt['totaal_boven_grens']
            # eerste_dag = resultaten_verwelkingspunt['eerste_dag_boven_grens']
            # laatste_dag = resultaten_verwelkingspunt['laatste_dag_boven_grens']
            # uitval = resultaten_verwelkingspunt['vegetatie_uitval']
            # cover = resultaten_verwelkingspunt['maatgevend_cover']
            # dag_cover = resultaten_verwelkingspunt['dag_maatgevend_cover']
            # cover_1april = resultaten_verwelkingspunt['cover_1_april']
            cover_1oktober = resultaten_verwelkingspunt['cover_1_oktober']
            
            # df_data['Locatie_ID'] = dv_id
            # df_data['Jaar'] = year
            # df_data['LivingPlantCover (%)'] = resultaten_verwelkingspunt["cover_per_dag"]
            # df_data['cover_1_10'] = resultaten_verwelkingspunt["cover_1_oktober"]
            # df_data['maatgevend_cover'] = resultaten_verwelkingspunt['maatgevend_cover']
            
            # overzicht_per_dag = df_data
            
            # if output_path:
            #     outfile_daily = os.path.join(
            #         output_path, f"{em_scen}_{grondsoort}_{cot}_LivingPlantCover_per_dag_{dv_id}.csv"
            #     )
            #     overzicht_per_dag.to_csv(outfile_daily, index=False, encoding="utf-8")
            

            print(f"Het percentage LPC op 1 oktober is {cover_1oktober:.2f} %")
        
        
            p = df_water['p'].sum()
            ETc = df_water['ETa'].sum()
            runoff = df_water['runoff'].sum()
            print(f'De totale neerslag in {jaar_ext} is {p:.2f} mm.')
            print(f'De totale verdamping in {jaar_ext} is {ETc:.2f} mm.')
            print(f'De totale runoff in {jaar_ext} is {runoff:.2f} mm.')
            
            #if output_path:
                #plot(df_water, TAW, grondsoort, dv_id, em_scen, cot, year, jaar=jaar_ext, ens=ens_ext)
                
            
