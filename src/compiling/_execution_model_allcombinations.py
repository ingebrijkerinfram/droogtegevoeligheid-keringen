# -*- coding: utf-8 -*-
"""
Draait het droogtemodel voor alle combinaties van scenario, helling en grondsoort.

Paar opties:
    
- TESTRUN = True  -> alleen een paar cellen en jaren, output in de submap 'Testrun'
- TESTRUN = False -> alle cellen en alle 30 jaren, output in de submap 'Volledige run'
- Voor het testen van een specifiek jaartal kan je alleen deze in een van de parameters in INSTELLINGEN laten staan.

"""

import os
import time
import traceback
import warnings
from datetime import datetime

from src.paden import bepaal_paden, klimaatbestand, VARIABELEN
from src.parsing.parsing_droogtegevoeligheid_model_vervolg import read_coordinates
from src.processing.processing_droogtegevoeligheid_model_vervolg import main_all_ensembles

warnings.filterwarnings("ignore")


# =============================== INSTELLINGEN ===============================
TESTRUN = True

SCENARIOS    = ["ref", "Hn_2050", "Hd_2050", "Ln_2100", "Ld_2100", "Hn_2100", "Hd_2100"]
HELLINGEN    = [1/2, 1/3, 1/4, 1/5]
GRONDSOORTEN = ["zand", "klei"]

# Instellingen voor als je een testrun wilt doen
TEST_AANTAL_CELLEN = 5    # eerste N gridcellen uit NL_grid
TEST_AANTAL_JAREN  = 2    # eerste N jaren van elk scenario
# ============================================================================


if TESTRUN:
    aantal_cellen = TEST_AANTAL_CELLEN
    aantal_jaren = TEST_AANTAL_JAREN
    submap = "Testrun"
else:
    aantal_cellen = None   # None = alles
    aantal_jaren = None
    submap = "Volledige run"

JARENRANGE = {"ref": range(1991, 2021), "2050": range(2036, 2066), "2100": range(2086, 2116)}


def jaren_voor(em_scen):
    """Geeft (jaar, jarenrange) zoals in het compiling-script."""
    if em_scen == "ref":
        jaar = "2005"
        jarenrange = JARENRANGE["ref"]
    else:
        jaar = em_scen.split("_")[1]
        jarenrange = JARENRANGE[jaar]
    if aantal_jaren is not None:
        jarenrange = jarenrange[:aantal_jaren]
    return jaar, jarenrange


def cot_voor(helling):
    """1/2 -> 'cot2', 1/3 -> 'cot3', enz."""
    return f"cot{round(1 / helling)}"


def al_klaar(output_path, em_scen, grondsoort, cot):
    """True als beide outputbestanden van deze combinatie al bestaan."""
    namen = [f"{em_scen}_{grondsoort}_{cot}_N_alle_jaren.csv",
             f"{em_scen}_{grondsoort}_{cot}_kans_LCP_per_locatie.csv"]
    return all(os.path.exists(os.path.join(output_path, n)) for n in namen)


def schrijf_logboek(output_path, tekst):
    with open(os.path.join(output_path, "logboek.txt"), "a", encoding="utf-8") as f:
        f.write(f"{datetime.now():%Y-%m-%d %H:%M:%S}  {tekst}\n")


# ------------------------------------------------------------------ voorbereiden
start_totaal = time.time()
paden = bepaal_paden()
output_path = paden["output"] / submap
os.makedirs(output_path, exist_ok=True)

# Eerst controleren of de bestanden van ALLE scenario's er zijn,
# zodat de run niet na uren alsnog stopt op een ontbrekend bestand.
if not paden["locaties"].exists():
    raise FileNotFoundError(f"Locatiebestand niet gevonden: {paden['locaties']}")
meldingen = []
for em_scen in SCENARIOS:
    missend = [v for v in VARIABELEN if not klimaatbestand(paden, v, em_scen).exists()]
    if len(missend) == len(VARIABELEN):
        meldingen.append(f"  - {em_scen}: alle {len(VARIABELEN)} bestanden ontbreken")
    elif missend:
        meldingen.append(f"  - {em_scen}: ontbreekt {', '.join(missend)}")
if meldingen:
    raise FileNotFoundError(
        f"Niet alle klimaatbestanden staan in {paden['klimaatdata']}:\n" + "\n".join(meldingen)
        + "\nVoeg ze toe, of haal die scenario's (tijdelijk) uit de lijst SCENARIOS."
    )
print(f"Klimaatbestanden voor alle {len(SCENARIOS)} scenario's gevonden.")

lat, lon, dijkvak_id, _ = read_coordinates(grid_path=paden["locaties"].parent)
if aantal_cellen is not None:
    lat, lon, dijkvak_id = lat[:aantal_cellen], lon[:aantal_cellen], dijkvak_id[:aantal_cellen]

combinaties = [(s, g, h) for s in SCENARIOS for g in GRONDSOORTEN for h in HELLINGEN]
print(f"\n{'TESTRUN' if TESTRUN else 'VOLLEDIGE RUN'}: {len(combinaties)} combinaties, "
      f"{len(dijkvak_id)} gridcellen, output naar: {output_path}\n")

# ------------------------------------------------------------------ draaien
gelukt, overgeslagen, mislukt = 0, 0, []
duur_per_run = []

for nr, (em_scen, grondsoort, helling) in enumerate(combinaties, start=1):
    cot = cot_voor(helling)
    naam = f"{em_scen} | {grondsoort} | {cot}"

    if al_klaar(output_path, em_scen, grondsoort, cot):
        print(f"[{nr}/{len(combinaties)}] {naam}: output bestaat al, overgeslagen.")
        overgeslagen += 1
        continue

    print(f"\n{'=' * 70}\n [{nr}/{len(combinaties)}] {naam}\n{'=' * 70}")
    jaar, jarenrange = jaren_voor(em_scen)
    t0 = time.time()
    try:
        main_all_ensembles(
            path=paden["klimaatdata"],
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
        )
    except Exception:
        # Eén mislukte combinatie stopt niet de hele run
        mislukt.append(naam)
        print(f"{naam} mislukt, de run gaat verder. Zie foutmeldingen.txt")
        with open(os.path.join(output_path, "foutmeldingen.txt"), "a", encoding="utf-8") as f:
            f.write(f"\n=== {datetime.now():%Y-%m-%d %H:%M:%S}  {naam} ===\n{traceback.format_exc()}")
        schrijf_logboek(output_path, f"MISLUKT   {naam}")
        continue

    duur = time.time() - t0
    duur_per_run.append(duur)
    gelukt += 1
    schrijf_logboek(output_path, f"klaar     {naam}  ({duur / 60:.1f} min)")

    nog_te_doen = len(combinaties) - nr
    schatting = sum(duur_per_run) / len(duur_per_run) * nog_te_doen
    print(f"\n{naam} klaar in {duur / 60:.1f} min. "
          f"Nog {nog_te_doen} combinaties, geschat nog {schatting / 3600:.1f} uur.")

# ------------------------------------------------------------------ samenvatting
print(f"\n{'=' * 70}")
print(f"Klaar in {(time.time() - start_totaal) / 60:.1f} min. "
      f"Gelukt: {gelukt}, overgeslagen: {overgeslagen}, mislukt: {len(mislukt)}.")
for naam in mislukt:
    print(f"{naam}")
print(f"Output: {output_path}")
