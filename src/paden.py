# -*- coding: utf-8 -*-
"""
Centraal padenbeheer voor het droogtegevoeligheidsmodel.

Er staan geen persoonlijke paden meer in de scripts. Paden worden zo bepaald:

1. config.ini in de root van de repo (persoonlijk, staat in .gitignore).
   Elke regel daarin overschrijft de standaard.
2. Anders zoekt het script zelf de gesynchroniseerde SharePoint-map van
   project 000687 onder je gebruikersmap, bijvoorbeeld
   C:\\Users\\<naam>\\Infram BV\\Infram Projecten - 000687 Vervolg stresstest droogte
   en gebruikt daarin de vaste indeling uit STANDAARD.

Gebruik in het hoofdscript:

    from src.paden import bepaal_paden, controleer_klimaatdata
    paden = bepaal_paden()
    controleer_klimaatdata(paden, "Hd_2050")
    paden["klimaatdata"], paden["locaties"], paden["output"]
"""
import configparser
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]          # src/paden.py -> repo-root
CONFIG = REPO / "config.ini"

PROJECTCODE = "000687"

# Vaste indeling binnen de projectmap (relatief t.o.v. de projectmap)
STANDAARD = {
    "klimaatdata": "Uitvoering/Methode/Input/KNMI Data Klimaatscenarios",
    "locaties": "Uitvoering/Methode/Input/Gridkaart Nederland/LatLong_grids.xlsx",
    "output": "Uitvoering/Methode/Output",
}

# KNMI'23-variabelen die het model per scenario inleest
VARIABELEN = ["hurs", "tas", "tasmax", "tasmin", "pr", "sfcwind", "pet", "rsds"]


def _lees_config():
    """Lees config.ini (als die bestaat). Lege waarden tellen als 'niet ingesteld'."""
    if not CONFIG.exists():
        return {}
    cp = configparser.ConfigParser(interpolation=None)
    cp.read(CONFIG, encoding="utf-8")
    if not cp.has_section("paden"):
        return {}
    return {k: v.strip().strip('"').strip("'") for k, v in cp["paden"].items() if v.strip()}


def _naar_pad(waarde):
    """Relatieve paden in config.ini zijn relatief t.o.v. de repo-root."""
    pad = Path(waarde).expanduser()
    return pad if pad.is_absolute() else (REPO / pad)


def _zoek_projectmap():
    """Zoek de gesynchroniseerde SharePoint-map van het project onder de gebruikersmap.

    Dekt zowel 'Synchroniseren' (~/Infram BV/...) als
    'Snelkoppeling toevoegen aan Mijn bestanden' (~/OneDrive - Infram BV/...).
    """
    home = Path.home()
    gevonden = sorted({p for p in home.glob(f"*Infram*/*{PROJECTCODE}*") if p.is_dir()})
    return gevonden[0] if gevonden else None


def bepaal_paden(verbose=True):
    """Geef een dict met de paden 'klimaatdata', 'locaties' en 'output'."""
    cfg = _lees_config()

    if "projectmap" in cfg:
        projectmap = _naar_pad(cfg["projectmap"])
        bron = "config.ini"
    else:
        projectmap = _zoek_projectmap()
        bron = "automatisch gevonden"

    paden = {}
    for naam, relatief in STANDAARD.items():
        if naam in cfg:
            paden[naam] = _naar_pad(cfg[naam])
        elif projectmap is not None:
            paden[naam] = projectmap / relatief
        else:
            raise FileNotFoundError(
                f"Geen pad voor '{naam}'. De projectmap {PROJECTCODE} is niet gevonden onder "
                f"{Path.home()}.\nSynchroniseer de SharePoint-projectmap, of maak een config.ini "
                f"aan (kopie van config.example.ini) en vul daar 'projectmap' of '{naam}' in."
            )

    if verbose:
        print("📁 Gebruikte paden:")
        if projectmap is not None:
            print(f"   projectmap  : {projectmap}  ({bron})")
        for naam, pad in paden.items():
            print(f"   {naam:<12}: {pad}")
    return paden


def klimaatbestand(paden, variabele, em_scen):
    """Pad naar één KNMI'23-bestand, bv. klimaatbestand(paden, 'hurs', 'Hd_2050')."""
    return paden["klimaatdata"] / f"{variabele}_{em_scen}_interp.nc"


def controleer_klimaatdata(paden, em_scen, variabelen=VARIABELEN):
    """Stop direct met een duidelijke melding als er KNMI-bestanden ontbreken."""
    if not paden["locaties"].exists():
        raise FileNotFoundError(f"Locatiebestand niet gevonden: {paden['locaties']}")

    ontbrekend = [klimaatbestand(paden, v, em_scen).name for v in variabelen
                  if not klimaatbestand(paden, v, em_scen).exists()]
    if ontbrekend:
        raise FileNotFoundError(
            f"In {paden['klimaatdata']} ontbreken {len(ontbrekend)} bestand(en) "
            f"voor scenario {em_scen}:\n  - " + "\n  - ".join(ontbrekend)
        )
    print(f"✅ Alle {len(variabelen)} klimaatbestanden voor {em_scen} gevonden.")
