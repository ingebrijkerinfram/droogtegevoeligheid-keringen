
import xarray as xr

# Laad het NetCDF-bestand
file_path = r'C:\Users\inge.brijker\OneDrive - Infram BV\Droogtegevoeligheid Keringen\KNMI files\tas_ref_interp.nc'
ds = xr.open_dataset(file_path)

# Print basisinformatie over het dataset
print(ds)

# Haal coördinaten op
lat = ds['lat']
lon = ds['lon']
time = ds['time']
ens = ds['ens']

# Toon basisinformatie over coördinaten
print("\nLatitude (lat):")
print(lat)

print("\nLongitude (lon):")
print(lon)

print("\nTijd (time):")
print(time)

print("\nEnsemble leden (ens):")
print(ens)

# Bekijk de eerste paar waarden van de temperatuurvariabele (tas)
print("\nTemperatuur (tas) voorbeeldwaarden:")
print(ds['tas'].isel(time=0, ens=0))  # eerste tijdstap en ensemble-lid