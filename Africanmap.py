import numpy as np
import sys
import regionmask
import pylab as plt
import cartopy.crs as ccrs      # map projections
import xarray as xr             # major tool to work with NetCDF data!
import netCDF4 as nc4
import pandas as pd

import warnings
warnings.filterwarnings('ignore')

#--------------------------------------------------------
# African country totals from Palmer et al, 2019
# For Jade Saunders
#
# Written by Paul Palmer with input from Mehliyat Sadiq
# 10/4/20
#--------------------------------------------------------

#Africa_extent = [-17, 51, -35, 38]
Africa_extent = [-17, 45, -35, 34] 

# read in a NetCDF file
fname = '/Users/ppalmer/Desktop/Africa_calculation/1x1files/L4_FLUXMAP_1x1_UoE_ln_v1.0.nc'


ds      = xr.open_dataset(fname)
drin    = ds['land'] # dataarray
latin   = ds['lat']  # dataarray
lonin   = ds['lon']  # dataarray
date    = ds['start_date']
nmonths = len(ds.n_months)

#--------------------------
# Get subset for a particular year
#--------------------------

instrument = 'OCO2_LN'
iyear = 2016

years = []
for ii in np.arange(nmonths): years.append(np.squeeze(date.data[ii][0]))
years = np.array(years)

daysinmonth = [31,28,31,30,31,30,31,\
               31,30,31,30,31,30,31]

ind = np.where(years == iyear)
ind = ind[0][:]

# subset of dataset, 12 months of global fluxes
dr_iyear = ds.land[ind,:,:] # use this to do calculations, dim = [12,180,360], unit: kgC/m2/s

#--------------------------
# Plot environment
#--------------------------
proj=ccrs.PlateCarree()
fig, (ax1,ax2,ax3) = plt.subplots(1,3,figsize=(15,4))

#--------------------------
# Plot mask over Africa
#--------------------------
ax1 = plt.subplot(131, projection=proj)

mask_countries = regionmask.defined_regions.natural_earth.countries_110.mask(dr_iyear)

mask_countries.plot(ax=ax1,transform=ccrs.PlateCarree(), cbar_kwargs={'shrink': 0.8,})

ax1.set_title("Countries 1:110m")
ax1.coastlines()
ax1.set_extent(Africa_extent, crs=ccrs.PlateCarree());

#--------------------------
# Make a smaller region
#--------------------------

mask_smaller     = mask_countries.sel(lat = slice(Africa_extent[2],Africa_extent[3]), \
                                      lon = slice(Africa_extent[0],Africa_extent[1]))
africa_countries = np.unique(mask_smaller.values) # unique values

ind_smaller = np.squeeze(np.where(~np.isnan(africa_countries)))
africa_countries_small = africa_countries[ind_smaller]

#print(africa_countries_small)
#print(len(ind_smaller))
#sys.exit()

#--------------------------
# Plot CO2 field over Africa
#--------------------------
#ax2 = plt.subplot(132, projection=proj)
#
#dr_iyear[5,:,:].plot(ax=ax2,transform=ccrs.PlateCarree(), cbar_kwargs={'shrink': 0.8})
#
#ax2.set_title("Countries 1:110m")
#ax2.coastlines()
#ax2.set_extent(Africa_extent, crs=ccrs.PlateCarree());

#--------------------------
# Area 
#--------------------------
infile = 'global_1x1_area.nc'
rootgrp = nc4.Dataset(infile, 'r')
aream2   = np.array(rootgrp.variables['area'])    
rootgrp.close()
aream2 = np.transpose(aream2)

#--------------------------
# Convert units of map
#--------------------------
for monthcount in np.arange(12):
    q = np.multiply(dr_iyear[monthcount,:,:],aream2) # kgC/s
    r = np.multiply(q,3600*24*1e3)                   #  gC/d
    s = np.multiply(r,daysinmonth[monthcount])       #  gC/month        
    dr_iyear[monthcount,:,:] = s

dr_annual = dr_iyear.sum(dim = 'n_months') * 1e-15 # unit: gC/month to PgC/year

#--------------------------
# Plot CO2 field over Africa
#--------------------------
ax2 = plt.subplot(132, projection=proj)

dr_annual.plot(ax=ax2,transform=ccrs.PlateCarree(), cbar_kwargs={'shrink': 0.8})

ax2.set_title("Countries 1:110m")
ax2.coastlines()
ax2.set_extent(Africa_extent, crs=ccrs.PlateCarree());
    
#--------------------------
# Quantify carbon budgets from individual countries
#--------------------------
# you can group over all integer values of the mask
# and compute mean (or total) values for each country
co2fluxbudget = dr_annual.groupby(mask_countries).sum('stacked_lat_lon')
country_names = regionmask.defined_regions.natural_earth.countries_110[co2fluxbudget.region.values].names
co2fluxbudget.coords['names'] = ('region', country_names)
co2flux_subset                = co2fluxbudget.sel(region = africa_countries_small) # get subset of the results

#--------------------------
# Prepare for plotting in alphabetical order
#--------------------------
d = {'country': co2flux_subset.names, 'co2': co2flux_subset.values}
df = pd.DataFrame(data=d)
df_sorted = df.sort_values(by=['country'])

# Remove from countries within coarsely defined domain that aren't based in Africa
not_african_countries = ['Iraq', 'Israel', 'Jordan', 'Lebanon', 'Saudi Arabia', 'Syria', 'Yemen']

for j in range(len(not_african_countries)):
    # Get names of indexes
    indexNames = df_sorted[df_sorted['country'] == not_african_countries[j]].index
    # Delete these row indexes from dataFrame
    df_sorted.drop(indexNames, inplace=True)

#--------------------------
# Plot estimates
#--------------------------

fig = plt.figure(figsize=[5,10])

plt.barh(df_sorted.country, df_sorted.co2)

for ii in np.arange(len(df_sorted.co2)):
    plt.plot([-0.2,0],[ii,ii],'--',color='red',alpha=0.1)

#plt.title('Number of countries: ' + str(len(df_sorted)))
plt.title('African carbon budgets: '+str(iyear))
plt.xlabel('CO$_2$ flux (PgC/yr)')
plt.xticks(rotation=90)
plt.xlim([-0.2,0.35])

CountryTotal = np.sum(df_sorted.co2)
OutString =  'Continental total: '+'{:4.2f}'.format(CountryTotal)+' PgC/yr'



plt.annotate(OutString, [0.0,51],ha='center')

fig.tight_layout()

plt.savefig('Palmer2019AfricanCO2budget_'+instrument+'_'+str(iyear)+'.png')

plt.show()
