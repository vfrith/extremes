"""
This script calculate the most extreme event for each region based on the moving mean method.

inputs
------

	extreme csv : csv file with columns: [year, month, day, 0,1, ..., 236]
	region netCDF : netCDF file of the just the region labels (not the raw netcdf from the Daithi Stone paper)

outputs
-------

	plot : global plot of the most extreme events

author : Emily Vosper 16.08.2021
"""

# impoty modules
import cartopy.crs as ccrs
import matplotlib.pyplot as plt
import xarray as xr
import pandas as pd
import seaborn as sns
# from cartopy.util import add_cyclic_point
sns.set_style("white")

# define filepaths
extreme_fp = '/scratch/al18709/heatwave/data/regional_monthly_data_236_regions_tasmax_alan.csv'
region_fp = '/scratch/al18709/heatwave/region_regrid_0.5.nc'

# open file
region_ds = xr.open_mfdataset(region_fp, parallel=True)
df = pd.read_csv(extreme_fp)

# define lats and lons for plotting
lats = region_ds.lat
lons = region_ds.lon

# extract only regional data and leave dates
regions = [str(i) for i in range(237)]
extremes  = df[regions]

# 1. Find the month with the highest Tmax from the whole record
max_idx = extremes.idxmax(axis=0)
max_month = list(df['month'][max_idx])

# 2. Subset those months and save to extremes
for column in extremes.columns:
	column_bool = list(df['month'] == max_month[int(column)])
	extremes[column] = extremes[column][column_bool]

# 3. For the years 1989-2021, calculate the mean and SD from the 10 years previous (i.e. 1989 is compared to 1979-1988)
# TODO: change this to 1959-2021
mean = pd.DataFrame()
sd = pd.DataFrame()
for year in range(1989,2022):
	# grab relevant decade
	bool = (df['year'].isin(range(year-10,year)))
	# calculate mean and standard deviation
	year_mean = extremes[bool].mean(skipna=True,axis=0)
	year_sd = extremes[bool].std(skipna=True,axis=0)
	# append to mean and sd arrays
	mean[str(year)] = year_mean.values
	sd[str(year)] = year_sd.values

# 4. Calculate the extremes metric for every day in the whole record
temperatures = pd.read_csv(extreme_fp)
extreme = pd.DataFrame()
for i in range(15539):
	print(i)
	yr = int(temperatures.loc[i]['year'])
	row = temperatures.iloc[i][regions]
	row.index = list(range(237))
	if yr in [1979,1980,1981,1982,1983,1984,1985,1986,1987,1988]:
		continue
	e = (row.sub(mean[str(yr)])).div(sd[str(yr)])
	extreme[str(i)] = e

# 5. Calculate what the most extreme day was/find the max
# extreme = (extremes.sub(mean)).div(sd)
extreme = extreme.transpose()
max_extreme = extreme.max(axis=0)
max_extreme.to_csv('region_extreme_means.csv')

# 6. plot 
region_extreme = region_ds.region
for region in range(237):
	print('region',region)
	max = max_extreme[region]
	print('max extreme: ',max)
	region_extreme = region_extreme.where(region_ds.region.values != region, max)

fig, ax = plt.subplots(subplot_kw={'projection': ccrs.PlateCarree()})
c = ax.contourf(lons,lats,region_extreme,20,cmap = sns.color_palette("flare", as_cmap=True), transform=ccrs.PlateCarree())

ax.outline_patch.set_linewidth(0.3)
cbar = plt.colorbar(c, shrink=0.54)
cbar.outline.set_linewidth(0.5)
cbar.ax.tick_params(labelsize=6,width=0.5)
plt.title('Extreme')

plt.savefig('extreme__236.png',dpi=600,bbox_inches='tight')
plt.clf()