'''
Editted 12 October 2021 
Plotting fig 2 - chance of exceeding observed extreme  
and supp fig 1 - distribution of CanESM5 and ERA5 for 1981-2010 & 2015-2024
@vikki.thompson
'''

# Load neccessary libraries
import iris
import iris.coord_categorisation as icc
from iris.coord_categorisation import add_season_membership
import numpy as np
import matplotlib.pyplot as plt
import iris.plot as iplt
import cartopy.crs as ccrs
import cartopy as cart
import glob
import matplotlib.cm as mpl_cm
import sys
from iris.experimental.equalise_cubes import equalise_attributes
import scipy.stats as sps

def load_ERA5_tasmax(constr):
    file_list = glob.glob('/bp1store/geog-tropical/data/ERA-5/day/tasmax/*')
    cubes = iris.load(file_list, constr)
    equalise_attributes(cubes) # to allow merge to one cube 
    return cubes.concatenate_cube()

def time_slice(cube, year1, year2):
    year_cons = iris.Constraint(time=lambda cell: year1 <= cell.point.year <= year2)
    return cube.extract(year_cons)

def cube_to_array(cube):
    return cube.data.reshape(np.shape(cube.data)[0]*np.shape(cube.data)[1])

############

### VARIABLES
#y1 = 1981
#y2 = 2010
#
#obs = iris.load_cube('tmp/WWA_ERA5_19792020.nc')
#obs = obs.collapsed(['longitude', 'latitude'], iris.analysis.MEAN)
#iris.coord_categorisation.add_year(obs, 'time') # add coord year
#obs = time_slice(obs, y1, y2) - 273.15
#
#mod1 = iris.load_cube('tmp/WWA_hist_CanESM.nc')
#mod1 = cube_to_array(time_slice(mod1, y1, y2)) - 273.15
#


## Exceedance
record = 39.5

mod_rcp85 = iris.load_cube('/home/hh21501/LargeEnsembles/USHeatwave/tmp/WWA_CanESM.nc') # JJA, 2015-2100
mod_rcp45 = iris.load_cube('/home/hh21501/LargeEnsembles/USHeatwave/tmp/WWA_ssp245_CanESM.nc') # JJA, 2015-2100
mod_rcp16 = iris.load_cube('/home/hh21501/LargeEnsembles/USHeatwave/tmp/WWA_ssp126_CanESM.nc') # JJA, 2015-2100

def chance(data):
    return ((len(data[data > record]))/len(data)) *120

hi = []; mid = []; lo = []
for yr in np.arange(2015, 2090):
    y1 = yr; y2 = yr+9
    print(y1)
    hi.append(chance(cube_to_array(time_slice(mod_rcp85, y1, y2)) - 273.15))
    mid.append(chance(cube_to_array(time_slice(mod_rcp45, y1, y2)) - 273.15))
    lo.append(chance(cube_to_array(time_slice(mod_rcp16, y1, y2)) - 273.15))

time = np.arange(2020, 2095)

## Figure 
fig = plt.figure(figsize=(5., 3.), dpi=80, num=None)
ax1 = plt.subplot2grid((1, 1), (0, 0), colspan=1, rowspan=1)
plt.plot(time[0::10], hi[0::10], 'r+', label='RCP 8.5')
plt.plot(time[0::10], mid[0::10], 'm+', label='RCP 4.5')
plt.plot(time[0::10], lo[0::10], 'c+', label='RCP 1.6')
#plt.xticks(time, labels, rotation=60, color = 'k')
plt.xlim([2018, 2092])
plt.ylim([0.005, 22])
ax1.spines['top'].set_color('grey')
ax1.spines['right'].set_color('grey')
ax1.spines['bottom'].set_color('grey')
ax1.spines['left'].set_color('grey') 
ax1.tick_params(axis='x', colors='k')
#plt.yticks([])
plt.ylabel('Chance, %')
#plt.title('Chance of exceeding the record each year')
plt.xlabel('Year')
plt.legend()
#ax1.yaxis.set_minor_formatter(matplotlib.ticker.ScalarFormatter())
#ax1.yaxis.set_minor_formatter(matplotlib.ticker.FormatStrFormatter("%.8f"))
plt.yscale('log')
plt.yticks(ticks=[0.001, 0.01, 0.1, 1, 5, 10, 20],labels=[0.001, 0.01, 0.1, 1, 5, 10, 20])
#plt.grid()
plt.tight_layout()
plt.savefig('future_chance_canesm.png', format="png", dpi=300)
plt.close()





### VARIABLES
y1 = 1981
y2 = 2010
obs = iris.load_cube('tmp/WWA_ERA5_19792020.nc')
obs = obs.collapsed(['longitude', 'latitude'], iris.analysis.MEAN)
iris.coord_categorisation.add_year(obs, 'time') # add coord year
obs = time_slice(obs, y1, y2) - 273.15
mod1 = iris.load_cube('tmp/WWA_hist_CanESM.nc')
mod1 = cube_to_array(time_slice(mod1, y1, y2)) - 273.15


' Dictionary of individual locations for assessing' 
lookup_lon = {'PacNW':[245, 250], 'Lytton':[236, 241], 'WWA':[237, 241]}
lookup_lat = {'PacNW':[40, 45], 'Lytton':[48, 53], 'WWA':[45, 52]}

name = 'WWA'
y1 = 2015
y2 = 2021

# Location constraint
loc_cons = iris.Constraint(latitude=lambda cell: lookup_lat[name][0] < cell < lookup_lat[name][1]) \
                    & iris.Constraint(longitude=lambda cell: lookup_lon[name][0] < cell < lookup_lon[name][1])
# Seasonal (JJA) constraint
season_cons = iris.Constraint(time=lambda cell: 6 <= cell.point.month <= 8)
# Years constraint
year_cons = iris.Constraint(time=lambda cell: y1 <= cell.point.year <= y2)

obs2 = load_ERA5_tasmax(loc_cons & season_cons & year_cons)
obs2 = obs2.collapsed(['longitude', 'latitude'], iris.analysis.MEAN)
#iris.coord_categorisation.add_year(obs, 'time') # add coord year
#obs = time_slice(obs, y1, y2) - 273.15
obs2 = obs2 - 273.15
mod_rcp85 = iris.load_cube('tmp/WWA_CanESM.nc') # JJA, 2015-2100
y1 = 2015; y2 = 2025
mod_rcp85_2020 = cube_to_array(time_slice(mod_rcp85, y1, y2)) - 273.15
mod = mod_rcp85_2020



plt.ion()

fig = plt.figure(figsize=(7., 5.), dpi=80, num=None)
ax1 = plt.subplot2grid((2, 1), (0, 0), colspan=1, rowspan=1)
x = 46
y = np.min(obs.data)-2
nbins = np.arange(0, x, x/100)
ax1.hist(mod1, bins=nbins, normed=True, histtype='stepfilled', \
             color='salmon', label='CanESM5')
ax1.hist(obs.data, bins=nbins, normed=True, histtype='stepfilled', \
             color='grey', alpha=0.5, label='ERA5')
plt.xticks([10, 20, 30, 40], color = 'k')
plt.xlim([0, x])
#plt.ylim([0, 0.012])
ax1.tick_params(axis='x', colors='k')
plt.yticks([])
#plt.xlabel('Daily maximum temperature, $^\circ$C')
plt.text(-3, .09, '(a)')
plt.text(1, .075, '1981-2010')
plt.legend(loc=1)


ax2 = plt.subplot2grid((2, 1), (1, 0), colspan=1, rowspan=1)
x = 46
y = np.min(obs.data)-2
nbins = np.arange(0, x, x/100)
ax2.hist(mod, bins=nbins, normed=True, histtype='stepfilled', \
             color='salmon', label='CanESM5')
ax2.hist(obs2.data, bins=nbins, normed=True, histtype='stepfilled', \
             color='grey', alpha=0.5, label='ERA5')
plt.xticks([10, 20, 30, 40], color = 'k')
plt.xlim([0, x])
#plt.ylim([0, 0.012])
ax2.tick_params(axis='x', colors='k')
plt.yticks([])
plt.xlabel('Daily maximum temperature, $^\circ$C')
plt.text(-3, .12, '(b)')
plt.text(1, .1, '2015-2024')
#plt.legend(loc=2)
sub_axes = plt.axes([.73, .3, .15, .15])
sub_axes.hist(mod, bins=nbins, normed=True, histtype='stepfilled', \
             color='salmon', label='CanESM5')
sub_axes.hist(obs2.data, bins=nbins, normed=True, histtype='stepfilled', \
             color='grey', alpha=0.5, label='ERA5')
sub_axes.set_xlim([37, 43])
sub_axes.set_xticks([38, 40, 42], color = 'k')
sub_axes.set_ylim([0, 0.001])
sub_axes.set_yticks([])

plt.tight_layout()
plt.savefig('bias_canesm_present.png', dpi=300)
plt.close()




