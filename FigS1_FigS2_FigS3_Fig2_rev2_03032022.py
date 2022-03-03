''' 
Last editted 03/03/2022 - checked for code repository

supp fig 1 - distribution of CanESM5, MIROC6, and ERA5 for 1981-2010 & 2015-2024
supp fig 2 - distribution mean and stdev for CanESM5 and ERA5
fig 2 - chance of exceeding extreme 
supp fig 3 - GEV of CanESM5 data for WWA region

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

plt.ion()
plt.show()

## Exceedance
record = 40.16
lookup_lon = {'PacNW':[245, 250], 'Lytton':[236, 241], 'WWA':[237, 241]}
lookup_lat = {'PacNW':[40, 45], 'Lytton':[48, 53], 'WWA':[45, 52]}
name = 'WWA'



## Plots of distribution - is obs within ensembles?
# Load ERA5 data
obs = iris.load_cube('era_19602020.nc') # wwa region era5 daily max temp jja
obs = time_slice(obs, 1981, 2010) - 273.15
loc_cons = iris.Constraint(latitude=lambda cell: lookup_lat[name][0] < cell < lookup_lat[name][1]) \
                    & iris.Constraint(longitude=lambda cell: lookup_lon[name][0] < cell < lookup_lon[name][1])
season_cons = iris.Constraint(time=lambda cell: 6 <= cell.point.month <= 8)
year_cons = iris.Constraint(time=lambda cell: 2015 <= cell.point.year <= 2021)
obs2 = load_ERA5_tasmax(loc_cons & season_cons & year_cons)
obs2 = obs2.collapsed(['longitude', 'latitude'], iris.analysis.MEAN)
obs2 = obs2 - 273.15

# Load CanESM5 data, two time periods
can1 = iris.load_cube('/user/home/hh21501/LargeEnsembles/USHeatwave/tmp/WWA_hist_CanESM.nc')
can1 = time_slice(can1, 1981, 2010) - 273.15
mod_rcp85 = iris.load_cube('/user/home/hh21501/LargeEnsembles/USHeatwave/tmp/WWA_CanESM.nc') # JJA, 2015-2100
can2 = time_slice(mod_rcp85, 2015, 2025) - 273.15
can2a = np.reshape((can2.data.data), (50*1012,))
can1a = np.reshape((can1.data.data), (50*2760,))

# Load MIROC6 data
mod1a = np.load('/user/home/hh21501/LargeEnsembles/USHeatwave/tmp/WWA_miroc_20002010.nc.npy')
mod1b = np.load('/user/home/hh21501/LargeEnsembles/USHeatwave/tmp/WWA_miroc_19901999.nc.npy')
mod1c = np.load('/user/home/hh21501/LargeEnsembles/USHeatwave/tmp/WWA_miroc_19811989.nc.npy')
mod1 = np.concatenate((mod1a, mod1b, mod1c))-273.15
mod2 = iris.load_cube('/user/home/hh21501/LargeEnsembles/USHeatwave/tmp/WWA_miroc_20152024.nc')
mod2a = np.reshape((mod2.data.data - 273.15).data, (50*920,)) 


# FIGS1
fig = plt.figure(figsize=(7., 5.), dpi=80, num=None)
ax1 = plt.subplot2grid((2, 1), (0, 0), colspan=1, rowspan=1)
x = 46
y = np.min(obs.data)-2
nbins = np.arange(0, x, x/100)
ax1.hist(can1a, bins=nbins, normed=True, histtype='stepfilled', \
             color='salmon', label='CanESM5') 
ax1.hist(mod1, bins=nbins, normed=True, histtype='stepfilled', \
             color='blue', alpha=.2, label='MIROC6')
ax1.hist(obs.data, bins=nbins, normed=True, histtype='stepfilled', \
             color='grey', alpha=0.5, label='ERA5')
plt.xticks([10, 20, 30, 40], color = 'k')
plt.xlim([0, x])
#plt.ylim([0, 0.012])
ax1.tick_params(axis='x', colors='k')
plt.yticks([])
#plt.xlabel('Daily maximum temperature, $^\circ$C')
plt.text(-3, .09, '(a)')
plt.text(1, .085, '1981-2010')
plt.legend(loc=1)

ax2 = plt.subplot2grid((2, 1), (1, 0), colspan=1, rowspan=1)
x = 46
y = np.min(obs.data)-2
nbins = np.arange(0, x, x/100)
ax2.hist(can2a, bins=nbins, normed=True, histtype='stepfilled', \
             color='salmon', label='CanESM5') 
ax2.hist(mod2a, bins=nbins, normed=True, histtype='stepfilled', \
             color='blue', alpha=.2,  label='MIROC6') 
ax2.hist(obs2.data, bins=nbins, normed=True, histtype='stepfilled', \
             color='grey', alpha=0.5, label='ERA5')
plt.xticks([10, 20, 30, 40], color = 'k')
plt.xlim([0, x])
ax2.tick_params(axis='x', colors='k')
plt.yticks([])
plt.xlabel('Daily maximum temperature, $^\circ$C')
plt.text(-3, .12, '(b)')
plt.text(1, .1, '2015-2024')

sub_axes = plt.axes([.73, .3, .15, .15])
sub_axes.hist(can2a, bins=nbins, normed=True, histtype='stepfilled', \
             color='salmon', label='CanESM5')
sub_axes.hist(mod2a, bins=nbins, normed=True, histtype='stepfilled', \
             color='blue', alpha=.5,  label='MIROC6') # doesn't show due to xlim
sub_axes.hist(obs2.data, bins=nbins, normed=True, histtype='stepfilled', \
             color='grey', alpha=0.5, label='ERA5')
sub_axes.set_xlim([37, 43])
sub_axes.set_xticks([38, 40, 42], color = 'k')
sub_axes.set_ylim([0, 0.001])
sub_axes.set_yticks([])

plt.tight_layout()
#plt.savefig('figs1_rev2.png', dpi=300) 
plt.close()

# FIGS2 - just CanESM5
mod1 = iris.load_cube('/user/home/hh21501/LargeEnsembles/USHeatwave/tmp/WWA_hist_CanESM.nc')
mod1 = time_slice(mod1, 1981, 2010) - 273.15

mod_mean = []
mod_std = []
for ens in np.arange(50):
    new_mod = mod1.data[ens,:]
    mod_mean.append(np.mean(new_mod))
    mod_std.append(np.std(new_mod))

mod_rcp85 = iris.load_cube('/user/home/hh21501/LargeEnsembles/USHeatwave/tmp/WWA_CanESM.nc') # JJA, 2015-2100
mod = time_slice(mod_rcp85, 2015, 2025) - 273.15

mod_mean2 = []
mod_std2 = []
for ens in np.arange(50):
    new_mod = mod.data[ens,:]
    mod_mean2.append(np.mean(new_mod))
    mod_std2.append(np.std(new_mod)) 

fig = plt.figure(figsize=(7., 5.), dpi=80, num=None)
ax1 = plt.subplot2grid((2, 2), (0, 0), colspan=1, rowspan=1)
y = np.min(obs.data)-2
nbins = np.arange(21, 23, 24/200)
ax1.hist(mod_mean, bins=nbins, normed=True, histtype='stepfilled', \
             color='salmon', label='CanESM5')
ax1.vlines(np.mean(obs.data), ymin=0, ymax=4, color='grey', alpha=0.5, label='ERA5')
ax1.set_xlim([21, 23])
ax1.set_ylim([0, 3])
ax1.set_xticks([22, 23])
ax1.set_yticks([])
ax1.text(20.8, 3.15, '(a)')
ax1.set_title('1981-2010')
ax1.set_ylabel('Mean')

ax2 = plt.subplot2grid((2, 2), (1, 0), colspan=1, rowspan=1)
nbins = np.arange(4, 6, 6/100)
ax2.hist(mod_std, bins=nbins, normed=True, histtype='stepfilled', \
             color='salmon', label='CanESM5')
ax2.vlines(np.std(obs.data), ymin=0, ymax=6, color='grey', alpha=0.5, label='ERA5')
ax2.set_xlim([4.5, 5.5])
ax2.set_ylim([0, 4])
ax2.set_yticks([])
ax2.text(4.4, 4.25, '(b)') 
ax2.set_ylabel('Std Dev')

ax3 = plt.subplot2grid((2, 2), (0, 1), colspan=1, rowspan=1)
y = np.min(obs.data)-2
nbins = np.arange(23.5, 25.5, 24/200)
ax3.hist(mod_mean2, bins=nbins, normed=True, histtype='stepfilled', \
             color='salmon', label='CanESM5')
ax3.vlines(np.mean(obs2.data), ymin=0, ymax=4, color='grey', alpha=0.5, label='ERA5')
ax3.set_xlim([23.5, 25.5])
ax3.set_ylim([0, 2])
ax3.set_xticks([24, 25])
ax3.set_yticks([])
ax3.text(23.3, 2.12, '(c)')
ax3.set_title('2015-2024')
plt.legend()

ax4 = plt.subplot2grid((2, 2), (1, 1), colspan=1, rowspan=1)
nbins = np.arange(4, 6, 6/100)
ax4.hist(mod_std2, bins=nbins, normed=True, histtype='stepfilled', \
             color='salmon')
ax4.vlines(np.std(obs2.data), ymin=0, ymax=6, color='grey', alpha=0.5)
ax4.set_xlim([4.5, 6])
ax4.set_ylim([0, 3])
ax4.set_yticks([])
ax4.set_xticks([5, 5.5])
ax4.text(4.35, 3.15, '(d)') 

#plt.savefig('figs2_rev.png', dpi=300)
plt.close()



# Fig2
# Generate data
'''
mod_rcp85 = iris.load_cube('/user/home/hh21501/LargeEnsembles/USHeatwave/tmp/WWA_CanESM.nc') # JJA, 2015-2100
mod_rcp45 = iris.load_cube('/user/home/hh21501/LargeEnsembles/USHeatwave/tmp/WWA_ssp245_CanESM.nc') # JJA, 2015-2100
mod_rcp16 = iris.load_cube('/user/home/hh21501/LargeEnsembles/USHeatwave/tmp/WWA_ssp126_CanESM.nc') # JJA, 2015-2100

def chance(data):
    return ((len(data[data > record]))/len(data)) *120

hi = []; mid = []; lo = []
for yr in np.arange(2015, 2090):
    y1 = yr; y2 = yr+9
    print(y1)
    hi.append(chance(cube_to_array(time_slice(mod_rcp85, y1, y2)) - 273.15))
    mid.append(chance(cube_to_array(time_slice(mod_rcp45, y1, y2)) - 273.15))
    lo.append(chance(cube_to_array(time_slice(mod_rcp16, y1, y2)) - 273.15))

np.save('/user/work/hh21501/canesm_ts/hi', hi)
np.save('/user/work/hh21501/canesm_ts/mid', mid)
np.save('/user/work/hh21501/canesm_ts/lo', lo)

# Range from ensemble members
hi_ens = []; mid_ens = []; lo_ens = []
for ens in np.arange(30, 49):
    hi = []; mid = []; lo = []
    for yr in np.arange(2015, 2090):
        y1 = yr; y2 = yr+9
        print(y1)
        hi.append(chance((time_slice(mod_rcp85[ens, :], y1, y2).data) - 273.15))
        mid.append(chance((time_slice(mod_rcp45[ens, :], y1, y2).data) - 273.15))
        lo.append(chance((time_slice(mod_rcp16[ens, :], y1, y2).data) - 273.15))
    hi_ens.append(hi)
    mid_ens.append(mid)
    lo_ens.append(lo)

np.save('/user/work/hh21501/canesm_ts/hi_ens3', hi_ens)
np.save('/user/work/hh21501/canesm_ts/mid_ens3', mid_ens)
np.save('/user/work/hh21501/canesm_ts/lo_ens3', lo_ens)
'''

hi = np.load('/user/work/hh21501/canesm_ts/hi.npy')
hi3 = np.load('/user/work/hh21501/canesm_ts/hi_ens3.npy')
hi1 = np.load('/user/work/hh21501/canesm_ts/hi_ens1.npy')
hi2 = np.load('/user/work/hh21501/canesm_ts/hi_ens2.npy')

mid = np.load('/user/work/hh21501/canesm_ts/mid.npy')
mid3 = np.load('/user/work/hh21501/canesm_ts/mid_ens3.npy')
mid1 = np.load('/user/work/hh21501/canesm_ts/mid_ens1.npy')
mid2 = np.load('/user/work/hh21501/canesm_ts/mid_ens2.npy')

lo = np.load('/user/work/hh21501/canesm_ts/lo.npy')
lo3 = np.load('/user/work/hh21501/canesm_ts/lo_ens3.npy')
lo1 = np.load('/user/work/hh21501/canesm_ts/lo_ens1.npy')
lo2 = np.load('/user/work/hh21501/canesm_ts/lo_ens2.npy')

hi_range = np.concatenate((hi1[:,0::10], hi2[:,0::10], hi3[:,0::10]), axis=0)
mid_range = np.concatenate((mid1[:,0::10], mid2[:,0::10], mid3[:,0::10]), axis=0)
lo_range = np.concatenate((lo1[:,0::10], lo2[:,0::10], lo3[:,0::10]), axis=0)

hi_max = np.max(hi_range, axis=0)
hi_min = np.min(hi_range, axis=0)
mid_max = np.max(mid_range, axis=0)
mid_min = np.min(mid_range, axis=0)
lo_max = np.max(lo_range, axis=0)
lo_min = np.min(lo_range, axis=0)

time = np.arange(2020, 2095)

fig = plt.figure(figsize=(5., 3.), dpi=80, num=None)
ax1 = plt.subplot2grid((1, 1), (0, 0), colspan=1, rowspan=1)
plt.plot(time[0::10]+.5, hi[0::10], 'ro', label='SSP585')
for i in np.arange(8):
    plt.vlines(time[0::10][i]+.5, hi_max[i], hi_min[i], 'r')

plt.plot(time[0::10], mid[0::10], 'mo', label='SSP245')
for i in np.arange(8):
    plt.vlines(time[0::10][i], mid_max[i], mid_min[i], 'm')

plt.plot(time[0::10]-.5, lo[0::10], 'co', label='SSP126')
for i in np.arange(8):
    plt.vlines(time[0::10][i]-.5, lo_max[i], lo_min[i], 'c')

plt.xlim([2018, 2092])
plt.ylim([0.005, 22])
ax1.spines['top'].set_color('grey')
ax1.spines['right'].set_color('grey')
ax1.spines['bottom'].set_color('grey')
ax1.spines['left'].set_color('grey') 
ax1.tick_params(axis='x', colors='k')
plt.ylabel('Chance, %')
plt.xlabel('Year')
plt.legend()
plt.yscale('log')
plt.yticks(ticks=[0.001, 0.01, 0.1, 1, 5, 10, 20],labels=[0.001, 0.01, 0.1, 1, 5, 10, 20])
plt.tight_layout()
#plt.savefig('fig2_rev_v3.png', format="png", dpi=300)
plt.close()


# SFig3
from scipy.stats import genextreme as gev

def return_levels_plot(distribution_pdf, x_values):
    '''
    Calculates probability of return levels 
    '''
    chance = []
    return_level = []
    for i, _ in enumerate(distribution_pdf):
        width = x_values[1] - x_values[0]
        P = []
        for a, b in zip(distribution_pdf[i:-1], distribution_pdf[i+1:]):
            P.append(((a+b) / 2) * width) 
        chance.append(sum(P)*100)
        return_level.append(x_values[i])
    return return_level, chance

def plot_points(data_array):
    obs_sort = np.sort(data_array)
    chance = []
    for i in range(len(obs_sort)):
        chance.append(((i+1) / len(obs_sort)) * 100)
    ret_per = []
    for each in chance:
        ret_per.append(100.*1/each)
    ret_per.reverse()
    plt.plot(ret_per, obs_sort, '+r', alpha=.5)

def plot_gev(data_array, min_val, max_val):
    shape, loc, scale = gev.fit(data_array)
    x_val = np.linspace(min_val, max_val, 1000)
    dist_pdf = gev.pdf(x_val, shape, loc, scale)
    ret_lev, chance = return_levels_plot(dist_pdf, x_val)
    chance2 = []
    for each in chance[:-1]:
        chance2.append(100.*1/each)
    plt.plot(chance2[10:], ret_lev[10:-1], 'r', label='GEV fit')

# CANESM PRESENT 2015-2024
mod_rcp85 = iris.load_cube('/user/home/hh21501/LargeEnsembles/USHeatwave/tmp/WWA_CanESM.nc') # JJA, 2015-2100
mod_rcp85 = time_slice(mod_rcp85, 2015, 2024)
mod_data2 = mod_rcp85.aggregated_by('year', iris.analysis.MAX)-273.15
mod = np.reshape(mod_data2.data, [50*10])

fig = plt.figure(figsize=(7., 5.), dpi=80, num=None)
plt.xscale('log')
plot_points(mod)
plot_gev(mod, 28, 42)
plt.hlines(40.1, 0, 10000, color='k', linestyles='--', label=('June 2021 event magnitude'))
plt.xlim([.5, 1000])
plt.xticks([2, 5, 10, 50, 100, 1000, 10000], [2, 5, 10, 50, 100, 1000, 10000])
plt.xlabel('Return Period, years')
plt.ylabel('Annual daily maximum temperature, celcius')
plt.legend()
plt.tight_layout()
#plt.savefig('figS3.png', dpi=300)
plt.close()

