'''
16/12/2021

Fig4 - plumes for CanESM5
FigS5 - plumes for CanESM5 and MIROC6
data generation scripts - 
 a) miroc_allreg_timeline_masked_101221.py
 b) miroc_run_plume_v2.py
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
from iris.analysis import geometry
import pandas as pd


# Fixed climatology masked
# load model data
sd1 = np.load('/user/work/hh21501/canesm_ts/sd1_clim_plume_mask.npy')
sd2 = np.load('/user/work/hh21501/canesm_ts/sd2_clim_plume_mask.npy')
sd3 = np.load('/user/work/hh21501/canesm_ts/sd3_clim_plume_mask.npy')
sd4 = np.load('/user/work/hh21501/canesm_ts/sd4_clim_plume_mask.npy')
sd5 = np.load('/user/work/hh21501/canesm_ts/sd5_clim_plume_mask.npy')

sd1_m = np.load('/user/work/hh21501/miroc_ts/sd1_clim_plume_mask.npy')
sd2_m = np.load('/user/work/hh21501/miroc_ts/sd2_clim_plume_mask.npy')
sd3_m = np.load('/user/work/hh21501/miroc_ts/sd3_clim_plume_mask.npy')
sd4_m = np.load('/user/work/hh21501/miroc_ts/sd4_clim_plume_mask.npy')
sd5_m = np.load('/user/work/hh21501/miroc_ts/sd5_clim_plume_mask.npy')


# load obs data
obs = pd.read_csv('era5_clim_mask.csv')
obs_year = obs['Year'].to_numpy()
obs_sd1 = obs['1SD'].to_numpy()/158*100
obs_sd2 = obs['2SD'].to_numpy()/158*100
obs_sd3 = obs['3SD'].to_numpy()/158*100
obs_sd4 = obs['4SD'].to_numpy()/158*100
obs_sd5 = obs['5SD'].to_numpy()/158*100

plt.ion()
plt.show()
# Plot as timeseries
fig = plt.figure(figsize=(12., 8.), dpi=80, num=None)
ax1 = plt.subplot2grid((2, 2), (0, 0), colspan=1, rowspan=1)
yr = np.arange(1950, 2101)


for ens in np.arange(49):
    plt.plot(yr, sd1[ens], color='cornflowerblue', alpha=.5, linewidth=.5)
    plt.plot(yr, sd2[ens], color='springgreen', alpha=.5, linewidth=.5)
    plt.plot(yr, sd3[ens], color='salmon', alpha=.5, linewidth=.5)
    plt.plot(yr, sd4[ens], color='darkturquoise', alpha=.5, linewidth=.5)
    plt.plot(yr, sd5[ens], color='plum', alpha=.5, linewidth=.5)


plt.plot(obs_year, obs_sd1, 'darkblue', linewidth=2, label='1SD')
plt.plot(obs_year, obs_sd2, 'g', linewidth=2, label='2SD')
plt.plot(obs_year, obs_sd3, 'r', linewidth=2, label='3SD')
plt.plot(obs_year, obs_sd4, 'teal', linewidth=2, label='4SD')
plt.plot(obs_year, obs_sd5, 'darkmagenta', linewidth=2, label='5SD')

#plt.legend()
plt.ylim([0, 100])
plt.xlim([1960, 2100])
#plt.title('Extremes above each threshold')
#plt.xlabel('year')
plt.text(2030, 88, 'Fixed climatology, 1981-2010')
plt.text(1950, 106, '(a)')
plt.ylabel('Percent of Regions')
plt.title('CanESM5')


    

ax2 = plt.subplot2grid((2, 2), (0, 1), colspan=1, rowspan=1)

for ens in np.arange(50):
    plt.plot(yr, sd1_m[ens]*1.068, color='cornflowerblue', alpha=.5, linewidth=.5)
    plt.plot(yr, sd2_m[ens]*1.068, color='springgreen', alpha=.5, linewidth=.5)
    plt.plot(yr, sd3_m[ens]*1.068, color='salmon', alpha=.5, linewidth=.5)
    plt.plot(yr, sd4_m[ens]*1.068, color='darkturquoise', alpha=.5, linewidth=.5)
    plt.plot(yr, sd5_m[ens]*1.068, color='plum', alpha=.5, linewidth=.5)

ax2.plot(yr, np.mean(sd1, axis=0), color='darkblue', linewidth=1)
ax2.plot(yr, np.mean(sd2, axis=0), color='g', linewidth=1)
ax2.plot(yr, np.mean(sd3, axis=0), color='r', linewidth=1)
ax2.plot(yr, np.mean(sd4, axis=0), color='teal', linewidth=1)
ax2.plot(yr, np.mean(sd5, axis=0), color='darkmagenta', linewidth=1)

ax2.plot(obs_year, obs_sd1, 'darkblue', linewidth=2)
ax2.plot(obs_year, obs_sd2, 'g', linewidth=2)
ax2.plot(obs_year, obs_sd3, 'r', linewidth=2)
ax2.plot(obs_year, obs_sd4, 'teal', linewidth=2)
ax2.plot(obs_year, obs_sd5, 'darkmagenta', linewidth=2)

#plt.legend()
plt.ylim([0, 100])
plt.xlim([1960, 2100])
#plt.title('Extremes above each threshold')
#plt.xlabel('year')
#plt.text(2030, 85, 'Fixed climatology, 1981-2010')
plt.text(1950, 106, '(b)')
plt.title('MIROC5')



# multiple ensembles lines
# Which regions to use, from Alan, keep the 1's (updated 21st Sep)
mask_regs = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,0,1,1,0,1,1,1,1,0,1,1,1,1,1,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,1,0,1,1,0,1,1,1,0,0,1,0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,1,1,1,1,1,1,0,0,1,1,1,1,1,1,0,0,1,0,0,0,1,0,1,1,1,1,1,0,1,0,1,0,0,0,1,1,1,0,1,0,0,1,1,0,0,0,1,0,1,1,0,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,0,0,0,1,1,0,1,1,0,1,0,1,1,0,0,1,1,0,0,1,1,0]


reg_data = []
for each in np.arange(237):
    # Region by regions
    print('region:')
    print(each)
    ann_max = []
    if mask_regs[each] == 0:
        pass # do nothing, regions not included in calculations
    if mask_regs[each] ==1:
        # Load data of maximum in each year in each ensemble (141 x 49)
        reg_data.append(np.load('/user/work/hh21501/canesm_ts/ann_max_index_reg'+str(each)+'.npy'))
        
# Calculate % of regions above each sd/year
sd1_list=[]; sd2_list = []; sd3_list = []; sd4_list = []; sd5_list = []
for ens in np.arange(49):
    sd1 = []; sd2 = []; sd3 = []; sd4 = []; sd5 = []
    for yr in np.arange(141): # number of years
        year_list = []
        for each in reg_data:
            year_list.append(each.data[yr, ens])
        year_list = np.array(year_list)
        sd1.append((((len(year_list[year_list > 1]))/len(year_list)) *100))
        sd2.append((((len(year_list[year_list > 2]))/len(year_list)) *100))
        sd3.append((((len(year_list[year_list > 3]))/len(year_list)) *100))
        sd4.append((((len(year_list[year_list > 4]))/len(year_list)) *100))
        sd5.append((((len(year_list[year_list > 5]))/len(year_list)) *100))
    sd1_list.append(sd1)
    sd2_list.append(sd2)
    sd3_list.append(sd3)
    sd4_list.append(sd4)
    sd5_list.append(sd5)

# load obs data
obs = pd.read_csv('era5_run_mask.csv')
obs_year = obs['Year'].to_numpy()
obs_sd1 = obs['1SD'].to_numpy()/158*100
obs_sd2 = obs['2SD'].to_numpy()/158*100
obs_sd3 = obs['3SD'].to_numpy()/158*100
obs_sd4 = obs['4SD'].to_numpy()/158*100
obs_sd5 = obs['5SD'].to_numpy()/158*100

# Plot as timeseries
#fig = plt.figure(figsize=(10., 5.), dpi=80, num=None)
ax3 = plt.subplot2grid((2, 2), (1, 0), colspan=1, rowspan=1)
yr = np.arange(1960, 2101)
for ens in np.arange(49):
    plt.plot(yr, sd1_list[ens], color='cornflowerblue', alpha=.5, linewidth=.5)
    plt.plot(yr, sd2_list[ens], color='springgreen', alpha=.5, linewidth=.5)
    plt.plot(yr, sd3_list[ens], color='salmon', alpha=.5, linewidth=.5)
    plt.plot(yr, sd4_list[ens], color='darkturquoise', alpha=.5, linewidth=.5)
    plt.plot(yr, sd5_list[ens], color='plum', alpha=.5, linewidth=.5)

plt.plot(obs_year, obs_sd1, 'darkblue', linewidth=2)
plt.plot(obs_year, obs_sd2, 'g', linewidth=2)
plt.plot(obs_year, obs_sd3, 'r', linewidth=2)
plt.plot(obs_year, obs_sd4, 'teal', linewidth=2)
plt.plot(obs_year, obs_sd5, 'darkmagenta', linewidth=2)

#fig.legend(loc='center right', bbox_to_anchor=(0.9, 0.5), borderaxespad=0)
#fig.legend()
plt.ylim([0, 100])
plt.xlim([1960, 2100])
#plt.title('Highest daily maximum temperature each year, relative to standard deviation')
plt.xlabel('Year')
plt.ylabel('Percent of Regions')
plt.text(2030, 88, 'Moving climatology')
plt.text(1950, 106, '(c)')


## MIROC RUNNING MEAN
reg_data2 = []
mask_regs[85] = 0; mask_regs[107]=0
for each in np.arange(109):
    # Region by regions
    print('region:')
    print(each)
    ann_max = []
    if mask_regs[each] == 0:
        pass # do nothing, regions not included in calculations
    if mask_regs[each] ==1:
        # Load data of maximum in each year in each ensemble (141 x 49)
        if each < 160:
            x=np.load('/user/work/hh21501/miroc_ts2/ann_max_index_reg'+str(each)+'.npy')
        else:
            x=np.load('/user/work/hh21501/miroc_ts/sd1_run_plume_annmax'+str(each)+'.npy')
        reg_data2.append(x[-141:, :])

# Calculate % of regions above each sd/year
sd1m_list=[]; sd2m_list = []; sd3m_list = []; sd4m_list = []; sd5m_list = []
for ens in np.arange(50):
    sd1 = []; sd2 = []; sd3 = []; sd4 = []; sd5 = []
    for yr in np.arange(141): # number of years
        year_list = []
        for each in reg_data2:
            year_list.append(each.data[yr, ens]) # single ens, single year, every region
        year_list = np.array(year_list)
        sd1.append((((len(year_list[year_list > 1]))/len(year_list)) *100))
        sd2.append((((len(year_list[year_list > 2]))/len(year_list)) *100))
        sd3.append((((len(year_list[year_list > 3]))/len(year_list)) *100))
        sd4.append((((len(year_list[year_list > 4]))/len(year_list)) *100))
        sd5.append((((len(year_list[year_list > 5]))/len(year_list)) *100))
    sd1m_list.append(sd1)
    sd2m_list.append(sd2)
    sd3m_list.append(sd3)
    sd4m_list.append(sd4)
    sd5m_list.append(sd5)

# Plot as timeseries
#fig = plt.figure(figsize=(10., 5.), dpi=80, num=None)
ax4 = plt.subplot2grid((2, 2), (1, 1), colspan=1, rowspan=1)
yr = np.arange(1960, 2101)
for ens in np.arange(50):
    plt.plot(yr, sd1m_list[ens], color='cornflowerblue', alpha=.5, linewidth=.5)
    plt.plot(yr, sd2m_list[ens], color='springgreen', alpha=.5, linewidth=.5)
    plt.plot(yr, sd3m_list[ens], color='salmon', alpha=.5, linewidth=.5)
    plt.plot(yr, sd4m_list[ens], color='darkturquoise', alpha=.5, linewidth=.5)
    plt.plot(yr, sd5m_list[ens], color='plum', alpha=.5, linewidth=.5)

plt.plot(yr, np.mean(sd1_list, axis=0), color='darkblue', linewidth=1)
plt.plot(yr, np.mean(sd2_list, axis=0), color='g', linewidth=1)
plt.plot(yr, np.mean(sd3_list, axis=0), color='r', linewidth=1)
plt.plot(yr, np.mean(sd4_list, axis=0), color='teal', linewidth=1)
plt.plot(yr, np.mean(sd5_list, axis=0), color='darkmagenta', linewidth=1)

plt.plot(obs_year, obs_sd1, 'darkblue', linewidth=2)
plt.plot(obs_year, obs_sd2, 'g', linewidth=2)
plt.plot(obs_year, obs_sd3, 'r', linewidth=2)
plt.plot(obs_year, obs_sd4, 'teal', linewidth=2)
plt.plot(obs_year, obs_sd5, 'darkmagenta', linewidth=2)

fig.legend(loc='center right', bbox_to_anchor=(0.9, 0.4), borderaxespad=0)
#fig.legend()
plt.ylim([0, 100])
plt.xlim([1960, 2100])
#plt.title('Highest daily maximum temperature each year, relative to standard deviation')
plt.xlabel('Year')
#plt.ylabel('Percent of Regions')
#plt.text(2030, 85, 'Moving climatology')
plt.text(1950, 106, '(d)')

plt.tight_layout
#plt.savefig('figS5.png', format="png", dpi=300)
plt.close()



