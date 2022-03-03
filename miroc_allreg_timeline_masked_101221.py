'''
10 Dec 2021 

Assessing all regions for MIROC rcp8.5 to make timeseries like for ERA-5

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


# From Alan, keeps the 1's
mask_regs = [1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,0,1,1,0,1,1,1,1,0,1,1,1,1,1,0,0,1,0,0,1,0,0,1,0,0,0,0,0,0,0,1,0,0,0,0,0,0,1,0,1,0,1,1,0,1,1,1,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,0,1,1,0,1,1,1,1,1,1,0,0,1,1,1,1,1,1,0,0,1,0,0,0,1,0,1,1,1,1,1,0,1,0,1,0,0,0,1,1,1,0,1,0,0,1,1,0,0,0,1,0,1,0,0,1,1,1,1,1,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,0,0,0,1,1,0,1,1,0,1,0,1,1,0,0,1,1,0,1,1,1,0]


# Set climatology
ann_max = []
for each in np.arange(237):  
    print(each)
    if mask_regs[each] == 0:
        pass
    if mask_regs[each] ==1:
        file_list = glob.glob('/user/work/hh21501/miroc_regions/reg'+str(each)+'_*')
        cubes = iris.load(file_list)
        mod = cubes.concatenate()[0]
        iris.coord_categorisation.add_year(mod, 'time') # add coord year
        iris.coord_categorisation.add_month_number(mod, 'time')
        # Calculate climatological mean and sd
        year_cons = iris.Constraint(year = lambda cell: 1981 <= cell  <= 2010) # climatology period
        mod_clim = mod.extract(year_cons)
        month_mean = mod_clim.aggregated_by('month_number', iris.analysis.MEAN)
        month_mean = month_mean.collapsed('realization', iris.analysis.MEAN)
        max_val = month_mean.collapsed('month_number', iris.analysis.MAX).data # monthly max
        month_number = np.where(month_mean.data == max_val)[0][0]+1 # month number of max
        if month_number == 1:
            mod_clim = mod_clim.extract(iris.Constraint(time = lambda cell: cell.point.month  == 1 or cell.point.month  == 2 or cell.point.month  == 12))
        else: 
            mod_clim = mod_clim.extract(iris.Constraint(time = lambda cell: month_number-1 <= cell.point.month  <= month_number+1))
            clim_mean = np.mean(mod_clim.data)
        clim_sd = np.std(mod_clim.data)
        # Calculate extreme index & identify max each year
        mod_index = (mod - clim_mean)/clim_sd
        ann_max.append(mod_index.aggregated_by('year', iris.analysis.MAX))

# Calculate % of regions above each sd/year
sd1_list=[]; sd2_list = []; sd3_list = []; sd4_list = []; sd5_list = []
for ens in np.arange(50):
    sd1 = []; sd2 = []; sd3 = []; sd4 = []; sd5 = []
    for yr in np.arange(151): # number of years
        year_list = []
        for each in ann_max:
            year_list.append(each.data[ens, yr])
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



# save data
np.save('/user/work/hh21501/miroc_ts/sd1_clim_plume_mask', sd1_list)
np.save('/user/work/hh21501/miroc_ts/sd2_clim_plume_mask', sd2_list)
np.save('/user/work/hh21501/miroc_ts/sd3_clim_plume_mask', sd3_list)
np.save('/user/work/hh21501/miroc_ts/sd4_clim_plume_mask', sd4_list)
np.save('/user/work/hh21501/miroc_ts/sd5_clim_plume_mask', sd5_list)

'''

# Running mean
ann_max = []
for each in np.arange(237):
    print(each)
    if mask_regs[each] == 0:
        pass
    if mask_regs[each] ==1:
        file_list = glob.glob('/user/work/hh21501/miroc_regions/reg'+str(each)+'*')
        cubes = iris.load(file_list)
        mod = cubes.concatenate()[0]
        iris.coord_categorisation.add_year(mod, 'time') # add coord year
        iris.coord_categorisation.add_month_number(mod, 'time')
        # Running climatology - fixed month
        year_cons = iris.Constraint(year = lambda cell: 1981 <= cell  <= 2010) # climatology period
        mod_clim = mod.extract(year_cons)
        month_mean = mod_clim.aggregated_by('month_number', iris.analysis.MEAN)
        month_mean = month_mean.collapsed('realization', iris.analysis.MEAN)
        max_val = month_mean.collapsed('month_number', iris.analysis.MAX).data # monthly max
        month_number = np.where(month_mean.data == max_val)[0][0]+1 # month number of max
        print(mod)
        if month_number == 1:
            mod_clim = mod.extract(iris.Constraint(time = lambda cell: cell.point.month  == 1 or cell.point.month  == 2 or cell.point.month  == 12))
        else: 
            mod_clim = mod.extract(iris.Constraint(time = lambda cell: month_number-1 <= cell.point.month  <= month_number+1))
        yr_extreme = []
        # running climatology
        for ens in np.arange(50):
            ann_max_ens = []
            for yr in np.arange(1960, 2101):
                year_cons = iris.Constraint(year = lambda cell: yr-10 <= cell  <= yr) # climatology period
                mod_new = mod.extract(year_cons)
                clim_mean = np.mean(mod_new.data)
                clim_sd = np.std(mod_new.data)
                # Calculate extreme index & identify max each year
                single_year = mod.extract(iris.Constraint(year = lambda cell: cell == yr))
                mod_index = (single_year - clim_mean)/clim_sd
                yr_extreme.append(np.max(mod_index.data))
                ann_max_ens.append(yr_extreme)
            print('Add to ann max', ens)
            ann_max.append(ann_max_ens)



# Calculate % of regions above each sd/year
sd1_list=[]; sd2_list = []; sd3_list = []; sd4_list = []; sd5_list = []
for ens in np.arange(49):
    sd1 = []; sd2 = []; sd3 = []; sd4 = []; sd5 = []
    for yr in np.arange(141): # number of years
        year_list = []
        for each in ann_max:
            year_list.append(each.data[ens, yr])
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


# save data
np.save('/user/work/hh21501/miroc_ts/sd1_run_plume_mask', sd1_list)
np.save('/user/work/hh21501/miroc_ts/sd2_run_plume_mask', sd2_list)
np.save('/user/work/hh21501/miroc_ts/sd3_run_plume_mask', sd3_list)
np.save('/user/work/hh21501/miroc_ts/sd4_run_plume_mask', sd4_list)
np.save('/user/work/hh21501/miroc_ts/sd5_run_plume_mask', sd5_list)



'''
