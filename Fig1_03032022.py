from __future__ import division
import glob
import iris
import iris.coord_categorisation as icc
from iris.time import PartialDateTime
import numpy as np
import numpy.ma as ma
from datetime import date, timedelta
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
import iris.plot as iplt
import cartopy.crs as ccrs
import ssl


''' This script plots the time series of T '''
''' and heat stress metrcis '''
''' in the NW Pacific region '''
''' from ERA5 data '''
''' Eunice Lo '''
''' Created 02/08/2021 '''
''' Updated 04/08/2021, with Daithi's mask and an added distribution panel '''
''' Updated 05/08/2021, with z500 '''
''' Updated 12/08/2021, with v500 '''
''' Updated 30/11/2021, with box on map '''


def roll_cube(cube):

    """Takes a cube which goes longitude -180 to 180 back to 0 to 360."""

    lon = cube.coord("longitude")

    new_cube = cube.copy()
    new_cube.data = np.roll(cube.data, len(lon.points) // 2)
    new_cube.coord("longitude").points = lon.points + 180.
    
    if new_cube.coord("longitude").bounds is not None:
        new_cube.coord("longitude").bounds = lon.bounds + 180.
    
    return new_cube


def roll_cube_180(cube):

    """Takes a cube which goes longitude 0 to 360 back to -180 to 180."""

    lon = cube.coord("longitude")

    new_cube = cube.copy()
    new_cube.data = np.roll(cube.data, len(lon.points) // 2)
    new_cube.coord("longitude").points = lon.points - 180.

    if new_cube.coord("longitude").bounds is not None:
        new_cube.coord("longitude").bounds = lon.bounds - 180.

    return new_cube


if __name__ == '__main__':
    
    # trouble shooting
    ssl._create_default_https_context = ssl._create_unverified_context
 
    # paths
    pera5 = "/bp1store/geog-tropical/data/ERA-5/"

    # 1979-present June and July data
    # tasmax, this is from 1950!
    fjun = glob.glob(pera5+"day/tasmax/*_*06.nc")
    fjul = glob.glob(pera5+"day/tasmax/*_*07.nc")
    # v500, june july 2021
    fj21_v = "/user/home/yl17544/ERA5_extras/ERA5_v500_hour_20210607.nc"
    # z500, from 1950
    fjun_z = glob.glob(pera5+"day/z/z500/*_*06.nc")
    fjul_z = glob.glob(pera5+"day/z/z500/*_*07.nc")
    # precip, from 1950
    fps_1 = glob.glob(pera5+"day/total_precipitation/*_*06.nc")
    #fps_2 = glob.glob(pera5+"day/total_precipitation/*_*07.nc")
    # evaporation, from 1950
    fep_1 = glob.glob(pera5+"day/evaporation/*_*06.nc")
    #fep_2 = glob.glob(pera5+"day/evaporation/*_*07.nc")

    # last week of June and first week of July?
    twoweek_cons = iris.Constraint(time=lambda cell: PartialDateTime(month=6, day=24) <= cell.point <= PartialDateTime(month=7, day=7))    
    
    # first week of 1/2/3 months prior
    #prior_cons = iris.Constraint(time=lambda cell: PartialDateTime(day=1) <= cell.point <= PartialDateTime(day=14))

    # load data
    # tasmax
    pnw_jun_cubes = iris.load(fjun, twoweek_cons)
    pnw_jul_cubes = iris.load(fjul, twoweek_cons)
    # v500
    pnw_junjul_cube_v = iris.load_cube(fj21_v, twoweek_cons)
    # turn into daily mean
    icc.add_day_of_year(pnw_junjul_cube_v, "time")
    pnw_daily_cube_v = pnw_junjul_cube_v.aggregated_by("day_of_year", iris.analysis.MEAN)
    # z500, at 1 deg res
    pnw_jun_cubes_z1 = iris.load(fjun_z, twoweek_cons)
    pnw_jul_cubes_z1 = iris.load(fjul_z, twoweek_cons)
    pnw_jun_cubes_z = iris.cube.CubeList()
    pnw_jul_cubes_z = iris.cube.CubeList()
    # some files have other pressure levels
    for c6, c7 in zip(pnw_jun_cubes_z1, pnw_jul_cubes_z1):
        if c6.ndim==4:
            c6_new = c6.extract(iris.Constraint(air_pressure=500.))
            c6_new.remove_coord("air_pressure")
            pnw_jun_cubes_z.append(c6_new)
        else:
            pnw_jun_cubes_z.append(c6)
        if c7.ndim==4:
            c7_new = c7.extract(iris.Constraint(air_pressure=500.))
            c7_new.remove_coord("air_pressure")
            pnw_jul_cubes_z.append(c7_new)
        else:
            pnw_jul_cubes_z.append(c7)
    # precip 
    pnw_ps_1 = iris.load(fps_1)
    #pnw_ps_2 = iris.load(fps_2, twoweek_cons)
    # evaporation
    pnw_ep_1 = iris.load(fep_1)
    #pnw_ep_2 = iris.load(fep_2, june_cons)

    # concatenate for each year, tasmax, z500, p-e
    pnw_junjul_cubes = iris.cube.CubeList()
    pnw_junjul_cubes_z = iris.cube.CubeList()
    pnw_cubes_pme = iris.cube.CubeList()
    # tasmax is from 1950
    for n in range(len(pnw_jun_cubes)):
        cubelist = iris.cube.CubeList([pnw_jun_cubes[n], pnw_jul_cubes[n]])
        iris.util.equalise_attributes(cubelist)
        cube = cubelist.concatenate_cube()
        pnw_junjul_cubes.append(cube)
    
    # other variables from 1950
    for n in range(len(pnw_jun_cubes_z)):
        # z500
        cubelist_z = iris.cube.CubeList([pnw_jun_cubes_z[n], pnw_jul_cubes_z[n]])
        iris.util.equalise_attributes(cubelist_z)
        cube_z = cubelist_z.concatenate_cube()
        # turn into m
        cube_z.data //= 9.80665
        cube_z.units = "m"
        pnw_junjul_cubes_z.append(cube_z)    
        
        # precip
        if pnw_ps_1[n].ndim==4:
            new_1 = pnw_ps_1[n][:,0,:,:].copy()
            new_1.remove_coord("expver")
            cubelist_p = iris.cube.CubeList([new_1])
            #cubelist_p = iris.cube.CubeList([new_1, pnw_ps_2[n]])
        else:
            cubelist_p = iris.cube.CubeList([pnw_ps_1[n]])
            #cubelist_p = iris.cube.CubeList([pnw_ps_1[n], pnw_ps_2[n]])
        iris.util.equalise_attributes(cubelist_p)
        cube_p = cubelist_p.concatenate_cube() 

        # evaporation
        if pnw_ep_1[n].ndim==4:
            new_e = pnw_ep_1[n][:,0,:,:].copy()
            new_e.remove_coord("expver")
            cubelist_e = iris.cube.CubeList([new_e])
            #cubelist_e = iris.cube.CubeList([new_e, pnw_ep_2[n]])
        else:
            cubelist_e = iris.cube.CubeList([pnw_ep_1[n]])
            #cubelist_e = iris.cube.CubeList([pnw_ep_1[n], pnw_ep_2[n]])
        iris.util.equalise_attributes(cubelist_e)
        cube_e = cubelist_e.concatenate_cube()
        cube_e.units = "m"

        # p-e
        pnw_cubes_pme.append(cube_p-cube_e)

    # area of interest
    # WWA box
    lat_cons = iris.Constraint(latitude=lambda cell: 45 <= cell <= 52) 
    lon_cons = iris.Constraint(longitude=lambda cell: 237 <= cell <= 241)
    
    # Daithi's region
    # WRAF0.5-v4.1
    #stone_reg = iris.load_cube(pera5+"masks/region_05_regrid.nc", "region_area_fraction")[8]
    # WRAF2-v4.1
    #stone_reg = iris.load_cube(pera5+"region_regrid.nc", "region_area_fraction")[2]
    # revert lats and change lons to 0->360
    #stone_reg2 = roll_cube(iris.util.reverse(stone_reg, 0))
    # 721 x 1440 
    #stone_mask = np.where(stone_reg2.data>0., 1, stone_reg2.data)         # turn regional area fracs into 0 and 1
    #stone_mask = ~(np.array(stone_mask, dtype=bool))                      # reverse 0 and 1, then turn to T(1) & F(0)
    #stone_mask_big = np.repeat(stone_mask[np.newaxis, :, :], 14, axis=0)  # broadcast in time dimension, 14 days
    # 1 deg res   
    #stone_reg2_1deg = stone_reg2.regrid(pnw_junjul_cubes_z[0][0], iris.analysis.Nearest())
    #stone_mask_1deg = np.where(stone_reg2_1deg.data>0., 1, stone_reg2_1deg.data)
    #stone_mask_1deg_big = np.repeat(stone_mask_1deg[np.newaxis, :, :], 14, axis=0)
    
    # Northern Hemisphere box
    pna_lat_cons = iris.Constraint(latitude=lambda cell: 0 <= cell <= 90) 
    pna_lon_cons = iris.Constraint(longitude=lambda cell: 0 <= cell <= 360)

    # extract area of interest
    # WWA
    # tasmax
    pnw_junjul_reg = pnw_junjul_cubes.extract(lat_cons & lon_cons)
     
    # Daithi's region, involves masking
    #pnw_junjul_reg = iris.cube.CubeList()
    # tasmax
    #for c in pnw_junjul_cubes:
        #new_c = c.copy()
        #new_c.data = ma.array(data=c.data, mask=stone_mask_big)
        #pnw_junjul_reg.append(new_c)
      
    # area weighted average, taxmax in deg C
    pnw_junjul_reg[0].coord("latitude").guess_bounds()
    pnw_junjul_reg[0].coord("longitude").guess_bounds()
    pnw_wghs = iris.analysis.cartography.area_weights(pnw_junjul_reg[0]) 
    pnw_junjul_means = [c.collapsed(["latitude", "longitude"], iris.analysis.MEAN, weights=pnw_wghs)-273.15 for c in pnw_junjul_reg]
     
    # v winds 500 hPa in the Northern Hemisphere, in m/s
    # 29th June
    #date_cons = iris.Constraint(time=lambda cell: cell.point == PartialDateTime(month=6, day=29)) 
    #pnw_plot_cube_v = pnw_daily_cube_v.extract(date_cons & pna_lat_cons & pna_lon_cons)
    # two-week average
    pnw_nh_cube_v = pnw_daily_cube_v.extract(pna_lat_cons & pna_lon_cons)
    pnw_plot_cube_v = pnw_nh_cube_v.collapsed("time", iris.analysis.MEAN)

    # z500 in WWA box
    pnw_z_reg = pnw_junjul_cubes_z.extract(lat_cons & lon_cons) 

    # z500 in Daithi's region
    #pnw_z_reg = iris.cube.CubeList()
    #for cz in pnw_junjul_cubes_z:
        #new_cz = cz.copy()
        #new_cz.data = ma.array(data=cz.data, mask=stone_mask_1deg_big)
        #pnw_z_reg.append(new_cz)
    
    # area-weighted average, z500 in m
    pnw_z_reg[0].coord("latitude").guess_bounds()
    pnw_z_reg[0].coord("longitude").guess_bounds()
    pnw_wghs_1deg = iris.analysis.cartography.area_weights(pnw_z_reg[0])
    pnw_z_means = [mzc.collapsed(["latitude", "longitude"], iris.analysis.MEAN, weights=pnw_wghs_1deg) for mzc in pnw_z_reg] 

    # P-E in WWA box, in m
    pnw_pme_reg = pnw_cubes_pme.extract(lat_cons & lon_cons)    

    # P-E in Daithi's region, in m
    #pnw_pme_reg = iris.cube.CubeList()
    #for cpe in pnw_cubes_pme:
        #new_cpe = cpe.copy()
        #new_cpe.data = ma.array(data=cpe.data, mask=stone_mask_big)
        #pnw_pme_reg.append(new_cpe)
    
    # area-weighted average, P-E in m
    pnw_pme_reg[0].coord("latitude").guess_bounds()
    pnw_pme_reg[0].coord("longitude").guess_bounds() 
    pnw_wghs_june = iris.analysis.cartography.area_weights(pnw_pme_reg[0]) 
    pnw_pme_means = [mpec.collapsed(["latitude", "longitude"], iris.analysis.MEAN, weights=pnw_wghs_june) for mpec in pnw_pme_reg]
    #pnw_pme_means = [mpec.collapsed(["latitude", "longitude"], iris.analysis.MEAN, weights=pnw_wghs) for mpec in pnw_pme_reg] 
    
    # plot graph
    fig, axs = plt.subplots(nrows=2, ncols=4, gridspec_kw={'width_ratios': [4, 1, 0.5, 5]}, figsize=(8,6))
    
    # time series plot, tasmax
    ax1 = axs[0,0]
    sdate = date(2021,6,24)
    edate = date(2021,7,8)
    x = [sdate+timedelta(days=x) for x in range((edate-sdate).days)]
    
    # loop years
    for y in range(len(pnw_junjul_means)-2):
        ax1.plot(x, pnw_junjul_means[y].data, color="black", alpha=0.6)
    ax1.plot(x, pnw_junjul_means[-2].data, color="black", alpha=0.6, label="1950-2020")    # for the sake of single labeling
    ax1.plot(x, pnw_junjul_means[-1].data, color="red", linewidth=3, label="2021")

    # y-axis
    ax1.set_ylabel("Daily max temp. (deg C)")
    
    # x-axis
    fmt_week = mdates.DayLocator(interval=5)
    fmt_day = mdates.DayLocator()
    ax1.xaxis.set_major_locator(fmt_week)
    ax1.xaxis.set_major_formatter(mdates.DateFormatter('%m-%d'))
    ax1.xaxis.set_minor_locator(fmt_day)
    ax1.set_xlabel("Date")
    #fig.autofmt_xdate()
    ax1.legend(loc=0, fontsize=8) 
    ax1.set_title("(a)", loc="left")

    # rotated distribution, tasmax
    ax2 = axs[0,1]
    ax2.sharey=ax1
    # convert into arrays
    pre2021_data = np.zeros((len(pnw_junjul_means[:-1]), 14))
    for yc in range(len(pnw_junjul_means[:-1])):
        pre2021_data[yc, :] = pnw_junjul_means[yc].data
    pre2021_data = pre2021_data.flatten()   
    yr2021_data = pnw_junjul_means[-1].data

    # plot dist
    ax2.hist(pre2021_data, color="black", alpha=0.6, histtype='bar', ec='black', orientation="horizontal")
    ax2.axhline(np.max(yr2021_data), color="red", linewidth=3)  
    ax2.set_yticklabels([]) 
    ax2.tick_params(axis='y', which='major', pad=-20)
    ax2.set_xlabel("Frequency") 
    ax2.set_title("(b)", loc="left")

    # share y axis
    ax1.get_shared_y_axes().join(ax1, ax2)
    ax2.autoscale()

    # invisible ax to add space 
    axv = axs[0,2]
    axv.set_visible(False)
    
    # drought plot    
    ax3 = axs[0,3]
    plt.sca(ax3)       # set current axis
    cm = plt.cm.get_cmap("Reds", 20)

    # z500 vs p-e, two-week mean in each year
    plot_pme = []
    plot_z = []
    for cpme, cz in zip(pnw_pme_means, pnw_z_means):
        twoweek_mean_pme = cpme.collapsed("time", iris.analysis.MEAN)
        twoweek_mean_z = cz.collapsed("time", iris.analysis.MEAN)
        plot_pme.append(twoweek_mean_pme.data)
        plot_z.append(twoweek_mean_z.data)

    # max tasmax each two-week period of each year, 1950-2021
    #tmax_short = pnw_junjul_means[-43:]   # for 1979-2021
    plot_tmax = np.array([ct.collapsed("time", iris.analysis.MAX).data for ct in pnw_junjul_means])

    # seperate top five max tasmax
    plot_tmax_top = []
    plot_pme_top = []
    plot_z_top = []
    for n in range(5):
        top = np.max(plot_tmax)
        pos = np.argmax(plot_tmax)
        # put in top arrays
        plot_tmax_top.append(top)
        plot_pme_top.append(plot_pme[pos])
        plot_z_top.append(plot_z[pos])
        # remove from existing arrays
        plot_tmax = np.delete(plot_tmax, pos)
        plot_pme.remove(plot_pme[pos])
        plot_z.remove(plot_z[pos]) 
    
    # scatter plots
    # vmin=15, vmax=35 if using Daithi's region
    sc = plt.scatter(plot_pme, plot_z, c=plot_tmax, \
                     vmin=20, vmax=40, cmap=cm, marker="o", alpha=0.7)
    sc2 = plt.scatter(plot_pme_top, plot_z_top, c=plot_tmax_top, \
                      vmin=20, vmax=40, cmap=cm, marker="^", alpha=0.7)
    plt.colorbar(sc, orientation="vertical", label="Max daily max temp. (deg C)")
    ax3.set_ylabel("Z500 (m)")
    ax3.set_xlabel("June P-E (m)")
    ax3.set_title("(c)", loc="left")

    # v500 synoptic map
    # combined axis
    gs = axs[1,0].get_gridspec()
    for ax in axs[1,:]:
        ax.remove()

    ax5 = fig.add_subplot(gs[1,:])
    plt.sca(ax5)       # set current axis
    cf = iplt.contourf(pnw_plot_cube_v, levels=np.linspace(-10,10,11), cmap="RdBu_r", extend="both")
    plt.gca().coastlines()
    # WWA box
    crs_latlon = ccrs.PlateCarree()
    plt.plot([237,241], [45,45], "dimgray", transform=crs_latlon, linewidth=1)
    plt.plot([241,241], [45,52], "dimgray", transform=crs_latlon, linewidth=1)
    plt.plot([237,241], [52,52], "dimgray", transform=crs_latlon, linewidth=1)
    plt.plot([237,237], [45,52], "dimgray", transform=crs_latlon, linewidth=1)
    plt.colorbar(cf, orientation="horizontal", label="V winds at 500 hPa (m/s)", aspect=40, shrink=0.9) 
    plt.title("(d)", loc="left") 
    #plt.title("2021-06-29", loc="center")

    plt.tight_layout(w_pad=-2, h_pad=1.5)
    plt.savefig("/user/home/yl17544/scripts/PNW_extremes_paper/plots/WWAbox_junjul_1950-2021_tasmax_1950-2021_z500_vs_juneP-E_2021_v500_withbox_dpi300.png", format="png", dpi=300)
    plt.close("all")
    print("Saved figure")