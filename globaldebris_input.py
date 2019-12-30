#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Model input for intercomparison experiment
"""
# Built-in libraries
import os
# External libraries
import numpy as np
#import pandas as pd
import xarray as xr


#%%
# Main directory
main_directory = os.getcwd()
output_fp = main_directory + '/../output/'
output_ostrem_fp = main_directory + '/../output/ostrem_curves/'
output_ostrem_fn_sample = 'XXXX_debris_melt_curve.nc'

# Region of Interest Data (lat, long, elevation, hr of satellite data acquisition)
roi = 'HMA'
roi_dict = {'HMA':[45, 25, 105, 65],
            '01':[71, 50, 233, 180]
            }

metdata_fp = main_directory + '/../climate_data/'
metdata_elev_fn = 'ERA5_elev.nc'
mb_binned_fp = main_directory + '/../mb_data/Shean_2019_0213/mb_combined_20190213_nmad_bins/'
mb_binned_fp_wdebris = main_directory + '/../mb_data/Shean_2019_0213/mb_combined_20190213_nmad_bins/_wdebris/'

# Simulation data
startyear = '2000'
endyear = '2018'
timezone = 0
metdata_fn_sample = roi + '_ERA5-metdata-XXXX' + str(startyear) + '_' + str(endyear) + '.nc'
debris_elevstats_fullfn = main_directory + '/../hma_data/' + roi + '_debris_elevstats.nc'

# Latitude and longitude index to run the model
#  Longitude must be 0 - 360 degrees
#latlon_list_raw = 'all'
latlon_list_raw = [(28.1,86.7)]
#latlon_list_raw = [(44.,83.25)]
if latlon_list_raw == 'all':
    ds_elevstats = xr.open_dataset(debris_elevstats_fullfn)
    latidx_list, lonidx_list = np.where(ds_elevstats['zmean'] > 0)
    lat_list = ds_elevstats.latitude[latidx_list].values
    lon_list = ds_elevstats.longitude[lonidx_list].values
    latlon_list = list(tuple(zip(list(lat_list), list(lon_list))))
else:
    ds_elevstats = xr.open_dataset(debris_elevstats_fullfn)
    lat_list_raw = np.array([x[0] for x in latlon_list_raw])
    lon_list_raw = np.array([x[1] for x in latlon_list_raw])
    #  argmin() finds the minimum distance between the glacier lat/lon and the GCM pixel
    lat_nearidx = np.abs(lat_list_raw[:,np.newaxis] - ds_elevstats['latitude'][:].values).argmin(axis=1)
    lon_nearidx = np.abs(lon_list_raw[:,np.newaxis] - ds_elevstats['longitude'][:].values).argmin(axis=1)
    latlon_nearidx = list(zip(lat_nearidx, lon_nearidx))
    latlon_nearidx_unique = sorted(list(set(latlon_nearidx)))
    latidx_list_unique = [x[0] for x in latlon_nearidx_unique]
    lonidx_list_unique= [x[1] for x in latlon_nearidx_unique]
    lat_list = ds_elevstats.latitude[latidx_list_unique].values
    lon_list = ds_elevstats.longitude[lonidx_list_unique].values
    latlon_list = list(tuple(zip(list(lat_list), list(lon_list))))

#latlon_list = latlon_list[0:5]

#%%
# Simulation data
start_date = '2000-05-28'   # start date for debris_ts_model.py
end_date = '2018-05-28'     # end date for debris_ts_model.py
fn_prefix = 'Rounce2015_' + roi + '-'
#elev_cns = ['zmean']
elev_cns = ['zmean', 'zstdlow', 'zstdhigh']

# Output info
output_option = 2           # 1: csv of all fluxes and internal temps, 2: netcdf of melt and ts
if output_option == 2:
    mc_stat_cns = ['mean']
elif output_option == 3:
    mc_stat_cns = ['mean', 'std', '25%', '75%']
    print('\nSTOP!!!!! NEED TO STORE ATTRIBUTES FOR STATISTICS!!!!\n\n')
date_start = '20191227'

# ===== Debris properties =====
experiment_no = 3
# Debris thickness
#debris_thickness_all = np.array([0.02])
#debris_thickness_all = np.array([0.2])
#debris_thickness_all = np.array([0.2, 0.3])
#debris_thickness_all = np.arange(0,5.001,0.05)
debris_thickness_all = np.arange(0,3.001,0.1)
debris_thickness_all[0] = 0.02
# Surface roughness, thermal conductivity, and albedo
debris_properties_fullfn = main_directory + '/../hma_data/hma_debris_properties.csv'
debris_properties = np.genfromtxt(debris_properties_fullfn, delimiter=',', skip_header=1)
if experiment_no == 3:
    z0_random = np.array([0.016])
    k_random = np.array([1.])
    albedo_random = np.array([0.3])
elif experiment_no == 4:    
    z0_random = debris_properties[:,1]
    k_random = debris_properties[:,2]
    albedo_random = np.array([debris_properties[:,5]])
    print('\n\nNEED TO MAKE DEBRIS PROPERTIES RANDOM FOR MC SIMULATIONS\n\n')

#met_data_netcdf_fp = '/Volumes/LaCie_Raid/ERA5_hrly/'
#orog_data_fullfn = '/Users/davidrounce/Documents/Dave_Rounce/HiMAT/Climate_data/ERA5/ERA5_geopotential_monthly.nc'

#pkl_fp = main_directory + '/hma_data/pkl/'
#pkl_fullfn_sample = pkl_fp + roi + '_ERA5-metadata_' + start_date + '_'+ end_date + '_XXXX.pkl' 

# Dates of satellite temperature data
#ts_dates = ['2015-09-30']      
#ts_hr = roi_dict[roi][3]            # hour of satellite temperature data acquisition  
    
# Extra 
debris_albedo = 0.3     # -, debris albedo
za = 2                  # m, height of air temperature instrument
zw = 10                 # m, height of wind instrument

# Snow model parameters
option_snow_fromAWS = 0 # Switch to use snow depth as opposed to snow fall
option_snow = 1         # Switch to use snow model (1) or not (0)
Tsnow_threshold = 273.15      # Snow temperature threshold [K]
snow_min = 0.0001       # minimum snowfall (mwe) to include snow on surface; since the density of falling snow is 
                        # much less (~50-100 kg m-3) 0.0001 m of snow w.e. will produce 0.001 - 0.002 m of snow
rain_min = 0.0001

#%%
#debris_elev = roi_dict[roi][0]       # m a.s.l. of the debris modeling
lr_gcm = -0.0065        # lapse rate from gcm to glacier [K m-1]    
delta_t = 60*60         # s, time step of AWS
slope_AWS_deg = 0       # assuming AWS flat on top of ridge or Sky View Factor = 1
aspect_AWS_deg = 0      # deg CW from N

#%% Constants
row_d = 2700                # Density of debris (kg m-3)
c_d = 750                   # Specific heat capacity of debris (J kg-1 K-1)
I0 = 1368                   # Solar constant (W/m2)
transmissivity=0.75         # Vertical atmospheric clear-sky transmissivity (assumed)
emissivity = 0.95           # Emissivity of debris surface (Nicholson and Benn, 2006)
P0 = 101325                 # Standard Atmospheric Pressure (Pa) 
density_air_0 = 1.29        # Density of air (kg/m3)
density_water = 1000        # Density of water (kg/m3)
density_ice = 900           # Density of ice (kg/m3) (Nicholson and Benn, 2006)
lapserate = 0.0065          # Temperature Lapse Rate (K/m)
Kvk = 0.41                  # Von Karman's Constant
Lv = 2.49e6                 # Latent Heat of Vaporation of water (J/kg)
Lf = 3.335e5                # Latent Heat of Fusion of Water (J/kg)
Ls = 2.834e6                # Latent Heat of Sublimation (J/kg) 
cA = 1005                   # Specific Heat Capacity of Air (J kg-1 K-1)
cW = 4.179e3                # Specific Heat Capacity of Water (J kg-1 K-1)
cSnow = 2.09e3              # Specific Heat Capacity of Snow (J kg-1 K-1)
R_const = 461               # Gas Constant for water vapor (J kg-1 K-1)
Rd = 287                    # Dry gas constant (J kg-1 K-1)
stefan_boltzmann = 5.67e-8  # Sefan-Bolzmann constant (W m-2 K-4)

# Snow melt model constants
snow_c_v = 0.2              # sensitivity of visible albedo to snow surface aging (Tarboten and Luce, 1996)
snow_c_ir = 0.5             # sensitivity of near infrared albedo to snow surface aging (Tarboten and Luce, 1996)
albedo_vo = 0.85            # fresh snow reflectance for visible band
albedo_iro = 0.65           # fresh snow reflectance for near infrared band
snow_tau_0 = 1e6            # non-dimensional snow surface age constant [s]
z0_snow = 0.002             # Brock etal (2006)
emissivity_snow = 0.99      # emissivity of snow (Tarboten and Luce, 1996); Collier etal (2014) use 0.97
eS_snow = 610.5             # Saturated vapor pressure of snow (Pa) (Colbeck 1990)
k_snow = 0.10               # Rahimi and Konrad (2012), Sturm etal (2002), Sturm etal (1997)
#density_snow = 150         # Density of snow (kg/m3) - Lejeune et al. (2007)
#albedo_snow = 0.75         # Collier etal (2014)

# Newton-Raphson Method constants
n_iter_max = 100