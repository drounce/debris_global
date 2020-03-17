#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Model input for intercomparison experiment
"""
# Built-in libraries
import os
# External libraries
import numpy as np
import pandas as pd
#import xarray as xr

#%%
# Main directory
main_directory = os.getcwd()
output_fp = main_directory + '/output/'

# Region of Interest Data (lat, long, elevation, hr of satellite data acquisition)
roi = 'Ngozumpa'
roi_dict = {'Ngozumpa':[28.1,86.7, 5000, 4.75, '15.03473'],
            'Langtang':[28.25,85.6, 4800, 4.75, '15.04121'],
            'Baltoro':[35.75,76.3, 4000, 5.5, '14.06794'],
            'Koxkar':[41.8,80.1, 3400, 5.5, '13.43232']
            }

# Simulation data
experiment_no = 3
start_date = '2000-05-28'   # start date for debris_ts_model.py
end_date = '2018-05-28'     # end date for debris_ts_model.py
fn_prefix = 'Rounce2015_' + roi + '_'
debris_properties_fullfn = main_directory + '/hma_data/hma_debris_properties.csv'
met_data_fullfn = main_directory + '/hma_data/hma_ERA5-metdata_2000_2018_' + roi + '.csv'
met_data_netcdf_fp = '/Volumes/LaCie_Raid/ERA5_hrly/'
orog_data_fullfn = '/Users/davidrounce/Documents/Dave_Rounce/HiMAT/Climate_data/ERA5/ERA5_geopotential_monthly.nc'

# Dates of satellite temperature data
ts_dates = ['2015-09-30']      
ts_hr = roi_dict[roi][3]            # hour of satellite temperature data acquisition  


lat_AWS_deg = roi_dict[roi][0]      # deg, latitude - NGOZUMPA
lon_AWS_deg = roi_dict[roi][1]      # deg, longitude - NGOZUMPA
timezone = 0                        # hr, from UTC

lon_AWS_deg360 = lon_AWS_deg
if lon_AWS_deg < 0:
    lon_AWS_deg360 = 360 + lon_AWS_deg  
    
# Extra 
debris_albedo = 0.3     # -, debris albedo
za = 2                  # m, height of air temperature instrument
zw = 10                 # m, height of wind instrument

# Snow model parameters
option_snow_fromAWS = 1 # Switch to use snow depth as opposed to snow fall
option_snow = 1         # Switch to use snow model (1) or not (0)
Tsnow_threshold = 273.15      # Snow temperature threshold [K]

#%%
# Meteorological data
# Create CSV for model
if os.path.exists(met_data_fullfn):
    met_data_all = pd.read_csv(met_data_fullfn)
else:
    #%%
    import xarray as xr
    met_data_all = None
    years = list(np.arange(int(start_date.split('-')[0]), int(end_date.split('-')[0])+1))
    for nyear, year in enumerate(years):
        
#        met_data_fullfn = main_directory + '/hma_data/hma_ERA5-metdata_' + str(year) + '.csv'
        
        
        for nmonth, month in enumerate(list(np.arange(1,12+1))):
            met_data_netcdf_fn = 'ERA5_' + str(year) + '-' + str(month).zfill(2) + '.nc'
            #%%
            met_data_netcdf_fn = 'ERA5_2000-01.nc'
            
            print('processing ' + met_data_netcdf_fn)
            
            
            ds_fullfn = met_data_netcdf_fp + met_data_netcdf_fn
            met_data = xr.open_dataset(ds_fullfn)
            
            # Extract lat/lon indices only once
            if nyear + nmonth == 0:
                orog_data = xr.open_dataset(orog_data_fullfn)
                # Lat/lon index
                lat_nearidx = (np.abs(lat_AWS_deg - met_data['latitude'][:].values).argmin(axis=0))
                lon_nearidx = (np.abs(lon_AWS_deg - met_data['longitude'][:].values).argmin(axis=0))
                orog_lat_nearidx = (np.abs(lat_AWS_deg - orog_data['latitude'][:].values).argmin(axis=0))
                orog_lon_nearidx = (np.abs(lon_AWS_deg360 - orog_data['longitude'][:].values).argmin(axis=0))    
                
            elev_timeseries = orog_data['z'].values[:,orog_lat_nearidx,orog_lon_nearidx] / 9.80665
            elev = elev_timeseries[0]
            # wind speed
            u10 = met_data['u10'][:,lat_nearidx,lon_nearidx].values
            v10 = met_data['v10'][:,lat_nearidx,lon_nearidx].values
            wind10 = (u10**2 + v10**2)**0.5
            # incoming solar radiation
            sin = met_data.ssrd[:,lat_nearidx,lon_nearidx].values / 3600
            sin[sin < 0.1] = 0
            lin = met_data.strd[:,lat_nearidx,lon_nearidx].values / 3600
            # relative humidity ('Arden Buck equation': approximation from Bogel modification)
            d2m = met_data.d2m[:,lat_nearidx,lon_nearidx].values
            t2m = met_data.t2m[:,lat_nearidx,lon_nearidx].values
            rh = (100 * np.exp((18.678*(d2m-273.15) / (257.14 + (d2m-273.15))) - 
                               ((18.678 - (t2m-273.15) / 234.5) * ((t2m-273.15) / (257.14 + (t2m-273.15))))))
            # local time
            local_time = met_data['time'].values + np.timedelta64(timezone,'h')

            # MET DATA CSV FILE            
            met_data_csv_cns = ['YEAR', 'MONTH', 'DAY', 'HOUR', 'MINUTE', 'T', 'RH', 'FF', 'DIR', 'P', 'SWIN', 'SWOUT', 
                                'LWIN',  'LWOUT', 'PP', 'SNOW', 'T_z', 'RH_z', 'FF_z']
            met_data_csv = pd.DataFrame(np.zeros((met_data.time.shape[0], len(met_data_csv_cns))), 
                                                  columns=met_data_csv_cns)
            met_data_csv['YEAR'] = pd.Series(local_time).apply(lambda x: x.year).values
            met_data_csv['MONTH'] = pd.Series(local_time).apply(lambda x: x.month).values
            met_data_csv['DAY'] = pd.Series(local_time).apply(lambda x: x.day).values
            met_data_csv['HOUR'] = pd.Series(local_time).apply(lambda x: x.hour).values
            met_data_csv['MINUTE'] = pd.Series(local_time).apply(lambda x: x.minute).values            
            met_data_csv['T'] = met_data['t2m'][:,lat_nearidx,lon_nearidx].values - 273.15
            met_data_csv['RH'] = rh
            met_data_csv['FF'] = wind10
            met_data_csv['DIR'] = 0 # wind direction needs to be calculated
#            met_data_csv['P'] = met_data['sp'].values[:,lat_nearidx,lon_nearidx]
            met_data_csv['SWIN'] = sin
            met_data_csv['SWOUT'] = met_data_csv['SWIN'].values * debris_albedo # computed from albedo
            met_data_csv['LWIN'] = lin
            met_data_csv['LWOUT'] = 0 # outgoing longwave radiation computed from surface temperature
            tp = met_data['tp'][:,lat_nearidx,lon_nearidx].values * 1000
            tp[tp<0] = 0
            tp_rain = tp.copy()
            tp_rain[met_data_csv['T'].values + 273.15 <= Tsnow_threshold] = 0
            tp_snow = tp.copy()
            tp_snow[met_data_csv['T'].values + 273.15 > Tsnow_threshold] = 0
            met_data_csv['PP'] = tp_rain
            met_data_csv['SNOW'] = tp_snow
            met_data_csv['T_z'] = za
            met_data_csv['RH_z'] = za
            met_data_csv['FF_z'] = zw
            met_data_csv['ELEV'] = elev
            
            if met_data_all is None:
                met_data_all = met_data_csv
            else:
                met_data_all = pd.concat([met_data_all, met_data_csv], sort=False)
                
            #%%
            ds = xr.open_dataset(main_directory + '/hma_data/hma_ERA5-metdata_2000_2018-t2m.nc')
            ds['t2m'][0:24,68,87] - 273.15
            print('\nDELETE ME AND OTHERS AND REPLACE EXPORT\n')
                #%%
            
    # Export to csv
#    met_data_all.to_csv(met_data_fullfn, index=False)

#%%
#Elev_AWS = 2200        # m
debris_elev = roi_dict[roi][0]       # m a.s.l. of the debris modeling
lr_gcm = -0.0065        # lapse rate from gcm to glacier [K m-1]    
delta_t = 60*60         # s, time step of AWS
slope_AWS_deg = 0       # assuming AWS flat on top of ridge or Sky View Factor = 1
aspect_AWS_deg = 0      # deg CW from N
#spinup_days = 5
        
met_data_all['YYYY-MM-DD'] = [str(met_data_all.loc[x,'YEAR']) + '-' + str(met_data_all.loc[x,'MONTH']).zfill(2) + '-' + 
                              str(met_data_all.loc[x,'DAY']).zfill(2) for x in met_data_all.index.values]
start_idx = list(met_data_all['YYYY-MM-DD']).index(start_date)
end_idx = list(met_data_all['YYYY-MM-DD']).index(end_date) + 23
Elev_AWS = met_data_all['ELEV'].values[0]

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



#%% ===== ADJUST THE TIMEZONE =====
#df_fullfn = ('/Users/davidrounce/Documents/Dave_Rounce/DebrisGlaciers_WG/Melt_Intercomparison/rounce_model/hma_data/' + 
#             'hma_ERA5-metdata_2000_2018.csv')
#df = pd.read_csv(df_fullfn)
#df_columns = list(df.columns)
#df['date'] = [str(df.loc[x,'YEAR']) + str(df.loc[x,'MONTH']).zfill(2) + str(df.loc[x,'DAY']).zfill(2) 
#              + str(df.loc[x,'HOUR']).zfill(2) + str(df.loc[x,'MINUTE']) + '00' for x in df.index.values]
#import datetime
#date_datetime = [datetime.datetime.strptime(x, '%Y%m%d%H%M%S') for x in df['date'].values]
#df['date_pd'] = pd.to_datetime(date_datetime) - np.timedelta64(5,'h') - np.timedelta64(45,'m')
#
#df['YEAR'] = pd.Series(df.date_pd.values).apply(lambda x: x.year).values
#df['MONTH'] = pd.Series(df.date_pd.values).apply(lambda x: x.month).values
#df['DAY'] = pd.Series(df.date_pd.values).apply(lambda x: x.day).values
#df['HOUR'] = pd.Series(df.date_pd.values).apply(lambda x: x.hour).values
#df['MINUTE'] = pd.Series(df.date_pd.values).apply(lambda x: x.minute).values
#
#df2 = df.loc[:,df_columns]
#df_fullfn2 = ('/Users/davidrounce/Documents/Dave_Rounce/DebrisGlaciers_WG/Melt_Intercomparison/rounce_model/hma_data/' + 
#                'hma_ERA5-metdata_2000_2018_v2.csv')
#df_new = df2.to_csv(df_fullfn2, index=False)