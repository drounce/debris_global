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
import xarray as xr

import globaldebris_input as input

#%%
print('Note: longitude is 0-360 degrees')

# Simulation data
output_metdata_sample = (input.metdata_fp + input.roi + '_ERA5-metdata-XXXX' + str(input.startyear) 
                         + '_' + str(input.endyear) + '.nc')
metdata_netcdf_fp = '/Volumes/LaCie_Raid/ERA5_hrly/'
orog_data_fullfn = '/Users/davidrounce/Documents/Dave_Rounce/HiMAT/Climate_data/ERA5/ERA5_geopotential_monthly.nc'
timezone = 0

#%%
ds = xr.open_dataset(input.metdata_fp +input.metdata_elev_fn)

lat_N = input.roi_latlon_dict[input.roi][0]
lat_S = input.roi_latlon_dict[input.roi][1]
lon_E = input.roi_latlon_dict[input.roi][2]
lon_W = input.roi_latlon_dict[input.roi][3]

lat_N_idx = 



#%% ===== MODIFY THIS TO EXTRACT SUBSET OF DATA ======================================================================
## Meteorological data
##for nlatlon, latlon in enumerate([input.latlon_list[0]]):
#for nlatlon, latlon in enumerate(input.latlon_list):
#    print(nlatlon, latlon)
#    
#    lat_deg = latlon[0]
#    lon_deg = latlon[1]
#    
#    output_metdata_fullfn = output_metdata_sample.replace('XXXX', str(int(lat_deg*100)) + 'N-' + str(int(lon_deg*100)) 
#                                                          + 'E-')
#    
#    if os.path.exists(output_metdata_fullfn):
#        print('\n\nFILE ALREADY EXISTS\n\n')
#    else:
#        # ===== Combine meteorological data =====
#        ds_all = None
#        years = list(np.arange(int(input.startyear), int(input.endyear)+1))
#        
#        for vn in ['t2m','tp', 'u10', 'v10', 'ssrd', 'strd', 'rh', 'z']:
##            print('\n',vn)
#        
#            for nyear, year in enumerate(years):
##            for nyear, year in enumerate([years[0]]):
##                print(year)   
#                
##                for nmonth, month in enumerate(list(np.arange(1,12+1))[0:2]):
#                for nmonth, month in enumerate(list(np.arange(1,12+1))):
#                    metdata_netcdf_fn = 'ERA5_' + str(year) + '-' + str(month).zfill(2) + '.nc'
#
#                    print('processing ' + vn + ' ', metdata_netcdf_fn)
#                    
#                    met_data = xr.open_dataset(metdata_netcdf_fp + metdata_netcdf_fn)
#        
#                    
#                    # Extract lat/lon indices only once
#                    if nyear + nmonth == 0:
#                        orog_data = xr.open_dataset(orog_data_fullfn)
#                        lat_idx = np.abs(lat_deg - orog_data['latitude'][:].values).argmin(axis=0)
#                        lon_idx = np.abs(lon_deg - orog_data['longitude'][:].values).argmin(axis=0)
#                        
#                        if ((abs(met_data.latitude.values[lat_idx] - orog_data.latitude.values[lat_idx]) > 0.01)
#                            or (abs(met_data.longitude.values[lon_idx] - orog_data.longitude.values[lon_idx]) > 0.01)):
#                            print('\n\nOROGRAPHY IS NOT LINED UP WITH CLIMATE DATA \n\n')
#                            print(met_data.latitude.values[lat_idx], met_data.longitude.values[lon_idx])
#                            print(orog_data.latitude.values[lat_idx], orog_data.longitude.values[lon_idx])
#                     
#                    if vn == 'rh':
#                        # relative humidity ('Arden Buck equation': approximation from Bogel modification)
#                        d2m = met_data.d2m[:,lat_idx,lon_idx].values
#                        t2m = met_data.t2m[:,lat_idx,lon_idx].values
#                        rh = (100 * np.exp((18.678*(d2m-273.15) / (257.14 + (d2m-273.15))) - 
#                                           ((18.678 - (t2m-273.15) / 234.5) * ((t2m-273.15) / (257.14 + (t2m-273.15))))))
#                        ds_roi_month = xr.Dataset({'rh': (['time'], rh)},
#                                                  coords={'time': met_data.time})
#                    elif vn == 'z':
#                        ds_roi_month = (orog_data[vn][0,lat_idx,lon_idx] / 9.80665).to_dataset()
#                    else:
#                        ds_roi_month = met_data[vn][:,lat_idx,lon_idx].to_dataset()
#                    
#                    # Concatenate dataset
#                    if nyear + nmonth == 0:
#                        ds_roi = ds_roi_month
#                    elif vn not in ['z']:
#                        ds_roi = xr.concat([ds_roi, ds_roi_month], 'time')
#                        
#            # Concatenate datasets                
#            if ds_all is None:
#                ds_all = ds_roi
#            else:
#                ds_all = ds_all.combine_first(ds_roi)
#            # Add attributes manually
#            if vn == 'tp':
#                ds_all.tp.attrs = {'units':'m', 'long_name':'Total precipitation'}
#            if vn == 'rh':
#                ds_all.rh.attrs = {'units':'%', 'long_name':'Relative humidity'}
#            elif vn == 'u10':
#                ds_all.u10.attrs = {'units': 'm s**-1', 'long_name': '10 metre U wind component'}
#            elif vn == 'v10':
#                ds_all.v10.attrs = {'units': 'm s**-1', 'long_name': '10 metre V wind component'}
#            elif vn == 'ssrd':
#                ds_all.ssrd.attrs = {'units':'J m**-2', 'long_name':'Surface solar radiation downwards', 
#                                   'standard_name': 'surface_downwelling_shortwave_flux_in_air'}
#            elif vn == 'strd':
#                ds_all.strd.attrs = {'units':'J m**-2', 'long_name':'Surface thermal radiation downwards'}
#            elif vn == 'z':
#                ds_all.z.attrs = {'units':'m a.s.l.', 'long_name':'Elevation', 'comment':'converted from geopotential'}
#               
#        
#        # Export array for each variable
#        if os.path.exists(input.metdata_fp) == False:
#            os.mkdir(input.metdata_fp)
#        ds_all.to_netcdf(output_metdata_fullfn)
#        
#        # =============================================================================================================


#%%
## Meteorological data
##for nlatlon, latlon in enumerate([input.latlon_list[0]]):
#for nlatlon, latlon in enumerate(input.latlon_list):
#    print(nlatlon, latlon)
#    
#    lat_deg = latlon[0]
#    lon_deg = latlon[1]
#    
#    output_metdata_fullfn = output_metdata_sample.replace('XXXX', str(int(lat_deg*100)) + 'N-' + str(int(lon_deg*100)) 
#                                                          + 'E-')
#    
#    if os.path.exists(output_metdata_fullfn):
#        print('\n\nFILE ALREADY EXISTS\n\n')
#    else:
#        # ===== Combine meteorological data =====
#        ds_all = None
#        years = list(np.arange(int(input.startyear), int(input.endyear)+1))
#        
#        for vn in ['t2m','tp', 'u10', 'v10', 'ssrd', 'strd', 'rh', 'z']:
##            print('\n',vn)
#        
#            for nyear, year in enumerate(years):
##            for nyear, year in enumerate([years[0]]):
##                print(year)   
#                
##                for nmonth, month in enumerate(list(np.arange(1,12+1))[0:2]):
#                for nmonth, month in enumerate(list(np.arange(1,12+1))):
#                    metdata_netcdf_fn = 'ERA5_' + str(year) + '-' + str(month).zfill(2) + '.nc'
#
#                    print('processing ' + vn + ' ', metdata_netcdf_fn)
#                    
#                    met_data = xr.open_dataset(metdata_netcdf_fp + metdata_netcdf_fn)
#        
#                    
#                    # Extract lat/lon indices only once
#                    if nyear + nmonth == 0:
#                        orog_data = xr.open_dataset(orog_data_fullfn)
#                        lat_idx = np.abs(lat_deg - orog_data['latitude'][:].values).argmin(axis=0)
#                        lon_idx = np.abs(lon_deg - orog_data['longitude'][:].values).argmin(axis=0)
#                        
#                        if ((abs(met_data.latitude.values[lat_idx] - orog_data.latitude.values[lat_idx]) > 0.01)
#                            or (abs(met_data.longitude.values[lon_idx] - orog_data.longitude.values[lon_idx]) > 0.01)):
#                            print('\n\nOROGRAPHY IS NOT LINED UP WITH CLIMATE DATA \n\n')
#                            print(met_data.latitude.values[lat_idx], met_data.longitude.values[lon_idx])
#                            print(orog_data.latitude.values[lat_idx], orog_data.longitude.values[lon_idx])
#                     
#                    if vn == 'rh':
#                        # relative humidity ('Arden Buck equation': approximation from Bogel modification)
#                        d2m = met_data.d2m[:,lat_idx,lon_idx].values
#                        t2m = met_data.t2m[:,lat_idx,lon_idx].values
#                        rh = (100 * np.exp((18.678*(d2m-273.15) / (257.14 + (d2m-273.15))) - 
#                                           ((18.678 - (t2m-273.15) / 234.5) * ((t2m-273.15) / (257.14 + (t2m-273.15))))))
#                        ds_roi_month = xr.Dataset({'rh': (['time'], rh)},
#                                                  coords={'time': met_data.time})
#                    elif vn == 'z':
#                        ds_roi_month = (orog_data[vn][0,lat_idx,lon_idx] / 9.80665).to_dataset()
#                    else:
#                        ds_roi_month = met_data[vn][:,lat_idx,lon_idx].to_dataset()
#                    
#                    # Concatenate dataset
#                    if nyear + nmonth == 0:
#                        ds_roi = ds_roi_month
#                    elif vn not in ['z']:
#                        ds_roi = xr.concat([ds_roi, ds_roi_month], 'time')
#                        
#            # Concatenate datasets                
#            if ds_all is None:
#                ds_all = ds_roi
#            else:
#                ds_all = ds_all.combine_first(ds_roi)
#            # Add attributes manually
#            if vn == 'tp':
#                ds_all.tp.attrs = {'units':'m', 'long_name':'Total precipitation'}
#            if vn == 'rh':
#                ds_all.rh.attrs = {'units':'%', 'long_name':'Relative humidity'}
#            elif vn == 'u10':
#                ds_all.u10.attrs = {'units': 'm s**-1', 'long_name': '10 metre U wind component'}
#            elif vn == 'v10':
#                ds_all.v10.attrs = {'units': 'm s**-1', 'long_name': '10 metre V wind component'}
#            elif vn == 'ssrd':
#                ds_all.ssrd.attrs = {'units':'J m**-2', 'long_name':'Surface solar radiation downwards', 
#                                   'standard_name': 'surface_downwelling_shortwave_flux_in_air'}
#            elif vn == 'strd':
#                ds_all.strd.attrs = {'units':'J m**-2', 'long_name':'Surface thermal radiation downwards'}
#            elif vn == 'z':
#                ds_all.z.attrs = {'units':'m a.s.l.', 'long_name':'Elevation', 'comment':'converted from geopotential'}
#               
#        
#        # Export array for each variable
#        if os.path.exists(input.metdata_fp) == False:
#            os.mkdir(input.metdata_fp)
#        ds_all.to_netcdf(output_metdata_fullfn)

#%%
#print('\nSHORTCUT FOR HMA WHICH IS ALREADY PROCESSED!\n')
## Meteorological data
##for nlatlon, latlon in enumerate([input.latlon_list[0]]):
#for nlatlon, latlon in enumerate(input.latlon_list):
#    print(nlatlon, latlon)
#    
#    lat_deg = latlon[0]
#    lon_deg = latlon[1]
#    
#    output_metdata_fullfn = output_metdata_sample.replace('XXXX', str(int(lat_deg*100)) + 'N-' + str(int(lon_deg*100)) 
#                                                          + 'E-')
#    
#    if os.path.exists(output_metdata_fullfn):
#        print('\n\nFILE ALREADY EXISTS\n\n')
#    else:
#        # ===== Combine meteorological data =====
#        ds_all = None
#        
#        for nvn, vn in enumerate(['t2m','tp', 'u10', 'v10', 'ssrd', 'strd', 'rh', 'z']):
##            print('\n',vn)
#        
#            metdata_netcdf_fp = input.main_directory + '/hma_data/'
#            metdata_netcdf_fn = 'HMA_ERA5-metdata_2000_2018-' + vn + '.nc'
#            orog_data_fullfn = input.main_directory + '/hma_data/HMA_ERA5-metdata_2000_2018-z.nc'
#
#            print('processing ' + vn + ' ', metdata_netcdf_fn)
#            
#            met_data = xr.open_dataset(metdata_netcdf_fp + metdata_netcdf_fn)
#
#            # Extract lat/lon indices only once
#            orog_data = xr.open_dataset(orog_data_fullfn)
#            lat_idx = np.abs(lat_deg - orog_data['latitude'][:].values).argmin(axis=0)
#            lon_idx = np.abs(lon_deg - orog_data['longitude'][:].values).argmin(axis=0)
#            
#            if ((abs(met_data.latitude.values[lat_idx] - orog_data.latitude.values[lat_idx]) > 0.01)
#                or (abs(met_data.longitude.values[lon_idx] - orog_data.longitude.values[lon_idx]) > 0.01)):
#                print('\n\nOROGRAPHY IS NOT LINED UP WITH CLIMATE DATA \n\n')
#                print(met_data.latitude.values[lat_idx], met_data.longitude.values[lon_idx])
#                print(orog_data.latitude.values[lat_idx], orog_data.longitude.values[lon_idx])
#                
#            if vn == 'z':
#                ds_roi_month = met_data[vn][lat_idx,lon_idx].to_dataset()
#            else:
#                ds_roi_month = met_data[vn][:,lat_idx,lon_idx].to_dataset()
#                        
#            # Concatenate datasets                
#            if ds_all is None:
#                ds_all = ds_roi_month
#            else:
#                ds_all = ds_all.combine_first(ds_roi_month)
#            # Add attributes manually
#            if vn == 'tp':
#                ds_all.tp.attrs = {'units':'m', 'long_name':'Total precipitation'}
#            if vn == 'rh':
#                ds_all.rh.attrs = {'units':'%', 'long_name':'Relative humidity'}
#            elif vn == 'u10':
#                ds_all.u10.attrs = {'units': 'm s**-1', 'long_name': '10 metre U wind component'}
#            elif vn == 'v10':
#                ds_all.v10.attrs = {'units': 'm s**-1', 'long_name': '10 metre V wind component'}
#            elif vn == 'ssrd':
#                ds_all.ssrd.attrs = {'units':'J m**-2', 'long_name':'Surface solar radiation downwards', 
#                                   'standard_name': 'surface_downwelling_shortwave_flux_in_air'}
#            elif vn == 'strd':
#                ds_all.strd.attrs = {'units':'J m**-2', 'long_name':'Surface thermal radiation downwards'}
#            elif vn == 'z':
#                ds_all.z.attrs = {'units':'m a.s.l.', 'long_name':'Elevation', 'comment':'converted from geopotential'}
#               
#        
#        # Export array for each variable
#        if os.path.exists(input.metdata_fp) == False:
#            os.mkdir(input.metdata_fp)
#        ds_all.to_netcdf(output_metdata_fullfn)    