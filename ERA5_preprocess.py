#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Model input for intercomparison experiment
"""
# Built-in libraries
import argparse
import os
# External libraries
import numpy as np
#import pandas as pd
import xarray as xr
# Local libraries
import globaldebris_input as input


def getparser():
    """
    Use argparse to add arguments from the command line

    Parameters
    ----------
    era5_fn (optional) : int
        batch number used to differentiate output on supercomputer
    debug (optional) : int
        Switch for turning debug printing on or off (default = 0 (off))

    Returns
    -------
    Object containing arguments and their respective values.
    """
    parser = argparse.ArgumentParser(description="run simulations from gcm list in parallel")
    # add arguments
    parser.add_argument('-process_era5_hrly_data', action='store', type=str, default=0,
                        help='switch to process era5 hrly data to regional subsets')
    parser.add_argument('-era5_fn', action='store', type=str, default=None,
                        help='ERA5 filename to process')
    parser.add_argument('-process_unique_latlon_data', action='store', type=str, default=0,
                        help='switch to process data of unique lat/lons')
    parser.add_argument('-lat_deg', action='store', type=str, default=0,
                        help='latitude * 100 in degrees')
    parser.add_argument('-lon_deg', action='store', type=str, default=0,
                        help='longitude * 100 in degrees')
    parser.add_argument('-debug', action='store', type=int, default=0,
                        help='Boolean for debugging to turn it on or off (default 0 is off')
    return parser 

# Simulation data
#orog_data_fullfn = '/Users/davidrounce/Documents/Dave_Rounce/HiMAT/Climate_data/ERA5/ERA5_geopotential_monthly.nc'
#timezone = 0
    
if __name__ == '__main__':
    parser = getparser()
    args = parser.parse_args()
    
    if args.debug == 1:
        debug = True
    else:
        debug = False


    #%% ===== SUBSET ERA5 HOURLY DATA TO REGIONAL EXTENTS TO REDUCE FILE SIZES TO MANAGEABLE SIZE =====
    if args.process_era5_hrly_data == '1':
        if args.era5_fn is None:
            era5_fns = []
            for i in os.listdir(input.era5_hrly_fp):
                if i.startswith('ERA5_') and i.endswith('.nc'):
                    era5_fns.append(i)
            era5_fns = sorted(era5_fns)
        else:
            era5_fns = [args.era5_fn]
    
        ds_elev = xr.open_dataset(input.metdata_fp +input.metdata_elev_fn)
        
        lat_N = input.roi_latlon_dict[input.roi][0]
        lat_S = input.roi_latlon_dict[input.roi][1]
        lon_E = input.roi_latlon_dict[input.roi][2]
        lon_W = input.roi_latlon_dict[input.roi][3]
        
        lat_N_idx = np.abs(lat_N - ds_elev['latitude'].values).argmin()
        lat_S_idx = np.abs(lat_S - ds_elev['latitude'].values).argmin()
        lon_E_idx = np.abs(lon_E - ds_elev['longitude'].values).argmin()
        lon_W_idx = np.abs(lon_W - ds_elev['longitude'].values).argmin()
        
        if debug:
            print(lat_N_idx, ds_elev['latitude'][lat_N_idx].values, lat_N)
            print(lat_S_idx, ds_elev['latitude'][lat_S_idx].values, lat_S)
            print(lon_E_idx, ds_elev['longitude'][lon_E_idx].values, lon_E)
            print(lon_W_idx, ds_elev['longitude'][lon_W_idx].values, lon_W)
    
        for n, era5_fn in enumerate(era5_fns):
            if debug:
                print(n, era5_fn)
            
            ds = xr.open_dataset(input.era5_hrly_fp + era5_fn)
            ds_out = ds.sel(latitude=slice(ds_elev['latitude'][lat_N_idx].values,ds_elev['latitude'][lat_S_idx].values), 
                            longitude=slice(ds_elev['longitude'][lon_W_idx].values,ds_elev['longitude'][lon_E_idx].values))
            
            # Export subset
            ds_out_fp = input.metdata_fp + input.roi + '/'
            ds_out_fn = input.roi + '-' + era5_fn
            if os.path.exists(ds_out_fp) == False:
                os.makedirs(ds_out_fp)
            
            if os.path.exists(ds_out_fp + ds_out_fn) == False:
                ds_out.to_netcdf(ds_out_fp + ds_out_fn)


    #%% ===== EXTRACT DATA FOR UNIQUE LAT/LONS =====
    if args.process_unique_latlon_data == '1':
        lat_deg = int(args.lat_deg) / 100
        lon_deg = int(args.lon_deg) / 100
        
        print(lat_deg, lon_deg)
        
        ds_out_fp = input.metdata_fp + input.roi + '/'
        era5_reg_fns = []
        for i in os.listdir(ds_out_fp):
            if i.startswith(input.roi + '-ERA5_') and i.endswith('.nc'):
                era5_reg_fns.append(i)
        era5_reg_fns = sorted(era5_reg_fns)
        
        output_metdata_fp = input.metdata_fp + input.roi + '/'
        output_metdata_fn = input.metdata_fn_sample.replace('XXXX', str(int(lat_deg*100)) + 'N-' + 
                                                            str(int(lon_deg*100)) + 'E-')
        if os.path.exists(output_metdata_fp) == False:
            os.makedirs(output_metdata_fp)
            
        
        if os.path.exists(output_metdata_fp + output_metdata_fn) == False:
            # ===== Combine meteorological data =====
            ds_all = None
            years = list(np.arange(int(input.startyear), int(input.endyear)+1))
            
            ds_elev = xr.open_dataset(input.metdata_fp +input.metdata_elev_fn)
            
            for vn in ['t2m','tp', 'u10', 'v10', 'ssrd', 'strd', 'rh', 'z']:
#            for vn in ['t2m', 'z']:
#                print('\n',vn)
            
                for nyear, year in enumerate(years):
#                for nyear, year in enumerate([years[0]]):
#                    print(year)   

                    for nmonth, month in enumerate(list(np.arange(1,12+1))):
#                    for nmonth, month in enumerate(list(np.arange(1,12+1))[0:2]):
                        metdata_netcdf_fn = ds_out_fn = input.roi + '-' + 'ERA5_' + str(year) + '-' + str(month).zfill(2) + '.nc'
                        
                        met_data = xr.open_dataset(output_metdata_fp + metdata_netcdf_fn)
            
                        
                        # Extract lat/lon indices only once
                        if nyear + nmonth == 0:
                            lat_idx_z = np.abs(lat_deg - ds_elev['latitude'][:].values).argmin(axis=0)
                            lon_idx_z = np.abs(lon_deg - ds_elev['longitude'][:].values).argmin(axis=0)
                            
                            lat_idx = np.abs(lat_deg - met_data['latitude'].values).argmin()
                            lon_idx = np.abs(lon_deg - met_data['longitude'].values).argmin()

                        if vn == 'rh':
                            # relative humidity ('Arden Buck equation': approximation from Bogel modification)
                            d2m = met_data.d2m[:,lat_idx,lon_idx].values
                            t2m = met_data.t2m[:,lat_idx,lon_idx].values
                            rh = (100 * np.exp((18.678*(d2m-273.15) / (257.14 + (d2m-273.15))) - 
                                               ((18.678 - (t2m-273.15) / 234.5) * ((t2m-273.15) / (257.14 + (t2m-273.15))))))
                            ds_roi_month = xr.Dataset({'rh': (['time'], rh)},
                                                      coords={'time': met_data.time})
                        elif vn == 'z':
                            ds_roi_month = ds_elev[vn][lat_idx_z,lon_idx_z].to_dataset()
                        else:
                            ds_roi_month = met_data[vn][:,lat_idx,lon_idx].to_dataset()
                        
                        # Concatenate dataset
                        if nyear + nmonth == 0:
                            ds_roi = ds_roi_month
                        elif vn not in ['z']:
                            ds_roi = xr.concat([ds_roi, ds_roi_month], 'time')
                            
                # Concatenate datasets                
                if ds_all is None:
                    ds_all = ds_roi
                else:
                    ds_all = ds_all.combine_first(ds_roi)
                # Add attributes manually
                if vn == 'tp':
                    ds_all.tp.attrs = {'units':'m', 'long_name':'Total precipitation'}
                if vn == 'rh':
                    ds_all.rh.attrs = {'units':'%', 'long_name':'Relative humidity'}
                elif vn == 'u10':
                    ds_all.u10.attrs = {'units': 'm s**-1', 'long_name': '10 metre U wind component'}
                elif vn == 'v10':
                    ds_all.v10.attrs = {'units': 'm s**-1', 'long_name': '10 metre V wind component'}
                elif vn == 'ssrd':
                    ds_all.ssrd.attrs = {'units':'J m**-2', 'long_name':'Surface solar radiation downwards', 
                                       'standard_name': 'surface_downwelling_shortwave_flux_in_air'}
                elif vn == 'strd':
                    ds_all.strd.attrs = {'units':'J m**-2', 'long_name':'Surface thermal radiation downwards'}
                elif vn == 'z':
                    ds_all.z.attrs = {'units':'m a.s.l.', 'long_name':'Elevation', 'comment':'converted from geopotential'}
                  
            # Export array for each variable
            ds_all.to_netcdf(output_metdata_fp + output_metdata_fn)    
        
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