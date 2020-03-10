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
#import globaldebris_input as input
import debrisglobal.globaldebris_input as debris_prms


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
    parser.add_argument('-roi', action='store', type=str, default=None,
                        help='region of interest')
    parser.add_argument('-fromexternal', action='store', type=str, default='0',
                        help='option to download data from external hard drive')
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
        
    if args.roi is None:
        roi = debris_prms.roi
    else:
        roi = args.roi
        
    option_fromexternal = int(args.fromexternal)

    #%% ===== SUBSET ERA5 HOURLY DATA TO REGIONAL EXTENTS TO REDUCE FILE SIZES TO MANAGEABLE SIZE =====
    if args.process_era5_hrly_data == '1':
        if args.era5_fn is None:
            era5_fns = []
            for i in os.listdir(debris_prms.era5_hrly_fp):
                if i.startswith('ERA5_') and i.endswith('.nc'):
                    era5_fns.append(i)
            era5_fns = sorted(era5_fns)
        else:
            era5_fns = [args.era5_fn]
    
        ds_elev = xr.open_dataset(debris_prms.metdata_fp + '../' + debris_prms.metdata_elev_fn)
        
        lat_N = debris_prms.roi_latlon_dict[roi][0]
        lat_S = debris_prms.roi_latlon_dict[roi][1]
        lon_E = debris_prms.roi_latlon_dict[roi][2]
        lon_W = debris_prms.roi_latlon_dict[roi][3]
        
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
            
            ds = xr.open_dataset(debris_prms.era5_hrly_fp + era5_fn)
            ds_out = ds.sel(latitude=slice(ds_elev['latitude'][lat_N_idx].values,ds_elev['latitude'][lat_S_idx].values), 
                            longitude=slice(ds_elev['longitude'][lon_W_idx].values,ds_elev['longitude'][lon_E_idx].values))
            
            # Export subset
            ds_out_fp = debris_prms.metdata_fp + '../' + roi + '/'
            ds_out_fn = roi + '-' + era5_fn
            if os.path.exists(ds_out_fp) == False:
                os.makedirs(ds_out_fp)
            
            if os.path.exists(ds_out_fp + ds_out_fn) == False:
                ds_out.to_netcdf(ds_out_fp + ds_out_fn)

    
    #%% ===== EXTRACT DATA FOR UNIQUE LAT/LONS =====
    if args.process_unique_latlon_data == '1':
        lat_deg = int(args.lat_deg) / 100
        lon_deg = int(args.lon_deg) / 100
        
        print(lat_deg, lon_deg)

        output_metdata_fp = debris_prms.metdata_fp + '../' + roi + '/'
        metdata_fn_sample = (roi + '_ERA5-metdata-XXXX' + str(debris_prms.roi_years[roi][0]) + '_' + 
                             str(debris_prms.roi_years[roi][1]) + '.nc')
        if os.path.exists(output_metdata_fp) == False:
            os.makedirs(output_metdata_fp)
        
        if option_fromexternal == 1:
            era5_fp = '/Volumes/LaCie_Raid/ERA5_hrly/'
        else:
            era5_fp = debris_prms.metdata_fp + '../' + roi + '/'
        era5_reg_fns = []
        for i in os.listdir(era5_fp):
            if option_fromexternal == 1:
                if i.startswith('ERA5_') and i.endswith('.nc'):
                    era5_reg_fns.append(i)
            else:
                if i.startswith(roi + '-ERA5_') and i.endswith('.nc'):
                    era5_reg_fns.append(i)
        era5_reg_fns = sorted(era5_reg_fns)
        
        # Met data filename
        if lat_deg < 0:
            lat_str = 'S-'
        else:
            lat_str = 'N-'   
        output_metdata_fn = (metdata_fn_sample.replace('XXXX', str(int(abs(lat_deg)*100)) + lat_str + 
                                                       str(int(lon_deg*100)) + 'E-'))
        
        lat_N = debris_prms.roi_latlon_dict[roi][0]
        lat_S = debris_prms.roi_latlon_dict[roi][1]
        lon_E = debris_prms.roi_latlon_dict[roi][2]
        lon_W = debris_prms.roi_latlon_dict[roi][3]
        
        if (os.path.exists(output_metdata_fp + output_metdata_fn) == False and 
            lat_deg <= lat_N and lat_deg >= lat_S and lon_deg >= lon_W and lon_deg <= lon_E):
            # ===== Combine meteorological data =====
            ds_all = None
            years = list(np.arange(int(debris_prms.roi_years[roi][0]), int(debris_prms.roi_years[roi][1])+1))
            
            ds_elev = xr.open_dataset(debris_prms.metdata_fp + '../' + debris_prms.metdata_elev_fn)
            
            for nyear, year in enumerate(years):
                print(year)

                for nmonth, month in enumerate(list(np.arange(1,12+1))):
                    print(year, month)
                    
                    if option_fromexternal == 1:
                        metdata_netcdf_fn = 'ERA5_' + str(year) + '-' + str(month).zfill(2) + '.nc'
                    else:
                        metdata_netcdf_fn = roi + '-' + 'ERA5_' + str(year) + '-' + str(month).zfill(2) + '.nc'
                    
                    met_data = xr.open_dataset(era5_fp + metdata_netcdf_fn)
                    
                    # Extract lat/lon indices only once
                    if nyear + nmonth == 0:
                        lat_idx_z = np.abs(lat_deg - ds_elev['latitude'][:].values).argmin(axis=0)
                        lon_idx_z = np.abs(lon_deg - ds_elev['longitude'][:].values).argmin(axis=0)
                        
                        lat_idx = np.abs(lat_deg - met_data['latitude'].values).argmin()
                        lon_idx = np.abs(lon_deg - met_data['longitude'].values).argmin()
                        
                    # Extract lat/lon data
                    met_data_latlon = met_data[dict(longitude=lon_idx, latitude=lat_idx)]
                    
                    # Add relative humidity
                    #   relative humidity ('Arden Buck equation': approximation from Bogel modification)
                    d2m = met_data_latlon.d2m.values
                    t2m = met_data_latlon.t2m.values
                    rh = (100 * np.exp((18.678*(d2m-273.15) / (257.14 + (d2m-273.15))) - 
                                       ((18.678 - (t2m-273.15) / 234.5) * ((t2m-273.15) / (257.14 + (t2m-273.15))))))
                    met_data_rh = xr.Dataset({'rh': (['time'], rh)},
                                              coords={'time': met_data.time,
                                                      'longitude': met_data_latlon.longitude,
                                                      'latitude': met_data_latlon.latitude})
                    met_data_rh.rh.attrs = {'units':'%', 'long_name':'Relative humidity'}
                    met_data_latlon = xr.merge([met_data_latlon, met_data_rh])
                    met_data_latlon = met_data_latlon.drop(['d2m'])
                    
                    # Add elevation
                    ds_z = ds_elev['z'][lat_idx_z,lon_idx_z].to_dataset()
                    ds_z.z.attrs = {'units':'m a.s.l.', 'long_name':'Elevation', 'comment':'converted from geopotential'}
                    met_data_latlon = xr.merge([met_data_latlon, ds_z])

                    try:
                        met_data.close()
                    except:
                        continue
                            
                    # Concatenate datasets                
                    if ds_all is None:
                        ds_all = met_data_latlon
                    else:
                        ds_all = ds_all.combine_first(met_data_latlon) 
                  
            # Export array for each variable
            print('exporting...' + output_metdata_fn)
            ds_all.to_netcdf(output_metdata_fp + output_metdata_fn)  
            
        
#%%
#print('\nSHORTCUT FOR HMA WHICH IS ALREADY PROCESSED!\n')
## Meteorological data
##for nlatlon, latlon in enumerate([debris_prms.latlon_list[0]]):
#for nlatlon, latlon in enumerate(debris_prms.latlon_list):
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
#            metdata_netcdf_fp = debris_prms.main_directory + '/hma_data/'
#            metdata_netcdf_fn = 'HMA_ERA5-metdata_2000_2018-' + vn + '.nc'
#            orog_data_fullfn = debris_prms.main_directory + '/hma_data/HMA_ERA5-metdata_2000_2018-z.nc'
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
#        if os.path.exists(debris_prms.metdata_fp) == False:
#            os.mkdir(debris_prms.metdata_fp)
#        ds_all.to_netcdf(output_metdata_fullfn)    