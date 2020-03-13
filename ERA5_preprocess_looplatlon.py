#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 2019
@author: davidrounce
"""
# Built-in libraries
import argparse
import os
import pickle
from subprocess import call
# External libraries
import numpy as np
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
    parser.add_argument('-process_era5_hrly_data', action='store', type=int, default=0,
                        help='switch to process era5 hrly data')
    parser.add_argument('-process_unique_latlon_data', action='store', type=int, default=0,
                        help='switch to process data of unique lat/lons')
    parser.add_argument('-roi', action='store', type=str, default=None,
                        help='region of interest')
    parser.add_argument('-fromexternal', action='store', type=str, default='0',
                        help='option to download data from external hard drive')
    parser.add_argument('-debug', action='store', type=int, default=0,
                        help='Boolean for debugging to turn it on or off (default 0 is off')
    return parser 

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

    if args.process_era5_hrly_data == 1:
        # Find ERA5 filenames to process
        era5_fns = []
        for i in os.listdir(debris_prms.era5_hrly_fp):
            if i.startswith('ERA5_') and i.endswith('.nc'):
                date_str = i.split('_')[1]
                year_str = date_str.split('-')[0]
                startyear = debris_prms.roi_years[roi][0]
                endyear = debris_prms.roi_years[roi][1]
                if int(year_str) >= startyear and int(year_str) <= endyear:
                    era5_fns.append(i)
        era5_fns = sorted(era5_fns)
        
        for n, era5_fn in enumerate(era5_fns):
#        for n, era5_fn in enumerate([era5_fns[0]]):
            
            if debug:
                print(n, era5_fn)
            
            # Export subset
#            ds_out_fp = debris_prms.metdata_fp
            ds_out_fp = debris_prms.main_directory + '/../climate_data/' + roi + '/'
            ds_out_fn = roi + '-' + era5_fn
            if os.path.exists(ds_out_fp) == False:
                os.makedirs(ds_out_fp)
            
            if os.path.exists(ds_out_fp + ds_out_fn) == False:
                # Append arguments to call list
                call_list = ["python", "ERA5_preprocess.py"]         
                call_list.append('-era5_fn=' + era5_fn)
                call_list.append('-process_era5_hrly_data=1')
                call_list.append('-roi=' + roi)
                
                # Run script
                call(call_list)
         
#%%
    if args.process_unique_latlon_data == 1:
        
        output_metdata_fp = debris_prms.metdata_fp + '../' + roi + '/'
        metdata_fn_sample = (roi + '_ERA5-metdata-XXXX' + str(debris_prms.roi_years[roi][0]) + '_' + 
                             str(debris_prms.roi_years[roi][1]) + '.nc')
        if os.path.exists(output_metdata_fp) == False:
                os.makedirs(output_metdata_fp)
                
        option_fromexternal = int(args.fromexternal)
        
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
        
        
        # Process unique lat/lons 
        #  (best to get pickle from the debris_stats.ipynb script)
        if debris_prms.latlon_list_raw == 'all':
            with open(debris_prms.latlon_unique_fp + debris_prms.latlon_unique_dict[roi], 'rb') as f:
                latlon_list = pickle.load(f)
                
        
        for nlatlon, latlon_idx in enumerate(latlon_list):
#        for nlatlon, latlon_idx in enumerate([latlon_list[0]]):
            
            lat_deg = latlon_idx[0]
            lon_deg = latlon_idx[1]
            
            # Met data filename
            if lat_deg < 0:
                lat_str = 'S-'
            else:
                lat_str = 'N-'   
            output_metdata_fn = (metdata_fn_sample.replace('XXXX', str(int(abs(lat_deg)*100)) + lat_str + 
                                                           str(int(lon_deg*100)) + 'E-'))
            
#            print(output_metdata_fn, 'exists?', os.path.exists(output_metdata_fp + output_metdata_fn))
            
            
            if os.path.exists(output_metdata_fp + output_metdata_fn) == False:
                print('  ', output_metdata_fn)    
                # Append arguments to call list
                call_list = ["python", "ERA5_preprocess.py"]         
                call_list.append('-process_unique_latlon_data=1')
                call_list.append('-lat_deg=' + "%04.0f"%(lat_deg*100))
                call_list.append('-lon_deg=' + "%04.0f"%(lon_deg*100))
                call_list.append('-roi=' + roi)
                call_list.append('-fromexternal=' + str(option_fromexternal))
                
                # Run script
                call(call_list)
        

        
#%% ===== OLD SCRIPT PRE-01/13/2020 =====
#        # Select glaciers and get unique lat/lons
#        main_glac_rgi = debris_prms.selectglaciersrgitable(rgi_regionsO1=debris_prms.roi_rgidict[debris_prms.roi],
#                                                     rgi_regionsO2='all',
#                                                     rgi_glac_number='all')
#        # Convert longitude to be 0-360
#        main_glac_rgi['CenLon_360'] = main_glac_rgi['CenLon']
#        main_glac_rgi.loc[main_glac_rgi['CenLon_360'] < 0, 'CenLon_360'] = (
#                360 + main_glac_rgi.loc[main_glac_rgi['CenLon_360'] < 0,'CenLon_360'])
#        
#        # Group datasets by nearest lat/lon
#        ds = xr.open_dataset(ds_out_fp + era5_reg_fns[0])
#        #  argmin() finds the minimum distance between the glacier lat/lon and the GCM pixel
#        lat_nearidx = (np.abs(main_glac_rgi['CenLat'].values[:,np.newaxis] - 
#                              ds['latitude'][:].values).argmin(axis=1))
#        lon_nearidx = (np.abs(main_glac_rgi['CenLon_360'].values[:,np.newaxis] - 
#                              ds['longitude'][:].values).argmin(axis=1))
#        
#        latlon_nearidx = list(zip(lat_nearidx, lon_nearidx))
#        latlon_nearidx_unique = sorted(list(set(latlon_nearidx)))
#        
#        main_glac_rgi['latlon_nearidx'] = latlon_nearidx
#        latlon_unique_dict = dict(zip(latlon_nearidx_unique,np.arange(0,len(latlon_nearidx_unique))))
#        latlon_unique_dict_reversed = dict(zip(np.arange(0,len(latlon_nearidx_unique)),latlon_nearidx_unique))
#        main_glac_rgi['latlon_unique_no'] = main_glac_rgi['latlon_nearidx'].map(latlon_unique_dict)
#        
#        print('\nunique lat/lons:', len(latlon_nearidx_unique),'\n')
#        
#        
#        for nlatlon, latlon_idx in enumerate(latlon_nearidx_unique):
##        for nlatlon, latlon_idx in enumerate([latlon_nearidx_unique[0]]):
#            
#            main_glac_rgi_subset = main_glac_rgi[main_glac_rgi['latlon_unique_no'] == nlatlon]
#            
#            lat_idx, lon_idx = latlon_unique_dict_reversed[nlatlon]
#            
#            # Meteorological data
#            lat_deg = ds.latitude[lat_idx].values
#            lon_deg = ds.longitude[lon_idx].values
#            
#            if debug:
#                print(nlatlon, lat_idx, lat_deg, lon_deg)
#            
#            output_metdata_fp = debris_prms.metdata_fp + debris_prms.roi + '/'
#            output_metdata_fn = debris_prms.metdata_fn_sample.replace('XXXX', str(int(lat_deg*100)) + 'N-' + 
#                                                                str(int(lon_deg*100)) + 'E-')
#            if os.path.exists(output_metdata_fp) == False:
#                os.makedirs(output_metdata_fp)
#                
#            
#            if os.path.exists(output_metdata_fp + output_metdata_fn) == False:
#                # Append arguments to call list
#                call_list = ["python", "ERA5_preprocess.py"]         
#                call_list.append('-process_unique_latlon_data=1')
#                call_list.append('-lat_deg=' + "%04.0f"%(lat_deg*100))
#                call_list.append('-lon_deg=' + "%04.0f"%(lon_deg*100))
#                
#                # Run script
#                call(call_list)