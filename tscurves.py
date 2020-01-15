# -*- coding: utf-8 -*-
"""
Created on Sat Aug 11 10:17:13 2018

@author: David
"""

# Built-in libraries
import argparse
import collections
import multiprocessing
import os
import pickle
import time

# External libraries
#import rasterio
#import gdal
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import xarray as xr

# Local libraries
import globaldebris_input as input
from spc_split_lists import split_list

#%% ===== SCRIPT SPECIFIC INFORMATION =====
#eb_fp = input.main_directory + '/../output/exp3/'
#eb_fn = input.fn_prefix + 'YYYYN-' + 'XXXXE-' + input.date_start + '.nc'
#
#elev_cns2analyze = ['zmean']
#elev_cns2analyze = ['zmean', 'zstdlow', 'zstdhigh']

#%% ===== FUNCTIONS =====
def getparser():
    """
    Use argparse to add arguments from the command line

    Parameters
    ----------
    batchno (optional) : int
        batch number used to differentiate output on supercomputer
    batches (optional) : int
        total number of batches based on supercomputer
    num_simultaneous_processes (optional) : int
        number of cores to use in parallels
    option_parallels (optional) : int
        switch to use parallels or not
    debug (optional) : int
        Switch for turning debug printing on or off (default = 0 (off))

    Returns
    -------
    Object containing arguments and their respective values.
    """
    parser = argparse.ArgumentParser(description="run simulations from gcm list in parallel")
    # add arguments
    parser.add_argument('-batchno', action='store', type=int, default=0,
                        help='Batch number used to differentiate output on supercomputer')
    parser.add_argument('-batches', action='store', type=int, default=1,
                        help='Total number of batches (nodes) for supercomputer')
    parser.add_argument('-latlon_fn', action='store', type=str, default=None,
                        help='Filename containing list of lat/lon tuples for running batches on spc')
    parser.add_argument('-num_simultaneous_processes', action='store', type=int, default=4,
                        help='number of simultaneous processes (cores) to use')
    parser.add_argument('-option_parallels', action='store', type=int, default=1,
                        help='Switch to use or not use parallels (1 - use parallels, 0 - do not)')
    parser.add_argument('-option_ordered', action='store', type=int, default=1,
                        help='switch to keep lists ordered or not')
    parser.add_argument('-debug', action='store', type=int, default=0,
                        help='Boolean for debugging to turn it on or off (default 0 is off')
    return parser


def create_xrdataset_ts(ds, time_values):
    """
    Create empty xarray dataset that will be used to record surface temperature from simulation runs.

    Parameters
    ----------
    ds : xarray dataset
        dataframe containing energy balance model runs

    Returns
    -------
    output_ds_all : xarray Dataset
        empty xarray dataset that contains variables and attributes to be filled in by simulation runs
    encoding : dictionary
        encoding used with exporting xarray dataset to netcdf
    """
    # Create empty datasets for each variable and merge them
    # Variable coordinates dictionary
    output_coords_dict = {
            'ts': collections.OrderedDict(
                    [('hd_cm', ds.hd_cm.values), ('time', time_values), ('stats', ds.stats.values), 
                     ('elev_cns', ds.elev_cns.values)]),
            'dsnow': collections.OrderedDict(
                    [('hd_cm', ds.hd_cm.values), ('time', time_values), ('stats', ds.stats.values), 
                     ('elev_cns', ds.elev_cns.values)]),
            'elev': collections.OrderedDict(
                    [('elev_cns', ds.elev_cns.values)])
            }
    # Attributes dictionary
    output_attrs_dict = {
            'hd_cm': {
                    'long_name': 'debris thickness in centimeters',
                    'comment': 'cm so values are integers'},
            'elev': {
                    'long_name': 'elevation',
                    'units': 'm a.s.l.',
                    'comment': 'elevation associated with the elevation column name (elev_cns)'},
            'stats': {
                    'long_name': 'variable statistics',
                    'comment': str(input.k_random.shape[0]) + ' simulations; % refers to percentiles'},
            'elev_cns': {
                    'long_name': 'elevation column names',
                    'comment': 'elevations used to run the simulations'},
            'ts': {
                    'long_name': 'debris surface temperature',
                    'units': 'K',
                    'temporal_resolution': 'hourly'},
            'dsnow': {
                    'long_name': 'snow depth',
                    'units': 'm',
                    'temporal_resolution': 'hourly'},
            'time': {
                    'long_name': 'time of satellite acquisition',
#                    'units': '-',
                    'temporal_resolution': 'daily'}
            }
    # Add variables to empty dataset and merge together
    count_vn = 0
    encoding = {}
    noencoding_vn = ['stats', 'hd_cm', 'elev_cns', 'elev']
    for vn in ['ts', 'dsnow', 'elev']:
        count_vn += 1
        empty_holder = np.zeros([len(output_coords_dict[vn][i]) for i in list(output_coords_dict[vn].keys())])
        output_ds = xr.Dataset({vn: (list(output_coords_dict[vn].keys()), empty_holder)},
                               coords=output_coords_dict[vn])
        if vn == 'elev':
            output_ds['elev'].values = ds.elev.values
        # Merge datasets of stats into one output
        if count_vn == 1:
            output_ds_all = output_ds
        else:
            output_ds_all = xr.merge((output_ds_all, output_ds))

    # Add attributes
    for vn in ['ts', 'dsnow', 'hd_cm', 'stats', 'elev_cns', 'elev', 'time']:
        try:
            output_ds_all[vn].attrs = output_attrs_dict[vn]
        except:
            pass
        # Encoding (specify _FillValue, offsets, etc.)
        if vn not in noencoding_vn:
            encoding[vn] = {'_FillValue': False}

    return output_ds_all, encoding


def main(list_packed_vars):
    """
    Model simulation

    Parameters
    ----------
    list_packed_vars : list
        list of packed variables that enable the use of parallels

    Returns
    -------
    netcdf files of the simulation output (specific output is dependent on the output option)
    """
    # Unpack variables
    count = list_packed_vars[0]
    latlon_list = list_packed_vars[1]
    
#    if debug:
#        print(count, latlon_list)
        
    # Surface temperature information (year, day of year, hour)
    ts_info_fullfn = input.ts_fp + input.roi + '_debris_tsinfo.nc'
    ds_ts_info = xr.open_dataset(ts_info_fullfn, decode_times=False)
    
    for nlatlon, latlon in enumerate(latlon_list):
#        if debug:
        print(nlatlon, latlon)
        
        lat_deg = latlon[0]
        lon_deg = latlon[1]
        
        # ===== Debris Thickness vs. Surface Lowering =====
        # stats column index (0=mean)
        stats_idx = 0
        
        # Filename
        plot_str = str(int(lat_deg*100)) + 'N-' + str(int(lon_deg*100)) + 'E'
        ds_tscurve_fn = input.output_ts_fn_sample.replace('XXXX', plot_str)
        
        # Melt model output fn                
        if lat_deg < 0:
            lat_str = 'S-'
        else:
            lat_str = 'N-'
        ds_meltmodel_fn = (input.fn_prefix + str(int(abs(lat_deg*100))) + lat_str + str(int(lon_deg*100)) + 'E-' + 
                           input.date_start + '.nc')
        
        if ((os.path.exists(input.tscurve_fp + ds_tscurve_fn) == False) and 
            (os.path.exists(input.eb_fp + ds_meltmodel_fn) == True)):
            
            # Dataset from energy balance modeling
            ds = xr.open_dataset(input.eb_fp + ds_meltmodel_fn)
            
            # Time information of surface temperature
            lat_idx = np.abs(lat_deg - ds_ts_info['latitude'][:].values).argmin(axis=0)
            lon_idx = np.abs(lon_deg - ds_ts_info['longitude'][:].values).argmin(axis=0)
            
            ts_year = np.round(ds_ts_info['year_mean'][lat_idx,lon_idx].values,0)
            ts_doy = np.round(ds_ts_info['doy_med'][lat_idx,lon_idx].values,0)
            ts_hr = ds_ts_info['dayfrac_mean'][lat_idx,lon_idx].values
            
#            if debug:            
#                print('ts_hr:', ts_hr, 'ts_doy:', ts_doy, 'ts_year', ts_year)
            
            if ts_year > 0 and ts_doy > 0:
                ts_str = str(int(ts_year)) + '-' + str(int(ts_doy))
                ts_date_pd = pd.to_datetime(pd.DataFrame(np.array([ts_str]))[0], format='%Y-%j')
                ts_date = ts_date_pd.values[0] + np.timedelta64(int(ts_hr),'h')
                
                # Index with model results      
                time_pd = pd.to_datetime(ds.time.values)  
                time_idx = np.where(ts_date == time_pd)[0][0]
                # index one month before and after to get statistics
                time_idx_all = np.arange(time_idx - 30*24, time_idx + 31*24, 24)
                time_all_interpolated = time_pd[time_idx_all] + np.timedelta64(int((ts_hr%1)*60),'m')
                
                # Output dataset
                ds_ts, encoding = create_xrdataset_ts(ds, time_all_interpolated)
                
                for nelev, elev_cn in enumerate(input.elev_cns):
#                    if debug:
#                        print(nelev, elev_cn)
                        
                    # Select data
                    #  each row is a debris thickness
                    ts_t1 = ds['ts'][:,time_idx_all,stats_idx,nelev].values
                    ts_t2 = ds['ts'][:,time_idx_all+1,stats_idx,nelev].values
                    ts_data = ts_t1 + ts_hr%1 * (ts_t2 - ts_t1)
                    
                    dsnow_t1 = ds['snow_depth'][:,time_idx_all,stats_idx,nelev].values
                    dsnow_t2 = ds['snow_depth'][:,time_idx_all+1,stats_idx,nelev].values
                    dsnow_data = dsnow_t1 + ts_hr%1 * (dsnow_t2 - dsnow_t1)
                    
                    ds_ts['ts'][:,:,stats_idx,nelev] = ts_data
                    ds_ts['dsnow'][:,:,stats_idx,nelev] = dsnow_data
                        
                # Export netcdf
                if os.path.exists(input.tscurve_fp) == False:
                    os.makedirs(input.tscurve_fp)
                ds_ts.to_netcdf(input.tscurve_fp + ds_tscurve_fn)

    if debug:
        return ds_ts
            

#%%
    
if __name__ == '__main__':
    time_start = time.time()
    parser = getparser()
    args = parser.parse_args()

    if args.debug == 1:
        debug = True
    else:
        debug = False

    time_start = time.time()
    
    # RGI glacier number
    if args.latlon_fn is not None:
        with open(args.latlon_fn, 'rb') as f:
            latlon_list = pickle.load(f)
    else:
        latlon_list = input.latlon_list
    
    # Number of cores for parallel processing
    if args.option_parallels != 0:
        num_cores = int(np.min([len(latlon_list), args.num_simultaneous_processes]))
    else:
        num_cores = 1

    # Glacier number lists to pass for parallel processing
    latlon_lsts = split_list(latlon_list, n=num_cores, option_ordered=args.option_ordered)

    # Pack variables for multiprocessing
    list_packed_vars = []
    for count, latlon_lst in enumerate(latlon_lsts):
        list_packed_vars.append([count, latlon_lst])

    # Parallel processing
    if args.option_parallels != 0:
        print('Processing in parallel with ' + str(args.num_simultaneous_processes) + ' cores...')
        with multiprocessing.Pool(args.num_simultaneous_processes) as p:
            p.map(main,list_packed_vars)
    # If not in parallel, then only should be one loop
    else:
        # Loop through the chunks and export bias adjustments
        for n in range(len(list_packed_vars)):
            if debug and num_cores == 1:
                ds_ts = main(list_packed_vars[n])
            else:
                main(list_packed_vars[n])
     
    print('\nProcessing time of :',time.time()-time_start, 's')
    
    
    #%%
    
