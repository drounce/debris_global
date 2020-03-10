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
import debrisglobal.globaldebris_input as debris_prms
from spc_split_lists import split_list

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
    parser.add_argument('-plotfigs', action='store', type=int, default=1,
                        help='Boolean for plotting figures or not (default 1 is to plot')
    return parser


# fit curve
# NOTE: two ways of writing the 2nd order reaction rate equations
#  1 / A = 1 / A0 + kt  (here we replaced t with h)
#def melt_fromdebris_func(h, a, k):
#    """ estimate melt from debris thickness (h is debris thickness, a and k are coefficients) """
#    return a / (1 + 2 * k * a * h)
#def debris_frommelt_func(b, a, k):
#    """ estimate debris thickness from melt (b is melt, a and k are coefficients) """
#    return (a - b) / (2*k*a*b)
def melt_fromdebris_func(h, a, k):
    """ Second order reaction rate equation used to estimate melt from debris thickness
    The standard form is 1/A = 1/A0 + kt
      1/b = 1/a + k*h  derived from the standard form: 
    where b is melt, h is debris thickness, and k and a are constants. """
    return 1 / (1 / a + k * h)


def debris_frommelt_func(b, a, k):
    """ estimate debris thickness from melt (b is melt, a and k are coefficients) """
    return 1 / k * (1 / b - 1 / a)


#def create_xrdataset_ostrem(ds):
#    """
#    Create empty xarray dataset that will be used to record melt data from simulation runs.
#
#    Parameters
#    ----------
#    ds : xarray dataset
#        dataframe containing energy balance model runs
#
#    Returns
#    -------
#    output_ds_all : xarray Dataset
#        empty xarray dataset that contains variables and attributes to be filled in by simulation runs
#    encoding : dictionary
#        encoding used with exporting xarray dataset to netcdf
#    """
#    # Create empty datasets for each variable and merge them
#    # Variable coordinates dictionary
#    output_coords_dict = {
#            'melt_mwea': collections.OrderedDict(
#                    [('hd_cm', ds.hd_cm.values), ('stats', ds.stats.values), ('elev_cns', ds.elev_cns.values)]),
#            'b0': collections.OrderedDict(
#                    [('stats', ds.stats.values), ('elev_cns', ds.elev_cns.values)]),
#            'k': collections.OrderedDict(
#                    [('stats', ds.stats.values), ('elev_cns', ds.elev_cns.values)]),
#            'elev': collections.OrderedDict(
#                    [('elev_cns', ds.elev_cns.values)])
#            }
#    # Attributes dictionary
#    output_attrs_dict = {
#            'hd_cm': {
#                    'long_name': 'debris thickness in centimeters',
#                    'comment': 'cm so values are integers'},
#            'elev': {
#                    'long_name': 'elevation',
#                    'units': 'm a.s.l.',
#                    'comment': 'elevation associated with the elevation column name (elev_cns)'},
#            'stats': {
#                    'long_name': 'variable statistics',
#                    'comment': str(debris_prms.k_random.shape[0]) + ' simulations; % refers to percentiles'},
#            'elev_cns': {
#                    'long_name': 'elevation column names',
#                    'comment': 'elevations used to run the simulations'},
#            'melt_mwea': {
#                    'long_name': 'annual sub-debris glacier melt',
#                    'units': 'meters water equivalent per year',
#                    'temporal_resolution': 'annual',
#                    'start_date': debris_prms.start_date,
#                    'end_date': debris_prms.end_date},
#            'b0': {
#                    'long_name': 'second order reaction rate initial concentration',
#                    'units': 'm w.e.'},
#            'k': {
#                    'long_name': 'second order reaction rate initial concentration',
#                    'units': 'm w.e.'}
#            }
#    # Add variables to empty dataset and merge together
#    count_vn = 0
#    encoding = {}
#    noencoding_vn = ['stats', 'hd_cm', 'elev_cns', 'elev']
#    for vn in ['melt_mwea', 'b0', 'k', 'elev']:
#        count_vn += 1
#        empty_holder = np.zeros([len(output_coords_dict[vn][i]) for i in list(output_coords_dict[vn].keys())])
#        output_ds = xr.Dataset({vn: (list(output_coords_dict[vn].keys()), empty_holder)},
#                               coords=output_coords_dict[vn])
#        if vn == 'elev':
#            output_ds['elev'].values = ds.elev.values
#        # Merge datasets of stats into one output
#        if count_vn == 1:
#            output_ds_all = output_ds
#        else:
#            output_ds_all = xr.merge((output_ds_all, output_ds))
#
#    # Add attributes
#    for vn in ['melt', 'hd_cm', 'stats', 'elev_cns', 'elev', 'b0', 'k']:
#        try:
#            output_ds_all[vn].attrs = output_attrs_dict[vn]
#        except:
#            pass
#        # Encoding (specify _FillValue, offsets, etc.)
#        if vn not in noencoding_vn:
#            encoding[vn] = {'_FillValue': False}
#
#    return output_ds_all, encoding
def export_ds_daily_melt(ds):
    """
    Create empty xarray dataset that will be used to record melt data from simulation runs.

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
    # Extract time values
    time_daily = pd.to_datetime(ds.time.values[0::24])
    
    # Melt daily
    ds_shape = ds.melt.values.shape
    melt_daily = ds.melt.values.reshape((ds_shape[0], int(ds_shape[1] / 24), 24, ds_shape[2], ds_shape[3])).sum(axis=2)
    
    # Variable coordinates dictionary
    output_coords_dict = {
            'melt': collections.OrderedDict(
                    [('hd_cm', ds.hd_cm.values), ('time', time_daily), 
                     ('stats', ds.stats.values), ('elev_cns', ds.elev_cns.values)]),
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
                    'comment': str(debris_prms.k_random.shape[0]) + ' simulations; % refers to percentiles'},
            'elev_cns': {
                    'long_name': 'elevation column names',
                    'comment': 'elevations used to run the simulations'},
            'melt': {
                    'long_name': 'glacier melt',
                    'units': 'meters water equivalent',
                    'temporal_resolution': 'daily'}
            }
    # Add variables to empty dataset and merge together
    count_vn = 0
    encoding = {}
    noencoding_vn = ['stats', 'hd_cm', 'elev_cns', 'elev']
    for vn in ['melt', 'elev']:
        count_vn += 1
        empty_holder = np.zeros([len(output_coords_dict[vn][i]) for i in list(output_coords_dict[vn].keys())])
        output_ds = xr.Dataset({vn: (list(output_coords_dict[vn].keys()), empty_holder)},
                               coords=output_coords_dict[vn])
        if vn == 'elev':
            output_ds['elev'].values = ds.elev.values
        elif vn == 'melt':
            output_ds['melt'].values = melt_daily
        # Merge datasets of stats into one output
        if count_vn == 1:
            output_ds_all = output_ds
        else:
            output_ds_all = xr.merge((output_ds_all, output_ds))

    # Add attributes
    for vn in ['melt', 'hd_cm', 'stats', 'elev_cns', 'elev', 'b0', 'k']:
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
    
    if debug:
        print(count, latlon_list)
    
    for nlatlon, latlon in enumerate(latlon_list):

        print(nlatlon, latlon)
        
        lat_deg = latlon[0]
        lon_deg = latlon[1]
        
        # ===== Debris Thickness vs. Surface Lowering =====        
        # Filename
        if os.path.exists(debris_prms.ostrem_fp) == False:
            os.makedirs(debris_prms.ostrem_fp)
        
        # Melt model output fn
        if lat_deg < 0:
            lat_str = 'S-'
        else:
            lat_str = 'N-'
        latlon_str = str(int(abs(lat_deg*100))) + lat_str + str(int(lon_deg*100)) + 'E-'
        # Raw meltmodel output filename
        ds_meltmodel_fn = debris_prms.fn_prefix + latlon_str + debris_prms.date_start + '.nc'
        # Processed "ostrem" filename, although this is really only daily data to enable each glacier to choose correct
        #  dates over which to sum the melt consistent with the DEM differencing
        ds_ostrem_fn = debris_prms.ostrem_fn_sample.replace('XXXX', latlon_str)
        
        
        if os.path.exists(debris_prms.ostrem_fp + ds_ostrem_fn) == False:
            
            # Debris thickness vs. melt dataset from energy balance modeling
            ds = xr.open_dataset(debris_prms.eb_fp + ds_meltmodel_fn)
            
            ds_ostrem, encoding = export_ds_daily_melt(ds)
            # Export netcdf
            ds_ostrem.to_netcdf(debris_prms.ostrem_fp + ds_ostrem_fn)
            

    if debug:
        return ds_ostrem  
            

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
        latlon_list = debris_prms.latlon_list    

        #%%
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
                ds_ostrem = main(list_packed_vars[n])
            else:
                main(list_packed_vars[n])
                
    print('\nProcessing time of :',time.time()-time_start, 's')
    














