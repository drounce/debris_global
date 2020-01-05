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
eb_fp = input.main_directory + '/../output/exp3/'
eb_fn = input.fn_prefix + 'YYYYN-' + 'XXXXE-' + input.date_start + '.nc'

#elev_cns2analyze = ['zmean']
elev_cns2analyze = ['zmean', 'zstdlow', 'zstdhigh']

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


def create_xrdataset_ostrem(ds):
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
    # Create empty datasets for each variable and merge them
    # Variable coordinates dictionary
    output_coords_dict = {
            'melt_mwea': collections.OrderedDict(
                    [('hd_cm', ds.hd_cm.values), ('stats', ds.stats.values), ('elev_cns', ds.elev_cns.values)]),
            'b0': collections.OrderedDict(
                    [('stats', ds.stats.values), ('elev_cns', ds.elev_cns.values)]),
            'k': collections.OrderedDict(
                    [('stats', ds.stats.values), ('elev_cns', ds.elev_cns.values)]),
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
            'melt_mwea': {
                    'long_name': 'annual sub-debris glacier melt',
                    'units': 'meters water equivalent per year',
                    'temporal_resolution': 'annual',
                    'start_date': input.start_date,
                    'end_date': input.end_date},
            'b0': {
                    'long_name': 'second order reaction rate initial concentration',
                    'units': 'm w.e.'},
            'k': {
                    'long_name': 'second order reaction rate initial concentration',
                    'units': 'm w.e.'}
            }
    # Add variables to empty dataset and merge together
    count_vn = 0
    encoding = {}
    noencoding_vn = ['stats', 'hd_cm', 'elev_cns', 'elev']
    for vn in ['melt_mwea', 'b0', 'k', 'elev']:
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
#        if debug:
        print(nlatlon, latlon)
        
        lat_deg = latlon[0]
        lon_deg = latlon[1]
        
        # ===== Debris Thickness vs. Surface Lowering =====
        # stats column index (0=mean)
        stats_idx = 0
        
        # Filename
        if os.path.exists(input.ostrem_fp) == False:
            os.makedirs(input.ostrem_fp)
        plot_str = str(int(lat_deg*100)) + 'N-' + str(int(lon_deg*100)) + 'E'
        ds_ostrem_fn = input.output_ostrem_fn_sample.replace('XXXX', plot_str)
        
        if os.path.exists(input.ostrem_fp + ds_ostrem_fn) == False:
            
            # Debris thickness vs. melt dataset from energy balance modeling
            ds_fn = eb_fn.replace('YYYY',str(int(lat_deg*100))).replace('XXXX',str(int(lon_deg*100)))
            ds = xr.open_dataset(eb_fp + ds_fn)
            
            ds_ostrem, encoding = create_xrdataset_ostrem(ds)
            
            # Debris thickness
            debris_thicknesses = ds.hd_cm.values
            
            debris_melt_df = pd.DataFrame(np.zeros((len(debris_thicknesses),2)), 
                                          columns=['debris_thickness', 'melt_mwea'])
        
            # Time information
            time_pd = pd.to_datetime(ds.time.values)    
            
            for nelev, elev_cn in enumerate(elev_cns2analyze):
#                if debug:
#                    print(nelev, elev_cn)
        
                for ndebris, debris_thickness in enumerate(debris_thicknesses):                    
                    melt_mwea = ds['melt'][ndebris,:,stats_idx,nelev].values.sum() / (len(time_pd)/24/365.25)
                                 
                    debris_melt_df.loc[ndebris] = debris_thickness / 100, melt_mwea
            
                # Fit curve
                fit_idx = list(np.where(debris_thicknesses >= 5)[0])            
                func_coeff, pcov = curve_fit(melt_fromdebris_func, 
                                             debris_melt_df.debris_thickness.values[fit_idx], 
                                             debris_melt_df.melt_mwea.values[fit_idx])
#                func_coeff_meltfromdebris = func_coeff.copy()
#                fit_melt = melt_fromdebris_func(debris_melt_df.debris_thickness.values, func_coeff[0], func_coeff[1])
                
                # Record ostrem curve information
                ds_ostrem['melt_mwea'][:,0,nelev] = debris_melt_df.melt_mwea.values
                ds_ostrem['b0'][0,nelev] = func_coeff[0]
                ds_ostrem['k'][0,nelev] = func_coeff[1]
            # Export netcdf
            ds_ostrem.to_netcdf(input.ostrem_fp + ds_ostrem_fn)
        
        # Plot debris vs. surface lowering
        if args.plotfigs == 1:
            
            # Open dataset
            ds_ostrem = xr.open_dataset(input.ostrem_fp + ds_ostrem_fn)
            
            fig, ax = plt.subplots(1, 1, squeeze=False, sharex=False, sharey=False, 
                                   gridspec_kw = {'wspace':0.4, 'hspace':0.15})
            
            elev_colordict = {'zmean':'k', 'zstdlow':'r', 'zstdhigh':'b'}
            elev_zorderdict = {'zmean':3, 'zstdlow':1, 'zstdhigh':1}
            elev_lwdict = {'zmean':1, 'zstdlow':0.5, 'zstdhigh':0.5}
            for nelev, elev_cn in enumerate(elev_cns2analyze):
        
#                if debug:
#                    print(nelev, elev_cn)
                debris_melt_df = pd.DataFrame(np.zeros((len(ds_ostrem.hd_cm.values),2)), 
                                          columns=['debris_thickness', 'melt_mwea'])
                debris_melt_df['debris_thickness'] = ds_ostrem.hd_cm.values / 100
                debris_melt_df['melt_mwea'] = ds_ostrem['melt_mwea'][:,stats_idx,nelev].values
                
                func_coeff = [ds_ostrem['b0'][stats_idx,nelev].values, ds_ostrem['k'][stats_idx,nelev].values]
                
                # Fitted curve
                debris_4curve = np.arange(0.02,5.01,0.01)
                melt_4curve = melt_fromdebris_func(debris_4curve, func_coeff[0], func_coeff[1])
                
                # Plot curve
                ax[0,0].plot(debris_melt_df['debris_thickness'], debris_melt_df['melt_mwea'], 'o', 
                             color=elev_colordict[elev_cn], markersize=3, markerfacecolor="None", markeredgewidth=0.75,
                             zorder=elev_zorderdict[elev_cn], label=elev_cn)
                ax[0,0].plot(debris_4curve, melt_4curve, 
                             color=elev_colordict[elev_cn], linewidth=elev_lwdict[elev_cn], linestyle='--', 
                             zorder=elev_zorderdict[elev_cn]+1)
                # text
                if nelev == 0:
                    ax[0,0].text(0.5, 1.05, plot_str, size=10, horizontalalignment='center', verticalalignment='top', 
                                 transform=ax[0,0].transAxes)
                    eqn_text = r'$b = \frac{b_{0}}{1 + kb_{0}h}$'
                    coeff1_text = r'$b_{0} = ' + str(np.round(func_coeff[0],2)) + '$' 
                    coeff2_text = r'$k = ' + str(np.round(func_coeff[1],2)) + '$' 
                    # coeff$\frac{b_{0}}{1 + 2kb_{0}h}$'
                    ax[0,0].text(0.9, 0.95, eqn_text, size=12, horizontalalignment='right', verticalalignment='top', 
                                 transform=ax[0,0].transAxes)
                    ax[0,0].text(0.615, 0.83, 'where', size=10, horizontalalignment='left', verticalalignment='top', 
                                 transform=ax[0,0].transAxes)
                    ax[0,0].text(0.66, 0.77, coeff1_text, size=10, horizontalalignment='left', verticalalignment='top', 
                                 transform=ax[0,0].transAxes)
                    ax[0,0].text(0.66, 0.7, coeff2_text, size=10, horizontalalignment='left', verticalalignment='top', 
                                 transform=ax[0,0].transAxes)
                    # X-label
                    ax[0,0].set_xlabel('Debris thickness(m)', size=12)
                    ax[0,0].set_xlim(0, 2.1)
                    #ax[0,0].set_xlim(0, debris_melt_df.debris_thickness.max())
                    ax[0,0].xaxis.set_tick_params(labelsize=12)
                    ax[0,0].xaxis.set_major_locator(plt.MultipleLocator(0.5))
                    ax[0,0].xaxis.set_minor_locator(plt.MultipleLocator(0.1))  
                    # Y-label
                    ax[0,0].set_ylabel('Melt (mwea)', size=12)
                    ax[0,0].set_ylim(0,(int(debris_melt_df.melt_mwea.values.max()/0.1)+3)*0.1)
                    ax[0,0].yaxis.set_major_locator(plt.MultipleLocator(1))
                    ax[0,0].yaxis.set_minor_locator(plt.MultipleLocator(0.1))
                    # Tick parameters
                    ax[0,0].yaxis.set_ticks_position('both')
                    ax[0,0].tick_params(axis='both', which='major', labelsize=12, direction='inout')
                    ax[0,0].tick_params(axis='both', which='minor', labelsize=10, direction='in') 
                      
            # Save plot
            ax[0,0].legend(loc=(0.65,0.45), fontsize=10, labelspacing=0.25, handlelength=1, handletextpad=0.25, 
                           borderpad=0, frameon=False)
            fig.set_size_inches(4, 4)
            figure_fn = input.output_ostrem_fn_sample.replace('XXXX',plot_str).replace('.nc','.png')
            if os.path.exists(input.ostrem_fp) == False:
                os.makedirs(input.ostrem_fp)
            fig.savefig(input.ostrem_fp + figure_fn, bbox_inches='tight', dpi=300) 

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
        latlon_list = input.latlon_list    

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