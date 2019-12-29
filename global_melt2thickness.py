# -*- coding: utf-8 -*-
"""
Created on Sat Aug 11 10:17:13 2018

@author: David
"""

# Built-in libraries
import collections
import os

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

#%% ===== SCRIPT SPECIFIC INFORMATION =====
debug = True
eb_fp = input.main_directory + '/../output/exp3/'
eb_fn = input.fn_prefix + 'YYYYN-' + 'XXXXE-' + input.date_start + '.nc'
ostrem_fp = input.main_directory + '/../output/ostrem_curves/'

#elev_cns2analyze = ['zmean']
elev_cns2analyze = ['zmean', 'zstdlow', 'zstdhigh']

#%% ===== FUNCTIONS =====
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
    Create empty xarray dataset that will be used to record simulation runs.

    Parameters
    ----------
    main_glac_rgi : pandas dataframe
        dataframe containing relevant rgi glacier information
    dates_table : pandas dataframe
        table of the dates, months, days in month, etc.
    sim_iters : int
        number of simulation runs included
    stat_cns : list
        list of strings containing statistics that will be used on simulations
    record_stats : int
        Switch to change from recording simulations to statistics

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

#%%
for nlatlon, latlon in enumerate(input.latlon_list):
    if debug:
        print(nlatlon, latlon)
    
    lat_deg = latlon[0]
    lon_deg = latlon[1]
    
        
    # Debris thickness vs. melt dataset from energy balance modeling
    ds_fn = eb_fn.replace('YYYY',str(int(lat_deg*100))).replace('XXXX',str(int(lon_deg*100)))
    ds = xr.open_dataset(eb_fp + ds_fn)
    
    ds_ostrem, encoding = create_xrdataset_ostrem(ds)
    
    # Debris thickness
    debris_thicknesses = ds.hd_cm.values

    # Time information
    time_pd = pd.to_datetime(ds.time.values)    
    
    # ===== Debris Thickness vs. Surface Lowering =====
    debris_melt_df = pd.DataFrame(np.zeros((len(debris_thicknesses),2)), columns=['debris_thickness', 'melt_mwea'])
    
    # stats column index (0=mean)
    stats_idx = 0
    
    # Plot
    plot_str = str(int(lat_deg*100)) + 'N-' + str(int(lon_deg*100)) + 'E'
    fig, ax = plt.subplots(1, 1, squeeze=False, sharex=False, sharey=False, 
                           gridspec_kw = {'wspace':0.4, 'hspace':0.15})
    
    elev_colordict = {'zmean':'k', 'zstdlow':'r', 'zstdhigh':'b'}
    elev_zorderdict = {'zmean':3, 'zstdlow':1, 'zstdhigh':1}
    elev_lwdict = {'zmean':1, 'zstdlow':0.5, 'zstdhigh':0.5}
    for nelev, elev_cn in enumerate(elev_cns2analyze):
#    for nelev, elev_cn in enumerate([ds.elev_cns.values[0]]):
        print(nelev, elev_cn)

        for ndebris, debris_thickness in enumerate(debris_thicknesses):
#            print(ndebris, debris_thickness)
            
            melt_mwea = ds['melt'][ndebris,:,stats_idx,nelev].values.sum() / (len(time_pd)/24/365.25)
                         
            debris_melt_df.loc[ndebris] = debris_thickness / 100, melt_mwea
    
        fit_idx = list(np.where(debris_thicknesses >= 5)[0])
#        fit_idx = list(debris_melt_df.index.values)
#        fit_idx = fit_idx[2:]
#        print('Fit to a subset of data (not thinnest values) to capture this portion of the curve well' + 
#              '\n  - any discrepancies at < 0.1 m are going to be very small because melt varies so drastically\n\n')
        
        # Fit curve
        func_coeff, pcov = curve_fit(melt_fromdebris_func, 
                                     debris_melt_df.debris_thickness.values[fit_idx], 
                                     debris_melt_df.melt_mwea.values[fit_idx])
        func_coeff_meltfromdebris = func_coeff.copy()
        fit_melt = melt_fromdebris_func(debris_melt_df.debris_thickness.values, func_coeff[0], func_coeff[1])
        debris_4curve = np.arange(0.02,5.01,0.01)
        melt_4curve = melt_fromdebris_func(debris_4curve, func_coeff[0], func_coeff[1])
        
        # Record ostrem curve information
        ds_ostrem['melt_mwea'][:,0,nelev] = debris_melt_df.melt_mwea.values
        ds_ostrem['b0'][0,nelev] = func_coeff[0]
        ds_ostrem['k'][0,nelev] = func_coeff[1]
        
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
    ax[0,0].legend(loc=(0.65,0.45), fontsize=10, labelspacing=0.25, handlelength=1, handletextpad=0.25, borderpad=0, 
                   frameon=False)
    fig.set_size_inches(4, 4)
    figure_fn = plot_str + '_debris_melt_curve.png'
    if os.path.exists(ostrem_fp) == False:
        os.makedirs(ostrem_fp)
    fig.savefig(ostrem_fp + figure_fn, bbox_inches='tight', dpi=300) 

        
#%% ===== Input Data =====
#glac_str = '15.03473'
##larsen_data_fullfn = '/Users/davidrounce/Documents/Dave_Rounce/HiMAT/DEMs/larsen/data/Kennicott.2000.2013.output.txt'
##larsen_ice_density = 850 # kg m-3
#mb_fullfn = ('/Users/davidrounce/Documents/Dave_Rounce/HiMAT/DEMs/Shean_2019_0213/' + 
#             'mb_combined_20190213_nmad_bins/' + glac_str + '_mb_bins.csv')
#glacier_name = 'Ngozumpa Glacier (RGI60-01.15645)'
#output_data_fp = input.main_directory + '/hma_data/output/change/'
##output_data_fp = input.main_directory + '/output/exp3/'
#output_fp = input.main_directory + '/hma_data/output/'
##anderson_debris_fullfn = input.main_directory + '/kennicott_data/anderson_debristhickness.csv'
##anderson_debris = pd.read_csv(anderson_debris_fullfn)

###%% ===== DERIVE DEBRIS THICKNESS FROM MASS BALANCE DATA =====
##mb_df = pd.read_csv(mb_fullfn)
##mb_df.loc[:,:] = mb_df.values.astype(np.float64)
##mb_df['debris_thickness'] = debris_frommelt_func(-1*mb_df[' mb_bin_mean_mwea'].values, func_coeff[0], func_coeff[1])
##mb_df['mb_bin_mean_mwea_1stdlow'] = mb_df[' mb_bin_mean_mwea'] - mb_df[' mb_bin_std_mwea']
##mb_df['mb_bin_mean_mwea_1stdhigh'] = mb_df[' mb_bin_mean_mwea'] + mb_df[' mb_bin_std_mwea']
##mb_df['debris_thickness_1stdlow'] = debris_frommelt_func(-1*mb_df['mb_bin_mean_mwea_1stdlow'].values, 
##                                                          func_coeff[0], func_coeff[1])
##mb_df['debris_thickness_1stdhigh'] = debris_frommelt_func(-1*mb_df['mb_bin_mean_mwea_1stdhigh'].values, 
##                                                          func_coeff[0], func_coeff[1])
##
### PLOT MASS BALANCE AND DEBRIS THICKNESS VS. ELEVATION WITHOUT EMERGENCE VELOCITIES
### Mass balance
##fig, ax = plt.subplots(1, 2, squeeze=False, sharex=False, sharey=True, gridspec_kw = {'wspace':0.1, 'hspace':0.15})
##ax[0,0].plot(mb_df[' mb_bin_mean_mwea'], mb_df['# bin_center_elev_m'], color='k', linewidth=1, zorder=2)
##ax[0,0].fill_betweenx(mb_df['# bin_center_elev_m'], 
##                      mb_df['mb_bin_mean_mwea_1stdlow'], mb_df['mb_bin_mean_mwea_1stdhigh'], 
##                      facecolor='k', alpha=0.2, zorder=1)
### text
##fig.text(0.5, 0.96, glacier_name, size=10, horizontalalignment='center', verticalalignment='top', 
##             transform=ax[0,0].transAxes)
### X-label
##ax[0,0].set_xlabel('Mass Balance\n(mwea)', size=12)
##ax[0,0].set_xlim(-3.25, 0.25)
###ax[0,0].set_xlim(0, debris_melt_df.debris_thickness.max())
##ax[0,0].xaxis.set_tick_params(labelsize=12)
##ax[0,0].xaxis.set_major_locator(plt.MultipleLocator(1))
##ax[0,0].xaxis.set_minor_locator(plt.MultipleLocator(0.1))  
### Y-label
##ylower = mb_df['# bin_center_elev_m'].min()
##yupper = ylower + 1000
##ax[0,0].set_ylabel('Elevation (m asl)', size=12)
##ax[0,0].set_ylim(ylower,yupper)
##ax[0,0].yaxis.set_major_locator(plt.MultipleLocator(100))
##ax[0,0].yaxis.set_minor_locator(plt.MultipleLocator(50))
### Tick parameters
##ax[0,0].yaxis.set_ticks_position('both')
##ax[0,0].tick_params(axis='both', which='major', labelsize=12, direction='inout')
##ax[0,0].tick_params(axis='both', which='minor', labelsize=10, direction='in')       
##
### Debris thickness
##ax[0,1].plot(mb_df['debris_thickness'], mb_df['# bin_center_elev_m'], color='k', linewidth=1, zorder=2)
###ax[0,0].plot(mb_df['MassBal_25'], mb_df['E'], color='k', linewidth=0.5, zorder=2)
###ax[0,0].plot(mb_df['MassBal_75'], mb_df['E'], color='k', linewidth=0.5, zorder=2)
##ax[0,1].fill_betweenx(mb_df['# bin_center_elev_m'], 
##                      mb_df['debris_thickness_1stdlow'], mb_df['debris_thickness_1stdhigh'], 
##                      facecolor='k', alpha=0.2, zorder=1)
### X-label
##ax[0,1].set_xlabel('Debris thickness\n(m)', size=12)
##ax[0,1].set_xlim(0, 1.5) 
##ax[0,1].xaxis.set_major_locator(plt.MultipleLocator(0.5))
##ax[0,1].xaxis.set_minor_locator(plt.MultipleLocator(0.1))  
### Y-label
##ax[0,1].set_ylim(ylower,yupper)
### Tick parameters
##ax[0,1].yaxis.set_ticks_position('both')
##ax[0,1].tick_params(axis='both', which='major', labelsize=12, direction='inout')
##ax[0,1].tick_params(axis='both', which='minor', labelsize=10, direction='in')       
##        
### Save plot
##fig.set_size_inches(4, 4)
##figure_fn = 'elev_mb_debris.png'
##fig.savefig(output_fp + figure_fn, bbox_inches='tight', dpi=300) 
##
###%% ===== EMERGENCE VELOCITIES ==========================================
### CURRENTLY PROCESSED IN IPYTHON NOTEBOOK USING SHEAN GMBTOOLS-TYPE FILE 'emergence_velocity.ipynb'
##emergence_fullfn = ('/Users/davidrounce/Documents/Dave_Rounce/HiMAT/DEMs/Shean_2019_0213/' + 
##                    'mb_combined_20190213_nmad_bins/csv/' + glac_str + '_mb_bins_wemvel.csv')
##emergence_df = pd.read_csv(emergence_fullfn)
##emergence_df['area_cumsum'] = np.cumsum(emergence_df['z1_bin_area_valid_km2'])
##
##binsize_mb = mb_df.loc[1,'# bin_center_elev_m'] - mb_df.loc[0,'# bin_center_elev_m']
##emvel_binsize = emergence_df['# bin_center_elev_m'].values[1] - emergence_df['# bin_center_elev_m'].values[0]
##emergence_shift_idx = np.where(emergence_df.area_cumsum.values < mb_df.loc[0,' z1_bin_area_valid_km2'])[0][-1]
##mb_offset = (mb_df.loc[0, '# bin_center_elev_m'] + binsize_mb/2 - emvel_binsize / 2 - 
##             emergence_df.loc[emergence_shift_idx, '# bin_center_elev_m'])
##                      
##emergence_df['# bin_center_elev_m'] = emergence_df['# bin_center_elev_m'] + mb_offset
##emergence_df['E_low'] = emergence_df['# bin_center_elev_m'] - emvel_binsize/2
##emergence_df['E_high'] = emergence_df['# bin_center_elev_m'] + emvel_binsize/2
##
### Get mean emergence velocity to coincide with elevation bins
##mb_df['E_low'] =  mb_df['# bin_center_elev_m'] - binsize_mb/2
##mb_df['E_high'] = mb_df['# bin_center_elev_m'] + binsize_mb/2
##
##mb_df['em_idx_low'] = np.nan
##mb_df['em_idx_high'] = np.nan
##for x in mb_df.index.values:
##    rows_low = np.where(mb_df.E_low.values[x] == emergence_df.E_low.values)[0]
##    if len(rows_low) > 0:
##        mb_df.loc[x,'em_idx_low'] = rows_low[0]
##    elif x == 0:
##        mb_df.loc[x,'em_idx_low'] = 0
##        
##    rows_high = np.where(mb_df.E_high.values[x] == emergence_df.E_high.values)[0]
##    if len(rows_high) > 0:
##        mb_df.loc[x,'em_idx_high'] = rows_high[0]
##    elif len(rows_high) == 0 and ~np.isnan(mb_df.loc[x,'em_idx_low']):
##        mb_df.loc[x,'em_idx_high'] = emergence_df.index.values[-1]
##
##emergence_df['emvel*area'] = emergence_df.emvel_mean * emergence_df.z1_bin_area_valid_km2
##emergence_df['emvel*area_1stdlow'] = (emergence_df.emvel_mean - emergence_df.emvel_std) * emergence_df.z1_bin_area_valid_km2
##emergence_df['emvel*area_1stdhigh'] = (emergence_df.emvel_mean + emergence_df.emvel_std) * emergence_df.z1_bin_area_valid_km2
##            
##mb_df['emvel_myr'] = np.nan
##for x in mb_df.index.values:
##    if ~np.isnan(mb_df.loc[x,'em_idx_low']):
##        mb_df.loc[x,'emvel_myr'] = (
##                emergence_df.loc[mb_df.loc[x,'em_idx_low']:mb_df.loc[x,'em_idx_high'], 'emvel*area'].sum() / 
##                emergence_df.loc[mb_df.loc[x,'em_idx_low']:mb_df.loc[x,'em_idx_high'], 'z1_bin_area_valid_km2'].sum())
###larsen_data['emvel_myr_1stdlow'] = (
###        [emvel_data.loc[larsen_data.loc[x,'em_idx_low']:larsen_data.loc[x,'em_idx_high'], 'emvel*area_1stdlow'].sum() / 
###         emvel_data.loc[larsen_data.loc[x,'em_idx_low']:larsen_data.loc[x,'em_idx_high'], 'z1_bin_area_valid_km2'].sum()
###         for x in larsen_data.index.values])
###larsen_data['emvel_myr_1stdhigh'] = (
###        [emvel_data.loc[larsen_data.loc[x,'em_idx_low']:larsen_data.loc[x,'em_idx_high'], 'emvel*area_1stdhigh'].sum() / 
###         emvel_data.loc[larsen_data.loc[x,'em_idx_low']:larsen_data.loc[x,'em_idx_high'], 'z1_bin_area_valid_km2'].sum()
###         for x in larsen_data.index.values])
##
##mb_df['mb_wem'] = mb_df[' mb_bin_mean_mwea'] - mb_df['emvel_myr']
####larsen_data['mb_wem_25'] = larsen_data['MassBal_25'] - larsen_data['emvel_myr_1stdhigh']
####larsen_data['mb_wem_75'] = larsen_data['MassBal_75'] - larsen_data['emvel_myr_1stdlow']
###larsen_data['mb_wem_25'] = larsen_data['MassBal_25'] - larsen_data['emvel_myr']
###larsen_data['mb_wem_75'] = larsen_data['MassBal_75'] - larsen_data['emvel_myr']
###larsen_data['debris_thickness_wem'] = debris_frommelt_func(-1*larsen_data['mb_wem'].values, func_coeff[0], func_coeff[1])
###larsen_data['debris_thickness_wem_25'] = debris_frommelt_func(-1*larsen_data['mb_wem_25'].values, 
###                                                          func_coeff[0], func_coeff[1])
###larsen_data['debris_thickness_wem_75'] = debris_frommelt_func(-1*larsen_data['mb_wem_75'].values, 
###                                                          func_coeff[0], func_coeff[1])
##
##mb_df['debris_thickness_wem'] = debris_frommelt_func(-1*mb_df['mb_wem'].values, func_coeff[0], func_coeff[1])
##
###%% ===== PLOT MASS BALANCE AND DEBRIS THICKNESS VS. ELEVATION WITH EMERGENCE VELOCITIES ======
### Mass balance
##fig, ax = plt.subplots(1, 2, squeeze=False, sharex=False, sharey=True, gridspec_kw = {'wspace':0.1, 'hspace':0.15})
##ax[0,0].plot(mb_df['mb_wem'], mb_df['# bin_center_elev_m'], color='k', linewidth=1, zorder=2)
###ax[0,0].fill_betweenx(mb_df['# bin_center_elev_m'], 
###                      mb_df['mb_bin_mean_mwea_1stdlow'], mb_df['mb_bin_mean_mwea_1stdhigh'], 
###                      facecolor='k', alpha=0.2, zorder=1)
### text
##fig.text(0.5, 0.96, glacier_name, size=10, horizontalalignment='center', verticalalignment='top', 
##             transform=ax[0,0].transAxes)
### X-label
##ax[0,0].set_xlabel('Mass Balance\n(mwea)', size=12)
##ax[0,0].set_xlim(-3.25, 0.25)
###ax[0,0].set_xlim(0, debris_melt_df.debris_thickness.max())
##ax[0,0].xaxis.set_tick_params(labelsize=12)
##ax[0,0].xaxis.set_major_locator(plt.MultipleLocator(1))
##ax[0,0].xaxis.set_minor_locator(plt.MultipleLocator(0.1))  
### Y-label
##ylower = mb_df['# bin_center_elev_m'].min()
##yupper = ylower + 1000
##ax[0,0].set_ylabel('Elevation (m asl)', size=12)
##ax[0,0].set_ylim(ylower,yupper)
##ax[0,0].yaxis.set_major_locator(plt.MultipleLocator(100))
##ax[0,0].yaxis.set_minor_locator(plt.MultipleLocator(50))
### Tick parameters
##ax[0,0].yaxis.set_ticks_position('both')
##ax[0,0].tick_params(axis='both', which='major', labelsize=12, direction='inout')
##ax[0,0].tick_params(axis='both', which='minor', labelsize=10, direction='in')       
##
### Debris thickness
##ax[0,1].plot(mb_df['debris_thickness_wem'], mb_df['# bin_center_elev_m'], color='k', linewidth=1, zorder=2)
###ax[0,0].plot(mb_df['MassBal_25'], mb_df['E'], color='k', linewidth=0.5, zorder=2)
###ax[0,0].plot(mb_df['MassBal_75'], mb_df['E'], color='k', linewidth=0.5, zorder=2)
###ax[0,1].fill_betweenx(mb_df['# bin_center_elev_m'], 
###                      mb_df['debris_thickness_1stdlow'], mb_df['debris_thickness_1stdhigh'], 
###                      facecolor='k', alpha=0.2, zorder=1)
### X-label
##ax[0,1].set_xlabel('Debris thickness\n(m)', size=12)
##ax[0,1].set_xlim(0, 1.5) 
##ax[0,1].xaxis.set_major_locator(plt.MultipleLocator(0.5))
##ax[0,1].xaxis.set_minor_locator(plt.MultipleLocator(0.1))  
### Y-label
##ax[0,1].set_ylim(ylower,yupper)
### Tick parameters
##ax[0,1].yaxis.set_ticks_position('both')
##ax[0,1].tick_params(axis='both', which='major', labelsize=12, direction='inout')
##ax[0,1].tick_params(axis='both', which='minor', labelsize=10, direction='in')       
##        
### Save plot
##fig.set_size_inches(4, 4)
##figure_fn = 'elev_mb_debris_wem.png'
##fig.savefig(output_fp + figure_fn, bbox_inches='tight', dpi=300) 