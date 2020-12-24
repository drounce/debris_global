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
from meltcurves import debris_frommelt_func

#%% ===== SCRIPT SPECIFIC INFORMATION =====
stats_idx = 0

#%%
# ===== SELECT GLACIERS WITH DATA =====
rgiid_list = []
rgiid_fn_list = []
for i in os.listdir(input.mb_binned_fp):
    if i.endswith('mb_bins.csv'):
        rgiid_list.append(i[0:8])
        rgiid_fn_list.append(i)
        
rgiid_list = sorted(rgiid_list)
rgiid_fn_list = sorted(rgiid_fn_list)

print(len(rgiid_list))

main_glac_rgi = input.selectglaciersrgitable(rgiid_list)
main_glac_rgi['bin_fn'] = rgiid_fn_list

# ==== DETERMINE NEAREST LAT/LON =====
ds_elev = xr.open_dataset(input.metdata_fp + input.metdata_elev_fn)
#  argmin() finds the minimum distance between the glacier lat/lon and the GCM pixel
lat_nearidx = (np.abs(main_glac_rgi['CenLat'].values[:,np.newaxis] - 
                      ds_elev['latitude'][:].values).argmin(axis=1))
lon_nearidx = (np.abs(main_glac_rgi['CenLon'].values[:,np.newaxis] - 
                      ds_elev['longitude'][:].values).argmin(axis=1))

latlon_nearidx = list(zip(lat_nearidx, lon_nearidx))
latlon_nearidx_unique = sorted(list(set(latlon_nearidx)))

main_glac_rgi['lat_nearest'] = ds_elev['latitude'][lat_nearidx].values 
main_glac_rgi['lon_nearest'] = ds_elev['longitude'][lon_nearidx].values 
main_glac_rgi['latlon_nearidx'] = latlon_nearidx
latlon_unique_dict = dict(zip(latlon_nearidx_unique,np.arange(0,len(latlon_nearidx_unique))))
latlon_unique_dict_reversed = dict(zip(np.arange(0,len(latlon_nearidx_unique)),latlon_nearidx_unique))
main_glac_rgi['latlon_unique_no'] = main_glac_rgi['latlon_nearidx'].map(latlon_unique_dict)

#%%
# ===== ESTIMATE DEBRIS THICKNESS =====
latlon_nodata = []
glac_nodata = []
#for nglac, glac_idx in enumerate(main_glac_rgi.index.values):
for nglac, glac_idx in enumerate([main_glac_rgi.index.values[6940]]):
#for nglac, glac_idx in enumerate([main_glac_rgi.index.values[3224]]):
    glac_str = main_glac_rgi.loc[glac_idx,'rgino_str']
    print(nglac, glac_idx, glac_str)
    
    # ===== Load Mass Balance Data =====
    mb_df_fn = main_glac_rgi.loc[glac_idx,'bin_fn']
    mb_df = pd.read_csv(input.mb_binned_fp + mb_df_fn)
    mb_df.loc[:,:] = mb_df.values.astype(np.float64)
    
    # ===== Load Ostrem data =====
    latlon_str = (str(int(main_glac_rgi.loc[glac_idx,'lat_nearest']*100)) + 'N-' + 
                  str(int(main_glac_rgi.loc[glac_idx,'lon_nearest']*100)) + 'E') 
    ostrem_fn = input.output_ostrem_fn_sample.replace('XXXX', latlon_str)
    
    
    if os.path.exists(input.output_ostrem_fp + ostrem_fn) and ' perc_debris' in mb_df.columns:
        ds_ostrem = xr.open_dataset(input.output_ostrem_fp + ostrem_fn)
        
        # Nearest elevation
        # Debris elevation (weighted mean)
        mb_area_total = main_glac_rgi.loc[glac_idx,'Area']
        mb_df['bin_area'] = mb_df[' z1_bin_area_perc'].values / 100 * mb_area_total
        mb_df['bin_area_debris'] = mb_df[' perc_debris'].values / 100 * mb_df['bin_area'].values
        debris_elev_mean = ((mb_df['# bin_center_elev_m'] * mb_df['bin_area_debris']).sum() 
                            /  mb_df['bin_area_debris'].sum()) 
        ds_elevidx = np.abs(debris_elev_mean - ds_ostrem.elev.values).argmin()
        
        # Melt function coefficients
        func_coeff = [float(ds_ostrem['b0'][stats_idx,ds_elevidx].values), 
                      float(ds_ostrem['k'][stats_idx,ds_elevidx].values)]
        
        # Estimate debris thickness
        mb_df['mb_bin_mean_mwea_1stdlow'] = mb_df[' mb_bin_mean_mwea'] - mb_df[' mb_bin_std_mwea']
        mb_df['mb_bin_mean_mwea_1stdhigh'] = mb_df[' mb_bin_mean_mwea'] + mb_df[' mb_bin_std_mwea']
        mb_df['hd'] = debris_frommelt_func(-1*mb_df[' mb_bin_mean_mwea'].values, 
                                                         func_coeff[0], func_coeff[1])
        mb_df['hd_1stdlow'] = debris_frommelt_func(-1*mb_df['mb_bin_mean_mwea_1stdlow'].values, 
                                                                  func_coeff[0], func_coeff[1])
        mb_df['hd_1stdhigh'] = debris_frommelt_func(-1*mb_df['mb_bin_mean_mwea_1stdhigh'].values, 
                                                                  func_coeff[0], func_coeff[1])

        # ==== Add Emergence Velocity =====
        # Load file
        emvel_fn = input.outdir_emvel_fp + glac_str + input.output_emvel_csv_ending
        emvel_df = pd.read_csv(emvel_fn)
        emvel_df.loc[:,:] = np.nan_to_num(emvel_df.values.astype(np.float64), 0)
        
        mb_df[' vm_med_itslive'] = emvel_df[' vm_med']
        mb_df[' vm_mad_itslive'] = emvel_df[' vm_mad']
        mb_df = pd.concat([mb_df, emvel_df[[' emvel_mean', ' emvel_std']]], axis=1)
        # Uncertainty with flux divergence from Farinotti et al. (2019)
        mb_df[' emvel_high'] = mb_df[' emvel_mean'] * 1.6
        mb_df[' emvel_low'] = mb_df[' emvel_mean'] * 0.8
    
        # Mass balance with emergence velocity and uncertainties
        mb_df['mb_wem'] = mb_df[' mb_bin_mean_mwea'] - mb_df[' emvel_mean']
        # Higher emergence --> more melt
        mb_df['mb_wemthick'] = mb_df[' mb_bin_mean_mwea'] - mb_df[' emvel_high'] - mb_df[' mb_bin_std_mwea']
        # Lower emergence --> less melt
        mb_df['mb_wemthin'] = mb_df[' mb_bin_mean_mwea'] - mb_df[' emvel_low'] + mb_df[' mb_bin_std_mwea']
        
        # Debris thickness with emergence velocity and uncertainties
        mb_df['hd_wem'] = debris_frommelt_func(-1*mb_df['mb_wem'].values, func_coeff[0], func_coeff[1])
        mb_df['hd_wemthick'] = debris_frommelt_func(-1*mb_df['mb_wemthick'].values, func_coeff[0], func_coeff[1])
        mb_df['hd_wemthin'] = debris_frommelt_func(-1*mb_df['mb_wemthin'].values, func_coeff[0], func_coeff[1])
    
        # Export dataset
        if os.path.exists(input.mb_binned_fp_wdebris) == False:
            os.makedirs(input.mb_binned_fp_wdebris)
        
        mb_df_fn_wdebris = mb_df_fn.replace('.csv','_wdebris.csv')
        mb_df.to_csv(input.mb_binned_fp_wdebris + mb_df_fn_wdebris, index=False)
        ds_ostrem.close()
    elif os.path.exists(input.output_ostrem_fp + ostrem_fn) == False:
        latlon_nodata.append(latlon_str)
        glac_nodata.append(main_glac_rgi.loc[glac_idx,'rgino_str'])

latlon_nodata = sorted(list(set(latlon_nodata)))
if len(latlon_nodata) > 0:
    print('\nMissing Ostrem data:\n', latlon_nodata,'\n')
    print('Glaciers without Ostrem data:\n', glac_nodata,'\n')
    print('\nTHESE GLACIERS DONT HAVE CONSIDERABLE DEBRIS AND WERE THEREFORE EXCLUDED\n')
        
#%% ===== Input Data =====
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
#
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