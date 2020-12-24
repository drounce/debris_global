# -*- coding: utf-8 -*-
"""
Created on Sat Aug 11 10:17:13 2018

@author: David
"""

# Built-in libraries
#import argparse
#import collections
#import multiprocessing
import os
#import pickle
#import time

# External libraries
#import rasterio
#import gdal
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
#from scipy.stats import linregress
from scipy.stats import median_absolute_deviation
import xarray as xr

# Local libraries
import debrisglobal.globaldebris_input as debris_prms
from meltcurves import melt_fromdebris_func
from meltcurves import debris_frommelt_func

#%% ===== SCRIPT OPTIONS =====

def clean_ice_melt(melt_fp, melt_fn, yearfracs):
    # Dataset of melt data
    ds_ostrem = xr.open_dataset(melt_fp + melt_fn)
    ds_ostrem = ds_ostrem.sortby('hd_cm')
    
    time_year = pd.to_datetime(ds_ostrem.time.values).year
    time_daysperyear = np.array([366 if x%4 == 0 else 365 for x in time_year])
    time_yearfrac = time_year + (pd.to_datetime(ds_ostrem.time.values).dayofyear-1) / time_daysperyear
    
    start_yearfrac = yearfracs[0]
    end_yearfrac = yearfracs[1]
    
    start_idx = np.where(abs(time_yearfrac - start_yearfrac) == abs(time_yearfrac - start_yearfrac).min())[0][0]
    end_idx = np.where(abs(time_yearfrac - end_yearfrac) == abs(time_yearfrac - end_yearfrac).min())[0][0]
    
    # Ostrem Curve
    debris_thicknesses = ds_ostrem.hd_cm.values
    debris_melt_df = pd.DataFrame(np.zeros((len(debris_thicknesses),2)), columns=['debris_thickness', 'melt_mmwed'])  
    
    nelev = 0
    for ndebris, debris_thickness in enumerate(debris_thicknesses):
        # Units: mm w.e. per day                
        melt_mmwed = (ds_ostrem['melt'][ndebris,start_idx:end_idx,nelev].values.sum() 
                      * 1000 / len(time_yearfrac[start_idx:end_idx]))
        debris_melt_df.loc[ndebris] = debris_thickness / 100, melt_mmwed
    
    print('  clean ice: ' + str(np.round(debris_melt_df.loc[0,'melt_mmwed'],2)) + ' mm w.e. d-1')
    return debris_melt_df

#%% ===== MODELED CLEAN ICE MELT =====
studies = ['Han etal 2010', 'Reid and Brock 2014', 'Brun etal 2016', 'Buri etal 2016', 'Watson etal 2017', 
           'Miles etal 2018']

if 'Han etal 2010' in studies:
    print('\nHan etal 2010:')
    melt_fp = debris_prms.output_fp + 'ostrem_curves/'
    melt_fn = '4175N-8000E-debris_melt_curve.nc'
    yearfracs = [2008 + 221/366, 2008 + 247/366]
    
    debris_melt_df = clean_ice_melt(melt_fp, melt_fn, yearfracs)


if 'Reid and Brock 2014' in studies:
    print('\nReid and Brock 2014:')
    melt_fp = debris_prms.output_fp + 'ostrem_curves/'
    melt_fn = '4650N-1050E-debris_melt_curve.nc'
    yearfracs_list = [[2010 + 156/365, 2010 + 165/365], [2011 + 156/365, 2011 + 166/365]]
    
    for yearfracs in yearfracs_list:    
        debris_melt_df = clean_ice_melt(melt_fp, melt_fn, yearfracs)
    
if 'Brun etal 2016' in studies:
    print('\nBrun etal 2016:')
    melt_fp = debris_prms.output_fp + 'ostrem_curves/'
    melt_fn = '2825N-8550E-debris_melt_curve.nc'
    yearfracs_list = [[2013 + 135/365, 2013 + 288/365], [2014 + 135/365, 2014 + 288/365]]
    
    for yearfracs in yearfracs_list:    
        debris_melt_df = clean_ice_melt(melt_fp, melt_fn, yearfracs) 


if 'Buri etal 2016' in studies:
    print('\nBuri etal 2016:')
    melt_fp = debris_prms.output_fp + 'ostrem_curves/'
    melt_fn = '2825N-8550E-debris_melt_curve.nc'
    yearfracs = [2013 + 128/365, 2013 + 140/365]
    
    debris_melt_df = clean_ice_melt(melt_fp, melt_fn, yearfracs)
    
if 'Watson etal 2017' in studies:
    print('\nWatson etal 2017:')
    melt_fp = debris_prms.output_fp + 'ostrem_curves/'
    melt_fn = '2800N-8700E-debris_melt_curve.nc'
    yearfracs = [2016 + 136/366, 2016 + 289/366]
    
    debris_melt_df = clean_ice_melt(melt_fp, melt_fn, yearfracs)

if 'Miles etal 2018' in studies:
    print('\nMiles etal 2018:')
    melt_fp = debris_prms.output_fp + 'ostrem_curves/'
    melt_fn = '2825N-8550E-debris_melt_curve.nc'
    yearfracs = [2006 + 207/365, 2015 + 43/365]
    
    debris_melt_df = clean_ice_melt(melt_fp, melt_fn, yearfracs)
        


#%%
#glaciers = [('15.04045', [0,4260], 'Immerzeel etal 2014'),
#            ('15.04045', [0,9000], 'Sakai etal 2000'),
#            ('15.04119', [0,9000], 'Miles etal 2018'),
#            ('15.04121', [0,9000], 'Miles etal 2018'),
#            ('15.04176', [0,9000], 'Miles etal 2018')]
#
#for glacier_info in glaciers:
#    print('\n')
#    
#    # Glacier info
#    glac_str = glacier_info[0]
#    elev_min = glacier_info[1][0]
#    elev_max = glacier_info[1][1]
#    reg = int(glac_str.split('.')[0])
#    
#    # Load binned data
#    try:
#        hdts_fp = debris_prms.mb_binned_fp_wdebris_hdts
#        hdts_fn = glac_str + '_mb_bins_hdts.csv'
#        hdts_df = pd.read_csv(hdts_fp + hdts_fn)
#    except:
#        hdts_fp = debris_prms.mb_binned_fp_wdebris_hdts + '../_wdebris_hdts_extrap/'
#        if reg < 10:
#            hdts_fn = (glac_str.split('.')[0].zfill(2) + '.' + glac_str.split('.')[1] + 
#                       '_mb_bins_hdts_extrap.csv')
#        else:
#            hdts_fn = glac_str + '_mb_bins_hdts_extrap.csv'
#        hdts_df = pd.read_csv(hdts_fp + hdts_fn)
#        
#        
#    hdts_df_subset = hdts_df[(hdts_df['bin_center_elev_m'] >= elev_min) & (hdts_df['bin_center_elev_m'] <= elev_max)]
#        
#    hd_glacwide = ((hdts_df_subset['dc_bin_area_valid_km2'] * hdts_df_subset['hd_ts_mean_m']).sum() / 
#                   hdts_df_subset['dc_bin_area_valid_km2'].sum())
#    mf_glacwide = ((hdts_df_subset['dc_bin_area_valid_km2'] * hdts_df_subset['mf_ts_mean']).sum() / 
#                   hdts_df_subset['dc_bin_area_valid_km2'].sum())
#    
#    print(glac_str, '  (' + glacier_info[2] + ')')
#    print('  hd: ' + str(np.round(hd_glacwide,2)) + '  mf: ' + str(np.round(mf_glacwide,2)))
            
            
