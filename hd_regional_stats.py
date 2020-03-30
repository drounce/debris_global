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
#import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
#from scipy.optimize import curve_fit
#from scipy.stats import linregress
#from scipy.stats import median_absolute_deviation
#import xarray as xr

import debrisglobal.globaldebris_input as debris_prms


#%%% ===== SCRIPT OPTIONS =====
#if not os.path.exists(melt_compare_fp):
#    os.makedirs(melt_compare_fp)
    
#rois = ['01','02','03','04','05','06','07','08','09','10','11','12','HMA','16','17','18']
rois = ['01','02','03','04','05','06','07','08','09','11','12','13','14', '15', '16','17','18']
#rois = ['18']

hd_cn = 'hd_ts_mean_m'
mf_cn = 'mf_ts_mean'

reg_stats_cns = ['roi', 'hd_mean', 'mf_mean']
reg_stats_df = pd.DataFrame(np.zeros((len(rois)+1,len(reg_stats_cns))), columns=reg_stats_cns)

## ===== REGIONAL MELT FACTOR STATISTICS =====
hdts_bin_fp = debris_prms.mb_binned_fp + '_wdebris_hdts/'
hdts_bin_fp_extrap = debris_prms.mb_binned_fp + '_wdebris_hdts_extrap/'

area_cumsum_all = 0
hdxarea_cumsum_all = 0
mfxarea_cumsum_all = 0
for nroi, roi in enumerate(rois):
    
    print(roi)
    
    # Glaciers optimized
    glac_hd_fullfns = []
    for i in os.listdir(hdts_bin_fp):
        if i.endswith('hdts.csv'):
            reg_str = str(int(i.split('.')[0])).zfill(2)
#            if region in debris_prms.roi_rgidict[roi]:
            if reg_str == roi:
                glac_hd_fullfns.append(hdts_bin_fp + i)
    
    # Glaciers extrapolated
    for i in os.listdir(hdts_bin_fp_extrap):
        if i.endswith('hdts_extrap.csv'):
#            region = int(i.split('.')[0])
#            if region in debris_prms.roi_rgidict[roi]:
            reg_str = str(int(i.split('.')[0])).zfill(2)
            if reg_str == roi:
                glac_hd_fullfns.append(hdts_bin_fp_extrap + i)
    glac_hd_fullfns = sorted(glac_hd_fullfns)
    
      
    area_cumsum = 0
    hdxarea_cumsum = 0
    mfxarea_cumsum = 0
    for nfn, fullfn in enumerate(glac_hd_fullfns):
#    for nfn, fullfn in enumerate([glac_hd_fullfns[0]]):
        if nfn%500 == 0:
            print('  ', nfn)
        df = pd.read_csv(fullfn)
        
    #     print(df.loc[:,['bin_center_elev_m', ' z1_bin_area_valid_km2', ' dc_bin_area_valid_km2', 
    #                     ' dc_bin_area_perc', 'debris_perc', 'hd_ts_med', 'mf_ts_med']])
        
        if 'hd_ts_mean_m' in list(df.columns):
            area_cumsum += np.sum(df['dc_bin_area_valid_km2'].values)
            hdxarea_cumsum += np.sum(df['dc_bin_area_valid_km2'] * df[hd_cn])
            mfxarea_cumsum += np.sum(df['dc_bin_area_valid_km2'] * df[mf_cn])
            
    #     print(area_cumsum)
    #     print(hdxarea_cumsum)
    #     print(mfxarea_cumsum)
    
    # Compute average
    hd_center = hdxarea_cumsum / area_cumsum    
    mf_center = mfxarea_cumsum / area_cumsum
    
    # All statistics
    area_cumsum_all += area_cumsum
    hdxarea_cumsum_all += hdxarea_cumsum
    mfxarea_cumsum_all += mfxarea_cumsum
    
    # Record stats
    reg_stats_df.loc[nroi,:] = [roi, hd_center, mf_center]
    
    print(roi, hd_center, mf_center)
    
# All regions statistic
hd_center_all = hdxarea_cumsum_all / area_cumsum_all
mf_center_all = mfxarea_cumsum_all / area_cumsum_all
nroi += 1
reg_stats_df.loc[nroi,:] = ['all', hd_center_all, mf_center_all]
    
# Export regional statistics
reg_stats_df.to_csv(debris_prms.output_fp + 'reg_stats_hd_mf.csv', index=False)
