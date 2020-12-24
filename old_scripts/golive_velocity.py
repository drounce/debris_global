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
import xarray as xr

# Local libraries
import globaldebris_input as input


#%%
golive_fp = '/Users/davidrounce/Documents/Dave_Rounce/Satellite_Images/GoLIVE/p064_r017-Kennicott/'
tdays = 32

golive_fns = []
for i in os.listdir(golive_fp):
    if i.endswith('.nc') and i.endswith('_nrt.nc') == False:
        if int(i.split('_')[3]) == tdays:
            golive_fns.append(i)
    
vx_stack = None
for nfn, fn in enumerate(golive_fns):
#for nfn, fn in enumerate([golive_fns[0]]):
#for nfn, fn in enumerate(golive_fns[0:2]):
    ds = xr.open_dataset(golive_fp + fn)
    
    print(fn, ds['vx'].values.shape)
#    if vx_stack is None:
#        vx_stack = ds['vx'].values
#    else:
#        print(vx_stack.shape, ds['vx'].values.shape)
#        vx_stack = np.dstack((vx_stack, ds['vx'].values))
        
    
    
    
    
#    n_nearest = 0
#    n_success = 0
#    min_n_nearest = 10
#    hd_ts_list = []
#    while n_nearest < n_glac_nearest and n_success < min_n_nearest:
#        rgi_str_nearest = rgiid_nearest_list[n_nearest]
#
#        # Load parameters
#        if int(region) < 10:
#            rgi_str_nearest = str(int(rgi_str_nearest.split('.')[0])) + '.' + rgi_str_nearest.split('.')[1]
#        df_opt_fn = rgi_str_nearest + '_hdts_opt_params.csv'
#        df_opt = pd.read_csv(ts_opt_fp + df_opt_fn)
#        func_coeff = [df_opt.loc[0,'a'], df_opt.loc[0,'b'], df_opt.loc[0,'c']]
#        ts_offset = df_opt.loc[0,'a_offset']
#        
#        # Estimate debris thickness
#        hd_array = debris_fromts_maskedarray(gf.ts, func_coeff[0]+ts_offset, func_coeff[1], func_coeff[2])
#        hd_array[hd_array>input.hd_max] = input.hd_max
#        hd_ma = np.ma.array(hd_array, mask=dc_mask)
#        hd_ma_median = np.median(hd_ma.compressed())
#        
#        # Only include estimates if they are plausible
#        if hd_ma_median > 0 and hd_ma_median < input.hd_max:
#            hd_ts_list.append(hd_ma)
#            n_success += 1
#        
#        n_nearest += 1
#        
#        
#        
#    # Debris thickness based on median of the plausible nearest values
#    hd_ts_all = np.ma.array(hd_ts_list)
#    hd_ts_med = np.median(hd_ts_all, axis=0)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

    
        