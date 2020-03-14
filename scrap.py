#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Mar 13 22:29:38 2020

@author: davidrounce
"""
import os
import xarray as xr

ds_fp = '/Users/davidrounce/Documents/Dave_Rounce/DebrisGlaciers_WG/Melt_Intercomparison/climate_data/HMA/'
mv_fns = []
keep_fns = []
for i in os.listdir(ds_fp):
    if i.endswith('.nc'):
        ds = xr.open_dataset(ds_fp + i)
        #%%
        if 'dc_zmean' not in list(ds.keys()):
            mv_fns.append(i)
        else:
            keep_fns.append(i)
        