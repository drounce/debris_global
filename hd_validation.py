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
#import xarray as xr

# Local libraries
import globaldebris_input as input


#%%
#glaciers = ['1.15645', '15.03473', '15.03733', '15.03734', '15.04121', '15.04045', '14.06794', '14.04477', 
#            '13.43232', '15.07886', '14.16042']
glaciers = ['1.15645', '15.04045', '15.03473', '15.03743', '14.06794']
#glaciers = ['15.04045']
glac_name_dict = {'1.15645':'Kennicott',
                  '15.04045':'Lirung',
                  '15.03473':'Ngozumpa',
                  '15.03743':'Imja',
                  '14.06794':'Baltoro'}

hd_obs_fp = input.main_directory + '/../hd_obs/'
hd_ts_fp = input.main_directory + '/../hd_obs/hd_ts_csv/'
hd_fn_dict = {'1.15645': 'kennicott_anderson_2019.csv',
              '15.04045': 'lirung_nicholson_gpr.csv',
              '15.03473': 'ngoz_mccarthy_gpr.csv',
              '15.03743': 'imja_rounce2014.csv',
              '14.06794': 'baltoro_mihalcea2006.csv'}

glac_symbol = {'1.15645':['o',20,'none'],
               '15.04045':['^',20,'none'],
               '15.03473':['x',20,'k'],
               '15.03743':['+',40,'k'],
               '14.06794':['D',15,'none']}
#   'Everest': [2009, 10, 5840, 320, '^', 'None', 30],
#   'West Nepal': [2009, 8, 5590, 138, '*', 'None', 50],
#   'Spiti Lahaul': [2002, 8, 5390, 140, 's', 'None', 25],
#   'Pamir': [2000, 7, 4580, 250, 'v', 'None', 30]}

hd_compare_all_array = None
for nglac, glac_str in enumerate(glaciers):
    hd_obs = pd.read_csv(hd_obs_fp + hd_fn_dict[glac_str])
    reg = int(glac_str.split('.')[0])
    glacno = int(glac_str.split('.')[1])
        
    for i in os.listdir(hd_ts_fp):
        if i.startswith(glac_str):
            mb_df = pd.read_csv(hd_ts_fp + i)
            mb_df.loc[:,:] = mb_df.values.astype(np.float64)
            
            # Bins
            zmin = hd_obs.elev.min()
            zmax = hd_obs.elev.max()
            zbin_width = mb_df.loc[1,'bin_center_elev_m'] - mb_df.loc[0,'bin_center_elev_m']
            zbincenter_min = mb_df.loc[0,'bin_center_elev_m']
            zbincenter_max = mb_df.loc[mb_df.shape[0]-1,'bin_center_elev_m']
            
            # Find minimum bin
            while zbincenter_min - zbin_width / 2 + zbin_width < zmin:
                zbincenter_min += zbin_width
            # Find maximum bin size
            while zbincenter_max - zbin_width /2 > zmax:
                zbincenter_max -= zbin_width
                
            # Statistics for each bin
            for nbin, zbincenter in enumerate(np.arange(zbincenter_min, zbincenter_max + zbin_width/2, zbin_width)):
                elev_idx = np.where((hd_obs.elev.values >= zbincenter - zbin_width/2) & 
                                    (hd_obs.elev.values < zbincenter + zbin_width/2))[0]
                if len(elev_idx) > 0:
                    # Observations
                    hd_obs_subset = hd_obs.loc[elev_idx,'hd_m']
                    hd_bin_mean = hd_obs_subset.mean()
                    hd_bin_std = hd_obs_subset.std()
                    hd_bin_med = np.median(hd_obs_subset)
                    hd_bin_mad = np.median(abs(hd_obs_subset - np.median(hd_obs_subset)))
                    
                    # Model
                    mb_df_idx = np.where(mb_df['bin_center_elev_m'] == zbincenter)[0][0]
                    if 'hd_ts_mean' in mb_df.columns:    
                        mb_df_mean = mb_df.loc[mb_df_idx, 'hd_ts_mean']
                        mb_df_std = mb_df.loc[mb_df_idx, 'hd_ts_std']
                        mb_df_med = mb_df.loc[mb_df_idx, 'hd_ts_med']
                        mb_df_mad = mb_df.loc[mb_df_idx, 'hd_ts_mad']
                    else:
                        mb_df_mean = mb_df.loc[mb_df_idx, 'debris_thick_ts_mean_m']
                        mb_df_std = mb_df.loc[mb_df_idx, 'debris_thick_ts_std_m']
                        mb_df_med = mb_df.loc[mb_df_idx, 'debris_thick_ts_med_m']
                        mb_df_mad = mb_df.loc[mb_df_idx, 'debris_thick_ts_mad_m']
                    
                    print(glac_str, zbincenter, 
                          '  ', np.round(hd_bin_mean,2), '+/-', np.round(hd_bin_std,2), 'vs',
                          '  ', np.round(mb_df_mean,2), '+/-', np.round(mb_df_std,2))
                    
                    bin_data = np.array([reg, glacno, zbincenter, hd_bin_mean, hd_bin_std, hd_bin_med, hd_bin_std,
                                         mb_df_mean, mb_df_std, mb_df_med, mb_df_mad]).reshape(1,11)
                    
                    if hd_compare_all_array is None:
                        hd_compare_all_array = bin_data
                    else:
                        hd_compare_all_array = np.concatenate((hd_compare_all_array, bin_data))
                        
hd_compare_all = pd.DataFrame(hd_compare_all_array, columns=['region', 'glacno', 'zbin', 
                                                             'hd_obs_mean', 'hd_obs_std', 'hd_obs_med', 'hd_obs_mad',
                                                             'hd_ts_mean', 'hd_ts_std', 'hd_ts_med', 'hd_ts_mad'])
        
#%% 
# ===== Plot comparison =====    
#cn_center_obs = 'hd_obs_mean'
#cn_center_ts = 'hd_ts_mean'
#cn_spread_obs = 'hd_obs_std'
#cn_spread_ts = 'hd_ts_std'
cn_center_obs = 'hd_obs_med'
cn_center_ts = 'hd_ts_med'
cn_spread_obs = 'hd_obs_mad'
cn_spread_ts = 'hd_ts_mad'

fig, ax = plt.subplots(1, 1, squeeze=False, figsize=(10,8), gridspec_kw = {'wspace':0, 'hspace':0})
glac_str_leg = []
for ndata, region in enumerate(hd_compare_all.region):
    glac_str = str(int(region)) + '.' + str(int(hd_compare_all.loc[ndata,'glacno'])).zfill(5)
    
    # Legend
    if glac_str not in glac_str_leg:
        glac_str_leg.append(glac_str)
        leg_label = glac_name_dict[glac_str]
    else:
        leg_label = ""
        
    ax[0,0].scatter(hd_compare_all.loc[ndata,cn_center_obs], hd_compare_all.loc[ndata,cn_center_ts], 
                    color='k', 
                    marker=glac_symbol[glac_str][0],
                    facecolor=glac_symbol[glac_str][2], 
                    linewidth=1,
                    s=glac_symbol[glac_str][1],
                    label=leg_label,
                    zorder=3)
    ax[0,0].errorbar(hd_compare_all.loc[ndata,cn_center_obs], hd_compare_all.loc[ndata,cn_center_ts], 
                     xerr=hd_compare_all.loc[ndata,cn_spread_obs], 
                     yerr=hd_compare_all.loc[ndata,cn_spread_ts], 
                     capsize=1, linewidth=0.5, 
                     color='darkgrey', 
                     zorder=2)    
ax[0,0].set_xlabel('Observed $h_d$ (m)', size=12)    
ax[0,0].set_ylabel('Modeled $h_d$ (m)', size=12)
ymin = 0
ymax = 2.
xmin = 0
xmax = 2.
ax[0,0].set_xlim(xmin,xmax)
ax[0,0].set_ylim(ymin,ymax)
ax[0,0].plot([np.min([xmin,ymin]),np.max([xmax,ymax])], [np.min([xmin,ymin]),np.max([xmax,ymax])], color='k', 
             linewidth=0.5, zorder=1)
# Ensure proper order for legend
handles, labels = ax[0,0].get_legend_handles_labels()
labels, handles = zip(*sorted(zip(labels, handles), key=lambda t:t[0]))
ax[0,0].legend(handles, labels, loc=(0.62,0.05), ncol=1, fontsize=10, frameon=True, handlelength=1, 
               handletextpad=0.15, columnspacing=0.5, borderpad=0, labelspacing=0)
fig.set_size_inches(3.45,3.45)
fig.savefig(hd_obs_fp + 'hd_bin_validation.png', bbox_inches='tight', dpi=300)
    

#%%
# ===== Individual comparisons =====
# ------ Kennicott -------
glac_str = '1.15645'
reg = int(glac_str.split('.')[0])
glacno = int(glac_str.split('.')[1])
hd_compare_all_subset = hd_compare_all[(hd_compare_all.region == reg) & (hd_compare_all.glacno == glacno)]
hd_compare_all_subset.reset_index(inplace=True, drop=True)
if glac_str in glaciers:
    fig, ax = plt.subplots(1, 1, squeeze=False, figsize=(10,8), gridspec_kw = {'wspace':0, 'hspace':0})
    for ndata, region in enumerate(hd_compare_all_subset.region):
        ax[0,0].scatter(hd_compare_all_subset.loc[ndata,cn_center_obs], hd_compare_all_subset.loc[ndata,cn_center_ts], 
                        color='k', marker='o', facecolor='none', s=30, zorder=3)
        ax[0,0].errorbar(hd_compare_all_subset.loc[ndata,cn_center_obs], hd_compare_all_subset.loc[ndata,cn_center_ts], 
                         xerr=hd_compare_all_subset.loc[ndata,cn_spread_obs], 
                         yerr=hd_compare_all_subset.loc[ndata,cn_spread_ts], 
                         capsize=1, linewidth=0.5, color='darkgrey', zorder=2)    
    ax[0,0].set_xlabel('Observed $h_d$ (m)', size=12)    
    ax[0,0].set_ylabel('Modeled $h_d$ (m)', size=12)
    ymin, ymax = 0, 0.5
    xmin, xmax = 0, 0.5
    ax[0,0].set_xlim(xmin,xmax)
    ax[0,0].set_ylim(ymin,ymax)
    ax[0,0].plot([np.min([xmin,ymin]),np.max([xmax,ymax])], [np.min([xmin,ymin]),np.max([xmax,ymax])], color='k', 
                 linewidth=0.5, zorder=1)
    fig.text(0.5, 0.9, glac_name_dict[glac_str] + ' (' + glac_str + ')', va='bottom', ha='center', size=12)
    fig.set_size_inches(3.45,3.45)
    fig.savefig(hd_obs_fp + glac_str + '-hd_bin_validation.png', bbox_inches='tight', dpi=300)
    
    
# ----- Lirung ------
glac_str = '15.04045'
reg = int(glac_str.split('.')[0])
glacno = int(glac_str.split('.')[1])
hd_compare_all_subset = hd_compare_all[(hd_compare_all.region == reg) & (hd_compare_all.glacno == glacno)]
hd_compare_all_subset.reset_index(inplace=True, drop=True)
if glac_str in glaciers:
    fig, ax = plt.subplots(1, 1, squeeze=False, figsize=(10,8), gridspec_kw = {'wspace':0, 'hspace':0})
    for ndata, region in enumerate(hd_compare_all_subset.region):
        ax[0,0].scatter(hd_compare_all_subset.loc[ndata,cn_center_obs], hd_compare_all_subset.loc[ndata,cn_center_ts], 
                        color='k', marker='o', facecolor='none', s=30, zorder=3)
        ax[0,0].errorbar(hd_compare_all_subset.loc[ndata,cn_center_obs], hd_compare_all_subset.loc[ndata,cn_center_ts], 
                         xerr=hd_compare_all_subset.loc[ndata,cn_spread_obs], 
                         yerr=hd_compare_all_subset.loc[ndata,cn_spread_ts], 
                         capsize=1, linewidth=0.5, color='darkgrey', zorder=2)    
    ax[0,0].set_xlabel('Observed $h_d$ (m)', size=12)    
    ax[0,0].set_ylabel('Modeled $h_d$ (m)', size=12)
    ymin, ymax = 0, 2
    xmin, xmax = 0, 2
    ax[0,0].set_xlim(xmin,xmax)
    ax[0,0].set_ylim(ymin,ymax)
    ax[0,0].plot([np.min([xmin,ymin]),np.max([xmax,ymax])], [np.min([xmin,ymin]),np.max([xmax,ymax])], color='k', 
                 linewidth=0.5, zorder=1)
    fig.text(0.5, 0.9, glac_name_dict[glac_str] + ' (' + glac_str + ')', va='bottom', ha='center', size=12)
    fig.set_size_inches(3.45,3.45)
    fig.savefig(hd_obs_fp + glac_str + '-hd_bin_validation.png', bbox_inches='tight', dpi=300)

# ----- Ngozumpa ------
glac_str = '15.03473'
reg = int(glac_str.split('.')[0])
glacno = int(glac_str.split('.')[1])
hd_compare_all_subset = hd_compare_all[(hd_compare_all.region == reg) & (hd_compare_all.glacno == glacno)]
hd_compare_all_subset.reset_index(inplace=True, drop=True)
if glac_str in glaciers:
    fig, ax = plt.subplots(1, 1, squeeze=False, figsize=(10,8), gridspec_kw = {'wspace':0, 'hspace':0})
    for ndata, region in enumerate(hd_compare_all_subset.region):
        ax[0,0].scatter(hd_compare_all_subset.loc[ndata,cn_center_obs], hd_compare_all_subset.loc[ndata,cn_center_ts], 
                        color='k', marker='o', facecolor='none', s=30, zorder=3)
        ax[0,0].errorbar(hd_compare_all_subset.loc[ndata,cn_center_obs], hd_compare_all_subset.loc[ndata,cn_center_ts], 
                         xerr=hd_compare_all_subset.loc[ndata,cn_spread_obs], 
                         yerr=hd_compare_all_subset.loc[ndata,cn_spread_ts], 
                         capsize=1, linewidth=0.5, color='darkgrey', zorder=2)    
    ax[0,0].set_xlabel('Observed $h_d$ (m)', size=12)    
    ax[0,0].set_ylabel('Modeled $h_d$ (m)', size=12)
    ymin, ymax = 0, 4
    xmin, xmax = 0, 4
    ax[0,0].set_xlim(xmin,xmax)
    ax[0,0].set_ylim(ymin,ymax)
    ax[0,0].plot([np.min([xmin,ymin]),np.max([xmax,ymax])], [np.min([xmin,ymin]),np.max([xmax,ymax])], color='k', 
                 linewidth=0.5, zorder=1)
    fig.text(0.5, 0.9, glac_name_dict[glac_str] + ' (' + glac_str + ')', va='bottom', ha='center', size=12)
    fig.set_size_inches(3.45,3.45)
    fig.savefig(hd_obs_fp + glac_str + '-hd_bin_validation.png', bbox_inches='tight', dpi=300)

# ----- Imja-Lhotse Shar ------
glac_str = '15.03743'
reg = int(glac_str.split('.')[0])
glacno = int(glac_str.split('.')[1])
hd_compare_all_subset = hd_compare_all[(hd_compare_all.region == reg) & (hd_compare_all.glacno == glacno)]
hd_compare_all_subset.reset_index(inplace=True, drop=True)
if glac_str in glaciers:
    fig, ax = plt.subplots(1, 1, squeeze=False, figsize=(10,8), gridspec_kw = {'wspace':0, 'hspace':0})
    for ndata, region in enumerate(hd_compare_all_subset.region):
        ax[0,0].scatter(hd_compare_all_subset.loc[ndata,cn_center_obs], hd_compare_all_subset.loc[ndata,cn_center_ts], 
                        color='k', marker='o', facecolor='none', s=30, zorder=3)
        ax[0,0].errorbar(hd_compare_all_subset.loc[ndata,cn_center_obs], hd_compare_all_subset.loc[ndata,cn_center_ts], 
                         xerr=hd_compare_all_subset.loc[ndata,cn_spread_obs], 
                         yerr=hd_compare_all_subset.loc[ndata,cn_spread_ts], 
                         capsize=1, linewidth=0.5, color='darkgrey', zorder=2)    
    ax[0,0].set_xlabel('Observed $h_d$ (m)', size=12)    
    ax[0,0].set_ylabel('Modeled $h_d$ (m)', size=12)
    ymin, ymax = 0, 1.5
    xmin, xmax = 0, 1.5
    ax[0,0].set_xlim(xmin,xmax)
    ax[0,0].set_ylim(ymin,ymax)
    ax[0,0].plot([np.min([xmin,ymin]),np.max([xmax,ymax])], [np.min([xmin,ymin]),np.max([xmax,ymax])], color='k', 
                 linewidth=0.5, zorder=1)
    fig.text(0.5, 0.9, glac_name_dict[glac_str] + ' (' + glac_str + ')', va='bottom', ha='center', size=12)
    fig.set_size_inches(3.45,3.45)
    fig.savefig(hd_obs_fp + glac_str + '-hd_bin_validation.png', bbox_inches='tight', dpi=300)

    
# ----- Imja-Lhotse Shar ------
glac_str = '14.06794'
reg = int(glac_str.split('.')[0])
glacno = int(glac_str.split('.')[1])
hd_compare_all_subset = hd_compare_all[(hd_compare_all.region == reg) & (hd_compare_all.glacno == glacno)]
hd_compare_all_subset.reset_index(inplace=True, drop=True)
if glac_str in glaciers:
    fig, ax = plt.subplots(1, 1, squeeze=False, figsize=(10,8), gridspec_kw = {'wspace':0, 'hspace':0})
    for ndata, region in enumerate(hd_compare_all_subset.region):
        ax[0,0].scatter(hd_compare_all_subset.loc[ndata,cn_center_obs], hd_compare_all_subset.loc[ndata,cn_center_ts], 
                        color='k', marker='o', facecolor='none', s=30, zorder=3)
        ax[0,0].errorbar(hd_compare_all_subset.loc[ndata,cn_center_obs], hd_compare_all_subset.loc[ndata,cn_center_ts], 
                         xerr=hd_compare_all_subset.loc[ndata,cn_spread_obs], 
                         yerr=hd_compare_all_subset.loc[ndata,cn_spread_ts], 
                         capsize=1, linewidth=0.5, color='darkgrey', zorder=2)    
    ax[0,0].set_xlabel('Observed $h_d$ (m)', size=12)    
    ax[0,0].set_ylabel('Modeled $h_d$ (m)', size=12)
    ymin, ymax = 0, 0.5
    xmin, xmax = 0, 0.5
    ax[0,0].set_xlim(xmin,xmax)
    ax[0,0].set_ylim(ymin,ymax)
    ax[0,0].plot([np.min([xmin,ymin]),np.max([xmax,ymax])], [np.min([xmin,ymin]),np.max([xmax,ymax])], color='k', 
                 linewidth=0.5, zorder=1)
    fig.text(0.5, 0.9, glac_name_dict[glac_str] + ' (' + glac_str + ')', va='bottom', ha='center', size=12)
    fig.set_size_inches(3.45,3.45)
    fig.savefig(hd_obs_fp + glac_str + '-hd_bin_validation.png', bbox_inches='tight', dpi=300)
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

    
        