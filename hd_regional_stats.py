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
from matplotlib.ticker import MultipleLocator
import numpy as np
import pandas as pd
#from scipy.optimize import curve_fit
#from scipy.stats import linregress
from scipy.stats import median_absolute_deviation
from scipy.stats import linregress
#import xarray as xr

import debrisglobal.globaldebris_input as debris_prms

option_regional_stats = False
option_regional_plots = True


#%%% ===== SCRIPT OPTIONS =====
#if not os.path.exists(melt_compare_fp):
#    os.makedirs(melt_compare_fp)
    
if option_regional_stats:
    #rois = ['01','02','03','04','05','06','07','08','09','10','11','12','HMA','16','17','18']
    rois = ['01','02','03','04','05','06','07','08','09','11','12','13','14', '15', '16','17','18']
    #rois = ['18']
    
    hd_cn = 'hd_ts_mean_m'
    mf_cn = 'mf_ts_mean'
    
    reg_stats_cns = ['roi', 'hd_mean', 'hd_std', 'hd_med', 'hd_mad', 'hd_25', 'hd_75', 
                     'mf_mean', 'mf_std', 'mf_med', 'mf_mad', 'mf_25', 'mf_75']
    reg_stats_df = pd.DataFrame(np.zeros((len(rois)+1,len(reg_stats_cns))), columns=reg_stats_cns)
    
    ## ===== REGIONAL MELT FACTOR STATISTICS =====
    hdts_bin_fp = debris_prms.mb_binned_fp + '_wdebris_hdts/'
    hdts_bin_fp_extrap = debris_prms.mb_binned_fp + '_wdebris_hdts_extrap/'
    
    hd_list_all = []
    mf_list_all = []
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
        
          
        area_reg_km2 = 0
        hd_list = []
        mf_list = []
        for nfn, fullfn in enumerate(glac_hd_fullfns):
    #    for nfn, fullfn in enumerate([glac_hd_fullfns[0]]):
            if nfn%500 == 0:
                print('  ', nfn)
            df = pd.read_csv(fullfn)
            
        #     print(df.loc[:,['bin_center_elev_m', ' z1_bin_area_valid_km2', ' dc_bin_area_valid_km2', 
        #                     ' dc_bin_area_perc', 'debris_perc', 'hd_ts_med', 'mf_ts_med']])
            
            if 'hd_ts_mean_m' in list(df.columns):
                # Can switch 1e5 to 1e6 for m accuracy, right now its 10 m, so 3x3 m pixel
                area_factor = 1e4
                n_values = np.round(df['dc_bin_area_valid_km2'].values*area_factor,0).astype(int)
                area_reg_km2 += n_values.sum() /area_factor
                hd_values = df[hd_cn].values
                mf_values = df[mf_cn].values
                for nidx, n_value in enumerate(n_values):
                    if n_value > 0 and not np.isnan(hd_values[nidx]):
                        hd_list.extend(np.repeat(hd_values[nidx],n_values[nidx]))
                        mf_list.extend(np.repeat(mf_values[nidx],n_values[nidx]))
        
        # Record stats
        hd_array = np.array(hd_list)
        hd_mean = hd_array.mean()
        hd_std = hd_array.std()
        hd_med = np.median(hd_array)
        hd_mad = median_absolute_deviation(hd_array)
        hd_25 = np.percentile(hd_array, 25)
        hd_75 = np.percentile(hd_array, 75)
        
        mf_array = np.array(mf_list)
        mf_mean = mf_array.mean()
        mf_std = mf_array.std()
        mf_med = np.median(mf_array)
        mf_mad = median_absolute_deviation(mf_array)
        mf_25 = np.percentile(mf_array, 25)
        mf_75 = np.percentile(mf_array, 75)
        
        reg_stats_df.loc[nroi,:] = [roi, hd_mean, hd_std, hd_med, hd_mad, hd_25, hd_75,
                                    mf_mean, mf_std, mf_med, mf_mad, mf_25, mf_75]
        
        print(roi, 'hd:', np.round(hd_med,2), '(' +  str(np.round(hd_25,2)) + ' - ', str(np.round(hd_75,2)) + ')', 
              '  mf:', np.round(mf_med,2), '(' +  str(np.round(mf_25,2)) + ' - ', str(np.round(mf_75,2)) + ')')
        
        #%%
        # ===== HISTOGRAM: regional debris thickness ======
        hd_bins = np.arange(0,2.01,0.05)
        label_frequency = 10
        hist_fn = roi + '_hd_hist.png'
        hist, bin_edges = np.histogram(hd_array,hd_bins) # make the histogram
        hist_area_km2 = hist / area_factor / area_reg_km2
        fig,ax = plt.subplots()    
        # Plot the histogram heights against integers on the x axis
        ax.bar(range(len(hist)),hist_area_km2, width=1, edgecolor='k', facecolor='limegreen', linewidth=0.5) 
        ax.set_xticks([i-0.5 for i,j in enumerate(hd_bins)], minor=True)
        bin_idx = np.arange(0,len(hd_bins),label_frequency)
        ax.set_xticks(bin_idx-0.5)
        ax.set_xticklabels([str(np.round(x,2)) for x in hd_bins[bin_idx]], rotation=0, ha='center')
        ax.set_xlabel('$h_{d}$ (m)', fontsize=12)
        ax.set_xlim(-0.5,len(hd_bins)-1.5)
        ax.set_ylim(0,0.2)
        ax.yaxis.set_major_locator(MultipleLocator(0.05))
        ax.yaxis.set_minor_locator(MultipleLocator(0.01))
        ax.set_ylabel('Area (-)', fontsize=12)
        text_str = str(int(np.round(area_reg_km2,0)))
        ax.text(0.98, 0.95, text_str, size=10, 
                horizontalalignment='right', verticalalignment='top', transform=ax.transAxes)
        # Save figure
        fig_fp = debris_prms.output_fp + 'figures/histograms/'
        if not os.path.exists(fig_fp):
            os.makedirs(fig_fp)
        fig.set_size_inches(3,2)
        fig.savefig(fig_fp + hist_fn, bbox_inches='tight', dpi=300, transparent=True)
        plt.close()
        
        # ===== HISTOGRAM: regional melt factors ======
        mf_bins = np.arange(0,1.51,0.05)
        label_frequency = 10
        hist_fn = roi + '_mf_hist.png'
        hist, bin_edges = np.histogram(mf_array,mf_bins) # make the histogram
        hist_area_km2 = hist / area_factor / area_reg_km2
        fig,ax = plt.subplots()    
        ax.bar(range(len(hist)),hist_area_km2, width=1, edgecolor='k', facecolor='deepskyblue', linewidth=0.5) 
        ax.set_xticks([i-0.5 for i,j in enumerate(hd_bins)], minor=True)
        bin_idx = np.arange(0,len(hd_bins),label_frequency)
        ax.set_xticks(bin_idx-0.5)
        ax.set_xticklabels([str(np.round(x,2)) for x in hd_bins[bin_idx]], rotation=0, ha='center')    
        ax.set_xlabel('Melt factor (-)', fontsize=12)
        ax.set_xlim(-0.5,len(mf_bins)-1.5)
        ax.set_ylim(0,0.125)
        ax.yaxis.set_major_locator(MultipleLocator(0.05))
        ax.yaxis.set_minor_locator(MultipleLocator(0.01))
        ax.set_ylabel('Area (-)', fontsize=12)
    #    text_str = str(int(np.round(area_reg_km2,0))) + ' $\mathregular{km^{2}}$'
    #    text_str = str(int(np.round(area_reg_km2,0)))
    #    ax.text(0.98, 0.95, text_str, size=10, 
    #            horizontalalignment='right', verticalalignment='top', transform=ax.transAxes)
        # Save figure
        fig_fp = debris_prms.output_fp + 'figures/histograms/'
        if not os.path.exists(fig_fp):
            os.makedirs(fig_fp)
        fig.set_size_inches(3,2)
        fig.savefig(fig_fp + hist_fn, bbox_inches='tight', dpi=300, transparent=True)
        plt.close()
            
        
        #%%
        
        hd_list_all.extend(hd_list)
        mf_list_all.extend(mf_list)
        
    # All regions statistic
    hd_array_all = np.array(hd_list_all)
    hd_mean_all = hd_array_all.mean()
    hd_std_all = hd_array_all.std()
    hd_med_all = np.median(hd_array_all)
    hd_mad_all = median_absolute_deviation(hd_array_all)
    hd_25_all = np.percentile(hd_array_all, 25)
    hd_75_all = np.percentile(hd_array_all, 75)
    
    mf_array_all = np.array(mf_list_all)
    mf_mean_all = mf_array_all.mean()
    mf_std_all = mf_array_all.std()
    mf_med_all = np.median(mf_array_all)
    mf_mad_all = median_absolute_deviation(mf_array_all)
    mf_25_all = np.percentile(mf_array_all, 25)
    mf_75_all = np.percentile(mf_array_all, 75)
    
    nroi += 1
    reg_stats_df.loc[nroi,:] = ['all', hd_mean_all, hd_std_all, hd_med_all, hd_mad_all, hd_25_all, hd_75_all,
                                mf_mean_all, mf_std_all, mf_med_all, mf_mad_all, mf_25_all, mf_75_all]
        
    # Export regional statistics
    reg_stats_df.to_csv(debris_prms.output_fp + 'reg_stats_hd_mf.csv', index=False)


#%%
if option_regional_plots:
    reg_stats_df = pd.read_csv(debris_prms.output_fp + 'reg_stats_hd_mf_wdc.csv')
    reg_stats_df.dropna(subset=['hd_med'], inplace=True)
    
    # ====== % debris cover vs. median debris thickness ===== 
    fig, ax = plt.subplots(1, 1, squeeze=False, sharex=False, sharey=False, gridspec_kw = {'wspace':0.35, 'hspace':0})
    ax[0,0].scatter(reg_stats_df['% DC'], reg_stats_df['hd_med'], color='k', marker='o', facecolor='none', s=30,
                    clip_on=False)
    
    slope, intercept, r_value, p_value, std_err = linregress(reg_stats_df['% DC'].values, reg_stats_df['hd_med'].values)
    print('r = ' + str(np.round(r_value,2)))
    
    lobf_x = np.arange(reg_stats_df['% DC'].min(),reg_stats_df['% DC'].max()+0.1,0.1)
    lobf_y = intercept + slope * lobf_x
    ax[0,0].plot(lobf_x, lobf_y, color='k', linewidth=1)
    
    # X-label
    ax[0,0].set_xlabel('Relative debris cover (%)', size=12)
    ax[0,0].set_xlim(0,30)
    ax[0,0].xaxis.set_major_locator(plt.MultipleLocator(10))
    ax[0,0].xaxis.set_minor_locator(plt.MultipleLocator(2))  
    # Y-label
    ax[0,0].set_ylabel('$h_{d}$ median (m)', size=12)
    ax[0,0].set_ylim(0,0.29)
    ax[0,0].yaxis.set_major_locator(plt.MultipleLocator(0.1))
    ax[0,0].yaxis.set_minor_locator(plt.MultipleLocator(0.02))
    # Tick parameters
    #ax[0,0].yaxis.set_ticks_position('both')
    ax[0,0].tick_params(axis='both', which='major', labelsize=10, direction='inout')
    ax[0,0].tick_params(axis='both', which='minor', labelsize=8, direction='in') 
    
    # Save plot
    fig.set_size_inches(3, 3)
    fig_fn = 'hd_med_vs_dc%.png'
    fig.savefig(debris_prms.output_fp + fig_fn, bbox_inches='tight', dpi=300, transparent=True)
    
    
    
    #%% ====== % debris cover vs. median debris thickness ===== 
    fig, ax = plt.subplots(1, 1, squeeze=False, sharex=False, sharey=False, gridspec_kw = {'wspace':0.35, 'hspace':0})
    ax[0,0].scatter(reg_stats_df['hd_med'], reg_stats_df['mf_med'], color='k', marker='o', facecolor='none', s=30,
                    clip_on=False)
    
    slope, intercept, r_value, p_value, std_err = linregress(reg_stats_df['hd_med'].values, 
                                                             reg_stats_df['mf_med'].values)
    print('r = ' + str(np.round(r_value,2)))
    
    lobf_x = np.arange(reg_stats_df['hd_med'].min(),reg_stats_df['mf_med'].max()+0.1,0.1)
    lobf_y = intercept + slope * lobf_x
    ax[0,0].plot(lobf_x, lobf_y, color='k', linewidth=1)
    
    # X-label
    ax[0,0].set_xlabel('$h_{d}$ median (m)', size=12)
    ax[0,0].set_xlim(0,0.29)
    ax[0,0].xaxis.set_major_locator(plt.MultipleLocator(0.1))
    ax[0,0].xaxis.set_minor_locator(plt.MultipleLocator(0.02))  
    # Y-label
    ax[0,0].set_ylabel('Melt factor (-)', size=12)
    ax[0,0].set_ylim(0,1)
    ax[0,0].yaxis.set_major_locator(plt.MultipleLocator(0.2))
    ax[0,0].yaxis.set_minor_locator(plt.MultipleLocator(0.05))
    # Tick parameters
    #ax[0,0].yaxis.set_ticks_position('both')
    ax[0,0].tick_params(axis='both', which='major', labelsize=10, direction='inout')
    ax[0,0].tick_params(axis='both', which='minor', labelsize=8, direction='in') 
    
    # Save plot
    fig.set_size_inches(3, 3)
    fig_fn = 'hd_med_vs_mf_med.png'
    fig.savefig(debris_prms.output_fp + fig_fn, bbox_inches='tight', dpi=300, transparent=True)
