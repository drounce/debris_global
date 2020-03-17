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

#%%% ===== SCRIPT OPTIONS =====
option_melt_comparison = False
option_hd_comparison = False
option_hd_spatial_compare = True

melt_compare_fp = debris_prms.main_directory + '/../hd_obs/figures/hd_melt_compare/'
hd_compare_fp = debris_prms.main_directory + '/../hd_obs/figures/hd_obs_compare/'

glac_name_dict = {'1.15645':'Kennicott',
                  '11.01604': 'Suldenferner',
                  '15.04045':'Lirung',
                  '15.03473':'Ngozumpa',
                  '15.03743':'Imja',
                  '14.06794':'Baltoro'}
hd_obs_fp = debris_prms.main_directory + '/../hd_obs/'
hd_ds_dict = {'1.15645': ['anderson2019'],
              '11.01604': ['nicholson-suld']}


#    hd_fn_dict = {
#                  '15.04045': [debris_prms.main_directory + '/../hd_obs/lirung_nicholson_gpr.csv'],
#                  '15.03473': [debris_prms.main_directory + '/../hd_obs/ngoz_mccarthy_gpr.csv'],
#                  '15.03743': [debris_prms.main_directory + '/../hd_obs/imja_rounce2014.csv'],
#                  '14.06794': [debris_prms.main_directory + '/../hd_obs/baltoro_mihalcea2006.csv']}


hd_ds_fn_dict = {
        'anderson2019': hd_obs_fp + 'kennicott_anderson_2019.csv',
        'nicholson-suld': hd_obs_fp + 'Nicholson_datasets/dz_SDF_all_whd_ts.csv',
        'mccarthy-ngoz': hd_obs_fp + 'ngoz_mccarthy_gpr_whdts.csv'} 

if os.path.exists(melt_compare_fp) == False:
    os.makedirs(melt_compare_fp)
if os.path.exists(hd_compare_fp) == False:
    os.makedirs(hd_compare_fp)

#%% ===== FUNCTIONS =====
def plot_hd_vs_melt_comparison(measured_hd, measured_melt, glac_name, fig_fn, melt_fn, start_yearfrac,
                               hd_min=0, hd_max=2, hd_tick_major=0.25, hd_tick_minor=0.05,
                               melt_min=0, melt_max=70, melt_tick_major=10, melt_tick_minor=5):
    """ Plot comparison of debris vs. melt for various sites """
    # Dataset of melt data
    melt_fp = debris_prms.ostrem_fp
    ds_ostrem = xr.open_dataset(melt_fp + melt_fn)
    
    # ===== Ostrem Curve =====
    time_year = pd.to_datetime(ds_ostrem.time.values).year
    time_daysperyear = np.array([366 if x%4 == 0 else 365 for x in time_year])
    time_yearfrac = time_year + (pd.to_datetime(ds_ostrem.time.values).dayofyear-1) / time_daysperyear

    start_idx = np.where(abs(time_yearfrac - start_yearfrac) == abs(time_yearfrac - start_yearfrac).min())[0][0]
    end_idx = np.where(abs(time_yearfrac - end_yearfrac) == abs(time_yearfrac - end_yearfrac).min())[0][0]

    # Debris thickness
    debris_thicknesses = ds_ostrem.hd_cm.values
    debris_melt_df = pd.DataFrame(np.zeros((len(debris_thicknesses),2)), 
                                  columns=['debris_thickness', 'melt_mwea'])  

    nelev = 0
    for ndebris, debris_thickness in enumerate(debris_thicknesses):    
        # Units: mm w.e. per day                
        melt_mmwed = (ds_ostrem['melt'][ndebris,start_idx:end_idx,nelev].values.sum() 
                      * 1000 / len(time_yearfrac[start_idx:end_idx]))
        debris_melt_df.loc[ndebris] = debris_thickness / 100, melt_mmwed
        
    # Fit curve
    fit_idx = list(np.where(debris_thicknesses >= 5)[0])            
    func_coeff, pcov = curve_fit(melt_fromdebris_func, 
                                 debris_melt_df.debris_thickness.values[fit_idx], 
                                 debris_melt_df.melt_mwea.values[fit_idx])
#    melt_cleanice = debris_melt_df.loc[0,'melt_mwea']
#    idx_2cm = np.where(debris_thicknesses == 2)[0][0]
#    melt_2cm = debris_melt_df.loc[idx_2cm, 'melt_mwea']
            
    # ===== PLOT DEBRIS VS. SURFACE LOWERING ===== 
    fig, ax = plt.subplots(1, 1, squeeze=False, sharex=False, sharey=False, 
                          gridspec_kw = {'wspace':0.4, 'hspace':0.15})

    # Fitted curve
    debris_4curve = np.arange(0.02,5.01,0.01)
    melt_4curve = melt_fromdebris_func(debris_4curve, func_coeff[0], func_coeff[1])

    # Plot curve
    ax[0,0].plot(measured_hd, measured_melt, 'D', 
                 color='k', markersize=3, markerfacecolor="None", markeredgewidth=0.75, zorder=1)
    ax[0,0].plot(debris_4curve, melt_4curve, 
                 color='k', linewidth=1, linestyle='--', zorder=2)
    # text
    ax[0,0].text(0.5, 1.05, glac_name, size=10, horizontalalignment='center', verticalalignment='top', 
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
    ax[0,0].set_xlim(hd_min, hd_max)
    ax[0,0].xaxis.set_tick_params(labelsize=12)
    ax[0,0].xaxis.set_major_locator(plt.MultipleLocator(hd_tick_major))
    ax[0,0].xaxis.set_minor_locator(plt.MultipleLocator(hd_tick_minor))  
    # Y-label
    ax[0,0].set_ylabel('Melt (mm w.e. d$^{-1}$)', size=12)
    ax[0,0].set_ylim(melt_min, melt_max)
    ax[0,0].yaxis.set_major_locator(plt.MultipleLocator(melt_tick_major))
    ax[0,0].yaxis.set_minor_locator(plt.MultipleLocator(melt_tick_minor))
    # Tick parameters
    ax[0,0].yaxis.set_ticks_position('both')
    ax[0,0].tick_params(axis='both', which='major', labelsize=12, direction='inout')
    ax[0,0].tick_params(axis='both', which='minor', labelsize=10, direction='in') 
    # Save plot
    fig.set_size_inches(4, 4)
    fig.savefig(melt_compare_fp + fig_fn, bbox_inches='tight', dpi=300)
    plt.close()


#%%
if option_melt_comparison:
    
    compare_khumbu = True
    compare_hailuogou = True
    compare_batal = True
    compare_miage = True
    
    # ===== KHUMBU ====
    if compare_khumbu:
        print('\nOstrem curve comparison with Kayastha et al (2000)\n')
        # Data from Kayastha et al. (2000)
        # units: m
        measured_hd = np.array([0, 0.3, 2, 5, 10, 20, 30, 40]) / 100
        # units: mm w.e. d-1
        measured_melt =  np.array([27.9, 57.6, 42.3, 29.7, 18, 12.6, 10.8, 9])

        # Plot and filename details
        glac_name = 'Khumbu Glacier (15.03733)'
        fig_fn = '15.03733_hd_melt_Kay2000.png'
        melt_fn = '2800N-8700E-debris_melt_curve.nc'
        
        # Manually estimate indices
        start_yearfrac = 2000 + 143/366
        end_yearfrac = 2000 + 153/366
        
        hd_min = 0
        hd_max = 0.45
        hd_tick_major = 0.1
        hd_tick_minor = 0.05
        
        melt_min = 0
        melt_max = 70
        melt_tick_major = 10
        melt_tick_minor = 5
        
        plot_hd_vs_melt_comparison(measured_hd, measured_melt, glac_name, fig_fn, melt_fn, start_yearfrac,
                                   hd_min=hd_min, hd_max=hd_max, 
                                   hd_tick_major=hd_tick_major, hd_tick_minor=hd_tick_minor,
                                   melt_min=melt_min, melt_max=melt_max, 
                                   melt_tick_major=melt_tick_major, melt_tick_minor=melt_tick_minor)
    
    # ===== HAILUOGOU (15.07866) ====
    if compare_hailuogou:
        print('\nOstrem curve comparison with Zhang et al (2011)\n')
        
        # Data from Zhang et al. (2011)
        # units: m
        measured_hd = np.array([2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 5, 5, 6, 7, 7, 10, 10, 11, 13]) / 100
        # units: mm w.e. d-1
        measured_melt =  np.array([65.2, 55.4, 52.8, 51.6, 47.0, 53.4, 44.4, 50.3, 58, 48.9, 58.4, 54.4, 44.8, 
                                   52.6, 43.7, 52.5, 38.5, 36.5, 34.2, 28.4])

        # Plot and filename details
        glac_name = 'Hailuogou Glacier (15.07866)'
        fig_fn = '15.07866_hd_melt_Zhang2011.png'
        melt_fn = '2950N-10175E-debris_melt_curve.nc'
        
        # Manually estimate indices
        start_yearfrac = 2008 + 134/366
        end_yearfrac = 2008 + 274/366
        
        hd_min = 0
        hd_max = 0.2
        hd_tick_major = 0.05
        hd_tick_minor = 0.01
        
        melt_min = 0
        melt_max = 70
        melt_tick_major = 10
        melt_tick_minor = 5
        
        plot_hd_vs_melt_comparison(measured_hd, measured_melt, glac_name, fig_fn, melt_fn, start_yearfrac,
                                   hd_min=hd_min, hd_max=hd_max, 
                                   hd_tick_major=hd_tick_major, hd_tick_minor=hd_tick_minor,
                                   melt_min=melt_min, melt_max=melt_max, 
                                   melt_tick_major=melt_tick_major, melt_tick_minor=melt_tick_minor)
        
        
    # ===== BATAL (14.16042) ====
    if compare_batal:
        print('\nOstrem curve comparison with Patel et al (2016)\n')
        
        # Data from Patel et al (2016)
        # units: m
        measured_hd = np.array([51, 54, 43, 53, 24, 27, 27, 50, 5, 2, 0, 5, 0, 2]) / 100
        # units: mm w.e. d-1
        measured_melt =  np.array([5.75, 6.125, 11.875, 7.875, 8.625, 10.375, 9.75, 7.375, 12.875, 19.875, 
                                   18.125, 14.375, 18.625, 10.75])

        # Plot and filename details
        glac_name = 'Batal Glacier (14.16042)'
        fig_fn = '14.16042_hd_melt_Patel2016.png'
        melt_fn = '3225N-7750E-debris_melt_curve.nc'
        
        # Manually estimate indices
        start_yearfrac = 2014 + 213/365
        end_yearfrac = 2014 + 288/365
        
        hd_min = 0
        hd_max = 0.6
        hd_tick_major = 0.1
        hd_tick_minor = 0.05
        
        melt_min = 0
        melt_max = 50
        melt_tick_major = 10
        melt_tick_minor = 5
        
        plot_hd_vs_melt_comparison(measured_hd, measured_melt, glac_name, fig_fn, melt_fn, start_yearfrac,
                                   hd_min=hd_min, hd_max=hd_max, 
                                   hd_tick_major=hd_tick_major, hd_tick_minor=hd_tick_minor,
                                   melt_min=melt_min, melt_max=melt_max, 
                                   melt_tick_major=melt_tick_major, melt_tick_minor=melt_tick_minor)
    
    # ===== MIAGE (11.03005) ====
    if compare_miage:
        print('\nOstrem curve comparison with Reid and Brock (2010)\n')
        
        # Data from Reid and Brock (2010)
        # units: m
        measured_hd = np.array([55, 26, 21, 19.5, 19.5, 19.5, 17, 16, 16, 16, 15, 12, 11, 10]) / 100
        # units: mm w.e. d-1
        measured_melt =  np.array([7,  17, 16, 15,   19,   22,   17, 19, 20, 22, 20, 20, 22, 21])

        # Plot and filename details
        glac_name = 'Miage Glacier (11.03005)'
        fig_fn = '11.03005_hd_melt_Reid2010.png'
        melt_fn = '4650N-1050E-debris_melt_curve.nc'
        
        # Manually estimate indices
        start_yearfrac = 2005 + 172/365
        end_yearfrac = 2005 + 247/365
        
        hd_min = 0
        hd_max = 0.6
        hd_tick_major = 0.1
        hd_tick_minor = 0.02
        
        melt_min = 0
        melt_max = 50
        melt_tick_major = 10
        melt_tick_minor = 5
        
        plot_hd_vs_melt_comparison(measured_hd, measured_melt, glac_name, fig_fn, melt_fn, start_yearfrac,
                                   hd_min=hd_min, hd_max=hd_max, 
                                   hd_tick_major=hd_tick_major, hd_tick_minor=hd_tick_minor,
                                   melt_min=melt_min, melt_max=melt_max, 
                                   melt_tick_major=melt_tick_major, melt_tick_minor=melt_tick_minor)

#%%
if option_hd_comparison: 
    
    glaciers = ['1.15645']
    glaciers = ['11.01604']
    
    
    #glaciers = ['1.15645', '15.03473', '15.03733', '15.03734', '15.04121', '15.04045', '14.06794', '14.04477', 
    #            '13.43232', '15.07886', '14.16042']
#    glaciers = ['1.15645', '15.04045', '15.03473', '15.03743', '14.06794']
    

    
    glac_symbol = {'1.15645':['o',20,'none'],
                   '15.04045':['^',20,'none'],
                   '15.03473':['x',20,'k'],
                   '15.03743':['+',40,'k'],
                   '14.06794':['D',15,'none']}

    bin_width = 10
    
    # ===== Plot comparison =====    
    #cn_center_obs = 'hd_obs_mean'
    #cn_center_ts = 'hd_ts_mean_m'
    #cn_spread_obs = 'hd_obs_std'
    #cn_spread_ts = 'hd_ts_std_m'
    cn_center_obs = 'hd_obs_med'
    cn_center_ts = 'hd_ts_med_m'
    cn_spread_obs = 'hd_obs_mad'
    cn_spread_ts = 'hd_ts_mad_m'
    
    
    
    #   'Everest': [2009, 10, 5840, 320, '^', 'None', 30],
    #   'West Nepal': [2009, 8, 5590, 138, '*', 'None', 50],
    #   'Spiti Lahaul': [2002, 8, 5390, 140, 's', 'None', 25],
    #   'Pamir': [2000, 7, 4580, 250, 'v', 'None', 30]}
    
    hd_compare_all_array = None
    ds_name_list = []
    for nglac, glac_str in enumerate(glaciers):
        
        for hd_ds in hd_ds_dict[glac_str]:
            hd_obs = pd.read_csv(hd_ds_fn_dict[hd_ds])
        
            reg = int(glac_str.split('.')[0])
            glacno = int(glac_str.split('.')[1])
                
            try:
                hdts_fp = debris_prms.mb_binned_fp_wdebris_hdts
                hdts_fn = glac_str + '_mb_bins_hdts.csv'
                mb_df = pd.read_csv(hdts_fp + hdts_fn)
            except:
                print('\n\nNEED TO ADD FILEPATH FOR EXTRAPOLATED GLACIERS\n\n')
                hdts_fp = debris_prms.mb_binned_fp_wdebris_hdts
                hdts_fn = glac_str + '_mb_bins_hdts.csv'
                mb_df = pd.read_csv(hdts_fp + hdts_fn)

            mb_df.loc[:,:] = mb_df.values.astype(np.float64)
        
            # Bins
            zmin = hd_obs.elev.min()
            zmax = hd_obs.elev.max()
            zbincenter_min = mb_df.loc[0,'bin_center_elev_m']
            zbincenter_max = mb_df.loc[mb_df.shape[0]-1,'bin_center_elev_m']
            
            # Find minimum bin
            while zbincenter_min - bin_width / 2 + bin_width < zmin:
                zbincenter_min += bin_width
            # Find maximum bin size
            while zbincenter_max - bin_width /2 > zmax:
                zbincenter_max -= bin_width
                
            # Statistics for each bin
            for nbin, zbincenter in enumerate(np.arange(zbincenter_min, zbincenter_max + bin_width/2, bin_width)):
    #            print('\n', zbincenter)
                zbincenter_min = zbincenter - bin_width/2
                zbincenter_max = zbincenter + bin_width/2
                elev_idx_obs = np.where((hd_obs['elev'].values >= zbincenter_min) & 
                                        (hd_obs['elev'].values < zbincenter_max))[0]
                if len(elev_idx_obs) > 1:
                    # Observations
                    hd_obs_subset = hd_obs.loc[elev_idx_obs,'hd_m']
                    hd_bin_mean = hd_obs_subset.mean()
                    hd_bin_std = hd_obs_subset.std()
                    hd_bin_med = np.median(hd_obs_subset)
                    hd_bin_mad = np.median(abs(hd_obs_subset - np.median(hd_obs_subset)))
                    
                    # Model
                    mb_df_elevs = mb_df['bin_center_elev_m'].values
                    mb_df_idxs = np.where((mb_df_elevs >= zbincenter_min) &
                                          (mb_df_elevs < zbincenter_max))[0]
                    
    #                print('  ', mb_df_idxs)
                    hd_list = []
                    for mb_df_idx in mb_df_idxs:
                        bin_hd_center = mb_df.loc[mb_df_idx,cn_center_ts]
                        bin_hd_spread = mb_df.loc[mb_df_idx,cn_spread_ts]
                        bin_hd_count = int(mb_df.loc[mb_df_idx,'dc_bin_count_valid'])
                        
                        # Randomly create thicknesses based on center and spread, but ensure mean is similar
                        hd_list_single = np.random.normal(loc=bin_hd_center, scale=bin_hd_spread, size=(bin_hd_count))
                        while (abs(np.mean(hd_list_single) - bin_hd_center) > 0.005 or 
                               abs(np.std(hd_list_single) - bin_hd_spread) > 0.01):
                            hd_list_single = np.random.normal(loc=bin_hd_center, scale=bin_hd_spread, size=(bin_hd_count))
                        hd_list.extend(hd_list_single)
                        
    #                    print(bin_hd_center, bin_hd_spread, bin_hd_count)
    #                    print(np.mean(hd_list_single), np.std(hd_list_single), len(hd_list_single))
                        
                    hd_array = np.array(hd_list)
                    mb_df_mean = hd_array.mean()
                    mb_df_std = hd_array.std()
                    mb_df_med = np.median(hd_array)
                    mb_df_mad = median_absolute_deviation(hd_array)
                    
                    print(glac_str, zbincenter, 
                          '  ', np.round(hd_bin_mean,2), '+/-', np.round(hd_bin_std,2), 'vs',
                          '  ', np.round(mb_df_mean,2), '+/-', np.round(mb_df_std,2))
                    
                    bin_data = np.array([reg, glacno, zbincenter, 
                                         hd_bin_mean, hd_bin_std, hd_bin_med, hd_bin_std,
                                         mb_df_mean, mb_df_std, mb_df_med, mb_df_mad]).reshape(1,11)
                    
                    if hd_compare_all_array is None:
                        hd_compare_all_array = bin_data
                    else:
                        hd_compare_all_array = np.concatenate((hd_compare_all_array, bin_data))
                    ds_name_list.append(hd_ds)
              
    hd_compare_all_cns = ['region', 'glacno', 'zbin', 
                          'hd_obs_mean', 'hd_obs_std', 'hd_obs_med', 'hd_obs_mad', 
                          'hd_ts_mean_m', 'hd_ts_std_m', 'hd_ts_med_m', 'hd_ts_mad_m']
    hd_compare_all = pd.DataFrame(hd_compare_all_array, columns=hd_compare_all_cns)
    hd_compare_all['hd_ds_name'] = ds_name_list
        
    #%% 
#    fig, ax = plt.subplots(1, 1, squeeze=False, figsize=(10,8), gridspec_kw = {'wspace':0, 'hspace':0})
#    glac_str_leg = []
#    for ndata, region in enumerate(hd_compare_all.region):
#        glac_str = str(int(region)) + '.' + str(int(hd_compare_all.loc[ndata,'glacno'])).zfill(5)
#        
#        # Legend
#        if glac_str not in glac_str_leg:
#            glac_str_leg.append(glac_str)
#            leg_label = glac_name_dict[glac_str]
#        else:
#            leg_label = ""
#            
#        ax[0,0].scatter(hd_compare_all.loc[ndata,cn_center_obs], hd_compare_all.loc[ndata,cn_center_ts], 
#                        color='k', 
#                        marker=glac_symbol[glac_str][0],
#                        facecolor=glac_symbol[glac_str][2], 
#                        linewidth=1,
#                        s=glac_symbol[glac_str][1],
#                        label=leg_label,
#                        zorder=3)
#        ax[0,0].errorbar(hd_compare_all.loc[ndata,cn_center_obs], hd_compare_all.loc[ndata,cn_center_ts], 
#                         xerr=hd_compare_all.loc[ndata,cn_spread_obs], 
#                         yerr=hd_compare_all.loc[ndata,cn_spread_ts], 
#                         capsize=1, linewidth=0.5, 
#                         color='darkgrey', 
#                         zorder=2)    
#    ax[0,0].set_xlabel('Observed $h_d$ (m)', size=12)    
#    ax[0,0].set_ylabel('Modeled $h_d$ (m)', size=12)
#    ymin = 0
#    ymax = 2.
#    xmin = 0
#    xmax = 2.
#    ax[0,0].set_xlim(xmin,xmax)
#    ax[0,0].set_ylim(ymin,ymax)
#    ax[0,0].plot([np.min([xmin,ymin]),np.max([xmax,ymax])], [np.min([xmin,ymin]),np.max([xmax,ymax])], color='k', 
#                 linewidth=0.5, zorder=1)
#    # Ensure proper order for legend
#    handles, labels = ax[0,0].get_legend_handles_labels()
#    labels, handles = zip(*sorted(zip(labels, handles), key=lambda t:t[0]))
#    ax[0,0].legend(handles, labels, loc=(0.62,0.05), ncol=1, fontsize=10, frameon=True, handlelength=1, 
#                   handletextpad=0.15, columnspacing=0.5, borderpad=0, labelspacing=0)
#    fig.set_size_inches(3.45,3.45)
#    fig.savefig(hd_compare_fp + 'hd_bin_validation.png', bbox_inches='tight', dpi=300)
        
    
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
        ymin, ymax = 0, 0.75
        xmin, xmax = 0, 0.75
        ax[0,0].set_xlim(xmin,xmax)
        ax[0,0].set_ylim(ymin,ymax)
        ax[0,0].plot([np.min([xmin,ymin]),np.max([xmax,ymax])], [np.min([xmin,ymin]),np.max([xmax,ymax])], color='k', 
                     linewidth=0.5, zorder=1)
        fig.text(0.5, 0.9, glac_name_dict[glac_str] + ' (' + glac_str + ')', va='bottom', ha='center', size=12)
        fig.set_size_inches(3.45,3.45)
        fig_fn = glac_str + '-hd_bin' + str(bin_width) + 'm_And2019.png'
        fig.savefig(hd_compare_fp + fig_fn, bbox_inches='tight', dpi=300)
        
        
    # ------ Suldenferner -------
    glac_str = '11.01604'
    reg = int(glac_str.split('.')[0])
    glacno = int(glac_str.split('.')[1])
    glac_idxs = np.where((hd_compare_all['region'].values == reg) & (hd_compare_all['glacno'] == glacno))[0]
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
        ymin, ymax = 0, 0.75
        xmin, xmax = 0, 0.75
        ax[0,0].set_xlim(xmin,xmax)
        ax[0,0].set_ylim(ymin,ymax)
        ax[0,0].plot([np.min([xmin,ymin]),np.max([xmax,ymax])], [np.min([xmin,ymin]),np.max([xmax,ymax])], color='k', 
                     linewidth=0.5, zorder=1)
        fig.text(0.5, 0.9, glac_name_dict[glac_str] + ' (' + glac_str + ')', va='bottom', ha='center', size=12)
        fig.set_size_inches(3.45,3.45)
        fig_fn = glac_str + '-hd_bin' + str(bin_width) + 'm_Nicholson.png'
        fig.savefig(hd_compare_fp + fig_fn, bbox_inches='tight', dpi=300)
        
        
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
        fig_fn = glac_str + '-hd_bin' + str(bin_width) + 'm.png'
        fig.savefig(hd_compare_fp + fig_fn, bbox_inches='tight', dpi=300)
    
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
        fig_fn = glac_str + '-hd_bin' + str(bin_width) + 'm.png'
        fig.savefig(hd_compare_fp + fig_fn, bbox_inches='tight', dpi=300)
    
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
        fig_fn = glac_str + '-hd_bin' + str(bin_width) + 'm_RM2014.png'
        fig.savefig(hd_compare_fp + fig_fn, bbox_inches='tight', dpi=300)
    
        
    # ----- Baltoro ------
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
        fig_fn = glac_str + '-hd_bin' + str(bin_width) + 'm_Mih2006.png'
        fig.savefig(hd_compare_fp + fig_fn, bbox_inches='tight', dpi=300)
    
#%%
if option_hd_spatial_compare:
    compare_suldenferner = False
    compare_ngozumpa = True
    
    # SULDENFERNER
    if compare_suldenferner:
        hd_ds = 'nicholson-suld'
    
        hd_obs = pd.read_csv(hd_ds_fn_dict[hd_ds])
        
        hd_obs[hd_obs['hd_ts_cal'] > 3] = np.nan
        hd_obs = hd_obs.dropna()
        hd_obs.reset_index(inplace=True, drop=True)
        
        glac_str = '11.01604'
        
        # Correlation
#        slope, intercept, r_value, p_value, std_err = linregress(hd_obs['hd_m'].values, hd_obs['hd_ts_cal'].values)
        
        # Plot
        fig, ax = plt.subplots(1, 1, squeeze=False, figsize=(10,8), gridspec_kw = {'wspace':0, 'hspace':0})

        ax[0,0].scatter(hd_obs['hd_m'], hd_obs['hd_ts_cal'], 
                        color='k', marker='o', facecolor='none', s=30, zorder=3)
 
        ax[0,0].set_xlabel('Observed $h_d$ (m)', size=12)    
        ax[0,0].set_ylabel('Modeled $h_d$ (m)', size=12)
        ymin, ymax = 0, 0.75
        xmin, xmax = 0, 0.75
        ax[0,0].set_xlim(xmin,xmax)
        ax[0,0].set_ylim(ymin,ymax)
        ax[0,0].plot([np.min([xmin,ymin]),np.max([xmax,ymax])], [np.min([xmin,ymin]),np.max([xmax,ymax])], color='k', 
                     linewidth=0.5, zorder=1)
        fig.text(0.5, 0.9, glac_name_dict[glac_str] + ' (' + glac_str + ')', va='bottom', ha='center', size=12)
        fig.set_size_inches(3.45,3.45)
        fig_fn = glac_str + '-hd_pts_Nicholson.png'
        fig.savefig(hd_compare_fp + fig_fn, bbox_inches='tight', dpi=300)
        
    
    # SULDENFERNER
    if compare_ngozumpa:
        hd_ds = 'mccarthy-ngoz'
    
        hd_obs = pd.read_csv(hd_ds_fn_dict[hd_ds])
        
        hd_obs[hd_obs['hd_ts_cal'] > 3] = np.nan
        hd_obs = hd_obs.dropna()
        hd_obs.reset_index(inplace=True, drop=True)
        
        glac_str = '15.03473'
        
        # Correlation
#        slope, intercept, r_value, p_value, std_err = linregress(hd_obs['hd_m'].values, hd_obs['hd_ts_cal'].values)
        
        # Plot
        fig, ax = plt.subplots(1, 1, squeeze=False, figsize=(10,8), gridspec_kw = {'wspace':0, 'hspace':0})

        ax[0,0].scatter(hd_obs['hd_m'], hd_obs['hd_ts_cal'], 
#                        color='k', marker='o', facecolor='none', s=30, zorder=3
                        color='k', marker='o', linewidth=0.1, facecolor='none', s=3, zorder=3)
 
        ax[0,0].set_xlabel('Observed $h_d$ (m)', size=12)    
        ax[0,0].set_ylabel('Modeled $h_d$ (m)', size=12)
        ymin, ymax = 0, 7.5
        xmin, xmax = 0, 7.5
        ax[0,0].set_xlim(xmin,xmax)
        ax[0,0].set_ylim(ymin,ymax)
        ax[0,0].plot([np.min([xmin,ymin]),np.max([xmax,ymax])], [np.min([xmin,ymin]),np.max([xmax,ymax])], color='k', 
                     linewidth=0.5, zorder=1)
        fig.text(0.5, 0.9, glac_name_dict[glac_str] + ' (' + glac_str + ')', va='bottom', ha='center', size=12)
        fig.set_size_inches(3.45,3.45)
        fig_fn = glac_str + '-hd_pts_McCarthyGPR.png'
        fig.savefig(hd_compare_fp + fig_fn, bbox_inches='tight', dpi=300)
        
        
        

    
        