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
option_hd_comparison = True
option_hd_spatial_compare = False

hd_obs_fp = debris_prms.main_directory + '/../hd_obs/'

melt_compare_fp = debris_prms.main_directory + '/../hd_obs/figures/hd_melt_compare/'
hd_compare_fp = debris_prms.main_directory + '/../hd_obs/figures/hd_obs_compare/'

if os.path.exists(melt_compare_fp) == False:
    os.makedirs(melt_compare_fp)
if os.path.exists(hd_compare_fp) == False:
    os.makedirs(hd_compare_fp)

#%% ===== FUNCTIONS =====
def plot_hd_vs_melt_comparison(measured_hd_list, measured_melt_list, yearfracs_list, glac_name, fig_fn, 
                               melt_fp, melt_fn, ds_names=None, 
                               hd_min=0, hd_max=2, hd_tick_major=0.25, hd_tick_minor=0.05, 
                               melt_min=0, melt_max=70, melt_tick_major=10, melt_tick_minor=5,
                               plot_meltfactor=False, z_value = 1.645, fontsize=11):
    #%%
    """ Plot comparison of debris vs. melt for various sites """
    # Dataset of melt data
    ds_ostrem = xr.open_dataset(melt_fp + melt_fn)
    ds_ostrem = ds_ostrem.sortby('hd_cm')
    
    time_year = pd.to_datetime(ds_ostrem.time.values).year
    time_daysperyear = np.array([366 if x%4 == 0 else 365 for x in time_year])
    time_yearfrac = time_year + (pd.to_datetime(ds_ostrem.time.values).dayofyear-1) / time_daysperyear

    color_dict = {0:'k', 1:'b', 2:'r'}
    symbol_dict = {0:'D', 1:'o', 2:'^'}
                    
    # ===== PLOT DEBRIS VS. SURFACE LOWERING ===== 
    fig, ax = plt.subplots(1, 1, squeeze=False, sharex=False, sharey=False, 
                          gridspec_kw = {'wspace':0.4, 'hspace':0.15})
    
#    #%%
#    # Seeing if there is melt over the winter - with the mean there is
#    ds_ostrem = xr.open_dataset('/Users/davidrounce/Documents/Dave_Rounce/DebrisGlaciers_WG/Melt_Intercomparison/' + 
#                                'output/ostrem_curves/2825N-8550E-debris_melt_curve.nc')
##    t1_idx = 5011   # Fall
##    t1_idx = 5080   # Winter
##    t1_idx = 5209   # spring
##    t2_idx = t1_idx + 12
#    t1_idx = 4000
#    t2_idx = 6000
#    fig, ax = plt.subplots(1, 1, squeeze=False, sharex=False, sharey=False, 
#                          gridspec_kw = {'wspace':0.4, 'hspace':0.15})
#    ax[0,0].plot(ds_ostrem.time.values[t1_idx:t2_idx], ds_ostrem.melt[0,t1_idx:t2_idx,0].values*1000)
#    
#    #%%
##    ds = xr.open_dataset('/Users/davidrounce/Documents/Dave_Rounce/DebrisGlaciers_WG/Melt_Intercomparison/output/' + 
##                         'exp4/HMA/Rounce2015_HMA-2825N-8550E-100MC_20200313.nc')
##    t1_idx = 125040 - 24 * 120
##    t2_idx = 125040+24*30 - 24 * 120
##    print(ds.melt[:,125040:125040+24*12,0].values.sum() * 1000 * 0.9 / 12)
##    print(ds.melt_med[:,125040:125040+24*12,0].values.sum() * 1000 * 0.9 / 12)
#    
#    ds = xr.open_dataset('/Users/davidrounce/Documents/Dave_Rounce/DebrisGlaciers_WG/Melt_Intercomparison/output/' + 
#                         'ostrem_curves/6150N-21700E-debris_melt_curve.nc')
#    t1_idx = 4000
#    t2_idx = 6000
#    
#    fig, ax = plt.subplots(1, 1, squeeze=False, sharex=False, sharey=False, 
#                          gridspec_kw = {'wspace':0.4, 'hspace':0.15})
#    ax[0,0].plot(ds.time.values[t1_idx:t2_idx], ds.melt[0,t1_idx:t2_idx,0].values*1000)
#    #%%

    for n in np.arange(0,len(measured_hd_list)):
        measured_hd = measured_hd_list[n]
        measured_melt = measured_melt_list[n]
        yearfracs = yearfracs_list[n]
        start_yearfrac = yearfracs[0]
        end_yearfrac = yearfracs[1]
        if ds_names is not None:
            ds_name = ds_names[n]
        else:
            ds_name = None
        
        start_idx = np.where(abs(time_yearfrac - start_yearfrac) == abs(time_yearfrac - start_yearfrac).min())[0][0]
        end_idx = np.where(abs(time_yearfrac - end_yearfrac) == abs(time_yearfrac - end_yearfrac).min())[0][0]
    
        # Ostrem Curve
        debris_thicknesses = ds_ostrem.hd_cm.values
        debris_melt_df = pd.DataFrame(np.zeros((len(debris_thicknesses),3)), 
                                      columns=['debris_thickness', 'melt_mmwed', 'melt_std_mmwed'])  
    
        nelev = 0
        for ndebris, debris_thickness in enumerate(debris_thicknesses):
            # Units: mm w.e. per day                
            melt_mmwed = (ds_ostrem['melt'][ndebris,start_idx:end_idx,nelev].values.sum() 
                          * 1000 / len(time_yearfrac[start_idx:end_idx]))
            melt_std_mmwed = (ds_ostrem['melt_std'][ndebris,start_idx:end_idx,nelev].values.sum() 
                              * 1000 / len(time_yearfrac[start_idx:end_idx]))
            debris_melt_df.loc[ndebris] = debris_thickness / 100, melt_mmwed, melt_std_mmwed
        debris_melt_df['melt_bndlow_mmwed'] = debris_melt_df['melt_mmwed'] - z_value * debris_melt_df['melt_std_mmwed']
        debris_melt_df['melt_bndhigh_mmwed'] = debris_melt_df['melt_mmwed'] + z_value * debris_melt_df['melt_std_mmwed']
            
        #%%
        # MEAN CURVE
        fit_idx = list(np.where(debris_thicknesses >= 5)[0])            
        func_coeff, pcov = curve_fit(melt_fromdebris_func, 
                                     debris_melt_df.debris_thickness.values[fit_idx], 
                                     debris_melt_df.melt_mmwed.values[fit_idx])
        melt_cleanice = debris_melt_df.loc[0,'melt_mmwed']
        # Fitted curve
        debris_4curve = np.arange(0.02,5.01,0.01)
        melt_4curve = melt_fromdebris_func(debris_4curve, func_coeff[0], func_coeff[1])
        # add clean ice
        debris_4curve = np.concatenate([[0.0], debris_4curve])
        melt_4curve = np.concatenate([[melt_cleanice], melt_4curve])
        
        # LOWER BOUND CURVE
        func_coeff_bndlow, pcov = curve_fit(melt_fromdebris_func, 
                                            debris_melt_df.debris_thickness.values[fit_idx], 
                                            debris_melt_df.melt_bndlow_mmwed.values[fit_idx])
        melt_cleanice_bndlow = debris_melt_df.loc[0,'melt_bndlow_mmwed']
        # Fitted curve
        debris_4curve = np.arange(0.02,5.01,0.01)
        melt_4curve_bndlow = melt_fromdebris_func(debris_4curve, func_coeff_bndlow[0], func_coeff_bndlow[1])
        # add clean ice
        debris_4curve = np.concatenate([[0.0], debris_4curve])
        melt_4curve_bndlow = np.concatenate([[melt_cleanice_bndlow], melt_4curve_bndlow])
        
        # UPPER BOUND CURVE
        func_coeff_bndhigh, pcov = curve_fit(melt_fromdebris_func, 
                                            debris_melt_df.debris_thickness.values[fit_idx], 
                                            debris_melt_df.melt_bndhigh_mmwed.values[fit_idx])
        melt_cleanice_bndhigh = debris_melt_df.loc[0,'melt_bndhigh_mmwed']
        # Fitted curve
        debris_4curve = np.arange(0.02,5.01,0.01)
        melt_4curve_bndhigh = melt_fromdebris_func(debris_4curve, func_coeff_bndhigh[0], func_coeff_bndhigh[1])
        # add clean ice
        debris_4curve = np.concatenate([[0.0], debris_4curve])
        melt_4curve_bndhigh = np.concatenate([[melt_cleanice_bndhigh], melt_4curve_bndhigh])
        
        if plot_meltfactor:
            melt_4curve = melt_4curve / melt_cleanice
            melt_4curve_bndlow = melt_4curve_bndlow / melt_cleanice
            melt_4curve_bndhigh = melt_4curve_bndhigh / melt_cleanice
        
        # Plot curve
        ax[0,0].plot(measured_hd, measured_melt, symbol_dict[n], color=color_dict[n], 
                     markersize=3, markerfacecolor="None", markeredgewidth=0.5, zorder=1, label=ds_name, clip_on=False)
        ax[0,0].plot(debris_4curve, melt_4curve, 
                     color=color_dict[n], linewidth=1, linestyle='--', zorder=2)
        ax[0,0].fill_between(debris_4curve, melt_4curve_bndlow, melt_4curve_bndhigh,
                             color=color_dict[n], linewidth=0, zorder=2, alpha=0.2)
#        ax[0,0].plot(debris_4curve, melt_4curve_bndhigh, 
#                     color=color_dict[n], linewidth=1, linestyle=':', zorder=2)
        
    # text
    ax[0,0].text(0.5, 1.09, glac_name, size=fontsize-2, horizontalalignment='center', verticalalignment='top', 
                 transform=ax[0,0].transAxes)
#    eqn_text = r'$b = \frac{b_{0}}{1 + kb_{0}h}$'
#    coeff1_text = r'$b_{0} = ' + str(np.round(func_coeff[0],2)) + '$' 
#    coeff2_text = r'$k = ' + str(np.round(func_coeff[1],2)) + '$' 
#    # coeff$\frac{b_{0}}{1 + 2kb_{0}h}$'
#    ax[0,0].text(0.9, 0.95, eqn_text, size=12, horizontalalignment='right', verticalalignment='top', 
#                 transform=ax[0,0].transAxes)
#    ax[0,0].text(0.615, 0.83, 'where', size=10, horizontalalignment='left', verticalalignment='top', 
#                 transform=ax[0,0].transAxes)
#    ax[0,0].text(0.66, 0.77, coeff1_text, size=10, horizontalalignment='left', verticalalignment='top', 
#                 transform=ax[0,0].transAxes)
#    ax[0,0].text(0.66, 0.7, coeff2_text, size=10, horizontalalignment='left', verticalalignment='top', 
#                 transform=ax[0,0].transAxes)
    # X-label
    ax[0,0].set_xlabel('Debris thickness (m)', size=fontsize)
    ax[0,0].set_xlim(hd_min, hd_max)
    ax[0,0].xaxis.set_major_locator(plt.MultipleLocator(hd_tick_major))
    ax[0,0].xaxis.set_minor_locator(plt.MultipleLocator(hd_tick_minor))  
    # Y-label
    if plot_meltfactor:
        ylabel_str = 'Melt (-)'
    else:
        ylabel_str = 'Melt (mm w.e. d$^{-1}$)'
    ax[0,0].set_ylabel(ylabel_str, size=fontsize)
    ax[0,0].set_ylim(melt_min, melt_max)
    ax[0,0].yaxis.set_major_locator(plt.MultipleLocator(melt_tick_major))
    ax[0,0].yaxis.set_minor_locator(plt.MultipleLocator(melt_tick_minor))
    # Tick parameters
    ax[0,0].yaxis.set_ticks_position('both')
    ax[0,0].tick_params(axis='both', which='major', labelsize=fontsize-2, direction='inout')
    ax[0,0].tick_params(axis='both', which='minor', labelsize=fontsize-4, direction='in') 
    # Legend
    ax[0,0].legend(ncol=1, fontsize=fontsize-3, frameon=True, handlelength=1, 
                   handletextpad=0.15, columnspacing=0.5, borderpad=0.25, labelspacing=0.5, framealpha=0.5)
    # Save plot
    fig.set_size_inches(2, 2)
    fig.savefig(melt_compare_fp + fig_fn, bbox_inches='tight', dpi=300, transparent=True)
    plt.close()


#%%
if option_melt_comparison:
#    glaciers = ['1.15645', '2.14297', '6.00474', '7.01044', '10.01732', '11.00719', '11.02810', '11.02858', '11.03005', 
#                '12.01012', '12.01132', '13.05000', '13.43232', '14.06794', '14.16042', '15.03733', '15.04045', 
#                '15.07886', '15.11758', '18.02397']
    glaciers = ['1.15645', '2.14297', '6.00474', '7.01044', '11.00719', '11.02472', '11.02810', '11.02858', '11.03005', 
                '12.01012', '12.01132', '13.05000', '13.43232', '14.06794', '14.16042', '15.03733', '15.04045', 
                '15.07886', '15.11758', '18.02397']
#    glaciers = ['10.01732']
    glaciers = ['13.43165']
    
    print('add comparison to Haidong etal 2006')
    
    z_value = 1.645
    
    # ===== KENNICOTT (1.15645) ====
    if '1.15645' in glaciers:
        print('\nmelt comparison with Anderson et al. 2019')
        # Data: debris thickness (m) and melt rate (mm w.e. d-1)
        mb_df = pd.read_csv(debris_prms.main_directory + '/../hd_obs/datasets/1.15645_kennicott_anderson_2019-melt.csv')
        measured_hd_list = [mb_df.hd_m.values]
        measured_melt_list = [mb_df.melt_mmwed.values]
        glac_name = "Kennicott (1.15645)"
        fig_fn = '1.15645_hd_melt_And2019.png'
#        ds_names = ['Anderson 2019\n(6/18/11$\u2009$-$\u2009$8/16/11)']
        ds_names = ['6/18/11$\u2009$-$\u2009$8/16/11']
        melt_fp = debris_prms.output_fp + 'ostrem_curves/exp4/'
        melt_fn = '6150N-21700E-debris_melt_curve.nc'
        yearfracs_list = [[2011 + 169/365, 2011 + 228/365]]

        hd_min, hd_max = 0, np.ceil(np.max([x.max() for x in measured_hd_list])/0.1)*0.1 + 0.05
        hd_tick_major, hd_tick_minor = 0.1, 0.02
        melt_min, melt_max = 0, np.round(np.ceil(np.max([x.max() for x in measured_melt_list]) / 10) * 10,0) + 5
        melt_tick_major, melt_tick_minor = 10, 5
    
        for n in np.arange(0,len(measured_hd_list)):
            assert len(measured_hd_list[n]) == len(measured_melt_list[n]), 'length of hd differs from melt'
        
        plot_hd_vs_melt_comparison(measured_hd_list, measured_melt_list, yearfracs_list, glac_name, fig_fn, 
                                   melt_fp, melt_fn,
                                   ds_names=ds_names, hd_min=hd_min, hd_max=hd_max, 
                                   hd_tick_major=hd_tick_major, hd_tick_minor=hd_tick_minor,
                                   melt_min=melt_min, melt_max=melt_max, 
                                   melt_tick_major=melt_tick_major, melt_tick_minor=melt_tick_minor)
        
    # ===== Emmons (2.14297) ====
    if '2.14297' in glaciers:
        print('\nmelt comparison with Moore et al. 2019')
        # Data: debris thickness (m) and melt rate (mm w.e. d-1)
        mb_df = pd.read_csv(debris_prms.main_directory + '/../hd_obs/datasets/2.14297_moore2019-melt.csv')
        measured_hd_list = [mb_df.hd_m.values]
        measured_melt_list = [mb_df.melt_mmwed.values]
        glac_name = "Emmons (2.14297)"
        fig_fn = '2.14297_hd_melt_Moo2019.png'
#        ds_names = ['Moore 2019\n(7/31/14$\u2009$-$\u2009$8/10/14)']
        ds_names = ['7/31/14$\u2009$-$\u2009$8/10/14']
        melt_fp = debris_prms.output_fp + 'ostrem_curves/exp4/'
        melt_fn = '4700N-23825E-debris_melt_curve.nc'
        yearfracs_list = [[2014 + 212/365, 2014 + 222/365]]

        hd_min, hd_max = 0, np.ceil(np.max([x.max() for x in measured_hd_list])/0.1)*0.1 + 0.05
        hd_tick_major, hd_tick_minor = 0.1, 0.02
        melt_min, melt_max = 0, np.round(np.ceil(np.max([x.max() for x in measured_melt_list]) / 10) * 10,0) + 5
        melt_tick_major, melt_tick_minor = 10, 5
    
        for n in np.arange(0,len(measured_hd_list)):
            assert len(measured_hd_list[n]) == len(measured_melt_list[n]), 'length of hd differs from melt'
        
        plot_hd_vs_melt_comparison(measured_hd_list, measured_melt_list, yearfracs_list, glac_name, fig_fn, 
                                   melt_fp, melt_fn,
                                   ds_names=ds_names, hd_min=hd_min, hd_max=hd_max, 
                                   hd_tick_major=hd_tick_major, hd_tick_minor=hd_tick_minor,
                                   melt_min=melt_min, melt_max=melt_max, 
                                   melt_tick_major=melt_tick_major, melt_tick_minor=melt_tick_minor)
    
    # ===== Svinafellsjokull (06.00474) ====
    if '6.00474' in glaciers:
        print('\nmelt comparison with Moller et al (2016)')
        # Data: debris thickness (m) and melt rate (mm w.e. d-1)
        mb_df = pd.read_csv(debris_prms.main_directory + '/../hd_obs/datasets/6.00474_moller2016-melt.csv')
        measured_hd_list = [mb_df.hd_m.values]
        measured_melt_list = [mb_df['melt_mf'].values]
        glac_name = "Svinafellsjokull (6.00474)"
        fig_fn = '6.00474_hd_melt_Moller2016.png'
#        ds_names = ['Moller 2016\n(5/17/13$\u2009$-$\u2009$5/30/13)']
        ds_names = ['5/17/13$\u2009$-$\u2009$5/30/13']
        melt_fp = debris_prms.output_fp + 'ostrem_curves/exp4/'
        melt_fn = '6400N-34325E-debris_melt_curve.nc'
        yearfracs_list = [[2013 + 137/365, 2013 + 150/365]]
        
        hd_min, hd_max = 0, np.ceil(np.max([x.max() for x in measured_hd_list])/0.1)*0.1 + 0.05
        hd_tick_major, hd_tick_minor = 0.1, 0.02
#        melt_min, melt_max = 0, np.round(np.ceil(np.max([x.max() for x in measured_melt_list]) / 10) * 10,0) + 5
#        melt_tick_major, melt_tick_minor = 10, 5
        melt_min, melt_max = 0, np.round(np.ceil(np.max([x.max() for x in measured_melt_list]) / 0.1) * 0.1,1) + 0.1
        melt_tick_major, melt_tick_minor = 0.5, 0.1
    
        for n in np.arange(0,len(measured_hd_list)):
            assert len(measured_hd_list[n]) == len(measured_melt_list[n]), 'length of hd differs from melt'
        
        plot_hd_vs_melt_comparison(measured_hd_list, measured_melt_list, yearfracs_list, glac_name, fig_fn, 
                                   melt_fp, melt_fn,
                                   ds_names=ds_names, hd_min=hd_min, hd_max=hd_max, 
                                   hd_tick_major=hd_tick_major, hd_tick_minor=hd_tick_minor,
                                   melt_min=melt_min, melt_max=melt_max, 
                                   melt_tick_major=melt_tick_major, melt_tick_minor=melt_tick_minor, 
                                   plot_meltfactor=True)
        
    # ===== Larsbreen (7.01044) ====
    if '7.01044' in glaciers:
        print('\nmelt comparison with Nicholson and Benn 2006')
        # Data: debris thickness (m) and melt rate (mm w.e. d-1)
        mb_df = pd.read_csv(debris_prms.main_directory + '/../hd_obs/datasets/7.01044_larsbreen_NB2006-melt.csv')
        measured_hd_list = [mb_df.hd_m.values]
        measured_melt_list = [mb_df.melt_mmwed.values]
        glac_name = "Larsbreen (7.01044)"
        fig_fn = '7.01044_hd_melt_NichBenn2006.png'
#        ds_names = ['Nicholson 2006\n(7/09/02$\u2009$-$\u2009$7/20/02)']
        ds_names = ['7/09/02$\u2009$-$\u2009$7/20/02']
        melt_fp = debris_prms.output_fp + 'ostrem_curves/exp4/'
        melt_fn = '7825N-1600E-debris_melt_curve.nc'
        yearfracs_list = [[2002 + 191/366, 2002 + 202/366]]

        hd_min, hd_max = 0, np.ceil(np.max([x.max() for x in measured_hd_list])/0.1)*0.1 + 0.05
        hd_tick_major, hd_tick_minor = 0.1, 0.02
        melt_min, melt_max = 0, np.round(np.ceil(np.max([x.max() for x in measured_melt_list]) / 10) * 10,0) + 5
        melt_tick_major, melt_tick_minor = 10, 5
    
        for n in np.arange(0,len(measured_hd_list)):
            assert len(measured_hd_list[n]) == len(measured_melt_list[n]), 'length of hd differs from melt'
        
        plot_hd_vs_melt_comparison(measured_hd_list, measured_melt_list, yearfracs_list, glac_name, fig_fn, 
                                   melt_fp, melt_fn,
                                   ds_names=ds_names, hd_min=hd_min, hd_max=hd_max, 
                                   hd_tick_major=hd_tick_major, hd_tick_minor=hd_tick_minor,
                                   melt_min=melt_min, melt_max=melt_max, 
                                   melt_tick_major=melt_tick_major, melt_tick_minor=melt_tick_minor)
        
    # ===== Maliy Aktru (10.01732) ====
    if '10.01732' in glaciers:
#        print('\nmelt comparison with Mayer et al (2011)')
        assert True == False, '10.01732 NEEDS TO DO THE MODELING FIRST!'
        # Data: debris thickness (m) and melt rate (mm w.e. d-1)
        mb_df = pd.read_csv(debris_prms.main_directory + '/../hd_obs/datasets/10.01732_mayer2011-melt.csv')
        measured_hd_list = [mb_df.hd_m.values]
        measured_melt_list = [mb_df['melt_mf'].values]
        glac_name = "Maliy Aktru (10.01732)"
        fig_fn = '10.01732_hd_melt_Mayer2011.png'
#        ds_names = ['Mayer 2011\n(7/11/07$\u2009$-$\u2009$7/30/07)']
        ds_names = ['7/11/07$\u2009$-$\u2009$7/30/07']
        melt_fp = debris_prms.output_fp + 'ostrem_curves/exp4/'
        melt_fn = '5000N-8775E-debris_melt_curve.nc'
        yearfracs_list = [[2007 + 192/365, 2007 + 211/365]]

        hd_min, hd_max = 0, np.ceil(np.max([x.max() for x in measured_hd_list])/0.1)*0.1 + 0.05
        hd_tick_major, hd_tick_minor = 0.1, 0.02
        melt_min, melt_max = 0, np.round(np.ceil(np.max([x.max() for x in measured_melt_list]) / 10) * 10,0) + 5
        melt_tick_major, melt_tick_minor = 10, 5
    
        for n in np.arange(0,len(measured_hd_list)):
            assert len(measured_hd_list[n]) == len(measured_melt_list[n]), 'length of hd differs from melt'
        
        plot_hd_vs_melt_comparison(measured_hd_list, measured_melt_list, yearfracs_list, glac_name, fig_fn, 
                                   melt_fp, melt_fn,
                                   ds_names=ds_names, hd_min=hd_min, hd_max=hd_max, 
                                   hd_tick_major=hd_tick_major, hd_tick_minor=hd_tick_minor,
                                   melt_min=melt_min, melt_max=melt_max, 
                                   melt_tick_major=melt_tick_major, melt_tick_minor=melt_tick_minor,
                                   plot_meltfactor=True)
    
    # ===== Vernagtferner (11.00719) ====
    if '11.00719' in glaciers:
        print('\nmelt comparison with Juen et al (2013)')
        # Data: debris thickness (m) and melt rate (mm w.e. d-1)
        mb_df = pd.read_csv(debris_prms.main_directory + '/../hd_obs/datasets/11.00719_vernagtferner_juen2013-melt.csv')
        measured_hd_list = [mb_df.hd_m.values]
        measured_melt_list = [mb_df.melt_mmwed.values]
        glac_name = "Vernagtferner (11.00719)"
        fig_fn = '11.00719_hd_melt_Juen2013.png'
#        ds_names = ['Juen 2013\n(6/25/10$\u2009$-$\u2009$7/10/10)']
        ds_names = ['6/25/10$\u2009$-$\u2009$7/10/10']
        melt_fp = debris_prms.output_fp + 'ostrem_curves/exp4/'
        melt_fn = '4700N-1075E-debris_melt_curve.nc'
        yearfracs_list = [[2010 + 176/365, 2010 + 191/365]]

        hd_min, hd_max = 0, np.ceil(np.max([x.max() for x in measured_hd_list])/0.1)*0.1 + 0.05
        hd_tick_major, hd_tick_minor = 0.1, 0.02
        melt_min, melt_max = 0, np.round(np.ceil(np.max([x.max() for x in measured_melt_list]) / 10) * 10,0) + 5
        melt_tick_major, melt_tick_minor = 10, 5
    
        for n in np.arange(0,len(measured_hd_list)):
            assert len(measured_hd_list[n]) == len(measured_melt_list[n]), 'length of hd differs from melt'
        
        plot_hd_vs_melt_comparison(measured_hd_list, measured_melt_list, yearfracs_list, glac_name, fig_fn, 
                                   melt_fp, melt_fn,
                                   ds_names=ds_names, hd_min=hd_min, hd_max=hd_max, 
                                   hd_tick_major=hd_tick_major, hd_tick_minor=hd_tick_minor,
                                   melt_min=melt_min, melt_max=melt_max, 
                                   melt_tick_major=melt_tick_major, melt_tick_minor=melt_tick_minor)
        
    # ===== Vernocolo (11.02472) =====
    if '11.02472' in glaciers:
        print('\nmelt comparison with bocchiola et al. (2015)')
        # Data: debris thickness (m) and melt rate (mm w.e. d-1)
        mb_df = pd.read_csv(debris_prms.main_directory + '/../hd_obs/datasets/11.02472_bocchiola2015-melt.csv')
        measured_hd_list = [mb_df.hd_m.values]
        measured_melt_list = [mb_df.melt_mmwed.values]
        glac_name = "Venerocolo (11.02472)"
        fig_fn = '11.02472_hd_melt_Boc2015.png'
        ds_names = ['8/10/07$\u2009$-$\u2009$9/13/07']
        melt_fp = debris_prms.output_fp + 'ostrem_curves/exp4/'
        melt_fn = '4625N-1050E-debris_melt_curve.nc'
        yearfracs_list = [[2007 + 222/365, 2007 + 256/365]]

        hd_min, hd_max = 0, np.ceil(np.max([x.max() for x in measured_hd_list])/0.1)*0.1 + 0.05
        hd_tick_major, hd_tick_minor = 0.1, 0.02
        melt_min, melt_max = 0, np.round(np.ceil(np.max([x.max() for x in measured_melt_list]) / 10) * 10,0) + 5
        melt_tick_major, melt_tick_minor = 10, 5
    
        for n in np.arange(0,len(measured_hd_list)):
            assert len(measured_hd_list[n]) == len(measured_melt_list[n]), 'length of hd differs from melt'
        
        plot_hd_vs_melt_comparison(measured_hd_list, measured_melt_list, yearfracs_list, glac_name, fig_fn, 
                                   melt_fp, melt_fn,
                                   ds_names=ds_names, hd_min=hd_min, hd_max=hd_max, 
                                   hd_tick_major=hd_tick_major, hd_tick_minor=hd_tick_minor,
                                   melt_min=melt_min, melt_max=melt_max, 
                                   melt_tick_major=melt_tick_major, melt_tick_minor=melt_tick_minor)
        
    # ===== Arolla (11.02810) ====
    if '11.02810' in glaciers:
        print('\nmelt comparison with Reid et al (2012)')
        # Data: debris thickness (m) and melt rate (mm w.e. d-1)
        measured_hd_list = [np.array([2, 3, 4, 6]) / 100]
        measured_melt_list = [np.array([1.22, 0.98, 0.72, 1.16]) * 1000 / (252-209)]
        glac_name = "Arolla (11.02810)"
        fig_fn = '11.02810_hd_melt_Reid2012.png'
#        ds_names = ['Reid 2012\n(7/28/10$\u2009$-$\u2009$9/09/10)']
        ds_names = ['7/28/10$\u2009$-$\u2009$9/09/10']
        melt_fp = debris_prms.output_fp + 'ostrem_curves/exp4/'
        melt_fn = '4600N-750E-debris_melt_curve.nc'
        yearfracs_list = [[2010 + 209/365, 2010 + 252/365]]

        hd_min, hd_max = 0, np.ceil(np.max([x.max() for x in measured_hd_list])/0.1)*0.1 + 0.05
        hd_tick_major, hd_tick_minor = 0.1, 0.02
        melt_min, melt_max = 0, np.round(np.ceil(np.max([x.max() for x in measured_melt_list]) / 10) * 10,0) + 5
        melt_tick_major, melt_tick_minor = 10, 5
    
        for n in np.arange(0,len(measured_hd_list)):
            assert len(measured_hd_list[n]) == len(measured_melt_list[n]), 'length of hd differs from melt'
        
        plot_hd_vs_melt_comparison(measured_hd_list, measured_melt_list, yearfracs_list, glac_name, fig_fn, 
                                   melt_fp, melt_fn,
                                   ds_names=ds_names, hd_min=hd_min, hd_max=hd_max, 
                                   hd_tick_major=hd_tick_major, hd_tick_minor=hd_tick_minor,
                                   melt_min=melt_min, melt_max=melt_max, 
                                   melt_tick_major=melt_tick_major, melt_tick_minor=melt_tick_minor)
        
    # ===== Belvedere (11.02858) ====
    if '11.02858' in glaciers:
        print('\nmelt comparison with Nicholson and Benn (2006)')
        # Data: debris thickness (m) and melt rate (mm w.e. d-1)
        mb_df = pd.read_csv(debris_prms.main_directory + '/../hd_obs/datasets/11.02858_belvedere_nb2006-melt.csv')
        measured_hd_list = [mb_df.hd_m.values]
        measured_melt_list = [mb_df.melt_mmwed.values]
        glac_name = "Belvedere (11.02858)"
        fig_fn = '11.02858_hd_melt_NB2006.png'
#        ds_names = ['Nicholson 2006\n(8/06/03$\u2009$-$\u2009$8/10/03)']
        ds_names = ['8/06/03$\u2009$-$\u2009$8/10/03']
        melt_fp = debris_prms.output_fp + 'ostrem_curves/exp4/'
        melt_fn = '4600N-800E-debris_melt_curve.nc'
        yearfracs_list = [[2003 + 218/365, 2003 + 222/365]]

        hd_min, hd_max = 0, np.ceil(np.max([x.max() for x in measured_hd_list])/0.1)*0.1 + 0.05
        hd_tick_major, hd_tick_minor = 0.1, 0.02
        melt_min, melt_max = 0, np.round(np.ceil(np.max([x.max() for x in measured_melt_list]) / 10) * 10,0) + 5
        melt_tick_major, melt_tick_minor = 10, 5
    
        for n in np.arange(0,len(measured_hd_list)):
            assert len(measured_hd_list[n]) == len(measured_melt_list[n]), 'length of hd differs from melt'
        
        plot_hd_vs_melt_comparison(measured_hd_list, measured_melt_list, yearfracs_list, glac_name, fig_fn, 
                                   melt_fp, melt_fn,
                                   ds_names=ds_names, hd_min=hd_min, hd_max=hd_max, 
                                   hd_tick_major=hd_tick_major, hd_tick_minor=hd_tick_minor,
                                   melt_min=melt_min, melt_max=melt_max, 
                                   melt_tick_major=melt_tick_major, melt_tick_minor=melt_tick_minor)
        
    # ===== MIAGE (11.03005) ====
    if '11.03005' in glaciers:
        print('\nmelt comparison with Reid and Brock (2010)')
        # Data: debris thickness (m) and melt rate (mm w.e. d-1)
        measured_hd_list = [np.array([55, 26, 21, 19.5, 19.5, 19.5, 17, 16, 16, 16, 15, 12, 11, 10]) / 100]
        measured_melt_list = [np.array([7, 17, 16, 15, 19, 22, 17, 19, 20, 22, 20, 20, 22, 21])]
        glac_name = 'Miage (11.03005)'
        fig_fn = '11.03005_hd_melt_Reid2010.png'
#        ds_names = ['Reid 2010\n(6/21/05$\u2009$-$\u2009$9/04/05)']
        ds_names = ['6/21/05$\u2009$-$\u2009$9/04/05']
        melt_fp = debris_prms.output_fp + 'ostrem_curves/exp4/'
        melt_fn = '4650N-1050E-debris_melt_curve.nc'
        yearfracs_list = [[2005 + 172/365, 2005 + 247/365]]

        hd_min, hd_max = 0, np.ceil(np.max([x.max() for x in measured_hd_list])/0.1)*0.1 + 0.05
        hd_tick_major, hd_tick_minor = 0.1, 0.02
        melt_min, melt_max = 0, np.round(np.ceil(np.max([x.max() for x in measured_melt_list]) / 10) * 10,0) + 5
        melt_tick_major, melt_tick_minor = 10, 5
    
        for n in np.arange(0,len(measured_hd_list)):
            assert len(measured_hd_list[n]) == len(measured_melt_list[n]), 'length of hd differs from melt'
        
        plot_hd_vs_melt_comparison(measured_hd_list, measured_melt_list, yearfracs_list, glac_name, fig_fn, 
                                   melt_fp, melt_fn,
                                   ds_names=ds_names, hd_min=hd_min, hd_max=hd_max, 
                                   hd_tick_major=hd_tick_major, hd_tick_minor=hd_tick_minor,
                                   melt_min=melt_min, melt_max=melt_max, 
                                   melt_tick_major=melt_tick_major, melt_tick_minor=melt_tick_minor)
        
    # ===== Zopkhito (12.01012) ====
    if '12.01012' in glaciers:
        print('\nmelt comparison with Lambrecht et al (2011)')
        # Data: debris thickness (m) and melt rate (mm w.e. d-1)
        mb_df = pd.read_csv(debris_prms.main_directory + '/../hd_obs/datasets/12.01012_lambrecht2011-melt2008.csv')
        mb_df2 = pd.read_csv(debris_prms.main_directory + '/../hd_obs/datasets/12.01012_lambrecht2011-melt2009.csv')
        measured_hd_list = [mb_df.hd_m.values, mb_df2.hd_m.values]
        measured_melt_list = [mb_df['melt_mf'].values, mb_df2['melt_mf'].values]
        glac_name = "Zopkhito (12.01012)"
        fig_fn = '12.01012_hd_melt_Lambrecht2011.png'
#        ds_names = ['Lambrecht 2011\n(6/20/08$\u2009$-$\u2009$6/27/08)', 
#                    'Lambrecht 2011\n(7/01/09$\u2009$-$\u2009$7/08/09)']
        ds_names = ['6/20/08$\u2009$-$\u2009$6/27/08', 
                    '7/01/09$\u2009$-$\u2009$7/08/09']
        melt_fp = debris_prms.output_fp + 'ostrem_curves/exp4/'
        melt_fn = '4300N-4350E-debris_melt_curve.nc'
        yearfracs_list = [[2008 + 172/366, 2008 + 179/366], [2009 + 182/365, 2009 + 189/365]]

        hd_min, hd_max = 0, np.ceil(np.max([x.max() for x in measured_hd_list])/0.1)*0.1 + 0.05
        hd_tick_major, hd_tick_minor = 0.1, 0.02
#        melt_min, melt_max = 0, np.round(np.ceil(np.max([x.max() for x in measured_melt_list]) / 10) * 10,0) + 5
#        melt_tick_major, melt_tick_minor = 10, 5
        melt_min, melt_max = 0, np.round(np.ceil(np.max([x.max() for x in measured_melt_list]) / 0.1) * 0.1,1) + 0.1
        melt_tick_major, melt_tick_minor = 0.5, 0.1
    
        for n in np.arange(0,len(measured_hd_list)):
            assert len(measured_hd_list[n]) == len(measured_melt_list[n]), 'length of hd differs from melt'
        
        plot_hd_vs_melt_comparison(measured_hd_list, measured_melt_list, yearfracs_list, glac_name, fig_fn, 
                                   melt_fp, melt_fn,
                                   ds_names=ds_names, hd_min=hd_min, hd_max=hd_max, 
                                   hd_tick_major=hd_tick_major, hd_tick_minor=hd_tick_minor,
                                   melt_min=melt_min, melt_max=melt_max, 
                                   melt_tick_major=melt_tick_major, melt_tick_minor=melt_tick_minor,
                                   plot_meltfactor=True)
        
    # ===== Djankuat (12.01132) ====
    if '12.01132' in glaciers:
        print('\nmelt comparison with Lambrecht et al (2011)')
        # Data: debris thickness (m) and melt rate (mm w.e. d-1)
        mb_df = pd.read_csv(debris_prms.main_directory + '/../hd_obs/datasets/12.01132_lambrecht2011-melt.csv')
        measured_hd_list = [mb_df.hd_m.values]
        measured_melt_list = [mb_df['melt_mf'].values]
        glac_name = "Djankuat (12.01132)"
        fig_fn = '12.01132_hd_melt_Lambrecht2011.png'
#        ds_names = ['Lambrecht 2011\n(6/20/08$\u2009$-$\u2009$9/02/08)']
        ds_names = ['6/20/08$\u2009$-$\u2009$9/02/08']
        melt_fp = debris_prms.output_fp + 'ostrem_curves/exp4/'
        melt_fn = '4325N-4275E-debris_melt_curve.nc'
        yearfracs_list = [[2008 + 172/366, 2008 + 246/366]]

        hd_min, hd_max = 0, np.ceil(np.max([x.max() for x in measured_hd_list])/0.1)*0.1 + 0.05
        hd_tick_major, hd_tick_minor = 0.1, 0.02
#        melt_min, melt_max = 0, np.round(np.ceil(np.max([x.max() for x in measured_melt_list]) / 10) * 10,0) + 5
#        melt_tick_major, melt_tick_minor = 10, 5
        melt_min, melt_max = 0, np.round(np.ceil(np.max([x.max() for x in measured_melt_list]) / 0.1) * 0.1,1) + 0.1
        melt_tick_major, melt_tick_minor = 0.5, 0.1
    
        for n in np.arange(0,len(measured_hd_list)):
            assert len(measured_hd_list[n]) == len(measured_melt_list[n]), 'length of hd differs from melt'
        
        plot_hd_vs_melt_comparison(measured_hd_list, measured_melt_list, yearfracs_list, glac_name, fig_fn, 
                                   melt_fp, melt_fn,
                                   ds_names=ds_names, hd_min=hd_min, hd_max=hd_max, 
                                   hd_tick_major=hd_tick_major, hd_tick_minor=hd_tick_minor,
                                   melt_min=melt_min, melt_max=melt_max, 
                                   melt_tick_major=melt_tick_major, melt_tick_minor=melt_tick_minor,
                                   plot_meltfactor=True)
        
    # ===== S Inylchek (13.05000) ====
    if '13.05000' in glaciers:
        print('\nmelt comparison with Hagg et al (2008)')
        # Data: debris thickness (m) and melt rate (mm w.e. d-1)
        mb_df = pd.read_csv(debris_prms.main_directory + '/../hd_obs/datasets/13.05000_hagg2008-melt.csv')
        measured_hd_list = [mb_df.hd_m.values]
        measured_melt_list = [mb_df.melt_mmwed.values]
        glac_name = "S Inylchek (13.05000)"
        fig_fn = '13.05000_hd_melt_Hagg2008.png'
#        ds_names = ['Hagg 2008\n(7/30/05$\u2009$-$\u2009$8/10/05)']
        ds_names = ['7/30/05$\u2009$-$\u2009$8/10/05']
        melt_fp = debris_prms.output_fp + 'ostrem_curves/exp4/'
        melt_fn = '4200N-8025E-debris_melt_curve.nc'
        yearfracs_list = [[2005 + 211/365, 2005 + 222/365]]

        hd_min, hd_max = 0, np.ceil(np.max([x.max() for x in measured_hd_list])/0.1)*0.1 + 0.05
        hd_tick_major, hd_tick_minor = 0.1, 0.02
        melt_min, melt_max = 0, np.round(np.ceil(np.max([x.max() for x in measured_melt_list]) / 10) * 10,0) + 5
        melt_tick_major, melt_tick_minor = 10, 5
    
        for n in np.arange(0,len(measured_hd_list)):
            assert len(measured_hd_list[n]) == len(measured_melt_list[n]), 'length of hd differs from melt'
        
        plot_hd_vs_melt_comparison(measured_hd_list, measured_melt_list, yearfracs_list, glac_name, fig_fn, 
                                   melt_fp, melt_fn,
                                   ds_names=ds_names, hd_min=hd_min, hd_max=hd_max, 
                                   hd_tick_major=hd_tick_major, hd_tick_minor=hd_tick_minor,
                                   melt_min=melt_min, melt_max=melt_max, 
                                   melt_tick_major=melt_tick_major, melt_tick_minor=melt_tick_minor)
        
    # ===== Koxkar (13.43232) ====
    if '13.43232' in glaciers:
        print('\nmelt comparison with Juen et al (2014)')
        # Data: debris thickness (m) and melt rate (mm w.e. d-1)
        mb_df = pd.read_csv(debris_prms.main_directory + '/../hd_obs/datasets/13.43232_juen2014-melt.csv')
        measured_hd_list = [mb_df.hd_m.values]
        measured_melt_list = [mb_df['melt_mf'].values]
        glac_name = "Koxkar (13.43232)"
        fig_fn = '13.43232_hd_melt_juen2014.png'
#        ds_names = ['Juen 2014\n(8/10/10$\u2009$-$\u2009$8/29/10)']
        ds_names = ['8/10/10$\u2009$-$\u2009$8/29/10']
        melt_fp = debris_prms.output_fp + 'ostrem_curves/exp4/'
        melt_fn = '4175N-8000E-debris_melt_curve.nc'
        yearfracs_list = [[2010 + 222/365, 2010 + 241/365]]

        hd_min, hd_max = 0, np.ceil(np.max([x.max() for x in measured_hd_list])/0.1)*0.1 + 0.05
        hd_tick_major, hd_tick_minor = 0.1, 0.02
#        melt_min, melt_max = 0, np.round(np.ceil(np.max([x.max() for x in measured_melt_list]) / 10) * 10,0) + 5
#        melt_tick_major, melt_tick_minor = 10, 5
        melt_min, melt_max = 0, np.round(np.ceil(np.max([x.max() for x in measured_melt_list]) / 0.1) * 0.1,1) + 0.1
        melt_tick_major, melt_tick_minor = 0.5, 0.1
    
        for n in np.arange(0,len(measured_hd_list)):
            assert len(measured_hd_list[n]) == len(measured_melt_list[n]), 'length of hd differs from melt'
        
        plot_hd_vs_melt_comparison(measured_hd_list, measured_melt_list, yearfracs_list, glac_name, fig_fn, 
                                   melt_fp, melt_fn,
                                   ds_names=ds_names, hd_min=hd_min, hd_max=hd_max, 
                                   hd_tick_major=hd_tick_major, hd_tick_minor=hd_tick_minor,
                                   melt_min=melt_min, melt_max=melt_max, 
                                   melt_tick_major=melt_tick_major, melt_tick_minor=melt_tick_minor,
                                   plot_meltfactor=True)
        
#    # ===== No 72 =====
#    if '13.43165' in glaciers:
#        print('\nmelt comparison with Wang et al (2011)')
#        # Data: debris thickness (m) and melt rate (mm w.e. d-1)
#        mb_df = pd.read_csv(debris_prms.main_directory + '/../hd_obs/datasets/13.43165_wang2011-melt.csv')
#        measured_hd_list = [mb_df.hd_m.values]
#        measured_melt_list = [mb_df['melt_mmwed'].values]
#        glac_name = "No 72 (13.43165)"
#        fig_fn = '13.43165_hd_melt_wang2011.png'
#        ds_names = ['8/10/10$\u2009$-$\u2009$8/29/10']
#        melt_fp = debris_prms.output_fp + 'ostrem_curves/exp4/'
#        melt_fn = '4175N-8000E-debris_melt_curve.nc'
#        yearfracs_list = [[2010 + 222/365, 2010 + 241/365]]
#
#        hd_min, hd_max = 0, np.ceil(np.max([x.max() for x in measured_hd_list])/0.1)*0.1 + 0.05
#        hd_tick_major, hd_tick_minor = 0.1, 0.02
##        melt_min, melt_max = 0, np.round(np.ceil(np.max([x.max() for x in measured_melt_list]) / 10) * 10,0) + 5
##        melt_tick_major, melt_tick_minor = 10, 5
#        melt_min, melt_max = 0, np.round(np.ceil(np.max([x.max() for x in measured_melt_list]) / 0.1) * 0.1,1) + 0.1
#        melt_tick_major, melt_tick_minor = 0.5, 0.1
#    
#        print('NEED THE DATES!')
#        for n in np.arange(0,len(measured_hd_list)):
#            assert len(measured_hd_list[n]) == len(measured_melt_list[n]), 'length of hd differs from melt'
#        
#        plot_hd_vs_melt_comparison(measured_hd_list, measured_melt_list, yearfracs_list, glac_name, fig_fn, 
#                                   melt_fp, melt_fn,
#                                   ds_names=ds_names, hd_min=hd_min, hd_max=hd_max, 
#                                   hd_tick_major=hd_tick_major, hd_tick_minor=hd_tick_minor,
#                                   melt_min=melt_min, melt_max=melt_max, 
#                                   melt_tick_major=melt_tick_major, melt_tick_minor=melt_tick_minor)
    
    # ===== Baltoro (14.06794) ====
    if '14.06794' in glaciers:
        print('\nmelt comparison with Mihalcea et al 2006 and Groos et al 2017')
        # Data: debris thickness (m) and melt rate (mm w.e. d-1)
        mb_df = pd.read_csv(debris_prms.main_directory + '/../hd_obs/datasets/14.06794_mihalcea2006-melt.csv')
        measured_hd_list = [mb_df.hd_m.values,
                np.array([0.095, 0., 0.035, 0.022, 0.02, 0.033, 0.088, 0.05, 0.035, 0.056, 0.126, 0.26, 0.13, 0.315, 
                          0.12, 0.375])]
        measured_melt_list = [mb_df['melt_mmwed'].values,
                np.array([2.1, 4.1, 3.8, 4.4, 4.3, 3.9, 3.8, 2.8, 4.6, 3.9, 2.5, 1.9, 2.5, 1.2, 2.9, 1.1]) * 10 * 0.9]
        glac_name = "Baltoro (14.06794)"
        fig_fn = '14.06794_hd_melt_Mih2006_Gro2017.png'
#        ds_names = ['Mihalcea 2006\n(7/04/04$\u2009$-$\u2009$7/14/04)', 'Groos 2017\n(7/22/11$\u2009$-$\u2009$8/10/11)']
        ds_names = ['7/04/04$\u2009$-$\u2009$7/14/04', '7/22/11$\u2009$-$\u2009$8/10/11']
        melt_fp = debris_prms.output_fp + 'ostrem_curves/exp4/'
        melt_fn = '3575N-7650E-debris_melt_curve.nc'
        yearfracs_list = [[2004 + 186/366, 2004 + 196/366],
                          [2011 + 203/365, 2011 + 222/365]]

        hd_min, hd_max = 0, np.ceil(np.max([x.max() for x in measured_hd_list])/0.1)*0.1 + 0.05
        hd_tick_major, hd_tick_minor = 0.1, 0.02
        melt_min, melt_max = 0, np.round(np.ceil(np.max([x.max() for x in measured_melt_list]) / 10) * 10,0) + 5
        melt_tick_major, melt_tick_minor = 10, 5
    
        for n in np.arange(0,len(measured_hd_list)):
            assert len(measured_hd_list[n]) == len(measured_melt_list[n]), 'length of hd differs from melt'
        plot_hd_vs_melt_comparison(measured_hd_list, measured_melt_list, yearfracs_list, glac_name, fig_fn, 
                                   melt_fp, melt_fn,
                                   ds_names=ds_names, hd_min=hd_min, hd_max=hd_max, 
                                   hd_tick_major=hd_tick_major, hd_tick_minor=hd_tick_minor,
                                   melt_min=melt_min, melt_max=melt_max, 
                                   melt_tick_major=melt_tick_major, melt_tick_minor=melt_tick_minor)
        
    # ===== Batal (14.16042) ====
    if '14.16042' in glaciers:
        print('\nmelt comparison with Patel et al (2016)')
        # Data: debris thickness (m) and melt rate (mm w.e. d-1)
        measured_hd_list = [np.array([51, 54, 43, 53, 24, 27, 27, 50, 5, 2, 0, 5, 0, 2]) / 100]
        measured_melt_list = [np.array([5.75, 6.125, 11.875, 7.875, 8.625, 10.375, 9.75, 7.375, 12.875, 19.875, 
                                        18.125, 14.375, 18.625, 10.75])]
        glac_name = 'Batal (14.16042)'
        fig_fn = '14.16042_hd_melt_Patel2016.png'
#        ds_names = ['Patel 2016\n(8/01/14$\u2009$-$\u2009$10/15/14)']
        ds_names = ['8/01/14$\u2009$-$\u2009$10/15/14']
        melt_fp = debris_prms.output_fp + 'ostrem_curves/exp4/'
        melt_fn = '3225N-7750E-debris_melt_curve.nc'
        yearfracs_list = [[2014 + 213/365, 2014 + 288/365]]

        hd_min, hd_max = 0, np.ceil(np.max([x.max() for x in measured_hd_list])/0.1)*0.1 + 0.05
        hd_tick_major, hd_tick_minor = 0.1, 0.02
        melt_min, melt_max = 0, np.round(np.ceil(np.max([x.max() for x in measured_melt_list]) / 10) * 10,0) + 5
        melt_tick_major, melt_tick_minor = 10, 5
    
        for n in np.arange(0,len(measured_hd_list)):
            assert len(measured_hd_list[n]) == len(measured_melt_list[n]), 'length of hd differs from melt'
        plot_hd_vs_melt_comparison(measured_hd_list, measured_melt_list, yearfracs_list, glac_name, fig_fn, 
                                   melt_fp, melt_fn,
                                   ds_names=ds_names, hd_min=hd_min, hd_max=hd_max, 
                                   hd_tick_major=hd_tick_major, hd_tick_minor=hd_tick_minor,
                                   melt_min=melt_min, melt_max=melt_max, 
                                   melt_tick_major=melt_tick_major, melt_tick_minor=melt_tick_minor)
        
    # ===== Khumbu (15.03733) ====
    if '15.03733' in glaciers:
        print('\nmelt comparison with Kayastha et al 2000')
        # Data: debris thickness (m) and melt rate (mm w.e. d-1)
        measured_hd_list = [np.array([0, 0.3, 2, 5, 10, 20, 30, 40]) / 100]
        measured_melt_list = [np.array([27.9, 57.6, 42.3, 29.7, 18, 12.6, 10.8, 9])]
        glac_name = 'Khumbu (15.03733)'
        fig_fn = '15.03733_hd_melt_Kay2000.png'
#        ds_names = ['Kayastha 2000\n(5/22/00$\u2009$-$\u2009$6/01/00)']
        ds_names = ['5/22/00$\u2009$-$\u2009$6/01/00']
        melt_fp = debris_prms.output_fp + 'ostrem_curves/exp4/'
        melt_fn = '2800N-8700E-debris_melt_curve.nc'
        yearfracs_list = [[2000 + 143/366, 2000 + 153/366]]

        hd_min, hd_max = 0, np.ceil(np.max([x.max() for x in measured_hd_list])/0.1)*0.1 + 0.05
        hd_tick_major, hd_tick_minor = 0.1, 0.02
        melt_min, melt_max = 0, np.round(np.ceil(np.max([x.max() for x in measured_melt_list]) / 10) * 10,0) + 5
        melt_tick_major, melt_tick_minor = 10, 5
    
        for n in np.arange(0,len(measured_hd_list)):
            assert len(measured_hd_list[n]) == len(measured_melt_list[n]), 'length of hd differs from melt'
        plot_hd_vs_melt_comparison(measured_hd_list, measured_melt_list, yearfracs_list, glac_name, fig_fn, 
                                   melt_fp, melt_fn,
                                   ds_names=ds_names, hd_min=hd_min, hd_max=hd_max, 
                                   hd_tick_major=hd_tick_major, hd_tick_minor=hd_tick_minor,
                                   melt_min=melt_min, melt_max=melt_max, 
                                   melt_tick_major=melt_tick_major, melt_tick_minor=melt_tick_minor)
        
    # ===== Lirung (15.04045) ====
    if '15.04045' in glaciers:
        print('\nmelt comparison with Chand et al 2015')
        # Data: debris thickness (m) and melt rate (mm w.e. d-1)
        mb_df1 = pd.read_csv(debris_prms.main_directory + '/../hd_obs/datasets/15.04045_chand2015_fall-melt.csv')
        mb_df2 = pd.read_csv(debris_prms.main_directory + '/../hd_obs/datasets/15.04045_chand2015_winter-melt.csv')
        mb_df3 = pd.read_csv(debris_prms.main_directory + '/../hd_obs/datasets/15.04045_chand2015_spring-melt.csv')
        measured_hd_list = [mb_df1.hd_m.values, mb_df2.hd_m.values, mb_df3.hd_m.values]
        measured_melt_list = [mb_df1.melt_mmwed.values, mb_df2.melt_mmwed.values, mb_df3.melt_mmwed.values]
        glac_name = 'Lirung (15.04045)'
        fig_fn = '15.04045_hd_melt_Cha2015.png'
#        ds_names = ['Chand 2000\n(9/22/13$\u2009$-$\u2009$10/03/13)', 'Chand 2000\n(11/29/13$\u2009$-$\u2009$12/12/13)', 
#                    'Chand 2000\n(4/07/14$\u2009$-$\u2009$4/19/14)']
        ds_names = ['9/22/13$\u2009$-$\u2009$10/03/13', '11/29/13$\u2009$-$\u2009$12/12/13', 
                    '4/07/14$\u2009$-$\u2009$4/19/14']
        melt_fp = debris_prms.output_fp + 'ostrem_curves/exp4/'
        melt_fn = '2825N-8550E-debris_melt_curve.nc'
        yearfracs_list = [[2013 + 265/365, 2013 + 276/365], [2013 + 333/365, 2013 + 346/365], 
                          [2014 + 97/365, 2014 + 109/365]]

        hd_min, hd_max = 0, np.ceil(np.max([x.max() for x in measured_hd_list])/0.1)*0.1 + 0.05
        hd_tick_major, hd_tick_minor = 0.1, 0.02
        melt_min, melt_max = 0, np.round(np.ceil(np.max([x.max() for x in measured_melt_list]) / 10) * 10,0) + 5
        melt_tick_major, melt_tick_minor = 10, 5
    
        for n in np.arange(0,len(measured_hd_list)):
            assert len(measured_hd_list[n]) == len(measured_melt_list[n]), 'length of hd differs from melt'
        plot_hd_vs_melt_comparison(measured_hd_list, measured_melt_list, yearfracs_list, glac_name, fig_fn, 
                                   melt_fp, melt_fn,
                                   ds_names=ds_names, hd_min=hd_min, hd_max=hd_max, 
                                   hd_tick_major=hd_tick_major, hd_tick_minor=hd_tick_minor,
                                   melt_min=melt_min, melt_max=melt_max, 
                                   melt_tick_major=melt_tick_major, melt_tick_minor=melt_tick_minor)
        
    # ===== Hailuogou (15.07886) ====
    if '15.07886' in glaciers:
        print('\nmelt comparison with Zhang et al (2011)')
        # Data: debris thickness (m) and melt rate (mm w.e. d-1)
        measured_hd_list = [np.array([2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 5, 5, 6, 7, 7, 10, 10, 11, 13]) / 100]
        measured_melt_list = [np.array([65.2, 55.4, 52.8, 51.6, 47.0, 53.4, 44.4, 50.3, 58, 48.9, 58.4, 54.4, 44.8, 
                                        52.6, 43.7, 52.5, 38.5, 36.5, 34.2, 28.4])]
        glac_name = 'Hailuogou (15.07886)'
        fig_fn = '15.07886_hd_melt_Zhang2011.png'
#        ds_names = ['Zhang 2011\n(7/02/08$\u2009$-$\u2009$9/30/08)']
        ds_names = ['7/02/08$\u2009$-$\u2009$9/30/08']
        melt_fp = debris_prms.output_fp + 'ostrem_curves/exp4/'
        melt_fn = '2950N-10200E-debris_melt_curve.nc'
        yearfracs_list = [[2008 + 184/366, 2008 + 274/366]]

        hd_min, hd_max = 0, np.ceil(np.max([x.max() for x in measured_hd_list])/0.1)*0.1 + 0.05
        hd_tick_major, hd_tick_minor = 0.1, 0.02
        melt_min, melt_max = 0, np.round(np.ceil(np.max([x.max() for x in measured_melt_list]) / 10) * 10,0) + 5
        melt_tick_major, melt_tick_minor = 10, 5
    
        for n in np.arange(0,len(measured_hd_list)):
            assert len(measured_hd_list[n]) == len(measured_melt_list[n]), 'length of hd differs from melt'
        plot_hd_vs_melt_comparison(measured_hd_list, measured_melt_list, yearfracs_list, glac_name, fig_fn, 
                                   melt_fp, melt_fn,
                                   ds_names=ds_names, hd_min=hd_min, hd_max=hd_max, 
                                   hd_tick_major=hd_tick_major, hd_tick_minor=hd_tick_minor,
                                   melt_min=melt_min, melt_max=melt_max, 
                                   melt_tick_major=melt_tick_major, melt_tick_minor=melt_tick_minor)
        
    # ===== 24K (15.11758) ====
    if '15.11758' in glaciers:
        print('\nmelt comparison with Wei et al (2010) and Yang et al (2017)')
        # Data: debris thickness (m) and melt rate (mm w.e. d-1)
        mb_df = pd.read_csv(debris_prms.main_directory + '/../hd_obs/datasets/24k_yang2017-melt.csv')
        mb_df2 = pd.read_csv(debris_prms.main_directory + '/../hd_obs/datasets/24k_wei2010-melt.csv')
        measured_hd_list = [mb_df.hd_m.values, mb_df2.hd_m.values]
        measured_melt_list = [mb_df.melt_mmwed.values, mb_df2.melt_mmwed.values]
        glac_name = "24K (15.11758)"
        fig_fn = '15.11758_hd_melt_Wei2010_Yang2017.png'
#        ds_names = ['Yang 2017\n(6/01/16$\u2009$-$\u2009$9/30/16)', 'Wei et al 2010\n(7/19/08$\u2009$-$\u2009$9/04/08)']
        ds_names = ['6/01/16$\u2009$-$\u2009$9/30/16', '7/19/08$\u2009$-$\u2009$9/04/08']
        melt_fp = debris_prms.output_fp + 'ostrem_curves/exp4/'
        melt_fn = '2975N-9575E-debris_melt_curve.nc'
        yearfracs_list = [[2016 + 153/366, 2016 + 274/366], [2008 + 201/366, 2008 + 248/366]]

        hd_min, hd_max = 0, np.ceil(np.max([x.max() for x in measured_hd_list])/0.1)*0.1 + 0.05
        hd_tick_major, hd_tick_minor = 0.1, 0.02
        melt_min, melt_max = 0, np.round(np.ceil(np.max([x.max() for x in measured_melt_list]) / 10) * 10,0) + 5
        melt_tick_major, melt_tick_minor = 10, 5
    
        for n in np.arange(0,len(measured_hd_list)):
            assert len(measured_hd_list[n]) == len(measured_melt_list[n]), 'length of hd differs from melt'
        
        plot_hd_vs_melt_comparison(measured_hd_list, measured_melt_list, yearfracs_list, glac_name, fig_fn, 
                                   melt_fp, melt_fn,
                                   ds_names=ds_names, hd_min=hd_min, hd_max=hd_max, 
                                   hd_tick_major=hd_tick_major, hd_tick_minor=hd_tick_minor,
                                   melt_min=melt_min, melt_max=melt_max, 
                                   melt_tick_major=melt_tick_major, melt_tick_minor=melt_tick_minor)
        
    # ===== Franz Josef (18.02397) ====
    if '18.02397' in glaciers:
        print('\nmelt comparison with Brook et al (2013)')
        # Data: debris thickness (m) and melt rate (mm w.e. d-1)
        mb_df = pd.read_csv(debris_prms.main_directory + '/../hd_obs/datasets/18.02397_brook2013-melt.csv')
        measured_hd_list = [mb_df.hd_m.values]
        measured_melt_list = [mb_df.melt_mmwed.values]
        glac_name = 'Franz Josef (18.02397)'
        fig_fn = '18.02397_hd_melt_Brook2013.png'
#        ds_names = ['Brook 2013\n(2/07/12$\u2009$-$\u2009$2/16/12)']
        ds_names = ['2/07/12$\u2009$-$\u2009$2/16/12']
        melt_fp = debris_prms.output_fp + 'ostrem_curves/exp4/'
        melt_fn = '4350S-17025E-debris_melt_curve.nc'
        yearfracs_list = [[2012 + 38/366, 2012 + 47/366]]

        hd_min, hd_max = 0, np.ceil(np.max([x.max() for x in measured_hd_list])/0.1)*0.1 + 0.05
        hd_tick_major, hd_tick_minor = 0.1, 0.02
        melt_min, melt_max = 0, np.round(np.ceil(np.max([x.max() for x in measured_melt_list]) / 10) * 10,0) + 5
        melt_tick_major, melt_tick_minor = 10, 5
    
        for n in np.arange(0,len(measured_hd_list)):
            assert len(measured_hd_list[n]) == len(measured_melt_list[n]), 'length of hd differs from melt'
        plot_hd_vs_melt_comparison(measured_hd_list, measured_melt_list, yearfracs_list, glac_name, fig_fn, 
                                   melt_fp, melt_fn,
                                   ds_names=ds_names, hd_min=hd_min, hd_max=hd_max, 
                                   hd_tick_major=hd_tick_major, hd_tick_minor=hd_tick_minor,
                                   melt_min=melt_min, melt_max=melt_max, 
                                   melt_tick_major=melt_tick_major, melt_tick_minor=melt_tick_minor)
    
    

#%%
if option_hd_comparison: 
    
    #%%
#    hd_list = list(np.repeat(0.025, 145))
#    hd_list.extend(list(np.repeat(0.075, 109)))
#    hd_list.extend(list(np.repeat(0.15, 110)))
#    hd_list.extend(list(np.repeat(0.25, 73)))
#    hd_list.extend(list(np.repeat(0.40, 182)))
#    hd_list.extend(list(np.repeat(0.625, 200)))
#    hd_list.extend(list(np.repeat(0.875, 127)))
#    hd_list.extend(list(np.repeat(1.000, 54)))
#    hd_array = np.array(hd_list)
#    print('hd_mean:', np.round(hd_array.mean(),3))
#    print('hd_std: ', np.round(hd_array.std(),3))
#    print('hd_med: ', np.round(np.median(hd_array), 3))
#    print('hd_mad: ', np.round(median_absolute_deviation(hd_array),3))
    #%%
        
    print('\n     !!!!! ADD KOXKAR !!!\n')
    print('\nNeed to add 10.01732 when North Asia has been run \n')
    
    glaciers = ['1.15645', '2.14297', '7.01044', '7.01107', '11.00106', '11.01604', '11.02472', '11.02810', '11.03005', 
                '12.01132', '14.06794', '14.16042', '15.03473', '15.04045']
#    glaciers = ['15.04045']
    

    process_files = True
    bin_width = 50

    if process_files:
        
        hd_datasets_fp = hd_obs_fp + 'datasets/'
        hd_ds_dict = {'1.15645': [(hd_datasets_fp + '1.15645_kennicott_anderson_2019-hd.csv', 'Anderson 2019')], 
                      '2.14297': [(hd_datasets_fp + '2.14297_moore2019-melt.csv', 'Moore 2019')], 
                      '7.01044': [(hd_datasets_fp + '7.01044_lukas2005-hd.csv', 'Lukas 2005')],
                      '7.01107': [(hd_datasets_fp + '7.01107_lukas2005-hd.csv', 'Lukas 2005')],
                      '11.00106': [(hd_datasets_fp + '11.00106_kellerer2008-hd.csv', 'Kellerer 2008')],
                      '11.01604': [(hd_datasets_fp + '11.01604_delGobbo2017-hd.csv', 'Del Gobbo 2017')],
                      '11.02472': [(hd_datasets_fp + '11.02472_bocchiola2015-melt.csv', 'Bocchiola 2015')],
                      '11.02810': [(hd_datasets_fp + '11.02810_reid2012-hd.csv', 'Reid 2012')],
                      '11.03005': [(hd_datasets_fp + '11.03005_foster2012-hd.csv', 'Foster 2012'),
                                   (hd_datasets_fp + '11.03005_mihalcea2008-hd.csv', 'Mihalcea 2008')],
                      '12.01132': [(hd_datasets_fp + '12.01132_popovnin2002_layers-hd.csv', 'Popovnin 2002')],
                      '14.06794': [(hd_datasets_fp + '14.06794_mihalcea2006-melt.csv', 'Mihalcea 2006'), 
                                   (hd_datasets_fp + '14.06794_minora2015-hd.csv', 'Minora 2015')],
                      '14.16042': [(hd_datasets_fp + '14.16042_patel2016-hd.csv', 'Patel 2016')],
                      '15.03473': [(hd_datasets_fp + '15.03473_nicholson2012-hd.csv', 'Nicholson 2012'), 
                                   (hd_datasets_fp + '15.03473_nicholson2017-hd.csv', 'Nicholson 2017'),
                                   (hd_datasets_fp + '15.03473_nicholson2018_gokyo-hd.csv', 'Nicholson 2018'),
                                   (hd_datasets_fp + '15.03473_nicholson2018_margin-hd.csv', 'Nicholson 2018')],
                      '15.04045': [(hd_datasets_fp + '15.04045_mccarthy2017-hd.csv', 'McCarthy 2017')],
                     }
            
        for glacier in glaciers:
            print('\n', glacier)
            for hd_ds_info in hd_ds_dict[glacier]:
                # Load observations
                hd_ds_fn, ds_name = hd_ds_info[0], hd_ds_info[1]
                hd_obs = pd.read_csv(hd_ds_fn)
                
                # remove clean ice values because only comparing to debris-covered areas
                hd_obs = hd_obs[hd_obs['hd_m'] > 0]
                hd_obs.reset_index(inplace=True, drop=True)
                
                # if lat/lon provided, remove "off-glacier" points, i.e., points not covered by debris cover maps
                if 'hd_ts_cal' in hd_obs.columns:
                    hd_obs = hd_obs.dropna(subset=['hd_ts_cal'])
                    hd_obs.reset_index(inplace=True, drop=True)
                
                glac_str = hd_ds_fn.split('/')[-1].split('_')[0]
                reg = int(glac_str.split('.')[0])
                glacno = int(glac_str.split('.')[1])
                
                # Load modeled debris thickness
                try:
                    hdts_fp = debris_prms.mb_binned_fp_wdebris_hdts
                    hdts_fn = glac_str + '_mb_bins_hdts.csv'
                    hdts_df = pd.read_csv(hdts_fp + hdts_fn)
                except:
                    hdts_fp = debris_prms.mb_binned_fp_wdebris_hdts + '../_wdebris_hdts_extrap/'
                    if reg < 10:
                        hdts_fn = (glac_str.split('.')[0].zfill(2) + '.' + glac_str.split('.')[1] + 
                                   '_mb_bins_hdts_extrap.csv')
                    else:
                        hdts_fn = glac_str + '_mb_bins_hdts_extrap.csv'
                    hdts_df = pd.read_csv(hdts_fp + hdts_fn)
    
                hdts_df.loc[:,:] = hdts_df.values.astype(np.float64)
                
                # ===== PROCESS BINNED DATA =====
                if not hd_obs['elev'].isnull().any():
                    # Bins
                    zmin = hd_obs.elev.min()
                    zmax = hd_obs.elev.max()
                    zbincenter_min = hdts_df.loc[0,'bin_center_elev_m']
                    zbincenter_max = hdts_df.loc[hdts_df.shape[0]-1,'bin_center_elev_m']
                    
                    # Find minimum bin
                    while zbincenter_min - bin_width / 2 + bin_width < zmin:
                        zbincenter_min += bin_width
                    # Find maximum bin size
                    while zbincenter_max - bin_width /2 > zmax:
                        zbincenter_max -= bin_width
                        
                    # Compute relevant statistics for each bin
                    hd_compare_all_array = None
                    for nbin, zbincenter in enumerate(np.arange(zbincenter_min, zbincenter_max + bin_width/2, 
                                                                bin_width)):
                        zbincenter_min = zbincenter - bin_width/2
                        zbincenter_max = zbincenter + bin_width/2
                        elev_idx_obs = np.where((hd_obs['elev'].values >= zbincenter_min) & 
                                                (hd_obs['elev'].values < zbincenter_max))[0]
                        obs_count = 0
                        
                        if len(elev_idx_obs) > 0:
                            # Observations
                            if 'hd_m_std' not in hd_obs.columns:
                                hd_obs_subset = hd_obs.loc[elev_idx_obs,'hd_m']
                                hd_bin_med = np.median(hd_obs_subset)
                                hd_bin_mad = np.median(abs(hd_obs_subset - np.median(hd_obs_subset)))
                                obs_count += hd_obs_subset.shape[0]
                            else:
                                hd_bin_med = hd_obs.loc[elev_idx_obs[0],'hd_m_med']
                                hd_bin_mad = hd_obs.loc[elev_idx_obs[0],'hd_m_mad']
                                obs_count = hd_obs.loc[elev_idx_obs[0],'n_obs']
                            
                            # Model
                            hdts_df_elevs = hdts_df['bin_center_elev_m'].values
                            hdts_df_idxs = np.where((hdts_df_elevs >= zbincenter_min) &
                                                    (hdts_df_elevs < zbincenter_max))[0]
                            hdts_list = []
                            for hdts_df_idx in hdts_df_idxs:
                                bin_hdts_center = hdts_df.loc[hdts_df_idx,'hd_ts_med_m']
                                bin_hdts_spread = hdts_df.loc[hdts_df_idx,'hd_ts_mad_m']
                                bin_hdts_count = int(hdts_df.loc[hdts_df_idx,'dc_bin_count_valid'])
                                
                                # Randomly create thicknesses based on center and spread, but ensure mean is similar
                                hdts_list_single = np.random.normal(loc=bin_hdts_center, scale=bin_hdts_spread, 
                                                                    size=(bin_hdts_count))
                                while (abs(np.mean(hdts_list_single) - bin_hdts_center) > 0.005 or 
                                       abs(np.std(hdts_list_single) - bin_hdts_spread) > 0.01):
                                    hdts_list_single = np.random.normal(loc=bin_hdts_center, scale=bin_hdts_spread, 
                                                                        size=(bin_hdts_count))
                                hdts_list.extend(hdts_list_single)
                                
                            hdts_array = np.array(hdts_list)
                            hdts_df_med = np.median(hdts_array)
                            hdts_df_mad = median_absolute_deviation(hdts_array)
                            
                            # Ensure the glacier has data as well
                            if not np.isnan(hdts_df_med):
                                print(glac_str, zbincenter, 
                                      '  ', np.round(hd_bin_med,2), '+/-', np.round(hd_bin_mad,2), 'vs',
                                      '  ', np.round(hdts_df_med,2), '+/-', np.round(hdts_df_mad,2))
                                
                                bin_data = np.array([reg, glacno, zbincenter, obs_count, hd_bin_med, hd_bin_mad, 
                                                     hdts_df_med, hdts_df_mad]).reshape(1,8)
                                
                                if hd_compare_all_array is None:
                                    hd_compare_all_array = bin_data
                                else:
                                    hd_compare_all_array = np.concatenate((hd_compare_all_array, bin_data))
                    
                # ===== PROCESS ENTIRE GLACIER =====
                else:
                    # Observations
                    if hd_obs.shape[0] > 1:
                        hd_bin_med = np.median(hd_obs['hd_m'])
                        hd_bin_mad = median_absolute_deviation(hd_obs['hd_m'].values)
                        obs_count = len(hd_obs['hd_m'].values)
                    else:
                        hd_bin_med = hd_obs.loc[0,'hd_m_med']
                        hd_bin_mad = hd_obs.loc[0,'hd_m_mad']
                        obs_count = hd_obs.loc[0,'n_obs']
                
                    # Modeled debris thickness
                    hdts_idxs = np.where(hdts_df['dc_bin_count_valid'].values > 0)[0]
                    hdts_list = []
                    for hdts_idx in hdts_idxs:
                        bin_hd_center = hdts_df.loc[hdts_idx,'hd_ts_med_m']
                        bin_hd_spread = hdts_df.loc[hdts_idx,'hd_ts_mad_m']
                        bin_hd_count = int(hdts_df.loc[hdts_idx,'dc_bin_count_valid'])
                        
                        # Randomly create thicknesses based on center and spread, but ensure mean is similar
                        hdts_list_single = np.random.normal(loc=bin_hd_center, scale=bin_hd_spread, size=(bin_hd_count))
                        while (abs(np.mean(hdts_list_single) - bin_hd_center) > 0.005 or 
                               abs(np.std(hdts_list_single) - bin_hd_spread) > 0.01):
                            hdts_list_single = np.random.normal(loc=bin_hd_center, scale=bin_hd_spread, 
                                                                size=(bin_hd_count))
                        hdts_list.extend(hdts_list_single)

                    hdts_array = np.array(hdts_list)
                    hdts_df_med = np.median(hdts_array)
                    hdts_df_mad = median_absolute_deviation(hdts_array)
                    
                    zbincenter = ((hdts_df['dc_bin_count_valid'] * hdts_df['bin_center_elev_m']).sum() / 
                                  hdts_df['dc_bin_count_valid'].sum())
                    
                    # Ensure the glacier has data as well
                    if not np.isnan(hdts_df_med):
                        print(glac_str, np.round(zbincenter), 
                              '  ', np.round(hd_bin_med,2), '+/-', np.round(hd_bin_mad,2), 'vs',
                              '  ', np.round(hdts_df_med,2), '+/-', np.round(hdts_df_mad,2))
                        
                        bin_data = np.array([reg, glacno, zbincenter, obs_count, hd_bin_med, hd_bin_mad, 
                                             hdts_df_med, hdts_df_mad]).reshape(1,8)
                        hd_compare_all_array = bin_data
            
                # EXPORT DATASET
                hd_compare_all_cns = ['region', 'glacno', 'zbin', 'obs_count', 'hd_obs_med', 'hd_obs_mad', 
                                      'hd_ts_med_m', 'hd_ts_mad_m']
                hd_compare_all = pd.DataFrame(hd_compare_all_array, columns=hd_compare_all_cns)
                hd_compare_all['hd_ds_name'] = ds_name
                hd_compare_all['glac_str'] = glac_str
                # Filename
                hd_processed_fp = hd_obs_fp + 'hd_processed/'
                hd_processed_fn = hd_ds_fn.split('/')[-1].replace('.csv','-processed.csv')
                if not os.path.exists(hd_processed_fp):
                    os.makedirs(hd_processed_fp)
                hd_compare_all.to_csv(hd_processed_fp + hd_processed_fn, index=False)

        
    #%% 
    # ===== INDIVIDUAL ELEVATION BIN COMPARISONS =====
    def plot_hd_obs_comparison(hd_compare_all_subset, glac_str, fig_fn, 
                               hd_min=0, hd_max=2, hd_tick_major=0.2, hd_tick_minor=0.05):
        
        glac_name_dict = {'1.15645': 'Kennicott',
                          '11.01604':'Suldenferner',
                          '11.02810':'Arolla',
                          '15.04045':'Lirung',
                          '15.03473':'Ngozumpa',
                          '15.03743':'Imja-Lhotse Shar',
                          '14.06794':'Baltoro'}
        
        fig, ax = plt.subplots(1, 1, squeeze=False, gridspec_kw = {'wspace':0, 'hspace':0})
        ds_list = []
        for ndata, region in enumerate(hd_compare_all_subset.region):
            ds_name = hd_compare_all_subset.loc[ndata,'hd_ds_name']
            if ds_name not in ds_list:
                label_str = ds_name + '\n(n = ' + str(int(hd_compare_all_subset['obs_count'].sum())) + ')'
                ds_list.append(ds_name)
            else:
                label_str = None
                
            # Size thresholds
            s_sizes = [10, 20, 40, 80]
            obs_count = hd_compare_all_subset.loc[ndata,'obs_count']
            if obs_count <= 5:
                s_plot = s_sizes[0]
            elif obs_count < 25:
                s_plot = s_sizes[1]
            elif obs_count < 100:
                s_plot = s_sizes[2]
            else:
                s_plot = s_sizes[3]
                
            ax[0,0].scatter(hd_compare_all_subset.loc[ndata,'hd_obs_med'], 
                            hd_compare_all_subset.loc[ndata,'hd_ts_med_m'], 
                            color='k', marker='o', facecolor='none', s=s_plot, zorder=3, label=label_str)
            ax[0,0].errorbar(hd_compare_all_subset.loc[ndata,'hd_obs_med'], 
                             hd_compare_all_subset.loc[ndata,'hd_ts_med_m'], 
                             xerr=1.483*hd_compare_all_subset.loc[ndata,'hd_obs_mad'], 
                             yerr=1.483*hd_compare_all_subset.loc[ndata,'hd_ts_mad_m'], 
                             capsize=1, linewidth=0.5, color='darkgrey', zorder=2)
        # Labels
        ax[0,0].set_xlabel('Observed $h_d$ (m)', size=12)    
        ax[0,0].set_ylabel('Modeled $h_d$ (m)', size=12)
        ax[0,0].set_xlim(hd_min,hd_max)
        ax[0,0].set_ylim(hd_min,hd_max)
        ax[0,0].plot([hd_min, hd_max], [hd_min, hd_max], color='k', linewidth=0.5, zorder=1)
#        ax[0,0].xaxis.set_tick_params(labelsize=12)
        ax[0,0].xaxis.set_major_locator(plt.MultipleLocator(hd_tick_major))
        ax[0,0].xaxis.set_minor_locator(plt.MultipleLocator(hd_tick_minor))  
#        ax[0,0].yaxis.set_tick_params(labelsize=12)
        ax[0,0].yaxis.set_major_locator(plt.MultipleLocator(hd_tick_major))
        ax[0,0].yaxis.set_minor_locator(plt.MultipleLocator(hd_tick_minor))
#        # Tick parameters
#        ax[0,0].yaxis.set_ticks_position('both')
        ax[0,0].tick_params(axis='both', which='major', labelsize=12, direction='inout')
        ax[0,0].tick_params(axis='both', which='minor', labelsize=10, direction='in') 
        # Legend
#        ax[0,0].legend(loc='lower right', ncol=1, fontsize=10, frameon=True, handlelength=1, 
#                       handletextpad=0.15, columnspacing=0.5, borderpad=0.25, labelspacing=0.5)
        fig.text(0.5, 0.9, glac_name_dict[glac_str] + ' (' + glac_str + ')', va='bottom', ha='center', size=12)
        fig.set_size_inches(3.45,3.45)
        fig.savefig(hd_compare_fp + fig_fn, bbox_inches='tight', dpi=300)
    
#    # ------ Kennicott (1.5645)-------
#    if '1.15645' in glaciers:
#        print('\nhd comparison with Anderson et al 2019')
#        glac_str = '1.15645'
#        fig_fn = glac_str + '-hd_bin' + str(bin_width) + 'm_And2019.png'
#        hd_min, hd_max = 0, 0.75
#        hd_tick_major, hd_tick_minor = 0.1, 0.05
#        reg = int(glac_str.split('.')[0])
#        glacno = int(glac_str.split('.')[1])
#        hd_compare_all_subset = hd_compare_all[(hd_compare_all.region == reg) & (hd_compare_all.glacno == glacno)]
#        hd_compare_all_subset.reset_index(inplace=True, drop=True)
#        plot_hd_obs_comparison(hd_compare_all_subset, glac_str, fig_fn, 
#                               hd_min=hd_min, hd_max=hd_max, hd_tick_major=hd_tick_major, hd_tick_minor=hd_tick_minor)

#    # ------ Suldenferner -------
#    glac_str = '11.01604'
##    reg = int(glac_str.split('.')[0])
##    glacno = int(glac_str.split('.')[1])
##    glac_idxs = np.where((hd_compare_all['region'].values == reg) & (hd_compare_all['glacno'] == glacno))[0]
##    hd_compare_all_subset = hd_compare_all[(hd_compare_all.region == reg) & (hd_compare_all.glacno == glacno)]
##    hd_compare_all_subset.reset_index(inplace=True, drop=True)
##    if glac_str in glaciers:
##        fig, ax = plt.subplots(1, 1, squeeze=False, figsize=(10,8), gridspec_kw = {'wspace':0, 'hspace':0})
##        for ndata, region in enumerate(hd_compare_all_subset.region):
##            ax[0,0].scatter(hd_compare_all_subset.loc[ndata,'hd_obs_med'], hd_compare_all_subset.loc[ndata,'hd_ts_med_m'], 
##                            color='k', marker='o', facecolor='none', s=30, zorder=3)
##            ax[0,0].errorbar(hd_compare_all_subset.loc[ndata,'hd_obs_med'], hd_compare_all_subset.loc[ndata,'hd_ts_med_m'], 
##                             xerr=hd_compare_all_subset.loc[ndata,'hd_obs_mad'], 
##                             yerr=hd_compare_all_subset.loc[ndata,'hd_ts_mad_m'], 
##                             capsize=1, linewidth=0.5, color='darkgrey', zorder=2)    
##        ax[0,0].set_xlabel('Observed $h_d$ (m)', size=12)    
##        ax[0,0].set_ylabel('Modeled $h_d$ (m)', size=12)
##        ymin, ymax = 0, 0.75
##        xmin, xmax = 0, 0.75
##        ax[0,0].set_xlim(xmin,xmax)
##        ax[0,0].set_ylim(ymin,ymax)
##        ax[0,0].plot([np.min([xmin,ymin]),np.max([xmax,ymax])], [np.min([xmin,ymin]),np.max([xmax,ymax])], color='k', 
##                     linewidth=0.5, zorder=1)
##        fig.text(0.5, 0.9, glac_name_dict[glac_str] + ' (' + glac_str + ')', va='bottom', ha='center', size=12)
##        fig.set_size_inches(3.45,3.45)
##        fig_fn = glac_str + '-hd_bin' + str(bin_width) + 'm_Nicholson.png'
##        fig.savefig(hd_compare_fp + fig_fn, bbox_inches='tight', dpi=300)
##        
##        
##    # ----- Lirung ------
##    glac_str = '15.04045'
##    reg = int(glac_str.split('.')[0])
##    glacno = int(glac_str.split('.')[1])
##    hd_compare_all_subset = hd_compare_all[(hd_compare_all.region == reg) & (hd_compare_all.glacno == glacno)]
##    hd_compare_all_subset.reset_index(inplace=True, drop=True)
##    if glac_str in glaciers:
##        fig, ax = plt.subplots(1, 1, squeeze=False, figsize=(10,8), gridspec_kw = {'wspace':0, 'hspace':0})
##        for ndata, region in enumerate(hd_compare_all_subset.region):
##            ax[0,0].scatter(hd_compare_all_subset.loc[ndata,'hd_obs_med'], hd_compare_all_subset.loc[ndata,'hd_ts_med_m'], 
##                            color='k', marker='o', facecolor='none', s=30, zorder=3)
##            ax[0,0].errorbar(hd_compare_all_subset.loc[ndata,'hd_obs_med'], hd_compare_all_subset.loc[ndata,'hd_ts_med_m'], 
##                             xerr=hd_compare_all_subset.loc[ndata,'hd_obs_mad'], 
##                             yerr=hd_compare_all_subset.loc[ndata,'hd_ts_mad_m'], 
##                             capsize=1, linewidth=0.5, color='darkgrey', zorder=2)    
##        ax[0,0].set_xlabel('Observed $h_d$ (m)', size=12)    
##        ax[0,0].set_ylabel('Modeled $h_d$ (m)', size=12)
##        ymin, ymax = 0, 2
##        xmin, xmax = 0, 2
##        ax[0,0].set_xlim(xmin,xmax)
##        ax[0,0].set_ylim(ymin,ymax)
##        ax[0,0].plot([np.min([xmin,ymin]),np.max([xmax,ymax])], [np.min([xmin,ymin]),np.max([xmax,ymax])], color='k', 
##                     linewidth=0.5, zorder=1)
##        fig.text(0.5, 0.9, glac_name_dict[glac_str] + ' (' + glac_str + ')', va='bottom', ha='center', size=12)
##        fig.set_size_inches(3.45,3.45)
##        fig_fn = glac_str + '-hd_bin' + str(bin_width) + 'm.png'
##        fig.savefig(hd_compare_fp + fig_fn, bbox_inches='tight', dpi=300)
##    
##    # ----- Ngozumpa ------
##    glac_str = '15.03473'
##    reg = int(glac_str.split('.')[0])
##    glacno = int(glac_str.split('.')[1])
##    hd_compare_all_subset = hd_compare_all[(hd_compare_all.region == reg) & (hd_compare_all.glacno == glacno)]
##    hd_compare_all_subset.reset_index(inplace=True, drop=True)
##    if glac_str in glaciers:
##        fig, ax = plt.subplots(1, 1, squeeze=False, figsize=(10,8), gridspec_kw = {'wspace':0, 'hspace':0})
##        for ndata, region in enumerate(hd_compare_all_subset.region):
##            ax[0,0].scatter(hd_compare_all_subset.loc[ndata,'hd_obs_med'], hd_compare_all_subset.loc[ndata,'hd_ts_med_m'], 
##                            color='k', marker='o', facecolor='none', s=30, zorder=3)
##            ax[0,0].errorbar(hd_compare_all_subset.loc[ndata,'hd_obs_med'], hd_compare_all_subset.loc[ndata,'hd_ts_med_m'], 
##                             xerr=hd_compare_all_subset.loc[ndata,'hd_obs_mad'], 
##                             yerr=hd_compare_all_subset.loc[ndata,'hd_ts_mad_m'], 
##                             capsize=1, linewidth=0.5, color='darkgrey', zorder=2)    
##        ax[0,0].set_xlabel('Observed $h_d$ (m)', size=12)    
##        ax[0,0].set_ylabel('Modeled $h_d$ (m)', size=12)
##        ymin, ymax = 0, 4
##        xmin, xmax = 0, 4
##        ax[0,0].set_xlim(xmin,xmax)
##        ax[0,0].set_ylim(ymin,ymax)
##        ax[0,0].plot([np.min([xmin,ymin]),np.max([xmax,ymax])], [np.min([xmin,ymin]),np.max([xmax,ymax])], color='k', 
##                     linewidth=0.5, zorder=1)
##        fig.text(0.5, 0.9, glac_name_dict[glac_str] + ' (' + glac_str + ')', va='bottom', ha='center', size=12)
##        fig.set_size_inches(3.45,3.45)
##        fig_fn = glac_str + '-hd_bin' + str(bin_width) + 'm.png'
##        fig.savefig(hd_compare_fp + fig_fn, bbox_inches='tight', dpi=300)
##    
##    # ----- Imja-Lhotse Shar ------
##    glac_str = '15.03743'
##    reg = int(glac_str.split('.')[0])
##    glacno = int(glac_str.split('.')[1])
##    hd_compare_all_subset = hd_compare_all[(hd_compare_all.region == reg) & (hd_compare_all.glacno == glacno)]
##    hd_compare_all_subset.reset_index(inplace=True, drop=True)
##    if glac_str in glaciers:
##        fig, ax = plt.subplots(1, 1, squeeze=False, figsize=(10,8), gridspec_kw = {'wspace':0, 'hspace':0})
##        for ndata, region in enumerate(hd_compare_all_subset.region):
##            ax[0,0].scatter(hd_compare_all_subset.loc[ndata,'hd_obs_med'], hd_compare_all_subset.loc[ndata,'hd_ts_med_m'], 
##                            color='k', marker='o', facecolor='none', s=30, zorder=3)
##            ax[0,0].errorbar(hd_compare_all_subset.loc[ndata,'hd_obs_med'], hd_compare_all_subset.loc[ndata,'hd_ts_med_m'], 
##                             xerr=hd_compare_all_subset.loc[ndata,'hd_obs_mad'], 
##                             yerr=hd_compare_all_subset.loc[ndata,'hd_ts_mad_m'], 
##                             capsize=1, linewidth=0.5, color='darkgrey', zorder=2)    
##        ax[0,0].set_xlabel('Observed $h_d$ (m)', size=12)    
##        ax[0,0].set_ylabel('Modeled $h_d$ (m)', size=12)
##        ymin, ymax = 0, 1.5
##        xmin, xmax = 0, 1.5
##        ax[0,0].set_xlim(xmin,xmax)
##        ax[0,0].set_ylim(ymin,ymax)
##        ax[0,0].plot([np.min([xmin,ymin]),np.max([xmax,ymax])], [np.min([xmin,ymin]),np.max([xmax,ymax])], color='k', 
##                     linewidth=0.5, zorder=1)
##        fig.text(0.5, 0.9, glac_name_dict[glac_str] + ' (' + glac_str + ')', va='bottom', ha='center', size=12)
##        fig.set_size_inches(3.45,3.45)
##        fig_fn = glac_str + '-hd_bin' + str(bin_width) + 'm_RM2014.png'
##        fig.savefig(hd_compare_fp + fig_fn, bbox_inches='tight', dpi=300)
##    
##        
##    # ----- Baltoro ------
##    glac_str = '14.06794'
##    reg = int(glac_str.split('.')[0])
##    glacno = int(glac_str.split('.')[1])
##    hd_compare_all_subset = hd_compare_all[(hd_compare_all.region == reg) & (hd_compare_all.glacno == glacno)]
##    hd_compare_all_subset.reset_index(inplace=True, drop=True)
##    if glac_str in glaciers:
##        fig, ax = plt.subplots(1, 1, squeeze=False, figsize=(10,8), gridspec_kw = {'wspace':0, 'hspace':0})
##        for ndata, region in enumerate(hd_compare_all_subset.region):
##            ax[0,0].scatter(hd_compare_all_subset.loc[ndata,'hd_obs_med'], hd_compare_all_subset.loc[ndata,'hd_ts_med_m'], 
##                            color='k', marker='o', facecolor='none', s=30, zorder=3)
##            ax[0,0].errorbar(hd_compare_all_subset.loc[ndata,'hd_obs_med'], hd_compare_all_subset.loc[ndata,'hd_ts_med_m'], 
##                             xerr=hd_compare_all_subset.loc[ndata,'hd_obs_mad'], 
##                             yerr=hd_compare_all_subset.loc[ndata,'hd_ts_mad_m'], 
##                             capsize=1, linewidth=0.5, color='darkgrey', zorder=2)    
##        ax[0,0].set_xlabel('Observed $h_d$ (m)', size=12)    
##        ax[0,0].set_ylabel('Modeled $h_d$ (m)', size=12)
##        ymin, ymax = 0, 0.5
##        xmin, xmax = 0, 0.5
##        ax[0,0].set_xlim(xmin,xmax)
##        ax[0,0].set_ylim(ymin,ymax)
##        ax[0,0].plot([np.min([xmin,ymin]),np.max([xmax,ymax])], [np.min([xmin,ymin]),np.max([xmax,ymax])], color='k', 
##                     linewidth=0.5, zorder=1)
##        fig.text(0.5, 0.9, glac_name_dict[glac_str] + ' (' + glac_str + ')', va='bottom', ha='center', size=12)
##        fig.set_size_inches(3.45,3.45)
##        fig_fn = glac_str + '-hd_bin' + str(bin_width) + 'm_Mih2006.png'
##        fig.savefig(hd_compare_fp + fig_fn, bbox_inches='tight', dpi=300)
#    
#    # ----- ALL TOGETHER -----
##    fig, ax = plt.subplots(1, 1, squeeze=False, figsize=(10,8), gridspec_kw = {'wspace':0, 'hspace':0})
##    glac_str_leg = []
##    for ndata, region in enumerate(hd_compare_all.region):
##        glac_str = str(int(region)) + '.' + str(int(hd_compare_all.loc[ndata,'glacno'])).zfill(5)
##        
##        # Legend
##        if glac_str not in glac_str_leg:
##            glac_str_leg.append(glac_str)
##            leg_label = glac_name_dict[glac_str]
##        else:
##            leg_label = ""
##            
##        ax[0,0].scatter(hd_compare_all.loc[ndata,'hd_obs_med'], hd_compare_all.loc[ndata,'hd_ts_med_m'], 
##                        color='k', 
##                        marker=glac_symbol[glac_str][0],
##                        facecolor=glac_symbol[glac_str][2], 
##                        linewidth=1,
##                        s=glac_symbol[glac_str][1],
##                        label=leg_label,
##                        zorder=3)
##        ax[0,0].errorbar(hd_compare_all.loc[ndata,'hd_obs_med'], hd_compare_all.loc[ndata,'hd_ts_med_m'], 
##                         xerr=hd_compare_all.loc[ndata,'hd_obs_mad'], 
##                         yerr=hd_compare_all.loc[ndata,'hd_ts_mad_m'], 
##                         capsize=1, linewidth=0.5, 
##                         color='darkgrey', 
##                         zorder=2)    
##    ax[0,0].set_xlabel('Observed $h_d$ (m)', size=12)    
##    ax[0,0].set_ylabel('Modeled $h_d$ (m)', size=12)
##    ymin = 0
##    ymax = 2.
##    xmin = 0
##    xmax = 2.
##    ax[0,0].set_xlim(xmin,xmax)
##    ax[0,0].set_ylim(ymin,ymax)
##    ax[0,0].plot([np.min([xmin,ymin]),np.max([xmax,ymax])], [np.min([xmin,ymin]),np.max([xmax,ymax])], color='k', 
##                 linewidth=0.5, zorder=1)
##    # Ensure proper order for legend
##    handles, labels = ax[0,0].get_legend_handles_labels()
##    labels, handles = zip(*sorted(zip(labels, handles), key=lambda t:t[0]))
##    ax[0,0].legend(handles, labels, loc=(0.62,0.05), ncol=1, fontsize=10, frameon=True, handlelength=1, 
##                   handletextpad=0.15, columnspacing=0.5, borderpad=0, labelspacing=0)
##    fig.set_size_inches(3.45,3.45)
##    fig.savefig(hd_compare_fp + 'hd_bin_validation.png', bbox_inches='tight', dpi=300)
#    
##%%
#if option_hd_spatial_compare:
#    compare_suldenferner = False
#    compare_ngozumpa = True
#    
#    # SULDENFERNER
#    if compare_suldenferner:
#        hd_ds = 'nicholson-suld'
#    
#        hd_obs = pd.read_csv(hd_ds_fn_dict[hd_ds])
#        
#        hd_obs[hd_obs['hd_ts_cal'] > 3] = np.nan
#        hd_obs = hd_obs.dropna()
#        hd_obs.reset_index(inplace=True, drop=True)
#        
#        glac_str = '11.01604'
#        
#        # Correlation
##        slope, intercept, r_value, p_value, std_err = linregress(hd_obs['hd_m'].values, hd_obs['hd_ts_cal'].values)
#        
#        # Plot
#        fig, ax = plt.subplots(1, 1, squeeze=False, figsize=(10,8), gridspec_kw = {'wspace':0, 'hspace':0})
#
#        ax[0,0].scatter(hd_obs['hd_m'], hd_obs['hd_ts_cal'], 
#                        color='k', marker='o', facecolor='none', s=30, zorder=3)
# 
#        ax[0,0].set_xlabel('Observed $h_d$ (m)', size=12)    
#        ax[0,0].set_ylabel('Modeled $h_d$ (m)', size=12)
#        ymin, ymax = 0, 0.75
#        xmin, xmax = 0, 0.75
#        ax[0,0].set_xlim(xmin,xmax)
#        ax[0,0].set_ylim(ymin,ymax)
#        ax[0,0].plot([np.min([xmin,ymin]),np.max([xmax,ymax])], [np.min([xmin,ymin]),np.max([xmax,ymax])], color='k', 
#                     linewidth=0.5, zorder=1)
#        fig.text(0.5, 0.9, glac_name_dict[glac_str] + ' (' + glac_str + ')', va='bottom', ha='center', size=12)
#        fig.set_size_inches(3.45,3.45)
#        fig_fn = glac_str + '-hd_pts_Nicholson.png'
#        fig.savefig(hd_compare_fp + fig_fn, bbox_inches='tight', dpi=300)
#        
#    
#    # SULDENFERNER
#    if compare_ngozumpa:
#        hd_ds = 'mccarthy-ngoz'
#    
#        hd_obs = pd.read_csv(hd_ds_fn_dict[hd_ds])
#        
#        hd_obs[hd_obs['hd_ts_cal'] > 3] = np.nan
#        hd_obs = hd_obs.dropna()
#        hd_obs.reset_index(inplace=True, drop=True)
#        
#        glac_str = '15.03473'
#        
#        # Correlation
##        slope, intercept, r_value, p_value, std_err = linregress(hd_obs['hd_m'].values, hd_obs['hd_ts_cal'].values)
#        
#        # Plot
#        fig, ax = plt.subplots(1, 1, squeeze=False, figsize=(10,8), gridspec_kw = {'wspace':0, 'hspace':0})
#
#        ax[0,0].scatter(hd_obs['hd_m'], hd_obs['hd_ts_cal'], 
##                        color='k', marker='o', facecolor='none', s=30, zorder=3
#                        color='k', marker='o', linewidth=0.1, facecolor='none', s=3, zorder=3)
# 
#        ax[0,0].set_xlabel('Observed $h_d$ (m)', size=12)    
#        ax[0,0].set_ylabel('Modeled $h_d$ (m)', size=12)
#        ymin, ymax = 0, 7.5
#        xmin, xmax = 0, 7.5
#        ax[0,0].set_xlim(xmin,xmax)
#        ax[0,0].set_ylim(ymin,ymax)
#        ax[0,0].plot([np.min([xmin,ymin]),np.max([xmax,ymax])], [np.min([xmin,ymin]),np.max([xmax,ymax])], color='k', 
#                     linewidth=0.5, zorder=1)
#        fig.text(0.5, 0.9, glac_name_dict[glac_str] + ' (' + glac_str + ')', va='bottom', ha='center', size=12)
#        fig.set_size_inches(3.45,3.45)
#        fig_fn = glac_str + '-hd_pts_McCarthyGPR.png'
#        fig.savefig(hd_compare_fp + fig_fn, bbox_inches='tight', dpi=300)
        