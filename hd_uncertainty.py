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

#%%% ===== SCRIPT OPTIONS =====
option_hd_melt_uncertainty = True
option_melt_diagram_template = False

#hd_obs_fp = debris_prms.main_directory + '/../hd_obs/'


#%% ===== FUNCTIONS =====
def hd_melt_uncertainty(measured_hd_list, yearfracs_list, melt_fp, melt_fn, z_value = 1.645):
    """ Calculate hd-melt relationship for uncertainty  """
    # Dataset of melt data
    ds_ostrem = xr.open_dataset(melt_fp + melt_fn)
    ds_ostrem = ds_ostrem.sortby('hd_cm')
    
    time_year = pd.to_datetime(ds_ostrem.time.values).year
    time_daysperyear = np.array([366 if x%4 == 0 else 365 for x in time_year])
    time_yearfrac = time_year + (pd.to_datetime(ds_ostrem.time.values).dayofyear-1) / time_daysperyear

    hd_wbnds_array_list = []
    for n in np.arange(0,len(measured_hd_list)):
        yearfracs = yearfracs_list[n]
        start_yearfrac = yearfracs[0]
        end_yearfrac = yearfracs[1]
#    # Hack to just have one curve per glacier from 2000 - 2015
#    for n in [0]:
#        start_yearfrac = 2000.6
#        end_yearfrac = 2018.6
        
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
            
        # MEAN CURVE
        fit_idx = list(np.where(debris_thicknesses >= 5)[0])            
        func_coeff, pcov = curve_fit(melt_fromdebris_func, 
                                     debris_melt_df.debris_thickness.values[fit_idx], 
                                     debris_melt_df.melt_mmwed.values[fit_idx])
#        melt_cleanice = debris_melt_df.loc[0,'melt_mmwed']
        
        # LOWER BOUND CURVE
        func_coeff_bndlow, pcov = curve_fit(melt_fromdebris_func, 
                                            debris_melt_df.debris_thickness.values[fit_idx], 
                                            debris_melt_df.melt_bndlow_mmwed.values[fit_idx])
#        melt_cleanice_bndlow = debris_melt_df.loc[0,'melt_bndlow_mmwed']
        
        # UPPER BOUND CURVE
        func_coeff_bndhigh, pcov = curve_fit(melt_fromdebris_func, 
                                            debris_melt_df.debris_thickness.values[fit_idx], 
                                            debris_melt_df.melt_bndhigh_mmwed.values[fit_idx])
#        melt_cleanice_bndhigh = debris_melt_df.loc[0,'melt_bndhigh_mmwed']
        
        debris_4curve = np.arange(0.02,3.01,0.01)
        # column 0 = hd
        # column 1 = melt
        # column 2 = melt bndlow
        # column 3 = melt bndhigh
        # column 4 = hd bndlow debris properties
        # column 5 = hd bndhigh debris properties
        # column 6 = hd bndlow elevchg
        # column 7 = hd bndhigh elevchg
        # column 8 = hd bndlow combined
        # column 9 = hd bndhigh combined
        hd_wbnds_array = np.zeros((len(debris_4curve), 10))
        for ndebris, hd in enumerate(debris_4curve):
            
            # Invert melt against bounded curves to get the uncertainty
            melt_mean = melt_fromdebris_func(hd, func_coeff[0], func_coeff[1])
            hd_low = debris_frommelt_func(melt_mean, func_coeff_bndlow[0], func_coeff_bndlow[1])
            if hd_low < 0:
                hd_low = 0
            hd_high = debris_frommelt_func(melt_mean, func_coeff_bndhigh[0], func_coeff_bndhigh[1])

            # Increase/decrease melt based on elevation change uncertainty and get change in debris thickness
            melt_bndlow = melt_mean + elevchg_mwea_zadj
            melt_bndhigh = melt_mean - elevchg_mwea_zadj
            hd_bndlow_elevchg = debris_frommelt_func(melt_bndlow, func_coeff[0], func_coeff[1])
            if hd_bndlow_elevchg < 0:
                hd_bndlow_elevchg = 0
            hd_bndhigh_elevchg = debris_frommelt_func(melt_bndhigh, func_coeff[0], func_coeff[1])
            
            # Combined root sum of squares of deviations
            hd_bndlow_both = hd - ((hd - hd_low)**2 + (hd - hd_bndlow_elevchg)**2)**0.5
            hd_bndhigh_both = hd + ((hd - hd_high)**2 + (hd - hd_bndhigh_elevchg)**2)**0.5
            
            # Max combined
#            hd_bndlow_max = debris_frommelt_func(melt_bndlow, func_coeff_bndlow[0], func_coeff_bndlow[1])
#            if hd_bndlow_max < 0:
#                hd_bndlow_max = 0
#            hd_bndhigh_max = debris_frommelt_func(melt_bndhigh, func_coeff_bndhigh[0], func_coeff_bndhigh[1])
            
            # Record data
            hd_wbnds_array[ndebris,:] = [hd, melt_mean, melt_bndlow, melt_bndhigh, hd_low, hd_high, 
                                         hd_bndlow_elevchg, hd_bndhigh_elevchg, hd_bndlow_both, hd_bndhigh_both]
            
#            print(np.round(hd,2), ' melt:', np.round(melt_mean,2), 
#                  'bnds:', str(np.round(hd_low,2)) + '-' + str(np.round(hd_high,2)),
#                  '  bndelev:', str(np.round(hd_bndlow_elevchg,2)) + '-' + str(np.round(hd_bndhigh_elevchg,2)), 
#                  '  bndboth:', str(np.round(hd_bndlow_both,2)) + '-' + str(np.round(hd_bndhigh_both,2)),
##                  '  bndmax:', str(np.round(hd_bndlow_max,2)) + '-' + str(np.round(hd_bndhigh_max,2))
#                  )
        hd_wbnds_array_list.append(hd_wbnds_array)
            
    return hd_wbnds_array_list


#%%
if option_hd_melt_uncertainty:
#    glaciers = ['1.15645', '2.14297', '6.00474', '7.01044', '10.01732', '11.00719', '11.02472', '11.02810', '11.02858',
#                '11.03005', '12.01012', '12.01132', '13.05000', '13.43232', '14.06794', '14.16042', '15.03733', 
#                '15.03743', '15.04045', '15.07886', '15.11758', '18.02397']
    glaciers = ['1.15645', '2.14297', '6.00474', '7.01044', '11.00719', '11.02472', '11.02810', '11.02858', '11.03005', 
                '12.01012', '12.01132', '13.05000', '13.43232', '14.06794', '14.16042', '15.03733', '15.03743', 
                '15.04045', '15.07886', '15.11758', '18.02397']
#    glaciers = ['10.01732']
    
    print('add comparison to Haidong etal 2006')

#    z_value = 1.645 # 5-95%
    z_value = 1     # 16-84%
#    z_value = 0.675 # 25-75%
    
    elevchg_mwea_std = 0.72
    elevchg_mwea_zadj = z_value * elevchg_mwea_std
    
    hd_wbnds_array_all = None
    # ===== KENNICOTT (1.15645) ====
    if '1.15645' in glaciers:
        # Data: debris thickness (m) and melt rate (mm w.e. d-1)
        mb_df = pd.read_csv(debris_prms.main_directory + '/../hd_obs/datasets/1.15645_kennicott_anderson_2019-melt.csv')
        measured_hd_list = [mb_df.hd_m.values]
        measured_melt_list = [mb_df.melt_mmwed.values]
        melt_fp = debris_prms.output_fp + 'ostrem_curves/exp4/'
        melt_fn = '6150N-21700E-debris_melt_curve.nc'
        yearfracs_list = [[2011 + 169/365, 2011 + 228/365]]
    
        for n in np.arange(0,len(measured_hd_list)):
            assert len(measured_hd_list[n]) == len(measured_melt_list[n]), 'length of hd differs from melt'

        hd_wbnds_array_list = hd_melt_uncertainty(measured_hd_list, yearfracs_list, melt_fp, melt_fn, z_value=z_value)
        
        for hd_wbnds_array in hd_wbnds_array_list:
            if hd_wbnds_array_all is None:
                hd_wbnds_array_all = hd_wbnds_array[:,:,np.newaxis]
            else:
                hd_wbnds_array_all = np.concatenate((hd_wbnds_array_all, hd_wbnds_array[:,:,np.newaxis]), axis=2)

    # ===== Emmons (2.14297) ====
    if '2.14297' in glaciers:
        # Data: debris thickness (m) and melt rate (mm w.e. d-1)
        mb_df = pd.read_csv(debris_prms.main_directory + '/../hd_obs/datasets/2.14297_moore2019-melt.csv')
        measured_hd_list = [mb_df.hd_m.values]
        measured_melt_list = [mb_df.melt_mmwed.values]
        melt_fp = debris_prms.output_fp + 'ostrem_curves/exp4/'
        melt_fn = '4700N-23825E-debris_melt_curve.nc'
        yearfracs_list = [[2014 + 212/365, 2014 + 222/365]]
    
        for n in np.arange(0,len(measured_hd_list)):
            assert len(measured_hd_list[n]) == len(measured_melt_list[n]), 'length of hd differs from melt'
        
        hd_wbnds_array_list = hd_melt_uncertainty(measured_hd_list, yearfracs_list, melt_fp, melt_fn, z_value=z_value)
        
        for hd_wbnds_array in hd_wbnds_array_list:
            if hd_wbnds_array_all is None:
                hd_wbnds_array_all = hd_wbnds_array[:,:,np.newaxis]
            else:
                hd_wbnds_array_all = np.concatenate((hd_wbnds_array_all, hd_wbnds_array[:,:,np.newaxis]), axis=2)
    
    # ===== Svinafellsjokull (06.00474) ====
    if '6.00474' in glaciers:
        # Data: debris thickness (m) and melt rate (mm w.e. d-1)
        mb_df = pd.read_csv(debris_prms.main_directory + '/../hd_obs/datasets/6.00474_moller2016-melt.csv')
        measured_hd_list = [mb_df.hd_m.values]
        measured_melt_list = [mb_df['melt_mf'].values]
        melt_fp = debris_prms.output_fp + 'ostrem_curves/exp4/'
        melt_fn = '6400N-34325E-debris_melt_curve.nc'
        yearfracs_list = [[2013 + 137/365, 2013 + 150/365]]
    
        for n in np.arange(0,len(measured_hd_list)):
            assert len(measured_hd_list[n]) == len(measured_melt_list[n]), 'length of hd differs from melt'
            
        hd_wbnds_array_list = hd_melt_uncertainty(measured_hd_list, yearfracs_list, melt_fp, melt_fn, z_value=z_value)
        
        for hd_wbnds_array in hd_wbnds_array_list:
            if hd_wbnds_array_all is None:
                hd_wbnds_array_all = hd_wbnds_array[:,:,np.newaxis]
            else:
                hd_wbnds_array_all = np.concatenate((hd_wbnds_array_all, hd_wbnds_array[:,:,np.newaxis]), axis=2)
        
    # ===== Larsbreen (7.01044) ====
    if '7.01044' in glaciers:
        # Data: debris thickness (m) and melt rate (mm w.e. d-1)
        mb_df = pd.read_csv(debris_prms.main_directory + '/../hd_obs/datasets/7.01044_larsbreen_NB2006-melt.csv')
        measured_hd_list = [mb_df.hd_m.values]
        measured_melt_list = [mb_df.melt_mmwed.values]
        melt_fp = debris_prms.output_fp + 'ostrem_curves/exp4/'
        melt_fn = '7825N-1600E-debris_melt_curve.nc'
        yearfracs_list = [[2002 + 191/366, 2002 + 202/366]]
    
        for n in np.arange(0,len(measured_hd_list)):
            assert len(measured_hd_list[n]) == len(measured_melt_list[n]), 'length of hd differs from melt'
        
        hd_wbnds_array_list = hd_melt_uncertainty(measured_hd_list, yearfracs_list, melt_fp, melt_fn, z_value=z_value)
        
        for hd_wbnds_array in hd_wbnds_array_list:
            if hd_wbnds_array_all is None:
                hd_wbnds_array_all = hd_wbnds_array[:,:,np.newaxis]
            else:
                hd_wbnds_array_all = np.concatenate((hd_wbnds_array_all, hd_wbnds_array[:,:,np.newaxis]), axis=2)
        
    # ===== Maliy Aktru (10.01732) ====
    if '10.01732' in glaciers:
#        print('\nmelt comparison with Mayer et al (2011)')
        assert True == False, '10.01732 NEEDS TO DO THE MODELING FIRST!'
        # Data: debris thickness (m) and melt rate (mm w.e. d-1)
        mb_df = pd.read_csv(debris_prms.main_directory + '/../hd_obs/datasets/10.01732_mayer2011-melt.csv')
        measured_hd_list = [mb_df.hd_m.values]
        measured_melt_list = [mb_df['melt_mf'].values]
        melt_fp = debris_prms.output_fp + 'ostrem_curves/exp4/'
        melt_fn = '5000N-8775E-debris_melt_curve.nc'
        yearfracs_list = [[2007 + 192/365, 2007 + 211/365]]
    
        for n in np.arange(0,len(measured_hd_list)):
            assert len(measured_hd_list[n]) == len(measured_melt_list[n]), 'length of hd differs from melt'
        
        hd_wbnds_array_list = hd_melt_uncertainty(measured_hd_list, yearfracs_list, melt_fp, melt_fn, z_value=z_value)
        
        for hd_wbnds_array in hd_wbnds_array_list:
            if hd_wbnds_array_all is None:
                hd_wbnds_array_all = hd_wbnds_array[:,:,np.newaxis]
            else:
                hd_wbnds_array_all = np.concatenate((hd_wbnds_array_all, hd_wbnds_array[:,:,np.newaxis]), axis=2)

    # ===== Vernagtferner (11.00719) ====
    if '11.00719' in glaciers:
        # Data: debris thickness (m) and melt rate (mm w.e. d-1)
        mb_df = pd.read_csv(debris_prms.main_directory + '/../hd_obs/datasets/11.00719_vernagtferner_juen2013-melt.csv')
        measured_hd_list = [mb_df.hd_m.values]
        measured_melt_list = [mb_df.melt_mmwed.values]
        melt_fp = debris_prms.output_fp + 'ostrem_curves/exp4/'
        melt_fn = '4700N-1075E-debris_melt_curve.nc'
        yearfracs_list = [[2010 + 176/365, 2010 + 191/365]]
    
        for n in np.arange(0,len(measured_hd_list)):
            assert len(measured_hd_list[n]) == len(measured_melt_list[n]), 'length of hd differs from melt'
        
        hd_wbnds_array_list = hd_melt_uncertainty(measured_hd_list, yearfracs_list, melt_fp, melt_fn, z_value=z_value)
        
        for hd_wbnds_array in hd_wbnds_array_list:
            if hd_wbnds_array_all is None:
                hd_wbnds_array_all = hd_wbnds_array[:,:,np.newaxis]
            else:
                hd_wbnds_array_all = np.concatenate((hd_wbnds_array_all, hd_wbnds_array[:,:,np.newaxis]), axis=2)
        
    # ===== Vernocolo (11.02472) =====
    if '11.02472' in glaciers:
        # Data: debris thickness (m) and melt rate (mm w.e. d-1)
        mb_df = pd.read_csv(debris_prms.main_directory + '/../hd_obs/datasets/11.02472_bocchiola2015-melt.csv')
        measured_hd_list = [mb_df.hd_m.values]
        measured_melt_list = [mb_df.melt_mmwed.values]
        melt_fp = debris_prms.output_fp + 'ostrem_curves/exp4/'
        melt_fn = '4625N-1050E-debris_melt_curve.nc'
        yearfracs_list = [[2007 + 222/365, 2007 + 256/365]]
    
        for n in np.arange(0,len(measured_hd_list)):
            assert len(measured_hd_list[n]) == len(measured_melt_list[n]), 'length of hd differs from melt'
        
        hd_wbnds_array_list = hd_melt_uncertainty(measured_hd_list, yearfracs_list, melt_fp, melt_fn, z_value=z_value)
        
        for hd_wbnds_array in hd_wbnds_array_list:
            if hd_wbnds_array_all is None:
                hd_wbnds_array_all = hd_wbnds_array[:,:,np.newaxis]
            else:
                hd_wbnds_array_all = np.concatenate((hd_wbnds_array_all, hd_wbnds_array[:,:,np.newaxis]), axis=2)
        
    # ===== Arolla (11.02810) ====
    if '11.02810' in glaciers:
        # Data: debris thickness (m) and melt rate (mm w.e. d-1)
        mb_df = pd.read_csv(debris_prms.main_directory + '/../hd_obs/datasets/11.02810_arolla_reid2012-melt.csv')
        measured_hd_list = [mb_df.hd_m.values]
        measured_melt_list = [mb_df.melt_mmwed.values]
        melt_fp = debris_prms.output_fp + 'ostrem_curves/exp4/'
        melt_fn = '4600N-750E-debris_melt_curve.nc'
        yearfracs_list = [[2010 + 209/365, 2010 + 252/365]]
    
        for n in np.arange(0,len(measured_hd_list)):
            assert len(measured_hd_list[n]) == len(measured_melt_list[n]), 'length of hd differs from melt'
        
        hd_wbnds_array_list = hd_melt_uncertainty(measured_hd_list, yearfracs_list, melt_fp, melt_fn, z_value=z_value)
        
        for hd_wbnds_array in hd_wbnds_array_list:
            if hd_wbnds_array_all is None:
                hd_wbnds_array_all = hd_wbnds_array[:,:,np.newaxis]
            else:
                hd_wbnds_array_all = np.concatenate((hd_wbnds_array_all, hd_wbnds_array[:,:,np.newaxis]), axis=2)
        
    # ===== Belvedere (11.02858) ====
    if '11.02858' in glaciers:
        # Data: debris thickness (m) and melt rate (mm w.e. d-1)
        mb_df = pd.read_csv(debris_prms.main_directory + '/../hd_obs/datasets/11.02858_belvedere_nb2006-melt.csv')
        measured_hd_list = [mb_df.hd_m.values]
        measured_melt_list = [mb_df.melt_mmwed.values]
        melt_fp = debris_prms.output_fp + 'ostrem_curves/exp4/'
        melt_fn = '4600N-800E-debris_melt_curve.nc'
        yearfracs_list = [[2003 + 218/365, 2003 + 222/365]]
    
        for n in np.arange(0,len(measured_hd_list)):
            assert len(measured_hd_list[n]) == len(measured_melt_list[n]), 'length of hd differs from melt'
        
        hd_wbnds_array_list = hd_melt_uncertainty(measured_hd_list, yearfracs_list, melt_fp, melt_fn, z_value=z_value)
        
        for hd_wbnds_array in hd_wbnds_array_list:
            if hd_wbnds_array_all is None:
                hd_wbnds_array_all = hd_wbnds_array[:,:,np.newaxis]
            else:
                hd_wbnds_array_all = np.concatenate((hd_wbnds_array_all, hd_wbnds_array[:,:,np.newaxis]), axis=2)
        
    # ===== MIAGE (11.03005) ====
    if '11.03005' in glaciers:
        # Data: debris thickness (m) and melt rate (mm w.e. d-1)
        mb_df = pd.read_csv(debris_prms.main_directory + '/../hd_obs/datasets/11.03005_reid2010-melt.csv')
        measured_hd_list = [mb_df.hd_m.values]
        measured_melt_list = [mb_df.melt_mmwed.values]
        melt_fp = debris_prms.output_fp + 'ostrem_curves/exp4/'
        melt_fn = '4650N-1050E-debris_melt_curve.nc'
        yearfracs_list = [[2005 + 172/365, 2005 + 247/365]]
    
        for n in np.arange(0,len(measured_hd_list)):
            assert len(measured_hd_list[n]) == len(measured_melt_list[n]), 'length of hd differs from melt'
        
        hd_wbnds_array_list = hd_melt_uncertainty(measured_hd_list, yearfracs_list, melt_fp, melt_fn, z_value=z_value)
        
        for hd_wbnds_array in hd_wbnds_array_list:
            if hd_wbnds_array_all is None:
                hd_wbnds_array_all = hd_wbnds_array[:,:,np.newaxis]
            else:
                hd_wbnds_array_all = np.concatenate((hd_wbnds_array_all, hd_wbnds_array[:,:,np.newaxis]), axis=2)
        
    # ===== Zopkhito (12.01012) ====
    if '12.01012' in glaciers:
        # Data: debris thickness (m) and melt rate (mm w.e. d-1)
        mb_df = pd.read_csv(debris_prms.main_directory + '/../hd_obs/datasets/12.01012_lambrecht2011-melt2008.csv')
        mb_df2 = pd.read_csv(debris_prms.main_directory + '/../hd_obs/datasets/12.01012_lambrecht2011-melt2009.csv')
        measured_hd_list = [mb_df.hd_m.values, mb_df2.hd_m.values]
        measured_melt_list = [mb_df.melt_mmwed.values, mb_df2.melt_mmwed.values]
#        measured_melt_list = [mb_df['melt_mf'].values, mb_df2['melt_mf'].values]
        melt_fp = debris_prms.output_fp + 'ostrem_curves/exp4/'
        melt_fn = '4300N-4350E-debris_melt_curve.nc'
        yearfracs_list = [[2008 + 172/366, 2008 + 179/366], [2009 + 182/365, 2009 + 189/365]]
    
        for n in np.arange(0,len(measured_hd_list)):
            assert len(measured_hd_list[n]) == len(measured_melt_list[n]), 'length of hd differs from melt'
        
        hd_wbnds_array_list = hd_melt_uncertainty(measured_hd_list, yearfracs_list, melt_fp, melt_fn, z_value=z_value)
        
        for hd_wbnds_array in hd_wbnds_array_list:
            if hd_wbnds_array_all is None:
                hd_wbnds_array_all = hd_wbnds_array[:,:,np.newaxis]
            else:
                hd_wbnds_array_all = np.concatenate((hd_wbnds_array_all, hd_wbnds_array[:,:,np.newaxis]), axis=2)
        
    # ===== Djankuat (12.01132) ====
    if '12.01132' in glaciers:
        # Data: debris thickness (m) and melt rate (mm w.e. d-1)
        mb_df = pd.read_csv(debris_prms.main_directory + '/../hd_obs/datasets/12.01132_lambrecht2011-melt.csv')
        measured_hd_list = [mb_df.hd_m.values]
        measured_melt_list = [mb_df['melt_mf'].values]
        melt_fp = debris_prms.output_fp + 'ostrem_curves/exp4/'
        melt_fn = '4325N-4275E-debris_melt_curve.nc'
        yearfracs_list = [[2008 + 172/366, 2008 + 246/366]]
    
        for n in np.arange(0,len(measured_hd_list)):
            assert len(measured_hd_list[n]) == len(measured_melt_list[n]), 'length of hd differs from melt'
        
        hd_wbnds_array_list = hd_melt_uncertainty(measured_hd_list, yearfracs_list, melt_fp, melt_fn, z_value=z_value)
        
        for hd_wbnds_array in hd_wbnds_array_list:
            if hd_wbnds_array_all is None:
                hd_wbnds_array_all = hd_wbnds_array[:,:,np.newaxis]
            else:
                hd_wbnds_array_all = np.concatenate((hd_wbnds_array_all, hd_wbnds_array[:,:,np.newaxis]), axis=2)
        
    # ===== S Inylchek (13.05000) ====
    if '13.05000' in glaciers:
        # Data: debris thickness (m) and melt rate (mm w.e. d-1)
        mb_df = pd.read_csv(debris_prms.main_directory + '/../hd_obs/datasets/13.05000_hagg2008-melt.csv')
        measured_hd_list = [mb_df.hd_m.values]
        measured_melt_list = [mb_df.melt_mmwed.values]
        melt_fp = debris_prms.output_fp + 'ostrem_curves/exp4/'
        melt_fn = '4200N-8025E-debris_melt_curve.nc'
        yearfracs_list = [[2005 + 211/365, 2005 + 222/365]]
    
        for n in np.arange(0,len(measured_hd_list)):
            assert len(measured_hd_list[n]) == len(measured_melt_list[n]), 'length of hd differs from melt'
        
        hd_wbnds_array_list = hd_melt_uncertainty(measured_hd_list, yearfracs_list, melt_fp, melt_fn, z_value=z_value)
        
        for hd_wbnds_array in hd_wbnds_array_list:
            if hd_wbnds_array_all is None:
                hd_wbnds_array_all = hd_wbnds_array[:,:,np.newaxis]
            else:
                hd_wbnds_array_all = np.concatenate((hd_wbnds_array_all, hd_wbnds_array[:,:,np.newaxis]), axis=2)
        
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
        
    # ===== Koxkar (13.43232) ====
    if '13.43232' in glaciers:
        # Data: debris thickness (m) and melt rate (mm w.e. d-1)
        mb_df = pd.read_csv(debris_prms.main_directory + '/../hd_obs/datasets/13.43232_juen2014-melt.csv')
        measured_hd_list = [mb_df.hd_m.values]
        measured_melt_list = [mb_df['melt_mf'].values]
        melt_fp = debris_prms.output_fp + 'ostrem_curves/exp4/'
        melt_fn = '4175N-8000E-debris_melt_curve.nc'
        yearfracs_list = [[2010 + 222/365, 2010 + 241/365]]
    
        for n in np.arange(0,len(measured_hd_list)):
            assert len(measured_hd_list[n]) == len(measured_melt_list[n]), 'length of hd differs from melt'
        
        hd_wbnds_array_list = hd_melt_uncertainty(measured_hd_list, yearfracs_list, melt_fp, melt_fn, z_value=z_value)
        
        for hd_wbnds_array in hd_wbnds_array_list:
            if hd_wbnds_array_all is None:
                hd_wbnds_array_all = hd_wbnds_array[:,:,np.newaxis]
            else:
                hd_wbnds_array_all = np.concatenate((hd_wbnds_array_all, hd_wbnds_array[:,:,np.newaxis]), axis=2)
    
    # ===== Baltoro (14.06794) ====
    if '14.06794' in glaciers:
        # Data: debris thickness (m) and melt rate (mm w.e. d-1)
        mb_df = pd.read_csv(debris_prms.main_directory + '/../hd_obs/datasets/14.06794_mihalcea2006-melt.csv')
        mb_df2 = pd.read_csv(debris_prms.main_directory + '/../hd_obs/datasets/14.06794_groos2017-melt.csv')
        measured_hd_list = [mb_df.hd_m.values, mb_df2.hd_m.values]
        measured_melt_list = [mb_df.melt_mmwed.values, mb_df2.melt_mmwed.values]
        melt_fp = debris_prms.output_fp + 'ostrem_curves/exp4/'
        melt_fn = '3575N-7650E-debris_melt_curve.nc'
        yearfracs_list = [[2004 + 186/366, 2004 + 196/366],
                          [2011 + 203/365, 2011 + 222/365]]

        for n in np.arange(0,len(measured_hd_list)):
            assert len(measured_hd_list[n]) == len(measured_melt_list[n]), 'length of hd differs from melt'
        
        hd_wbnds_array_list = hd_melt_uncertainty(measured_hd_list, yearfracs_list, melt_fp, melt_fn, z_value=z_value)
        
        for hd_wbnds_array in hd_wbnds_array_list:
            if hd_wbnds_array_all is None:
                hd_wbnds_array_all = hd_wbnds_array[:,:,np.newaxis]
            else:
                hd_wbnds_array_all = np.concatenate((hd_wbnds_array_all, hd_wbnds_array[:,:,np.newaxis]), axis=2)
        
    # ===== Batal (14.16042) ====
    if '14.16042' in glaciers:
        # Data: debris thickness (m) and melt rate (mm w.e. d-1)
        mb_df = pd.read_csv(debris_prms.main_directory + '/../hd_obs/datasets/14.16042_patel2016-melt.csv')
        measured_hd_list = [mb_df.hd_m.values]
        measured_melt_list = [mb_df.melt_mmwed.values]
        melt_fp = debris_prms.output_fp + 'ostrem_curves/exp4/'
        melt_fn = '3225N-7750E-debris_melt_curve.nc'
        yearfracs_list = [[2014 + 213/365, 2014 + 288/365]]
    
        for n in np.arange(0,len(measured_hd_list)):
            assert len(measured_hd_list[n]) == len(measured_melt_list[n]), 'length of hd differs from melt'
        
        hd_wbnds_array_list = hd_melt_uncertainty(measured_hd_list, yearfracs_list, melt_fp, melt_fn, z_value=z_value)
        
        for hd_wbnds_array in hd_wbnds_array_list:
            if hd_wbnds_array_all is None:
                hd_wbnds_array_all = hd_wbnds_array[:,:,np.newaxis]
            else:
                hd_wbnds_array_all = np.concatenate((hd_wbnds_array_all, hd_wbnds_array[:,:,np.newaxis]), axis=2)
        
    # ===== Khumbu (15.03733) ====
    if '15.03733' in glaciers:
        # Data: debris thickness (m) and melt rate (mm w.e. d-1)
        mb_df = pd.read_csv(debris_prms.main_directory + '/../hd_obs/datasets/15.03733_kayastha2000-melt.csv')
        measured_hd_list = [mb_df.hd_m.values]
        measured_melt_list = [mb_df.melt_mmwed.values]
        melt_fp = debris_prms.output_fp + 'ostrem_curves/exp4/'
        melt_fn = '2800N-8700E-debris_melt_curve.nc'
        yearfracs_list = [[2000 + 143/366, 2000 + 153/366]]
    
        for n in np.arange(0,len(measured_hd_list)):
            assert len(measured_hd_list[n]) == len(measured_melt_list[n]), 'length of hd differs from melt'
        
        hd_wbnds_array_list = hd_melt_uncertainty(measured_hd_list, yearfracs_list, melt_fp, melt_fn, z_value=z_value)
        
        for hd_wbnds_array in hd_wbnds_array_list:
            if hd_wbnds_array_all is None:
                hd_wbnds_array_all = hd_wbnds_array[:,:,np.newaxis]
            else:
                hd_wbnds_array_all = np.concatenate((hd_wbnds_array_all, hd_wbnds_array[:,:,np.newaxis]), axis=2)

    # ===== Imja-Lhotse Shar (15.03743) ====
    if '15.03743' in glaciers:
        # Data: debris thickness (m) and melt rate (mm w.e. d-1)
        mb_df = pd.read_csv(debris_prms.main_directory + '/../hd_obs/datasets/15.03743_rounce2015-melt.csv')
        measured_hd_list = [mb_df.hd_m.values]
        measured_melt_list = [mb_df.melt_mmwed.values]
        melt_fp = debris_prms.output_fp + 'ostrem_curves/exp4/'
        melt_fn = '2800N-8700E-debris_melt_curve.nc'
        yearfracs_list = [[2014 + 138/365, 2014 + 315/365]]
    
        for n in np.arange(0,len(measured_hd_list)):
            assert len(measured_hd_list[n]) == len(measured_melt_list[n]), 'length of hd differs from melt'
        
        hd_wbnds_array_list = hd_melt_uncertainty(measured_hd_list, yearfracs_list, melt_fp, melt_fn, z_value=z_value)
        
        for hd_wbnds_array in hd_wbnds_array_list:
            if hd_wbnds_array_all is None:
                hd_wbnds_array_all = hd_wbnds_array[:,:,np.newaxis]
            else:
                hd_wbnds_array_all = np.concatenate((hd_wbnds_array_all, hd_wbnds_array[:,:,np.newaxis]), axis=2)
        
    # ===== Lirung (15.04045) ====
    if '15.04045' in glaciers:
        # Data: debris thickness (m) and melt rate (mm w.e. d-1)
        mb_df1 = pd.read_csv(debris_prms.main_directory + '/../hd_obs/datasets/15.04045_chand2015_fall-melt.csv')
        mb_df2 = pd.read_csv(debris_prms.main_directory + '/../hd_obs/datasets/15.04045_chand2015_winter-melt.csv')
        mb_df3 = pd.read_csv(debris_prms.main_directory + '/../hd_obs/datasets/15.04045_chand2015_spring-melt.csv')
        measured_hd_list = [mb_df1.hd_m.values, mb_df2.hd_m.values, mb_df3.hd_m.values]
        measured_melt_list = [mb_df1.melt_mmwed.values, mb_df2.melt_mmwed.values, mb_df3.melt_mmwed.values]
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
        
        hd_wbnds_array_list = hd_melt_uncertainty(measured_hd_list, yearfracs_list, melt_fp, melt_fn, z_value=z_value)
        
        for hd_wbnds_array in hd_wbnds_array_list:
            if hd_wbnds_array_all is None:
                hd_wbnds_array_all = hd_wbnds_array[:,:,np.newaxis]
            else:
                hd_wbnds_array_all = np.concatenate((hd_wbnds_array_all, hd_wbnds_array[:,:,np.newaxis]), axis=2)
        
    # ===== Hailuogou (15.07886) ====
    if '15.07886' in glaciers:
        # Data: debris thickness (m) and melt rate (mm w.e. d-1)
        measured_hd_list = [np.array([2, 2, 2, 2, 2, 2, 3, 3, 3, 3, 4, 5, 5, 6, 7, 7, 10, 10, 11, 13]) / 100]
        measured_melt_list = [np.array([65.2, 55.4, 52.8, 51.6, 47.0, 53.4, 44.4, 50.3, 58, 48.9, 58.4, 54.4, 44.8, 
                                        52.6, 43.7, 52.5, 38.5, 36.5, 34.2, 28.4])]
        melt_fp = debris_prms.output_fp + 'ostrem_curves/exp4/'
        melt_fn = '2950N-10200E-debris_melt_curve.nc'
        yearfracs_list = [[2008 + 184/366, 2008 + 274/366]]
    
        for n in np.arange(0,len(measured_hd_list)):
            assert len(measured_hd_list[n]) == len(measured_melt_list[n]), 'length of hd differs from melt'
        
        hd_wbnds_array_list = hd_melt_uncertainty(measured_hd_list, yearfracs_list, melt_fp, melt_fn, z_value=z_value)
        
        for hd_wbnds_array in hd_wbnds_array_list:
            if hd_wbnds_array_all is None:
                hd_wbnds_array_all = hd_wbnds_array[:,:,np.newaxis]
            else:
                hd_wbnds_array_all = np.concatenate((hd_wbnds_array_all, hd_wbnds_array[:,:,np.newaxis]), axis=2)
        
    # ===== 24K (15.11758) ====
    if '15.11758' in glaciers:
        # Data: debris thickness (m) and melt rate (mm w.e. d-1)
        mb_df = pd.read_csv(debris_prms.main_directory + '/../hd_obs/datasets/15.11758_yang2017-melt.csv')
        mb_df2 = pd.read_csv(debris_prms.main_directory + '/../hd_obs/datasets/15.11758_wei2010-melt.csv')
        measured_hd_list = [mb_df.hd_m.values, mb_df2.hd_m.values]
        measured_melt_list = [mb_df.melt_mmwed.values, mb_df2.melt_mmwed.values]
        melt_fp = debris_prms.output_fp + 'ostrem_curves/exp4/'
        melt_fn = '2975N-9575E-debris_melt_curve.nc'
        yearfracs_list = [[2016 + 153/366, 2016 + 274/366], [2008 + 201/366, 2008 + 248/366]]
    
        for n in np.arange(0,len(measured_hd_list)):
            assert len(measured_hd_list[n]) == len(measured_melt_list[n]), 'length of hd differs from melt'
        
        hd_wbnds_array_list = hd_melt_uncertainty(measured_hd_list, yearfracs_list, melt_fp, melt_fn, z_value=z_value)
        
        for hd_wbnds_array in hd_wbnds_array_list:
            if hd_wbnds_array_all is None:
                hd_wbnds_array_all = hd_wbnds_array[:,:,np.newaxis]
            else:
                hd_wbnds_array_all = np.concatenate((hd_wbnds_array_all, hd_wbnds_array[:,:,np.newaxis]), axis=2)
        
    # ===== Fox (18.02375) ====
    if '18.02375' in glaciers:
        # Data: debris thickness (m) and melt rate (mm w.e. d-1)
        mb_df = pd.read_csv(debris_prms.main_directory + '/../hd_obs/datasets/18.02375_brook2012-melt.csv')
        measured_hd_list = [mb_df.hd_m.values]
        measured_melt_list = [mb_df.melt_mmwed.values]
        melt_fp = debris_prms.output_fp + 'ostrem_curves/exp4/'
        melt_fn = '4350S-17025E-debris_melt_curve.nc'
        yearfracs_list = [[2007 + 327/365, 2007 + 337/365]]
    
        for n in np.arange(0,len(measured_hd_list)):
            assert len(measured_hd_list[n]) == len(measured_melt_list[n]), 'length of hd differs from melt'
        
        hd_wbnds_array_list = hd_melt_uncertainty(measured_hd_list, yearfracs_list, melt_fp, melt_fn, z_value=z_value)
        
        for hd_wbnds_array in hd_wbnds_array_list:
            if hd_wbnds_array_all is None:
                hd_wbnds_array_all = hd_wbnds_array[:,:,np.newaxis]
            else:
                hd_wbnds_array_all = np.concatenate((hd_wbnds_array_all, hd_wbnds_array[:,:,np.newaxis]), axis=2)
        
    # ===== Franz Josef (18.02397) ====
    if '18.02397' in glaciers:
        # Data: debris thickness (m) and melt rate (mm w.e. d-1)
        mb_df = pd.read_csv(debris_prms.main_directory + '/../hd_obs/datasets/18.02397_brook2013-melt.csv')
        measured_hd_list = [mb_df.hd_m.values]
        measured_melt_list = [mb_df.melt_mmwed.values]
        melt_fp = debris_prms.output_fp + 'ostrem_curves/exp4/'
        melt_fn = '4350S-17025E-debris_melt_curve.nc'
        yearfracs_list = [[2012 + 38/366, 2012 + 47/366]]
    
        for n in np.arange(0,len(measured_hd_list)):
            assert len(measured_hd_list[n]) == len(measured_melt_list[n]), 'length of hd differs from melt'
        
        hd_wbnds_array_list = hd_melt_uncertainty(measured_hd_list, yearfracs_list, melt_fp, melt_fn, z_value=z_value)
        
        for hd_wbnds_array in hd_wbnds_array_list:
            if hd_wbnds_array_all is None:
                hd_wbnds_array_all = hd_wbnds_array[:,:,np.newaxis]
            else:
                hd_wbnds_array_all = np.concatenate((hd_wbnds_array_all, hd_wbnds_array[:,:,np.newaxis]), axis=2)
                
    #%%
    # ===== EXPORT THE MEAN RELATIONSHIPS! =====
    hd_wbnds_array_mean = np.mean(hd_wbnds_array_all, axis=2)
    hd_wbnds_array_med = np.median(hd_wbnds_array_all, axis=2)
    hd_cns = ['hd_m', 'melt', 'melt_bndlow', 'melt_bndhigh', 'hd_bndlow_debris', 'hd_bndhigh_debris', 
              'hd_bndlow_elevchg', 'hd_bndhigh_elevchg', 'hd_bndlow_both', 'hd_bndhigh_both']
    hd_wbnds_df = pd.DataFrame(hd_wbnds_array_med, columns=hd_cns)
    
    # Export regional statistics
    hd_wbnds_df.to_csv(debris_prms.output_fp + 'hd_uncertainty_bnds.csv', index=False)
    
    # ===== Plot the relationship ======
    fontsize = 12
    hd_major = 1
    hd_minor = 0.1
    fig, ax = plt.subplots(1, 2, squeeze=False, sharex=False, sharey=False, gridspec_kw = {'wspace':0.35, 'hspace':0})
    for n in [0,1]:
        ax[0,n].plot(hd_wbnds_df['hd_m'], hd_wbnds_df['hd_bndlow_both'], zorder=3, color='k')
        ax[0,n].plot(hd_wbnds_df['hd_m'], hd_wbnds_df['hd_bndhigh_both'], zorder=3, color='k')
        ax[0,n].fill_between(hd_wbnds_df['hd_m'], hd_wbnds_df['hd_bndlow_both'], hd_wbnds_df['hd_bndhigh_both'],
                             color='k', linewidth=0, zorder=3, alpha=0.2)
        ax[0,n].plot([0,5],[0,5], zorder=1, color='k', linewidth=1)
    
        # X-label
        ax[0,n].set_xlabel('$h_{d}$ (m)', size=fontsize)
        ax[0,n].set_xlim(0,5)
        ax[0,n].xaxis.set_major_locator(plt.MultipleLocator(hd_major))
        ax[0,n].xaxis.set_minor_locator(plt.MultipleLocator(hd_minor))  
        # Y-label
        ax[0,n].set_ylabel('$h_{d}$ bounds (m)', size=fontsize)
        ax[0,n].set_ylim(0,5)
        ax[0,n].yaxis.set_major_locator(plt.MultipleLocator(hd_major))
        ax[0,n].yaxis.set_minor_locator(plt.MultipleLocator(hd_minor))
        # Tick parameters
        ax[0,n].yaxis.set_ticks_position('both')
        ax[0,n].tick_params(axis='both', which='major', labelsize=fontsize-2, direction='inout')
        ax[0,n].tick_params(axis='both', which='minor', labelsize=fontsize-4, direction='in') 
        
    # Inset ticks
    ax[0,1].set_xlim(0,0.5)
    ax[0,1].xaxis.set_major_locator(plt.MultipleLocator(0.1))
    ax[0,1].xaxis.set_minor_locator(plt.MultipleLocator(0.02))  
    ax[0,1].set_ylim(0,0.5)
    ax[0,1].yaxis.set_major_locator(plt.MultipleLocator(0.1))
    ax[0,1].yaxis.set_minor_locator(plt.MultipleLocator(0.02))  
    
    # Labels
    ax[0,0].text(0.1, 0.98, 'a', size=12, fontweight='bold',
                horizontalalignment='right', verticalalignment='top', transform=ax[0,0].transAxes)
    ax[0,1].text(0.1, 0.98, 'b', size=12, fontweight='bold',
                horizontalalignment='right', verticalalignment='top', transform=ax[0,1].transAxes)
    
    # Save plot
    fig.set_size_inches(6, 3)
    fig_fn = 'hd_uncertainty_bnds.png'
    fig.savefig(debris_prms.output_fp + fig_fn, bbox_inches='tight', dpi=300, transparent=True)
    
#%%
if (os.path.exists(debris_prms.output_fp + 'hd_uncertainty_bnds-IQR.csv') and 
    os.path.exists(debris_prms.output_fp + 'hd_uncertainty_bnds-90.csv')):
    
    # Export regional statistics
    hd_wbnds_df = pd.read_csv(debris_prms.output_fp + 'hd_uncertainty_bnds-IQR.csv')
    hd_wbnds_df_90 = pd.read_csv(debris_prms.output_fp + 'hd_uncertainty_bnds-90.csv')
    
    # ===== Plot the relationship ======
    fontsize = 12
    hd_major = 1
    hd_minor = 0.1
    lw_both = 1
    lw_component = 0.5
    fig, ax = plt.subplots(1, 3, squeeze=False, sharex=False, sharey=False, gridspec_kw = {'wspace':0.4, 'hspace':0})
    ax[0,2].axis('off')
    for n in [0,1]:
#        ax[0,n].plot(hd_wbnds_df['hd_m'], hd_wbnds_df['hd_bndlow_both'], zorder=3, color='k', linewidth=lw_both)
#        ax[0,n].plot(hd_wbnds_df['hd_m'], hd_wbnds_df['hd_bndhigh_both'], zorder=3, color='k', linewidth=lw_both)
        ax[0,n].plot(hd_wbnds_df['hd_m'], hd_wbnds_df['hd_bndlow_debris'], zorder=3, 
                     color='b', linewidth=lw_component, linestyle=':', label='IQR\ndebris\nproperties')
        ax[0,n].plot(hd_wbnds_df['hd_m'], hd_wbnds_df['hd_bndhigh_debris'], zorder=3, 
                     color='b', linewidth=lw_component, linestyle=':', label='')
        ax[0,n].plot(hd_wbnds_df['hd_m'], hd_wbnds_df['hd_bndlow_elevchg'], zorder=3, 
                     color='r', linewidth=lw_component, linestyle=':', label='IQR\nobserved\nmelt')
        ax[0,n].plot(hd_wbnds_df['hd_m'], hd_wbnds_df['hd_bndhigh_elevchg'], zorder=3, 
                     color='r', linewidth=lw_component, linestyle=':', label='')
        ax[0,n].fill_between(hd_wbnds_df['hd_m'], hd_wbnds_df['hd_bndlow_both'], hd_wbnds_df['hd_bndhigh_both'],
                             color='k', linewidth=0, zorder=3, alpha=0.2, label='IQR')
        ax[0,n].plot(hd_wbnds_df_90['hd_m'], hd_wbnds_df_90['hd_bndlow_both'], zorder=1, color='k', linestyle='--', 
                      linewidth=lw_both, label='90%')
        ax[0,n].plot(hd_wbnds_df_90['hd_m'], hd_wbnds_df_90['hd_bndhigh_both'], zorder=1, color='k', linestyle='--',
                      linewidth=lw_both, label='')
        ax[0,n].plot([0,5],[0,5], zorder=1, color='k', linewidth=1)
    
        # X-label
        ax[0,n].set_xlabel('$h_{d}$ (m)', size=fontsize)
        ax[0,n].set_xlim(0,5)
        ax[0,n].xaxis.set_major_locator(plt.MultipleLocator(hd_major))
        ax[0,n].xaxis.set_minor_locator(plt.MultipleLocator(hd_minor))  
        # Y-label
        ax[0,n].set_ylabel('$h_{d}$ bounds (m)', size=fontsize)
        ax[0,n].set_ylim(0,5)
        ax[0,n].yaxis.set_major_locator(plt.MultipleLocator(hd_major))
        ax[0,n].yaxis.set_minor_locator(plt.MultipleLocator(hd_minor))
        # Tick parameters
        ax[0,n].yaxis.set_ticks_position('both')
        ax[0,n].tick_params(axis='both', which='major', labelsize=fontsize-2, direction='inout')
        ax[0,n].tick_params(axis='both', which='minor', labelsize=fontsize-4, direction='in') 
        
        if n == 0:
            # Legend
            lgd = fig.legend(bbox_to_anchor=(0.69, 0.9), ncol=1, fontsize=10, frameon=False, handlelength=1, 
                             handletextpad=0.15, columnspacing=0.5, borderpad=0.25, labelspacing=1)

    # Second figure
    ax[0,1].set_xlim(0,0.5)
    ax[0,1].xaxis.set_major_locator(plt.MultipleLocator(0.1))
    ax[0,1].xaxis.set_minor_locator(plt.MultipleLocator(0.02))  
    ax[0,1].set_ylim(0,0.5)
    ax[0,1].yaxis.set_major_locator(plt.MultipleLocator(0.1))
    ax[0,1].yaxis.set_minor_locator(plt.MultipleLocator(0.02))  
    
    # Labels
    ax[0,0].text(0.1, 0.98, 'a', size=12, fontweight='bold',
                horizontalalignment='right', verticalalignment='top', transform=ax[0,0].transAxes)
    ax[0,1].text(0.1, 0.98, 'b', size=12, fontweight='bold',
                horizontalalignment='right', verticalalignment='top', transform=ax[0,1].transAxes)
    
    # Save plot
    fig.set_size_inches(8, 3)
    fig_fn = 'hd_uncertainty_bnds.png'
    fig.savefig(debris_prms.output_fp + fig_fn, bbox_inches='tight', dpi=300, transparent=True)
    

#%%
if option_melt_diagram_template:
    fig, ax = plt.subplots(1, 1, squeeze=False, sharex=False, sharey=False, gridspec_kw = {'wspace':0.35, 'hspace':0})
    xmin = 0.05
    xmax = 0.5
    ymin = 0
    ymax = 26.19 / (1 + 1.34 * 26.19 * xmin)
    x_values = np.arange(xmin+0.01,xmax-0.07,0.01)
    y_values = 26.19 / (1 + 1.34 * 26.19 * x_values)
    ax[0,0].plot(x_values, y_values, color='k', linewidth=1)
    
    x2_values = np.arange(xmin-0.01, xmax-0.07, 0.01)
    y2_values = 26.19 / (1 + 2.5 * 26.19 * x2_values)
    ax[0,0].plot(x2_values, y2_values, color='r', linewidth=1, linestyle='--')
    
    # Remove spines
    ax[0,0].spines['right'].set_visible(False)
    ax[0,0].spines['top'].set_visible(False)
    ax[0,0].spines['bottom'].set_visible(False)
    ax[0,0].spines['left'].set_visible(False)
    
    # Remove ticks
    plt.xticks([]) # labels 
    plt.yticks([])
    ax[0,0].xaxis.set_ticks_position('none') # tick markers
    ax[0,0].yaxis.set_ticks_position('none')
    
    # Automatic arrows
    dps = fig.dpi_scale_trans.inverted()
    bbox = ax[0,0].get_window_extent().transformed(dps)
    width, height = bbox.width, bbox.height
    
    # manual arrowhead width and length
    hw = 1./20.*(ymax-ymin) 
    hl = 1./20.*(xmax-xmin)
    lw = 1. # axis line width
    ohg = 0.3 # arrow overhang

    # compute matching arrowhead length and width
    yhw = hw/(ymax-ymin)*(xmax-0)* height/width 
    yhl = hl/(xmax-0)*(ymax-0)* width/height
    
    # draw x axis
    ax[0,0].arrow(0, 0, xmax-xmin, 0., fc='k', ec='k', lw = lw, 
             head_width=hw, head_length=hl, overhang = ohg, 
             length_includes_head= True, clip_on = False) 
    # draw x axis
    ax[0,0].arrow(0, 0, 0., ymax, fc='k', ec='k', lw = lw, 
             head_width=yhw, head_length=yhl, overhang = ohg, 
             length_includes_head= True, clip_on = False)
    
    ax[0,0].set_xlim(0,0.5)
    
    # X-label
#    ax[0,0].text(0.75, 0.9, 'a', size=12, fontweight='bold',
#                horizontalalignment='right', verticalalignment='top', transform=ax[0,0].transAxes)
#    ax[0,0].text(0.9, -0.1, '$h_{d}$', size=12,
#                horizontalalignment='right', verticalalignment='top', transform=ax[0,0].transAxes)
#    ax[0,0].text(-0.05, 0.97, 'Melt', size=12, rotation=90,
#                horizontalalignment='right', verticalalignment='top', transform=ax[0,0].transAxes)
    
    # Save plot
    fig.set_size_inches(3, 3)
    fig_fn = 'hd_melt_uncertainty_template.png'
    fig.savefig(debris_prms.output_fp + fig_fn, bbox_inches='tight', dpi=300, transparent=True)