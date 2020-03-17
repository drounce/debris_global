# -*- coding: utf-8 -*-
"""
Created on Sat Aug 11 10:17:13 2018

@author: David
"""


import os
#import rasterio
#import gdal
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

# Local libraries
import hma_input as input

#%% ===== Input Data =====
glac_str = '15.03473'
#debris_thicknesses = [1,2,4,6,8,10,12,15,20,30,50,100,200,300,400,500]
output_sample_fn = 'Rounce2015_hma_debris_XXcm_20191206.csv'
#larsen_data_fullfn = '/Users/davidrounce/Documents/Dave_Rounce/HiMAT/DEMs/larsen/data/Kennicott.2000.2013.output.txt'
#larsen_ice_density = 850 # kg m-3
mb_fullfn = ('/Users/davidrounce/Documents/Dave_Rounce/HiMAT/DEMs/Shean_2019_0213/' + 
             'mb_combined_20190213_nmad_bins/' + glac_str + '_mb_bins.csv')
glacier_name = 'Ngozumpa Glacier (RGI60-01.15645)'
output_data_fp = input.main_directory + '/hma_data/output/change/'
#output_data_fp = input.main_directory + '/output/exp3/'
output_fp = input.main_directory + '/hma_data/output/'
#anderson_debris_fullfn = input.main_directory + '/kennicott_data/anderson_debristhickness.csv'
#anderson_debris = pd.read_csv(anderson_debris_fullfn)
if os.path.exists(output_fp) == False:
    os.makedirs(output_fp)
    

#%% ===== Debris Thickness vs. Surface Lowering =====
#debris_thicknesses = [2,15,30,50,100,200,300,400,500]
debris_thicknesses = np.arange(0.0,5,0.1)
debris_thicknesses[0] = 0.02
debris_thicknesses = [int(x*100) for x in debris_thicknesses]
debris_melt_df = pd.DataFrame(np.zeros((len(debris_thicknesses),2)), columns=['debris_thickness', 'melt_mwea'])
#for ndebris, debris_thickness in enumerate(debris_thicknesses):
for i in os.listdir(output_data_fp):
    if i.endswith('.csv') and 'hma' in i:
#        output_fullfn = output_data_fp + output_sample_fn.replace('XX',str(debris_thickness))
        output_fullfn = output_data_fp + i
        print(output_fullfn)
        output_ds = pd.read_csv(output_fullfn)
        output_ds.to_csv(output_data_fp + i.replace('hma','Ngozumpa'))
#    output_ds['Time'] = [np.datetime64(x) for x in output_ds['Time'].values]
#    
#    melt_mwea = (output_ds['Melt [mwe]'].sum() / 
#                 ((output_ds['Time'].values[-1] - output_ds['Time'].values[0]).astype('timedelta64[D]').astype(int) 
#                  / 365.25))
#    debris_melt_df.loc[ndebris] = debris_thickness / 100, melt_mwea
#    print(str(np.round(debris_thickness/100,2)) + ' m: ' + str(np.round(melt_mwea,2)) + ' mwea')

## PLOT CURVE
## fit curve
#def melt_fromdebris_func(h, a, k):
#    """ estimate melt from debris thickness (h is debris thickness, a and k are coefficients) """
#    return a / (1 + 2 * k * a * h)
#def debris_frommelt_func(b, a, k):
#    """ estimate debris thickness from melt (b is melt, a and k are coefficients) """
#    return (a - b) / (2*k*a*b)
#func_coeff, pcov = curve_fit(melt_fromdebris_func, 
#                             debris_melt_df.debris_thickness.values, debris_melt_df.melt_mwea.values)
#fit_melt = melt_fromdebris_func(debris_melt_df.debris_thickness.values, func_coeff[0], func_coeff[1])
#debris_4curve = np.arange(0.01,5.01,0.01)
#melt_4curve = melt_fromdebris_func(debris_4curve, func_coeff[0], func_coeff[1])
#
## plot curve
#fig, ax = plt.subplots(1, 1, squeeze=False, sharex=False, sharey=False, gridspec_kw = {'wspace':0.4, 'hspace':0.15})
#ax[0,0].plot(debris_melt_df['debris_thickness'], debris_melt_df['melt_mwea'], 'o', 
#             color='k', markersize=3, zorder=1)
#ax[0,0].plot(debris_4curve, melt_4curve, 
#             color='b', linewidth=0.5, linestyle='--', zorder=2, label='plot1')
## text
#ax[0,0].text(0.5, 0.99, glacier_name, size=10, horizontalalignment='center', verticalalignment='top', 
#             transform=ax[0,0].transAxes)
#eqn_text = r'$b = \frac{b_{0}}{1 + 2kb_{0}h}$'
#coeff1_text = r'$b_{0} = ' + str(np.round(func_coeff[0],2)) + '$' 
#coeff2_text = r'$k = ' + str(np.round(func_coeff[1],2)) + '$' 
## coeff$\frac{b_{0}}{1 + 2kb_{0}h}$'
#ax[0,0].text(0.9, 0.85, eqn_text, size=12, horizontalalignment='right', verticalalignment='top', 
#             transform=ax[0,0].transAxes)
#ax[0,0].text(0.615, 0.73, 'where', size=10, horizontalalignment='left', verticalalignment='top', 
#             transform=ax[0,0].transAxes)
#ax[0,0].text(0.66, 0.67, coeff1_text, size=10, horizontalalignment='left', verticalalignment='top', 
#             transform=ax[0,0].transAxes)
#ax[0,0].text(0.66, 0.6, coeff2_text, size=10, horizontalalignment='left', verticalalignment='top', 
#             transform=ax[0,0].transAxes)
## X-label
#ax[0,0].set_xlabel('Debris thickness(m)', size=12)
#ax[0,0].set_xlim(0, 2.1)
##ax[0,0].set_xlim(0, debris_melt_df.debris_thickness.max())
#ax[0,0].xaxis.set_tick_params(labelsize=12)
#ax[0,0].xaxis.set_major_locator(plt.MultipleLocator(0.5))
#ax[0,0].xaxis.set_minor_locator(plt.MultipleLocator(0.1))  
## Y-label
#ax[0,0].set_ylabel('Melt (mwea)', size=12)
#ax[0,0].set_ylim(0,(int(debris_melt_df.melt_mwea.values.max()/0.1)+3)*0.1)
#ax[0,0].yaxis.set_major_locator(plt.MultipleLocator(1))
#ax[0,0].yaxis.set_minor_locator(plt.MultipleLocator(0.1))
## Tick parameters
#ax[0,0].yaxis.set_ticks_position('both')
#ax[0,0].tick_params(axis='both', which='major', labelsize=12, direction='inout')
#ax[0,0].tick_params(axis='both', which='minor', labelsize=10, direction='in')               
## Save plot
#fig.set_size_inches(4, 4)
#figure_fn = 'debris_melt_curve.png'
#fig.savefig(output_fp + figure_fn, bbox_inches='tight', dpi=300)  
#
##%% ===== DERIVE DEBRIS THICKNESS FROM MASS BALANCE DATA =====
#mb_df = pd.read_csv(mb_fullfn)
#mb_df.loc[:,:] = mb_df.values.astype(np.float64)
#mb_df['debris_thickness'] = debris_frommelt_func(-1*mb_df[' mb_bin_mean_mwea'].values, func_coeff[0], func_coeff[1])
#mb_df['mb_bin_mean_mwea_1stdlow'] = mb_df[' mb_bin_mean_mwea'] - mb_df[' mb_bin_std_mwea']
#mb_df['mb_bin_mean_mwea_1stdhigh'] = mb_df[' mb_bin_mean_mwea'] + mb_df[' mb_bin_std_mwea']
#mb_df['debris_thickness_1stdlow'] = debris_frommelt_func(-1*mb_df['mb_bin_mean_mwea_1stdlow'].values, 
#                                                          func_coeff[0], func_coeff[1])
#mb_df['debris_thickness_1stdhigh'] = debris_frommelt_func(-1*mb_df['mb_bin_mean_mwea_1stdhigh'].values, 
#                                                          func_coeff[0], func_coeff[1])
#
## PLOT MASS BALANCE AND DEBRIS THICKNESS VS. ELEVATION WITHOUT EMERGENCE VELOCITIES
## Mass balance
#fig, ax = plt.subplots(1, 2, squeeze=False, sharex=False, sharey=True, gridspec_kw = {'wspace':0.1, 'hspace':0.15})
#ax[0,0].plot(mb_df[' mb_bin_mean_mwea'], mb_df['# bin_center_elev_m'], color='k', linewidth=1, zorder=2)
#ax[0,0].fill_betweenx(mb_df['# bin_center_elev_m'], 
#                      mb_df['mb_bin_mean_mwea_1stdlow'], mb_df['mb_bin_mean_mwea_1stdhigh'], 
#                      facecolor='k', alpha=0.2, zorder=1)
## text
#fig.text(0.5, 0.96, glacier_name, size=10, horizontalalignment='center', verticalalignment='top', 
#             transform=ax[0,0].transAxes)
## X-label
#ax[0,0].set_xlabel('Mass Balance\n(mwea)', size=12)
#ax[0,0].set_xlim(-3.25, 0.25)
##ax[0,0].set_xlim(0, debris_melt_df.debris_thickness.max())
#ax[0,0].xaxis.set_tick_params(labelsize=12)
#ax[0,0].xaxis.set_major_locator(plt.MultipleLocator(1))
#ax[0,0].xaxis.set_minor_locator(plt.MultipleLocator(0.1))  
## Y-label
#ylower = mb_df['# bin_center_elev_m'].min()
#yupper = ylower + 1000
#ax[0,0].set_ylabel('Elevation (m asl)', size=12)
#ax[0,0].set_ylim(ylower,yupper)
#ax[0,0].yaxis.set_major_locator(plt.MultipleLocator(100))
#ax[0,0].yaxis.set_minor_locator(plt.MultipleLocator(50))
## Tick parameters
#ax[0,0].yaxis.set_ticks_position('both')
#ax[0,0].tick_params(axis='both', which='major', labelsize=12, direction='inout')
#ax[0,0].tick_params(axis='both', which='minor', labelsize=10, direction='in')       
#
## Debris thickness
#ax[0,1].plot(mb_df['debris_thickness'], mb_df['# bin_center_elev_m'], color='k', linewidth=1, zorder=2)
##ax[0,0].plot(mb_df['MassBal_25'], mb_df['E'], color='k', linewidth=0.5, zorder=2)
##ax[0,0].plot(mb_df['MassBal_75'], mb_df['E'], color='k', linewidth=0.5, zorder=2)
#ax[0,1].fill_betweenx(mb_df['# bin_center_elev_m'], 
#                      mb_df['debris_thickness_1stdlow'], mb_df['debris_thickness_1stdhigh'], 
#                      facecolor='k', alpha=0.2, zorder=1)
## X-label
#ax[0,1].set_xlabel('Debris thickness\n(m)', size=12)
#ax[0,1].set_xlim(0, 1.5) 
#ax[0,1].xaxis.set_major_locator(plt.MultipleLocator(0.5))
#ax[0,1].xaxis.set_minor_locator(plt.MultipleLocator(0.1))  
## Y-label
#ax[0,1].set_ylim(ylower,yupper)
## Tick parameters
#ax[0,1].yaxis.set_ticks_position('both')
#ax[0,1].tick_params(axis='both', which='major', labelsize=12, direction='inout')
#ax[0,1].tick_params(axis='both', which='minor', labelsize=10, direction='in')       
#        
## Save plot
#fig.set_size_inches(4, 4)
#figure_fn = 'elev_mb_debris.png'
#fig.savefig(output_fp + figure_fn, bbox_inches='tight', dpi=300) 
#
##%% ===== EMERGENCE VELOCITIES ==========================================
## CURRENTLY PROCESSED IN IPYTHON NOTEBOOK USING SHEAN GMBTOOLS-TYPE FILE 'emergence_velocity.ipynb'
#emergence_fullfn = ('/Users/davidrounce/Documents/Dave_Rounce/HiMAT/DEMs/Shean_2019_0213/' + 
#                    'mb_combined_20190213_nmad_bins/csv/' + glac_str + '_mb_bins_wemvel.csv')
#emergence_df = pd.read_csv(emergence_fullfn)
#emergence_df['area_cumsum'] = np.cumsum(emergence_df['z1_bin_area_valid_km2'])
#
#binsize_mb = mb_df.loc[1,'# bin_center_elev_m'] - mb_df.loc[0,'# bin_center_elev_m']
#emvel_binsize = emergence_df['# bin_center_elev_m'].values[1] - emergence_df['# bin_center_elev_m'].values[0]
#emergence_shift_idx = np.where(emergence_df.area_cumsum.values < mb_df.loc[0,' z1_bin_area_valid_km2'])[0][-1]
#mb_offset = (mb_df.loc[0, '# bin_center_elev_m'] + binsize_mb/2 - emvel_binsize / 2 - 
#             emergence_df.loc[emergence_shift_idx, '# bin_center_elev_m'])
#                      
#emergence_df['# bin_center_elev_m'] = emergence_df['# bin_center_elev_m'] + mb_offset
#emergence_df['E_low'] = emergence_df['# bin_center_elev_m'] - emvel_binsize/2
#emergence_df['E_high'] = emergence_df['# bin_center_elev_m'] + emvel_binsize/2
#
## Get mean emergence velocity to coincide with elevation bins
#mb_df['E_low'] =  mb_df['# bin_center_elev_m'] - binsize_mb/2
#mb_df['E_high'] = mb_df['# bin_center_elev_m'] + binsize_mb/2
#
#mb_df['em_idx_low'] = np.nan
#mb_df['em_idx_high'] = np.nan
#for x in mb_df.index.values:
#    rows_low = np.where(mb_df.E_low.values[x] == emergence_df.E_low.values)[0]
#    if len(rows_low) > 0:
#        mb_df.loc[x,'em_idx_low'] = rows_low[0]
#    elif x == 0:
#        mb_df.loc[x,'em_idx_low'] = 0
#        
#    rows_high = np.where(mb_df.E_high.values[x] == emergence_df.E_high.values)[0]
#    if len(rows_high) > 0:
#        mb_df.loc[x,'em_idx_high'] = rows_high[0]
#    elif len(rows_high) == 0 and ~np.isnan(mb_df.loc[x,'em_idx_low']):
#        mb_df.loc[x,'em_idx_high'] = emergence_df.index.values[-1]
#
#emergence_df['emvel*area'] = emergence_df.emvel_mean * emergence_df.z1_bin_area_valid_km2
#emergence_df['emvel*area_1stdlow'] = (emergence_df.emvel_mean - emergence_df.emvel_std) * emergence_df.z1_bin_area_valid_km2
#emergence_df['emvel*area_1stdhigh'] = (emergence_df.emvel_mean + emergence_df.emvel_std) * emergence_df.z1_bin_area_valid_km2
#            
#mb_df['emvel_myr'] = np.nan
#for x in mb_df.index.values:
#    if ~np.isnan(mb_df.loc[x,'em_idx_low']):
#        mb_df.loc[x,'emvel_myr'] = (
#                emergence_df.loc[mb_df.loc[x,'em_idx_low']:mb_df.loc[x,'em_idx_high'], 'emvel*area'].sum() / 
#                emergence_df.loc[mb_df.loc[x,'em_idx_low']:mb_df.loc[x,'em_idx_high'], 'z1_bin_area_valid_km2'].sum())
##larsen_data['emvel_myr_1stdlow'] = (
##        [emvel_data.loc[larsen_data.loc[x,'em_idx_low']:larsen_data.loc[x,'em_idx_high'], 'emvel*area_1stdlow'].sum() / 
##         emvel_data.loc[larsen_data.loc[x,'em_idx_low']:larsen_data.loc[x,'em_idx_high'], 'z1_bin_area_valid_km2'].sum()
##         for x in larsen_data.index.values])
##larsen_data['emvel_myr_1stdhigh'] = (
##        [emvel_data.loc[larsen_data.loc[x,'em_idx_low']:larsen_data.loc[x,'em_idx_high'], 'emvel*area_1stdhigh'].sum() / 
##         emvel_data.loc[larsen_data.loc[x,'em_idx_low']:larsen_data.loc[x,'em_idx_high'], 'z1_bin_area_valid_km2'].sum()
##         for x in larsen_data.index.values])
#
#mb_df['mb_wem'] = mb_df[' mb_bin_mean_mwea'] - mb_df['emvel_myr']
###larsen_data['mb_wem_25'] = larsen_data['MassBal_25'] - larsen_data['emvel_myr_1stdhigh']
###larsen_data['mb_wem_75'] = larsen_data['MassBal_75'] - larsen_data['emvel_myr_1stdlow']
##larsen_data['mb_wem_25'] = larsen_data['MassBal_25'] - larsen_data['emvel_myr']
##larsen_data['mb_wem_75'] = larsen_data['MassBal_75'] - larsen_data['emvel_myr']
##larsen_data['debris_thickness_wem'] = debris_frommelt_func(-1*larsen_data['mb_wem'].values, func_coeff[0], func_coeff[1])
##larsen_data['debris_thickness_wem_25'] = debris_frommelt_func(-1*larsen_data['mb_wem_25'].values, 
##                                                          func_coeff[0], func_coeff[1])
##larsen_data['debris_thickness_wem_75'] = debris_frommelt_func(-1*larsen_data['mb_wem_75'].values, 
##                                                          func_coeff[0], func_coeff[1])
#
#mb_df['debris_thickness_wem'] = debris_frommelt_func(-1*mb_df['mb_wem'].values, func_coeff[0], func_coeff[1])
#
##%% ===== PLOT MASS BALANCE AND DEBRIS THICKNESS VS. ELEVATION WITH EMERGENCE VELOCITIES ======
## Mass balance
#fig, ax = plt.subplots(1, 2, squeeze=False, sharex=False, sharey=True, gridspec_kw = {'wspace':0.1, 'hspace':0.15})
#ax[0,0].plot(mb_df['mb_wem'], mb_df['# bin_center_elev_m'], color='k', linewidth=1, zorder=2)
##ax[0,0].fill_betweenx(mb_df['# bin_center_elev_m'], 
##                      mb_df['mb_bin_mean_mwea_1stdlow'], mb_df['mb_bin_mean_mwea_1stdhigh'], 
##                      facecolor='k', alpha=0.2, zorder=1)
## text
#fig.text(0.5, 0.96, glacier_name, size=10, horizontalalignment='center', verticalalignment='top', 
#             transform=ax[0,0].transAxes)
## X-label
#ax[0,0].set_xlabel('Mass Balance\n(mwea)', size=12)
#ax[0,0].set_xlim(-3.25, 0.25)
##ax[0,0].set_xlim(0, debris_melt_df.debris_thickness.max())
#ax[0,0].xaxis.set_tick_params(labelsize=12)
#ax[0,0].xaxis.set_major_locator(plt.MultipleLocator(1))
#ax[0,0].xaxis.set_minor_locator(plt.MultipleLocator(0.1))  
## Y-label
#ylower = mb_df['# bin_center_elev_m'].min()
#yupper = ylower + 1000
#ax[0,0].set_ylabel('Elevation (m asl)', size=12)
#ax[0,0].set_ylim(ylower,yupper)
#ax[0,0].yaxis.set_major_locator(plt.MultipleLocator(100))
#ax[0,0].yaxis.set_minor_locator(plt.MultipleLocator(50))
## Tick parameters
#ax[0,0].yaxis.set_ticks_position('both')
#ax[0,0].tick_params(axis='both', which='major', labelsize=12, direction='inout')
#ax[0,0].tick_params(axis='both', which='minor', labelsize=10, direction='in')       
#
## Debris thickness
#ax[0,1].plot(mb_df['debris_thickness_wem'], mb_df['# bin_center_elev_m'], color='k', linewidth=1, zorder=2)
##ax[0,0].plot(mb_df['MassBal_25'], mb_df['E'], color='k', linewidth=0.5, zorder=2)
##ax[0,0].plot(mb_df['MassBal_75'], mb_df['E'], color='k', linewidth=0.5, zorder=2)
##ax[0,1].fill_betweenx(mb_df['# bin_center_elev_m'], 
##                      mb_df['debris_thickness_1stdlow'], mb_df['debris_thickness_1stdhigh'], 
##                      facecolor='k', alpha=0.2, zorder=1)
## X-label
#ax[0,1].set_xlabel('Debris thickness\n(m)', size=12)
#ax[0,1].set_xlim(0, 1.5) 
#ax[0,1].xaxis.set_major_locator(plt.MultipleLocator(0.5))
#ax[0,1].xaxis.set_minor_locator(plt.MultipleLocator(0.1))  
## Y-label
#ax[0,1].set_ylim(ylower,yupper)
## Tick parameters
#ax[0,1].yaxis.set_ticks_position('both')
#ax[0,1].tick_params(axis='both', which='major', labelsize=12, direction='inout')
#ax[0,1].tick_params(axis='both', which='minor', labelsize=10, direction='in')       
#        
## Save plot
#fig.set_size_inches(4, 4)
#figure_fn = 'elev_mb_debris_wem.png'
#fig.savefig(output_fp + figure_fn, bbox_inches='tight', dpi=300) 