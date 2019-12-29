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
debris_thicknesses = [1,2,4,6,8,10,12,15,20,30,50,100,200,300,400,500]
output_sample_fn = 'Rounce2015_hma_debris_XXcm_20191205.csv'
glacier_name = 'Ngozumpa Glacier (RGI60-01.15645)'
output_data_fp = input.main_directory + '/hma_data/output/exp3/'
output_fp = input.main_directory + '/hma_data/output/'
if os.path.exists(output_fp) == False:
    os.makedirs(output_fp)
    
#output_ds['YYYY-MM-DD'] = [str(met_data_all.loc[x,'YEAR']) + '-' + str(met_data_all.loc[x,'MONTH']).zfill(2) + '-' + 
#                              str(met_data_all.loc[x,'DAY']).zfill(2) for x in met_data_all.index.values]
#start_idx = list(met_data_all['YYYY-MM-DD']).index(start_date)
#end_idx = list(met_data_all['YYYY-MM-DD']).index(end_date) + 23

#%% ===== Debris Thickness vs. Surface Temperature FOR SPECIFIC DATES=====
ts_cns = ['debris_thickness']
for ts_date in input.ts_dates:
    ts_cns.append('ts_degC_' + ts_date)
debris_ts_df = pd.DataFrame(np.zeros((len(debris_thicknesses),len(ts_cns))), columns=ts_cns)
debris_ts_df['debris_thickness'] = np.array(debris_thicknesses) / 100
for ts_date in input.ts_dates:
    print(ts_date)
    ts_cn = 'ts_degC_' + ts_date 
    for ndebris, debris_thickness in enumerate(debris_thicknesses):
        output_fullfn = output_data_fp + output_sample_fn.replace('XX',str(debris_thickness))
        output_ds = pd.read_csv(output_fullfn)
        output_ds['Time_raw'] = output_ds['Time'].values
        output_ds['YYYY-MM-DD'] = [output_ds['Time_raw'].values[x].split(' ')[0] for x in output_ds.index.values]
        output_ds['Time'] = [np.datetime64(x) for x in output_ds['Time'].values]
        ts_idx = list(output_ds['YYYY-MM-DD']).index(ts_date) + int(input.ts_hr)
        ts_hr_adj = input.ts_hr%1   # amount to interpolate between values
        
        ts_degC = (output_ds['Ts [K]'].values[ts_idx] + ts_hr_adj * (output_ds['Ts [K]'].values[ts_idx+1] - 
                   output_ds['Ts [K]'].values[ts_idx]) - 273.15)
        
        debris_ts_df.loc[ndebris, ts_cn] = ts_degC
        print('  ', "{:.2f}".format(np.round(debris_thickness/100,2)) + ' m: ' + str(np.round(ts_degC,1)) + ' degC')

# PLOT CURVES
fig, ax = plt.subplots(1, 1, squeeze=False, sharex=False, sharey=False, gridspec_kw = {'wspace':0.4, 'hspace':0.15})
for ndate, ts_date in enumerate(input.ts_dates):
    ts_cn = 'ts_degC_' + ts_date 
    ax[0,0].plot(debris_ts_df['debris_thickness'], debris_ts_df[ts_cn], '-', zorder=1, label=ts_date)
# text
ax[0,0].text(0.5, 0.99, glacier_name, size=10, horizontalalignment='center', verticalalignment='top', 
             transform=ax[0,0].transAxes)
# X-label
ax[0,0].set_xlabel('Debris thickness(m)', size=12)
ax[0,0].set_xlim(0, 5)
#ax[0,0].set_xlim(0, debris_melt_df.debris_thickness.max())
ax[0,0].xaxis.set_tick_params(labelsize=12)
ax[0,0].xaxis.set_major_locator(plt.MultipleLocator(0.5))
ax[0,0].xaxis.set_minor_locator(plt.MultipleLocator(0.1))  
# Y-label
ax[0,0].set_ylabel('Ts (degC)', size=12)
#ax[0,0].set_ylim(0,(int(debris_ts_df.melt_mwea.values.max()/0.1)+3)*0.1)
ax[0,0].yaxis.set_major_locator(plt.MultipleLocator(10))
ax[0,0].yaxis.set_minor_locator(plt.MultipleLocator(2))
# Tick parameters
ax[0,0].yaxis.set_ticks_position('both')
ax[0,0].tick_params(axis='both', which='major', labelsize=12, direction='inout')
ax[0,0].tick_params(axis='both', which='minor', labelsize=10, direction='in')     
# Legend
ax[0,0].legend(fontsize=8, labelspacing=0.25, handlelength=1, handletextpad=0.25, borderpad=0, frameon=False)          
# Save plot
fig.set_size_inches(4, 4)
figure_fn = 'debris_Ts_curve_specificdates.png'
fig.savefig(output_fp + figure_fn, bbox_inches='tight', dpi=300)  
        
#%% ===== Debris Thickness vs. Surface Temperature - STATS OVER MELT SEASON=====
month_start = 8
month_end = 9
ts_cns = ['debris_thickness', 'ts_degC_mean', 'ts_degC_std', 'ts_degC_med', 'ts_degC_nmad']
debris_ts_df_stats = pd.DataFrame(np.zeros((len(debris_thicknesses),len(ts_cns))), columns=ts_cns)
debris_ts_df_stats['debris_thickness'] = np.array(debris_thicknesses) / 100

for ndebris, debris_thickness in enumerate(debris_thicknesses):
    output_fullfn = output_data_fp + output_sample_fn.replace('XX',str(debris_thickness))
    output_ds = pd.read_csv(output_fullfn)
    output_ds['Time'] = [np.datetime64(x) for x in output_ds['Time'].values]

    # Calculate the Ts at the acquisition time for all dates during melt season
    ts_11 = output_ds.loc[int(input.ts_hr)::24,['Time', 'Melt [mwe]', 'Ts [K]', 'Rn [W m2]', 'LE [W m2]', 
                                                'H [W m2]', 'P [W m2]', 'Qc [W m2]', 'd_snow [m]']]
    ts_12 = output_ds.loc[(int(input.ts_hr)+1)::24,['Time', 'Melt [mwe]', 'Ts [K]', 'Rn [W m2]', 'LE [W m2]', 
                                                    'H [W m2]', 'P [W m2]', 'Qc [W m2]', 'd_snow [m]']]
    ts_acq = ts_11.copy()
    ts_acq.loc[:,:] = (ts_11.values + ts_hr_adj * (ts_12.values - ts_11.values))
    ts_acq['Month'] = [ts_acq.loc[x,'Time'].month for x in ts_acq.index.values]
    ts_acq_month = ts_acq['Month'].values        
    ts_acq_subset = ts_acq.loc[(ts_acq['Month'] >= month_start) & (ts_acq['Month'] <= month_end),:]
    ts_acq_subset.reset_index(inplace=True, drop=True)
    
    debris_ts_df_stats.loc[ndebris,'ts_degC_mean'] = ts_acq_subset['Ts [K]'].mean() - 273.15
    debris_ts_df_stats.loc[ndebris,'ts_degC_std'] = ts_acq_subset['Ts [K]'].std()
    debris_ts_df_stats.loc[ndebris,'ts_degC_med'] = ts_acq_subset['Ts [K]'].median() - 273.15
    debris_ts_df_stats.loc[ndebris,'ts_degC_nmad'] = (
            1.483 * np.median(np.absolute(ts_acq_subset['Ts [K]'] - ts_acq_subset['Ts [K]'].median()).values))

    
#%%
# PLOT CURVE
# fit curve
    
## Kraaijenbrink et al. (2017)
##def ts_fromdebris_func(h, hmax, T95):
##    """ estimate surface temperature from debris thickness (h is debris thickness, a and k are coefficients) 
##    Kraaijenbrink et al. (2017) equation"""
##    return 0 + np.log(h*100) / np.log(hmax*100) * (T95 - 0)
##def debris_fromts_func(ts, hmax, T95):
##    """ estimate debris thickness from Ts (b is melt, a and k are coefficients) """
##    return np.exp((ts - 0) / (T95 - 0) * np.log(hmax)) / 100
#def ts_fromdebris_func(h, T95):
#    """ estimate surface temperature from debris thickness (h is debris thickness, a and k are coefficients) 
#    Kraaijenbrink et al. (2017) equation"""
#    return 0 + np.log(h*100) / np.log(1*100) * (T95 - 0)
#def debris_fromts_func(ts, T95):
#    """ estimate debris thickness from Ts (b is melt, a and k are coefficients) """
#    return np.exp((ts - 0) / (T95 - 0) * np.log(1)) / 100
#fit_idx = debris_ts_df_stats.index.values
#func_coeff, pcov = curve_fit(ts_fromdebris_func, 
#                             debris_ts_df_stats.debris_thickness.values[fit_idx], debris_ts_df_stats.ts_degC_med.values[fit_idx],
##                             p0=[30,0.3]
#                             )
#fit_melt = ts_fromdebris_func(debris_ts_df_stats.debris_thickness.values, func_coeff[0])
#debris_4curve = np.arange(0.01,5.01,0.01)
#ts_4curve = ts_fromdebris_func(debris_4curve, func_coeff[0])
#ts_4curve_setmax = ts_fromdebris_func(debris_4curve, func_coeff[0])
#label_4tscurve = 'Kra2017'
#figure_fn = 'debris_Ts_curve_Kra2017.png'

#def ts_fromdebris_func(h, a, k):
#    """ estimate surface temperature from debris thickness (h is debris thickness, a and k are coefficients) 
#    Michaelis-Menten Equation from single-substrate reactions"""
#    return a - a*np.exp(-k*h)
def ts_fromdebris_func(h, a, k):
    """ estimate surface temperature from debris thickness (h is debris thickness, a and k are coefficients) 
    hyperbolic fit"""
    return a * h / (k + h)
def debris_fromts_func(ts, a, k):
    """ estimate debris thickness from Ts (b is melt, a and k are coefficients) """
    return k * ts / (a - ts)

#fit_idx = debris_ts_df_stats.index.values
fit_idx = [2,4,7,8,9,10,11,12,13,14]
print('\nFIT ONLY TO DEBRIS THICKNESSES:', debris_ts_df_stats.debris_thickness.values[fit_idx], 
      '\nTO ENSURE CAPTURES CURVE \n(plot with equal spacing to avoid this issue, but computationally expensive)')
func_coeff, pcov = curve_fit(ts_fromdebris_func, 
                             debris_ts_df_stats.debris_thickness.values[fit_idx], debris_ts_df_stats.ts_degC_med.values[fit_idx],
#                             p0=[30,0.3]
                             )
fit_melt = ts_fromdebris_func(debris_ts_df_stats.debris_thickness.values, func_coeff[0], func_coeff[1])
debris_4curve = np.arange(0.02,5.01,0.01)
ts_4curve = ts_fromdebris_func(debris_4curve, func_coeff[0], func_coeff[1])
ts_max = 25
ts_4curve_setmax = ts_fromdebris_func(debris_4curve, ts_max, func_coeff[1])
label_4tscurve = 'set max ('+ str(ts_max) + ' degC)'

if len(np.where(ts_4curve < 0)[0]) > 0:
    print('\nNEGATIVE FIT!\n')

figure_fn = 'debris_Ts_curve_avg_hyperbolic.png'

# plot curve
fig, ax = plt.subplots(1, 1, squeeze=False, sharex=False, sharey=False, gridspec_kw = {'wspace':0.4, 'hspace':0.15})
#ax[0,0].plot(debris_ts_df_stats['debris_thickness'], debris_ts_df_stats['ts_degC_mean'], '-', zorder=1, label='mean')
ax[0,0].plot(debris_ts_df_stats['debris_thickness'], debris_ts_df_stats['ts_degC_med'], 
             '-', color='k', zorder=1, label='median')
ax[0,0].fill_between(debris_ts_df_stats['debris_thickness'], 
                     debris_ts_df_stats['ts_degC_med'] - debris_ts_df_stats['ts_degC_nmad'], 
                     debris_ts_df_stats['ts_degC_med'] + debris_ts_df_stats['ts_degC_nmad'],
                      facecolor='k', alpha=0.2, zorder=1)
ax[0,0].plot(debris_4curve, ts_4curve, 
             color='b', linewidth=0.5, linestyle='--', zorder=2, label='best fit')
ax[0,0].plot(debris_4curve, ts_4curve_setmax, 
             color='r', linewidth=0.5, linestyle='--', zorder=2, label=label_4tscurve)
# text
ax[0,0].text(0.5, 0.99, glacier_name, size=10, horizontalalignment='center', verticalalignment='top', 
             transform=ax[0,0].transAxes)
#eqn_text = r'$T_{s} = \frac{T_{s,max} h}{k + h}$'
#coeff1_text = r'$T_{s,max} = ' + str(np.round(func_coeff[0],2)) + '$' 
#coeff2_text = r'$k = ' + str(np.round(func_coeff[1],2)) + '$' 
## coeff$\frac{b_{0}}{1 + 2kb_{0}h}$'
#ax[0,0].text(0.05, 0.9, eqn_text, size=12, horizontalalignment='left', verticalalignment='top', 
#             transform=ax[0,0].transAxes)
#ax[0,0].text(0.06, 0.79, 'where', size=10, horizontalalignment='left', verticalalignment='top', 
#             transform=ax[0,0].transAxes)
#ax[0,0].text(0.09, 0.74, coeff1_text, size=10, horizontalalignment='left', verticalalignment='top', 
#             transform=ax[0,0].transAxes)
#ax[0,0].text(0.09, 0.67, coeff2_text, size=10, horizontalalignment='left', verticalalignment='top', 
#             transform=ax[0,0].transAxes)
# X-label
ax[0,0].set_xlabel('Debris thickness(m)', size=12)
ax[0,0].set_xlim(0, 2)
#ax[0,0].set_xlim(0, debris_melt_df.debris_thickness.max())
ax[0,0].xaxis.set_tick_params(labelsize=12)
#ax[0,0].xaxis.set_major_locator(plt.MultipleLocator(0.5))
#ax[0,0].xaxis.set_minor_locator(plt.MultipleLocator(0.1))  
# Y-label
ax[0,0].set_ylabel('Ts (degC)', size=12)
#ax[0,0].set_ylim(0,(int(debris_ts_df.melt_mwea.values.max()/0.1)+3)*0.1)
#ax[0,0].yaxis.set_major_locator(plt.MultipleLocator(10))
#ax[0,0].yaxis.set_minor_locator(plt.MultipleLocator(2))
# Tick parameters
ax[0,0].yaxis.set_ticks_position('both')
ax[0,0].tick_params(axis='both', which='major', labelsize=12, direction='inout')
ax[0,0].tick_params(axis='both', which='minor', labelsize=10, direction='in')     
# Legend
ax[0,0].legend(loc=(0.05, 0.05), 
               fontsize=8, labelspacing=0.25, handlelength=1, handletextpad=0.25, borderpad=0, 
               frameon=False)          
# Save plot
fig.set_size_inches(4, 4)
fig.savefig(output_fp + figure_fn, bbox_inches='tight', dpi=300)  


#%% ===== THE "HILL" EQUATION FITS NGOZUMPA NICELY ======
def ts_fromdebris_func(h, a, b, c):
    """ estimate surface temperature from debris thickness (h is debris thickness, a and k are coefficients) 
    nonlinear"""
    return a * h**c / (b**c + h**c)

#fit_idx = debris_ts_df_stats.index.values
#fit_idx = [1,4,7,8,9,10,12]
fit_idx = [2,4,7,8,9,10,11,12,13,14]
print('\nFIT ONLY TO DEBRIS THICKNESSES:', debris_ts_df_stats.debris_thickness.values[fit_idx], 
      '\nTO ENSURE CAPTURES CURVE \n(plot with equal spacing to avoid this issue, but computationally expensive)')
func_coeff, pcov = curve_fit(ts_fromdebris_func, 
                             debris_ts_df_stats.debris_thickness.values[fit_idx], debris_ts_df_stats.ts_degC_med.values[fit_idx],
                             p0=[25,1,1]
                             )
fit_melt = ts_fromdebris_func(debris_ts_df_stats.debris_thickness.values, func_coeff[0], func_coeff[1], func_coeff[2])
debris_4curve = np.arange(0.02,5.01,0.01)
ts_4curve = ts_fromdebris_func(debris_4curve, func_coeff[0], func_coeff[1], func_coeff[2])
ts_max = 18
ts_4curve_setmax = ts_fromdebris_func(debris_4curve, ts_max, func_coeff[1], func_coeff[2])
label_4tscurve = 'set max ('+ str(ts_max) + ' degC)'

if len(np.where(ts_4curve < 0)[0]) > 0:
    print('\nNEGATIVE FIT!\n')

figure_fn = 'debris_Ts_curve_avg_Hill.png'

# plot curve
fig, ax = plt.subplots(1, 1, squeeze=False, sharex=False, sharey=False, gridspec_kw = {'wspace':0.4, 'hspace':0.15})
#ax[0,0].plot(debris_ts_df_stats['debris_thickness'], debris_ts_df_stats['ts_degC_mean'], '-', zorder=1, label='mean')
ax[0,0].plot(debris_ts_df_stats['debris_thickness'], debris_ts_df_stats['ts_degC_med'], 
             '-', color='k', zorder=1, label='median')
ax[0,0].fill_between(debris_ts_df_stats['debris_thickness'], 
                     debris_ts_df_stats['ts_degC_med'] - debris_ts_df_stats['ts_degC_nmad'], 
                     debris_ts_df_stats['ts_degC_med'] + debris_ts_df_stats['ts_degC_nmad'],
                      facecolor='k', alpha=0.2, zorder=1)
ax[0,0].plot(debris_4curve, ts_4curve, 
             color='b', linewidth=0.5, linestyle='--', zorder=2, label='best fit')
ax[0,0].plot(debris_4curve, ts_4curve_setmax, 
             color='r', linewidth=0.5, linestyle='--', zorder=2, label=label_4tscurve)
# text
ax[0,0].text(0.5, 0.99, glacier_name, size=10, horizontalalignment='center', verticalalignment='top', 
             transform=ax[0,0].transAxes)
#eqn_text = r'$T_{s} = \frac{T_{s,max} h^}{k + h}$'
#coeff1_text = r'$T_{s,max} = ' + str(np.round(func_coeff[0],2)) + '$' 
#coeff2_text = r'$k = ' + str(np.round(func_coeff[1],2)) + '$' 
# coeff$\frac{b_{0}}{1 + 2kb_{0}h}$'
#ax[0,0].text(0.05, 0.9, eqn_text, size=12, horizontalalignment='left', verticalalignment='top', 
#             transform=ax[0,0].transAxes)
#ax[0,0].text(0.06, 0.79, 'where', size=10, horizontalalignment='left', verticalalignment='top', 
#             transform=ax[0,0].transAxes)
#ax[0,0].text(0.09, 0.74, coeff1_text, size=10, horizontalalignment='left', verticalalignment='top', 
#             transform=ax[0,0].transAxes)
#ax[0,0].text(0.09, 0.67, coeff2_text, size=10, horizontalalignment='left', verticalalignment='top', 
#             transform=ax[0,0].transAxes)
# X-label
ax[0,0].set_xlabel('Debris thickness(m)', size=12)
ax[0,0].set_xlim(0, 2)
#ax[0,0].set_xlim(0, debris_melt_df.debris_thickness.max())
ax[0,0].xaxis.set_tick_params(labelsize=12)
#ax[0,0].xaxis.set_major_locator(plt.MultipleLocator(0.5))
#ax[0,0].xaxis.set_minor_locator(plt.MultipleLocator(0.1))  
# Y-label
ax[0,0].set_ylabel('Ts (degC)', size=12)
#ax[0,0].set_ylim(0,(int(debris_ts_df.melt_mwea.values.max()/0.1)+3)*0.1)
#ax[0,0].yaxis.set_major_locator(plt.MultipleLocator(10))
#ax[0,0].yaxis.set_minor_locator(plt.MultipleLocator(2))
# Tick parameters
ax[0,0].yaxis.set_ticks_position('both')
ax[0,0].tick_params(axis='both', which='major', labelsize=12, direction='inout')
ax[0,0].tick_params(axis='both', which='minor', labelsize=10, direction='in')     
# Legend
ax[0,0].legend(loc=(0.05, 0.05), 
               fontsize=8, labelspacing=0.25, handlelength=1, handletextpad=0.25, borderpad=0, 
               frameon=False)          
# Save plot
fig.set_size_inches(4, 4)
fig.savefig(output_fp + figure_fn, bbox_inches='tight', dpi=300)  
print('hill coeffs:', func_coeff)
#%%


##%% ===== THIS EQUATION CAN'T BE INVERTED EASILY BUT FITS NGOZUMPA NICELY ======
#def ts_fromdebris_func(h, a, b, c):
#    """ estimate surface temperature from debris thickness (h is debris thickness, a and k are coefficients) 
#    nonlinear"""
#    return 1 / (a + b * h**c)
#
#fit_idx = debris_ts_df_stats.index.values
##fit_idx = [2,4,7,8,9,10,11,12,13,14]
#print('\nFIT ONLY TO DEBRIS THICKNESSES:', debris_ts_df_stats.debris_thickness.values[fit_idx], 
#      '\nTO ENSURE CAPTURES CURVE \n(plot with equal spacing to avoid this issue, but computationally expensive)')
#func_coeff, pcov = curve_fit(ts_fromdebris_func, 
#                             debris_ts_df_stats.debris_thickness.values[fit_idx], debris_ts_df_stats.ts_degC_med.values[fit_idx]
#                             )
#fit_melt = ts_fromdebris_func(debris_ts_df_stats.debris_thickness.values, func_coeff[0], func_coeff[1], func_coeff[2])
#debris_4curve = np.arange(0.02,5.01,0.01)
#ts_4curve = ts_fromdebris_func(debris_4curve, func_coeff[0], func_coeff[1], func_coeff[2])
#ts_max = 25
#ts_4curve_setmax = ts_fromdebris_func(debris_4curve, ts_max, func_coeff[1], func_coeff[2])
#label_4tscurve = 'set max ('+ str(ts_max) + ' degC)'
#
#if len(np.where(ts_4curve < 0)[0]) > 0:
#    print('\nNEGATIVE FIT!\n')
#
#figure_fn = 'debris_Ts_curve_avg.png'
#
## plot curve
#fig, ax = plt.subplots(1, 1, squeeze=False, sharex=False, sharey=False, gridspec_kw = {'wspace':0.4, 'hspace':0.15})
##ax[0,0].plot(debris_ts_df_stats['debris_thickness'], debris_ts_df_stats['ts_degC_mean'], '-', zorder=1, label='mean')
#ax[0,0].plot(debris_ts_df_stats['debris_thickness'], debris_ts_df_stats['ts_degC_med'], 
#             '-', color='k', zorder=1, label='median')
#ax[0,0].fill_between(debris_ts_df_stats['debris_thickness'], 
#                     debris_ts_df_stats['ts_degC_med'] - debris_ts_df_stats['ts_degC_nmad'], 
#                     debris_ts_df_stats['ts_degC_med'] + debris_ts_df_stats['ts_degC_nmad'],
#                      facecolor='k', alpha=0.2, zorder=1)
#ax[0,0].plot(debris_4curve, ts_4curve, 
#             color='b', linewidth=0.5, linestyle='--', zorder=2, label='best fit')
#ax[0,0].plot(debris_4curve, ts_4curve_setmax, 
#             color='r', linewidth=0.5, linestyle='--', zorder=2, label=label_4tscurve)
## text
#ax[0,0].text(0.5, 0.99, glacier_name, size=10, horizontalalignment='center', verticalalignment='top', 
#             transform=ax[0,0].transAxes)
##eqn_text = r'$T_{s} = \frac{T_{s,max} h}{k + h}$'
##coeff1_text = r'$T_{s,max} = ' + str(np.round(func_coeff[0],2)) + '$' 
##coeff2_text = r'$k = ' + str(np.round(func_coeff[1],2)) + '$' 
### coeff$\frac{b_{0}}{1 + 2kb_{0}h}$'
##ax[0,0].text(0.05, 0.9, eqn_text, size=12, horizontalalignment='left', verticalalignment='top', 
##             transform=ax[0,0].transAxes)
##ax[0,0].text(0.06, 0.79, 'where', size=10, horizontalalignment='left', verticalalignment='top', 
##             transform=ax[0,0].transAxes)
##ax[0,0].text(0.09, 0.74, coeff1_text, size=10, horizontalalignment='left', verticalalignment='top', 
##             transform=ax[0,0].transAxes)
##ax[0,0].text(0.09, 0.67, coeff2_text, size=10, horizontalalignment='left', verticalalignment='top', 
##             transform=ax[0,0].transAxes)
## X-label
#ax[0,0].set_xlabel('Debris thickness(m)', size=12)
#ax[0,0].set_xlim(0, 2)
##ax[0,0].set_xlim(0, debris_melt_df.debris_thickness.max())
#ax[0,0].xaxis.set_tick_params(labelsize=12)
##ax[0,0].xaxis.set_major_locator(plt.MultipleLocator(0.5))
##ax[0,0].xaxis.set_minor_locator(plt.MultipleLocator(0.1))  
## Y-label
#ax[0,0].set_ylabel('Ts (degC)', size=12)
##ax[0,0].set_ylim(0,(int(debris_ts_df.melt_mwea.values.max()/0.1)+3)*0.1)
##ax[0,0].yaxis.set_major_locator(plt.MultipleLocator(10))
##ax[0,0].yaxis.set_minor_locator(plt.MultipleLocator(2))
## Tick parameters
#ax[0,0].yaxis.set_ticks_position('both')
#ax[0,0].tick_params(axis='both', which='major', labelsize=12, direction='inout')
#ax[0,0].tick_params(axis='both', which='minor', labelsize=10, direction='in')     
## Legend
#ax[0,0].legend(loc=(0.05, 0.05), 
#               fontsize=8, labelspacing=0.25, handlelength=1, handletextpad=0.25, borderpad=0, 
#               frameon=False)          
## Save plot
#fig.set_size_inches(4, 4)
#fig.savefig(output_fp + figure_fn, bbox_inches='tight', dpi=300)  



#%% ===== DERIVE DEBRIS THICKNESS FROM THERMAL DATA =====
print('TAKE EQUATION AND MOVE TO QGIS FOR MAP - SHOULD MAKE WORKFLOW BETTER AND HAVE IT ALL IN PYTHON')