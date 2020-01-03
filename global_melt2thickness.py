# -*- coding: utf-8 -*-
"""
Created on Sat Aug 11 10:17:13 2018

@author: David
"""

# Built-in libraries
import collections
import os

# External libraries
#import rasterio
#import gdal
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit
import xarray as xr

# Local libraries
import globaldebris_input as input
from meltcurves import debris_frommelt_func

#%% ===== SCRIPT SPECIFIC INFORMATION =====
stats_idx = 0

#%% ===== FUNCTIONS ====
def selectglaciersrgitable(glac_no=None,
                           rgi_regionsO1=None,
                           rgi_regionsO2=None,
                           rgi_glac_number=None,
#                            rgi_fp=input.rgi_fp,
                           rgi_fp = '/Users/davidrounce/Documents/Dave_Rounce/HiMAT/RGI/rgi60/00_rgi60_attribs/',
                           rgi_cols_drop=['GLIMSId','BgnDate','EndDate','Status','Connect','Linkages','Name'],
                           rgi_O1Id_colname='glacno',
                           rgi_glacno_float_colname='RGIId_float',
                           indexname='GlacNo'):
    """
    Select all glaciers to be used in the model run according to the regions and glacier numbers defined by the RGI
    glacier inventory. This function returns the rgi table associated with all of these glaciers.

    glac_no : list of strings
        list of strings of RGI glacier numbers (e.g., ['1.00001', '13.00001'])
    rgi_regionsO1 : list of integers
        list of integers of RGI order 1 regions (e.g., [1, 13])
    rgi_regionsO2 : list of integers or 'all'
        list of integers of RGI order 2 regions or simply 'all' for all the order 2 regions
    rgi_glac_number : list of strings
        list of RGI glacier numbers without the region (e.g., ['00001', '00002'])

    Output: Pandas DataFrame of the glacier statistics for each glacier in the model run
    (rows = GlacNo, columns = glacier statistics)
    """
    if glac_no is not None:
        glac_no_byregion = {}
        rgi_regionsO1 = [int(i.split('.')[0]) for i in glac_no]
        rgi_regionsO1 = list(set(rgi_regionsO1))
        for region in rgi_regionsO1:
            glac_no_byregion[region] = []
        for i in glac_no:
            region = i.split('.')[0]
            glac_no_only = i.split('.')[1]
            glac_no_byregion[int(region)].append(glac_no_only)

        for region in rgi_regionsO1:
            glac_no_byregion[region] = sorted(glac_no_byregion[region])

    # Create an empty dataframe
    rgi_regionsO1 = sorted(rgi_regionsO1)
    glacier_table = pd.DataFrame()
    for region in rgi_regionsO1:

        if glac_no is not None:
            rgi_glac_number = glac_no_byregion[region]

#        if len(rgi_glac_number) < 50:

        for i in os.listdir(rgi_fp):
            if i.startswith(str(region).zfill(2)) and i.endswith('.csv'):
                rgi_fn = i
        try:
            csv_regionO1 = pd.read_csv(rgi_fp + rgi_fn)
        except:
            csv_regionO1 = pd.read_csv(rgi_fp + rgi_fn, encoding='latin1')
        
        # Populate glacer_table with the glaciers of interest
        if rgi_regionsO2 == 'all' and rgi_glac_number == 'all':
            print("All glaciers within region(s) %s are included in this model run." % (region))
            if glacier_table.empty:
                glacier_table = csv_regionO1
            else:
                glacier_table = pd.concat([glacier_table, csv_regionO1], axis=0)
        elif rgi_regionsO2 != 'all' and rgi_glac_number == 'all':
            print("All glaciers within subregion(s) %s in region %s are included in this model run." %
                  (rgi_regionsO2, region))
            for regionO2 in rgi_regionsO2:
                if glacier_table.empty:
                    glacier_table = csv_regionO1.loc[csv_regionO1['O2Region'] == regionO2]
                else:
                    glacier_table = (pd.concat([glacier_table, csv_regionO1.loc[csv_regionO1['O2Region'] ==
                                                                                regionO2]], axis=0))
        else:
            if len(rgi_glac_number) < 20:
                print("%s glaciers in region %s are included in this model run: %s" % (len(rgi_glac_number), region,
                                                                                       rgi_glac_number))
            else:
                print("%s glaciers in region %s are included in this model run: %s and more" %
                      (len(rgi_glac_number), region, rgi_glac_number[0:50]))
                
            rgiid_subset = ['RGI60-' + str(region).zfill(2) + '.' + x for x in rgi_glac_number] 
            rgiid_all = list(csv_regionO1.RGIId.values)
            rgi_idx = [rgiid_all.index(x) for x in rgiid_subset]
            if glacier_table.empty:
                glacier_table = csv_regionO1.loc[rgi_idx]
            else:
                glacier_table = (pd.concat([glacier_table, csv_regionO1.loc[rgi_idx]],
                                           axis=0))
                    
    glacier_table = glacier_table.copy()
    # reset the index so that it is in sequential order (0, 1, 2, etc.)
    glacier_table.reset_index(inplace=True)
    # change old index to 'O1Index' to be easier to recall what it is
    glacier_table.rename(columns={'index': 'O1Index'}, inplace=True)
    # Record the reference date
    glacier_table['RefDate'] = glacier_table['BgnDate']
    # if there is an end date, then roughly average the year
    enddate_idx = glacier_table.loc[(glacier_table['EndDate'] > 0), 'EndDate'].index.values
    glacier_table.loc[enddate_idx,'RefDate'] = (
            np.mean((glacier_table.loc[enddate_idx,['BgnDate', 'EndDate']].values / 10**4).astype(int),
                    axis=1).astype(int) * 10**4 + 9999)
    # drop columns of data that is not being used
    glacier_table.drop(rgi_cols_drop, axis=1, inplace=True)
    # add column with the O1 glacier numbers
    glacier_table[rgi_O1Id_colname] = (
            glacier_table['RGIId'].str.split('.').apply(pd.Series).loc[:,1].astype(int))
    glacier_table['rgino_str'] = [x.split('-')[1] for x in glacier_table.RGIId.values]
    glacier_table[rgi_glacno_float_colname] = (np.array([np.str.split(glacier_table['RGIId'][x],'-')[1]
                                                    for x in range(glacier_table.shape[0])]).astype(float))
    # set index name
    glacier_table.index.name = indexname

    print("This study is focusing on %s glaciers in region %s" % (glacier_table.shape[0], rgi_regionsO1))

    return glacier_table


#%% ==== TO-DO LIST =====
print('\n%TO-DO LIST:\n-----------')
print('1. estimate debris with and without flux divergence')
print('  a. quality control emergence velocities ')
print('  b. total uncertainty should be based on root sum of squares from MB and emvel')
print('2. identify glaciers with relatively stable velocities where confident in debris thickness')
print('  a. consider a mixture of velocities and % of the terminus?')
print('3. assess max debris thickness for each region and extrapolate to regions without data')
print('  a. assess any regional trends based on Bolch et al. regions')
print('  b. quality control: ex. 13.00147 has positive mass balance - surging glacier?\n\n')

#%%
# ===== SELECT GLACIERS WITH DATA =====
rgiid_list = []
rgiid_fn_list = []
for i in os.listdir(input.mb_binned_fp):
    if i.endswith('mb_bins.csv'):
        rgiid_list.append(i[0:8])
        rgiid_fn_list.append(i)
        
rgiid_list = sorted(rgiid_list)
rgiid_fn_list = sorted(rgiid_fn_list)

print(len(rgiid_list))

main_glac_rgi = selectglaciersrgitable(rgiid_list)
main_glac_rgi['bin_fn'] = rgiid_fn_list

# ==== DETERMINE NEAREST LAT/LON =====
ds_elev = xr.open_dataset(input.metdata_fp + input.metdata_elev_fn)
#  argmin() finds the minimum distance between the glacier lat/lon and the GCM pixel
lat_nearidx = (np.abs(main_glac_rgi['CenLat'].values[:,np.newaxis] - 
                      ds_elev['latitude'][:].values).argmin(axis=1))
lon_nearidx = (np.abs(main_glac_rgi['CenLon'].values[:,np.newaxis] - 
                      ds_elev['longitude'][:].values).argmin(axis=1))

latlon_nearidx = list(zip(lat_nearidx, lon_nearidx))
latlon_nearidx_unique = sorted(list(set(latlon_nearidx)))

main_glac_rgi['lat_nearest'] = ds_elev['latitude'][lat_nearidx].values 
main_glac_rgi['lon_nearest'] = ds_elev['longitude'][lon_nearidx].values 
main_glac_rgi['latlon_nearidx'] = latlon_nearidx
latlon_unique_dict = dict(zip(latlon_nearidx_unique,np.arange(0,len(latlon_nearidx_unique))))
latlon_unique_dict_reversed = dict(zip(np.arange(0,len(latlon_nearidx_unique)),latlon_nearidx_unique))
main_glac_rgi['latlon_unique_no'] = main_glac_rgi['latlon_nearidx'].map(latlon_unique_dict)

#%%
# ===== ESTIMATE DEBRIS THICKNESS =====
latlon_nodata = []
glac_nodata = []
#for nglac, glac_idx in enumerate(main_glac_rgi.index.values):
#for nglac, glac_idx in enumerate([main_glac_rgi.index.values[6940]]):
for nglac, glac_idx in enumerate([main_glac_rgi.index.values[4957]]):
    glac_str = main_glac_rgi.loc[glac_idx,'rgino_str']
    print(nglac, glac_idx, glac_str)
    
    # ===== Load Mass Balance Data =====
    mb_df_fn = main_glac_rgi.loc[glac_idx,'bin_fn']
    mb_df = pd.read_csv(input.mb_binned_fp + mb_df_fn)
    mb_df.loc[:,:] = mb_df.values.astype(np.float64)
    
    # ===== Load Ostrem data =====
    latlon_str = (str(int(main_glac_rgi.loc[glac_idx,'lat_nearest']*100)) + 'N-' + 
                  str(int(main_glac_rgi.loc[glac_idx,'lon_nearest']*100)) + 'E') 
    ostrem_fn = input.output_ostrem_fn_sample.replace('XXXX', latlon_str)
    
    
    if os.path.exists(input.output_ostrem_fp + ostrem_fn) and ' perc_debris' in mb_df.columns:
        ds_ostrem = xr.open_dataset(input.output_ostrem_fp + ostrem_fn)
        
        # Nearest elevation
        # Debris elevation (weighted mean)
        mb_area_total = main_glac_rgi.loc[glac_idx,'Area']
        mb_df['bin_area'] = mb_df[' z1_bin_area_perc'].values / 100 * mb_area_total
        mb_df['bin_area_debris'] = mb_df[' perc_debris'].values / 100 * mb_df['bin_area'].values
        debris_elev_mean = ((mb_df['# bin_center_elev_m'] * mb_df['bin_area_debris']).sum() 
                            /  mb_df['bin_area_debris'].sum()) 
        ds_elevidx = np.abs(debris_elev_mean - ds_ostrem.elev.values).argmin()
        
        # Melt function coefficients
        func_coeff = [float(ds_ostrem['b0'][stats_idx,ds_elevidx].values), 
                      float(ds_ostrem['k'][stats_idx,ds_elevidx].values)]
        
        # Estimate debris thickness
        mb_df['mb_bin_mean_mwea_1stdlow'] = mb_df[' mb_bin_mean_mwea'] - mb_df[' mb_bin_std_mwea']
        mb_df['mb_bin_mean_mwea_1stdhigh'] = mb_df[' mb_bin_mean_mwea'] + mb_df[' mb_bin_std_mwea']
        mb_df['hd'] = debris_frommelt_func(-1*mb_df[' mb_bin_mean_mwea'].values, 
                                                         func_coeff[0], func_coeff[1])
        mb_df['hd_1stdlow'] = debris_frommelt_func(-1*mb_df['mb_bin_mean_mwea_1stdlow'].values, 
                                                                  func_coeff[0], func_coeff[1])
        mb_df['hd_1stdhigh'] = debris_frommelt_func(-1*mb_df['mb_bin_mean_mwea_1stdhigh'].values, 
                                                                  func_coeff[0], func_coeff[1])

        # ==== Add Emergence Velocity =====
        # Load file
        emvel_fn = input.outdir_emvel_fp + glac_str + input.output_emvel_csv_ending
        emvel_df = pd.read_csv(emvel_fn)
        emvel_df.loc[:,:] = np.nan_to_num(emvel_df.values.astype(np.float64), 0)
        
        mb_df[' vm_med_itslive'] = emvel_df[' vm_med']
        mb_df[' vm_mad_itslive'] = emvel_df[' vm_mad']
        mb_df = pd.concat([mb_df, emvel_df[[' emvel_mean', ' emvel_std']]], axis=1)
        # Uncertainty with flux divergence from Farinotti et al. (2019)
        mb_df[' emvel_high'] = mb_df[' emvel_mean'] * 1.6
        mb_df[' emvel_low'] = mb_df[' emvel_mean'] * 0.8
    
        # Mass balance with emergence velocity and uncertainties
        mb_df['mb_wem'] = mb_df[' mb_bin_mean_mwea'] - mb_df[' emvel_mean']
        # Higher emergence --> more melt
        mb_df['mb_wemthick'] = mb_df[' mb_bin_mean_mwea'] - mb_df[' emvel_high'] - mb_df[' mb_bin_std_mwea']
        # Lower emergence --> less melt
        mb_df['mb_wemthin'] = mb_df[' mb_bin_mean_mwea'] - mb_df[' emvel_low'] + mb_df[' mb_bin_std_mwea']
        
        # Debris thickness with emergence velocity and uncertainties
        mb_df['hd_wem'] = debris_frommelt_func(-1*mb_df['mb_wem'].values, func_coeff[0], func_coeff[1])
        mb_df['hd_wemthick'] = debris_frommelt_func(-1*mb_df['mb_wemthick'].values, func_coeff[0], func_coeff[1])
        mb_df['hd_wemthin'] = debris_frommelt_func(-1*mb_df['mb_wemthin'].values, func_coeff[0], func_coeff[1])
    
        # Export dataset
        if os.path.exists(input.mb_binned_fp_wdebris) == False:
            os.makedirs(input.mb_binned_fp_wdebris)
        
        mb_df_fn_wdebris = mb_df_fn.replace('.csv','_wdebris.csv')
        mb_df.to_csv(input.mb_binned_fp_wdebris + mb_df_fn_wdebris, index=False)
        ds_ostrem.close()
    elif os.path.exists(input.output_ostrem_fp + ostrem_fn) == False:
        latlon_nodata.append(latlon_str)
        glac_nodata.append(main_glac_rgi.loc[glac_idx,'rgino_str'])

latlon_nodata = sorted(list(set(latlon_nodata)))
if len(latlon_nodata) > 0:
    print('\nMissing Ostrem data:\n', latlon_nodata,'\n')
    print('Glaciers without Ostrem data:\n', glac_nodata,'\n')
    print('\nTHESE GLACIERS DONT HAVE CONSIDERABLE DEBRIS AND WERE THEREFORE EXCLUDED\n')
        
#%% ===== Input Data =====
### PLOT MASS BALANCE AND DEBRIS THICKNESS VS. ELEVATION WITHOUT EMERGENCE VELOCITIES
### Mass balance
##fig, ax = plt.subplots(1, 2, squeeze=False, sharex=False, sharey=True, gridspec_kw = {'wspace':0.1, 'hspace':0.15})
##ax[0,0].plot(mb_df[' mb_bin_mean_mwea'], mb_df['# bin_center_elev_m'], color='k', linewidth=1, zorder=2)
##ax[0,0].fill_betweenx(mb_df['# bin_center_elev_m'], 
##                      mb_df['mb_bin_mean_mwea_1stdlow'], mb_df['mb_bin_mean_mwea_1stdhigh'], 
##                      facecolor='k', alpha=0.2, zorder=1)
### text
##fig.text(0.5, 0.96, glacier_name, size=10, horizontalalignment='center', verticalalignment='top', 
##             transform=ax[0,0].transAxes)
### X-label
##ax[0,0].set_xlabel('Mass Balance\n(mwea)', size=12)
##ax[0,0].set_xlim(-3.25, 0.25)
###ax[0,0].set_xlim(0, debris_melt_df.debris_thickness.max())
##ax[0,0].xaxis.set_tick_params(labelsize=12)
##ax[0,0].xaxis.set_major_locator(plt.MultipleLocator(1))
##ax[0,0].xaxis.set_minor_locator(plt.MultipleLocator(0.1))  
### Y-label
##ylower = mb_df['# bin_center_elev_m'].min()
##yupper = ylower + 1000
##ax[0,0].set_ylabel('Elevation (m asl)', size=12)
##ax[0,0].set_ylim(ylower,yupper)
##ax[0,0].yaxis.set_major_locator(plt.MultipleLocator(100))
##ax[0,0].yaxis.set_minor_locator(plt.MultipleLocator(50))
### Tick parameters
##ax[0,0].yaxis.set_ticks_position('both')
##ax[0,0].tick_params(axis='both', which='major', labelsize=12, direction='inout')
##ax[0,0].tick_params(axis='both', which='minor', labelsize=10, direction='in')       
##
### Debris thickness
##ax[0,1].plot(mb_df['debris_thickness'], mb_df['# bin_center_elev_m'], color='k', linewidth=1, zorder=2)
###ax[0,0].plot(mb_df['MassBal_25'], mb_df['E'], color='k', linewidth=0.5, zorder=2)
###ax[0,0].plot(mb_df['MassBal_75'], mb_df['E'], color='k', linewidth=0.5, zorder=2)
##ax[0,1].fill_betweenx(mb_df['# bin_center_elev_m'], 
##                      mb_df['debris_thickness_1stdlow'], mb_df['debris_thickness_1stdhigh'], 
##                      facecolor='k', alpha=0.2, zorder=1)
### X-label
##ax[0,1].set_xlabel('Debris thickness\n(m)', size=12)
##ax[0,1].set_xlim(0, 1.5) 
##ax[0,1].xaxis.set_major_locator(plt.MultipleLocator(0.5))
##ax[0,1].xaxis.set_minor_locator(plt.MultipleLocator(0.1))  
### Y-label
##ax[0,1].set_ylim(ylower,yupper)
### Tick parameters
##ax[0,1].yaxis.set_ticks_position('both')
##ax[0,1].tick_params(axis='both', which='major', labelsize=12, direction='inout')
##ax[0,1].tick_params(axis='both', which='minor', labelsize=10, direction='in')       
##        
### Save plot
##fig.set_size_inches(4, 4)
##figure_fn = 'elev_mb_debris.png'
##fig.savefig(output_fp + figure_fn, bbox_inches='tight', dpi=300) 
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

###%% ===== PLOT MASS BALANCE AND DEBRIS THICKNESS VS. ELEVATION WITH EMERGENCE VELOCITIES ======
### Mass balance
##fig, ax = plt.subplots(1, 2, squeeze=False, sharex=False, sharey=True, gridspec_kw = {'wspace':0.1, 'hspace':0.15})
##ax[0,0].plot(mb_df['mb_wem'], mb_df['# bin_center_elev_m'], color='k', linewidth=1, zorder=2)
###ax[0,0].fill_betweenx(mb_df['# bin_center_elev_m'], 
###                      mb_df['mb_bin_mean_mwea_1stdlow'], mb_df['mb_bin_mean_mwea_1stdhigh'], 
###                      facecolor='k', alpha=0.2, zorder=1)
### text
##fig.text(0.5, 0.96, glacier_name, size=10, horizontalalignment='center', verticalalignment='top', 
##             transform=ax[0,0].transAxes)
### X-label
##ax[0,0].set_xlabel('Mass Balance\n(mwea)', size=12)
##ax[0,0].set_xlim(-3.25, 0.25)
###ax[0,0].set_xlim(0, debris_melt_df.debris_thickness.max())
##ax[0,0].xaxis.set_tick_params(labelsize=12)
##ax[0,0].xaxis.set_major_locator(plt.MultipleLocator(1))
##ax[0,0].xaxis.set_minor_locator(plt.MultipleLocator(0.1))  
### Y-label
##ylower = mb_df['# bin_center_elev_m'].min()
##yupper = ylower + 1000
##ax[0,0].set_ylabel('Elevation (m asl)', size=12)
##ax[0,0].set_ylim(ylower,yupper)
##ax[0,0].yaxis.set_major_locator(plt.MultipleLocator(100))
##ax[0,0].yaxis.set_minor_locator(plt.MultipleLocator(50))
### Tick parameters
##ax[0,0].yaxis.set_ticks_position('both')
##ax[0,0].tick_params(axis='both', which='major', labelsize=12, direction='inout')
##ax[0,0].tick_params(axis='both', which='minor', labelsize=10, direction='in')       
##
### Debris thickness
##ax[0,1].plot(mb_df['debris_thickness_wem'], mb_df['# bin_center_elev_m'], color='k', linewidth=1, zorder=2)
###ax[0,0].plot(mb_df['MassBal_25'], mb_df['E'], color='k', linewidth=0.5, zorder=2)
###ax[0,0].plot(mb_df['MassBal_75'], mb_df['E'], color='k', linewidth=0.5, zorder=2)
###ax[0,1].fill_betweenx(mb_df['# bin_center_elev_m'], 
###                      mb_df['debris_thickness_1stdlow'], mb_df['debris_thickness_1stdhigh'], 
###                      facecolor='k', alpha=0.2, zorder=1)
### X-label
##ax[0,1].set_xlabel('Debris thickness\n(m)', size=12)
##ax[0,1].set_xlim(0, 1.5) 
##ax[0,1].xaxis.set_major_locator(plt.MultipleLocator(0.5))
##ax[0,1].xaxis.set_minor_locator(plt.MultipleLocator(0.1))  
### Y-label
##ax[0,1].set_ylim(ylower,yupper)
### Tick parameters
##ax[0,1].yaxis.set_ticks_position('both')
##ax[0,1].tick_params(axis='both', which='major', labelsize=12, direction='inout')
##ax[0,1].tick_params(axis='both', which='minor', labelsize=10, direction='in')       
##        
### Save plot
##fig.set_size_inches(4, 4)
##figure_fn = 'elev_mb_debris_wem.png'
##fig.savefig(output_fp + figure_fn, bbox_inches='tight', dpi=300) 