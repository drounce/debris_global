#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Model input for intercomparison experiment
"""
# Built-in libraries
import os
import pickle
# External libraries
import numpy as np
import pandas as pd

#%% ===== FREQUENTLY CHANGED PARAMETERS (at top for convenience) =====
# Main directory
main_directory = os.getcwd()
date_start = '20200313'
rgi_fp = main_directory + '/../00_rgi60_attribs/'
output_fp = main_directory + '/../output/'
ostrem_fp = main_directory + '/../output/ostrem_curves/'
ostrem_fn_sample = 'XXXXdebris_melt_curve.nc'

# Region of Interest Data (lat, long, elevation, hr of satellite data acquisition)
#roi = '01'
#roi = '02'
#roi = '03'
#roi = '04'
#roi = '05'
#roi = '06'
#roi = '07'
#roi = '08'
#roi = '09'
#roi = '10'
roi = '11'
#roi = '12'
#roi = 'HMA'
#roi = '16'
#roi = '17'
#roi = '18'

# ===== Debris thickness =====
#debris_thickness_all = np.array([0])
debris_thickness_all = np.array([0, 0.02])
#debris_thickness_all = np.concatenate((np.array([0]), np.arange(0,3.001,0.05)))
#debris_thickness_all[1] = 0.02

# Experiment number 3 is single run, 4 is Monte Carlo simulations 
experiment_no = 4
if experiment_no == 4:
    mc_simulations = 100
    mc_stat_cns = ['mean', 'std']
else:
    mc_simulations = 1
    mc_stat_cns = ['mean']

#eb_fp = output_fp + 'exp' + str(experiment_no) + '/' + roi + '/'
eb_fp = '/Volumes/LaCie/debris_output/exp3-20200313/' + roi + '/'

# Latitude and longitude index to run the model
#  Longitude must be 0 - 360 degrees
latlon_list = [(46.5, 10.5)] # Miage (11.03005)
#latlon_list = None
if latlon_list is None:
    latlon_unique_fp = output_fp + 'latlon_unique/'
    latlon_unique_dict = {'01':'01_latlon_unique.pkl',
                          '02':'02_latlon_unique.pkl',
                          '03':'03_latlon_unique.pkl',
                          '04':'04_latlon_unique.pkl',
                          '05':'05_latlon_unique.pkl',
                          '06':'06_latlon_unique.pkl',
                          '07':'07_latlon_unique.pkl',
                          '08':'08_latlon_unique.pkl',
                          '09':'09_latlon_unique.pkl',
                          '10':'10_latlon_unique.pkl',
                          '11':'11_latlon_unique.pkl',
                          '12':'12_latlon_unique.pkl',
                          'HMA':'HMA_latlon_unique.pkl',
                          '16':'16_latlon_unique.pkl',
                          '17':'17_latlon_unique.pkl',
                          '18':'18_latlon_unique.pkl'}
    with open(latlon_unique_fp + latlon_unique_dict[roi], 'rb') as f:
        latlon_list = pickle.load(f)


#%% ===== OTHER PARAMETERS =====

roi_latlon_dict = {'01':[71, 50, 233, 180],
                   '02':[66, 36, 256, 225],
                   '03':[84, 74, 300, 237],
                   '04':[74, 58, 300, 272],
                   '05':[84, 59, 347, 287],
                   '06':[67, 63, 347, 336],
                   '07':[81, 71, 33, 347],
                   '08':[71, 59, 34, 5],
                   '09':[82, 72, 105, 36],
                   '10':[78, 46, 188, 59],
                   '11':[48, 42, 20, -1],
                   '12':[44, 30, 53, 35],
                   'HMA':[45, 25, 105, 65],
                   '16':[12, -26, 294, 279],
                   '17':[-26, -56, 293, 285],
                   '18':[-39, -46, 176, 167]}
roi_rgidict = {'01': [1],
               '02': [2],
               '03': [3],
               '04': [4],
               '05': [5],
               '06': [6],
               '07': [7],
               '08': [8],
               '09': [9],
               '10': [10],
               '11': [11],
               '12': [12],
               'HMA':[13,14,15],
               '16':[16],
               '17':[17],
               '18':[18]}
roi_years = {'01':[1994,2018],
             '02':[2000,2018],
             '03':[2000,2018],
             '04':[2000,2018],
             '05':[2000,2018],
             '06':[2000,2018],
             '07':[2000,2018],
             '08':[2000,2018],
             '09':[2000,2018],
             '10':[2000,2018],
             '11':[2000,2018],
             '12':[2000,2018],
             'HMA':[2000,2018],
             '16':[2000,2018],
             '17':[2000,2018],
             '18':[2000,2018]}

braun_fp =  main_directory + '/../mb_data/Braun/'
mcnabb_fp = main_directory + '/../mb_data/McNabb/'
shean_fp =  main_directory + '/../mb_data/Shean/'
mb_fp_list_roi = {'01': [mcnabb_fp + '01/'],
                  '02': [braun_fp + '02/'],
                  '03': [mcnabb_fp + '03/'],
                  '04': [mcnabb_fp + '04/'],
                  '05': [mcnabb_fp + '05/'],
                  '06': [mcnabb_fp + '06/'],
                  '07': [mcnabb_fp + '07/'],
                  '08': [mcnabb_fp + '08/'],
                  '09': [mcnabb_fp + '09/'],
                  '10': [braun_fp + '10/'],
                  '11': [braun_fp + '11/'],
                  '12': [braun_fp + '12/'],
                  'HMA': [shean_fp + 'HMA/'],
                  '16': [braun_fp + 'SouthAmerica/'],
                  '17': [braun_fp + 'SouthAmerica/'],
                  '18': [braun_fp + '18/']}
mb_yrfrac_dict = {'01': [2000.6, 2018.6],
                  '02': [2000.128, 2012],
                  '03': [2000.6, 2018.6],
                  '04': [2000.6, 2018.6],
                  '05': [2000.6, 2018.6],
                  '06': [2000.6, 2018.6],
                  '07': [2000.6, 2018.6],
                  '08': [2000.6, 2018.6],
                  '09': [2000.6, 2018.6],
                  '10': [2000.128, 2012],
                  '11': [2000.128, 2013],
                  '12': [2000.128, 2012],
                  'HMA': [2000.6, 2018.6],
                  '16': None,
                  '17': None,
                  '18': [2000.128, 2013]}
#dhdt_fn_dict = {'01': main_directory + '/../mb_data/McNabb/01_rgi60_Alaska_ls_dh.tif',
#                '02': None,
#                '03': main_directory + '/../mb_data/McNabb/03_rgi60_ArcticCanadaNorth_ls_dh.tif',
#                '11': main_directory + '/../mb_data/Braun/region_11_Europe_dh_dt_on_ice.tif',
#                'HMA': main_directory + '/../mb_data/Shean/' + 
#                       'dem_align_ASTER_WV_index_2000-2018_aea_trend_3px_filt_mos_retile.tif',
#                '18': main_directory + '/../mb_data/Braun/region_18_NZ_dh_dt_on_ice.tif'}

oggm_fp = main_directory + '/../oggm_project/debris_project/'
width_min_dict = {'01': 240,
                  '02': 240,
                  '03': 240,
                  '04': 240,
                  '05': 240,
                  '06': 240,
                  '07': 240,
                  '08': 100,
                  '09': 240,
                  '10': 100,
                  '11': 100,
                  '12': 100,
                  'HMA': 240,
                  '16': 100,
                  '17': 100,
                  '18': 100}
min_bin_samp_count = 0

# ===== CLIMATE DATA =====
metdata_fp = main_directory + '/../climate_data/' + roi + '/'
#metdata_fp = '/Volumes/LaCie/ERA5/hourly/' + roi + '/'
metdata_elev_fn = 'ERA5_elev.nc'
metdata_lr_fullfn = main_directory + '/../climate_data/ERA5_lapserates_monthly.nc'
mb_binned_fp = main_directory + '/../output/mb_bins/csv/'
mb_bin_size = 10
output_fig_fp = main_directory + '/../output/figures/'
mb_binned_fp_wdebris = main_directory + '/../output/mb_bins/csv/_wdebris/'
mb_binned_fp_wdebris_hdts = main_directory + '/../output/mb_bins/csv/_wdebris_hdts/'
era5_hrly_fp = '/Volumes/LaCie_Raid/ERA5_hrly/'

oggm_ts_fp = oggm_fp + 'ts_dc/'
oggm_ts_info_fp = oggm_fp + 'ts_info_dc/'

ts_fp = main_directory + '/../output/ts_tif/' + roi + '_ts_data/'
ts_fns_fn = roi + '-ts_fns.csv'
#ts_fn = roi + '_debris_tsurfC.tif'
#ts_dayfrac_fn = roi + '_debris_dayfrac.tif'
#ts_year_fn = roi + '_debris_year.tif'
#ts_doy_fn = roi + '_debris_doy.tif'

ts_fullfns_dict = {'01': [ts_fp + '01_debris_tsurfC.tif'],
                   '02': [ts_fp + '02_debris_tsurfC.tif'],
                   '03': [ts_fp + '03_debris_tsurfC.tif'],
                   '04': [ts_fp + '04_debris_tsurfC.tif'],
                   '05': [ts_fp + '05_debris_tsurfC.tif'],
                   '06': [ts_fp + '06_debris_tsurfC.tif'],
                   '07': [ts_fp + '07_debris_tsurfC.tif'],
                   '08': [ts_fp + '08_debris_tsurfC.tif'],
                   '09': [ts_fp + '09_debris_tsurfC.tif'],
                   '10': [ts_fp + '10-01_debris_tsurfC.tif',
                          ts_fp + '10-02_debris_tsurfC.tif',
                          ts_fp + '10-03_debris_tsurfC.tif',
                          ts_fp + '10-04_debris_tsurfC.tif',
                          ts_fp + '10-05_debris_tsurfC.tif',
                          ts_fp + '10-06_debris_tsurfC.tif',
                          ts_fp + '10-07_debris_tsurfC.tif'],
                   '11': [ts_fp + '11_debris_tsurfC.tif'],
                   '12': [ts_fp + '12_debris_tsurfC.tif'],
                   'HMA': [ts_fp + 'HMA_debris_tsurfC.tif'],
                   '16': [ts_fp + '16-01_debris_tsurfC.tif',
                          ts_fp + '16-02_debris_tsurfC.tif',
                          ts_fp + '16-03_debris_tsurfC.tif'],
                   '17': [ts_fp + '17_debris_tsurfC.tif'],
                   '18': [ts_fp + '18_debris_tsurfC.tif'],}
ts_year_fullfns_dict = {'01': [ts_fp + '01_debris_year.tif'],
                        '02': [ts_fp + '02_debris_year.tif'],
                        '03': [ts_fp + '03_debris_year.tif'],
                        '04': [ts_fp + '04_debris_year.tif'],
                        '05': [ts_fp + '05_debris_year.tif'],
                        '06': [ts_fp + '06_debris_year.tif'],
                        '07': [ts_fp + '07_debris_year.tif'],
                        '08': [ts_fp + '08_debris_year.tif'],
                        '09': [ts_fp + '09_debris_year.tif'],
                        '10': [ts_fp + '10-01_debris_year.tif',
                               ts_fp + '10-02_debris_year.tif',
                               ts_fp + '10-03_debris_year.tif',
                               ts_fp + '10-04_debris_year.tif',
                               ts_fp + '10-05_debris_year.tif',
                               ts_fp + '10-06_debris_year.tif',
                               ts_fp + '10-07_debris_year.tif'],
                        '11': [ts_fp + '11_debris_year.tif'],
                        '12': [ts_fp + '12_debris_year.tif'],
                        'HMA': [ts_fp + 'HMA_debris_year.tif'],
                        '16': [ts_fp + '16-01_debris_year.tif',
                               ts_fp + '16-02_debris_year.tif',
                               ts_fp + '16-03_debris_year.tif'],
                        '17': [ts_fp + '17_debris_year.tif'],
                        '18': [ts_fp + '18_debris_year.tif'],}
ts_doy_fullfns_dict = {'01': [ts_fp + '01_debris_doy.tif'],
                       '02': [ts_fp + '02_debris_doy.tif'],
                       '03': [ts_fp + '03_debris_doy.tif'],
                       '04': [ts_fp + '04_debris_doy.tif'],
                       '05': [ts_fp + '05_debris_doy.tif'],
                       '06': [ts_fp + '06_debris_doy.tif'],
                       '07': [ts_fp + '07_debris_doy.tif'],
                       '08': [ts_fp + '08_debris_doy.tif'],
                       '09': [ts_fp + '09_debris_doy.tif'],
                       '10': [ts_fp + '10-01_debris_doy.tif',
                              ts_fp + '10-02_debris_doy.tif',
                              ts_fp + '10-03_debris_doy.tif',
                              ts_fp + '10-04_debris_doy.tif',
                              ts_fp + '10-05_debris_doy.tif',
                              ts_fp + '10-06_debris_doy.tif',
                              ts_fp + '10-07_debris_doy.tif'],
                       '11': [ts_fp + '11_debris_doy.tif'],
                       '12': [ts_fp + '12_debris_doy.tif'],
                       'HMA': [ts_fp + 'HMA_debris_doy.tif'],
                       '16': [ts_fp + '16-01_debris_doy.tif',
                              ts_fp + '16-02_debris_doy.tif',
                              ts_fp + '16-03_debris_doy.tif'],
                       '17': [ts_fp + '17_debris_doy.tif'],
                       '18': [ts_fp + '18_debris_doy.tif'],}
ts_dayfrac_fullfns_dict = {'01': [ts_fp + '01_debris_dayfrac.tif'],
                           '02': [ts_fp + '02_debris_dayfrac.tif'],
                           '03': [ts_fp + '03_debris_dayfrac.tif'],
                           '04': [ts_fp + '04_debris_dayfrac.tif'],
                           '05': [ts_fp + '05_debris_dayfrac.tif'],
                           '06': [ts_fp + '06_debris_dayfrac.tif'],
                           '07': [ts_fp + '07_debris_dayfrac.tif'],
                           '08': [ts_fp + '08_debris_dayfrac.tif'],
                           '09': [ts_fp + '09_debris_dayfrac.tif'],
                           '10': [ts_fp + '10-01_debris_dayfrac.tif',
                                  ts_fp + '10-02_debris_dayfrac.tif',
                                  ts_fp + '10-03_debris_dayfrac.tif',
                                  ts_fp + '10-04_debris_dayfrac.tif',
                                  ts_fp + '10-05_debris_dayfrac.tif',
                                  ts_fp + '10-06_debris_dayfrac.tif',
                                  ts_fp + '10-07_debris_dayfrac.tif'],
                           '11': [ts_fp + '11_debris_dayfrac.tif'],
                           '12': [ts_fp + '12_debris_dayfrac.tif'],
                           'HMA': [ts_fp + 'HMA_debris_dayfrac.tif'],
                           '16': [ts_fp + '16-01_debris_dayfrac.tif',
                                  ts_fp + '16-02_debris_dayfrac.tif',
                                  ts_fp + '16-03_debris_dayfrac.tif'],
                           '17': [ts_fp + '17_debris_dayfrac.tif'],
                           '18': [ts_fp + '18_debris_dayfrac.tif'],}
ts_stats_res = 50 # common resolution needed such that resolution does not interfere with regional stats
output_ts_csv_ending = '_ts_hd_opt.csv'
tscurve_fp = ts_fp + '../ts_curves/'
output_ts_fn_sample = 'XXXXdebris_ts_curve.nc'
hd_fp = ts_fp + '../hd_tifs/' + roi + '/'
hd_fn_sample = 'XXXX_hdts_m.tif'
mf_fp = ts_fp + 'hd_tifs/_meltfactor/'
mf_fn_sample = 'XXXX_meltfactor.tif'
hd_max = 3
vel_threshold = 7.5
debrisperc_threshold = 50

# Debris datasets
debriscover_fp = main_directory + '/../scherler_debris/LS8_2013-2017_RATIO/fixed_v2/'
debriscover_fn_dict = {'01':'01_rgi60_L8ratio_fixed_v2.shp',
                       '02':'02_rgi60_L8ratio_fixed_v2.shp',
                       '03':'03_rgi60_L8ratio_fixed_v2.shp',
                       '04':'04_rgi60_L8ratio_fixed_v2.shp',
                       '05':'05_rgi60_L8ratio_fixed_v2.shp',
                       '06':'06_rgi60_L8ratio_fixed_v2.shp',
                       '07':'07_rgi60_L8ratio_fixed_v2.shp',
                       '08':'08_rgi60_L8ratio_fixed_v2.shp',
                       '09':'09_rgi60_L8ratio_fixed_v2.shp',
                       '10':'10_rgi60_L8ratio_fixed_v2.shp',
                       '11':'11_rgi60_L8ratio_fixed_v2.shp',
                       '12':'12_rgi60_L8ratio_fixed_v2.shp',
                       'HMA':'HMA_rgi60_L8ratio_fixed_v2.shp',
                       '16':'16_rgi60_L8ratio_fixed_v2.shp',
                       '17':'17_rgi60_L8ratio_fixed_v2.shp',
                       '18':'18_rgi60_L8ratio_fixed_v2.shp'}
dc_percarea_threshold = 5   # percent area threshold (%)
dc_area_threshold = 1       # debris-covered area threshold (km2) for large glaciers with low % but significant debris
min_glac_area = 2           # minimum glacier area (only work with large glaciers)
term_elevrange_perc = 0.1   # Terminus elevation range (ex. 0.1 = lower 10% of glacier)
term_area_perc = 10         # Terminus area percentage (ex. 10 = lower 10% of glacier by area)

# Glacier data
glac_shp_fn_dict = {
        '01': main_directory + '/../../../HiMAT/RGI/rgi60/01_rgi60_Alaska/01_rgi60_Alaska.shp',
        '02': main_directory + '/../../../HiMAT/RGI/rgi60/02_rgi60_WesternCanadaUS/02_rgi60_WesternCanadaUS.shp',
        '03': main_directory + '/../../../HiMAT/RGI/rgi60/03_rgi60_ArcticCanadaNorth/03_rgi60_ArcticCanadaNorth.shp',
        '04': main_directory + '/../../../HiMAT/RGI/rgi60/04_rgi60_ArcticCanadaSouth/04_rgi60_ArcticCanadaSouth.shp',
        '05': main_directory + '/../../../HiMAT/RGI/rgi60/05_rgi60_GreenlandPeriphery/05_rgi60_GreenlandPeriphery.shp',
        '06': main_directory + '/../../../HiMAT/RGI/rgi60/06_rgi60_Iceland/06_rgi60_Iceland.shp',
        '07': main_directory + '/../../../HiMAT/RGI/rgi60/07_rgi60_Svalbard/07_rgi60_Svalbard.shp',
        '08': main_directory + '/../../../HiMAT/RGI/rgi60/08_rgi60_Scandinavia/08_rgi60_Scandinavia.shp',
        '09': main_directory + '/../../../HiMAT/RGI/rgi60/09_rgi60_RussianArctic/09_rgi60_RussianArctic.shp',
        '10': main_directory + '/../../../HiMAT/RGI/rgi60/10_rgi60_NorthAsia/10_rgi60_NorthAsia.shp',
        '11': main_directory + '/../../../HiMAT/RGI/rgi60/11_rgi60_CentralEurope/11_rgi60_CentralEurope.shp',
        '12': main_directory + '/../../../HiMAT/RGI/rgi60/12_rgi60_CaucasusMiddleEast/12_rgi60_CaucasusMiddleEast.shp',
        '13': main_directory + '/../../../HiMAT/RGI/rgi60/13_rgi60_CentralAsia/13_rgi60_CentralAsia.shp',
        '14': main_directory + '/../../../HiMAT/RGI/rgi60/14_rgi60_SouthAsiaWest/14_rgi60_SouthAsiaWest.shp',
        '15': main_directory + '/../../../HiMAT/RGI/rgi60/15_rgi60_SouthAsiaEast/15_rgi60_SouthAsiaEast.shp',
        '16': main_directory + '/../../../HiMAT/RGI/rgi60/16_rgi60_LowLatitudes/16_rgi60_LowLatitudes.shp',
        '17': main_directory + '/../../../HiMAT/RGI/rgi60/17_rgi60_SouthernAndes/17_rgi60_SouthernAndes.shp',
        '18': main_directory + '/../../../HiMAT/RGI/rgi60/18_rgi60_NewZealand/18_rgi60_NewZealand.shp'}
glac_shp_proj_fp = output_fp + 'glac_shp_proj/'
if os.path.exists(glac_shp_proj_fp) == False:
    os.makedirs(glac_shp_proj_fp)
#DEM
z1_dir_sample = main_directory + '/../oggm_dems/dem_qc/RGI60-XXXX/'
z1_fn_sample = 'RGI60-XXXX-dem.tif'
# Surface velocity
sat_img_dir = main_directory + '/../../../Satellite_Images/'
vx_dir_dict_list = {'01': [sat_img_dir + 'ITS_Live/ALA_G0120_0000_vx.tif'],
                    '02': [sat_img_dir + 'ITS_Live/ALA_G0120_0000_vx.tif'],
                    '03': [sat_img_dir + 'ITS_Live/CAN_G0120_0000_vx.tif'],
                    '04': [sat_img_dir + 'ITS_Live/CAN_G0120_0000_vx.tif'],
                    '05': [sat_img_dir + 'ITS_Live/GRE_G0120_0000_vx.tif'],
                    '06': [sat_img_dir + 'ITS_Live/ICE_G0120_0000_vx.tif'],
                    '07': [sat_img_dir + 'ITS_Live/SRA_G0120_0000_vx.tif'],
                    '08': [sat_img_dir + 'romain_velocity/NORWAY/Mosaic__vxo.tif'],
                    '09': [sat_img_dir + 'ITS_Live/SRA_G0120_0000_vx.tif'],
                    '10': [sat_img_dir + 'romain_velocity/Asia_north/Mosaic__vxo.tif'],
                    '11': [sat_img_dir + 'romain_velocity/ALPS/Mosaic__vxo.tif'],
                    '12': [sat_img_dir + 'romain_velocity/Caucasus/Mosaic__vxo.tif'],
                    'HMA': [sat_img_dir + 'ITS_Live/HMA_G0120_0000_vx.tif'],
                    '16': [sat_img_dir + 'romain_velocity/AndesBlanca/Mosaic__vxo.tif'],
                    '17': [sat_img_dir + 'romain_velocity/Cordillera_South/Mosaic__vxo.tif'],
                    '18': [sat_img_dir + 'romain_velocity/NEWZEALAND/Mosaic__vxo.tif']
                    }
dhdt_vel_fns_fp = output_fp + 'dhdt_vel_fns/'
dhdt_vel_fns_fn = 'XXXX-dhdt_vel_fns.csv'

# Emergence Velocity data
min_glac_area_writeout=0
min_valid_area_perc = 0
buff_dist = 1000
#emvel_bin_width = 50
emvel_filter_pixsize = 3
#Surface to column average velocity scaling
v_col_f = 0.8
output_emvel_csv_ending = '_emvel_stats_woffset.csv'
outdir_emvel_fp = output_fp + 'csv/'

# Simulation data
startyear = roi_years[roi][0]
endyear = roi_years[roi][1]
timezone = 0
metdata_fn_sample = roi + '_ERA5-metdata-XXXX' + str(startyear) + '_' + str(endyear) + '.nc'

#%%
# Simulation data
roi_datedict = {'01': ['1994-01-01', '2018-12-31'],
                '02': ['2000-01-01', '2018-12-31'],
                '03': ['2000-01-01', '2018-12-31'],
                '04': ['2000-01-01', '2018-12-31'],
                '05': ['2000-01-01', '2018-12-31'],
                '06': ['2000-01-01', '2018-12-31'],
                '07': ['2000-01-01', '2018-12-31'],
                '08': ['2000-01-01', '2018-12-31'],
                '09': ['2000-01-01', '2018-12-31'],
                '10': ['2000-01-01', '2018-12-31'],
                '11': ['2000-01-01', '2018-12-31'],
                '12': ['2000-01-01', '2018-12-31'],
                'HMA': ['2000-01-01', '2018-12-31'],
                '16': ['2000-01-01', '2018-12-31'],
                '17': ['2000-01-01', '2018-12-31'],
                '18': ['2000-01-01', '2018-12-31']}
start_date = roi_datedict[roi][0]  # start date for debris_ts_model.py
end_date = roi_datedict[roi][1]     # end date for debris_ts_model.py
fn_prefix = 'Rounce2015_' + roi + '-'
elev_cns = ['zmean']
#elev_cns = ['zmean', 'zstdlow', 'zstdhigh']

# ===== Debris properties =====
# Surface roughness, thermal conductivity, and albedo
if experiment_no == 3:
    z0_random = np.array([0.016])
    k_random = np.array([1.])
    albedo_random = np.array([0.2])
    z0_random_ice = np.array([0.002])             # Brock etal (2006)
    z0_random_snow = np.array([z0_random_ice])    # Hock and Holmgren (2005) - 0.0001 to 0.0027, z0_ice = z0_snow
    albedo_random_ice = np.array([0.4])           # Gardner and Sharp (2010); Hock (2005)
    sin_factor_random = np.array([1.])            # multiplicative factor
    
elif experiment_no == 4:
    debris_properties_fp = output_fp + 'debris_properties/'
    debris_properties_fn = 'debris_properties_global.csv'
    if os.path.exists(debris_properties_fp + debris_properties_fn):
        debris_properties_df = pd.read_csv(debris_properties_fp + debris_properties_fn)
        albedo_random = debris_properties_df['albedo'].values
        z0_random = debris_properties_df['z0'].values
        k_random = debris_properties_df['k'].values
        albedo_random_ice = debris_properties_df['albedo_ice'].values
        z0_random_ice = debris_properties_df['z0_ice'].values
        z0_random_snow = z0_random_ice
        sin_factor_random = debris_properties_df['Sin_factor'].values
    else:
        # Albedo (uniform distribution 0.1 - 0.3)
        albedo_random = np.random.uniform(low=0.1, high=0.3, size=mc_simulations)
        # Surface roughness (m)
        z0_random = np.random.uniform(low=0.008, high=0.024, size=mc_simulations)
        # Thermal conductivity (uniform distribution 0.5 - 1.5 W m-1 K-1)
        k_random = np.random.uniform(low=0.5, high=1.5, size=mc_simulations)
        # Clean ice albedo (uniform distribution 0.3 - 0.5)
        albedo_random_ice = np.random.uniform(low=0.3, high=0.5, size=mc_simulations)
        # Clean ice surface roughness (m)
        z0_random_ice = np.random.uniform(low=0.0001, high=0.004, size=mc_simulations)
        z0_random_snow = z0_random_ice
        # Sin multiplicative factor to adjust Sin for topography, etc. (uniform distribution 0.8 - 1.2)
        sin_factor_random = np.random.uniform(low=0.8, high=1.2, size=mc_simulations)
        
        debris_properties_values = np.column_stack((albedo_random, k_random, z0_random, 
                                                    albedo_random_ice, z0_random_ice, sin_factor_random))
        # Export properties
        debris_properties_cns = ['albedo', 'k', 'z0', 'albedo_ice', 'z0_ice', 'Sin_factor']
        debris_properties_df = pd.DataFrame(debris_properties_values, columns=debris_properties_cns)
        if not os.path.exists(debris_properties_fp):
            os.makedirs(debris_properties_fp)
        debris_properties_df.to_csv(debris_properties_fp + debris_properties_fn, index=False)

# Extra
#debris_albedo = 0.2     # -, debris albedo
za = 2                  # m, height of air temperature instrument
zw = 10                 # m, height of wind instrument

# Snow model parameters
option_snow_fromAWS = 0 # Switch to use snow depth as opposed to snow fall
option_snow = 1         # Switch to use snow model (1) or not (0)
Tsnow_threshold = 274.15      # Snow temperature threshold [K] - Regine get source
snow_min = 0.0001       # minimum snowfall (mwe) to include snow on surface; since the density of falling snow is
                        # much less (~50-100 kg m-3) 0.0001 m of snow w.e. will produce 0.001 - 0.002 m of snow
rain_min = 0.0001
option_lr_fromdata = 1  # Switch to use lapse rate from data (e.g., ERA5) or specified value

#%%
delta_t = 60*60         # s, time step of AWS
slope_AWS_deg = 0       # assuming AWS flat on top of ridge or Sky View Factor = 1
aspect_AWS_deg = 0      # deg CW from N

#%% Constants
row_d = 2700                # Density of debris (kg m-3)
c_d = 750                   # Specific heat capacity of debris (J kg-1 K-1)
I0 = 1368                   # Solar constant (W/m2)
transmissivity=0.75         # Vertical atmospheric clear-sky transmissivity (assumed)
emissivity = 0.95           # Emissivity of debris surface (Nicholson and Benn, 2006)
P0 = 101325                 # Standard Atmospheric Pressure (Pa)
density_air_0 = 1.29        # Density of air (kg/m3)
density_water = 1000        # Density of water (kg/m3)
density_ice = 900           # Density of ice (kg/m3) (Nicholson and Benn, 2006)
lapserate = -0.0055          # Temperature Lapse Rate (K/m)
Kvk = 0.41                  # Von Karman's Constant
Lv = 2.49e6                 # Latent Heat of Vaporation of water (J/kg)
Lf = 3.335e5                # Latent Heat of Fusion of Water (J/kg)
Ls = 2.834e6                # Latent Heat of Sublimation (J/kg)
cA = 1005                   # Specific Heat Capacity of Air (J kg-1 K-1)
cW = 4.179e3                # Specific Heat Capacity of Water (J kg-1 K-1)
cSnow = 2.09e3              # Specific Heat Capacity of Snow (J kg-1 K-1)
R_const = 461               # Gas Constant for water vapor (J kg-1 K-1)
Rd = 287                    # Dry gas constant (J kg-1 K-1)
stefan_boltzmann = 5.67e-8  # Sefan-Bolzmann constant (W m-2 K-4)

# Snow melt model constants
snow_c_v = 0.2              # sensitivity of visible albedo to snow surface aging (Tarboten and Luce, 1996)
snow_c_ir = 0.5             # sensitivity of near infrared albedo to snow surface aging (Tarboten and Luce, 1996)
albedo_vo = 0.85            # fresh snow reflectance for visible band
albedo_iro = 0.65           # fresh snow reflectance for near infrared band
snow_tau_0 = 1e6            # non-dimensional snow surface age constant [s]
emissivity_snow = 0.99      # emissivity of snow (Tarboten and Luce, 1996); Collier etal (2014) use 0.97
eS_snow = 610.5             # Saturated vapor pressure of snow (Pa) (Colbeck 1990)
k_snow = 0.10               # Rahimi and Konrad (2012), Sturm etal (2002), Sturm etal (1997)
#density_snow = 150         # Density of snow (kg/m3) - Lejeune et al. (2007)
#albedo_snow = 0.75         # Collier etal (2014)

#k_ice = 1

# Newton-Raphson Method constants
n_iter_max = 100

#%% FUNCTIONS
def selectglaciersrgitable(glac_no=None,
                           rgi_regionsO1=None,
                           rgi_regionsO2=None,
                           rgi_glac_number=None,
                           rgi_fp=rgi_fp,
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
