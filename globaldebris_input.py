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
import xarray as xr


#%%
# Main directory
main_directory = os.getcwd()
rgi_fp = main_directory + '/../00_rgi60_attribs/'
output_fp = main_directory + '/../output/'
ostrem_fp = main_directory + '/../output/ostrem_curves/'
ostrem_fn_sample = 'XXXX_debris_melt_curve.nc'

# Region of Interest Data (lat, long, elevation, hr of satellite data acquisition)
roi = '01'
#roi = '11'
#roi = 'HMA'
#roi = '17'
roi_latlon_dict = {'01':[71, 50, 233, 180],
                   '11':[48, 42, 20, -1],
                   'HMA':[45, 25, 105, 65],
                   '17':[-26, -56, 293, 285]}
roi_rgidict = {'01': [1],
               '11': [11],
               'HMA':[13,14,15],
               '17':[17]}
roi_years = {'01':[1994,2018],
             '11':[2000,2018],
             'HMA':[2000,2018],
             '17':[2000,2018]}
latlon_unique_fp = output_fp + 'latlon_unique/'
latlon_unique_dict = {'01':'01_latlon_unique.pkl',
                      '11':'11_latlon_unique.pkl',
                      'HMA':'HMA_latlon_unique.pkl',
                      '17':'17_latlon_unique.pkl'}

#mb_datasets = ['larsen', 'braun']
mb_datasets = ['shean']

# Climate data
metdata_fp = main_directory + '/../climate_data/' + roi + '/'
metdata_elev_fn = 'ERA5_elev.nc'
mb_binned_fp = main_directory + '/../output/mb_bins/csv/'
#mb_binned_fp_wdebris = main_directory + '/../mb_data/Shean_2019_0213/mb_combined_20190213_nmad_bins/_wdebris/'
mb_binned_fp_wdebris = main_directory + '/../output/mb_bins/csv/_wdebris/'
era5_hrly_fp = '/Volumes/LaCie_Raid/ERA5_hrly/'

ts_fp = main_directory + '/../output/ts_tif/'
ts_fn_dict = {'01':'01_debris_tsurfC.tif',
              'HMA':'hma_debris_tsurfC_wbuffer.tif'}
ts_dayfrac_fn_dict = {'01':'01_debris_dayfrac.tif',
                      'HMA':'hma_debris_dayfrac.tif'}
ts_year_fn_dict = {'01':'01_debris_year.tif',
                   'HMA':'hma_debris_year.tif'}
ts_doy_fn_dict = {'01':'01_debris_doy.tif',
                  'HMA':'hma_debris_doy.tif'}
ts_stats_res = 50 # common resolution needed such that resolution does not interfere with regional stats
#ts_fn = ts_fn_dict[roi]
output_ts_csv_ending = '_ts_hd_opt.csv'
tscurve_fp = ts_fp + 'ts_curves/'
output_ts_fn_sample = 'XXXX_debris_ts_curve.nc'
hd_fp = ts_fp + 'hd_tifs/'
hd_fn_sample = 'XXXX_hdts_m.tif'
hd_max = 2.5
vel_threshold = 7.5
debrisperc_threshold = 50

# Debris datasets
debriscover_fp = main_directory + '/../scherler_debris/LS8_2013-2017_RATIO/fixed_v2/'
debriscover_fn_dict = {'01':'01_rgi60_L8ratio_fixed_v2.shp',
                       '11':'11_rgi60_L8ratio_fixed_v2.shp',
                       'HMA':'HMA_rgi60_L8ratio_fixed_v2.shp',
                       '17':'17_rgi60_L8ratio_fixed_v2.shp'}
dc_percarea_threshold = 5   # percent area threshold (%)
dc_area_threshold = 1       # debris-covered area threshold (km2) for large glaciers with low % but significant debris
min_glac_area = 2                # minimum glacier area (only work with large glaciers)

# Glacier data
glac_shp_fn_dict = {
        '01': main_directory + '/../../../HiMAT/RGI/rgi60/01_rgi60_Alaska/01_rgi60_Alaska.shp',
        '11': main_directory + '/../../../HiMAT/RGI/rgi60/11_rgi60_CentralEurope/11_rgi60_CentralEurope.shp',
        '13': main_directory + '/../../../HiMAT/RGI/rgi60/13_rgi60_CentralAsia/13_rgi60_CentralAsia.shp',
        '14': main_directory + '/../../../HiMAT/RGI/rgi60/14_rgi60_SouthAsiaWest/14_rgi60_SouthAsiaWest.shp',
        '15': main_directory + '/../../../HiMAT/RGI/rgi60/15_rgi60_SouthAsiaEast/15_rgi60_SouthAsiaEast.shp',
        '17': main_directory + '/../../../HiMAT/RGI/rgi60/17_rgi60_SouthernAndes/17_rgi60_SouthernAndes.shp'}
glac_shp_proj_fp = output_fp + 'glac_shp_proj/'
if os.path.exists(glac_shp_proj_fp) == False:
    os.makedirs(glac_shp_proj_fp)
#DEM
z1_dir_sample = main_directory + '/../../../HiMAT/IceThickness_Farinotti/surface_DEMs_RGI60/surface_DEMs_RGI60-XXXX/'
z1_fn_sample = 'surface_DEM_RGI60-XXXX.tif'
z1_backup_dict = {'01': main_directory + '/../../../Satellite_Images/Alaska_albers_V3_mac/Alaska_albers_V3.tif',
                  '11': None,
                  'HMA': None,
                  '17': None}
# Ice thickness
huss_dir_sample = (
        main_directory + '/../../../HiMAT/IceThickness_Farinotti/composite_thickness_RGI60-all_regions/RGI60-XXXX/')
huss_fn_sample = 'RGI60-XXXX_thickness.tif'
# Surface velocity
v_dir = main_directory + '/../../../Satellite_Images/ITS_Live/'
vx_fn_dict = {'01': 'ALA_G0120_0000_vx.tif',
              '11': None,
              'HMA': 'HMA_G0120_0000_vx.tif',
              '17': None
              }
vy_fn_dict = {'01': 'ALA_G0120_0000_vy.tif',
              '11': None,
              'HMA': 'HMA_G0120_0000_vy.tif',
              '17': None}

# Mass balance data
larsen_fp = main_directory + '/../mb_data/Larsen/'
larsen_binned_fp = larsen_fp + 'binned_data/'
larsen_fn = 'larsen2015_supplementdata_wRGIIds_v3.csv'

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
debris_elevstats_fullfn = main_directory + '/../hma_data/' + roi + '_debris_elevstats.nc'

# Latitude and longitude index to run the model
#  Longitude must be 0 - 360 degrees
latlon_list_raw = 'all'
#latlon_list_raw = None
if latlon_list_raw == 'all':
    with open(latlon_unique_fp + latlon_unique_dict[roi], 'rb') as f:
        latlon_list = pickle.load(f)
elif latlon_list_raw is not None:
    ds_elevstats = xr.open_dataset(debris_elevstats_fullfn)
    lat_list_raw = np.array([x[0] for x in latlon_list_raw])
    lon_list_raw = np.array([x[1] for x in latlon_list_raw])
    #  argmin() finds the minimum distance between the glacier lat/lon and the GCM pixel
    lat_nearidx = np.abs(lat_list_raw[:,np.newaxis] - ds_elevstats['latitude'][:].values).argmin(axis=1)
    lon_nearidx = np.abs(lon_list_raw[:,np.newaxis] - ds_elevstats['longitude'][:].values).argmin(axis=1)
    latlon_nearidx = list(zip(lat_nearidx, lon_nearidx))
    latlon_nearidx_unique = sorted(list(set(latlon_nearidx)))
    latidx_list_unique = [x[0] for x in latlon_nearidx_unique]
    lonidx_list_unique= [x[1] for x in latlon_nearidx_unique]
    lat_list = ds_elevstats.latitude[latidx_list_unique].values
    lon_list = ds_elevstats.longitude[lonidx_list_unique].values
    latlon_list = list(tuple(zip(list(lat_list), list(lon_list))))

#latlon_list = latlon_list[0:5]
#latlon_list = [latlon_list[0]]
#latlon_list = [(61.5, 217.0)]
#latlon_list = [(55.25, 230.50)]
#latlon_list = [(28.0, 86.75)]

#%%
# Simulation data
roi_datedict = {'01': ['1994-01-01', '2018-12-31'],
                '11': ['2000-01-01', '2018-12-31'],
                'HMA': ['2000-05-28', '2018-05-28'],
                '17': ['2000-01-01', '2017-12-31']}
start_date = roi_datedict[roi][0]  # start date for debris_ts_model.py
end_date = roi_datedict[roi][1]     # end date for debris_ts_model.py
#start_date = '2000-05-28'   # start date for debris_ts_model.py
#end_date = '2018-05-28'     # end date for debris_ts_model.py
fn_prefix = 'Rounce2015_' + roi + '-'
elev_cns = ['zmean']
#elev_cns = ['zmean', 'zstdlow', 'zstdhigh']

# Output info
output_option = 2           # 1: csv of all fluxes and internal temps, 2: netcdf of melt and ts
if output_option == 2:
    mc_stat_cns = ['mean']
elif output_option == 3:
    mc_stat_cns = ['mean', 'std', '25%', '75%']
    print('\nSTOP!!!!! NEED TO STORE ATTRIBUTES FOR STATISTICS!!!!\n\n')
date_start = '20200113'

# ===== Debris properties =====
experiment_no = 3
# Debris thickness
debris_thickness_all = np.array([0.2])
#debris_thickness_all = np.array([0, 0.02])
#debris_thickness_all = np.concatenate((np.array([0]), np.arange(0,3.001,0.05)))
#debris_thickness_all[1] = 0.02

# Surface roughness, thermal conductivity, and albedo
debris_properties_fullfn = main_directory + '/../hma_data/hma_debris_properties.csv'
debris_properties = np.genfromtxt(debris_properties_fullfn, delimiter=',', skip_header=1)
if experiment_no == 3:
    z0_random = np.array([0.016])
    k_random = np.array([1.])
    albedo_random = np.array([0.2])
elif experiment_no == 4:    
    z0_random = debris_properties[:,1]
    k_random = debris_properties[:,2]
    albedo_random = np.array([debris_properties[:,5]])
    print('\n\nNEED TO MAKE DEBRIS PROPERTIES RANDOM FOR MC SIMULATIONS\n\n')

#met_data_netcdf_fp = '/Volumes/LaCie_Raid/ERA5_hrly/'
#orog_data_fullfn = '/Users/davidrounce/Documents/Dave_Rounce/HiMAT/Climate_data/ERA5/ERA5_geopotential_monthly.nc'

#pkl_fp = main_directory + '/hma_data/pkl/'
#pkl_fullfn_sample = pkl_fp + roi + '_ERA5-metadata_' + start_date + '_'+ end_date + '_XXXX.pkl' 

# Dates of satellite temperature data
#ts_dates = ['2015-09-30']      
#ts_hr = roi_dict[roi][3]            # hour of satellite temperature data acquisition  
    
# Extra 
debris_albedo = 0.2     # -, debris albedo
za = 2                  # m, height of air temperature instrument
zw = 10                 # m, height of wind instrument

# Snow model parameters
option_snow_fromAWS = 0 # Switch to use snow depth as opposed to snow fall
option_snow = 1         # Switch to use snow model (1) or not (0)
Tsnow_threshold = 274.15      # Snow temperature threshold [K] - Regine get source
snow_min = 0.0001       # minimum snowfall (mwe) to include snow on surface; since the density of falling snow is 
                        # much less (~50-100 kg m-3) 0.0001 m of snow w.e. will produce 0.001 - 0.002 m of snow
rain_min = 0.0001



#%%
#debris_elev = roi_dict[roi][0]       # m a.s.l. of the debris modeling
lr_gcm = -0.0065        # lapse rate from gcm to glacier [K m-1]    
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
lapserate = 0.0065          # Temperature Lapse Rate (K/m)
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
z0_snow = 0.002             # Brock etal (2006)
emissivity_snow = 0.99      # emissivity of snow (Tarboten and Luce, 1996); Collier etal (2014) use 0.97
eS_snow = 610.5             # Saturated vapor pressure of snow (Pa) (Colbeck 1990)
k_snow = 0.10               # Rahimi and Konrad (2012), Sturm etal (2002), Sturm etal (1997)
#density_snow = 150         # Density of snow (kg/m3) - Lejeune et al. (2007)
#albedo_snow = 0.75         # Collier etal (2014)

# Clean ice melt model constants
albedo_ice = 0.4            # Gardner and Sharp (2010); Hock (2005)
z0_ice = z0_snow            # Hock and Holmgren (2005) - vary from 0.0001 to 0.0027, assume z0_ice = z0_snow
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