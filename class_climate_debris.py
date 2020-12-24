"""class of climate data and functions associated with manipulating the dataset to be in the proper format"""

# External libraries
import pandas as pd
import numpy as np
import xarray as xr
# Local libraries
import debrisglobal.globaldebris_input as debris_prms

class GCM():
    """
    Global climate model data properties and functions used to automatically retrieve data.
    
    Attributes
    ----------
    name : str
        name of climate dataset.
    rcp_scenario : str
        rcp scenario (example: 'rcp26')
    """
    def __init__(self, 
                 name=str(),
                 rcp_scenario=str()):
        """
        Add variable name and specific properties associated with each gcm.
        """
        
        # Source of climate data
        self.name = name
        # Set parameters for ERA5, ERA-Interim, and CMIP5 netcdf files
        if self.name == 'ERA5':
            # Variable names
            self.temp_vn = 't2m'
            self.tempstd_vn = 't2m_std'
            self.prec_vn = 'tp'
            self.elev_vn = 'z'
            self.lat_vn = 'latitude'
            self.lon_vn = 'longitude'
            self.time_vn = 'time'
            self.lr_vn = 'lapserate'
            # Variable filenames
            self.temp_fn = debris_prms.era5_temp_fn
            self.tempstd_fn = debris_prms.era5_tempstd_fn
            self.prec_fn = debris_prms.era5_prec_fn
            self.elev_fn = debris_prms.era5_elev_fn
            self.lr_fn = debris_prms.era5_lr_fn
            # Variable filepaths
            self.var_fp = debris_prms.era5_fp
            self.fx_fp = debris_prms.era5_fp
            # Extra information
            self.timestep = debris_prms.timestep
            self.rgi_lat_colname=debris_prms.rgi_lat_colname
            self.rgi_lon_colname=debris_prms.rgi_lon_colname
            
            
    def importGCMfxnearestneighbor_xarray(self, filename, vn, main_glac_rgi):
        """
        Import time invariant (constant) variables and extract nearest neighbor.
        
        Note: cmip5 data used surface height, while ERA-Interim data is geopotential
        
        Parameters
        ----------
        filename : str
            filename of variable
        variablename : str
            variable name
        main_glac_rgi : pandas dataframe
            dataframe containing relevant rgi glacier information
        
        Returns
        -------
        glac_variable : numpy array
            array of nearest neighbor values for all the glaciers in model run (rows=glaciers, column=variable)
        """
        # Import netcdf file
        data = xr.open_dataset(self.fx_fp + filename)
        glac_variable = np.zeros(main_glac_rgi.shape[0])
        # If time dimension included, then set the time index (required for ERA Interim, but not for CMIP5 or COAWST)
        if 'time' in data[vn].coords:
            time_idx = 0
            #  ERA Interim has only 1 value of time, so index is 0
        # Find Nearest Neighbor
        #  argmin() finds the minimum distance between the glacier lat/lon and the GCM pixel
        lat_nearidx = (np.abs(main_glac_rgi[self.rgi_lat_colname].values[:,np.newaxis] - 
                              data.variables[self.lat_vn][:].values).argmin(axis=1))
        lon_nearidx = (np.abs(main_glac_rgi[self.rgi_lon_colname].values[:,np.newaxis] - 
                              data.variables[self.lon_vn][:].values).argmin(axis=1))
        
        latlon_nearidx = list(zip(lat_nearidx, lon_nearidx))
        latlon_nearidx_unique = list(set(latlon_nearidx))
        
        glac_variable_dict = {}
        for latlon in latlon_nearidx_unique:
            try:
                glac_variable_dict[latlon] = data[vn][time_idx, latlon[0], latlon[1]].values
            except:
                glac_variable_dict[latlon] = data[vn][latlon[0], latlon[1]].values
        
        glac_variable = np.array([glac_variable_dict[x] for x in latlon_nearidx])    
            
        
        # Correct units if necessary (CMIP5 already in m a.s.l., ERA Interim is geopotential [m2 s-2])
        if vn == self.elev_vn:
            # If the variable has units associated with geopotential, then convert to m.a.s.l (ERA Interim)
            if 'units' in data[vn].attrs and (data[vn].attrs['units'] == 'm**2 s**-2'):  
                # Convert m2 s-2 to m by dividing by gravity (ERA Interim states to use 9.80665)
                glac_variable = glac_variable / 9.80665
            # Elseif units already in m.a.s.l., then continue
            elif 'units' in data[vn].attrs and data[vn].attrs['units'] == 'm':
                pass
            # Otherwise, provide warning
            else:
                print('Check units of elevation from GCM is m.')
        return glac_variable

    
    def importGCMvarnearestneighbor_xarray(self, filename, vn, main_glac_rgi, dates_table):
        """
        Import time series of variables and extract nearest neighbor.
        
        Note: "NG" refers to a homogenized "new generation" of products from ETH-Zurich.
              The function is setup to select netcdf data using the dimensions: time, latitude, longitude (in that 
              order). Prior to running the script, the user must check that this is the correct order of the dimensions 
              and the user should open the netcdf file to determine the names of each dimension as they may vary.
        
        Parameters
        ----------
        filename : str
            filename of variable
        vn : str
            variable name
        main_glac_rgi : pandas dataframe
            dataframe containing relevant rgi glacier information
        dates_table: pandas dataframe
            dataframe containing dates of model run
        
        Returns
        -------
        glac_variable_series : numpy array
            array of nearest neighbor values for all the glaciers in model run (rows=glaciers, columns=time series)
        time_series : numpy array
            array of dates associated with the meteorological data (may differ slightly from those in table for monthly
            timestep, i.e., be from the beginning/middle/end of month)
        """
        # Import netcdf file
        data = xr.open_dataset(self.var_fp + filename)
        glac_variable_series = np.zeros((main_glac_rgi.shape[0],dates_table.shape[0]))
        # Determine the correct time indices
        if self.timestep == 'monthly':
            start_idx = (np.where(pd.Series(data[self.time_vn]).apply(lambda x: x.strftime('%Y-%m')) == 
                                  dates_table['date'].apply(lambda x: x.strftime('%Y-%m'))[0]))[0][0]
            end_idx = (np.where(pd.Series(data[self.time_vn]).apply(lambda x: x.strftime('%Y-%m')) == 
                                dates_table['date']
                                .apply(lambda x: x.strftime('%Y-%m'))[dates_table.shape[0] - 1]))[0][0]              
            #  np.where finds the index position where to values are equal
            #  pd.Series(data.variables[gcm_time_varname]) creates a pandas series of the time variable associated with
            #  the netcdf
            #  .apply(lambda x: x.strftime('%Y-%m')) converts the timestamp to a string with YYYY-MM to enable the 
            #  comparison
            #    > different climate dta can have different date formats, so this standardization for comparison is 
            #      important
            #      ex. monthly data may provide date on 1st of month or middle of month, so YYYY-MM-DD would not work
            #  The same processing is done for the dates_table['date'] to facilitate the comparison
            #  [0] is used to access the first date
            #  dates_table.shape[0] - 1 is used to access the last date
            #  The final indexing [0][0] is used to access the value, which is inside of an array containing extraneous 
            #  information

        # Extract the time series
        time_series = pd.Series(data[self.time_vn][start_idx:end_idx+1])
        # Find Nearest Neighbor
        #  argmin() finds the minimum distance between the glacier lat/lon and the GCM pixel; .values is used to 
        #  extract the position's value as opposed to having an array
        lat_nearidx = (np.abs(main_glac_rgi[self.rgi_lat_colname].values[:,np.newaxis] - 
                              data.variables[self.lat_vn][:].values).argmin(axis=1))
        lon_nearidx = (np.abs(main_glac_rgi[self.rgi_lon_colname].values[:,np.newaxis] - 
                              data.variables[self.lon_vn][:].values).argmin(axis=1))
        # Find unique latitude/longitudes
        latlon_nearidx = list(zip(lat_nearidx, lon_nearidx))
        latlon_nearidx_unique = list(set(latlon_nearidx))
        # Create dictionary of time series for each unique latitude/longitude
        glac_variable_dict = {}
        for latlon in latlon_nearidx_unique:
            glac_variable_dict[latlon] = data[vn][start_idx:end_idx+1, latlon[0], latlon[1]].values
        # Convert to series
        glac_variable_series = np.array([glac_variable_dict[x] for x in latlon_nearidx])    
                
        # Perform corrections to the data if necessary
        # Surface air temperature corrections
        if vn in ['tas', 't2m', 'T2']:
            if 'units' in data[vn].attrs and data[vn].attrs['units'] == 'K':
                # Convert from K to deg C
                glac_variable_series = glac_variable_series - 273.15
            else:
                print('Check units of air temperature from GCM is degrees C.')
        elif vn in ['t2m_std']:
            if 'units' in data[vn].attrs and data[vn].attrs['units'] not in ['C', 'K']:
                print('Check units of air temperature standard deviation from GCM is degrees C or K')
        # Precipitation corrections
        # If the variable is precipitation
        elif vn in ['pr', 'tp', 'TOTPRECIP']:
            # If the variable has units and those units are meters (ERA Interim)
            if 'units' in data[vn].attrs and data[vn].attrs['units'] == 'm':
                pass
            # Elseif the variable has units and those units are kg m-2 s-1 (CMIP5)
            elif 'units' in data[vn].attrs and data[vn].attrs['units'] == 'kg m-2 s-1':  
                # Convert from kg m-2 s-1 to m day-1
                glac_variable_series = glac_variable_series/1000*3600*24
                #   (1 kg m-2 s-1) * (1 m3/1000 kg) * (3600 s / hr) * (24 hr / day) = (m day-1)
            # Elseif the variable has units and those units are mm (COAWST)
            elif 'units' in data[vn].attrs and data[vn].attrs['units'] == 'mm':
                glac_variable_series = glac_variable_series/1000
            # Else check the variables units
            else:
                print('Check units of precipitation from GCM is meters per day.')
            if self.timestep == 'monthly' and self.name != 'COAWST':
                # Convert from meters per day to meters per month (COAWST data already 'monthly accumulated precipitation')
                if 'daysinmonth' in dates_table.columns:
                    glac_variable_series = glac_variable_series * dates_table['daysinmonth'].values[np.newaxis,:]
        elif vn != self.lr_vn:
            print('Check units of air temperature or precipitation')
        return glac_variable_series, time_series
    

#%% Testing
if __name__ == '__main__':
##    gcm = GCM(name='CanESM2', rcp_scenario='rcp85')
    gcm = GCM(name='ERA5')
##    gcm = GCM(name='ERA-Interim')
#    
#    main_glac_rgi = modelsetup.selectglaciersrgitable(rgi_regionsO1=pygem_prms.rgi_regionsO1, rgi_regionsO2 = 'all',
#                                                      rgi_glac_number=pygem_prms.rgi_glac_number)
#    dates_table = modelsetup.datesmodelrun(startyear=1980, endyear=2017, spinupyears=0, 
#                                           option_wateryear=pygem_prms.gcm_wateryear)
#
#    # Air temperature [degC], Precipitation [m], Elevation [masl], Lapse rate [K m-1]
#    gcm_temp, gcm_dates = gcm.importGCMvarnearestneighbor_xarray(gcm.temp_fn, gcm.temp_vn, main_glac_rgi, dates_table)
#    gcm_prec, gcm_dates = gcm.importGCMvarnearestneighbor_xarray(gcm.prec_fn, gcm.prec_vn, main_glac_rgi, dates_table)
#    gcm_elev = gcm.importGCMfxnearestneighbor_xarray(gcm.elev_fn, gcm.elev_vn, main_glac_rgi)
#    if gcm.name == 'ERA-Interim' or gcm.name == 'ERA5':
#        gcm_lr, gcm_dates = gcm.importGCMvarnearestneighbor_xarray(gcm.lr_fn, gcm.lr_vn, main_glac_rgi, dates_table)
#    if gcm.name == 'ERA5':
#        gcm_tempstd, gcm_dates = gcm.importGCMvarnearestneighbor_xarray(gcm.tempstd_fn, gcm.tempstd_vn, main_glac_rgi, 
#                                                                        dates_table)
#    else:
#        gcm_lr = np.tile(ref_lr_monthly_avg, int(gcm_temp.shape[1]/12))