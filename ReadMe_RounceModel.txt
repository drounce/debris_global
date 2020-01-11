Debris-covered Glacier Energy Balance Model from Rounce et al. (2015)


Credits: If using the energy balance model, please cite the following:
Rounce, D.R., Quincey, D.J., and McKinney, D.C. (2015). Debris-covered glacier energy balance model for Imja-Lhotse Shar Glacier in the Everest region of Nepal, The Cryosphere, 9:2295-2310, doi:10.5194/tc-9-2295-2015.


Overview: the meltmodel.py script runs the code from Rounce et al. (2015), which is a debris-covered glacier energy balance model that requires input concerning the debris properties and meteorological data. Model output includes sub-debris melt rates, debris temperature, and energy fluxes.


Model details: the meltmodel.py script runs the code from Rounce et al. (2015) with a few modifications as discussed below:
- The code has been converted from Matlab to Python
- The code uses the LE_Rain option from Rounce et al. (2015), although other options (LE_RH100 and LE_Dry can easily be added to the code if precipitation data is unavailable).
- There are now two options for handling snow accumulation: (1) assuming that the snow depth is prescribed by the automatic weather station, so snow melt does not need to be computed and instead the net radiation, latent heat flux, and sensible heat flux for the debris is set to zero, and (2) explicitly accounting for snow accumulation and melt based on a modified version of Tarboten and Luce (1996) described below.
- A snow model has been added based on a modified version of Tarboten and Luce (1996) to compute fluxes with snow ison the surface. Specifically, the incoming shortwave radiation is computed using the scheme already in the melt model as opposed to using their corrections for local slope and illumination. The ground heat flux is computed from the debris as opposed to using their estimates from diurnal soil temperature. The model uses a set temperature threshold for snowfall, which can be specified in the input, instead of their threshold. We do not calculate the heat flux from precipitation for rain on snow because it alters the cold content of the snow pack and therefore we do not want to double count this energy. We do not account for wind redistribution, since it is site specific. We do not iterate to solve for snow temperature, but instead simply calculate a depth average. The thermal conductivity at the debris/ice interface is estimated using the depth of snow and debris height.
- Two bugs were fixed: (1) dealing with the specification of layers in the internal debris structure and (2) in counting the number of iterations for which to cut off the Newton-Raphson method if agreement cannot be reached. Both bugs were found to cause minor/no changes in the results (i.e., the total melt over an ablation season in Nepal for various debris thicknesses was changed by less than 1 cm).
- The model was significantly cleaned up during the transition and large portions of model code were converted to functions in order to make the code easier to read and easier to modify in the future.


# ===== TO-DO LIST / FUTURE WORK ===========================================================================================================
- MULTI-GLACIER CALIBRATION - all glaciers that share ERA5 data should be calibrated together
- AREA-WEIGHTED BINS FOR GLACIER CALIBRATION & MULTI-GLACIER CALIBRATION!
  --> this is important to ensure we place most weight on glaciers/bins that have the most data
- IMPROVE FIT BASED ON VARIATIONS IN THE STANDARD DEVIATION AS OPPOSED TO THE TOFFSET

- consider how to use hd_max: should these bins be avoided during calibration, since debris thickness estimates are very sensitive to them and  they are the most uncertain?


- ADD IN ABILITY TO MODEL CLEAN ICE MELT!
- ADD SNOW AS AN OUTPUT to be able to discern whether there is snow on the surface and mask these values from the Ts-debris curve fitting
- add back in user-specified options for LE_RH100 and LE_Dry
- add back in user-specified options for instability corrections
- add option for wind speed height to vary over time (at present if it varies, the wind speed height is assumed to be the average over the entire study period)
- test/validate snow model (as of August 2019, the snow model has only been tested on Miage Glacier to ensure snow melt is "reasonable" during ablation season)
- test/validate model performance during transition and winter seasons (as of August 2019, the model has only been run during ablation seasons and into transition seasons)


# ===== MODEL WORKFLOW =====================================================================================================================
- process Scherler debris cover extents in QGIS:
  - fix shape file
  - remove isolated pixels (< 2 surface temperature pixels):
     - multipart to single part
     - open attribute -> new field 'area_single'
     - select features area_single > min threshold (20000 m2)
     - dissolve using RGIId
     - open attribute -> new fields 'DC_Area_v2' and 'DC_Area_v2_%' - make sure multiply by 100 and convert km2 to m2 when divided by area
  - remove holes (if desirable)

- process MB data:
  - Shean sent binned data
  - Braun processed using some PyGEM scripts; need to clean up

- process ERA5 orography data to elevation, which will be used to get lat/lon in future scripts

- process ERA5 data: download and pre-process into netcdf files for each site with a glacier

- debris_stats.ipynb:
   ---> environment: debris_thickness_global environment
   ---> identifies lat/lon of all glaciers with data
   ---> processes emergence velocity, debris cover area, and mass balance data

   - HMA: done
   - Alaska: done

- meltmodel_global.py: run model to get melt/Ts/snow at every timestep

- HMA_emergence_velocity.ipynb: estimate binned emergence velocities
    

- meltcurves.py: processes meltmodel_global.py output to develop Ostrem curves for each glacier

- global_melt2thickness.py: sub-debris melt inversion for debris thickness based on mass balance data and emergence velocities
    --> pygem_v2 environment because debris_thickness_global has no array

- Kraaijenbrink2017-ts-hma Google Earth Engine: surface temperature, etc.
    --> Modifications: export year, day of year and time; 
			  added buffer which is needed to ensure full coverage of glaciers by debris thickness estimates

- ts_datetime_stats.ipynb: calculate information associated with the composite surface temperature raster from GEE
    --> debris_thickness_global environment


# ===== Validation =========================================================================================================================
'15.03473' # Ngozumpa	- good
'15.03733' # Khumbu 		- good
'15.03734' # Changri Nup	- good
'15.04121' # Langtang	- good
'15.04045' # Lirung
'14.06794' # Baltoro 	- good
'14.04477' # Hispar 		- good
'13.43232' # Koxkar		- good (positive first bin)

# ===== UNCERTAINTY WITH MASS BALANCE DATA ================================================================================================
- Some glaciers have data from Braun and Larsen: can use to quantify uncertainty due to data source!


# ===== HMA MISSING DATA IN MELT2THICKNESS ================================================================================================
Missing Ostrem data:
 ['2850N-9000E', '2850N-9050E', '2900N-9325E', '3050N-8150E', '3150N-8025E', '3250N-8800E', '3350N-8650E', '3375N-8250E', '3400N-7975E', '3425N-8050E', '3425N-8650E', '3450N-8025E', '3450N-8150E', '3450N-8975E', '3525N-8225E', '3525N-8250E', '3550N-8025E', '3550N-8200E', '3550N-8250E', '3550N-9475E', '3575N-7800E', '3575N-7925E', '3575N-8100E', '3575N-8200E', '3575N-8550E', '3575N-8575E', '3575N-8850E', '3575N-8875E', '3600N-7725E', '3600N-7825E', '3600N-7875E', '3600N-7900E', '3625N-7800E', '3625N-7825E', '3625N-8725E', '3625N-8750E', '3650N-8700E', '3775N-8925E', '3825N-7075E', '4525N-8050E', '4525N-8075E', '4525N-8100E'] 

Glaciers without Ostrem data:
 ['13.00062', '13.04185', '13.04832', '13.04837', '13.05185', '13.05197', '13.05455', '13.05465', '13.05479', '13.05480', '13.05481', '13.20715', '13.26255', '13.26432', '13.33582', '13.34649', '13.34658', '13.36206', '13.36211', '13.36213', '13.36966', '13.36969', '13.36973', '13.36977', '13.36982', '13.36983', '13.36986', '13.37061', '13.37079', '13.37086', '13.37494', '13.37496', '13.38210', '13.38220', '13.38289', '13.38307', '13.38313', '13.38315', '13.38316', '13.39135', '13.39476', '13.39623', '13.39693', '13.39702', '13.39705', '13.39722', '13.39814', '13.39819', '13.39864', '13.39894', '13.39900', '13.39901', '13.39942', '13.39947', '13.39948', '13.39959', '13.39977', '13.39981', '13.39988', '13.39995', '13.40006', '13.40030', '13.40068', '13.40075', '13.49017', '13.49041', '13.49043', '13.49044', '13.49045', '13.49046', '13.49052', '13.49061', '13.49095', '13.49113', '13.49974', '13.51409', '13.51418', '13.51437', '13.51540', '13.51692', '13.51851', '13.52478', '13.53312', '13.53334', '13.53578', '13.53643', '13.53645', '13.53649', '13.53651', '13.53655', '13.53660', '13.53662', '13.53800', '13.53876', '13.53945', '13.53946', '13.53947', '13.53951', '13.53953', '13.54051', '13.54052', '13.54059', '13.54067', '13.54068', '13.54072', '13.54074', '13.54075', '13.54077', '13.54086', '13.54090', '13.54145', '13.54270', '13.54339', '13.54345', '13.54355', '14.06570', '14.26244', '14.27801'] 

