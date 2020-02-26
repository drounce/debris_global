Debris-covered Glacier Energy Balance Model from Rounce et al. (2015)


Credits: If using the energy balance model, please cite the following:
Rounce, D.R., Quincey, D.J., and McKinney, D.C. (2015). Debris-covered glacier energy balance model for Imja-Lhotse Shar Glacier in the Everest region of Nepal, The Cryosphere, 9:2295-2310, doi:10.5194/tc-9-2295-2015.


Overview: the meltmodel.py script runs the code from Rounce et al. (2015), which is a debris-covered glacier energy balance model that requires input concerning the debris properties and meteorological data. Model output includes sub-debris melt rates, debris temperature, and energy fluxes.


Model details: the meltmodel.py script runs the code from Rounce et al. (2015) with a few modifications as discussed below:
- The code has been converted from Matlab to Python
- The code uses the LE_Rain option from Rounce et al. (2015), although other options (LE_RH100 and LE_Dry can easily be added to the code if precipitation data is unavailable).
- There are now two options for handling snow accumulation: (1) assuming that the snow depth is prescribed by the automatic weather station, so snow melt does not need to be computed and instead the net radiation, latent heat flux, and sensible heat flux for the debris is set to zero, and (2) explicitly accounting for snow accumulation and melt based on a modified version of Tarboten and Luce (1996) described below.
- A snow model has been added based on a modified version of Tarboten and Luce (1996) to compute fluxes with snow on the surface. Specifically, the incoming shortwave radiation is computed using the scheme already in the melt model as opposed to using their corrections for local slope and illumination. The ground heat flux is computed from the debris as opposed to using their estimates from diurnal soil temperature. The model uses a set temperature threshold for snowfall, which can be specified in the input, instead of their threshold. We do not calculate the heat flux from precipitation for rain on snow because it alters the cold content of the snow pack and therefore we do not want to double count this energy. We do not account for wind redistribution, since it is site specific. We do not iterate to solve for snow temperature, but instead simply calculate a depth average. The thermal conductivity at the debris/ice interface is estimated using the depth of snow and debris height. To avoid cooling the snow pack excessively (e.g, when the snow depth is < 1 cm and the net energy is highly negative), the change in snow pack temperature for a given time step is limited to 1 degree. Any remaining energy is used to cool the debris beneath the snow pack; however, this energy is also limited to cool the upper layer by 1 degree, which is necessary because the turbulent heat fluxes are estimated by the snow, yet there is no feedback between the turbulent heat fluxes and net energy as the debris cools thus the debris would cool to unrealistically low temperatures (< -100 degrees C). 

- Two bugs were fixed: (1) dealing with the specification of layers in the internal debris structure and (2) in counting the number of iterations for which to cut off the Newton-Raphson method if agreement cannot be reached. Both bugs were found to cause minor/no changes in the results (i.e., the total melt over an ablation season in Nepal for various debris thicknesses was changed by less than 1 cm).
- The model was significantly cleaned up during the transition and large portions of model code were converted to functions in order to make the code easier to read and easier to modify in the future.


# ===== TO-DO LIST / FUTURE WORK ===========================================================================================================
- add back in user-specified options for LE_RH100 and LE_Dry
- add back in user-specified options for instability corrections
- add option for wind speed height to vary over time (at present if it varies, the wind speed height is assumed to be the average over the entire study period)
- test/validate snow model (as of August 2019, the snow model has only been tested on Miage Glacier to ensure snow melt is "reasonable" during ablation season)
- test/validate model performance during transition and winter seasons (as of August 2019, the model has only been run during ablation seasons and into transition seasons)


# ===== MODEL WORKFLOW =====================================================================================================================

# ----- INPUT DATA -----
0. DEMs:
   - download OGGM's best DEM with ice thickness and widths:
	 https://cluster.klima.uni-bremen.de/~fmaussion/output/debris_project/
     Use wget:
       wget -r -nH --cut-dirs=2 -np -R "index.html*" https://cluster.klima.uni-bremen.de/~fmaussion/output/debris_project/

1. Debris cover extents:
   - process Scherler debris cover extents in QGIS:
     a. fix shape file
     b. multipart to single part
     c. open attribute -> new field 'area_single'
     d. select features area_single > min threshold (20000 m2): removes isolated pixels (< 2 surface temperature pixels)
     e. dissolve using RGIId
     f. open attribute -> new fields 'DC_Area_v2' and 'DC_Area_v2_%' - make sure multiply by 100 and convert km2 to m2 when divided by area
     - remove holes (if desirable)

2. ERA5 data
   - process ERA5 data: 
       a. python ERA5_preprocess.py -process_era5_hrly_data=1 -debug=1
            set latlon_list_raw = None in globaldebris_input.py if want to download before computing the latlon
       b. debris_stats.ipynb
	    run this until get to first unique lat/lons such that you can process the relevant data
       c. python ERA5_preprocess_looplatlon.py -process_unique_latlon_data=1
            download and pre-process into netcdf files for each site with a glacier

3. Mass Balance data:
   - process MB data: (this is done as part of the workflow now)

4. Surface temperature data:
   - Kraaijenbrink2017-ts-hma Google Earth Engine: surface temperature, etc.
       --> Modifications: export year, day of year and time; 
			  added buffer which is needed to ensure full coverage of glaciers by debris thickness estimates

   - ts_datetime_stats.ipynb: calculate information associated with the composite surface temperature raster from GEE
       --> debris_thickness_global environment

5. Velocity data:
   - ITS-Live: High Mountain Asia, Alaska, Arctic, part of South America
   - GoLive: Europe, New Zealand, Caucasus, part of South America, continental US
     --> ftp://dtn.rc.colorado.edu/work/nsidc0710/nsidc0710_landsat8_golive_ice_velocity_v1.1/p193_r028/


# ----- WORKFLOW -----
0. process_mb_bin.ipynb:
    ---> environment: debris_thickness_global environment
    ---> processes mass balance data, debris-covered areas, emergence velocities, etc.

1. debris_stats.ipynb:
    ---> environment: debris_thickness_global environment
    ---> identifies lat/lon of all glaciers with data (THIS SHOULD BE RUN FIRST SO CAN SELECT INDIVIDUAL LAT/LON)
    ---> adds debris cover elevation statistics for each lat/lon to the ERA5 dataset

      --> (now in debris_stats) HMA_emergence_velocity.ipynb: estimate binned emergence velocities

2. meltmodel_global.py: 
    ---> run model to get melt/Ts/snow at every timestep




3. meltcurves.py: 
    ---> processes meltmodel_global.py output to develop Ostrem curves for each glacier

4. ts_datetime_stats.ipynb: pull stats of surface temperature composite image for each lat/lon
    ---> environment: debris_thickness_global

5. tscurves.py: processes meltmodel_global.py output to develop Ts-hd curves for each lat/lon

6. global_melt2thickness.py: sub-debris melt inversion for debris thickness based on mass balance data and emergence velocities
    --> pygem_v2 environment because debris_thickness_global has no array
   #### MOVED THIS INTO HMA DEBRIS THICKNESS #####


7. HMA_debris_thickness.ipynb: calibrate surface temperature inversion with sub-debris melt inversion


# ===== SCRIPTS ACTUALLY USED (for clean-up of the repository) =====
  - debris_stats.ipynb
  - ERA5_preprocess.py
  - ERA5_preprocess_looplatlon.py
  - process_mb_bin.ipynb
  - globaldebris_input.py
  - meltmodel_global.py
  - spc_run_meltmodel.sh
  - spc_split_lists.py


# ===== Validation =========================================================================================================================
Alaska (01):
'01.15645' # Kennicott  - good

Western Canada (02):
'02.12438' # Dome Glacier -    Mattson (2000) - 

Europe (11):
'11.01604' # Suldenferner (Lindsey has datasets)
'11.02810' # Haut d'Arolla (Carenzo etal 2016)
'11.03005' # Miage (Mihalcea et al. 2008; Foster et al. 2012; Ablation stake data from Reid and Brock 2010)

HMA (13,14,15):
'13.43232' # Koxkar	- good (positive first bin)
'14.04477' # Hispar 	- good
'14.06794' # Baltoro 	- good
'14.16042' # Batal      - good (Patel et al. 2016 - with altitude dependence)
'14.15447' # Bara Shigri- good (Schauwecker etal 2015; no measurements?--> our analysis suggest debris is thicker)
'15.03473' # Ngozumpa	- good
'15.03733' # Khumbu 	- good (Gades et al 2000 referencing Nakawo et al 1986 "less than 0.1 m below the icefall to more than 2 m near the terminus"
'15.03734' # Changri Nup- good
'15.04045' # Lirung	- good (McCarthy et al. 2017) Gades et al. (2000) "0.5 m below the Rockwell to 3 m near the terminus"
'15.04121' # Langtang	- good
'15.07886' # Hailuogou  - HORRIBLE (Zhang et al. 2011, Figure 2)


# ===== UNCERTAINTIES ====================================================================================================================
- Model parameter uncertainty
   --> select glaciers in each region to come up with uncertainty associated with model parameters

- Mass balance data uncertainty 
    TO-DO: compare debris thickness estimates for select glaciers derived with different data sources
      - Alaska (McNabb vs. Braun)
      - S America (Braun vs. Berthier)
    TO-DO: add uncertainty with emergence velocity into this analysis

Are debris thickness estimates within the error bars of one another when derived from different sources?
  - Are the main discrepancies within the accumulation error due to the handling of the penetration corrections? 
      If so, then good for debris estimates.

Is uncertainty associated with the model parameters parameters comparable to the uncertainty associated with the mass balance data?


# ===== TO-DO LIST ========================================================================================================================

- PROCESS ERA5 DATA FOR ALL OTHER REGIONS

- SEE NOTES FROM .IPYNB FILES
  - convert .ipynb 'debris_stats' and 'HMA_debris_thickness' to .py scripts

- Ts for each region from GEE

- 15.07886: Hailuogou Glacier is HORRIBLE - the reason is emissivity - now fixed
  --> still the velocity is super low likely because the glacier is very
  --> add max terminus cells (e.g., lowest 10% based on elevation?)

  --> 10% of the elevation as opposed to the jump cells


-LIMIT MELTFACTOR BASED ON 2 CM


- Validation figure
  a) using the full x-axis needed to show the full ranges and in  that case make  in  inlet figure with 0-5 m so that information does not get too  squeezed.
  b) indicate the number of samples per dot by using different symbol size. Would need to be the same  symbol  for  comparability  in which case you would need  to use colors for the different sites.

