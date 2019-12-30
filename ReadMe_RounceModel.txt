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


Future work:
- ADD IN ABILITY TO MODEL CLEAN ICE MELT!
- add back in user-specified options for LE_RH100 and LE_Dry
- add back in user-specified options for instability corrections
- add option for wind speed height to vary over time (at present if it varies, the wind speed height is assumed to be the average over the entire study period)
- test/validate snow model (as of August 2019, the snow model has only been tested on Miage Glacier to ensure snow melt is "reasonable" during ablation season)
- test/validate model performance during transition and winter seasons (as of August 2019, the model has only been run during ablation seasons and into transition seasons)


Model workflow:
- process ERA5 orography data to elevation, which will be used to get lat/lon in future scripts
- debris_elevstats: identify lat/lon of all glaciers
- process ERA5 data: download and pre-process into netcdf files for each site with a glacier
- meltmodel_global.py: run model to get melt/Ts at every tilmestep
- meltcurves.py: processes meltmodel_global.py output to develop Ostrem curves for each glacier

# ===== Validation ======
'15.03473' # Ngozumpa	- good
'15.03733' # Khumbu 		- good
'15.03734' # Changri Nup	- good
'15.04121' # Langtang	- good
'15.04045' # Lirung
'14.06794' # Baltoro 	- good
'14.04477' # Hispar 		- good
'13.43232' # Koxkar		- good (positive first bin)


