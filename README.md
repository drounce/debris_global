# debris_global
Code to calculate the debris thickness of debris-covered glaciers globally using remote sensing data.

This project is currently under active development and under review.  Please contact David Rounce (drounce@cmu.edu) for any questions.

## CREDITS: 
If using the model/methods for estimating debris thickness, please contact David Rounce for the appropriate citation, which will soon replace:
<br> - Rounce, D.R., Hock, r., McNabb, R.W., Millan, R., Sommer, C., Braun, M.H., Malz, P., Maussion, F., Mouginot, J., Seehaus, T.C., and Shean, D.E. (in review). Distributed global debris thickness estimates reveal debris significantly impacts glacier mass balance.

If an emphasis is also placed on the debris-covered glacier energy balance model, please cite the following:
<br> - Rounce, D.R., Quincey, D.J., and McKinney, D.C. (2015). Debris-covered glacier energy balance model for Imja-Lhotse Shar Glacier in the Everest region of Nepal, The Cryosphere, 9:2295-2310, doi:10.5194/tc-9-2295-2015.

## OVERVIEW:
The methods and model in this repository have been used in the studies above.  Unfortunately due to time constraints, the model and methods are in various formats including python scripts and jupyter notebooks.  I have attempted to document the code extensively within those scripts such that they are easy to use and follow.  The following is a brief overview of the most relevant elements and scripts to get started.

### Melt Model
*meltmodel_global.py* and *globaldebris_input.py* can be used to run the code from Rounce et al. (2015), which is a debris-covered glacier energy balance model that requires input concerning the debris properties and meteorological data. Model output includes sub-debris melt rates, debris temperature, and energy fluxes.  *globaldebris_input.py* contains all the input for the debris thickness inversion and therefore all fields do not need to be filled to run just the melt model.  I have attemped create sections within this file such that this is easy to navigate. 
<br><br>There are also some important distinctions/updates from Rounce et al. (2015):
   - The code has been converted from Matlab to Python
   - The code uses the LE_Rain option from Rounce et al. (2015), although other options (LE_RH100 and LE_Dry can easily be added to the code if precipitation data is unavailable).
   - There are now two options for handling snow accumulation: (1) assuming that the snow depth is prescribed by the automatic weather station, so snow melt does not need to be computed and instead the net radiation, latent heat flux, and sensible heat flux for the debris is set to zero when snow is on the surface, and (2) explicitly accounting for snow accumulation and melt based on a modified version of Tarboten and Luce (1996) described below.
   - A *snow model* has been added based on a modified version of Tarboten and Luce (1996) to compute fluxes with snow on the surface. Specifically, the incoming shortwave radiation is computed using the scheme already in the melt model as opposed to using their corrections for local slope and illumination. The ground heat flux is computed from the debris as opposed to using their estimates from diurnal soil temperature. The model uses a set temperature threshold for snowfall, which can be specified in the input, instead of their threshold. We do not calculate the heat flux from precipitation for rain on snow because it alters the cold content of the snow pack and therefore we do not want to double count this energy. We do not account for wind redistribution, since it is site specific. We do not iterate to solve for snow temperature, but instead simply calculate a depth average. The thermal conductivity at the debris/ice interface is estimated using the depth of snow and debris height. To avoid cooling the snow pack excessively (e.g, when the snow depth is < 1 cm and the net energy is highly negative), the change in snow pack temperature for a given time step is limited to 1 degree. Any remaining energy is used to cool the debris beneath the snow pack; however, this energy is also limited to cool the upper layer by 1 degree, which is necessary because the turbulent heat fluxes are estimated by the snow, yet there is no feedback between the turbulent heat fluxes and net energy as the debris cools thus the debris would cool to unrealistically low temperatures (< -100 degrees C). 
   - Two *bugs* were fixed: (1) dealing with the specification of layers in the internal debris structure and (2) in counting the number of iterations for which to cut off the Newton-Raphson method if agreement cannot be reached. Both bugs were found to cause minor/no changes in the results (i.e., the total melt over an ablation season in Nepal for various debris thicknesses was changed by less than 1 cm).
   - The model was significantly cleaned up during the transition and large portions of model code were converted to functions in order to make the code easier to read and easier to modify in the future.
   
### Python Environment
All the code is run using the *debris_thickness_global* environment, which can be installed on your computer using the .yml file with anaconda.  For instructions, please use see the following: https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html.  

<br> To run the scripts in the code, clone this repository and run on your local computer.  The workflow described below takes you through the general steps for estimating the debris thickness.

### Input Data
- DEMs
  - download OGGM's best DEM with ice thickness and widths (https://oggm.org/)

- Debris cover extents:
  - process Scherler et al. (2018) debris cover extents in QGIS:
  <br> (a) fix shape file
  <br> (b) multipart to single part
  <br> (c) open attribute -> new field 'area_single'
  <br> (d) select features area_single > min threshold (20000 m2): removes isolated pixels (< 2 surface temperature pixels)
  <br> (e) dissolve using RGIId
  <br> (f) open attribute -> new fields 'DC_Area_v2' and 'DC_Area_v2_%', make sure multiply by 100 and convert km2 to m2 when divided by area.
  - remove holes (if desired)

- ERA5 data
  - process ERA5 data: 
  <br> (a) $ *python ERA5_preprocess.py -process_era5_hrly_data=1 -debug=1*
  <br> &nbsp;&nbsp;&nbsp;&nbsp; set latlon_list_raw = None in *globaldebris_input.py* if want to download before computing the latlon
  <br> (b) *debris_stats.ipynb* 
  <br> &nbsp;&nbsp;&nbsp;&nbsp; run until get to first unique lat/lons (i.e., only part of the notebook) such that you can process the relevant data
  <br> (c) $ *python ERA5_preprocess_looplatlon.py -process_unique_latlon_data=1*
  <br> &nbsp;&nbsp;&nbsp;&nbsp; download and pre-process into netcdf files for each site with a glacier
	<br> &nbsp;&nbsp;&nbsp;&nbsp; (or from external directly: python ERA5_preprocess_looplatlon.py -process_unique_latlon_data=1 -roi='12' -fromexternal='1')

- Mass Balance (elevation change) data:
  - the workflow is set up to use .tif files of elevation change rates

- Surface temperature data:
  - Kraaijenbrink2017-ts-hma Google Earth Engine: surface temperature, etc.
  <br> Modifications include exporting the year, day of year and time.  Also, added buffer which is needed to ensure full coverage of glaciers by debris thickness estimates
  - *ts_datetime_stats.ipynb*
  <br> used to calculate information associated with the composite surface temperature raster from GEE

- Velocity data:
  - ITS-Live: High Mountain Asia, Alaska, Arctic, part of South America
  - Other regions were provided by Romain Millan (see study)
 
- Debris thickness and sub-debris melt data:
  - Used for validation of the modeled sub-debris melt rates and debris thickness estimates
  - Many studies do not provide the data, so the following tool was used to automatically create .csv files from figures: https://apps.automeris.io/wpd/


### Debris Thickness Workflow
- *process_mb_bin.ipynb*
  - processes mass balance data, debris-covered areas, emergence velocities, etc.

- *debris_elev_stats.ipynb*
  - identifies lat/lon of all glaciers with data (this should be run in pre-processing input data, so you can select individual lat/lons)
  - adds debris cover elevation statistics for each lat/lon to the ERA5 dataset

- *meltmodel_global.py*
  - run model to get melt/Ts/snow at every timestep

- *meltcurves.py*
  - processes *meltmodel_global.py* output to develop Ostrem curves for each glacier

- *ts_datetime_stats.ipynb*
  - pull stats of surface temperature composite image for each lat/lon

- *tscurves.py*
  - processes meltmodel_global.py output to develop Ts-hd curves for each lat/lon

- *debris_thickness-calibration.ipynb*
  - invert ostrem curves to estimate the binned debris thickness
  - calibrate surface temperature inversion with sub-debris melt inversion

- *debris_thickness-extrapolation.ipynb*
  - extrapolates the debris thickness from nearest neighbors

Note that some debris-covered glaciers may be missing due to a lack of surface temperature data or other issues.  We have processed some of these sites as an additional step, but did not distribute them via NSIDC to avoid individuals using them without thinking of the quality control.  If you would like a specific debris-covered glacier that is not provided, please contact David Rounce.

### Directory / File Description
Below is a brief description of each directory and file.

<br> Directories:
- *debrisglobal/*: directory containing the *glacfeat.py* and *globaldebris_input.py*.  Ideally, this directory structure will continue to be developed to enable pip install in the future.
- *old_scripts/*: directory containing old scripts that were used in model development that I've kept around just in case there are useful elements.

<br> Scripts and jupyter notebooks:
- *bin_hd_add_uncertainty.ipynb*: jupyter notebook that contains code to estimate the uncertainty associated with the debris thickness and enhancement factors for each elevation bin.  The output is useful for dealing with uncertainty in large-scale glacier evolution model simulations.
- *class_climate_debris.py*: python script containing the GCM class that is used to load ERA5 data into the energy balance model.
- *debris_elev_stats.ipynb*: jupyter notebook used to compute elevation statistics for the debris-covered areas in each latitude and longitude.
- *debris_thickness-calibration.ipynb*: jupyter notebook used to compute the debris thickness using sub-debris and temperature inversion methods.
- *debris_thickness-extrapolation.ipynb*: jupyter notebook used to extrapolate the debris thickness from the nearest calibrated glaciers
- *debris_thickness-extrapolation-crossvalidation.ipynb*: jupyter notebook used to perform a leave-one-out cross validation for the extrapolation of calibrated glaciers.
- *debris_thickness-extrapolation-missing.ipynb*: jupyter notebook used to extrapolate the debris thickness for any glaciers that were initially missing (note this may require adjusting the input data, e.g., if surface temperature data was not available).
- *ERA5_preprocess.py*: python script to (i) subset ERA5 hourly data into regional files to make the file sizes manageable or (ii) extract the data for a given lat/lon that includes all the required meteorological data for the debris-covered glacier energy balance model.
- *ERA5_preprocess_looplatlon.py*: python script to loop through and process all available hourly datasets or just those associated with a specific lat/lon.
- *hd_pygem_mbdif.ipynb*: jupyter notebook used to compute the impact accounting for debris has on regional mass balance after the simulations have been run.
- *hd_resampling_wreg_stats.ipynb*: jupyter notebook used to resample the debris thickness to a common resolution and calculate the regional statistics for the debris thickness and enhancement factors.
- *hd_resampling_wreg_stats.ipynb*: jupyter notebook used above but ignoring any melt enhancement in the statistics (at the request of a reviewer).
- *hd_uncertainty.py*: python script with a few options to create figures such as the debris thickness vs. melt uncertainty figure, the schematic of how debris thickness uncertainty was calculated, and a methods diagram figure.
- *hd_validation.py*: python script with various options to compare modeled and observed sub-debris melt rates or debris thicknesses at both the bin and point scale. This was used to create some of the figures in the study.
- *hist_emvel_calibrated_glaciers.ipynb*: jupyter notebook containing histograms and scatter plots of the emergence velocities of the calibrated glaciers.
- *mb_gradient_model*: modification of the mass balance gradient model used by Kraaijenbrink et al. (2017) to estimate the impact of debris on regional mass balance. Unfortunately, due to the apparent assumptions used in this model (e.g., the uppermost bin has a mass balance equal to the total precipitation, i.e., there's no melt in that bin), we did not use it in this study.  However, it's potentially useful to have as reference in the future.
- *meltcurves.py*: python script used to process the debris-covered glacier energy balance output into debris thickness versus sub-debris melt curves that were then used for the sub-debris melt inversion methods.
- *meltmodel_global.py*: python script used to run the debris-covered glacier energy balance model
- *nsidc_process_mb_bins.ipynb*: jupyter notebook used to process the calibrated mass balance bins into a suitable .csv file for sharing via NSIDC.
- *plots_calibrated_prms_hist.ipynb*: jupyter notebook used to plot histograms of the calibrated model parameters in response to a reviewer comment.
- *plots_extrapolation-nearestneighbor.ipynb*: jupyter notebook containing histograms and calculations of the distance to the nearest neighbor for each extrapolated glacier.
- *process_Anderson_hd_data*: jupyter notebook used to process the debris thickness data provided by Leif Anderson (https://zenodo.org/record/4317470#.X-TlbOlKhTa) into a suitable format for comparison
- *process_mb_bin.ipynb*: jupyter notebook used to process all of the debris-covered glaciers' datasets into elevation bins including mass balance, surface temperature, and emergence velocities
- *process_mb_bin-all.ipynb*: jupyter notebook used to process all of the elevation change datasets into elevation bins
- *ts_datetime_stats.ipynb*: jupyter notebook used to process all of the surface temperature data and pull out the date and time statistics associated with a given lat/lon, i.e, the date/time you'd want to pull the corresonding modeled surface temperatures from the debris-covered glacier energy balance.
- *ts_preprocess.ipynb*: jupyter notebook used to determine the filename associated with the surface temperature data for each glacier to assist the processing scripts.
- *ts_preprocess-missing.ipynb*: jupyter notebook used to preprocess the missing surface temperature data for glaciers that didn't have data in the first round.
- *tscurves.py*: jupyter notebook used to process the surface temperature output from the debris-covered glacier energy balance model into a format used by the surface temperature inversion methods.

<br> Bash scripts for UAF's supercomputer:
- *spc_run_meltmodel.sh*: bash script used to run the melt model on UAF's supercomputer.
- *spc_run_ostremcurves.sh*: bash script used to process the ostrem curves (i.e., *meltcurves.py*) on UAF's supercomputer.
- *spc_run_tsdata.sh*: bash script used to process the modeled surface temperature data (i.e., *tscurves.py*) on UAF's supercomputer.
- *spc_split_lists.py*: python script used to split lists for processing on UAF's supercomputer.
