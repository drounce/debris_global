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
     Use wget and run from the oggm_project directory (i.e., change directory to oggm_project)
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
	  (or from external directly: python ERA5_preprocess_looplatlon.py -process_unique_latlon_data=1 -roi='12' -fromexternal='1')

3. Mass Balance data:
   - process MB data:

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

2. meltmodel_global.py: 
    ---> run model to get melt/Ts/snow at every timestep

3. meltcurves.py: 
    ---> processes meltmodel_global.py output to develop Ostrem curves for each glacier

4. ts_datetime_stats.ipynb: pull stats of surface temperature composite image for each lat/lon
    ---> environment: debris_thickness_global

5. tscurves.py: processes meltmodel_global.py output to develop Ts-hd curves for each lat/lon


6. debris_thickness.ipynb: 
    ---> invert ostrem curves to estimate the binned debris thickness
    ---> calibrate surface temperature inversion with sub-debris melt inversion


# ===== SCRIPTS ACTUALLY USED (for clean-up of the repository) =====
  - debris_stats.ipynb
  - ERA5_preprocess.py
  - ERA5_preprocess_looplatlon.py
  - process_mb_bin.ipynb
  - globaldebris_input.py
  - meltcurves.py
  - meltmodel_global.py
  - spc_run_meltmodel.sh
  - spc_split_lists.py
  - tscurves.py


# ===== Validation =========================================================================================================================
* = melt rates measurements
^ = debris thickness measurements

Alaska (01):
'01.15645' # Kennicott*^       - (Anderson etal 2019a, debris thickness and melt rates [Fig 7])

Western Canada (02):
'02.12438' # Dome              - (Mattson 2000, debris up to 0.5 m thick along the margins - no measurements)
'02.14297' # Emmons*^          - (Moore etal 2019, ablation measurements July 31 - Aug 10 2014 [Table 1, Figure 4]])
'02.18792' # Eliot	       - (Lundstrom 1993, thickness and ablation measurements - MAP with some values - TOO OLD)
'' Sierra Nevada region        - (Clark etal 1994)

Iceland (06):
'06.00474' # Svinafellsjokull* - (Moller etal 2016, melt rates from May 17-30 2013, but only can compare melt factor curves)

Svalbard (07):
'7.01044' # Larsbreen*^        - (Nicholson and Benn 2006, melt rates from July 9-20 2002 from 0.0005 to 0.10 m)
				 (Lukas etal 2005, hd all 1.33 +/- 0.72 m, 12 measurements)
'7.01107' # Longyearbreen^     - (Lukas etal 2005, hd all 1.84 +/- 1.32 m, 10 measurements)
'7.01106' # Nordenskioldtoppenbreen^ - (Lukas etal 2005, hd all 0.38 +/- 0.34 m, 32 measurements; NOT IN SCHERLER)

North Asia (10):
'10.01732' # Maliy Aktru       - (Mayer etal 2011 TCD discussion, 2007 ablation stakes, but only DDFs shown, could compare the melt factor curves)
			         (10 measurements available at ablation stakes with location on map)

Europe (11):
'11.00106' # Pasterze^	       - (Kellerer-Pirklbauer et al. 2008, thickness measurements: "Mean debris thickness along a longitudinal profile follows a power function decreasing from 				  47.3 cm close to the glacier terminus to 7.5 cm 3.8 km up-glacier"; SLL - mean 15 cm (2300 masl), VPL - mean 40 cm (2140 masl))
'11.00719' # Vernagtferner*    - (Juen etal 2013, melt rates June 25 - July 10 2010 [Fig 3])
'11.01604' # Suldenferner*^    - (del Gobbo 2017, 0.32 +/- 0.10)
'11.02472' # Venerocolo*^      - (Bocchiola etal 2015, melt rates 8/10/2007 - 9/13/2007 for various debris thicknesses)
'11.02810' # Haut d'Arolla*^   - (Reid et al. 2012; debris thicknesses [0.05 +/- 0.12 m, Table 4 and Fig 3] and ablation stakes [Fig 6])
'11.02858' # Belvedere*        - (Nicholson and Benn 2006 stake data)
'11.03005' # Miage*^           - (Foster et al. 2012, Table 2 [E: 0.23 +/- 0.16 m at 2030 mass, C: 0.32 +/- 0.13 m at 2060 masl], Figure 8) 
                   	         (Others: Mihalcea et al. 2008; Ablation stake data from Reid and Brock 2010, good fit)

Caucasus and Middle East (12):
'12.01012' # Zopkhito*         - (Lambrecht etal 2011, ablation stakes - data not in paper, could compare the melt factor curves)
'12.01132' # Djankuat*^        - (Popovnin and Rozova 2002, maps of debris thickness at various elevation bins from 1994 [Table 1])
		 	         (Lambrecht etal 2011, ablation stakes - data not in paper only DDF, could compare the melt factor curves)

HMA (13,14,15):
'13.05000' # S Inylcheck*      - good (Hagg etal 2008, melt rates from July 30 - August 10[Fig 2])
'13.43165' # No. 72*^  	       - good (Wang etal 2011, ablation rates and debris thickness; Wang et al. 2017)
'13.43174' # No. 74^	       - good (Wang etal 2011, debris thickness)
'13.43207' # Tuomuer	       - good (Wang etal 2011, debris thickness)
'13.43232' # Koxkar*^	       - good (Juen etal 2014, better map of When and Shiyin 2012, no ablation stake only degree-day factors, so compare could be melt factors; 
				      (Han etal 2006, ablation rate beneath 3 stakes from heat conduction method)
				      (Han etal 2006, debris thickness at 3 sites further up
				      (When and Shiyin 2012, debris thickness at terminus - see Juen etal 2014 for better representation)
				      (Chengwei etal 2007, hd vs elev)

'14.04477' # Hispar 	       - good
'14.06794' # Baltoro*^ 	       - good (Mihalcea et al. 2006, melt rate [Fig 7, Table 1, July 4-14 2004])
		 		      (Groos etal 2017, melt rate July 22 - Aug 10 2011 [Table 3, Fig 4])
				      (Minora etal 2015, melt rates SAME DATA AS GROOS July 22 - Aug 10 2011 and debris thickness with elevation, Table 2)
'14.16042' # Batal*^           - good (Patel etal 2016, debris thickness data [fig with altitude dependence; melt rate too - Table 1 good fit)
'14.15447' # Bara Shigri       - good (Schauwecker etal 2015; no measurements?--> our analysis suggest debris is thicker)
'14.20154' # Rakhiot	       - (Mattson etal 1989 - June 22 - August 8 1986, 10 stakes 0 - 0.4 m)
'15.03473' # Ngozumpa^	       - good (Nicholson and Benn 2012, Nicholson and Mertes 2017, Nicholson et al 2018)
'15.03733' # Khumbu* 	       - good (Gades etal 2000 referencing Nakawo et al 1986 "less than 0.1 m below the icefall to more than 2 m near the terminus")
			 	      (Kayastha etal 2000, melt rates - Table 1, excellent fit)
				      (Soncini etal 2016 up to 3.0 m with a mean value of 0.35 m (n = 64))
				      (Soncini etal 2016 have ablation stake data as well)
				      (Buizza 2014 - ablation rates, but not really usable, very thin debris)
'15.03734' # Changri Nup       - good
'15.04045' # Lirung*^	       - good (McCarthy etal 2017 - thickness, Chand etal 2015 melt rates seasonally) 
				      (Dahal 2015 Thesis  [April 6-16, 2004; Table 3.2 for heights Fig 4.11 - m w.e.d])
				      (Gades et al. (2000) "0.5 m below the Rockwell to 3 m near the terminus")
				      (Ragettli et al. 2015, debris thickness > 0.5 "too thick" to measure in 93% of the measurements)
'15.04121' # Langtang	       - good
'15.07886' # Hailuogou*^       - good (Zhang etal 2011, debris thickness data [fig 2] and melt data [fig 5c] underestimated)
				      (Zhang et al. 2016 - Hailougou glacier and a couple others (Dagongba and Xiaogongba) measured in 1982) 
				       debris thicknesses at various elevations [Fig 4] from Li and Su (1996, in Chinese))
'15.11758' # 24K  	       - Yang etal 2010 debris thickness and stake data [fig 2]


South America (17):
'17.13720' # Piramide	    - (Ayala etal 2016): use QGIS and map to get the distances roughly, then center over the DCG

New Zealand (18):
'18.02397' # Franz Josef*^  - (Brook etal 2013, debris thickness [Fig 4] and melt rates from February 7-16 2012 [Table 1])
'18.02375' # Fox*	    - (Brook etal 2012, ablation from 11/23/2007 - 12/03/20017 from 11 stakes, Figure 4a)
'18.02342' # Tasman	    - (Kirkbride 1989, 150 measurements of debris over lower 10 km [fig 2.18]; melt rate estimated from field data [Table 
			       2.20, fig 2.17, fig 2.2 shows stake locations])
			      (Kirkbride and Warren 1999, 15 debris thickness measurements (and estimated melt rate) from 1986)
			      (Purdie and Fitzharris 1999, one estimate from conduction method and some clean ice from 1995)



Fujita and Ageta 2000 - Xiao Dongkemadi glacier from 1992-1993


Mattson etal 1989 other ostrem curves older
Lundstrom metal 1993 - western North America



Automatically digitize figures by clicking on points (very fast):
   https://apps.automeris.io/wpd/

"Sample Raster Values" for direct comparison between bins and rasters

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
- Validation figure
  a) using the full x-axis needed to show the full ranges and in  that case make  in  inlet figure with 0-5 m so that information does not get too  squeezed.
  b) indicate the number of samples per dot by using different symbol size. Would need to be the same  symbol  for  comparability  in which case you would need  to use colors for the different sites.

  --> may be better to only compare over regions that are known to avoid errors with moraines increasing the debris thickness (or we just mention this)
  --> shall we add uncertainty to the model results based on the curves?



Limitations/issues:
  - 13.00604 - debris cover is up-glacier (this is a rock crop, not a DCG!): showcases issues with debris outlines and glacier outlines
  - Jan Mayen doesn't have surface temperature data from GEE for some reason
  - Composite Ts image, so two different dates comes into play
     --> Focused studies should apply method over known clear sky image
  - glacial lakes show up as very thin debris due to cool temperatures
  - Errors with debris cover prevent application to known debris-covered regions (e.g., Svalbard comparing Scherler etal 2018 with Lukas etal 2005 - N'breen is not even identified as debris-covered)
  - Errors with debris cover extent, e.g., Suldenferner prevent robust comparison with del Gobbo (2017)
      Suldenferner (RGI60-11.01604): clear case of poor glacier boundary
         --> good reason not to clip the debris thickness maps
      Hailougou (RGI60-15.07886): clear case of rock outcrops in accumulation being debris
         --> hence whey calibration is performed "without jumps"
  - Errors associated with different DEM sources
  - Errors associated with different time periods (elevation may have changed or the debris thickness too)
  - Most debris thicknesses from center of glacier, while binned data includes entire width including moraine which are typically thicker


Discussion:
  - additional benefit of temperature inversion over the sub-debris melt inversion method is that the temperature inversion can have debris thickness < 0.02,
    while this is an indeterminate problem for the sub-debris melt inversion method since there are two correct answers 
    - fortunately, this doesn't affect us because the terminus bin inversions are typically thickness > 0.02 m


# ===== Ideas =====
print('IDEA FOR PROCESSING:\n  set a larger minimum DC area for processing (remove discontinuous DC sections' + 
      '  --> then when clipping the debris thickness maps, clip using the original Scherler extents...')
print('  right now there are still way too many small bins messing up binned signals')
print('BETTER: SELECT ONLY THE BIGGEST POLYGON AND DISCARD THE REST?')




Test extrapolation method by comparing Ngozumpa and some other well-known glaciers (assuming they didn't have data...)



DEM citations
-  the correct citations for each source is there: https://github.com/OGGM/oggm/blob/master/oggm/data/dem_sources.txt
-  What we used for you is:
     - SRTM everywhere where available
     - ARCTICDEM and REMA for high latitudes
     - As a fallback for ARCTICDEM and REMA (when more that 5% of DEM data is missing on glacier) we use RAMP and GIMP
     - If everything above fails we fallback on DEM3
     - For ALASKA we use the product you gave me (the one from Christian Kienholz)


===== FOR NEXT ROUND ======
- use gaussian filter instead of mean filter (idea from David Shean)
- use ITS-Live grids to mask bins that have highly uncertain velocities!  This might avoid the terminus criteria.
- Thorsten has data for Peru and Bolivia from 2000-2016 that may be better than Braun et al. (2019) from 2000-2013.



===== Debris cover extent difference =====
- Fix geometries
- difference (input: RGI, overlay: new outlines)
- difference (input: scherler, overlay: difference of RGI new outlines)
    --> this is the debris cover extent that is similar for both


# ===== FIGURE 3 =====
- the region boundaries still pop out better on the SROCC figure. They don’t merge to something coherent without thinking hard.
7) Perhaps give the following at least a try: put a horizontal bar (perhaps only covering half the horizontal distance) in the red bar at value 0.56. This make help a lot to extract visually a bit more quantitative information.


# ===== JUSTIFICATON FOR OGGM RES =====
There are so many datasets I don't think we can list all native resolutions here. The Farinotti ice thickness product has its own arbitrary resolution based on glacier size (10, 20, 50, 100) in the idea very similar to what OGGM does (small glacier -> hig resolution) but not exact same of course.


# ===== 
Debris cover extent:
Similarly, we can use the binned dh/dt to identify areas where sub-debris melt is less than that associated with 3m (which we discard for this very reason!) as a way to estimate the uncertainty due to glacier retreat.  However, this could only be done for glaciers with data.

Additional analyses:
Regional debris thickness as a function of elevation (interquartile range can be shown in grey)
Maps of the spatial distribution of the calibrated parameters (a,b,c): are they similar across regions??
issue: it’ll be difficult to plot all three at the same time, so likely need 3 subplots
could assess spatial variability by performing autocorrelation analysis (semivariogram)

