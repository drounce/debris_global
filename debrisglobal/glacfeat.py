#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 08:41:09 2020

@author: davidrounce
"""
import os
import re
from collections import OrderedDict

import geopandas as gpd
import numpy as np
import pandas as pd
from osgeo import gdal, ogr, osr
from scipy import ndimage

from pygeotools.lib import malib, warplib, geolib, iolib, timelib

#import globaldebris_input as input
import debrisglobal.globaldebris_input as debris_prms


class GlacFeat:
    def __init__(self, feat, glacname_fieldname, glacnum_fieldname):

        self.glacname = feat.GetField(glacname_fieldname)
        if self.glacname is None:
            self.glacname = ""
        else:
            #RGI has some nonstandard characters
            #self.glacname = self.glacname.decode('unicode_escape').encode('ascii','ignore')
            #glacname = re.sub(r'[^\x00-\x7f]',r'', glacname)
            self.glacname = re.sub(r'\W+', '', self.glacname)
            self.glacname = self.glacname.replace(" ", "")
            self.glacname = self.glacname.replace("_", "")
            self.glacname = self.glacname.replace("/", "")

        self.glacnum = feat.GetField(glacnum_fieldname)
#        fn = feat.GetDefnRef().GetName()
        #RGIId (String) = RGI50-01.00004
        self.glacnum = '%0.5f' % float(self.glacnum.split('-')[-1])

        if self.glacname:
            self.feat_fn = "%s_%s" % (self.glacnum, self.glacname)
        else:
            self.feat_fn = str(self.glacnum)

        self.glac_geom_orig = geolib.geom_dup(feat.GetGeometryRef())
        self.glac_geom = geolib.geom_dup(self.glac_geom_orig)
        #Hack to deal with fact that this is not preserved in geom when loaded from pickle on disk
        self.glac_geom_srs_wkt = self.glac_geom.GetSpatialReference().ExportToWkt()

        #Attributes written by mb_calc
        self.z1 = None
        self.z1_hs = None
        self.z1_stats = None
        self.z1_ela = None
        self.z2 = None
        self.z2_hs = None
        self.z2_stats = None
        self.z2_ela = None
        self.z2_aspect = None
        self.z2_aspect_stats = None
        self.z2_slope = None
        self.z2_slope_stats = None
        self.res = None
        self.dhdt = None
        self.dc_dhdt = None
        self.mb = None
        self.dc_mb = None
        self.mb_mean = None
        self.t1 = None
        self.t2 = None
        self.dt = None
        self.t1_mean = None
        self.t2_mean = None
        self.dt_mean = None

        self.H = None
        self.H_mean = np.nan
        self.vx = None
        self.vy = None
        self.vm = None
        self.vm_mean = np.nan
        self.divQ = None
        self.emvel = None
#         self.debris_class = None
        self.debris_thick = None
        self.debris_thick_mean = np.nan
        self.perc_clean = np.nan
        self.perc_debris = np.nan
        self.perc_pond = np.nan

    def geom_srs_update(self, srs=None):
        if self.glac_geom.GetSpatialReference() is None:
            if srs is None:
                srs = osr.SpatialReference()
                srs.ImportFromWkt(self.glac_geom_srs_wkt)
            self.glac_geom.AssignSpatialReference(srs)

    def geom_attributes(self, srs=None):
        self.geom_srs_update()
        if srs is not None:
            #Should reproject here to equal area, before geom_attributes
            #self.glac_geom.AssignSpatialReference(glac_shp_srs)
            #self.glac_geom_local = geolib.geom2localortho(self.glac_geom)
            geolib.geom_transform(self.glac_geom, srs)

        self.glac_geom_extent = geolib.geom_extent(self.glac_geom)
        self.glac_area = self.glac_geom.GetArea()
        self.glac_area_km2 = self.glac_area / 1E6
        self.cx, self.cy = self.glac_geom.Centroid().GetPoint_2D()
        

    #%%
    def emergence_pixels(self, vel_x_raw, vel_y_raw, icethickness_raw, xres, yres, 
                         vel_min=0, max_velocity=600, vel_depth_avg_factor=0.8, option_border=1,
                         positive_is_east=True, positive_is_north=True, constant_icethickness=False, debug=True):
        """ Compute the emergence velocity using an ice flux approach
        """
        
        # Glacier mask
        glac_mask = np.zeros(vel_x_raw.shape) + 1
        glac_mask[self.z1.mask] = 0
        
        # Replace nan with 0
        vel_x_raw = np.nan_to_num(vel_x_raw,0)
        vel_y_raw = np.nan_to_num(vel_y_raw,0)
        
        # Modify vel_y by multiplying velocity by -1 such that matrix operations agree with flow direction
        #    Specifically, a negative y velocity means the pixel is flowing south.
        #    However, if you were to subtract that value from the rows, it would head north in the matrix.
        #    This is due to the fact that the number of rows start at 0 at the top.
        #    Therefore, multipylying by -1 aligns the matrix operations with the flow direction
        if positive_is_north:
            vel_y = -1*vel_y_raw * vel_depth_avg_factor
        else:
            vel_y = vel_y_raw * vel_depth_avg_factor
        if positive_is_east:
            vel_x = vel_x_raw * vel_depth_avg_factor
        else:
            vel_x = -1*vel_x_raw * vel_depth_avg_factor
        vel_total = (vel_y**2 + vel_x**2)**0.5
        # Ice thickness
        icethickness = icethickness_raw.copy()
        if constant_icethickness:
            icethickness[:,:] = 1
            icethickness = icethickness * glac_mask
    #     print('mean ice thickness:', np.round(icethickness.mean(),0), 'm')
        # Compute the initial volume
        volume_initial = icethickness * (xres * yres)
        pix_maxres = xres
        if yres > pix_maxres:
            pix_maxres = yres
        # Quality control options:
        # Apply a border based on the max specified velocity to prevent errors associated with pixels going out of bounds
        if option_border == 1:
            border = int(max_velocity / pix_maxres) + 1
            for r in range(vel_x.shape[0]):
                for c in range(vel_x.shape[1]):
                    if (r < border) | (r >= vel_x.shape[0] - border) | (c < border) | (c >= vel_x.shape[1] - border):
                        vel_x[r,c] = 0
                        vel_y[r,c] = 0
        # Minimum/maximum velocity bounds
        vel_x[vel_total < vel_min] = 0
        vel_y[vel_total < vel_min] = 0
        vel_x[vel_total > max_velocity] = 0
        vel_y[vel_total > max_velocity] = 0
    #     # Remove clusters of high velocity on stagnant portions of glaciers due to feature tracking of cliffs and ponds
    #     if option_stagnantbands == 1:
    #         vel_x[bands <= stagnant_band] = 0
    #         vel_y[bands <= stagnant_band] = 0        
        # Compute displacement in units of pixels
        vel_x_pix = vel_x / xres
        vel_y_pix = vel_y / yres
        # Compute the displacement and fraction of pixels moved for all columns (x-axis)
        # col_x1 is the number of columns to the closest pixel receiving ice [ex. 2.6 returns 2, -2.6 returns -2]
        #    int() automatically rounds towards zero
        col_x1 = vel_x_pix.astype(int)
        # col_x2 is the number of columns to the further pixel receiving ice [ex. 2.6 returns 3, -2.6 returns -3]
        #    np.sign() returns a value of 1 or -1, so it's adding 1 pixel away from zero
        col_x2 = (vel_x_pix + np.sign(vel_x_pix)).astype(int)
        # rem_x2 is the fraction of the pixel that remains in the further pixel (col_x2) 
        #    [ex. 2.6 returns 0.6, -2.6 returns 0.6]
        #    np.sign() returns a value of 1 or -1, so multiplying by that ensures you have a positive value
        #    then when you take the remainder using "% 1", you obtain the desired fraction
        rem_x2 = np.multiply(np.sign(vel_x_pix), vel_x_pix) % 1
        # rem_x1 is the fraction of the pixel that remains in the closer pixel (col_x1) 
        #    [ex. 2.6 returns 0.4, -2.6 returns 0.4]
        rem_x1 = 1 - rem_x2
        # Repeat the displacement and fraction computations for all rows (y-axis)
        row_y1 = vel_y_pix.astype(int)
        row_y2 = (vel_y_pix + np.sign(vel_y_pix)).astype(int)
        rem_y2 = np.multiply(np.sign(vel_y_pix), vel_y_pix) % 1
        rem_y1 = 1 - rem_y2
              
        # Compute the mass flux for each pixel
        volume_final = np.zeros(volume_initial.shape)
        for r in range(vel_x.shape[0]):
            for c in range(vel_x.shape[1]):
                volume_final[r+row_y1[r,c], c+col_x1[r,c]] = (
                    volume_final[r+row_y1[r,c], c+col_x1[r,c]] + rem_y1[r,c]*rem_x1[r,c]*volume_initial[r,c]
                    )
                volume_final[r+row_y2[r,c], c+col_x1[r,c]] = (
                    volume_final[r+row_y2[r,c], c+col_x1[r,c]] + rem_y2[r,c]*rem_x1[r,c]*volume_initial[r,c]
                    )
                volume_final[r+row_y1[r,c], c+col_x2[r,c]] = (
                    volume_final[r+row_y1[r,c], c+col_x2[r,c]] + rem_y1[r,c]*rem_x2[r,c]*volume_initial[r,c]
                    )
                volume_final[r+row_y2[r,c], c+col_x2[r,c]] = (
                    volume_final[r+row_y2[r,c], c+col_x2[r,c]] + rem_y2[r,c]*rem_x2[r,c]*volume_initial[r,c]
                    )
             
        # Redistribute off-glacier volume back onto the nearest pixel on the glacier
        offglac_row, offglac_col = np.where((glac_mask == 0) & (volume_final > 0))
        for nidx in range(0,len(offglac_row)):
            nrow = offglac_row[nidx]
            ncol = offglac_col[nidx]
            ridx, cidx = nearest_nonzero_idx(glac_mask, nrow, ncol)
            # Add off-glacier volume back onto nearest pixel on glacier
            volume_final[ridx,cidx] += volume_final[nrow,ncol]
            volume_final[nrow,ncol] = 0
                
        # Check that mass is conserved (threshold = 0.1 m x pixel_size**2)
        if debug:
            print('Mass is conserved?', np.absolute(volume_final.sum() - volume_initial.sum()) / 
                  volume_initial.sum() < 0.01)
            print('volume_final - volume initial:', 
                  np.round(np.absolute(volume_final.sum() - volume_initial.sum()),1), 'm3 or',
                  np.absolute(volume_final.sum() - volume_initial.sum()) / volume_initial.sum() * 100, '%')
            
        if np.absolute(volume_final.sum() - volume_initial.sum()) / volume_initial.sum() > 0.01:
            print('MASS NOT CONSERVED FOR EMERGENCE VELOCITY')
        # Final ice thickness
        icethickness_final = volume_final / (xres * yres)
        # Emergence velocity
        emergence_velocity = icethickness_final - icethickness
        return emergence_velocity
        
        
    #%%
#    def add_elev_data(self, thick_dir, thick_fn, verbose=False):
    def add_layers(self, dc_shp_lyr, gf_add_dhdt=True, gf_add_vel=True, gf_add_ts=True, 
                   gf_add_ts_info=False, gf_add_slope_aspect=False, calc_emergence=False,
                   verbose=False, debug_emergence=False):
        
        glac_str = self.glacnum
        region = glac_str.split('.')[0]
        rgiid = 'RGI60-' + self.glacnum
        
        # =====FILENAMES =====
        # Add the filenames
        fn_dict = OrderedDict()
        # DEM
        z1_fp = debris_prms.oggm_fp + 'dems/RGI60-' + str(region.zfill(2)) + '/'
        z1_fn = 'RGI60-' + str(region.zfill(2)) + '.' + rgiid.split('.')[1] + '_dem.tif'
        fn_dict['z1'] = z1_fp + z1_fn
        # Ice thickness
        thick_dir = debris_prms.oggm_fp + 'thickness/RGI60-' + str(region.zfill(2)) + '/'
        thick_fn = 'RGI60-' + str(region.zfill(2)) + '.' + rgiid.split('.')[1] + '_thickness.tif'
        fn_dict['ice_thick'] = thick_dir + thick_fn
        # dh/dt
        if gf_add_dhdt:
            dhdt_fn = debris_prms.dhdt_fn_dict[debris_prms.roi]
            if dhdt_fn is not None:
                fn_dict['dhdt'] = dhdt_fn
        # Velocity
        if gf_add_vel:
            if os.path.exists(debris_prms.v_dir + debris_prms.vx_fn_dict[debris_prms.roi]):
                fn_dict['vx'] = debris_prms.v_dir + debris_prms.vx_fn_dict[debris_prms.roi]
                fn_dict['vy'] = debris_prms.v_dir + debris_prms.vy_fn_dict[debris_prms.roi]
        # Surface temperature
        if gf_add_ts:
            if os.path.exists(debris_prms.ts_fp + debris_prms.ts_fn_dict[debris_prms.roi]):
                fn_dict['ts'] = debris_prms.ts_fp + debris_prms.ts_fn_dict[debris_prms.roi]

        if gf_add_ts_info:
            if os.path.exists(debris_prms.ts_fp + debris_prms.ts_dayfrac_fn_dict[debris_prms.roi]):
                fn_dict['ts_dayfrac'] = debris_prms.ts_fp + debris_prms.ts_dayfrac_fn_dict[debris_prms.roi]
            if os.path.exists(debris_prms.ts_fp + debris_prms.ts_year_fn_dict[debris_prms.roi]):
                fn_dict['ts_year'] = debris_prms.ts_fp + debris_prms.ts_year_fn_dict[debris_prms.roi]
            if os.path.exists(debris_prms.ts_fp + debris_prms.ts_doy_fn_dict[debris_prms.roi]):
                fn_dict['ts_doy'] = debris_prms.ts_fp + debris_prms.ts_doy_fn_dict[debris_prms.roi]


        # ===== PROCESS THE DATA =====
        #Expand extent to include buffered region around glacier polygon
        warp_extent = geolib.pad_extent(self.glac_geom_extent, width=debris_prms.buff_dist)
        if verbose:
            print("Expanding extent")
            print(self.glac_geom_extent)
            print(warp_extent)
            print(self.aea_srs)

        #Warp everything to common res/extent/proj
        z1_gt = gdal.Open(fn_dict['z1']).GetGeoTransform()
        z1_res = np.min([z1_gt[1], -z1_gt[5]])
        ds_list = warplib.memwarp_multi_fn(fn_dict.values(), res=z1_res, extent=warp_extent, 
                                           t_srs=self.aea_srs, verbose=verbose, r='cubic')
        ds_dict = dict(zip(fn_dict.keys(), ds_list))
        self.ds_dict = ds_dict

        if verbose:
            print(ds_list)
            print(fn_dict.keys())

        #Get global glacier mask
        #Want this to be True over ALL glacier surfaces, not just the current polygon
#        glac_shp_lyr_mask = geolib.lyr2mask(glac_shp_lyr, ds_dict['ice_thick'])
        dc_shp_lyr_mask = geolib.lyr2mask(dc_shp_lyr, ds_dict['ice_thick'])

        if 'z1' in ds_dict:
            #This is False over glacier polygon surface, True elsewhere - can be applied directly
            glac_geom_mask = geolib.geom2mask(self.glac_geom, ds_dict['z1'])
            self.z1 = np.ma.array(iolib.ds_getma(ds_dict['z1']), mask=glac_geom_mask)

            self.res = geolib.get_res(ds_dict['z1'])

            # Debris cover
            self.dc_mask = np.ma.mask_or(dc_shp_lyr_mask, glac_geom_mask)
            self.dc_area = np.ma.array(iolib.ds_getma(ds_dict['z1']), mask=self.dc_mask)

            if verbose:
                print('\n\n# z1 pixels:', self.z1.count(), '\n')

        # ===== ADD VARIOUS LAYERS TO gf =====
        if gf_add_slope_aspect:
            #Caluclate stats for aspect and slope
            self.z2_aspect = np.ma.array(geolib.gdaldem_mem_ds(ds_dict['z1'], processing='aspect', returnma=True), 
                                       mask=glac_geom_mask)
            self.z2_aspect_stats = malib.get_stats(self.z2_aspect)
            self.z2_slope = np.ma.array(geolib.gdaldem_mem_ds(ds_dict['z1'], processing='slope', returnma=True), 
                                      mask=glac_geom_mask)
            self.z2_slope_stats = malib.get_stats(self.z2_slope)

        # copy for Ts because it changes the mask otherwise and messes up binned statistics for whole glacier
        glac_geom_mask_copy = glac_geom_mask.copy()

        # ==== ADD LAYERS =====
        if 'ice_thick' in ds_dict:
            #Load ice thickness
            self.H = np.ma.array(iolib.ds_getma(ds_dict['ice_thick']), mask=glac_geom_mask)
            self.H_mean = self.H.mean()
            if verbose:
                print('mean ice thickness [m]:', self.H_mean)

        if 'vx' in ds_dict and 'vy' in ds_dict:
            #Load surface velocity maps
            self.vx = np.ma.array(iolib.ds_getma(ds_dict['vx']), mask=glac_geom_mask)
            self.vy = np.ma.array(iolib.ds_getma(ds_dict['vy']), mask=glac_geom_mask)
            self.vm = np.ma.sqrt(self.vx**2 + self.vy**2)
            self.vm_mean = self.vm.mean()
            if verbose:
                print('mean velocity [m/s]:', self.vm_mean)

        if 'ts_dayfrac' in ds_dict:
            #Load surface temperature maps
            self.ts_dayfrac = np.ma.array(iolib.ds_getma(ds_dict['ts_dayfrac']), mask=glac_geom_mask)
            self.ts_dayfrac.mask = np.ma.mask_or(
                    glac_geom_mask, np.ma.getmask(np.ma.masked_array(self.ts_dayfrac.data, 
                                                                     np.isnan(self.ts_dayfrac.data))))
        else:
            self.ts_dayfrac = None

        if 'ts_year' in ds_dict:
            #Load surface temperature maps
            self.ts_year = np.ma.array(iolib.ds_getma(ds_dict['ts_year']), mask=glac_geom_mask)
            self.ts_year.mask = np.ma.mask_or(
                    glac_geom_mask, np.ma.getmask(np.ma.masked_array(self.ts_year.data, np.isnan(self.ts_year.data))))
        else:
            self.ts_year = None

        if 'ts_doy' in ds_dict:
            #Load surface temperature maps
            self.ts_doy = np.ma.array(iolib.ds_getma(ds_dict['ts_doy']), mask=glac_geom_mask)
            self.ts_doy.mask = np.ma.mask_or(
                    glac_geom_mask, np.ma.getmask(np.ma.masked_array(self.ts_doy.data, np.isnan(self.ts_doy.data))))
        else:
            self.ts_doy = None

        # Emergence velocity
        if calc_emergence and 'vx' in ds_dict and 'vy' in ds_dict and self.H is not None:
            vx = np.ma.filled(self.vx,0)
            vy = np.ma.filled(self.vy,0)
            H = np.ma.filled(self.H,0)
            vx[self.z1 > self.z1.max()] = 0
            vy[self.z1 > self.z1.max()] = 0
            H[self.z1 > self.z1.max()] = 0
            vmax = np.nanmax((vx**2 + vy**2)**0.5)

            # Emergence computation
            emvel = self.emergence_pixels(vx, vy, H, self.res[0], self.res[1], 
                                          positive_is_east=True, positive_is_north=True, 
                                          constant_icethickness=False, max_velocity=vmax, vel_min=0, 
                                          debug=debug_emergence)
            # 3x3 filter to reduce
            if debris_prms.emvel_filter_pixsize > 0:
                emvel = ndimage.filters.convolve(
                    emvel, weights=np.full((debris_prms.emvel_filter_pixsize, debris_prms.emvel_filter_pixsize), 
                                            1.0/debris_prms.emvel_filter_pixsize**2))
            # Add to glacier feature
            self.emvel = np.ma.masked_array(emvel, mask=np.ma.getmask(self.z1))

        if 'ts' in ds_dict:
            #Load surface temperature maps
            self.ts = np.ma.array(iolib.ds_getma(ds_dict['ts']), mask=glac_geom_mask_copy)
            self.ts.mask = np.ma.mask_or(glac_geom_mask, 
                                       np.ma.getmask(np.ma.masked_array(self.ts.data, np.isnan(self.ts.data))))
            # Debris only
            self.dc_ts = self.ts.copy()
            self.dc_ts.mask = self.dc_mask
        else:
            self.ts = None
            self.dc_ts = None

        if 'dhdt' in ds_dict:
            self.dhdt = np.ma.array(iolib.ds_getma(ds_dict['dhdt']), mask=glac_geom_mask_copy)
            self.dhdt.mask = np.ma.mask_or(
                glac_geom_mask, np.ma.getmask(np.ma.masked_array(self.dhdt.data, np.isnan(self.dhdt.data))))
            self.mb = self.dhdt.copy() * debris_prms.density_ice / debris_prms.density_water

            # Debris only
            self.dc_dhdt = np.ma.array(iolib.ds_getma(ds_dict['dhdt']), mask=glac_geom_mask_copy)
            self.dc_dhdt.mask = self.dc_mask
            self.dc_mb = self.dc_dhdt.copy() * debris_prms.density_ice / debris_prms.density_water

        if 'debris_thick_ts' in ds_dict:
            # Load debris thickness map
            self.debris_thick_ts = np.ma.array(
                iolib.ds_getma(ds_dict['debris_thick_ts']), mask=glac_geom_mask_copy)
            self.meltfactor_ts = None
        else:
            self.debris_thick_ts = None
            self.meltfactor_ts = None

        if verbose:
            print('Area [km2]:', self.glac_area / 1e6)
            print('-------------------------------')
        
    
        
    #%%
    def hist_plot(self, bin_width=50.0, dz_clim=(-2.0, 2.0), exportcsv=False, csv_ending='', mb_df=None, 
                  outdir_csv=None):
        #print("Generating histograms")
        #Create bins for full range of input data and specified bin width
    
        #NOTE: these counts/areas are for valid pixels only
        #Not necessarily a true representation of actual glacier hypsometry
        #Need a void-filled DEM for this
        if mb_df is not None:
            # Align bins with mass balance data
            bin_center_min = mb_df.loc[0,'bin_center_elev_m']
            while bin_center_min > self.z1.min() + bin_width/2:
                bin_center_min -= debris_prms.mb_bin_size
            bin_center_max = mb_df['bin_center_elev_m'].values[-1]
            while bin_center_max < self.z1.max():
                bin_center_max += debris_prms.mb_bin_size    
            z_bin_centers = np.arange(bin_center_min, bin_center_max + debris_prms.mb_bin_size/2, debris_prms.mb_bin_size)
            z_bin_edges = np.arange(bin_center_min - debris_prms.mb_bin_size / 2, 
                                    debris_prms.bin_center_max + debris_prms.mb_bin_size, debris_prms.mb_bin_size)
        else:
            z_bin_edges, z_bin_centers = malib.get_bins(self.z1, bin_width)
            
        #Need to compress here, otherwise histogram uses masked values!
        z1_bin_counts, z1_bin_edges = np.histogram(self.z1.compressed(), bins=z_bin_edges)
        z1_bin_areas = z1_bin_counts * self.res[0] * self.res[1] / 1E6
        #RGI standard is integer thousandths of glaciers total area
        #Should check to make sure sum of bin areas equals total area
        #z1_bin_areas_perc = 100. * z1_bin_areas / np.sum(z1_bin_areas)
        z1_bin_areas_perc = 100. * (z1_bin_areas / self.glac_area_km2)
        z1_bin_areas_perc_cum = np.cumsum(z1_bin_areas) / z1_bin_areas.sum() * 100
    
        #If we only have one elevation grid with dhdt
        if self.z2 is not None:
            z2_bin_counts, z2_bin_edges = np.histogram(self.z2.compressed(), bins=z_bin_edges)
            z2_bin_areas = z2_bin_counts * self.res[0] * self.res[1] / 1E6
            #z2_bin_areas_perc = 100. * z2_bin_areas / np.sum(z2_bin_areas)
            z2_bin_areas_perc = 100. * (z1_bin_areas / self.glac_area_km2)
        else:
            z2_bin_counts = z1_bin_counts
#            z2_bin_edges = z1_bin_edges
            z2_bin_areas = z1_bin_areas
            z2_bin_areas_perc = z1_bin_areas_perc
            
        if self.dc_area is not None:
            dc_bin_counts, dc_bin_edges = np.histogram(self.dc_area.compressed(), bins=z_bin_edges)
            dc_bin_areas = dc_bin_counts * self.res[0] * self.res[1] / 1E6
            dc_bin_areas_perc = np.zeros(dc_bin_counts.shape)
            dc_bin_areas_perc[dc_bin_counts > 0] = (100. * (dc_bin_areas[dc_bin_counts > 0] / 
                                                            z1_bin_areas[dc_bin_counts > 0]))
    #         outbins_df['bin_debris_perc'] = outbins_df['dc_bin_count_valid'] / outbins_df['z1_bin_count_valid'] * 100
            # Cumulative debris cover
            dc_bin_area_cumsum = np.cumsum(dc_bin_areas)
            dc_bin_areas_perc_cum = dc_bin_area_cumsum / dc_bin_areas.sum() * 100
    
        #Create arrays to store output
        slope_bin_med = np.ma.masked_all_like(z1_bin_areas)
        slope_bin_mad = np.ma.masked_all_like(z1_bin_areas)
        aspect_bin_med = np.ma.masked_all_like(z1_bin_areas)
        aspect_bin_mad = np.ma.masked_all_like(z1_bin_areas)
        if self.dhdt is not None:
            mb_bin_med = np.ma.masked_all_like(z1_bin_areas)
            np.ma.set_fill_value(mb_bin_med, np.nan)
            mb_bin_mad = np.ma.masked_all_like(mb_bin_med)
            mb_bin_mean = np.ma.masked_all_like(mb_bin_med)
            mb_bin_std = np.ma.masked_all_like(mb_bin_med)
            dhdt_bin_med = np.ma.masked_all_like(mb_bin_med)
            dhdt_bin_mad = np.ma.masked_all_like(mb_bin_med)
            dhdt_bin_mean = np.ma.masked_all_like(mb_bin_med)
            dhdt_bin_std = np.ma.masked_all_like(mb_bin_med)
            dhdt_bin_count = np.ma.masked_all_like(mb_bin_med)
        if self.dc_dhdt is not None:
            dc_mb_bin_med = np.ma.masked_all_like(z1_bin_areas)
            np.ma.set_fill_value(dc_mb_bin_med, np.nan)
            dc_mb_bin_mad = np.ma.masked_all_like(dc_mb_bin_med)
            dc_mb_bin_mean = np.ma.masked_all_like(dc_mb_bin_med)
            dc_mb_bin_std = np.ma.masked_all_like(dc_mb_bin_med)
            dc_dhdt_bin_med = np.ma.masked_all_like(dc_mb_bin_med)
            dc_dhdt_bin_mad = np.ma.masked_all_like(dc_mb_bin_med)
            dc_dhdt_bin_mean = np.ma.masked_all_like(dc_mb_bin_med)
            dc_dhdt_bin_std = np.ma.masked_all_like(dc_mb_bin_med)
            dc_dhdt_bin_count = np.ma.masked_all_like(dc_mb_bin_med)
        if self.vm is not None:
            vm_bin_med = np.ma.masked_all_like(z1_bin_areas)
            vm_bin_mad = np.ma.masked_all_like(z1_bin_areas)
        if self.H is not None:
            H_bin_mean = np.ma.masked_all_like(z1_bin_areas)
            H_bin_std = np.ma.masked_all_like(z1_bin_areas)
        if self.emvel is not None:
            emvel_bin_mean = np.ma.masked_all_like(z1_bin_areas)
            emvel_bin_std = np.ma.masked_all_like(z1_bin_areas)
            emvel_bin_med = np.ma.masked_all_like(z1_bin_areas)
            emvel_bin_mad = np.ma.masked_all_like(z1_bin_areas)
    #     if self.debris_class is not None:
    #         debris_thick_med = np.ma.masked_all_like(z1_bin_areas)
    #         debris_thick_mad = np.ma.masked_all_like(z1_bin_areas)
    
        if self.ts is not None:
            ts_mean = np.ma.masked_all_like(z1_bin_areas)
            ts_std = np.ma.masked_all_like(z1_bin_areas)
            ts_med = np.ma.masked_all_like(z1_bin_areas)
            ts_mad = np.ma.masked_all_like(z1_bin_areas)
            dc_ts_mean = np.ma.masked_all_like(z1_bin_areas)
            dc_ts_std = np.ma.masked_all_like(z1_bin_areas)
            dc_ts_med = np.ma.masked_all_like(z1_bin_areas)
            dc_ts_mad = np.ma.masked_all_like(z1_bin_areas)
        if self.debris_thick_ts is not None:
            debris_thick_ts_mean = np.ma.masked_all_like(z1_bin_areas)
            debris_thick_ts_std = np.ma.masked_all_like(z1_bin_areas)
            debris_thick_ts_med = np.ma.masked_all_like(z1_bin_areas)
            debris_thick_ts_mad = np.ma.masked_all_like(z1_bin_areas)
        if self.meltfactor_ts is not None:
            meltfactor_ts_mean = np.ma.masked_all_like(z1_bin_areas)
            meltfactor_ts_std = np.ma.masked_all_like(z1_bin_areas)
            meltfactor_ts_med = np.ma.masked_all_like(z1_bin_areas)
            meltfactor_ts_mad = np.ma.masked_all_like(z1_bin_areas)
    
        #Bin sample count must be greater than this value
        min_bin_samp_count = debris_prms.min_bin_samp_count
    
        #Loop through each bin and extract stats
        idx = np.digitize(self.z1, z_bin_edges)
        for bin_n in range(z_bin_centers.size):
            if self.dhdt is not None:
                mb_bin_samp = self.mb[(idx == bin_n+1)]
                if mb_bin_samp.count() > min_bin_samp_count:
                    mb_bin_med[bin_n] = malib.fast_median(mb_bin_samp)
                    mb_bin_mad[bin_n] = malib.mad(mb_bin_samp)
                    mb_bin_mean[bin_n] = mb_bin_samp.mean()
                    mb_bin_std[bin_n] = mb_bin_samp.std()
                    
                dhdt_bin_samp = self.dhdt[(idx == bin_n+1)]
                if dhdt_bin_samp.count() > min_bin_samp_count:
                    dhdt_bin_med[bin_n] = malib.fast_median(dhdt_bin_samp)
                    dhdt_bin_mad[bin_n] = malib.mad(dhdt_bin_samp)
                    dhdt_bin_mean[bin_n] = dhdt_bin_samp.mean()
                    dhdt_bin_std[bin_n] = dhdt_bin_samp.std()
                    dhdt_bin_count[bin_n] = dhdt_bin_samp.count()
            
            if self.dc_dhdt is not None:
                dc_mb_bin_samp = self.dc_mb[(idx == bin_n+1)]
                if dc_mb_bin_samp.count() > min_bin_samp_count:
                    dc_mb_bin_med[bin_n] = malib.fast_median(dc_mb_bin_samp)
                    dc_mb_bin_mad[bin_n] = malib.mad(dc_mb_bin_samp)
                    dc_mb_bin_mean[bin_n] = dc_mb_bin_samp.mean()
                    dc_mb_bin_std[bin_n] = dc_mb_bin_samp.std()
                dc_dhdt_bin_samp = self.dc_dhdt[(idx == bin_n+1)]
                if dc_dhdt_bin_samp.count() > min_bin_samp_count:
                    dc_dhdt_bin_med[bin_n] = malib.fast_median(dc_dhdt_bin_samp)
                    dc_dhdt_bin_mad[bin_n] = malib.mad(dc_dhdt_bin_samp)
                    dc_dhdt_bin_mean[bin_n] = dc_dhdt_bin_samp.mean()
                    dc_dhdt_bin_std[bin_n] = dc_dhdt_bin_samp.std()
                    dc_dhdt_bin_count[bin_n] = dc_dhdt_bin_samp.count()
                    
    #        if self.debris_thick is not None:
    #            debris_thick_bin_samp = self.debris_thick[(idx == bin_n+1)]
    #            if debris_thick_bin_samp.size > min_bin_samp_count:
    #                debris_thick_med[bin_n] = malib.fast_median(debris_thick_bin_samp)
    #                debris_thick_mad[bin_n] = malib.mad(debris_thick_bin_samp)
            
            if self.ts is not None:
                ts_bin_samp = self.ts[(idx == bin_n+1)]
                if ts_bin_samp.size > min_bin_samp_count:
                    ts_mean[bin_n] = ts_bin_samp.mean()
                    ts_std[bin_n] = ts_bin_samp.std()
                    ts_med[bin_n] = malib.fast_median(ts_bin_samp)
                    ts_mad[bin_n] = malib.mad(ts_bin_samp)
                dc_ts_bin_samp = self.dc_ts[(idx == bin_n+1)]
                if dc_ts_bin_samp.size > min_bin_samp_count:
                    dc_ts_mean[bin_n] = dc_ts_bin_samp.mean()
                    dc_ts_std[bin_n] = dc_ts_bin_samp.std()
                    dc_ts_med[bin_n] = malib.fast_median(dc_ts_bin_samp)
                    dc_ts_mad[bin_n] = malib.mad(dc_ts_bin_samp)
            if self.debris_thick_ts is not None:
                debris_thick_ts_bin_samp = self.debris_thick_ts[(idx == bin_n+1)]
                if debris_thick_ts_bin_samp.size > min_bin_samp_count:
                    debris_thick_ts_mean[bin_n] = debris_thick_ts_bin_samp.mean()
                    debris_thick_ts_std[bin_n] = debris_thick_ts_bin_samp.std()
                    debris_thick_ts_med[bin_n] = malib.fast_median(debris_thick_ts_bin_samp)
                    debris_thick_ts_mad[bin_n] = malib.mad(debris_thick_ts_bin_samp)
            if self.meltfactor_ts is not None:
                meltfactor_ts_bin_samp = self.meltfactor_ts[(idx == bin_n+1)]
                if meltfactor_ts_bin_samp.size > min_bin_samp_count:
                    meltfactor_ts_mean[bin_n] = meltfactor_ts_bin_samp.mean()
                    meltfactor_ts_std[bin_n] = meltfactor_ts_bin_samp.std()
                    meltfactor_ts_med[bin_n] = malib.fast_median(meltfactor_ts_bin_samp)
                    meltfactor_ts_mad[bin_n] = malib.mad(meltfactor_ts_bin_samp)

    #         if self.debris_class is not None:
    #             debris_class_bin_samp = self.debris_class[(idx == bin_n+1)]
    #             dhdt_clean_bin_samp = self.dhdt_clean[(idx == bin_n+1)]
    #             dhdt_debris_bin_samp = self.dhdt_debris[(idx == bin_n+1)]
    #             dhdt_pond_bin_samp = self.dhdt_pond[(idx == bin_n+1)]
    #             if debris_class_bin_samp.count() > min_bin_samp_count:
    #                 perc_clean[bin_n] = 100. * (debris_class_bin_samp == 1).sum()/debris_class_bin_samp.count()
    #                 perc_debris[bin_n] = 100. * (debris_class_bin_samp == 2).sum()/debris_class_bin_samp.count()
    #                 perc_pond[bin_n] = 100. * (debris_class_bin_samp == 3).sum()/debris_class_bin_samp.count()
    #             if dhdt_clean_bin_samp.count() > min_bin_samp_count:
    #                 dhdt_clean_bin_med[bin_n] = malib.fast_median(dhdt_clean_bin_samp)
    #             if dhdt_debris_bin_samp.count() > min_bin_samp_count:
    #                 dhdt_debris_bin_med[bin_n] = malib.fast_median(dhdt_debris_bin_samp)
    #             if dhdt_pond_bin_samp.count() > min_bin_samp_count:
    #                 dhdt_pond_bin_med[bin_n] = malib.fast_median(dhdt_pond_bin_samp)
            if self.vm is not None:
                vm_bin_samp = self.vm[(idx == bin_n+1)]
                if vm_bin_samp.size > min_bin_samp_count:
                    vm_bin_med[bin_n] = malib.fast_median(vm_bin_samp)
                    vm_bin_mad[bin_n] = malib.mad(vm_bin_samp)
            if self.H is not None:
                H_bin_samp = self.H[(idx == bin_n+1)]
                if H_bin_samp.size > min_bin_samp_count:
                    H_bin_mean[bin_n] = H_bin_samp.mean()
                    H_bin_std[bin_n] = H_bin_samp.std()
            if self.emvel is not None:
                emvel_bin_samp = self.emvel[(idx == bin_n+1)]
                if emvel_bin_samp.size > min_bin_samp_count:
                    emvel_bin_mean[bin_n] = emvel_bin_samp.mean()
                    emvel_bin_std[bin_n] = emvel_bin_samp.std()
                    emvel_bin_med[bin_n] = malib.fast_median(emvel_bin_samp)
                    emvel_bin_mad[bin_n] = malib.mad(emvel_bin_samp)
            if self.z2_slope is not None:
                slope_bin_samp = self.z2_slope[(idx == bin_n+1)]
                if slope_bin_samp.size > min_bin_samp_count:
                    slope_bin_med[bin_n] = malib.fast_median(slope_bin_samp)
                    slope_bin_mad[bin_n] = malib.mad(slope_bin_samp)
            if self.z2_aspect is not None:
                aspect_bin_samp = self.z2_aspect[(idx == bin_n+1)]
                if aspect_bin_samp.size > min_bin_samp_count:
                    aspect_bin_med[bin_n] = malib.fast_median(aspect_bin_samp)
                    aspect_bin_mad[bin_n] = malib.mad(aspect_bin_samp)
    
        if self.dhdt is not None:
            dhdt_bin_areas = dhdt_bin_count * self.res[0] * self.res[1] / 1E6
            #dhdt_bin_areas_perc = 100. * dhdt_bin_areas / np.sum(dhdt_bin_areas)
            dhdt_bin_areas_perc = 100. * (dhdt_bin_areas / self.glac_area_km2)
    
        outbins_header = ('bin_center_elev_m,z1_bin_count_valid,z1_bin_area_valid_km2,z1_bin_area_perc,' + 
                          'z1_bin_areas_perc_cum,z2_bin_count_valid,z2_bin_area_valid_km2,z2_bin_area_perc,' + 
                          'slope_bin_med,aspect_bin_med')
        fmt = '%0.1f,%0.0f,%0.3f,%0.2f,%0.2f,%0.0f,%0.3f,%0.2f,%0.2f,%0.2f'
        outbins = [z_bin_centers, z1_bin_counts, z1_bin_areas, z1_bin_areas_perc, z1_bin_areas_perc_cum, 
                   z2_bin_counts, z2_bin_areas, z2_bin_areas_perc, slope_bin_med, aspect_bin_med]
        if self.dhdt is not None:
            outbins_header = (outbins_header + ',dhdt_bin_count,dhdt_bin_area_valid_km2,dhdt_bin_area_perc,' + 
                              'dhdt_bin_mean_ma,dhdt_bin_std_ma,dhdt_bin_med_ma,dhdt_bin_mad_ma,mb_bin_mean_mwea,' + 
                              'mb_bin_std_mwea,mb_bin_med_mwea,mb_bin_mad_mwea')
            fmt += ',%0.0f,%0.3f,%0.2f,%0.2f,%0.2f,%0.2f,%0.2f,%0.2f,%0.2f,%0.2f,%0.2f'
            outbins.extend([dhdt_bin_count, dhdt_bin_areas, dhdt_bin_areas_perc, dhdt_bin_mean, dhdt_bin_std,
                            dhdt_bin_med, dhdt_bin_mad, mb_bin_mean, mb_bin_std, mb_bin_med, mb_bin_mad, ])
        if self.dc_dhdt is not None:
            outbins_header = (outbins_header + ',dc_dhdt_bin_count,dc_dhdt_bin_mean_ma,dc_dhdt_bin_std_ma,' + 
                              'dc_dhdt_bin_med_ma,dc_dhdt_bin_mad_ma,dc_mb_bin_mean_mwea,dc_mb_bin_std_mwea,' + 
                              'dc_mb_bin_med_mwea,dc_mb_bin_mad_mwea')
            fmt += ',%0.0f,%0.2f,%0.2f,%0.2f,%0.2f,%0.2f,%0.2f,%0.2f,%0.2f'
            outbins.extend([dc_dhdt_bin_count, dc_dhdt_bin_mean, dc_dhdt_bin_std, dc_dhdt_bin_med, dc_dhdt_bin_mad,  
                            dc_mb_bin_mean, dc_mb_bin_std, dc_mb_bin_med, dc_mb_bin_mad])
        if self.dc_area is not None:
            outbins_header += ',dc_bin_count_valid,dc_bin_area_valid_km2,dc_bin_area_perc,dc_bin_area_perc_cum'
            fmt += ',%0.0f,%0.3f,%0.2f,%0.2f'
            outbins.extend([dc_bin_counts, dc_bin_areas, dc_bin_areas_perc, dc_bin_areas_perc_cum])
    #         outbins.extend([z1_bin_counts, z1_bin_areas, z1_bin_areas_perc])
            
            
    #    if self.debris_thick is not None:
    #        outbins_header += ',hd_med_m,hd_mad_m'
    #        fmt += ',%0.2f,%0.2f'
    #        debris_thick_med[debris_thick_med == -(np.inf)] = 0.00
    #        debris_thick_mad[debris_thick_mad == -(np.inf)] = 0.00
    #        outbins.extend([debris_thick_med, debris_thick_mad])
        if self.ts is not None:
            outbins_header += ',ts_mean,ts_std,ts_med,ts_mad'
            fmt += ',%0.2f,%0.2f,%0.2f,%0.2f'
            outbins.extend([ts_mean, ts_std, ts_med, ts_mad])
            # debris cover only
            outbins_header += ',dc_ts_mean,dc_ts_std,dc_ts_med,dc_ts_mad'
            fmt += ',%0.2f,%0.2f,%0.2f,%0.2f'
            outbins.extend([dc_ts_mean, dc_ts_std, dc_ts_med, dc_ts_mad])
        if self.debris_thick_ts is not None:
            outbins_header += ',hd_ts_mean_m,hd_ts_std_m,hd_ts_med_m,hd_ts_mad_m'
            fmt += ',%0.2f,%0.2f,%0.2f,%0.2f'
            debris_thick_ts_mean[debris_thick_ts_mean == -(np.inf)] = 0.00
            debris_thick_ts_std[debris_thick_ts_std == -(np.inf)] = 0.00
            debris_thick_ts_med[debris_thick_ts_med == -(np.inf)] = 0.00
            debris_thick_ts_mad[debris_thick_ts_mad == -(np.inf)] = 0.00
            outbins.extend([debris_thick_ts_mean, debris_thick_ts_std, debris_thick_ts_med, debris_thick_ts_mad])
        if self.meltfactor_ts is not None:
            outbins_header += ',mf_ts_mean,mf_ts_std,mf_ts_med,mf_ts_mad'
            fmt += ',%0.2f,%0.2f,%0.2f,%0.2f'
            meltfactor_ts_mean = np.ma.filled(meltfactor_ts_mean,1)
            meltfactor_ts_std = np.ma.filled(meltfactor_ts_std,0)
            meltfactor_ts_med = np.ma.filled(meltfactor_ts_med,1)
            meltfactor_ts_mad = np.ma.filled(meltfactor_ts_mad,0)
            meltfactor_ts_mean[meltfactor_ts_mean == -(np.inf)] = 1
            meltfactor_ts_mean[meltfactor_ts_mean <= 0] = 1
            meltfactor_ts_std[meltfactor_ts_std == -(np.inf)] = 0
            meltfactor_ts_std[meltfactor_ts_std <= 0] = 0
            meltfactor_ts_med[meltfactor_ts_med == -(np.inf)] = 1
            meltfactor_ts_med[meltfactor_ts_med <= 0] = 1
            meltfactor_ts_mad[meltfactor_ts_mad == -(np.inf)] = 0
            meltfactor_ts_mad[meltfactor_ts_mad <= 0] = 0
            outbins.extend([meltfactor_ts_mean, meltfactor_ts_std, meltfactor_ts_med, meltfactor_ts_mad])
        
        if self.vm is not None:
            outbins_header += ',vm_med,vm_mad'
            fmt += ',%0.2f,%0.2f'
            outbins.extend([vm_bin_med, vm_bin_mad])
        if self.H is not None:
            outbins_header += ',H_mean,H_std'
            fmt += ',%0.2f,%0.2f'
            outbins.extend([H_bin_mean, H_bin_std])
    
        if self.emvel is not None:
            outbins_header += ',emvel_mean,emvel_std,emvel_med,emvel_mad'
            fmt += ',%0.3f,%0.3f,%0.3f,%0.3f'
            outbins.extend([emvel_bin_mean, emvel_bin_std, emvel_bin_med, emvel_bin_mad])
        
        outbins = np.ma.array(outbins).T.astype('float32')
        np.ma.set_fill_value(outbins, np.nan)
        outbins = outbins.filled(np.nan)
        
        outbins_df = pd.DataFrame(outbins, columns=outbins_header.split(','))
        
#        if mb_df is not None:
#            # ADD MASS BALANCE DATA
#            mb_df = mb_df[np.isfinite(mb_df['bin_center_elev_m'])]
#            mb_df.reset_index(inplace=True, drop=True)
#            # start index for merge
#            if mb_df.loc[0,'bin_center_elev_m'] >= outbins_df.loc[0,'bin_center_elev_m']:
#                mb_df_idx1 = 0
#                outbins_idx1 = np.where(outbins_df['bin_center_elev_m'] == mb_df.loc[0,'bin_center_elev_m'])[0][0]
#            else:
#                outbins_idx1 = 0
#                mb_df_idx1 = np.where(outbins_df.loc[0,'bin_center_elev_m'] == mb_df['bin_center_elev_m'])[0][0]
#        #     print('idx1:', 
#        #           '\noutbins:', outbins_idx1, outbins_df.loc[outbins_idx1,'bin_center_elev_m'],
#        #           '\ndfbins:', mb_df_idx1, mb_df.loc[mb_df_idx1,'# bin_center_elev_m'])
#            # end index for merge
#            if (outbins_df.loc[outbins_df.shape[0]-1,'bin_center_elev_m'] >= 
#                mb_df.loc[mb_df.shape[0]-1,'bin_center_elev_m']):
#                outbins_idx2 = np.where(
#                    outbins_df['bin_center_elev_m'] == mb_df.loc[mb_df.shape[0]-1,'bin_center_elev_m'])[0][0]
#                mb_df_idx2 = mb_df.shape[0]-1
#            else:
#                outbins_idx2 = outbins_df.shape[0]-1
#                mb_df_idx2 = np.where(
#                    outbins_df.loc[outbins_df.shape[0]-1,'bin_center_elev_m'] == mb_df['bin_center_elev_m'])[0][0]
#        #     print('idx2:', 
#        #           '\noutbins:', outbins_idx2, outbins_df.loc[outbins_idx2,'bin_center_elev_m'],
#        #           '\ndfbins:', mb_df_idx2, mb_df.loc[mb_df_idx2,'# bin_center_elev_m'])
#            outbins_df[' mb_bin_mean_mwea'] = np.nan
#            outbins_df[' mb_bin_std_mwea'] = np.nan
#            outbins_df[' mb_bin_area_valid_km2'] = np.nan
#            outbins_df.loc[outbins_idx1:outbins_idx2+1,' mb_bin_mean_mwea'] = (
#                mb_df.loc[mb_df_idx1:mb_df_idx2+1,' mb_bin_mean_mwea'])
#            outbins_df.loc[outbins_idx1:outbins_idx2+1,' mb_bin_std_mwea'] = (
#                mb_df.loc[mb_df_idx1:mb_df_idx2+1,' mb_bin_std_mwea'])
#            outbins_df.loc[outbins_idx1:outbins_idx2+1,' mb_bin_area_valid_km2'] = (
#                mb_df.loc[mb_df_idx1:mb_df_idx2+1,' z1_bin_area_valid_km2'])
#            try:
#                outbins_df['startyear'] = mb_df.loc[mb_df_idx1,'startyear']
#                outbins_df['endyear'] = mb_df.loc[mb_df_idx1,'endyear']
#            except:
#                outbins_df['startyear'] = 2000
#                outbins_df['endyear'] = 2012
        
        if exportcsv:
            if int(self.feat_fn.split('.')[0]) < 10:
                outbins_fullfn = os.path.join(outdir_csv, self.feat_fn[0:7] + csv_ending)
            else:
                outbins_fullfn = os.path.join(outdir_csv, self.feat_fn[0:8] + csv_ending)
            outbins_df.to_csv(outbins_fullfn, index=False)
        
        outbins_df = pd.DataFrame(outbins, columns=outbins_header.split(','))
        
        return outbins_df, z_bin_edges
    
    
def create_glacfeat(thick_dir, thick_fn, verbose=False):
    """ Create the glacier feature"""
    
    glac_str = thick_fn.split('-')[-1].split('_')[0]
    region = thick_fn.split('-')[-1].split('.')[0]
    rgiid = thick_fn.split('_')[0]
    
    proj_fn = os.path.join(thick_dir, thick_fn) # THIS PROJECTION IS KEY!
    ds = gdal.Open(proj_fn)
    prj = ds.GetProjection()
    srs = osr.SpatialReference(wkt=prj)
    aea_srs = srs

    # Shape layer processing
    # If projected shapefile already exists, then skip projection
    glac_shp_proj_fn = (debris_prms.glac_shp_proj_fp + glac_str + '_crs' + 
                        str(aea_srs.GetAttrValue("AUTHORITY", 1)) + '.shp')
    if os.path.exists(glac_shp_proj_fn) == False:
        glac_shp_init = gpd.read_file(debris_prms.glac_shp_fn_dict[region])
        if verbose:
            print('Shp init crs:', glac_shp_init.crs)
        glac_shp_single = glac_shp_init[glac_shp_init['RGIId'] == rgiid]
        glac_shp_single = glac_shp_single.reset_index()
        glac_shp_proj = glac_shp_single.to_crs({'init': 'epsg:' + str(aea_srs.GetAttrValue("AUTHORITY", 1))})
        glac_shp_proj.to_file(glac_shp_proj_fn)
    glac_shp_ds = ogr.Open(glac_shp_proj_fn, 0)
    glac_shp_lyr = glac_shp_ds.GetLayer()
    #This should be contained in features
#    glac_shp_srs = glac_shp_lyr.GetSpatialRef()
    feat_count = glac_shp_lyr.GetFeatureCount()
    if verbose:
        print("Input glacier polygon count: %i" % feat_count)

    # ===== CREATE GLACIER FEATURE =====
    glacname_fieldname = "Name"
    glacnum_fieldname = "RGIId"
#    glacnum_fmt = '%08.5f'
    for n, feat in enumerate(glac_shp_lyr):
        gf = GlacFeat(feat, glacname_fieldname, glacnum_fieldname)
        if verbose:
            print("%i of %i: %s" % (n+1, feat_count, gf.feat_fn))
        #NOTE: Input must be in projected coordinate system, ideally equal area
        gf.geom_attributes(srs=aea_srs)
        gf.aea_srs = srs
    if verbose:
        print(gf.feat_fn)

    return gf


def nearest_nonzero_idx(a,x,y):
    r,c = np.nonzero(a)
    min_idx = ((r - x)**2 + (c - y)**2).argmin()
    return r[min_idx], c[min_idx]