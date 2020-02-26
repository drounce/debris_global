#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Feb 26 08:41:09 2020

@author: davidrounce
"""
import os
import re

import numpy as np
import pandas as pd
from osgeo import osr

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
        

    #RGI uses 50 m bins
    def hist_plot(self, bin_width=50.0, dz_clim=(-2.0, 2.0), exportcsv=True, csv_ending='', mb_df=None, outdir_csv=None):
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
    
        #If we only have one elevation grid with dhdt
        if self.z2 is not None:
            z2_bin_counts, z2_bin_edges = np.histogram(self.z2.compressed(), bins=z_bin_edges)
            z2_bin_areas = z2_bin_counts * self.res[0] * self.res[1] / 1E6
            #z2_bin_areas_perc = 100. * z2_bin_areas / np.sum(z2_bin_areas)
            z2_bin_areas_perc = 100. * (z1_bin_areas / self.glac_area_km2)
        else:
            z2_bin_counts = z1_bin_counts
            z2_bin_edges = z1_bin_edges
            z2_bin_areas = z1_bin_areas
            z2_bin_areas_perc = z1_bin_areas_perc
            
        if self.dc_area is not None:
            dc_bin_counts, dc_bin_edges = np.histogram(self.dc_area.compressed(), bins=z_bin_edges)
            dc_bin_areas = dc_bin_counts * self.res[0] * self.res[1] / 1E6
    #         dc_bin_areas_perc = 100. * (dc_bin_areas / self.glac_area_km2)
            dc_bin_areas_perc = 100. * (dc_bin_areas / z1_bin_areas)
    #         outbins_df['bin_debris_perc'] = outbins_df['dc_bin_count_valid'] / outbins_df['z1_bin_count_valid'] * 100
    
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
            slope_bin_samp = self.z1_slope[(idx == bin_n+1)]
            if slope_bin_samp.size > min_bin_samp_count:
                slope_bin_med[bin_n] = malib.fast_median(slope_bin_samp)
                slope_bin_mad[bin_n] = malib.mad(slope_bin_samp)
            aspect_bin_samp = self.z1_aspect[(idx == bin_n+1)]
            if aspect_bin_samp.size > min_bin_samp_count:
                aspect_bin_med[bin_n] = malib.fast_median(aspect_bin_samp)
                aspect_bin_mad[bin_n] = malib.mad(aspect_bin_samp)
    
        if self.dhdt is not None:
            dhdt_bin_areas = dhdt_bin_count * self.res[0] * self.res[1] / 1E6
            #dhdt_bin_areas_perc = 100. * dhdt_bin_areas / np.sum(dhdt_bin_areas)
            dhdt_bin_areas_perc = 100. * (dhdt_bin_areas / self.glac_area_km2)
    
        outbins_header = ('bin_center_elev_m,z1_bin_count_valid,z1_bin_area_valid_km2,z1_bin_area_perc,z2_bin_count_valid'
                          + ',z2_bin_area_valid_km2,z2_bin_area_perc,slope_bin_med,aspect_bin_med')
        fmt = '%0.1f,%0.0f,%0.3f,%0.2f,%0.0f,%0.3f,%0.2f,%0.2f,%0.2f'
        outbins = [z_bin_centers, z1_bin_counts, z1_bin_areas, z1_bin_areas_perc, z2_bin_counts, z2_bin_areas, 
                   z2_bin_areas_perc, slope_bin_med, aspect_bin_med]
        if self.dhdt is not None:
            outbins_header = (outbins_header + ',dhdt_bin_count,dhdt_bin_area_valid_km2,dhdt_bin_area_perc,' + 
                              'dhdt_bin_med_ma,dhdt_bin_mad_ma,dhdt_bin_mean_ma,dhdt_bin_std_ma,mb_bin_med_mwea,' + 
                              'mb_bin_mad_mwea,mb_bin_mean_mwea,mb_bin_std_mwea')
            fmt += ',%0.0f,%0.3f,%0.2f,%0.2f,%0.2f,%0.2f,%0.2f,%0.2f,%0.2f,%0.2f,%0.2f'
            outbins.extend([dhdt_bin_count, dhdt_bin_areas, dhdt_bin_areas_perc, dhdt_bin_med, dhdt_bin_mad, 
                            dhdt_bin_mean, dhdt_bin_std, mb_bin_med, mb_bin_mad, mb_bin_mean, mb_bin_std])
        if self.dc_dhdt is not None:
            outbins_header = (outbins_header + ',dc_dhdt_bin_count,dc_dhdt_bin_med_ma,dc_dhdt_bin_mad_ma,' +
                              'dc_dhdt_bin_mean_ma,dc_dhdt_bin_std_ma,dc_mb_bin_med_mwea,dc_mb_bin_mad_mwea,' + 
                              'dc_mb_bin_mean_mwea,dc_mb_bin_std_mwea')
            fmt += ',%0.0f,%0.2f,%0.2f,%0.2f,%0.2f,%0.2f,%0.2f,%0.2f,%0.2f'
            outbins.extend([dc_dhdt_bin_count, dc_dhdt_bin_med, dc_dhdt_bin_mad, dc_dhdt_bin_mean, dc_dhdt_bin_std, 
                            dc_mb_bin_med, dc_mb_bin_mad, dc_mb_bin_mean, dc_mb_bin_std])
        if self.dc_area is not None:
            outbins_header += ',dc_bin_count_valid,dc_bin_area_valid_km2,dc_bin_area_perc'
            fmt += ',%0.0f,%0.3f,%0.2f'
            outbins.extend([dc_bin_counts, dc_bin_areas, dc_bin_areas_perc])
    #         outbins.extend([z1_bin_counts, z1_bin_areas, z1_bin_areas_perc])
            
            
    #    if self.debris_thick is not None:
    #        outbins_header += ',hd_med_m,hd_mad_m'
    #        fmt += ',%0.2f,%0.2f'
    #        debris_thick_med[debris_thick_med == -(np.inf)] = 0.00
    #        debris_thick_mad[debris_thick_mad == -(np.inf)] = 0.00
    #        outbins.extend([debris_thick_med, debris_thick_mad])
        
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
        
        if mb_df is not None:
            # ADD MASS BALANCE DATA
            mb_df = mb_df[np.isfinite(mb_df['bin_center_elev_m'])]
            mb_df.reset_index(inplace=True, drop=True)
            # start index for merge
            if mb_df.loc[0,'bin_center_elev_m'] >= outbins_df.loc[0,'bin_center_elev_m']:
                mb_df_idx1 = 0
                outbins_idx1 = np.where(outbins_df['bin_center_elev_m'] == mb_df.loc[0,'bin_center_elev_m'])[0][0]
            else:
                outbins_idx1 = 0
                mb_df_idx1 = np.where(outbins_df.loc[0,'bin_center_elev_m'] == mb_df['bin_center_elev_m'])[0][0]
        #     print('idx1:', 
        #           '\noutbins:', outbins_idx1, outbins_df.loc[outbins_idx1,'bin_center_elev_m'],
        #           '\ndfbins:', mb_df_idx1, mb_df.loc[mb_df_idx1,'# bin_center_elev_m'])
            # end index for merge
            if (outbins_df.loc[outbins_df.shape[0]-1,'bin_center_elev_m'] >= 
                mb_df.loc[mb_df.shape[0]-1,'bin_center_elev_m']):
                outbins_idx2 = np.where(
                    outbins_df['bin_center_elev_m'] == mb_df.loc[mb_df.shape[0]-1,'bin_center_elev_m'])[0][0]
                mb_df_idx2 = mb_df.shape[0]-1
            else:
                outbins_idx2 = outbins_df.shape[0]-1
                mb_df_idx2 = np.where(
                    outbins_df.loc[outbins_df.shape[0]-1,'bin_center_elev_m'] == mb_df['bin_center_elev_m'])[0][0]
        #     print('idx2:', 
        #           '\noutbins:', outbins_idx2, outbins_df.loc[outbins_idx2,'bin_center_elev_m'],
        #           '\ndfbins:', mb_df_idx2, mb_df.loc[mb_df_idx2,'# bin_center_elev_m'])
            outbins_df[' mb_bin_mean_mwea'] = np.nan
            outbins_df[' mb_bin_std_mwea'] = np.nan
            outbins_df[' mb_bin_area_valid_km2'] = np.nan
            outbins_df.loc[outbins_idx1:outbins_idx2+1,' mb_bin_mean_mwea'] = (
                mb_df.loc[mb_df_idx1:mb_df_idx2+1,' mb_bin_mean_mwea'])
            outbins_df.loc[outbins_idx1:outbins_idx2+1,' mb_bin_std_mwea'] = (
                mb_df.loc[mb_df_idx1:mb_df_idx2+1,' mb_bin_std_mwea'])
            outbins_df.loc[outbins_idx1:outbins_idx2+1,' mb_bin_area_valid_km2'] = (
                mb_df.loc[mb_df_idx1:mb_df_idx2+1,' z1_bin_area_valid_km2'])
            try:
                outbins_df['startyear'] = mb_df.loc[mb_df_idx1,'startyear']
                outbins_df['endyear'] = mb_df.loc[mb_df_idx1,'endyear']
            except:
                outbins_df['startyear'] = 2000
                outbins_df['endyear'] = 2012
        
        if exportcsv:
            if int(self.feat_fn.split('.')[0]) < 10:
                outbins_fullfn = os.path.join(outdir_csv, self.feat_fn[0:7] + csv_ending)
            else:
                outbins_fullfn = os.path.join(outdir_csv, self.feat_fn[0:8] + csv_ending)
            outbins_df.to_csv(outbins_fullfn, index=False)
        
        outbins_df = pd.DataFrame(outbins, columns=outbins_header.split(','))
        
        return outbins_df, z_bin_edges