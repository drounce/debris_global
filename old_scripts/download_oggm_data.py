#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
OVERVIEW: Download DEMs from OGGM
environment: pygem_v3_viapip
"""

## Built-in libaries
#import argparse
#import collections
##import datetime
#import multiprocessing
#import os
#import pickle
import time
## External libraries
##import matplotlib.pyplot as plt
#import numpy as np
#import pandas as pd
#import xarray as xr
## Local libraries
#import globaldebris_input as input
#from spc_split_lists import split_list


import pandas as pd
import numpy as np
from oggm import cfg, utils, workflow, tasks, graphics, GlacierDirectory
import xarray as xr
import geopandas as gpd
import salem
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import AxesGrid
import itertools

from oggm.utils import DEM_SOURCES
from oggm.workflow import init_glacier_regions

# Make sure the plot directory exists
utils.mkdir(plot_dir);
# Use OGGM to download the data
cfg.initialize()
cfg.PATHS['working_dir'] = utils.gettempdir(dirname='OGGM-DEMS', reset=True)
cfg.PARAMS['use_intersects'] = False

# URL of the preprocessed GDirs
gdir_url = 'https://cluster.klima.uni-bremen.de/data/gdirs/dems_v0/'
# We use OGGM to download the data
gdir = init_glacier_regions([rgi_id], from_prepro_level=1, prepro_border=10, 
                            prepro_rgi_version='62', prepro_base_url=gdir_url)[0]

if sources is None:
    sources = [src for src in os.listdir(gdir.dir) if src in utils.DEM_SOURCES]



print('RGI ID:', rgi_id)
print('Available DEM sources:', sources)
print('Plotting directory:', plot_dir)
                
#%%
if __name__ == '__main__':
    time_start = time.time()
    
    print('Download oggm glacier data:\n - widths \n - DEMs')
    
    
    
    print('\nProcessing time of :',time.time()-time_start, 's')
