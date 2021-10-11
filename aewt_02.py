# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 10:08:19 2021

@author: mar886

Basing this script on my PhD -- Baby Hope "BH_05"
"""
    
import sys
import os
import shutil
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import flopy
from flopy.utils.reference import SpatialReference
import flopy.utils.binaryfile as bf
import subprocess
from flopy import mt3d
import pandas as pd
import fiona
import geopandas as gpd
from rasterio.features import rasterize
from rasterio import Affine
import shapefile as sf
from matplotlib.patches import Ellipse, Polygon
from itertools import groupby
import contextily
import shapely
from shapely.geometry import Point


print(sys.version)
print('numpy version: {}'.format(np.__version__))
print('matplotlib version: {}'.format(mpl.__version__))
print('flopy version: {}'.format(flopy.__version__))

#=====DIRECTORIES==============================================================

overallName = "AE_01"
#modelname = overallName
#modelname_mt = modelname + "_MT3D" 

#project_folder = r'C:\workspace\Proj4_BabyHope'

#MODFLOW_folder      = r'C:\workspace\modflow_dbl\mfnwtdbl.exe' # double precision MODFLOW NWT
#mfnwt_exe_name = "mfnwtdbl.exe" 

#MT3D_USGS_folder      = r'C:\workspace\modflow_dbl\mt3dusgsdbl.exe' # Double precision MT3D-USGS
#mt3d_version = 'mt3d-usgs' 
#mt3d_exe_name = "mt3dusgsdbl.exe"

#spatialFolder = r'C:\SMarshall_PhD\BabyHope\Sample_results\final_compiled'
#spatialFolder2 = r'C:\SMarshall_PhD\Spatial\from_Saskia\East_Pilbara_Concept'

#arc_output_fldr = r'C:\SMarshall_PhD\Spatial\Baby_Hope_Sampling_Location'

#------------------------------------------------------------------------------
# FOLDER STRUCTURES

os.getcwd() # Something weird is going on with this, unless I run using green triangle at the top, it is giving the wrong cwd

# get current directory
path = os.getcwd()
print("Current Directory", path)

# Make folders for data and figures
if not os.path.exists("input_data"):
    os.makedirs("input_data")
       
if not os.path.exists("figures"):
    os.makedirs("figures")
    
if not os.path.exists("output_data"):
    os.makedirs("output_data")
   
#------------------------------------------------------------------------------
# DATA SET UP

wt_contour_fldr = os.path.join("input_data", "GAB_WT", "Watertable_Contour")
 
# epsg code specifying coordinate reference system
model_epsg = 32755 # Got zone from here: https://spatialreference.org/ref/epsg/wgs-84-utm-zone-50s/
# WGS 84 / UTM zone 55S - this is an estimate, need to check


studyarea = shapely.geometry.box(-800000, 7100000, -200000, 7450000) # minx, miny, maxx, maxy

studyarea_df = gpd.GeoDataFrame(geometry=studyarea)

# See this for making polygon https://gis.stackexchange.com/questions/208069/make-a-shapefile-with-geopandas-from-polygons

# Shape files -----------------------------------------------------------
# Convert all to same projection (UTM)

wt_contour_1 = gpd.read_file(os.path.join(wt_contour_fldr, 'wt_contour.shp'))



wt_contour_1.crs # Check existing coordinate reference system
#? It doesn't have an existing coordinate reference system?
wt_contour_1_utm = wt_contour_1.to_crs(epsg=32755)

wt_contour_1_utm.crs # It should now be in UTM

print(wt_contour_1_utm.type)
print(type(wt_contour_1_utm))
print(type(wt_contour_1_utm.geometry))

print(wt_contour_1_utm.head())
print(wt_contour_1_utm.area)
# Prev issue with conversion solved, try this:https://stackoverflow.com/questions/55390492/runtimeerror-bno-arguments-in-initialization-list

# Checking out the data ------------------------------------------------------

print(wt_contour_1_utm.head())

# Make a plot
fig, ax = plt.subplots()
wt_contour_1_utm.plot(ax=ax, color="k", alpha=0.5, linewidth=2)


#contextily.add_basemap(ax, crs=wt_contour_1_utm.crs) # The basemap is in the complete wrong area???

# Adding DEM ------------------------------------------------------------------

























