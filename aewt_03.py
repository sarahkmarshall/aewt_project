# -*- coding: utf-8 -*-
"""
Created on Tue Oct  5 10:08:19 2021

@author: mar886

Basing this script on my PhD -- Baby Hope "BH_05
Using luk's jupyter notebook: https://ppts-bracewell-i1.hpc.csiro.au:9503/notebooks/Methodology/methodology/GKIS_Load_DataSets.ipynb#"
"""
    
import contextily
import flopy
from flopy.utils.reference import SpatialReference
import flopy.utils.binaryfile as bf
from flopy import mt3d
import fiona
import geopandas as gpd
from itertools import groupby
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
import os
import numpy as np
import pyproj
import pandas as pd
import rasterio as rio
from rasterio.features import rasterize
from rasterio import Affine
from rasterio.plot import show
from rasterio.transform import from_bounds
import shutil
import shapefile as sf
import subprocess
import sys
import shapely
from shapely.geometry import Point
from shapely.geometry import mapping, Polygon
from owslib.wcs import WebCoverageService

print('python version: {}' .format(sys.version))
print('numpy version: {}'.format(np.__version__))
print('matplotlib version: {}'.format(mpl.__version__))
print('flopy version: {}'.format(flopy.__version__))
print('pyproj version: {}'.format(pyproj.__version__))

# GDAL 2.3.3, released 2018/12/14. Got this by typing into conda prompt: gdalinfo --version
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
os.chdir(r'C:\Users\mar886\WaterTableProject\aewt_project')

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
# DATA SET UP & PROJECTIONS

wt_contour_fldr = os.path.join("input_data", "GAB_WT", "Watertable_Contour")
 
# epsg code specifying coordinate reference system
model_epsg = 32755 # Got zone from here: https://spatialreference.org/ref/epsg/wgs-84-utm-zone-50s/
# WGS 84 / UTM zone 55S - this is an estimate, need to check

wgs84 = pyproj.CRS('epsg:4326')
utm = pyproj.CRS('epsg:32755')
project = pyproj.Transformer.from_crs(wgs84, utm, always_xy=True).transform

gda94 = pyproj.CRS('epsg:4283')

# DEFINING STUDY AREA ---------------------------------------------------------

# This is currently in UTM
minx = -800000
miny = 7100000
maxx = -200000
maxy = 7450000


studyarea = shapely.geometry.box(minx, miny, maxx, maxy) # minx, miny, maxx, maxy
#sa_polygon = Polygon([(minx, miny), (maxx, miny), (maxx, maxy), (minx, maxy)])
print(type(studyarea))

studyarea.area
studyarea.length

sa_df = pd.DataFrame()
sa_df["geometry"] = [studyarea]

sa_df.head()
sa_gdf = gpd.GeoDataFrame(sa_df, geometry='geometry',crs="epsg:32755")
sa_gdf.head()
sa_gdf.columns
print(sa_gdf.crs)


sa_gdf_wgs84 = sa_gdf.to_crs(wgs84)
sa_gdf_wgs84.crs
print(sa_gdf_wgs84.geometry[0])

sa_gdf_gda94 = sa_gdf.to_crs(gda94)
sa_gdf_gda94.crs
print(sa_gdf_gda94.geometry[0])


# -----------------------------------------------------------------------------
# Getting the bounding box/study area in wgs84 coordinates so I can use it as the 
# bounding box to extract the elevation data

which_sa = sa_gdf_wgs84

which_sa.geometry[0] # It's a wonky shape because of the transformation

coordinates_array = np.asarray(which_sa.geometry[0].exterior.coords)[:-1]# Note that the first and last points are the same


latitudes = []
longitudes = []
for i, crd in enumerate(coordinates_array):
    latitudes.append(coordinates_array[i][1])
    longitudes.append(coordinates_array[i][0])


#flat_coords = coordinates_array.flatten() # This is to easily get max and 
#min values to get the maximum extent for the dem (because it's a bit warped and not a proper rectangle once you change coordinate systems)
# Issue with this is that it doesn't work for the northern hemisphere.
#latitudes = [coord for coord in flat_coords if coord < 0]
#longitudes = [coord for coord in flat_coords if coord > 0]

north_extent  =   max(latitudes)
south_extent  =   min(latitudes)
west_extent   =   min(longitudes)
east_extent   =   max(longitudes)

print("north extent:   %2.2f" %north_extent)
print("south extent:   %2.2f" %south_extent)
print("west extent:   %2.2f" %west_extent)
print("east extent:   %2.2f" %east_extent)

ew_extent = abs(east_extent - west_extent)
sn_extent = abs(south_extent - north_extent)

size_of_box = ew_extent*sn_extent # Degrees squared
# For elvis the size of the downloaded area is supposed to be 1.5 degrees squared?

stepsize = 0.01 # 1.2 gets me the size I need for elvis
long_intervals = np.arange(west_extent, east_extent, stepsize)
lat_intervals  = np.arange(south_extent, north_extent, stepsize)

# Elvis online portal is WGS84, epsg 4326


# Start with the first square, south eastern most part of the map

n1 = lat_intervals[1]
s1 = lat_intervals[0]
w1 = long_intervals[-2]
e1 = long_intervals[-1]

print("north extent:   %2.2f" %n1)
print("south extent:   %2.2f" %s1)
print("west extent:   %2.2f" %w1)
print("east extent:   %2.2f" %e1)

# Create bounding box (bbox) -
# this can be used for automatic extraction of tiff files

bbox_sa = '%f,%f,%f,%f' % (w1,s1,e1,n1)

#------------------------------------------------------------------------------
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# DEM Data

# I downloaded as geotiff or as esri ascii data
dem_fldr = os.path.join("input_data", "Elvis_dl_25Oct21_geotiff")
dem_rstr_nm = 'Hydro_Enforced_1_Second_DEM.tif'

#dem_rstr = rio.open(os.path.join(dem_fldr, dem_rstr_nm))

#plt.figure()
#plt.imshow(dem_rstr,cmap='terrain')
#plt.colorbar()
#plt.axis('off')

# Getting error saying: 'RasterioIOError: This is a BigTIFF file.  BigTIFF is not supported by this version of GDAL and libtiff.'
# ??? Not sure what I can do

#------------------------------------------------------------------------------
# Practicing with a smaller tiff
dem_fldr_small = os.path.join("input_data", "Elvis_dl_27Oct21_geotiff")

fname_dem_small = os.path.join(dem_fldr_small, dem_rstr_nm)
#with rio.open(fname_dem_small) as grd:
#    hsu = grd.read()[0,:,:]
#    grid_meta = grd.profile   # Gets the metadata automatically
#    bounds = grd.bounds
#    res = grd.res
# Still says that it is a BigTIFF file, even though it is a small file.

#------------------------------------------------------------------------------

# Import ascii file
dem_fldr_asc = os.path.join("input_data", "Elvis_dl_25Oct21")
dem_rstr_nm_asc = 'Hydro_Enforced_1_Second_DEM.asc'

dem_rstr_as = rio.open(os.path.join(dem_fldr_asc, dem_rstr_nm_asc))
dem_rstr_as.crs # 4326

## A different way to import, Luk's jupyter nb
with rio.open(os.path.join(dem_fldr_asc, dem_rstr_nm_asc)) as grd:
    hsu = grd.read()[0,:,:]
    grid_meta = grd.profile   # Gets the metadata automatically
    bounds = grd.bounds
    res = grd.res

mask = np.zeros_like(hsu)
nrows,ncols = np.shape(mask)

print("nrows: %s, ncols: %s" %(nrows,ncols))

###############################################################################
# bounds: BoundingBox(left=138.92013888913038, bottom=-26.049861111151657, right=140.1198611113643, top=-24.849861111139962)
# res: (0.000277777777780485, 0.000277777777780485)
###############################################################################

# Try to download dem directly from web services

wcs = WebCoverageService('https://services.ga.gov.au/site_9/services/DEM_SRTM_1Second_Hydro_Enforced/MapServer/WCSServer',
                         version='1.0.0')

print(wcs.contents.keys)



bounds2 = (140.167619, 26.049768, 140.177619, -26.039768)

dat = wcs.getCoverage(identifier='1', bbox=bounds2, format='GeoTIFF',
                     crs='4283',resx=res[0],resy=res[0])


fname = 'DEM_1s_HE.tif'

with open(os.path.join("input_data", fname),'wb') as file:
    file.write(dat.read())
    
with rio.open(os.path.join("input_data", fname)) as grd:
    dem = grd.read()[0,:,:]
    
# I'm getting an error: RasterioIOError: 'input_data\DEM_1s_HE.tif' not recognized as a supported file format.

plt.figure()
plt.imshow(dem,cmap='terrain')
plt.colorbar()
plt.axis('off')




#------------------------------------------------------------------------------
# Shape files -----------------------------------------------------------------
# Convert all to same projection (UTM)

wt_contour_1 = gpd.read_file(os.path.join(wt_contour_fldr, 'wt_contour.shp'))
wt_contour_1.crs # Check existing coordinate reference system
#? It doesn't have an existing coordinate reference system?
wt_contour_1_utm = wt_contour_1.to_crs(epsg=model_epsg)

# Define study area crs

wt_contour_1_utm.crs # It should now be in UTM

print(wt_contour_1_utm.type)
print(type(wt_contour_1_utm))
print(type(wt_contour_1_utm.geometry))

print(wt_contour_1_utm.head())
print(wt_contour_1_utm.area)
# Prev issue with conversion solved, try this:https://stackoverflow.com/questions/55390492/runtimeerror-bno-arguments-in-initialization-list

# Convert study area polygon to the same projection

#sa_gdf.to_crs({'init': 'epsg:32755', 'no_defs': True}) # sa_gdf.to_crs(wt_contour_1_utm.crs)

# Checking out the data ------------------------------------------------------

print(wt_contour_1_utm.head())

# Make a plot
fig, ax = plt.subplots()
wt_contour_1_utm.plot(ax=ax, color="k", alpha=0.5, linewidth=2)
#plt.plot(*studyarea.exterior.xy, color="m")

sa_gdf.plot(ax=ax, color="gold", alpha=0.5, label="Study extent", legend=True)
plt.xlabel("Longitude")
plt.ylabel("Latitude")

# Make a legend

pmarks = []
color = "gold"
label = "Study extent"
pmarks.append(Patch(facecolor=color, label=label))

handles, _ = ax.get_legend_handles_labels()
ax.legend(handles=[*handles,*pmarks], loc='lower right')

# Save figure
plt.savefig(os.path.join("figures", "original_contours"), dpi=300)


# Crop the data to the study area ---------------------------------------------
print(sa_gdf.crs)
print(wt_contour_1_utm.crs)


#contours_icpt = wt_contour_1_utm[wt_contour_1_utm.geometry.intersects(studyarea)]
#contours_icpt = wt_contour_1_utm.intersection(studyarea) # This does work also

# AHA this works! It needs to be just the polygon to do intersection....that's good to know.
contours_icpt = wt_contour_1_utm.intersection(sa_gdf.geometry[0])

fig, ax = plt.subplots()
contours_icpt.plot(ax=ax, color="k", alpha=0.5, linewidth=2)
sa_gdf.plot(ax=ax,color="gold", alpha=0.5)

# Save figure
plt.savefig(os.path.join("figures", "cropped_contours"), dpi=300)

#plt.plot(*studyarea.exterior.xy, color="m")
#contextily.add_basemap(ax, crs=wt_contour_1_utm.crs) # The basemap is in the complete wrong area???

# convert back to wgs84

contours_sa_wgs84 = contours_icpt.to_crs('epsg:4326')
contours_sa_wgs84.crs

# Adding DEM ------------------------------------------------------------------
# Using this website to help: https://rasterio.readthedocs.io/en/latest/quickstart.html

dem_rstr_as.name
dem_rstr_as.closed
dem_rstr_as.count
dem_rstr_as.width
dem_rstr_as.height
dem_bounds = dem_rstr_as.bounds

# Values are relative to the origin of the dataset's crs
upperleft = dem_rstr_as.transform * (0,0)
lowerright = dem_rstr_as.transform * (dem_rstr_as.width, dem_rstr_as.height)

# Reading the raster data
band1 = dem_rstr_as.read(1) # Read from 1 as the first

# To get the spatial coordinates of a pixel:
# This gets the coordinates of the centre of the image
dem_rstr_as.xy(dem_rstr_as.height // 2, dem_rstr_as.width // 2)

### Plotting the raster

# This line below doesn't seem to be plotting on the same figure
plt.imshow(dem_rstr_as.read(1), cmap='pink')

# Plot DEM with contours
fig1, ax1 = plt.subplots()
show(dem_rstr_as, ax=ax1, alpha=0.5)
contours_sa_wgs84.plot(ax=ax1, color="k", alpha=0.5, linewidth=1)

#------------------------------------------------------------------------------
# Re-crop the contours to the new study area - just within the dem extent

# NEW study area

studyarea2 = shapely.geometry.box(dem_bounds[0], # minx
                                 dem_bounds[1],  # miny, 
                                 dem_bounds[2],  # maxx, 
                                 dem_bounds[3])  # maxy
print(type(studyarea2))

sa_df2 = pd.DataFrame()
sa_df2["geometry"] = [studyarea2]

sa_gdf2 = gpd.GeoDataFrame(sa_df2, geometry='geometry',crs=wgs84)

# Now crop the contours again

contours_icpt2 = contours_sa_wgs84.intersection(sa_gdf2.geometry[0])
contours_sa_wgs842 = contours_icpt2.to_crs('epsg:4326')
contours_sa_wgs842.crs

# Plot DEM with contours
fig1, ax1 = plt.subplots()
show(dem_rstr_as, ax=ax1, alpha=0.5)
contours_sa_wgs842.plot(ax=ax1, color="k", alpha=0.5, linewidth=1)

#------------------------------------------------------------------------------
#------------------------------------------------------------------------------
# Turn my water table contours into a raster file
df= wt_contour_1 # contours_icpt2
shape = 1000, 1000

# I don't really know what this transform line does
transform = rio.transform.from_bounds(*df['geometry'].total_bounds, *shape)

# Why has cropping the wt df converted it into a geoseries?

raster_wt = rasterize(
    ((s, h) for s, h in zip(df['geometry'], df['height'])),
    out_shape=shape,
    transform=transform,
    fill = 0,
    all_touched = True,
    default_value = 0,
    dtype=rio.int32)

raster_wt = rasterize(
    [(s, 1) for s in df['geometry']],
    out_shape=shape,
    transform=transform,
    fill=0,
    all_touched=True,
    dtype=rio.uint8)

plt.figure()
plt.imshow(raster_wt, cmap='cool')
plt.colorbar()
plt.axis('off')

type(raster_wt)
np.size(raster_wt)
# Get some stats
raster_wt.min()
raster_wt.max()

raster_wt_contours = rasterize(
    [(shape, 1) for shape in df['geometry']],
    out_shape=shape,
    transform=transform,
    fill=0,
    all_touched=True,
    dtype=rio.uint8)

type(raster_wt_contours)
np.size(raster_wt_contours)
    
# Plot rasterised contours
fig1, ax1 = plt.subplots()
show(raster_wt_contours, ax=ax1, alpha=0.5)

# Get some stats
raster_wt_contours.min()
raster_wt_contours.max()
# It's just 1 where there is a contour, no attribute assigned


# New plot
plt.figure()
plt.imshow(raster_wt_contours, cmap='terrain')
plt.colorbar()
plt.axis('off')


# Write it to a file
with rio.open(
    os.path.join('output_data', 'rasterise_wt.tif'), 'w',
    driver='GTiff',
    dtype=rio.uint8,
    count=1,
    width=shape[0],
    height=shape[1],
    transform=transform
) as dst:
    dst.write(raster_wt_contours, indexes=1)
    
    #### Trying to change the default value
    
shapes = ((geom,value) for geom, value in zip(df.geometry, df["height"]))

    
varray = rasterize(shapes,
                   out_shape = shape,
                   transform = transform,
                   fill = 0,
                   all_touched = True,
                   dtype = rio.int32)    
    
raster_wt = rasterize(
    [(s, 1) for s in df['geometry']],
    out_shape=shape,
    transform=transform,
    fill=0,
    all_touched=True,
    dtype=rio.uint8)