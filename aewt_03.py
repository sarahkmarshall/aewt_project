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
from requests import Request
import shutil
import shapefile as sf
import subprocess
from scipy.interpolate import interpn
from scipy.interpolate import griddata
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
 
wgs84 = pyproj.CRS('epsg:4326')
utm   = pyproj.CRS('epsg:32755')
gda94 = pyproj.CRS('epsg:4283')

#------------------------------------------------------------------------------
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# DEM Data --------------------------------------------------------------------

# Start with the DEM data to make the grid, this sets the shape for the other 
# datasets to work with/crop to

# Import ascii file (because I am having trouble working with tiffs)
dem_fldr_asc = os.path.join("input_data", "Elvis_dl_25Oct21")
dem_rstr_nm_asc = 'Hydro_Enforced_1_Second_DEM.asc'

# Import data with the crs included
dem_projected = rio.open(os.path.join(dem_fldr_asc, dem_rstr_nm_asc))
dem_projected.crs

# Now import the data only as a raster object (i.e. as an array)
with rio.open(os.path.join(dem_fldr_asc, dem_rstr_nm_asc)) as grd:
    dem = grd.read()[0,:,:]
    grid_meta = grd.profile   # Gets the metadata automatically
    bounds = grd.bounds
    res = grd.res

#Use these parameters from the download to set up my extent for future data analysis
dem.shape
mask = np.zeros_like(dem)
nrows,ncols = np.shape(mask)

print("nrows: %s, ncols: %s" %(nrows,ncols))

#------------------------------------------------------------------------------
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# Water table contours (shapefiles) -------------------------------------------

wt_contour = gpd.read_file(os.path.join(wt_contour_fldr, 'wt_contour.shp'))
wt_contour.crs # Check existing coordinate reference system

wt_contour_wgs84 = wt_contour.to_crs(wgs84)
wt_contour_wgs84.crs

# Checking out the data ------------------------------------------------------

print(wt_contour_wgs84.head())

# Make a plot
fig, ax = plt.subplots()
wt_contour_wgs84.plot(ax=ax, color="k", alpha=0.5, linewidth=0.5)


ax.add_patch(mpl.patches.Polygon([[bounds[0], bounds[1]],   # bottom left
                                  [bounds[2], bounds[1]],   # bottom right
                                  [bounds[2], bounds[3]],   # top right
                                  [bounds[0], bounds[3]]], # top left

                                  facecolor="m", edgecolor="m", closed=True, 
                                  lw=0, fill=True))

plt.xlabel("Longitude")
plt.ylabel("Latitude")

# Make a legend

pmarks = []
color = "m"
label = "Study extent"
pmarks.append(Patch(facecolor=color, label=label))

handles, _ = ax.get_legend_handles_labels()
ax.legend(handles=[*handles,*pmarks], loc='lower right')

# Save figure
plt.savefig(os.path.join("figures", "original_contours"), dpi=300)


# Make a geodataframe of the bounding box to crop with ------------------------

studyarea = shapely.geometry.box(bounds[0], bounds[1], bounds[2], bounds[3]) # minx, miny, maxx, maxy

studyarea.area
studyarea.length

sa_df = pd.DataFrame()
sa_df["geometry"] = [studyarea]

sa_gdf = gpd.GeoDataFrame(sa_df, geometry='geometry',crs=wgs84)

# Crop the data to the study area ---------------------------------------------

# First get rid of all contours that aren't contours (i.e. the gab boundary)
wt_contour_wgs84_filtr = wt_contour_wgs84[wt_contour_wgs84['height']>-9999]

print(sa_gdf.crs)
print(wt_contour_wgs84_filtr.crs)


type(sa_gdf)
type(wt_contour_wgs84_filtr)
len(wt_contour_wgs84_filtr) # 1122 rows

# Make a copy of the dataframe so I can add the intersection geoseries over the top
contours_icpt_gdf = wt_contour_wgs84_filtr.copy()
len(contours_icpt_gdf)

# Crop contours based on the study area
contours_icpt_gdf['contours_icpt'] = wt_contour_wgs84_filtr.intersection(sa_gdf.geometry[0])

type(contours_icpt_gdf)
len(contours_icpt_gdf) # 1122 rows
type(contours_icpt_gdf['contours_icpt'])
type(contours_icpt_gdf['contours_icpt'][0])

# Get only values where the geometry is not empty
print(contours_icpt_gdf['contours_icpt'][0])


is_empty_list = []
for i in range(len(contours_icpt_gdf)):
    is_empty_list.append(contours_icpt_gdf.iloc[i, 4].is_empty)
len(is_empty_list)


contours_icpt_gdf['is_empty'] = is_empty_list
contours_icpt_gdf_notempty = contours_icpt_gdf[contours_icpt_gdf['is_empty'] == False]

# Plot this version of the data
fig, ax = plt.subplots()
contours_icpt_gdf_notempty.plot(ax=ax, color="k", alpha=0.5, linewidth=2)

# Now it is filtered make a new one where the name geometry is changed

contours_sa = contours_icpt_gdf_notempty.copy()
contours_sa = contours_sa.drop(labels="geometry", axis=1)
contours_sa.columns
contours_sa = contours_sa.rename(columns={"contours_icpt":"geometry"})

########################################
fig, ax = plt.subplots()
contours_sa.plot(ax=ax, color="k", alpha=0.5, linewidth=1)

# Now it is a geopandas opject and is the right shape, cropped to the study area

sa_gdf.plot(ax=ax,color="gold", alpha=0.5)

plt.xlabel("Longitude")
plt.ylabel("Latitude")
# Save figure
plt.savefig(os.path.join("figures", "cropped_contours"), dpi=300)
#######################################

######################################
# Plot DEM with the contours

fig1, ax1 = plt.subplots()
show(dem_projected, ax=ax1, alpha=0.5)
contours_sa.plot(ax=ax1, color="k", alpha=0.5, linewidth=1)
plt.xlabel("Longitude")
plt.ylabel("Latitude")
ymin, ymax = ax1.get_ylim()
xmin, xmax = ax1.get_xlim()
axes_extent = ax1.axis()

plt.savefig(os.path.join("figures", "dem_with_contours"), dpi=300)
######################################


#------------------------------------------------------------------------------
# Turn my water table contours into a raster file
df = contours_sa # contours_icpt2
shape = dem.shape # 5000, 5000

# I don't really know what this transform line does
transform = grid_meta['transform'] # rio.transform.from_bounds(*df['geometry'].total_bounds, *shape)

# Why has cropping the wt df converted it into a geoseries?

raster_wt = rasterize(
    ((s, h) for s, h in zip(df['geometry'], df['height'])),
    out_shape=shape,
    transform=transform,
    fill = 0,
    all_touched = True,
    default_value = 0,
    dtype=rio.int32)

plt.figure()
plt.imshow(raster_wt, cmap='gist_ncar_r')
plt.colorbar()
plt.xlabel("x")
plt.ylabel("y")

plt.savefig(os.path.join("figures", "raster_contours"), dpi=300)


#------------------------------------------------------------------------------
# Download river network

# url to Geofabric Surface Hydrology Cartography
WfsUrl = 'http://geofabric.bom.gov.au/simplefeatures/ahgf_shcarto/ows'
params = dict(service='WFS', 
              version='1.0.0', 
              request='GetFeature',
              typename='ahgf_shcarto:AHGFMappedStream', 
              bbox='%f,%f,%f,%f' % (bounds[0],bounds[1],bounds[2],bounds[3]),
              outputFormat='json')

q = Request('GET',WfsUrl,params=params).prepare().url
riv_wgs = gpd.read_file(q)
riv_wgs.crs # Already in wgs84
river_gdf = riv_wgs.to_crs(wgs84)


# PLot the river network

fig1, ax1 = plt.subplots()
show(dem_projected, ax=ax1, alpha=0.5)
contours_sa.plot(ax=ax1, color="k", alpha=0.5, linewidth=1)
river_gdf.plot(ax=ax1, edgecolor='b')

plt.xlabel("Longitude")
plt.ylabel("Latitude")

ax1.set_xlim([xmin, xmax])
ax1.set_ylim([ymin, ymax])

plt.savefig(os.path.join("figures", "dem_riv_contours"), dpi=300)


#------------------------------------------------------------------------------
# Vectorise river network

# Select major rivers
inds = river_gdf.hierarchy=='Major'

# Create features for rasterio
riv_shp = ((geom,value) for geom, value in zip(river_gdf[inds].geometry, 
                                               river_gdf[inds].hydroid))

# burn in features in raster
riv_brn = rasterize(shapes=riv_shp, out_shape=dem.shape, fill=0, 
                    transform=grid_meta['transform'])

# change grid values so river cells are 1
riv_brn[riv_brn > 0] = 1

# Plot the rasterised river network

plt.figure()
plt.imshow(riv_brn, cmap='gist_stern', interpolation=None)
plt.savefig(os.path.join("figures", "riv_raster"), dpi=300)

#------------------------------------------------------------------------------
n_elements = raster_wt.shape[0]*raster_wt.shape[1]

# This is what I want to interpolate:
raster_wt
type(raster_wt)
raster_wt.shape

# Define a 3d np array based on this raster, i.e. change to [x, y, z] values 

x = np.arange(0, ncols)
y = np.arange(0, nrows)

xi, yi = np.meshgrid(x, y) # , sparse=True

print("number of rows: %i should be equal to y: %i" %(nrows, len(y)))
print("number of cols: %i should be equal to x: %i" %(ncols, len(x)))

# Find the points where the water table exists
points_wt = (raster_wt).nonzero()   

# Find the points where the water table doesn't exist
points_no_wt = np.where(raster_wt == 0) 

# Now find the values of the water table at each of these points
values_wt = raster_wt[points_wt]

# BUT I need to preserve the structure of the array

# Now do an interpolation over every other point in the raster    

raster_wt_interp = raster_wt.copy()


# Get interpolated value for a point

grid_wt_z0 = griddata(points_wt, values_wt, (xi, yi), method='nearest')
grid_wt_z1 = griddata(points_wt, values_wt, (xi, yi), method='linear')
grid_wt_z2 = griddata(points_wt, values_wt, (xi, yi), method='cubic')

###########################################
### Plot the results
fig = plt.figure()

ax1 = plt.subplot(2, 2, 1)
contours_sa.plot(ax=ax1, column=contours_sa["height"], alpha=0.5, linewidth=1)
ax1.set_xlim([xmin, xmax])
ax1.set_ylim([ymin, ymax])
ax1.text(139, -25, "Original contours")


ax2 = plt.subplot(2, 2, 2)
plt.imshow(np.flipud(grid_wt_z0.T), extent=(xmin,xmax,ymin,ymax), origin='lower')
ax2.text(139, -25, "Nearest interp")


ax3 = plt.subplot(2, 2, 3)
plt.imshow(np.flipud(grid_wt_z1.T), extent=(xmin,xmax,ymin,ymax), origin='lower')
ax3.text(139, -25, "Linear interp")


ax4 = plt.subplot(2, 2, 4)
plt.imshow(np.flipud(grid_wt_z2.T), extent=(xmin,xmax,ymin,ymax), origin='lower')
ax4.text(139, -25, "Cubic interp")

plt.savefig(os.path.join("figures", "wt_interps"), dpi=300)
###########################################

###########################################
### Exactly the same plot but with contours
fig = plt.figure()

ax1 = plt.subplot(2, 2, 1)
contours_sa.plot(ax=ax1, column=contours_sa["height"], alpha=0.5, linewidth=1)
ax1.set_xlim([xmin, xmax])
ax1.set_ylim([ymin, ymax])
#ax1.text(139, -25, "Original contours")


ax2 = plt.subplot(2, 2, 2)
plt.imshow(np.flipud(grid_wt_z0.T), extent=(xmin,xmax,ymin,ymax), origin='lower')
#ax2.text(139, -25, "Nearest interp")
contours_sa.plot(ax=ax2, color='k', linewidth=1)


ax3 = plt.subplot(2, 2, 3)
plt.imshow(np.flipud(grid_wt_z1.T), extent=(xmin,xmax,ymin,ymax), origin='lower')
contours_sa.plot(ax=ax3, color='k', linewidth=1)
#ax3.text(139, -25, "Linear interp")


ax4 = plt.subplot(2, 2, 4)
plt.imshow(np.flipud(grid_wt_z2.T), extent=(xmin,xmax,ymin,ymax), origin='lower')
contours_sa.plot(ax=ax4, color='k', linewidth=1)
#ax4.text(139, -25, "Cubic interp")

plt.savefig(os.path.join("figures", "wt_interps_cont"), dpi=300)
###########################################
