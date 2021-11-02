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
from scipy import stats
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

# Convert to utm to get the size in metres

sa_gdf_utm = sa_gdf.to_crs(utm)
sa_gdf_utm.head()


coordinates_array = np.asarray(sa_gdf_utm.geometry[0].exterior.coords)[:-1]# Note that the first and last points are the same

# Check out the coordinates and see the warping
sa_gdf_utm.plot()
plt.plot(coordinates_array[0][0], coordinates_array[0][1], "m*", label="index=0")
plt.plot(coordinates_array[1][0], coordinates_array[1][1], "k.", label="index=1")
plt.plot(coordinates_array[2][0], coordinates_array[2][1], "g^", label="index=2")
plt.plot(coordinates_array[3][0], coordinates_array[3][1], "yP", label="index=3")
plt.legend()
axes = plt.gca()
plt.xlabel("Longitude")
plt.ylabel("Latitude")

latitudes = []
longitudes = []
for i, crd in enumerate(coordinates_array):
    latitudes.append(coordinates_array[i][1])
    longitudes.append(coordinates_array[i][0])
    

# What is the error from the warped shape box??
lat_error_1_km = abs(latitudes[0] - latitudes[3])/1000    
lat_error_2_km = abs(latitudes[1] - latitudes[2])/1000

long_error_1_km = abs(longitudes[2] - longitudes[3])/1000    
long_error_2_km = abs(longitudes[0] - longitudes[1])/1000     

# There is substantial error...ok so it is an estimate? Use the largest values

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
    
# Calculate error as proportion of size (max error from conversion to utm)

error_lat_scaled_percent = max(lat_error_1_km, lat_error_2_km)*100000/sn_extent
error_long_scaled_percent = max(long_error_1_km, long_error_2_km)*100000/ew_extent

# Get the size of each cell

ns_extent = sa_gdf_utm[3] - sa_gdf_utm[1] # north-south
ew_extent = sa_gdf_utm[2] - sa_gdf_utm[0]

# This is in degrees, I need it in metres.


delr = sn_extent/nrows # Spacing along rows, i.e. this is the width of the columns
delc = ew_extent/ncols # Spacing along columns, i.e. this is the width of the rows

print("nrows: %s, ncols: %s" %(nrows,ncols))

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
show(dem_projected, ax=ax1, cmap="terrain", alpha=0.5)
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

#### ! Double check these are the right way around
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

#------------------------------------------------------------------------------
# Making comparisons between water table contours and the DEM

# But are they in different formats? One has lat long one has grid based
# Start with just comparing to z0 then add in others

wt_raster_optns = [grid_wt_z0,
                   grid_wt_z1,
                   grid_wt_z2]

wt_raster_idx = 0

for wt_raster_idx in range(len(wt_raster_optns)):
    
    #choose wt raster to work with
    wt_raster = wt_raster_optns[wt_raster_idx].T #grid_wt_z0.T
    
    print("Starting loop with wt raster: z%i" %wt_raster_idx)
    
    print(type(dem))       # numpy.ndarray
    print(type(wt_raster)) # numpy.ndarray
    
    print(dem.shape)
    print(wt_raster.shape)
    
    print("DEM max: %2.2f, Water table max: %2.2f " %(np.nanmax(dem), np.nanmax(wt_raster)))
    print("DEM min: %2.2f, Water table min: %2.2f " %(np.nanmin(dem), np.nanmin(wt_raster)))
    print("DEM average: %2.2f, Water table average: %2.2f " %(np.nanmean(dem), np.nanmean(wt_raster)))
    print("DEM st deviation: %2.2f, Water table st deviation: %2.2f " %(np.nanstd(dem), np.nanstd(wt_raster)))
    
    ###########################################
    plt.figure()
    ax1 = plt.subplot(2, 1, 1)
    plt.imshow(wt_raster)
    
    ax2 = plt.subplot(2, 1, 2)
    plt.imshow(dem, cmap="terrain", alpha=0.5)
    plt.savefig(os.path.join("figures", "dem_vs_wt_raster_z%i"%wt_raster_idx), 
                dpi=300)
    ###########################################
    
    #-----------------------------------------------------------------------------
    # Check if there are any places where the height of the DEM is below wt
    
    residual_dem_wt = dem - wt_raster
    
    print("Max of residuals: %2.2f" %np.nanmax(residual_dem_wt))
    print("Min of residuals: %2.2f" %np.nanmin(residual_dem_wt))
    print("Mean of residuals: %2.2f" %np.nanmean(residual_dem_wt))
    print("Standard deviation of residuals: %2.2f" %np.nanstd(residual_dem_wt))
    
    max_diff = max(abs(np.nanmin(residual_dem_wt)), 
                   abs(np.nanmax(residual_dem_wt)))
    
    ###########################################
    plt.figure()
    ax1 = plt.subplot(1, 1, 1)
    img1 = plt.imshow(residual_dem_wt, cmap="PiYG", vmin=-max_diff, vmax=max_diff)
    
    cb1 = plt.colorbar(img1)  
    cb1.set_label('DEM - wt (z%i) [m]' %wt_raster_idx)
    
    plt.savefig(os.path.join("figures", "diff_dem_wt_%i" %wt_raster_idx), dpi=300)
    ###########################################
    
    # Places where residuals are < 0, i.e. the DEM is higher than the wt
    
    wt_abv_dem_idx = np.where(residual_dem_wt > 0) 
    wt_abv_dem_values = residual_dem_wt[wt_abv_dem_idx] # This gets me the values
    wt_abv_dem = residual_dem_wt > 0 # This shows results as a Boolean array
    
    ###########################################
    plt.figure()
    plt.suptitle("Is wt above DEM? wt: z%i" %wt_raster_idx)
    ax1 = plt.subplot(1, 1, 1)
    img1 = plt.imshow(wt_abv_dem, cmap="binary")
    
    plt.savefig(os.path.join("figures", "wt_abv_dem_z%i" %wt_raster_idx), dpi=300)
    ###########################################

#------------------------------------------------------------------------------  
# Regression analysis between the elevation and the water table

    # Flatten arrays
    dem_flat        = dem.flatten()
    wt_raster_flat  = wt_raster.flatten()
    
    dem_flat.shape
    wt_raster_flat.shape
    
    # Create mask of locations with nans. ~ means "is not" (only in np)
    mask = ~np.isnan(dem_flat) & ~np.isnan(wt_raster_flat)
    
    res = stats.linregress(dem_flat[mask], wt_raster_flat[mask])
    res.slope    
    r2 = (res.rvalue)**2
    
    ###########################################
    plt.figure()
    plt.plot(dem_flat, wt_raster_flat, 'o', label='dem vs wt: z%i' %wt_raster_idx)
    plt.plot(dem_flat, res.intercept + res.slope*dem_flat, 'r', label='fitted line')
    plt.legend()
    plt.text(110, 70, "R2 = %2.2f" %r2)
    
    plt.savefig(os.path.join("figures", "wt_dem_linearregr_z%i" %wt_raster_idx), dpi=300)
    ###########################################
    
    #------------------------------------------------------------------------------
    # Calculating the gradient of the wter table map - to determine a flux
    
for wt_raster_idx in range(len(wt_raster_optns)):
    
    #choose wt raster to work with
    wt_raster = wt_raster_optns[wt_raster_idx].T #grid_wt_z0.T
    
    print("Starting loop with wt raster: z%i" %wt_raster_idx)
    
    gradient_wt = np.gradient(wt_raster)  
    
    max_grad = max(abs(np.nanmax(gradient_wt[0])), abs(np.nanmax(gradient_wt[1])))
    min_grad = min(np.nanmin(gradient_wt[0]), np.nanmin(gradient_wt[1]))
    vvalue = max(abs(max_grad), abs(min_grad))
    
    ###########################################
    # Plot the gradients as a map
    plt.figure()
    plt.suptitle("Gradients of wt: z%i" %wt_raster_idx)
    ax1 = plt.subplot(1, 2, 1)
    img1 = plt.imshow(gradient_wt[0], cmap="gist_ncar") # , vmin = vvalue, vmax=-vvalue, 
    ax1.text(500, 500, "x direction")
    cb1 = plt.colorbar(img1)  
    ax2 = plt.subplot(1, 2, 2)
    img2 = plt.imshow(gradient_wt[1], cmap="gist_ncar") # , vmin = vvalue, vmax=-vvalue
    ax2.text(500, 500, "y direction")
    cb2 = plt.colorbar(img2)  
    ax1.axis("off")
    ax2.axis("off")
    
    plt.savefig(os.path.join("figures", "wt_grad_z%i" %wt_raster_idx), dpi=300)
    ###########################################
    # there are lots of very low values, some larger values obscuring data
    # Check out the data more with a histogram
    
    mask_x = ~np.isnan(gradient_wt[0].flatten())
    mask_y = ~np.isnan(gradient_wt[1].flatten())
    
    ###########################################
#    # Box plots of the gradients
#    plt.figure()
#    a = plt.boxplot(gradient_wt[0].flatten()[mask_x], vert=True)
#    b = plt.boxplot(gradient_wt[1].flatten()[mask_y], vert=True)
#    
    ###########################################
    # There are a lot of "outliers"
    # Get the interquartlile range and just look at this
    
    iqr_x = stats.iqr(gradient_wt[0].flatten()[mask_x])
    iqr_y = stats.iqr(gradient_wt[1].flatten()[mask_y])

    p1 = 45
    p2 = 55

    pcile_1_x = np.percentile(gradient_wt[0].flatten()[mask_x], p1)
    pcile_2_x = np.percentile(gradient_wt[0].flatten()[mask_x], p1)
    
    pcile_1_y = np.percentile(gradient_wt[1].flatten()[mask_y], p2)
    pcile_2_y = np.percentile(gradient_wt[1].flatten()[mask_y], p2)
    
    ###########################################
    plt.figure()
    plt.suptitle("Histograms of wt gradients, x (b) & y(m) wt: z%i" %wt_raster_idx)
    
    ax1 = plt.subplot(2,1,1)
    plt.hist(gradient_wt[0].flatten()[mask_x], 
             density=True, range=(pcile_1_x, pcile_2_x),
             bins=100, color="b")  # density=False would make counts
    plt.ylabel('Probability')

    ax2 = plt.subplot(2,1,2)
    plt.hist(gradient_wt[1].flatten()[mask_y], 
             density=True, range=(pcile_1_y, pcile_2_y),
             bins=100, color="m")  # density=False would make counts
    plt.ylabel('Probability')
    plt.xlabel('Data')  
    plt.savefig(os.path.join("figures", "wt_grad_hists_z%i" %wt_raster_idx), dpi=300)

    ###########################################
    # Re-plot the gradients as a map
    plt.figure()
    plt.suptitle("Gradients of wt (narrow range): z%i" %wt_raster_idx)
    ax1 = plt.subplot(1, 2, 1)
    img1 = plt.imshow(gradient_wt[0], cmap="gist_ncar", vmin=pcile_1_x, vmax=pcile_2_x) 
    ax1.text(500, 500, "x direction")
    cb1 = plt.colorbar(img1)  
    ax2 = plt.subplot(1, 2, 2)
    img2 = plt.imshow(gradient_wt[1], cmap="gist_ncar", vmin=pcile_1_y, vmax=pcile_2_y)
    ax2.text(500, 500, "y direction")
    cb2 = plt.colorbar(img2)  
    ax1.axis("off")
    ax2.axis("off")
    
    plt.savefig(os.path.join("figures", "wt_grad_narrow_z%i" %wt_raster_idx), dpi=300)
    
    ###########################################
    # Trying to plot up flow lines
wt_raster_idx=0    
for wt_raster_idx in range(len(wt_raster_optns)):
    
    #choose wt raster to work with
    wt_raster = wt_raster_optns[wt_raster_idx].T #grid_wt_z0.T
    
    print("Starting loop with wt raster: z%i" %wt_raster_idx)
    
    gradient_wt = np.gradient(wt_raster)      
    # Get the x and y locations of the arrows
    # xi, yi --> These are the locations of every cell. Do I want every cell?
    len(xi)
    
    skipcells = 100
    X = xi[::skipcells, ::skipcells]
    Y = yi[::skipcells, ::skipcells]
    
    len(X)
    
    x_mult = 0.1
    y_mult = 0.1
    
    U = x_mult*(gradient_wt[0][::skipcells, ::skipcells])
    V = y_mult*(gradient_wt[1][::skipcells, ::skipcells])

    X.shape
    Y.shape
    U.shape
    V.shape 
    type(U)
    type(V)

    ### 
    plt.figure()
    plt.suptitle("Quiver map from gradients, wt: z%i" %wt_raster_idx)
    plt.imshow(wt_raster, alpha=0.5)

    plt.quiver(X, Y, U, V, 
               color='k', scale=5, headwidth=3,
               headlength=2,headaxislength=2, width=0.0025)
        
    plt.savefig(os.path.join("figures", "wt_grad_quiver_z%i" 
                             %wt_raster_idx), dpi=300)
    ###########################################
    
    # Make a streamplot from gradients (instead of quiver)
#    fig1, ax1 = plt.subplots()
#    show(dem, ax=ax1, cmap="terrain", alpha=0.5)
#    plt.suptitle("Streamplot map from gradients, wt: z%i" %wt_raster_idx)
#
#    ax1.streamplot(X, Y, U, V, density=0.5)
#    plt.savefig(os.path.join("figures", "wt_grad_streamplot_z%i" 
#                                 %wt_raster_idx), dpi=300)
    ###########################################
wt_raster_idx=0    
for wt_raster_idx in range(len(wt_raster_optns)):
    
    #choose wt raster to work with
    wt_raster = wt_raster_optns[wt_raster_idx].T #grid_wt_z0.T
    
    print("Starting loop with wt raster: z%i" %wt_raster_idx)
    
    gradient_wt = np.gradient(wt_raster)      
    
    # Make a streamplot directly from the wt raster data (not the gradient data)
    # For streamplot to work the U and V values are actually the 
    # VELOCITY at each point on the grid...so this is not correct
    
    X = xi
    Y = yi
    
    K_value = 1
    B = 100 # Aquifer thickness
    prsty = 0.1 # Porosity estimate    
    
    A_x = delc*B # Cross sectional area of flow in x direction
    A_y = delr*B # Cross sectional area of flow in y direction
    
    X_v = ((gradient_wt[0]/delr)*K_value*A_x)/prsty # Estimate of darcy's flux
    Y_v = ((gradient_wt[1]/delc)*K_value*A_y)/prsty # Estimate of darcy's flux

    xi.shape
    yi.shape
    X_v.shape
    Y_v.shape
    
    
    
    ###########################################
    # Plot the Darcy's velocity (q)
    fig = plt.figure()
    plt.suptitle("q map from wt: z%i, K: %2.2f, B: %i, porosity: %2.2f" 
                 %(wt_raster_idx, K_value, B, prsty))
    ax1 = plt.subplot(2,1,1)

    img1 = plt.imshow(X_v, alpha=0.5, vmin=-2, vmax=2)
    plt.colorbar(img1)  
    ax2 = plt.subplot(2,1,2)
    

    img2 = plt.imshow(Y_v, alpha=0.5, vmin=-2, vmax=2) 
    plt.colorbar(img2)
    plt.savefig(os.path.join("figures", "wt_q_z%i" 
                                 %wt_raster_idx), dpi=300)
    
    ###########################################
    fig1, ax1 = plt.subplots()
    plt.suptitle("Streamplot map from wt: z%i" %wt_raster_idx)
    img1 = ax1.imshow(wt_raster, alpha=0.5)
    
    ax1.streamplot(xi, yi, X_v, Y_v, density=0.8)
    cbar = plt.colorbar(img1)
    cbar.set_label('Hydraulic head [m]')

    plt.savefig(os.path.join("figures", "wt_streamplot_z%i" 
                                 %wt_raster_idx), dpi=300)
    ###########################################
#------------------------------------------------------------------------------
# Making comparisons between water table contours and the river network

# River raster is just 1 where there is a river, 0 where there isn't one...
# can I get stage data for the river and add to the raster?
    
riv_raster = riv_brn
wt_raster_idx = 0

for wt_raster_idx in range(len(wt_raster_optns)):
    
    #choose wt raster to work with 
    wt_raster = wt_raster_optns[wt_raster_idx].T #grid_wt_z0.T
    
    print("Starting loop with wt raster: z%i" %wt_raster_idx)
    
    print(type(riv_raster))       # numpy.ndarray
    print(type(wt_raster)) # numpy.ndarray
    
    print(riv_raster.shape)
    print(wt_raster.shape)
    
    print("Riv max: %2.2f, Water table max: %2.2f " %(np.nanmax(riv_raster), np.nanmax(wt_raster)))
    print("Riv min: %2.2f, Water table min: %2.2f " %(np.nanmin(riv_raster), np.nanmin(wt_raster)))


#------------------------------------------------------------------------------
    
# Adding streamlines to a vector field
    
#------------------------------------------------------------------------------
# Adding streamlines to the wt contour map

for wt_raster_idx in range(len(wt_raster_optns)):
    
