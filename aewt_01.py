# -*- coding: utf-8 -*-
"""
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

print(sys.version)
print('numpy version: {}'.format(np.__version__))
print('matplotlib version: {}'.format(mpl.__version__))
print('flopy version: {}'.format(flopy.__version__))

#=====DIRECTORIES==============================================================

overallName = "BH_05"
modelname = overallName
modelname_mt = modelname + "_MT3D" 

project_folder = r'C:\workspace\Proj4_BabyHope'

MODFLOW_folder      = r'C:\workspace\modflow_dbl\mfnwtdbl.exe' # double precision MODFLOW NWT
mfnwt_exe_name = "mfnwtdbl.exe" 

MT3D_USGS_folder      = r'C:\workspace\modflow_dbl\mt3dusgsdbl.exe' # Double precision MT3D-USGS
mt3d_version = 'mt3d-usgs' 
mt3d_exe_name = "mt3dusgsdbl.exe"

spatialFolder = r'C:\SMarshall_PhD\BabyHope\Sample_results\final_compiled'
spatialFolder2 = r'C:\SMarshall_PhD\Spatial\from_Saskia\East_Pilbara_Concept'

arc_output_fldr = r'C:\SMarshall_PhD\Spatial\Baby_Hope_Sampling_Location'

#------------------------------------------------------------------------------
# SET UP OVER-ARCHING FOLDER FOR SAVING OBSERVATION WELLS

os.getcwd()

this_model_folder = os.path.join(project_folder, modelname)

if not os.path.exists(this_model_folder):
    os.makedirs(this_model_folder)

# Set working directory
os.chdir(this_model_folder)

# outpath = 'data'
# if os.path.isdir(outpath):
#     shutil.rmtree(outpath)
# os.mkdir(outpath)

dataDirectory = os.path.join(this_model_folder, 'Data')
if not os.path.exists(dataDirectory):
    os.makedirs(dataDirectory)
       
figureDirectory = os.path.join(this_model_folder, 'Figures')
if not os.path.exists(figureDirectory):
    os.makedirs(figureDirectory)
    
shapeFilesDirectory = os.path.join(this_model_folder, 'shapeFiles')
if not os.path.exists(shapeFilesDirectory):
    os.makedirs(shapeFilesDirectory)
   
    
# Copy modflow to folder
this_model_modflow = os.path.join(this_model_folder, mfnwt_exe_name)
   
if not os.path.isfile(this_model_modflow):
    print("modflow nwt dbl not yet in folder")
    shutil.copy2(MODFLOW_folder, this_model_modflow)
else:
       print("modflow nwt dbl is already in folder")
       
# Copy mt3d to folder
this_model_mt3d = os.path.join(this_model_folder, mt3d_exe_name)
   
if not os.path.isfile(this_model_mt3d):
    print("mt3d usgs dbl not yet in folder")
    shutil.copy2(MT3D_USGS_folder, this_model_mt3d)
else:
    print("mt3d usgs dbl is already in folder")

#------------------------------------------------------------------------------
# DATA SET UP

# I am using the "Site_ID" as the index for the first dataFrame, the site_details
# But NO index for the other dataFrames, because I can have multiple samples from
# one site (well/bore). Ideally these should have a sample ID.

df_bh_site_details      = pd.read_excel(os.path.join(spatialFolder, 
                                        'BH_site_details.xlsx'),
                              sheet_name=0, header=0, index_col=0)

# =============================================================================
# MODEL GRID SET UP

nrow, ncol = 25, 100
xll, yll = 698948.82, 7447852.71 # origin of the model (lower left corner)
dxdy = 250 # grid spacing (in model units)
rot = 10 # rotation (positive counterclockwise)

# epsg code specifying coordinate reference system
model_epsg = 32750 # Got zone from here: https://spatialreference.org/ref/epsg/wgs-84-utm-zone-50s/

# row and column spacings
# (note that delc is column spacings along a row; delr the row spacings along a column)
delc = np.ones(nrow, dtype=float) * dxdy
delr = np.ones(ncol, dtype=float) * dxdy

sr = SpatialReference(delr=delr, delc=delc, xll=xll, yll=yll, rotation=rot, epsg=model_epsg)
sr

# Get information about grid coordinates
sr.xcentergrid, sr.ycentergrid

# Get cell vertices
sr.vertices

# Get grid bounds
sr.bounds

# Write shapefile of the grid
sr.write_shapefile(os.path.join(dataDirectory, 'grid.shp'))

# Write grid spec file (for PEST)
sr.write_gridSpec(os.path.join(dataDirectory, 'grid.spc'))

#The Jupyter notebook has info on how to interpolate head at some points, i.e. to 
# get well data (go through when I need this)

# Set up where rivers etc are?

# Other Shape files -----------------------------------------------------------
# Convert all to same projection (UTM)

surface_water = gpd.read_file(os.path.join(spatialFolder2, 'Hydrology_WA.shp'))
surface_water.crs # Check existing coordinate reference system
surface_water_utm = surface_water.to_crs(epsg=32750)
surface_water_utm.crs # It should now be in UTM

catchment = gpd.read_file(os.path.join(spatialFolder2, 'hydro_subcat.shp'))
catchment.crs
catchment_utm = catchment.to_crs(epsg=32750)

dykes = gpd.read_file(os.path.join(spatialFolder2, 'Hydraulic_Barriers.shp'))
dykes.crs
dykes_utm = dykes # Already in UTM

aquifer_boundaries = gpd.read_file(os.path.join(spatialFolder2, 'Aquifer_polyline.shp'))
aquifer_boundaries.crs
aquifer_boundaries_utm = aquifer_boundaries # Already in UTM

aquifer_bound_fill = gpd.read_file(os.path.join(spatialFolder2, 'Boundaries.shp'))
aquifer_bound_fill.crs
aquifer_bound_fill_utm = aquifer_bound_fill # Already in UTM

#=====MODEL DISCRETISATION=====================================================

## FIXED PARAMETERS --> Use an external script where all of this information is set...

Lx    = ncol*dxdy  # Length of the model sides in metres.
Ly    = nrow*dxdy  # Length of the model sides in metres.
delr  = 100.  # Spacing length across rows in metres.
delc  = 100.  # Spacing length across columns in metres.

maxThickness = 300
sHead = 300  # Starting head across aquifer.
sy    = 0.1 ## Specific yield
prsity = sy
ss    = sy/maxThickness ## Specific storage
laytyp = 1 ## Layer type (1 = convertible; 0 = confined)
length_simulation = 1
recharge_flux = 1.3699e-5 # m/d HIGH This is equivalent ot 5 mm per year.
hk_aquifer = 1.  # Hydraulic conductvity along rows (m/day) It worked ar 0.09
vka   = hk_aquifer/10. # Vertical hydraulic conductivity

 #=====IBOUND PACKAGE==========================================================
    
    
# Setting up the real iBound array

grid_df         = pd.read_excel(os.path.join(arc_output_fldr, 'grid_from_gis.xls'))
boundr_grid_df  = pd.read_excel(os.path.join(arc_output_fldr, 'grid_boundary_from_gis.xls'))

"""
Note that 'grid_boundary_from_gis.xls' has two columns for "columns", I have combined these
into one column - there was no overlap. The old one unaltered is saves as: 
    "grid_boundary_from_gis_unaltered.xls" 
"""

grid_df.index = grid_df.node
boundr_grid_df.index = boundr_grid_df.node

# IBOUND
# I then imported grid file into ArcGIS and manually selected boundary layers
# (details here in my OneNote --> Useful Tips --> ArcGIS Making Grid)

grid_df["boundary_cell"] = ["N"]*len(grid_df.index)

# See which nodes are in both the grid df and boundary grid df.
for idx in grid_df.index:
    if (idx in boundr_grid_df.index) == True:
        grid_df.loc[idx, "boundary_cell"] = "Y"
    else:
        pass

for idx in grid_df.index:
    print(grid_df.loc[idx, "boundary_cell"])

# 0 = inactive cell 

ibound_with_bounds = np.ones((nrow, ncol), dtype=np.int32)

countibound = 1   # Use this to count through each cell - should correspond to the node in the dataframe
for ir in range(nrow):
    for ic in range(ncol):
        yes_or_no = grid_df.loc[countibound, "boundary_cell"]
        
        if yes_or_no == "Y":
            ibound_with_bounds[ir, ic] = 0
        else:
            pass
        countibound += 1

extent = (delr/2., Lx - delr/2., delc/2., Ly - delc/2.)
        
# Plot new ibound array
plt.figure()
plt.subplot(1, 1, 1, aspect='equal')    

plt.imshow((ibound_with_bounds[:, :]), extent=extent, cmap="PiYG")
cbr = plt.colorbar()
cbr.set_label('iBound New - With Boundaries')

#=====ARRAY WITH BOUNDARIES AND RIVER==========================================

river_grid_df = pd.read_excel(os.path.join(arc_output_fldr, 'river_grid_selection.xls'))
river_grid_df.index = river_grid_df.node

# - - - - 
# This won't be used for anything but I just want to make a new array with the river also

grid_df["river_cell"] = ["N"]*len(grid_df.index)

# See which nodes are in both the grid df and boundary grid df.
for idx in grid_df.index:
    if (idx in river_grid_df.index) == True:
        grid_df.loc[idx, "river_cell"] = "Y"
    else:
        pass
    
rivercell_value = 5

ibound_with_riv = np.copy(ibound_with_bounds)

countibound = 1   # Use this to count through each cell - should correspond to the node in the dataframe
for ir in range(nrow):
    for ic in range(ncol):
        yes_or_no = grid_df.loc[countibound, "river_cell"]
        
        if yes_or_no == "Y":
            ibound_with_riv[ir, ic] = rivercell_value
        else:
            pass
        countibound += 1

extent = (delr/2., Lx - delr/2., delc/2., Ly - delc/2.)
        
# Plot new ibound array
plt.figure()
plt.suptitle("iBound - With Boundaries & River cells")
plt.subplot(1, 1, 1, aspect='equal')    

plt.imshow((ibound_with_riv[:, :]), extent=extent, cmap="jet_r") #cubehelix_r 

#=====VALLEY SHAPED AQUIFER====================================================

# Assign different number of layers depending on the location along the rows, this
# is to emulate the valley shape depenging on the position of the river.

edgelayers = 10 # Each layer is 10 m each, 100 m is side depth
riverlayers = 30 # 300 is total depth


icol = 0
# Firstly group the array so it's separated into what is a boundary, valley or river
grouped_cols_v = {}

for icol in range(ncol):

    ar_col = ibound_with_riv[:, icol]
    
    grouped_cols_v[icol] = [(k, sum(1 for i in g)) for k,g in groupby(ar_col)]  

key = 0
# Now loop through this and make the row layers
all_row_layers_v = {}
all_row_tog_v = {} # Put them together

for key in grouped_cols_v:    
    
    all_row_layers_v[key] = []
    all_row_tog_v[key]    = []
    
    print(key)
    grouped_col = grouped_cols_v[key]
    
    # i = 4
    for i in range(len(grouped_col)):
        count = 0

        group = grouped_col[i]
        
        if (i == 0 and len(grouped_col) == 1): # Only one section - not sure what to do here?
            lays = []
                
            for ig in range(group[1]):
                lays.append(-50)
                
            all_row_layers_v[key].append(lays)   
            all_row_tog_v[key] = all_row_tog_v[key] + lays
        
        elif (i == 0 and len(grouped_col) > 1): # It's the first one but there is more than one
        
            if (group[0] == 1) and (grouped_col[i+1][0] == 0): # Before a boundary
                
                lays = []
                
                # Even number in group
                if group[1] % 2 == 0:
                    for ig in range(int(group[1]/2)):
                        depth = edgelayers+(ig*10)
                        lays.append(depth)
                    maxdepth=depth
                    for ig in range(int(group[1]/2)):
                        depth = maxdepth-(ig*10)
                        lays.append(depth)
                        
                # Odd number in group
                else:
                    for ig in range(int(group[1]/2)):
                        depth = edgelayers+(ig*10)
                        lays.append(depth)
                    maxdepth=depth
                    lays.append(maxdepth)
                    for ig in range(int(group[1]/2)):
                        depth = maxdepth-(ig*10)
                        lays.append(depth)
                
                all_row_layers_v[key].append(lays)
                all_row_tog_v[key] = all_row_tog_v[key] + lays
                    
    
            elif group[0] == 0: # Boundary cells, not sure what to do for depth here??
                lays = []
                
                for ig in range(group[1]):
                    lays.append(-10)
                
                all_row_layers_v[key].append(lays)
                all_row_tog_v[key] = all_row_tog_v[key] + lays
            
            elif (group[0] == 1) and (grouped_col[i+1][0] == 5): # Leading to a river 
                lays = []
                
                slope1 = (riverlayers-edgelayers)/group[1]
                
                rowadd1 = 0
            
                for ir in range(0, group[1]):
                
                    lays.append(edgelayers+rowadd1)
                
                    if (ir % 2) == 0:
                        rowadd1 = rowadd1 + np.floor(slope1)
                    else:
                        rowadd1 = rowadd1 + np.ceil(slope1)
                        
                all_row_layers_v[key].append(lays)
                all_row_tog_v[key] = all_row_tog_v[key] + lays
                
            elif (grouped_col[i][0] == 5): # At a river 
                lays = []
                
                for ig in range(group[1]):
                    lays.append(riverlayers)
        
                all_row_layers_v[key].append(lays)
                all_row_tog_v[key] = all_row_tog_v[key] + lays
    
                
            else:
                print("Check it out, there is some other category I haven't thought of for Key: %d; i: %d" %(key, i))
                
        
        elif i != (len(grouped_col)-1): # It's not the last one, but also not the first one also
        
            if group[0] == 0: # Boundary cells, not sure what to do for depth here??
                lays = []
                
                for ig in range(group[1]):
                    lays.append(-10)
                
                all_row_layers_v[key].append(lays)
                all_row_tog_v[key] = all_row_tog_v[key] + lays
            
            elif (group[0] == 1) and (grouped_col[i+1][0] == 5): # Leading to a river 
                lays = []
                
                slope1 = (riverlayers-edgelayers)/group[1]
                
                rowadd1 = 0
            
                for ir in range(0, group[1]):
                
                    lays.append(edgelayers+rowadd1)
                
                    if (ir % 2) == 0:
                        rowadd1 = rowadd1 + np.floor(slope1)
                    else:
                        rowadd1 = rowadd1 + np.ceil(slope1)
                        
                all_row_layers_v[key].append(lays)
                all_row_tog_v[key] = all_row_tog_v[key] + lays
    
                        
            elif (grouped_col[i][0] == 5): # At a river 
                lays = []
                
                for ig in range(group[1]):
                    lays.append(riverlayers)
        
                all_row_layers_v[key].append(lays)
                all_row_tog_v[key] = all_row_tog_v[key] + lays
    
        
            elif (group[0] == 1) and (grouped_col[i-1][0] == 5): # Leading away from a river 
                lays = []
                
                slope2 = (riverlayers-edgelayers)/group[1]
                
                rowtake2 = 1
            
                for ir in range(0, group[1]):
                
                    lays.append(riverlayers-rowtake2)
                
                    if (ir % 2) == 0:
                        rowtake2 = rowtake2 + np.floor(slope2)
                    else:
                        rowtake2 = rowtake2 + np.ceil(slope2)
                        
                all_row_layers_v[key].append(lays)
                all_row_tog_v[key] = all_row_tog_v[key] + lays
                
            elif (grouped_col[i-1][0] == 5 and grouped_col[i+1][0] == 5): # Between a river
                lays = []
                
                for ig in range(group[1]):
                    lays.append(riverlayers)
        
                all_row_layers_v[key].append(lays)
                all_row_tog_v[key] = all_row_tog_v[key] + lays
                
            elif (group[0] == 1 and grouped_col[i+1][0] == 0 and 
                  grouped_col[i-1][0] != 5): # Before a boundary, but not after a river
                
                lays = []
                
                # Even number in group
                if group[1] % 2 == 0:
                    for ig in range(int(group[1]/2)):
                        depth = edgelayers+(ig*10)
                        lays.append(depth)
                    maxdepth=depth
                    for ig in range(int(group[1]/2)):
                        depth = maxdepth-(ig*10)
                        lays.append(depth)
                        
                # Odd number in group
                else:
                    for ig in range(int(group[1]/2)):
                        depth = edgelayers+(ig*10)
                        lays.append(depth)
                    maxdepth=depth
                    lays.append(maxdepth)
                    for ig in range(int(group[1]/2)):
                        depth = maxdepth-(ig*10)
                        lays.append(depth)
                
                all_row_layers_v[key].append(lays)
                all_row_tog_v[key] = all_row_tog_v[key] + lays
                
                
            elif (group[0] == 1 and grouped_col[i+1][0] == 5 and 
                  grouped_col[i-1][0] != 0): # After a boundary, but not before a river
                
                lays = []
                
                # Even number in group
                if group[1] % 2 == 0:
                    for ig in range(int(group[1]/2)):
                        depth = edgelayers+(ig*10)
                        lays.append(depth)
                    maxdepth=depth
                    for ig in range(int(group[1]/2)):
                        depth = maxdepth-(ig*10)
                        lays.append(depth)
                        
                # Odd number in group
                else:
                    for ig in range(int(group[1]/2)):
                        depth = edgelayers+(ig*10)
                        lays.append(depth)
                    maxdepth=depth
                    lays.append(maxdepth)
                    for ig in range(int(group[1]/2)):
                        depth = maxdepth-(ig*10)
                        lays.append(depth)
                
                all_row_layers_v[key].append(lays)
                all_row_tog_v[key] = all_row_tog_v[key] + lays
                
            else:
                print("Check it out, there is some other category I haven't thought of for Key: %d; i: %d" %(key, i))
                
        else: # It is the last one
        
            if (group[0] == 1) and (grouped_col[i-1][0] == 0): # After a boundary
                
                lays = []
                
                # Even number in group
                if group[1] % 2 == 0:
                    for ig in range(int(group[1]/2)):
                        depth = edgelayers+(ig*10)
                        lays.append(depth)
                    maxdepth=depth
                    for ig in range(int(group[1]/2)):
                        depth = maxdepth-(ig*10)
                        lays.append(depth)
                        
                # Odd number in group
                else:
                    for ig in range(int(group[1]/2)):
                        depth = edgelayers+(ig*10)
                        lays.append(depth)
                    maxdepth=depth
                    lays.append(maxdepth)
                    for ig in range(int(group[1]/2)):
                        depth = maxdepth-(ig*10)
                        lays.append(depth)
                
                all_row_layers_v[key].append(lays)
                all_row_tog_v[key] = all_row_tog_v[key] + lays
                
            elif (group[0] == 1) and (grouped_col[i-1][0] == 5): # Leading away from a river 
                lays = []
                
                slope2 = (riverlayers-edgelayers)/group[1]
                
                rowtake2 = 1
            
                for ir in range(0, group[1]):
                
                    lays.append(riverlayers-rowtake2)
                
                    if (ir % 2) == 0:
                        rowtake2 = rowtake2 + np.floor(slope2)
                    else:
                        rowtake2 = rowtake2 + np.ceil(slope2)
                        
                all_row_layers_v[key].append(lays)
                all_row_tog_v[key] = all_row_tog_v[key] + lays
                
            elif group[0] == 0: # Boundary cells, not sure what to do for depth here??
                lays = []
                
                for ig in range(group[1]):
                    lays.append(-10)
                
                all_row_layers_v[key].append(lays)
                all_row_tog_v[key] = all_row_tog_v[key] + lays

                
            else:
                print("Check it out, there is some other category I haven't thought of for Key: %d; i: %d" %(key, i))

'''

I think some have layers less than 10 -- > I need to look through that.

'''

#=====VALLEY SHAPED AQUIFER====================================================

# Assign different number of layers depending on the location along the rows, this
# is to emulate the valley shape depenging on the position of the river.

edgelayers = 10 # Each layer is 10 m each, 100 m is side depth
riverlayers = 30 # 300 is total depth

all_row_layers = []

ic = 0
for ic in range(ncol):
    
    colval = ic+1 # ArcGIS has 1-based numbering
    row_layers = []
    
# Is there a river depth for this? I.e. what is the max column for the river because it then veers north
    if colval <= (np.max(river_grid_df["column"])):
        
        # Now make the river cell the deepest one. rr = river row
        df  = river_grid_df[river_grid_df["column"]==colval]
        
        nRiverRows = len(df)
        
        rr1 = df.iloc[0, 2]   # Get the first river row
        rr2 = df.iloc[-1, 2]  # Get the last river row
        
        sidelen1 = rr1 - 0  # In number of cells
        sidelen2 = 25 - rr2 # In number of cells
        
        slope1 = (riverlayers-edgelayers)/(sidelen1-1)
        slope2 = (riverlayers-edgelayers)/(sidelen2-1)
        
        # Now I should make all of the row layers for one column
        
        # Do one side leading up to the river
        rowadd1 = 0
        
        for ir in range(0, (rr1-1)):
            
            row_layers.append(edgelayers+rowadd1)
            
            if (ir % 2) == 0:
                rowadd1 = rowadd1 + np.floor(slope1)
            else:
                rowadd1 = rowadd1 + np.ceil(slope1)
        
        # Add river rows
        for i in range(nRiverRows):
            row_layers.append(riverlayers)
        
        # Now do the other side leading away from the river
        rowtake2 = 1
        
        for ir in range(rr2, nrow):
            
            row_layers.append(riverlayers-rowtake2)
            
            if (ir % 2) == 0:
                rowtake2 = rowtake2 + np.floor(slope2)
            else:
                rowtake2 = rowtake2 + np.ceil(slope2)
                
        # Check the length, it should be the same as the number of rows
        if len(row_layers) == nrow:
            print("Correct same as number of layers")
        else:
            print("You have an issue it is not the same as the number of layers")
            if len(row_layers) > nrow:
                row_layers = row_layers[0:nrow]
            elif len(row_layers) < nrow:
                
                while len(row_layers) < nrow:
                    for i in range(nrow-len(row_layers)):
                        row_layers.append(edgelayers)
        
        
        # Check none are more than 30 
        for i in range(len(row_layers)):
            if row_layers[i] > 30:
                row_layers[i] = 30
            elif row_layers[i] < 10:
                row_layers[i] = 10
            else:
                print("All row layers are between 10 and 30")
    
    else: # I.e. there is no river anymore
    
        for ir in range(nrow):
            row_layers.append(riverlayers)
                
    # Done making this row layers for the column
    all_row_layers.append(row_layers)
        
        
        
# side1_height = float(ztop)

nlay  = 1     # Number of layers.
ztop  = 300.  # Top of model.
zbot  = 0.    # Base of model.
delv  = (ztop - zbot) / nlay # Length of cells in z direction.



    



botm  = np.linspace(ztop, zbot, nlay + 1) 

#=====BOUNDARY CONDITIONS======================================================
    
stageleft = float(sHead) # to have a head gradient of 0.001
stageright = float(sHead-100)
bound_sp1 = []
for il in range(nlay):
    condleft  = hk_aquifer * (stageleft - zbot) * delc
    condright = hk_aquifer * (stageright - zbot) * delc
    for ir in range(nrow):
        bound_sp1.append([il, ir, 0, stageleft, condleft])
        bound_sp1.append([il, ir, ncol - 1, stageright, condright])
print('Adding ', len(bound_sp1), 'GHBs for stress period 1.')
ghb_spd = {0: bound_sp1}    

#=====STARTING HEADS===========================================================


if nlay == 1:
    # ibound = np.ones((nrow, ncol), dtype=np.int32)
    strt   = sHead * np.ones((nrow, ncol), dtype=np.float32)
else:
    # ibound = np.ones((nlay, nrow, ncol), dtype=np.int32)
    strt   = sHead * np.ones((nlay, nrow, ncol), dtype=np.float32)
    
#=====MAKE HYDRAULIC CONDUCTIVITY VECTOR=======================================
       
hk_array_barrier = hk_aquifer * np.ones((nlay, nrow, ncol), dtype=np.float32)   

#==== RECHARGE CHARACTERISTICS ================================================

number_of_recharge_scenarios = 1

# Recharge array 1 = diffuse recharge across the whole aquifer in equal amounts.
recharge_array = recharge_flux * np.ones((nrow, 
                 ncol), dtype=np.float32)

#==============================================================================
#=====RUNNING THE FLOW MODEL===================================================  
#==============================================================================
print(modelname)
modelname_mt = modelname_mt
print(modelname_mt)

mf = flopy.modflow.Modflow(modelname, exe_name=mfnwt_exe_name, version='mfnwt') # MODFLOW-NWT
#mf = flopy.modflow.Modflow(modelname, exe_name=MODFLOW_folder) # MODFLOW-2005

nper = 1 # Number of model stress periods
perlen = [length_simulation] 

nstp = [1] # Number of time steps in each stress period
tsmult = [1.0]
steady = [True] # True = steady state

dis = flopy.modflow.ModflowDis(mf, nlay, nrow, 
                               ncol, delr=delr, 
                               delc=delc, top=ztop, 
                               botm=botm[1:],
                               tsmult = tsmult,
                               nper=nper, perlen=perlen, nstp=nstp, 
                               steady=steady)

print("Discretisation Module Set Up")

# Spatial reference

mf.sr = SpatialReference(delr=mf.dis.delr, delc=mf.dis.delc, 
                         xll=xll, yll=yll, rotation=rot, epsg=model_epsg)
mf.sr

#        # UPSTREAM WEIGHTING PACKAGE

hk_array = hk_array_barrier

uwp = flopy.modflow.mfupw.ModflowUpw(mf, hk=hk_array, vka=vka, 
                                     sy=sy, 
                                     ss=ss, 
                                     laytyp=laytyp) # MODFLOW- NWT

print("Upstream Weighting Package Set Up")

## NEWTON SOLVER
nwt = flopy.modflow.mfnwt.ModflowNwt(mf, headtol=1e-03, fluxtol=0.05, maxiterout=100, 
            thickfact=1e-07, linmeth=2, iprnwt=1, ibotav=1, options='SPECIFIED', 
            Continue=True, dbdtheta=0.9, dbdkappa=1e-05, dbdgamma=0.0, momfact=0.1, 
            backflag=1, maxbackiter=30, backtol=1.05, backreduce=0.9, maxitinner=50, 
            ilumethod=2, levfill=5, stoptol=1e-10, msdr=15, iacl=2, norder=1, level=3, 
            north=7, iredsys=1, rrctols=0.0, idroptol=1, epsrn=0.0001, 
            hclosexmd=0.0001, mxiterxmd=200) # MODFLOW- NWT

print("Newton Solver Set Up")

# --- Recharge --- #
        
# RECHARGE (RCH) PACKAGE

rch = flopy.modflow.ModflowRch(mf, rech=recharge_array)

# GENERAL HEAD BOUNDARY (GHB) PACKAGE 
ghb = flopy.modflow.ModflowGhb(mf, stress_period_data=ghb_spd)

print("General Head Boundary Package Set Up")   
             
# OUTPUT CONTROL (OC) PACKAGE
                       
spd = {}
for strsp in range(nper):
    tmstp = nstp[strsp]
    for time_ in range(tmstp):
        spd[strsp, time_] = ['save head', 'print budget', 
                               'save budget'] 
                               
oc = flopy.modflow.ModflowOc(mf, stress_period_data=spd)

print("Output Control Set Up")

#bas = flopy.modflow.ModflowBas(mf, ibound=ibound, strt=300)   

bas = flopy.modflow.ModflowBas(mf, ibound=ibound_with_bounds, strt=300., 
                               ifrefm=True, ixsec=False, 
                               hnoflo=-999.99)     

print("Basic Package Set Up") 

# --- LINKING FILE --- #

# Set up the linking file for output from flow model to be used in transport model
lmt = flopy.modflow.ModflowLmt(mf, output_file_name= (modelname + str('_mt3dLink.ftl')),
                               output_file_header='extended',
                               output_file_format='unformatted')   
                                
print("Link-Mass Transport Package Set Up") 
  
mf.write_input() # Write the model input files
print("input written")

# What's the best way to run the command?
subprocess.run([mfnwt_exe_name, modelname]) # Run command 
   
#==============================================================================
#=====RUNNING THE TRANSPORT MODEL==============================================  
#==============================================================================

#==== TRACER PARAMETERS =======================================================

# Instantiate MT3D-USGS object in flopy

mt = mt3d.Mt3dms(modflowmodel=mf, modelname=modelname_mt, model_ws=this_model_folder, 
                 version=mt3d_version, namefile_ext='mtnam', exe_name=mt3d_exe_name,  
                 ftlfilename=(modelname + str('_mt3dLink.ftl'))) 

# ftlfree=True. This was on the online code but I can't find anywhere what it means.

#------------------------------------------------------------------------------
# BASIC TRANSPORT PACKAGE

ncomp = 1 # Total number of chemical species in simulation  
mcomp = 1 # Total number of "mobile" species 
sconc= np.zeros((nlay, nrow, ncol), dtype=np.float32)  # initial conc. = 0 yrs   
nprs = 0 # Flag indicating frequency of output. = 0: results not saved except 
# at end of simulation; > 0: saved at times specified by timprs; < 0: saved
# whenever the  number of transport steps is an even multiple of nprs.

# Instantiate basic transport (BTN) package for MT3D-USGS
btn = mt3d.Mt3dBtn(mt, nlay=nlay, nrow=nrow, ncol=ncol, nper=nper, perlen=perlen, 
                   ncomp=ncomp, mcomp=mcomp, sconc=sconc, prsity=prsity,
                   delr=delr, delc=delc, 
                   icbund=1, ifmtcn=-1, savucn=True,
                   nprs=nprs, nprobs=1, cinact=0, ssflag='SState', laycon=1) 

#------------------------------------------------------------------------------
# ADVECTION PACKAGE

# Instantiate advection (ADV) package for MT3D-USGS
adv = mt3d.Mt3dAdv(mt, mixelm=0)  

# mixelm = 0 means that I am using the standard finite-difference method with upstream or 
# central-in-stream weighting, depending on the value of NADVFD. 
# Other options, MOC (1); MMOC (2); HMOC (3); TVD (-1). 

#------------------------------------------------------------------------------
# GENERALISED CONJUGATE GRADIENT SOLVER for MT3D-USGS

# Instatiate generalized conjugate gradient solver (GCG) package for MT3D-USGS
gcg = mt3d.Mt3dGcg(mt, mxiter=30, iter1=50, isolve=2, accl=1, cclose=1e-06)

#------------------------------------------------------------------------------
# REACTION PACKAGE

rc1 = np.zeros((nlay, nrow, ncol), 
               dtype=np.float32)

rc1[:, :, :] = -1/365.25

isothm = 0      # 0 = no sorption
ireact = 100    # 100 = zeroth-order reaction option
rc1 = rc1       # First order reaction rate for diffolved phase of first species.

# Setting up the Reaction package
rct= mt3d.Mt3dRct(mt, isothm=isothm, ireact=ireact, igetsc=0, rc1=rc1) 

#------------------------------------------------------------------------------
# SOURCE-SINK MIXING PACKAGE

crch = np.zeros((nrow, ncol), 
                dtype=np.float32) # The age of recharging water is 0

itype = mt3d.Mt3dSsm.itype_dict()
#ssm_data = {}
#ssm_data[0] = [(9, 0, -1, 0., itype['CHD'])]
#
#for i in range(nrow):
#    ssm_data[0].append((9, i, -1, 0., itype['CHD']))
#

#Instantiate source-sink mixing (SSM) package for MT3D-USGS
ssm = mt3d.Mt3dSsm(mt, crch=crch) #, stress_period_data=ssm_data) 
      
#------------------------------------------------------------------------------
# DISPERSION PACKAGE 
al = 1.5 # The longitudinal dispersivity, default = 0.01
#        trpt = 0.1  #ratio of horizontal transverse: longitudinal dispersivity, default = 0.1
#        trpv = 0.01 #ratio of vertical transverse dispersivity: longitudinal dispersivity, default = 0.01
dmcoef = 1E-4 # Effective molecular diffusion coefficient (for water in my model), default = 1E-9
                # 9E-6 is a value that Lenny used --> but I think this is in m^2/s not m^2/d!

# Instantiate up the Dispersion package for MT3D-USGS
dsp = mt3d.Mt3dDsp(mt, al=al, dmcoef=dmcoef) 

#----- WRITING transport MODEL ------------------------------------------------
 
mt.write_input()
mt.write_name_file()

# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- -- --
# "Manual" changes to the input files
conc_filename_dissolved = str(modelname_mt) + ".ucn"
conc_filename_sorbed = str(modelname_mt) + "_S.ucn"
mass_filename = str(modelname_mt) + ".mas"
cnf_filename = str(modelname_mt) + ".cnf"

##Add a line to the MT3D .mtnam file re-naming output UCN file + MAS,CNF
mt_name_file = modelname_mt + str('.mtnam')

namfile = open(mt_name_file, 'a')
namfile.write('data(binary) 201 %s \n' %(conc_filename_dissolved))
namfile.close()

namfile = open(mt_name_file, 'a')
namfile.write('data(binary) 301 %s \n' %(conc_filename_sorbed))
namfile.close()
      
namfile = open(mt_name_file, 'a')
namfile.write('data 601 %s \n' %(mass_filename))
namfile.close()

namfile = open(mt_name_file, 'a')
namfile.write('data 17 %s \n' %(cnf_filename))
namfile.close()
        
###For USGS need to add DRYCELL keyword to the .btn file on line 3
mt_btn_file = modelname_mt + str('.btn')


btnfile = open(mt_btn_file, "r")
contents = btnfile.readlines()
btnfile.close()
contents.insert(2, "DRYCELL \n")
btnfile = open(mt_btn_file, "w")
contents = "".join(contents)
btnfile.write(contents)
btnfile.close()

# ADD SSTATE FLAG manually in correct spot

fpath = os.path.join(this_model_folder, mt_btn_file)
with open(mt_btn_file, 'r+') as f:
    lines = f.readlines()
    f.seek(0)
    f.truncate()

    a = length_simulation # '%.1E' % Decimal(length_simulation)
    b = str(a) + "         1         1 S"
    c = str(a) + "         1         1         SState"
 
    ##

    for line in lines:
        if b in line:
                        line = line.replace(b, c)
        f.write(line)

mt3d_namefile = modelname_mt + ".mtnam"

# Running the transport model
subprocess.run([mt3d_exe_name, mt3d_namefile]) # Run command 

#==============================================================================
#=====CREATE DATA OUTPUTS======================================================  
#==============================================================================
 
# CREATE HEADFILE OBJECTS  
    
headobj = bf.HeadFile(modelname + '.hds')
data = headobj.get_alldata()
times = headobj.get_times()
head_raster = headobj.get_data(totim=times[0]) # One time period: steady-state

# CONCENTRATION/AGE FILES - FROM TRANSPORT MODEL

directage_file = modelname_mt + str('.ucn') 

ucnobjc_direct_age = flopy.utils.binaryfile.UcnFile(directage_file, precision="double")

age_times = ucnobjc_direct_age.get_times()

age_raster = ucnobjc_direct_age.get_data(totim=age_times[0])

# Plotting options
head_cmap = 'coolwarm'
age_cmap = 'jet'

#==============================================================================
#=====PLOTTING HEAD DATA=======================================================  
#==============================================================================

# PLOT HEAD IN GENERAL FORM

hds = headobj.get_data(totim=times[0]) # One time period: steady-state

levels = np.linspace(0, 400, 51)

# Change head data in inactive cells to equal NaN instead of -999.9 which is output
# from MODFLOW

hds[hds == -999.99] = np.nan



# Plot steady state hydraulic heads

plt.figure()
plt.subplot(1, 1, 1, aspect='equal')    

plt.imshow(np.flipud(hds[0, :, :]), extent=extent, cmap=head_cmap) 

cbr = plt.colorbar()
cbr.ax.tick_params(labelsize=14) 
cbr.set_label('Head in mAHD', fontsize = 14)
cs = plt.contour((hds[0, :, :]), levels=levels, extent=extent, 
                           colors = "k") 
plt.clabel(cs, inline=1, fontsize=14, fmt='%1i')

# Plot location of observation wells
# plt.savefig(figureDirectory + "\head_modelled")

# Save the head data so I can look at them
np.savetxt("heads.csv", hds[0, :, :], delimiter=",")

#==============================================================================
# WELL LOCATIONS
#==============================================================================

mf.sr
sr.bounds
mf.sr.get_extent()

# Test one well location
well =     df_bh_site_details.index[1]
print(well)
easting =  df_bh_site_details.loc[well, "Easting"]
northing = df_bh_site_details.loc[well, "Northing"]

# Get interpolated head at well location
# ?? This doesn't seem to be working
head_well = mf.sr.interpolate(hds[0], ([northing, easting], [northing, easting]))
print(head_well)

# Plot well location
well_ij = sr.get_ij(easting, northing)
sr.get_vertices(well_ij[0], well_ij[1])

plt.figure()
plt.plot([well_ij[1]*delc], [well_ij[0]*delr], "k*")

# PLOT HEAD SPATIALLY REFERENCED
plt.figure()
ax1 = sr.plot_array(hds[0], cmap=head_cmap)
plt.plot(easting, northing, "k*")
cbr = plt.colorbar(ax1)
cbr.ax.tick_params(labelsize=14) 
cbr.set_label('Head in mAHD', fontsize = 14)

#==============================================================================
# PLOT AGE WITHOUT BARRIER
#==============================================================================

# Save the age data so I can look at them
np.savetxt("ages.csv", age_raster[0, :, :], delimiter=",")

# Remove inactive age values
age_raster[age_raster == -0.0] = np.nan



plt.figure()
plt.subplot(1, 1, 1, aspect='equal')    

plt.imshow((age_raster[0, :, :]), extent=extent, cmap=age_cmap)
cbr = plt.colorbar()
cbr.set_label('Age (years)')

# plt.savefig(figureDirectory + "\plotage_modelled")

# PLOT AGE SPATIALLY REFERENCED
plt.figure()
ax1 = sr.plot_array(age_raster[0], cmap=age_cmap)
plt.plot(easting, northing, "m*")
cbr = plt.colorbar(ax1)
cbr.ax.tick_params(labelsize=14) 
cbr.set_label('Age in yrs', fontsize = 14)

bl = [sr.bounds[0], sr.bounds[1]]
br = [sr.bounds[2], sr.bounds[1]]
tl = [sr.bounds[0], sr.bounds[3]]
tr = [sr.bounds[2], sr.bounds[3]]

fig = plt.figure()
axes = plt.subplot(1, 1, 1)
plt.plot(easting, northing, "m*")

axes.add_patch(Polygon([bl, br, tr, tl], 
                      facecolor="k", edgecolor="k", closed=True, lw=0,
                      fill=True))
plt.show()

# PLot model grid --> this is far too small I don't know why
mf.dis.sr
mf.dis.sr.epsg
mf.modelgrid.set_coord_info(xoff=sr.bounds[0], yoff=sr.bounds[1], angrot=rot, epsg=model_epsg)

fig = plt.figure(figsize=(8,8))
ax = fig.add_subplot(1, 1, 1, aspect='equal')
modelmap = flopy.plot.PlotMapView(model=mf)
linecollection = modelmap.plot_grid(linewidth=0.5, color='pink') # royalblue

for well in df_bh_site_details.index:
    plt.plot(df_bh_site_details.loc[well, "Easting"], df_bh_site_details.loc[well, "Northing"], "k*")

#Representation of model geometry with all the boundary conditions
fig = plt.figure(figsize=(8,8))

modelmap = flopy.plot.PlotMapView(model=mf)
linecollection = modelmap.plot_grid(linewidth=0.5, color='cyan')
# shpRiver = flopy.plot.plot_shapefile(surface_water_utm, facecolor='none')    


#==============================================================================
# Make new iBound Array

test_ibound_array = np.ones((nrow, ncol), dtype=np.int32)

countibound = 0
for ir in range(nrow):
    for ic in range(ncol):
        test_ibound_array[ir, ic] = countibound
        countibound += 1

# Plot original ibound array
plt.figure()
plt.subplot(1, 1, 1, aspect='equal')    

plt.imshow((test_ibound_array[:, :]), extent=extent, cmap="PiYG")
cbr = plt.colorbar()
cbr.set_label('iBound Test - coloured by cells')

#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#_#
