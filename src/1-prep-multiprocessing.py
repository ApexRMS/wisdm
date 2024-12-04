## ---------------------------------
## wisdm - prep multiprocessing
## ApexRMS, March 2024
## ---------------------------------

# built under Python version 3.11.0 & SyncroSim version 3.0.0
# This script generates a spatial multiprocessing (smp) grid from a 
# template raster -- here specificaly set up to generaet smp grids 
# for conus map zones

#%%
## Source dependencies --------------------------------------------------------

# Set up environment and load helper functions
import os
import sys

# import modules
import rasterio
import pysyncrosim as ps
from rasterio import enums
import rioxarray
import xarray
import numpy as np
import math

import dask
from dask.distributed import Client

#%% Set progress bar ---------------------------------------------------------

steps = 6
ps.environment.progress_bar(report_type = "begin", total_steps = steps)

#%% Connect to SyncroSim library --------------------------------------------

# Load current scenario
myScenario = ps.Scenario()  
myLibrary = myScenario.library

# Create a temporary folder for storing rasters
# ssimTempDir = myLibrary.info["Value"][myLibrary.info.Property == "Temporary files:"].item() 
ssimTempDir = ps.runtime_temp_folder(os.path.join("DataTransfer", "Scenario-" + str(myScenario.sid)))

# Load datasheets
networkSheet = myScenario.library.datasheets("wisdm_Network")
templateRasterSheet = myScenario.datasheets("wisdm_TemplateRaster", show_full_paths=True)
multiprocessingSheet = myScenario.datasheets("core_Multiprocessing")

spatialMultiprocessingSheet = myScenario.datasheets("core_SpatialMultiprocessing")

# update progress bar
ps.environment.progress_bar()

#%% Set up dask client ------------------------------------------------------
if multiprocessingSheet.EnableMultiprocessing.item() == "Yes":
    num_threads = multiprocessingSheet.MaximumJobs.item()
else:
    num_threads = 1

# Note: Follow link in output to view progress
dask.config.set(**{'temporary-directory': os.path.join(ssimTempDir, 'dask-worker-space')})
client = Client(threads_per_worker = num_threads, n_workers = 1, processes=False)
# client

#%%
## Functions -----------------------------------------------------------------
# Resample a raster
def resample_grid(input_grid_filepath, template_grid_filepath, output_filename, output_directory):
    with rasterio.open(input_grid_filepath) as src:
      with rasterio.open(template_grid_filepath) as template:

        # Resample data to target shape
        updated_grid = src.read(
            out_shape = (
                src.count,
                template.shape[0],
                template.shape[1]
            ),
            resampling = enums.Resampling.nearest
        )

        out_meta = src.meta
        out_meta.update({
            "height": updated_grid.shape[1],
            "width": updated_grid.shape[2],
            "transform": template.transform})

        resampled_grid_filepath = os.path.join(output_directory, output_filename)

      with rasterio.open(resampled_grid_filepath, "w", **out_meta) as dst:
          dst.write(updated_grid)
      
      return resampled_grid_filepath

#%% Load data ---------------------------------------------------------------

# Load template raster
templatePath = templateRasterSheet.RasterFilePath.item()
templateRaster = rioxarray.open_rasterio(templatePath)

# Set desired number of cells per tile
if templateRasterSheet.TileCount.isnull().item():
    if templateRaster.size <= 10000:
        raise ValueError("Template raster has only " + str(templateRaster.size) + " pixels. Multiprocessing is not required.")
    elif templateRaster.size <= 10000:
        numTiles = templateRaster.size/10000
    elif templateRaster.size <= 100000:
        numTiles = templateRaster.size/100000
    elif templateRaster.size <= 1e6:
        numTiles = templateRaster.size/1e6
    elif templateRaster.size <= 1e7:
        numTiles = templateRaster.size/1e7
    elif templateRaster.size > 1e7:
        numTiles = templateRaster.size/1e7
        if numTiles > 20:
            numTiles = templateRaster.size/5e7
else:
    numTiles = templateRasterSheet.TileCount.item()   
    
tile_size = templateRaster.size/numTiles # 7500000

# Create contiguous tiles
# contig = True

# update progress bar
ps.environment.progress_bar()

#%% Build tiling grid --------------------------------------------------------

# Calculate ncol and nrow for one tile 
tile_dimension = math.ceil(math.sqrt((templateRaster.size / tile_size)))

# Define coords of smp grid ----
x_left = templateRaster.coords['x'][0].item()
x_right = templateRaster.coords['x'][-1].item()
x_coords = np.linspace(x_left, x_right, num = tile_dimension)

y_left = templateRaster.coords['y'][0].item()
y_right = templateRaster.coords['y'][-1].item()
y_coords = np.linspace(y_left, y_right, num = tile_dimension)

coords = {'band': [1], 'x': x_coords, 'y': y_coords}

#%% Create smp grid ----------------------------------------------------------
small_grid = xarray.DataArray(
    data = np.linspace(1, pow(tile_dimension, 2), num = pow(tile_dimension, 2)).reshape(1, tile_dimension, tile_dimension),
    dims = templateRaster.dims,
    coords = coords,
    attrs = templateRaster.attrs)

small_grid.rio.write_crs(templateRaster.rio.crs)

# Write to disk
with rasterio.open(templatePath) as src:

  # Scale image transform
  transform = src.transform * src.transform.scale(
      (src.width / small_grid.shape[-1]),
      (src.height / small_grid.shape[-2])
  )

  out_meta = src.meta
  out_meta.update({
     "dtype": 'int16',
     "nodata": -9999,
     "height": small_grid.shape[1],
     "width": small_grid.shape[2],
     "transform": transform})
  
  # profile = src.profile

  grid_filepath = os.path.join(ssimTempDir, "smp-grid-temp.tif")

  with rasterio.open(grid_filepath, "w", **out_meta) as dst:
      dst.write(small_grid)

# Resample to match primary stratum resolution
resample_grid(grid_filepath, templatePath, "smp-grid.tif", ssimTempDir)

# update progress bar
ps.environment.progress_bar()

#%% Mask to analysis area ----------------------------------------------------

mask = templateRaster
nodata_value = templateRaster.rio.nodata
if np.isnan(nodata_value):
   mask = ~np.isnan(mask)*1
else:
   mask = (mask != nodata_value)*1

smp_grid = rioxarray.open_rasterio(os.path.join(ssimTempDir, "smp-grid.tif"))

smp_grid = smp_grid.astype(np.int16) * mask

# Set output filename of smp grid
oldVals = np.unique(smp_grid)
num_tiles = len(oldVals)-1
output_filename = "smpGrid-" + str(num_tiles) + "-" + str(int(tile_size/1e3)) + "k.tif"

# update template raster datasheet
templateRasterSheet.TileCount = [num_tiles]
myScenario.save_datasheet(name="wisdm_TemplateRaster", data=templateRasterSheet)

# update progress bar
ps.environment.progress_bar()

#%%  Reclassify tiles --------------------------------------------------------

newVals = list(range(0,num_tiles+1))

lookups = list(zip(oldVals, newVals))
idx, val = np.asarray(lookups).T
lookup_array = np.zeros(idx.max() + 1)
lookup_array[idx] = val

smp_grid = lookup_array[smp_grid]

# Update 0 to NA Value
smp_grid[smp_grid == 0] = -9999

# Convert to int
smp_grid = smp_grid.astype(np.int16)

#%% Save sampling grid -------------------------------------------------------

with rasterio.open(os.path.join(ssimTempDir, "smp-grid.tif")) as src:
    out_meta = src.meta    
with rasterio.open(os.path.join(ssimTempDir, output_filename), "w", compress='deflate', **out_meta) as dst:
    dst.write(smp_grid)

# Save to sample grid to scenario
# add a row to the datasheet
spatialMultiprocessingSheet["MaskFileName"] = [os.path.join(ssimTempDir, output_filename)]
myScenario.save_datasheet(name="core_SpatialMultiprocessing", data=spatialMultiprocessingSheet) 

# update progress bar
ps.environment.progress_bar(report_type = "end")
# %%
