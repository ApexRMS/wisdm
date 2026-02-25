# ---------------------------------
# wisdm - prep multiprocessing
# ApexRMS, March 2024
# ---------------------------------

# built under Python version 3.11.0 & SyncroSim version 3.0.0
# This script generates a spatial multiprocessing (smp) grid from a
# template raster -- here specificaly set up to generaet smp grids
# for conus map zones

# Source dependencies --------------------------------------------------------

# IMPORTANT: setup_functions must be imported before any non-stdlib packages.
# It removes user site-packages on import and setupCondaEnv() configures
# conda DLL paths. Linters must not reorder these imports.  # noqa: E402
import os
import sys
import platform
from setup_functions import setupCondaEnv, checkGdalVersion, setupGdalProj  # noqa: E402

setupCondaEnv()
checkGdalVersion()

# Non-stdlib imports (must follow setup calls above)
import math  # noqa: E402
import numpy as np  # noqa: E402
import pyproj  # noqa: E402
import rasterio  # noqa: E402
import rioxarray  # noqa: E402
import xarray  # noqa: E402
from rasterio import enums  # noqa: E402
import pysyncrosim as ps  # noqa: E402

# Set progress bar ---------------------------------------------------------

steps = 6
ps.environment.update_run_log('1 - Prepare Multiprocessing => Begin')
ps.environment.progress_bar(report_type="begin", total_steps=steps)

# Connect to SyncroSim library --------------------------------------------

# Load current scenario
myScenario = ps.Scenario()
myLibrary = myScenario.library
mySession = ps.Session()

# Set GDAL/PROJ paths if using conda environment ---------------------------

setupGdalProj(myLibrary)

# Create a temporary folder for storing rasters
# ssimTempDir = myLibrary.info["Value"][myLibrary.info.Property == "Temporary files:"].item()
ssimTempDir = ps.runtime_temp_folder(os.path.join(
    "DataTransfer", "Scenario-" + str(myScenario.sid)))

# Load datasheets
networkSheet = myScenario.library.datasheets("wisdm_Network")
templateRasterSheet = myScenario.datasheets(
    "wisdm_TemplateRaster", show_full_paths=True)
multiprocessingSheet = myScenario.datasheets("core_Multiprocessing")

spatialMultiprocessingSheet = myScenario.datasheets(
    "core_SpatialMultiprocessing")

# update progress bar
ps.environment.progress_bar()

# Functions -----------------------------------------------------------------
# Resample a raster


def resample_grid(input_grid_filepath, template_grid_filepath, output_filename, output_directory):
    with rasterio.open(input_grid_filepath) as src:
        with rasterio.open(template_grid_filepath) as template:

            # Resample data to target shape
            updated_grid = src.read(
                out_shape=(
                    src.count,
                    template.shape[0],
                    template.shape[1]
                ),
                resampling=enums.Resampling.nearest
            )

            out_meta = src.meta
            out_meta.update({
                "height": updated_grid.shape[1],
                "width": updated_grid.shape[2],
                "transform": template.transform})

            resampled_grid_filepath = os.path.join(
                output_directory, output_filename)

        with rasterio.open(resampled_grid_filepath, "w", **out_meta) as dst:
            dst.write(updated_grid)

        return resampled_grid_filepath

# Error handling ----------------------------------------------------------


# check that template raster was provided and valid
if templateRasterSheet.RasterFilePath.isnull().item():
    raise ValueError(
        "No template raster provided for spatial multiprocessing.")
# check that template is a tif file
if not templateRasterSheet.RasterFilePath.item().lower().endswith(".tif"):
    raise ValueError("Template raster must be a .tif file.")

# stop if tile count is set to 1
if templateRasterSheet.TileCount.item() == 1:
    # inform user and skip rest of script
    ps.environment.update_run_log(
        "The tile count for spatial multiprocessing is set to 1. Spatial multiprocessing grid not created.")
    ps.environment.progress_bar(report_type="end")
else:
    # Load data ---------------------------------------------------------------

    # Load template raster
    templatePath = templateRasterSheet.RasterFilePath.item()
    templateRaster = rioxarray.open_rasterio(templatePath)

    # Set defaults -------------------------------------------------------------

    # Set desired number of cells per tile
    # # Note: max cell count is 100 million cells per tile
    # This is approximate, the actual number of cells per tile may vary
    skipProcessing = False
    if templateRasterSheet.TileCount.isnull().item():
        if templateRaster.size <= 1e5:
            ps.environment.update_run_log("Template raster has only " + str(
                templateRaster.size) + " pixels. Spatial multiprocessing is not required.")
            ps.environment.progress_bar(report_type="end")
            skipProcessing = True
        else:
            # Compute number of tiles based on max tile size
            for maxTileSize in [1e5, 5e5, 1e6, 2.5e6, 5e6, 7.5e6, 1e7, 2e7, 3e7, 4e7, 5e7, 6e7, 7e7, 8e7, 9e7, 1e8, 5e8]:
                if templateRaster.size >= maxTileSize:
                    numTiles = math.ceil(templateRaster.size / maxTileSize)
                    if maxTileSize <= 5e6 and numTiles < 10:
                        break
                    if maxTileSize >= 5e6 and numTiles < 60:
                        break

            if numTiles < 2:
                numTiles = 2
    else:
        numTiles = templateRasterSheet.TileCount.item()

    if not skipProcessing:
        tile_size = templateRaster.size/numTiles  # 7500000

        # update progress bar
        ps.environment.progress_bar()

        # Build tiling grid --------------------------------------------------------

        # Calculate ncol and nrow for one tile
        tile_dimension = math.ceil(
            math.sqrt((templateRaster.size / tile_size)))

        # Define coords of smp grid ----
        x_left = templateRaster.coords['x'][0].item()
        x_right = templateRaster.coords['x'][-1].item()
        x_coords = np.linspace(x_left, x_right, num=tile_dimension)

        y_left = templateRaster.coords['y'][0].item()
        y_right = templateRaster.coords['y'][-1].item()
        y_coords = np.linspace(y_left, y_right, num=tile_dimension)

        coords = {'band': [1], 'x': x_coords, 'y': y_coords}

        # Create smp grid ----------------------------------------------------------
        small_grid = xarray.DataArray(
            data=np.linspace(1, pow(tile_dimension, 2), num=pow(
                tile_dimension, 2)).reshape(1, tile_dimension, tile_dimension),
            dims=templateRaster.dims,
            coords=coords,
            attrs=templateRaster.attrs)

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

        # Mask to analysis area ----------------------------------------------------

        mask = templateRaster
        nodata_value = templateRaster.rio.nodata
        if np.isnan(nodata_value):
            mask = ~np.isnan(mask)*1
        else:
            mask = (mask != nodata_value)*1

        smp_grid = rioxarray.open_rasterio(
            os.path.join(ssimTempDir, "smp-grid.tif"))

        smp_grid = smp_grid.astype(np.int16) * mask

        # Set output filename of smp grid
        oldVals = np.unique(smp_grid)
        num_tiles = len(oldVals)-1
        output_filename = "smpGrid-" + \
            str(num_tiles) + "-" + str(int(tile_size/1e3)) + "k.tif"

        # update template raster datasheet
        templateRasterSheet.TileCount = [num_tiles]
        myScenario.save_datasheet(
            name="wisdm_TemplateRaster", data=templateRasterSheet)

        # update progress bar
        ps.environment.progress_bar()

        # Reclassify tiles --------------------------------------------------------

        newVals = list(range(0, num_tiles+1))

        lookups = list(zip(oldVals, newVals))
        idx, val = np.asarray(lookups).T
        lookup_array = np.zeros(idx.max() + 1)
        lookup_array[idx] = val

        smp_grid = lookup_array[smp_grid]

        # Update 0 to NA Value
        smp_grid[smp_grid == 0] = -9999

        # Convert to int
        smp_grid = smp_grid.astype(np.int16)

        # Save sampling grid -------------------------------------------------------

        with rasterio.open(os.path.join(ssimTempDir, "smp-grid.tif")) as src:
            out_meta = src.meta
        with rasterio.open(os.path.join(ssimTempDir, output_filename), "w", compress='deflate', **out_meta) as dst:
            dst.write(smp_grid)

        # Save to sample grid to scenario
        # add a row to the datasheet
        spatialMultiprocessingSheet["MaskFileName"] = [
            os.path.join(ssimTempDir, output_filename)]
        myScenario.save_datasheet(
            name="core_SpatialMultiprocessing", data=spatialMultiprocessingSheet)

        # update progress bar
        ps.environment.progress_bar(report_type="end")
