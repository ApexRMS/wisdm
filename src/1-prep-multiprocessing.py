## ---------------------------------
## wisdm - prep multiprocessing
## ApexRMS, March 2024
## ---------------------------------

# built under Python version 3.11.0 & SyncroSim version 3.0.0
# This script generates a spatial multiprocessing (smp) grid from a 
# template raster -- here specificaly set up to generaet smp grids 
# for conus map zones

## Source dependencies --------------------------------------------------------

# Set up environment 
import os
import sys
import platform

# Ensure conda environment packages take priority over system Python packages
# Detect conda prefix from CONDA_PREFIX env var OR derive from sys.executable
conda_prefix = os.environ.get("CONDA_PREFIX")
print(f"[DEBUG] Initial CONDA_PREFIX: {conda_prefix}", file=sys.stderr)
print(f"[DEBUG] sys.executable: {sys.executable}", file=sys.stderr)

if not conda_prefix:
    # CONDA_PREFIX not set (happens when Python is called directly without activation)
    # Derive it from sys.executable path
    python_path = os.path.dirname(sys.executable)
    print(f"[DEBUG] Derived path from executable: {python_path}", file=sys.stderr)
    if "envs" in python_path or "conda" in python_path.lower():
        conda_prefix = python_path
        print(f"[DEBUG] Set conda_prefix to: {conda_prefix}", file=sys.stderr)

if conda_prefix and os.path.exists(conda_prefix):
    print(f"[DEBUG] conda_prefix exists: {conda_prefix}", file=sys.stderr)
    # CRITICAL: Set up DLL search paths for Windows BEFORE importing any packages
    if platform.system() == "Windows":
        # Add conda environment's Library\bin to PATH for DLL loading
        library_bin = os.path.join(conda_prefix, "Library", "bin")
        print(f"[DEBUG] library_bin path: {library_bin}", file=sys.stderr)
        if os.path.exists(library_bin):
            print(f"[DEBUG] library_bin exists! Adding to PATH", file=sys.stderr)
            # Prepend to PATH so conda env DLLs take priority
            os.environ["PATH"] = library_bin + os.pathsep + os.environ.get("PATH", "")

            # Set CONDA_PREFIX if not already set (some packages expect it)
            if "CONDA_PREFIX" not in os.environ:
                os.environ["CONDA_PREFIX"] = conda_prefix
                print(f"[DEBUG] Set CONDA_PREFIX env var", file=sys.stderr)

            # Also use os.add_dll_directory for Python 3.8+ (more reliable on Windows)
            if hasattr(os, "add_dll_directory"):
                os.add_dll_directory(library_bin)
                print(f"[DEBUG] Called os.add_dll_directory", file=sys.stderr)
        else:
            print(f"[DEBUG] ERROR: library_bin does NOT exist!", file=sys.stderr)

    # Remove user site-packages from sys.path to prevent conflicts
    sys.path = [p for p in sys.path if not ("AppData\\Roaming\\Python" in p or "AppData/Roaming/Python" in p)]

    # Try Windows path structure (Lib)
    conda_site_packages = os.path.join(conda_prefix, "Lib", "site-packages")

    # If Windows path doesn't exist, try Unix structure (lib)
    if not os.path.exists(conda_site_packages):
        python_version = f"python{sys.version_info.major}.{sys.version_info.minor}"
        conda_site_packages = os.path.join(conda_prefix, "lib", python_version, "site-packages")

    # Only modify sys.path if we found a valid site-packages directory
    if os.path.exists(conda_site_packages) and conda_site_packages not in sys.path:
        sys.path.insert(0, conda_site_packages)

# import modules
import rasterio
import pysyncrosim as ps
from rasterio import enums
import rioxarray
import xarray
import numpy as np
import math
import pyproj

import dask
from dask.distributed import Client

# Set progress bar ---------------------------------------------------------

steps = 6
ps.environment.update_run_log('1 - Prepare Multiprocessing => Begin')
ps.environment.progress_bar(report_type = "begin", total_steps = steps)

# Connect to SyncroSim library --------------------------------------------

# Load current scenario
myScenario = ps.Scenario()  
myLibrary = myScenario.library
mySession = ps.Session()

# Set GDAL/PROJ paths if using conda environment ---------------------------

if myLibrary.datasheets("core_Option").UseConda.item() == "Yes":
    conda_env_path = os.environ.get("CONDA_PREFIX") or sys.prefix

    # Platform-specific paths: Windows uses Library subdirectory, Unix doesn't
    if platform.system() == "Windows":
        library_folder = os.path.join(conda_env_path, "Library")
    else:
        library_folder = conda_env_path

    gdal_folder = os.path.join(library_folder, "share", "gdal")
    proj_folder = os.path.join(library_folder, "share", "proj")
    certifi_folder = os.path.join(library_folder, "ssl", "cacert.pem")
    ps.environment.update_run_log("GDAL path: " + gdal_folder)
    ps.environment.update_run_log("PROJ path: " + proj_folder)
    os.environ['GDAL_DATA'] = gdal_folder
    os.environ['GDAL_CURL_CA_BUNDLE'] = certifi_folder
    os.environ["PROJ_DATA"] = proj_folder
    os.environ['PROJ_CURL_CA_BUNDLE'] = certifi_folder
    os.environ["PROJ_LIB"] = proj_folder
    pyproj.datadir.set_data_dir(proj_folder)
    pyproj.network.set_ca_bundle_path(certifi_folder)
    ps.environment.update_run_log("pyproj data directory: " + pyproj.datadir.get_data_dir())    

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

# Set up dask client ------------------------------------------------------
if multiprocessingSheet.EnableMultiprocessing.item() == "Yes":
    num_threads = multiprocessingSheet.MaximumJobs.item()
else:
    num_threads = 1

# Note: Follow link in output to view progress
dask.config.set(**{'temporary-directory': os.path.join(ssimTempDir, 'dask-worker-space')})
client = Client(threads_per_worker = num_threads, n_workers = 1, processes=False)
# client

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

# Error handling ----------------------------------------------------------

# check that template raster was provided and valid
if templateRasterSheet.RasterFilePath.isnull().item():
    raise ValueError("No template raster provided for spatial multiprocessing.")
# check that template is a tif file
if not templateRasterSheet.RasterFilePath.item().lower().endswith(".tif"):
    raise ValueError("Template raster must be a .tif file.")

# stop if tile count is set to 1
if templateRasterSheet.TileCount.item() == 1:
    # inform user and skip rest of script
    ps.environment.update_run_log("The tile count for spatial multiprocessing is set to 1. Spatial multiprocessing grid not created.")
    ps.environment.progress_bar(report_type = "end")
else:
    # Load data ---------------------------------------------------------------

    # Load template raster
    templatePath = templateRasterSheet.RasterFilePath.item()
    templateRaster = rioxarray.open_rasterio(templatePath)

    # Set defaults -------------------------------------------------------------

    # Set desired number of cells per tile 
    # # Note: max cell count is 100 million cells per tile
    # This is approximate, the actual number of cells per tile may vary
    if templateRasterSheet.TileCount.isnull().item():
        if templateRaster.size <= 1e5:
            ps.environment.update_run_log("Template raster has only " + str(templateRaster.size) + " pixels. Spatial multiprocessing is not required.")
            ps.environment.progress_bar(report_type = "end")
            sys.exit(0)
        elif templateRaster.size > 1e5:
            # Compute number of tiles based on max tile size
            for maxTileSize in [1e5, 5e5, 1e6, 2.5e6, 5e6, 7.5e6, 1e7, 2e7, 3e7, 4e7, 5e7, 6e7, 7e7, 8e7, 9e7, 1e8, 5e8]:
                if templateRaster.size >= maxTileSize:
                    numTiles = math.ceil(templateRaster.size/ maxTileSize)
                    if maxTileSize <= 5e6 and numTiles < 10:
                        break
                    if maxTileSize >= 5e6 and numTiles < 60:
                        break

            if numTiles < 2:
                numTiles = 2
    else:
        numTiles = templateRasterSheet.TileCount.item()   
        
    tile_size = templateRaster.size/numTiles # 7500000

    # Create contiguous tiles
    # contig = True

    # update progress bar
    ps.environment.progress_bar()

    # Build tiling grid --------------------------------------------------------

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

    # Create smp grid ----------------------------------------------------------
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

    # Mask to analysis area ----------------------------------------------------

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

    # Reclassify tiles --------------------------------------------------------

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

    # Save sampling grid -------------------------------------------------------

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

