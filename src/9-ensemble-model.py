## ---------------------------------
## wisdm - ensemble model
## ApexRMS, December 2023
## ---------------------------------

# built under Python version 3.11.0 & SyncroSim version 2.4.42
# Script pulls in template, probability, and binary rasters and outputs and
# summmarizes the data based on user-defined options (e.g., sum or mean) and 
# outputs the ensemble model results

#%% Source dependencies ----------------------------------------------------------

## Modify os path if multiple GDAL installations ----
import os
import glob
from win32api import GetFileVersionInfo, LOWORD, HIWORD

gdal_installations = []
if "PATH" in os.environ:
  for p in os.environ["PATH"].split(os.pathsep):
    if p and glob.glob(os.path.join(p, "gdal*.dll")):
      gdal_installations.append(os.path.abspath(p))

if len(gdal_installations) > 1:
    for folder in gdal_installations:
        filenames = [f for f in os.listdir(folder) if f.startswith("gdal") & f.endswith(".dll")]

        for filename in filenames:
            filename = os.path.join(folder, filename)
        
            if not os.path.exists(filename):           
                print("no gdal dlls found in " + folder)
                os.environ['PATH'] = os.pathsep.join(
                        [p for p in os.environ['PATH'].split(os.pathsep) if folder not in p])
                continue
            try:
                info = GetFileVersionInfo (filename, "\\")
            except:
                continue
            
            major_version = HIWORD (info['FileVersionMS'])
            minor_version = LOWORD (info['FileVersionMS'])

            if (major_version < 3) | (minor_version < 6):
                os.environ['PATH'] = os.pathsep.join(
                    [p for p in os.environ['PATH'].split(os.pathsep) if folder not in p])

## dependencies
import pysyncrosim as ps 
import numpy as np
import rioxarray
import xarray as xr
import dask
from dask.distributed import Client, Lock
import pandas as pd   

## Modify the os PROJ path (when running with Conda) ----
myLibrary = ps.Library()
mySession = ps.Session()

result = mySession._Session__call_console(["--conda", "--config"])
conda_fpath = result.stdout.decode('utf-8').strip().split(": ")[1]
if myLibrary.datasheets("core_Options").UseConda.item() == "Yes":
    os.environ["PROJ_DATA"] = os.path.join(conda_fpath , "envs\\wisdm\\wisdm-py-conda\\Library\\share\\proj")
    os.environ['PROJ_CURL_CA_BUNDLE'] = os.path.join(conda_fpath , "envs\\wisdm\\wisdm-py-conda\\Library\\ssl\\cacert.pem")

#%% Connect to SyncroSim library ------------------------------------------------

# Load current scenario
myScenario = ps.Scenario()  

# Create a temporary folder for storing rasters
ssimTempDir = ps.runtime_temp_folder("Data")

# Get path to scnario inputs 
ssimInputDir = myScenario.library.location + ".input\Scenario-" + str(myScenario.sid)

# Load datasheets
# inputs
networkSheet = myScenario.library.datasheets("Network")
spatialOutputSheet = myScenario.datasheets("SpatialOutputs", show_full_paths=True)
ensembleOptionsSheet = myScenario.datasheets("EnsembleOptions")
multiprocessingSheet = myScenario.datasheets("core_Multiprocessing")

# outputs
ensembleOutputSheet = myScenario.datasheets("EnsembleOutputs", empty = True)

#%% Set progress bar ---------------------------------------------------------

steps = 4
ps.environment.progress_bar(report_type = "begin", total_steps = steps)

#%% Set up dask client -------------------------------------------------------
if multiprocessingSheet.EnableMultiprocessing.item() == "Yes":
    num_threads = multiprocessingSheet.MaximumJobs.item()
else:
    num_threads = 1

# Note: Follow link in output to view progress
dask.config.set(**{'temporary-directory': os.path.join(ssimTempDir, 'dask-worker-space')})
client = Client(threads_per_worker = num_threads, n_workers = 1, processes=False)
# client

#%% Check inputs and set defaults ---------------------------------------------

# Set enemble defaults if not provided
if len(ensembleOptionsSheet) == 0:
    ensembleOptionsSheet = ensembleOptionsSheet.append({'MakeProbabilityEnsemble': "Yes", 'ProbabilityMethod': "Mean", 'MakeBinaryEnsemble': "No"}, ignore_index=True)
if ensembleOptionsSheet.MakeProbabilityEnsemble.item() == "Yes":
    if pd.isnull(ensembleOptionsSheet.ProbabilityMethod.item()):
        ensembleOptionsSheet['ProbabilityMethod'] = "Mean"
if ensembleOptionsSheet.MakeBinaryEnsemble.item() == "Yes":
    if pd.isnull(ensembleOptionsSheet.BinaryMethod.item()):
        ensembleOptionsSheet['BinaryMethod'] = "Mean"

myScenario.save_datasheet(name="EnsembleOptions", data=ensembleOptionsSheet)

# update progress bar
ps.environment.progress_bar()

#%% Load maps ---------------------------------------------------------------------------

# set defualt chunk dimensions
chunkDims = 4096 #1024

# mean function
def mean(dx, axis=0):
    dn = dx.to_numpy()
    dn = np.where(dn == -9999, np.nan, dn)
    dnOut = np.nanmean(dn, axis=axis)
    dnOut = np.where(np.isnan(dnOut), -9999, dnOut)
    dnOut = np.expand_dims(dnOut, axis=0)
    dxOut = xr.DataArray(dnOut, dims=dx[range(1)].dims, coords=dx[range(1)].coords)
    return dxOut

# sum function
def sum(dx, axis=0):
    dn = dx.to_numpy()
    dn = np.where(dn == -9999, np.nan, dn)
    dnOut = np.sum(dn, axis=axis)
    dnOut = np.where(np.isnan(dnOut), -9999, dnOut)
    # dnOut = np.nansum(dn, axis=axis) # nansum will treat nan values as zero - output will have zeros where all inputs are nan
    # dnOut = np.where(dnOut == 0, -9999, dnOut) 
    dnOut = np.expand_dims(dnOut, axis=0)
    dxOut = xr.DataArray(dnOut, dims=dx[range(1)].dims, coords=dx[range(1)].coords)
    return dxOut

#%% if ensemble probability is yes, load probability raster
if ensembleOptionsSheet.MakeProbabilityEnsemble.item() == "Yes":
    
    inputFiles = spatialOutputSheet.ProbabilityRaster.tolist()
    inputStack = xr.concat([rioxarray.open_rasterio(f, chunks = {'y': chunkDims, 'x': chunkDims}, lock=False) for f in inputFiles], dim = "band").chunk({'band': -1, 'y': chunkDims, 'x': chunkDims})
        
    if ensembleOptionsSheet.ProbabilityMethod.item() == "Mean":
        
        # Calculate mean values block-by-block
        outputStack = inputStack.map_blocks(mean, template=inputStack[range(1)])
        
        # Set nodata flag
        outputStack.rio.write_nodata(-9999, inplace = True)
        
        outputStack.rio.to_raster(os.path.join(ssimTempDir, 'prob_mean.tif'), tiled=True, lock=Lock("rio", client=client), windowed = True, overwrite= True, compress = 'lzw')

        if len(ensembleOutputSheet) == 0:
            ensembleOutputSheet = ensembleOutputSheet.append({'ProbabilityRasterMean': os.path.join(ssimTempDir, 'prob_mean.tif')}, ignore_index=True)
        else:
            ensembleOutputSheet['ProbabilityRasterMean'] = os.path.join(ssimTempDir, 'prob_mean.tif')
    
    if ensembleOptionsSheet.ProbabilityMethod.item() == "Sum":
        
        # Calculate sum values block-by-block
        outputStack = inputStack.map_blocks(sum, template=inputStack[range(1)])
        
        # Set nodata flag
        outputStack.rio.write_nodata(-9999, inplace = True)
        
        outputStack.rio.to_raster(os.path.join(ssimTempDir, 'prob_sum.tif'), tiled=True, lock=Lock("rio", client=client), windowed = True, overwrite= True, compress = 'lzw')

        if len(ensembleOutputSheet) == 0:
            ensembleOutputSheet = ensembleOutputSheet.append({'ProbabilityRasterSum': os.path.join(ssimTempDir, 'prob_sum.tif')}, ignore_index=True)
        else:
            ensembleOutputSheet['ProbabilityRasterSum'] = os.path.join(ssimTempDir, 'prob_sum.tif')

# update progress bar
ps.environment.progress_bar()

#%% if ensemble binary is yes, load binary raster
if ensembleOptionsSheet.MakeBinaryEnsemble.item() == "Yes":
        
        inputFiles = spatialOutputSheet.BinaryRaster.tolist()
        inputStack = xr.concat([rioxarray.open_rasterio(f, chunks = {'y': chunkDims, 'x': chunkDims}, lock=False) for f in inputFiles], dim = "band").chunk({'band': -1, 'y': chunkDims, 'x': chunkDims})
            
        if ensembleOptionsSheet.BinaryMethod.item() == "Mean":
            
            # Calculate mean values block-by-block
            outputStack = inputStack.map_blocks(mean, template=inputStack[range(1)])
            
            # Set nodata flag
            outputStack.rio.write_nodata(-9999, inplace = True)
            
            outputStack.rio.to_raster(os.path.join(ssimTempDir, 'bin_mean.tif'), tiled=True, lock=Lock("rio", client=client), windowed = True, overwrite= True, compress = 'lzw')
    
            if len(ensembleOutputSheet) == 0:
                ensembleOutputSheet = ensembleOutputSheet.append({'BinaryRasterMean': os.path.join(ssimTempDir, 'bin_mean.tif')}, ignore_index=True)
            else:
                ensembleOutputSheet['BinaryRasterMean'] = os.path.join(ssimTempDir, 'bin_mean.tif')
        
        if ensembleOptionsSheet.BinaryMethod.item() == "Sum":
            
            # Calculate sum values block-by-block
            outputStack = inputStack.map_blocks(sum, template=inputStack[range(1)])
            
            # Set nodata flag
            outputStack.rio.write_nodata(-9999, inplace = True)
            
            outputStack.rio.to_raster(os.path.join(ssimTempDir, 'bin_sum.tif'), tiled=True, lock=Lock("rio", client=client), windowed = True, overwrite= True, compress = 'lzw')
    
            if len(ensembleOutputSheet) == 0:
                ensembleOutputSheet = ensembleOutputSheet.append({'BinaryRasterSum': os.path.join(ssimTempDir, 'bin_sum.tif')}, ignore_index=True)
            else:
                ensembleOutputSheet['BinaryRasterSum'] = os.path.join(ssimTempDir, 'bin_sum.tif')

# update progress bar
ps.environment.progress_bar()

#%% Save ensemble maps -------------------------------------------------------------------
    
# Save data to scenario 
myScenario.save_datasheet(name="EnsembleOutputs", data=ensembleOutputSheet)   

# update progress bar
ps.environment.progress_bar(report_type="end")
# %%
