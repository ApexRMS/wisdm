## ---------------------------------
## wisdm - ensemble model
## ApexRMS, March 2024
## ---------------------------------

# built under Python version 3.11.0 & SyncroSim version 3.0.0
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
import rasterio
import pysyncrosim as ps 
import numpy as np
import rioxarray
import xarray as xr
import dask
from dask.distributed import Client, Lock
import pandas as pd   
import pyproj

ps.environment.update_run_log('10 - Ensemble Model => Begin')

## Modify the os PROJ path (when running with Conda) ----
myLibrary = ps.Library()
mySession = ps.Session()

result = mySession._Session__call_console(["--conda", "--config"])
conda_fpath = result.stdout.decode('utf-8').strip().split("Conda path is currently: ")[1]
conda_fpath = os.path.normpath(conda_fpath)
if myLibrary.datasheets("core_Option").UseConda.item() == "Yes":
    library_folder = os.path.join(conda_fpath, "envs", "wisdm-2", "wisdm-conda-s3", "Library")
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

#%% Connect to SyncroSim library ------------------------------------------------

# Load current scenario
myScenario = ps.Scenario()  

# Create a temporary folder for storing rasters
ssimTempDir = ps.runtime_temp_folder("DataTransfer\\Scenario-" + str(myScenario.sid))
# ssimTempDir = myLibrary.info["Value"][myLibrary.info.Property == "Temporary files:"].item()

# Get path to scnario inputs 
ssimInputDir = myScenario.library.location + ".input\\Scenario-" + str(myScenario.sid)

# Load datasheets
# inputs
networkSheet = myScenario.library.datasheets("wisdm_Network")
spatialOutputSheet = myScenario.datasheets("wisdm_OutputSpatial", show_full_paths=True)
ensembleOptionsSheet = myScenario.datasheets("wisdm_EnsembleOptions")
multiprocessingSheet = myScenario.datasheets("core_Multiprocessing")

# outputs
ensembleOutputSheet = myScenario.datasheets("wisdm_OutputEnsemble", empty = True)

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

# Set ensemble defaults if not provided
if len(ensembleOptionsSheet) == 0:
    ensembleOptionsSheet = pd.DataFrame({'MakeProbabilityEnsemble': "Yes", 
                                         'ProbabilityMethod': "Mean", 
                                         "NormalizeProbability": "No", 
                                         'MakeBinaryEnsemble': "No", 
                                         'IgnoreNA': "Yes"})
if ensembleOptionsSheet.MakeProbabilityEnsemble.item() == "Yes":
    if pd.isnull(ensembleOptionsSheet.ProbabilityMethod.item()):
        ensembleOptionsSheet['ProbabilityMethod'] = "Mean"
    if pd.isnull(ensembleOptionsSheet.NormalizeProbability.item()):
        ensembleOptionsSheet['NormalizeProbability'] = "No"
if ensembleOptionsSheet.MakeBinaryEnsemble.item() == "Yes":
    if pd.isnull(ensembleOptionsSheet.BinaryMethod.item()):
        ensembleOptionsSheet['BinaryMethod'] = "Mean"
if pd.isnull(ensembleOptionsSheet.IgnoreNA.item()):
    ensembleOptionsSheet['IgnoreNA'] = "Yes"

myScenario.save_datasheet(name="wisdm_EnsembleOptions", data=ensembleOptionsSheet)

# update progress bar
ps.environment.progress_bar()

#%% Load maps ---------------------------------------------------------------------------
# set defualt chunk dimensions
chunkDims = 4096 #1024

# if normalze probability is yes, normalize probability rasters
if ensembleOptionsSheet.NormalizeProbability.item() == "Yes":
    # normalize function
    def norm(dx, min, max):
        dn = dx.to_numpy()
        dn = np.where(dn == -9999, np.nan, dn)
        dn = dn/100
        dnOut = (dn - min) / (max - min)
        dnOut = dnOut * 100
        dnOut = np.where(np.isnan(dnOut), -9999, dnOut)
        dxOut = xr.DataArray(dnOut, dims=dx[range(1)].dims, coords=dx[range(1)].coords)
        return dxOut

    inputFiles = spatialOutputSheet.ProbabilityRaster.tolist()
    for i in range(len(inputFiles)):
        # get min and max value of raster
        with rioxarray.open_rasterio(inputFiles[i]) as src:
            dn = src.to_numpy()
            dn = np.where(dn == -9999, np.nan, dn)
            minVal = np.nanmin(dn)/100
            maxVal = np.nanmax(dn)/100
            
        # normalize raster
        r = rioxarray.open_rasterio(inputFiles[i], chunks = {'y': chunkDims, 'x': chunkDims}, lock=False)
        rOut = r.map_blocks(norm, args=[minVal, maxVal], template=r)
        
        # save normalized raster to temp folder
        fname = "norm_map_" + str(i) + ".tif"
        rOut.rio.to_raster(os.path.join(ssimTempDir, fname), tiled=True, lock=Lock("rio", client=client), windowed = True, overwrite= True, compress = 'lzw')

#%% Set ensemble functions 
if ensembleOptionsSheet.IgnoreNA.item() == "Yes":
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
        dnNan = np.sum(dn, axis=axis)
        dnNan = dnNan == len(dn)*-9999
        dnNan = dnNan.astype(int)
        dnNan = np.where(dnNan == 1, np.nan, dnNan)
        
        dn = np.where(dn == -9999, np.nan, dn)
        dnOut = np.nansum(dn, axis=axis) # nansum will treat nan values as zero - output will have zeros where all inputs are nan
        
        dnOut= np.sum([dnNan,dnOut], axis=0)
        dnOut = np.where(np.isnan(dnOut), -9999, dnOut)
        dnOut = np.expand_dims(dnOut, axis=0)
        dxOut = xr.DataArray(dnOut, dims=dx[range(1)].dims, coords=dx[range(1)].coords)
        return dxOut
else:
    # mean function
    def mean(dx, axis=0):
        dn = dx.to_numpy()
        dn = np.where(dn == -9999, np.nan, dn)
        dnOut = np.mean(dn, axis=axis)
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
        dnOut = np.expand_dims(dnOut, axis=0)
        dxOut = xr.DataArray(dnOut, dims=dx[range(1)].dims, coords=dx[range(1)].coords)
        return dxOut

#%% if ensemble probability is yes, load probability rasters
if ensembleOptionsSheet.MakeProbabilityEnsemble.item() == "Yes":
    
    if ensembleOptionsSheet.NormalizeProbability.item() == "Yes":
        inputFiles = [os.path.join(ssimTempDir, f) for f in os.listdir(ssimTempDir) if f.startswith("norm_map_")]
    else:
        inputFiles = spatialOutputSheet.ProbabilityRaster.tolist()
    
    inputStack = xr.concat([rioxarray.open_rasterio(f, chunks = {'y': chunkDims, 'x': chunkDims}, lock=False) for f in inputFiles], dim = "band").chunk({'band': -1, 'y': chunkDims, 'x': chunkDims})
        
    if ensembleOptionsSheet.ProbabilityMethod.item() == "Mean":
        
        # Calculate mean values block-by-block
        outputStack = inputStack.map_blocks(mean, template=inputStack[range(1)])
        
        # Set nodata flag
        outputStack.rio.write_nodata(-9999, inplace = True)
        
        outputStack.rio.to_raster(os.path.join(ssimTempDir, 'prob_mean.tif'), tiled=True, lock=Lock("rio", client=client), windowed = True, overwrite= True, compress = 'lzw')

        if len(ensembleOutputSheet) == 0:
            newRow = pd.DataFrame({'ProbabilityRasterMean': [os.path.join(ssimTempDir, 'prob_mean.tif')]})
            ensembleOutputSheet = pd.concat([ensembleOutputSheet, newRow], ignore_index=True)
        else:
            ensembleOutputSheet['ProbabilityRasterMean'] = os.path.join(ssimTempDir, 'prob_mean.tif')
    
    if ensembleOptionsSheet.ProbabilityMethod.item() == "Sum":
        
        # Calculate sum values block-by-block
        outputStack = inputStack.map_blocks(sum, template=inputStack[range(1)])
        
        # Set nodata flag
        outputStack.rio.write_nodata(-9999, inplace = True)
        
        outputStack.rio.to_raster(os.path.join(ssimTempDir, 'prob_sum.tif'), tiled=True, lock=Lock("rio", client=client), windowed = True, overwrite= True, compress = 'lzw')

        if len(ensembleOutputSheet) == 0:
            newRow = pd.DataFrame({'ProbabilityRasterSum': [os.path.join(ssimTempDir, 'prob_sum.tif')]})
            ensembleOutputSheet = pd.concat([ensembleOutputSheet, newRow], ignore_index=True)
        else:
            ensembleOutputSheet['ProbabilityRasterSum'] = os.path.join(ssimTempDir, 'prob_sum.tif')

# update progress bar
ps.environment.progress_bar()

#%% if ensemble binary is yes, load binary rasters
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
                newRow = pd.DataFrame({'BinaryRasterMean': [os.path.join(ssimTempDir, 'bin_mean.tif')]})
                ensembleOutputSheet = pd.concat([ensembleOutputSheet, newRow], ignore_index=True)
            else:
                ensembleOutputSheet['BinaryRasterMean'] = os.path.join(ssimTempDir, 'bin_mean.tif')
        
        if ensembleOptionsSheet.BinaryMethod.item() == "Sum":
            
            # Calculate sum values block-by-block
            outputStack = inputStack.map_blocks(sum, template=inputStack[range(1)])
            
            # Set nodata flag
            outputStack.rio.write_nodata(-9999, inplace = True)
            
            outputStack.rio.to_raster(os.path.join(ssimTempDir, 'bin_sum.tif'), tiled=True, lock=Lock("rio", client=client), windowed = True, overwrite= True, compress = 'lzw')
    
            if len(ensembleOutputSheet) == 0:
                newRow = pd.DataFrame({'BinaryRasterSum': [os.path.join(ssimTempDir, 'bin_sum.tif')]})
                ensembleOutputSheet = pd.concat([ensembleOutputSheet, newRow], ignore_index=True)
            else:
                ensembleOutputSheet['BinaryRasterSum'] = os.path.join(ssimTempDir, 'bin_sum.tif')

# update progress bar
ps.environment.progress_bar()

#%% Save ensemble maps -------------------------------------------------------------------
    
# Save data to scenario 
myScenario.save_datasheet(name="wisdm_OutputEnsemble", data=ensembleOutputSheet)   

# update progress bar
ps.environment.progress_bar(report_type="end")
# %%
