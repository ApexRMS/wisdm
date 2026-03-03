# ---------------------------------
# wisdm - ensemble model
# ApexRMS, March 2024
# ---------------------------------

# built under Python version 3.11.0 & SyncroSim version 3.0.0
# Script pulls in template, probability, and binary rasters and outputs and
# summmarizes the data based on user-defined options (e.g., sum or mean) and
# outputs the ensemble model results

# %% Source dependencies ----------------------------------------------------------

# IMPORTANT: setup_functions must be imported before any non-stdlib packages.
# It removes user site-packages on import and setupCondaEnv() configures
# conda DLL paths. Linters must not reorder these imports.  # noqa: E402
import os  # noqa: E402
import sys  # noqa: E402
import time  # noqa: E402
import warnings  # noqa: E402
from setup_functions import setupCondaEnv, checkGdalVersion, setupGdalProj  # noqa: E402

setupCondaEnv()
checkGdalVersion()

# dependencies
import rasterio  # noqa: E402
import pysyncrosim as ps  # noqa: E402
import numpy as np  # noqa: E402
import rioxarray  # noqa: E402
import xarray as xr  # noqa: E402
import dask  # noqa: E402
from dask.distributed import Client, LocalCluster  # noqa: E402
import pandas as pd  # noqa: E402
import pyproj  # noqa: E402


def run():

    ps.environment.update_run_log('10 - Ensemble Model => Begin')

    # Modify the os PROJ path (when running with Conda) ----
    myLibrary = ps.Library()
    worker_env = setupGdalProj(myLibrary)

    # %% Connect to SyncroSim library ------------------------------------------------

    # Load current scenario
    myScenario = ps.Scenario()

    # Create a temporary folder for storing rasters
    ssimTempDir = ps.runtime_temp_folder(
        "DataTransfer\\Scenario-" + str(myScenario.sid))
    # ssimTempDir = myLibrary.info["Value"][myLibrary.info.Property == "Temporary files:"].item()

    # Get path to scnario inputs
    ssimInputDir = myScenario.library.location + \
        ".input\\Scenario-" + str(myScenario.sid)

    # Load datasheets
    # inputs
    networkSheet = myScenario.library.datasheets("wisdm_Network")
    spatialOutputSheet = myScenario.datasheets(
        "wisdm_OutputSpatial", show_full_paths=True)
    ensembleOptionsSheet = myScenario.datasheets("wisdm_EnsembleOptions")
    multiprocessingSheet = myScenario.datasheets("core_Multiprocessing")

    # outputs
    ensembleOutputSheet = myScenario.datasheets(
        "wisdm_OutputEnsemble", empty=True)

    # %% Set progress bar ---------------------------------------------------------

    steps = 4
    ps.environment.progress_bar(report_type="begin", total_steps=steps)

    # %% Set up dask client -------------------------------------------------------
    if multiprocessingSheet.EnableMultiprocessing.item() == "Yes":
        num_threads = multiprocessingSheet.MaximumJobs.item()
    else:
        num_threads = 1

    # Set up Dask temporary directory
    dask_temp = os.path.join(ssimTempDir, 'dask-worker-space')
    os.makedirs(dask_temp, exist_ok=True)

    conf = {
        'temporary-directory': dask_temp,
        'distributed.scheduler.worker-ttl': None,
    }
    if worker_env:
        conf['distributed.nanny.environ'] = worker_env

    dask.config.set(conf)

    ps.environment.update_run_log(
        f'Initializing Dask cluster with {num_threads} workers...')

    cluster = LocalCluster(
        n_workers=num_threads,
        memory_limit='auto',
        processes=True)

    ps.environment.update_run_log('Cluster created, connecting client...')

    # Connect client with timeout
    client = Client(cluster, timeout='60s')

    ps.environment.update_run_log(
        f'Dask client ready: {num_threads} workers')

    # %% Check inputs and set defaults ---------------------------------------------

    # Set ensemble defaults if not provided
    if len(ensembleOptionsSheet) == 0:
        ensembleOptionsSheet = pd.DataFrame({'MakeProbabilityEnsemble': ["Yes"],
                                             'ProbabilityMethod': ["Mean"],
                                             "NormalizeProbability": ["No"],
                                             'MakeBinaryEnsemble': ["No"],
                                             'BinaryMethod': [None],
                                             'IgnoreNA': ["Yes"]})
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

    myScenario.save_datasheet(name="wisdm_EnsembleOptions",
                              data=ensembleOptionsSheet)

    # update progress bar
    ps.environment.progress_bar()

    # %% Load maps ---------------------------------------------------------------------------
    # Use auto chunking to respect native tiling from input rasters
    chunkDims = 'auto'

    ps.environment.update_run_log('Starting raster processing...')

    # if normalze probability is yes, normalize probability rasters
    if ensembleOptionsSheet.NormalizeProbability.item() == "Yes":
        ps.environment.update_run_log('Normalizing probability rasters...')
        # normalize function

        def norm(dx, min, max):
            dn = dx.to_numpy()
            dn = np.where(dn == -9999, np.nan, dn)
            dn = dn/100
            dnOut = (dn - min) / (max - min)
            dnOut = dnOut * 100
            dnOut = np.where(np.isnan(dnOut), -9999, dnOut)
            dxOut = xr.DataArray(
                dnOut, dims=dx[range(1)].dims, coords=dx[range(1)].coords)
            return dxOut

        inputFiles = spatialOutputSheet.ProbabilityRaster.tolist()
        ps.environment.update_run_log(
            f'Normalizing {len(inputFiles)} rasters...')

        for i in range(len(inputFiles)):
            ps.environment.update_run_log(
                f'  Processing raster {i+1}/{len(inputFiles)}...')

            # Get min and max value of raster using chunked computation (memory efficient)
            ps.environment.update_run_log(f'  Computing min/max statistics...')
            r = rioxarray.open_rasterio(inputFiles[i], chunks={
                                        'y': chunkDims, 'x': chunkDims}, lock=False)
            r_masked = r.where(r != -9999)  # Mask nodata values
            minVal = float(r_masked.min().compute()) / 100
            maxVal = float(r_masked.max().compute()) / 100
            ps.environment.update_run_log(
                f'    Min: {minVal:.4f}, Max: {maxVal:.4f}')

            # Normalize raster (reuse already-opened raster)
            ps.environment.update_run_log(f'  Normalizing raster {i+1}...')
            rOut = r.map_blocks(norm, args=[minVal, maxVal], template=r)

            # Save normalized raster to temp folder
            fname = "norm_map_" + str(i) + ".tif"
            ps.environment.update_run_log(
                f'  Writing normalized raster {i+1}...')
            start_time = time.time()
            rOut.rio.to_raster(os.path.join(ssimTempDir, fname),
                               tiled=True, overwrite=True, compress='lzw')
            elapsed = (time.time() - start_time) / 60
            ps.environment.update_run_log(
                f'  Completed raster {i+1} in {elapsed:.1f} minutes')

    # %% Set ensemble functions
    if ensembleOptionsSheet.IgnoreNA.item() == "Yes":
        # mean function
        def mean(dx, axis=0):
            dn = dx.to_numpy()
            dn = np.where(dn == -9999, np.nan, dn)
            # Suppress "mean of empty slice" warnings - we handle the NaN result correctly
            with warnings.catch_warnings():
                warnings.filterwarnings(
                    'ignore', message='Mean of empty slice')
                dnOut = np.nanmean(dn, axis=axis)
            dnOut = np.where(np.isnan(dnOut), -9999, dnOut)
            dnOut = np.expand_dims(dnOut, axis=0)
            dxOut = xr.DataArray(
                dnOut, dims=dx[range(1)].dims, coords=dx[range(1)].coords)
            return dxOut

        # sum function
        def sum(dx, axis=0):
            dn = dx.to_numpy()
            dnNan = np.sum(dn, axis=axis)
            dnNan = dnNan == len(dn)*-9999
            dnNan = dnNan.astype(int)
            dnNan = np.where(dnNan == 1, np.nan, dnNan)

            dn = np.where(dn == -9999, np.nan, dn)
            # nansum will treat nan values as zero - output will have zeros where all inputs are nan
            dnOut = np.nansum(dn, axis=axis)

            dnOut = np.sum([dnNan, dnOut], axis=0)
            dnOut = np.where(np.isnan(dnOut), -9999, dnOut)
            dnOut = np.expand_dims(dnOut, axis=0)
            dxOut = xr.DataArray(
                dnOut, dims=dx[range(1)].dims, coords=dx[range(1)].coords)
            return dxOut
    else:
        # mean function
        def mean(dx, axis=0):
            dn = dx.to_numpy()
            dn = np.where(dn == -9999, np.nan, dn)
            dnOut = np.mean(dn, axis=axis)
            dnOut = np.where(np.isnan(dnOut), -9999, dnOut)
            dnOut = np.expand_dims(dnOut, axis=0)
            dxOut = xr.DataArray(
                dnOut, dims=dx[range(1)].dims, coords=dx[range(1)].coords)
            return dxOut

        # sum function
        def sum(dx, axis=0):
            dn = dx.to_numpy()
            dn = np.where(dn == -9999, np.nan, dn)
            dnOut = np.sum(dn, axis=axis)
            dnOut = np.where(np.isnan(dnOut), -9999, dnOut)
            dnOut = np.expand_dims(dnOut, axis=0)
            dxOut = xr.DataArray(
                dnOut, dims=dx[range(1)].dims, coords=dx[range(1)].coords)
            return dxOut

    # %% if ensemble probability is yes, load probability rasters
    if ensembleOptionsSheet.MakeProbabilityEnsemble.item() == "Yes":
        ps.environment.update_run_log('Creating probability ensemble...')

        if ensembleOptionsSheet.NormalizeProbability.item() == "Yes":
            inputFiles = [os.path.join(ssimTempDir, f) for f in os.listdir(
                ssimTempDir) if f.startswith("norm_map_")]
        else:
            inputFiles = spatialOutputSheet.ProbabilityRaster.tolist()

        ps.environment.update_run_log(
            f'Loading {len(inputFiles)} probability rasters...')
        inputStack = xr.concat([rioxarray.open_rasterio(f, chunks={'y': chunkDims, 'x': chunkDims}, lock=False)
                               for f in inputFiles], dim="band").chunk({'band': -1, 'y': chunkDims, 'x': chunkDims})

        if ensembleOptionsSheet.ProbabilityMethod.item() == "Mean":
            ps.environment.update_run_log('Calculating mean ensemble...')

            # Calculate mean values block-by-block
            outputStack = inputStack.map_blocks(
                mean, template=inputStack[range(1)])

            # Set nodata flag
            outputStack.rio.write_nodata(-9999, inplace=True)

            ps.environment.update_run_log('Writing mean ensemble raster...')
            write_start = time.time()
            outputStack.rio.to_raster(os.path.join(
                ssimTempDir, 'prob_mean.tif'), tiled=True, overwrite=True, compress='lzw')
            ps.environment.update_run_log(
                f'Mean ensemble complete in {(time.time() - write_start)/60:.1f} minutes')

            if len(ensembleOutputSheet) == 0:
                newRow = pd.DataFrame(
                    {'ProbabilityRasterMean': [os.path.join(ssimTempDir, 'prob_mean.tif')]})
                ensembleOutputSheet = pd.concat(
                    [ensembleOutputSheet, newRow], ignore_index=True)
            else:
                ensembleOutputSheet['ProbabilityRasterMean'] = os.path.join(
                    ssimTempDir, 'prob_mean.tif')

        if ensembleOptionsSheet.ProbabilityMethod.item() == "Sum":
            ps.environment.update_run_log('Calculating sum ensemble...')

            # Calculate sum values block-by-block
            outputStack = inputStack.map_blocks(
                sum, template=inputStack[range(1)])

            # Set nodata flag
            outputStack.rio.write_nodata(-9999, inplace=True)

            ps.environment.update_run_log('Writing sum ensemble raster...')
            write_start = time.time()
            outputStack.rio.to_raster(os.path.join(
                ssimTempDir, 'prob_sum.tif'), tiled=True, overwrite=True, compress='lzw')
            ps.environment.update_run_log(
                f'Sum ensemble complete in {(time.time() - write_start)/60:.1f} minutes')

            if len(ensembleOutputSheet) == 0:
                newRow = pd.DataFrame(
                    {'ProbabilityRasterSum': [os.path.join(ssimTempDir, 'prob_sum.tif')]})
                ensembleOutputSheet = pd.concat(
                    [ensembleOutputSheet, newRow], ignore_index=True)
            else:
                ensembleOutputSheet['ProbabilityRasterSum'] = os.path.join(
                    ssimTempDir, 'prob_sum.tif')

    # update progress bar
    ps.environment.progress_bar()

    # %% if ensemble binary is yes, load binary rasters
    if ensembleOptionsSheet.MakeBinaryEnsemble.item() == "Yes":
        ps.environment.update_run_log('Creating binary ensemble...')

        inputFiles = spatialOutputSheet.BinaryRaster.tolist()
        ps.environment.update_run_log(
            f'Loading {len(inputFiles)} binary rasters...')
        inputStack = xr.concat([rioxarray.open_rasterio(f, chunks={'y': chunkDims, 'x': chunkDims}, lock=False)
                               for f in inputFiles], dim="band").chunk({'band': -1, 'y': chunkDims, 'x': chunkDims})

        if ensembleOptionsSheet.BinaryMethod.item() == "Mean":
            ps.environment.update_run_log(
                'Calculating binary mean ensemble...')

            # Calculate mean values block-by-block
            # Cast template to float so the 0-1 probability values aren't truncated to integer
            outputStack = inputStack.map_blocks(
                mean, template=inputStack[range(1)].astype(float))

            # Set nodata flag
            outputStack.rio.write_nodata(-9999, inplace=True)

            ps.environment.update_run_log(
                'Writing binary mean ensemble raster...')
            write_start = time.time()
            outputStack.rio.to_raster(os.path.join(
                ssimTempDir, 'bin_mean.tif'), tiled=True, overwrite=True, compress='lzw')
            ps.environment.update_run_log(
                f'Binary mean ensemble complete in {(time.time() - write_start)/60:.1f} minutes')

            if len(ensembleOutputSheet) == 0:
                newRow = pd.DataFrame(
                    {'BinaryRasterMean': [os.path.join(ssimTempDir, 'bin_mean.tif')]})
                ensembleOutputSheet = pd.concat(
                    [ensembleOutputSheet, newRow], ignore_index=True)
            else:
                ensembleOutputSheet['BinaryRasterMean'] = os.path.join(
                    ssimTempDir, 'bin_mean.tif')

        if ensembleOptionsSheet.BinaryMethod.item() == "Sum":
            ps.environment.update_run_log('Calculating binary sum ensemble...')

            # Calculate sum values block-by-block
            outputStack = inputStack.map_blocks(
                sum, template=inputStack[range(1)])

            # Set nodata flag
            outputStack.rio.write_nodata(-9999, inplace=True)

            ps.environment.update_run_log(
                'Writing binary sum ensemble raster...')
            write_start = time.time()
            outputStack.rio.to_raster(os.path.join(
                ssimTempDir, 'bin_sum.tif'), tiled=True, overwrite=True, compress='lzw')
            ps.environment.update_run_log(
                f'Binary sum ensemble complete in {(time.time() - write_start)/60:.1f} minutes')

            if len(ensembleOutputSheet) == 0:
                newRow = pd.DataFrame(
                    {'BinaryRasterSum': [os.path.join(ssimTempDir, 'bin_sum.tif')]})
                ensembleOutputSheet = pd.concat(
                    [ensembleOutputSheet, newRow], ignore_index=True)
            else:
                ensembleOutputSheet['BinaryRasterSum'] = os.path.join(
                    ssimTempDir, 'bin_sum.tif')

    # update progress bar
    ps.environment.progress_bar()

    # %% Save ensemble maps -------------------------------------------------------------------

    # Save data to scenario
    myScenario.save_datasheet(name="wisdm_OutputEnsemble",
                              data=ensembleOutputSheet)

    # Close Dask resources
    client.close()
    cluster.close()

    # update progress bar
    ps.environment.progress_bar(report_type="end")


if __name__ == "__main__":
    run()
# %%
