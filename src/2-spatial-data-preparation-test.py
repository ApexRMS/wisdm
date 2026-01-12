# ---------------------------------
# wisdm - spatial data preparation
# ApexRMS, March 2024
# ---------------------------------

# built under Python version 3.11.0 & SyncroSim version 3.0.0
# Script pulls in template and covariate rasters and projects, aggragates/resamples, and
# clips (PARC) covariate rasters to match template; processes field data to ensure sites
# are in the template CRS and extent; aggregates or weights sites by spatial distribution;
# extracts site-specific covaraite data and generates datasheet of covariate site data.

# Source dependencies ----------------------------------------------------------

# Modify os path if multiple GDAL installations ----
def prep_spatial_data():

    import os, sys
    import glob
    from win32api import GetFileVersionInfo, LOWORD, HIWORD

    gdal_installations = []
    if "PATH" in os.environ:
        for p in os.environ["PATH"].split(os.pathsep):
            if p and glob.glob(os.path.join(p, "gdal*.dll")):
                gdal_installations.append(os.path.abspath(p))

    if len(gdal_installations) > 0:
        for folder in gdal_installations:
            filenames = [f for f in os.listdir(
                folder) if f.startswith("gdal") & f.endswith(".dll")]

            for filename in filenames:
                filename = os.path.join(folder, filename)

                if not os.path.exists(filename):
                    print("no gdal dlls found in " + folder)
                    os.environ['PATH'] = os.pathsep.join(
                        [p for p in os.environ['PATH'].split(os.pathsep) if folder not in p])
                    continue
                try:
                    info = GetFileVersionInfo(filename, "\\")
                except:
                    continue

                major_version = HIWORD(info['FileVersionMS'])
                minor_version = LOWORD(info['FileVersionMS'])

                if (major_version < 3) | (minor_version < 6):
                    os.environ['PATH'] = os.pathsep.join(
                        [p for p in os.environ['PATH'].split(os.pathsep) if folder not in p])

    # dependencies
    import rasterio
    import pysyncrosim as ps
    import numpy as np
    import pandas as pd
    import rioxarray
    import xarray as xr
    import geopandas as gpd
    import shapely
    import dask
    import pyproj
    from datetime import datetime
    import threading

    from uuid import uuid4
    import time
    from contextlib import suppress
    import certifi

    from dask import delayed
    from dask.distributed import Client
    from dask.utils import SerializableLock
    
    from shapely.geometry import Point  # , shape
    from rasterio.enums import Resampling  # , MergeAlg
    from rasterio.vrt import WarpedVRT

    ps.environment.update_run_log('2 - Spatial Data Preparation => Begin')
    ps.environment.progress_bar(
        "message", message="Preparing inputs for spatial data prep...")

    # Suppress pandas Setting with Copy warning
    pd.options.mode.chained_assignment = None

    # Modify the os PROJ path (when running with Conda) ----
    myLibrary = ps.Library()
    
    worker_env = None
    if myLibrary.datasheets("core_Option").UseConda.item() == "Yes":
        conda_env_path = os.environ.get("CONDA_PREFIX") or sys.prefix
        if os.name == "nt": # Windows
            library_folder = os.path.join(conda_env_path, "Library")
            gdal_folder = os.path.join(library_folder, "share", "gdal")
            proj_folder = os.path.join(library_folder, "share", "proj")
            certifi_folder = os.path.join(library_folder, "ssl", "cacert.pem")
        else:  # macOS / Linux
            gdal_folder = os.path.join(conda_env_path, "share", "gdal")
            proj_folder = os.path.join(conda_env_path, "share", "proj")
            try: 
                certifi_folder = certifi.where()
            except Exception:
                certifi_folder = ""
        ps.environment.update_run_log("GDAL path: " + gdal_folder)
        ps.environment.update_run_log("PROJ path: " + proj_folder)
        os.environ['GDAL_DATA'] = gdal_folder
        os.environ["PROJ_LIB"] = proj_folder
        os.environ["PROJ_DATA"] = proj_folder
        pyproj.datadir.set_data_dir(proj_folder)
        if certifi_folder:
            os.environ['GDAL_CURL_CA_BUNDLE'] = certifi_folder     
            os.environ['PROJ_CURL_CA_BUNDLE'] = certifi_folder
            pyproj.network.set_ca_bundle_path(certifi_folder)
        ps.environment.update_run_log(
            "pyproj data directory: " + pyproj.datadir.get_data_dir())
        
        # What future Dask workers should inherit:
        worker_env = {
            "GDAL_DATA": gdal_folder,
            "PROJ_LIB": proj_folder,
        }
        if certifi_folder:
            worker_env["GDAL_CURL_CA_BUNDLE"] = certifi_folder
            worker_env["PROJ_CURL_CA_BUNDLE"] = certifi_folder

    # Helper functions ------------------------------------------------------------
    
    # Function to choose Dask mode
    def chooseDaskMode(template_width: int,
                       template_height: int,
                       num_covariates: int,
                       max_jobs: int) -> dict:
        """
        Minimal heuristic:
        - If template pixels >= 200M OR covariates >= 20 → multi-process
        - Else → single process with multiple threads

        Returns dict for Dask Client + GDAL threads.
        """
        pixelCount = int(template_width) * int(template_height)

        manyPixels = 200_000_000    # ~200M px
        manyCovs  = 20

        useWorkers = (pixelCount >= manyPixels) or (num_covariates >= manyCovs)

        if useWorkers:
            # several worker processes, 1 thread each
            n_workers = max_jobs
            threads_per_worker = 1
            gdal_threads = "2"
        else:
            # one worker, several threads
            n_workers = 1
            threads_per_worker = max_jobs
            gdal_threads = str(max_jobs)

        return {
            "processes": useWorkers,
            "n_workers": n_workers,
            "threads_per_worker": threads_per_worker,
            "gdal_threads": gdal_threads
        }

    # Masking function
    def mask(block, input_nodata, template_nodata):
        # Convert to numpy for faster masking
        npblock = block.to_numpy()

        # Mask and convert all nodata to common value
        temp = np.where(npblock[[0]] == input_nodata,
                        nodata_value, npblock[[0]])
        masked = np.where(npblock[[1]] == template_nodata, nodata_value, temp)

        # Copy data back into appropriate xarray
        return block[[0]].copy(data=masked)
    
    # Function to write raster with retries
    def writeReplaceWithRetry(src_tmp, dst_path, retries=5, delay=0.25):
        # ensure destination is free
        if os.path.exists(dst_path):
            safeUnlink(dst_path)
        backoff = delay
        for i in range(retries):
            try:
                os.replace(src_tmp, dst_path)
                return
            except PermissionError:
                if i == retries - 1:
                    raise
                time.sleep(backoff)
                backoff *= 1.8

    # Function to generate unique temporary file path
    def uniqueTmp(base_dir: str, prefix: str) -> str:
        # Keep trying until we find a non-existent name (Windows-safe)
        for _ in range(50):
            p = os.path.join(base_dir, f"{prefix}_{uuid4().hex[:8]}.tif")
            if not os.path.exists(p):
                return p
        # last-resort slug
        return os.path.join(base_dir, f"{prefix}_{uuid4().hex}.tif")
    
    # Function to safely delete files with retries
    def safeUnlink(path: str, retries: int = 6, delay: float = 0.25):
        # Best-effort, with backoff; ignore if missing
        for i in range(retries):
            with suppress(FileNotFoundError, PermissionError):
                os.remove(path)
            if not os.path.exists(path):
                break
            time.sleep(delay * (1.8 ** i))
        # also clean common GDAL sidecars
        for sfx in (".aux.xml", ".ovr", ".msk"):
            sp = path + sfx
            for i in range(retries):
                with suppress(FileNotFoundError, PermissionError):
                    os.remove(sp)
                if not os.path.exists(sp):
                    break
                time.sleep(delay * (1.8 ** i))

    def predictor_for_dtype(dtype_str: str):
        # GTiff PREDICTOR=2 (horizontal differencing) for integers, =3 for floats
        if dtype_str.startswith("float"):
            return 3
        if dtype_str.startswith("int") or dtype_str.startswith("uint"):
            return 2
        return None

    # Connect to SyncroSim library ------------------------------------------------

    # Load current scenario
    myScenario = ps.Scenario()

    # Create a temporary folder for storing rasters
    ssimTempDir = ps.runtime_temp_folder(os.path.join(
        "DataTransfer", "Scenario-" + str(myScenario.sid)))

    # Load datasheets
    # inputs
    networkSheet = myScenario.library.datasheets("wisdm_Network")
    covariatesSheet = myScenario.project.datasheets("wisdm_Covariates")
    covariateDataSheet = myScenario.datasheets(
        "wisdm_CovariateData", show_full_paths=True)
    templateRasterSheet = myScenario.datasheets(
        "wisdm_TemplateRaster", show_full_paths=True)
    restrictionRasterSheet = myScenario.datasheets(
        "wisdm_RestrictionRaster", show_full_paths=True)
    multiprocessingSheet = myScenario.datasheets("core_Multiprocessing")

    # outputs
    outputCovariateSheet = myScenario.datasheets(
        "wisdm_CovariateData", empty=True)

    # Set progress bar ---------------------------------------------------------

    steps = 3 + len(covariateDataSheet.CovariatesID)
    ps.environment.progress_bar(report_type="begin", total_steps=steps)

    # Check inputs and set defaults ---------------------------------------------

    # Set PROJ network connection
    if networkSheet.NetworkEnabled.item() == "No":
        pyproj.network.set_network_enabled(active=False)

    # Check that a template raster was provided
    if templateRasterSheet.empty or pd.isnull(templateRasterSheet.RasterFilePath.iloc[0]):
        raise ValueError("Template raster is missing.")

    # if provided, ensure restriction raster data is between 0-1
    def check_raster_range(filepath, min_allowed=0.0, max_allowed=1.0, epsilon=1e-8):
        with rasterio.open(filepath) as src:
            for _, window in src.block_windows(1):  # band 1
                data = src.read(1, window=window, masked=True)
                # Skip if all values are nodata
                if np.ma.count(data) == 0:
                    continue
                # Scalarize masked mins/maxes and compare with tolerance
                block_min = float(np.ma.min(data))
                block_max = float(np.ma.max(data))
                if block_min < (min_allowed - epsilon) or block_max > (max_allowed + epsilon):
                    return False  # Out of range found
        return True  # All values are within range

    restrictionPath = None
    if not restrictionRasterSheet.empty:
        restrictionPath = restrictionRasterSheet.RasterFilePath.iloc[0]
        if pd.notnull(restrictionPath):
            if not check_raster_range(restrictionPath, 0, 1):
                raise ValueError(
                    "Restriction raster values must range from 0 to 1.")

    # Identify categorical variables
    catCovs = covariatesSheet.query(
        'IsCategorical == "Yes"').CovariateName.tolist()

    # Ensure that covariate names do not contain spaces
    if covariateDataSheet.CovariatesID.str.contains(r"\s", na=False).any():
        raise ValueError("Covariate names cannot contain spaces.")

    # Convert Resample and Aggregation method columns to strings
    covariateDataSheet = covariateDataSheet.astype({'ResampleMethod': 'string',
                                                    "AggregationMethod": "string"})

    # Set resample defaults
    for i in range(len(covariateDataSheet.ResampleMethod)):
        if pd.isnull(covariateDataSheet.ResampleMethod[i]):
            if covariateDataSheet.CovariatesID[i] in catCovs:
                covariateDataSheet.loc[i,
                                       "ResampleMethod"] = "Nearest Neighbor"
            else:
                covariateDataSheet.loc[i, "ResampleMethod"] = "Bilinear"

    # Set aggregate defaults
    for i in range(len(covariateDataSheet.AggregationMethod)):
        if pd.isnull(covariateDataSheet.AggregationMethod[i]):
            if covariateDataSheet.CovariatesID[i] in catCovs:
                covariateDataSheet.loc[i, "AggregationMethod"] = "Majority"
            else:
                covariateDataSheet.loc[i, "AggregationMethod"] = "Mean"

    # Define rasterio up/down-sample method
    resampleAggregateMethodName = ["Nearest Neighbor", "Bilinear",
                                   "Cubic", "Cubic Spline", "Lanczos", "Mean", "Min", "Max", "Majority"]
    resampleAggregateMethodCodes = ["nearest", "bilinear", "cubic",
                                    "cubic_spline", "lanczos", "average", "min", "max", "mode"]
    covariateDataSheet["rioResample"] = covariateDataSheet.ResampleMethod.replace(
        resampleAggregateMethodName, resampleAggregateMethodCodes)
    covariateDataSheet["rioAggregate"] = covariateDataSheet.AggregationMethod.replace(
        resampleAggregateMethodName, resampleAggregateMethodCodes)

    # Load template raster ----------------------------------------------------------------

    # set defualt chunk dimensions and no data value
    warpMemoryLimit = 8192  # MB
    chunkDims = 1024  # 1024
    nodata_value = -9999

    templatePath = templateRasterSheet.RasterFilePath.item()
    templateRaster = rioxarray.open_rasterio(
        templatePath, chunks={'x': chunkDims, 'y': chunkDims}, lock=False)

    # Get information about template
    templateCRS = templateRaster.rio.crs
    if templateCRS.is_valid == False:
        raise ValueError(
            "Template has an invalid CRS (authority code). See {documention} for a list of accepted authority codes.")
    templateTransform = templateRaster.rio.transform()
    templatePixelSize = np.abs(templateRaster.rio.resolution()).prod()
    templateBounds = templateRaster.rio.bounds()
    templateHeight = templateRaster.rio.height
    templateWidth = templateRaster.rio.width

    # update progress bar
    ps.environment.progress_bar()

    # Set up Dask client (single, distributed, process workers) -----------------
    
    mode = chooseDaskMode(
        template_width=templateWidth,
        template_height=templateHeight,
        num_covariates=len(covariateDataSheet.CovariatesID),
        max_jobs=int(multiprocessingSheet.MaximumJobs.item()) if multiprocessingSheet.EnableMultiprocessing.item() == "Yes" else 1
    )    

    # keep worker temp in scenario temp dir (atomic on same volume)
    conf = {
        'temporary-directory': os.path.join(ssimTempDir, 'dask-worker-space'),
        'distributed.scheduler.worker-ttl': None,
    }
    if worker_env:
        conf['distributed.worker.environ'] = worker_env
    
    dask.config.set(conf)

    client = Client(processes=mode["processes"], n_workers=int(mode["n_workers"]), threads_per_worker=mode["threads_per_worker"])
    rio_lock = SerializableLock()
    ps.environment.update_run_log(f"Dask client: {mode["n_workers"]} worker(s) | {getattr(client, 'dashboard_link', '')}")

    # Loop through and "PARC" covariate rasters -------------------------------------------

    # First check that all rasters have a valid crs and have valid dimensions
    invalidCRS = []
    for i in range(len(covariateDataSheet.CovariatesID)):
        # Load covariate rasters
        covariatePath = covariateDataSheet.RasterFilePath[i]
        covariateRaster = rioxarray.open_rasterio(covariatePath, chunks=True)
        if pd.isnull(covariateRaster.rio.crs):
            invalidCRS.append(covariateDataSheet.CovariatesID[i])
        elif covariateRaster.rio.crs.is_valid == False:
            invalidCRS.append(covariateDataSheet.CovariatesID[i])
        # elif "unnamed" in str(covariateRaster.rio.crs):
        #     invalidCRS.append(covariateDataSheet.CovariatesID[i])
    if len(invalidCRS) > 0:
        msg = "The following covariate rasters have an invalid or unknown CRS:\n"
        for i in invalidCRS:
            msg += i + "\n"
        msg += "\nEnsure that the covariate rasters have a valid CRS before continuing."
        ps.environment.update_run_log(msg, type="status")
        raise ValueError("Covariate rasters have invalid CRS.")

    # define processing function for each covariate
    def process_one_covariate_worker_fast(
        covariatePath, covId, tmpFinal, 
        templatePath, templateCRS, templateTransform, templateHeight, templateWidth,
        rioResampleMethod,  warpMemoryLimit, 
        signed_dtype, nodata_value,
        blockDims=512, compress="lzw", predictor=None, bigtiff="IF_NEEDED"
    ):
        with rasterio.open(templatePath) as templ:
            profile = templ.profile
            profile.update(dtype=signed_dtype, 
                           nodata=nodata_value,
                           compress=compress, 
                           tiled=True,
                           blockxsize=blockDims, blockysize=blockDims, 
                           bigtiff=bigtiff
            )
            if predictor is not None:
                profile["predictor"] = predictor

            with rasterio.open(tmpFinal, "w", **profile) as dst:
                with rasterio.open(covariatePath) as src:
                    with WarpedVRT(
                        src,
                        resampling=rioResampleMethod,
                        crs=templateCRS,
                        transform=templateTransform,
                        height=templateHeight,
                        width=templateWidth,
                        warp_mem_limit=warpMemoryLimit,
                        warp_extras={'NUM_THREADS': mode['gdal_threads']}
                    ) as vrt:
                        # quick bounds sanity
                        tb = templateBounds
                        cb = vrt.bounds
                        if cb.left > tb[2] or cb.right < tb[0] or cb.bottom > tb[3] or cb.top < tb[1]:
                            raise ValueError(f"The extent of the {covId} raster does not overlap the template raster.")

                        # iterate destination by block
                        for _, w in dst.block_windows(1):
                            src_block  = vrt.read(1, window=w, out_dtype=dst.dtypes[0], masked=False)
                            tmpl_block = templ.read(1, window=w, masked=False)

                            if vrt.nodata is not None:
                                src_block = np.where(src_block == vrt.nodata, nodata_value, src_block)
                            if templ.nodata is not None:
                                src_block = np.where(tmpl_block == templ.nodata, nodata_value, src_block)

                            dst.write(src_block, 1, window=w)
        return True

    def process_one_covariate_worker(
        covariatePath, covId, tmpUnmasked, tmpFinal,
        templatePath, templateBounds, templateCRS, templateTransform, templateHeight, templateWidth,
        rioResampleMethod, chunkDims, warpMemoryLimit,
        signed_dtype, nodata_value
    ):

        # 1) Reproject/clip to uncompressed temp via VRT
        with rasterio.open(covariatePath) as src:
            with WarpedVRT(
                src,
                resampling=rioResampleMethod,
                crs=templateCRS,
                transform=templateTransform,
                height=templateHeight,
                width=templateWidth,
                warp_mem_limit=warpMemoryLimit,
                warp_extras={'NUM_THREADS': mode['gdal_threads']},
            ) as vrt:
                with rioxarray.open_rasterio(vrt, chunks={'x': chunkDims, 'y': chunkDims}, lock=False) as covariateRaster:
                    covariateRaster.rio.to_raster(
                        tmpUnmasked,
                        tiled=True,
                        windowed=True,
                        overwrite=False  # never trigger delete-on-open
                    )

        # 2) Mask with template + nodata + compressed final
        with rioxarray.open_rasterio(templatePath, chunks={'x': chunkDims, 'y': chunkDims}, lock=False) as templateRaster:
            with rioxarray.open_rasterio(tmpUnmasked, chunks={'x': chunkDims, 'y': chunkDims}, lock=False) as covariateRaster:
                tb = templateBounds
                cb = covariateRaster.rio.bounds()
                if cb[0] > tb[2] or cb[2] < tb[0] or cb[1] > tb[3] or cb[3] < tb[1]:
                    raise ValueError(f"The extent of the {covId} raster does not overlap the template raster.")

                covariate_nodata  = covariateRaster.rio.nodata
                template_nodata = templateRaster.rio.nodata

                # Replace source nodata → common nodata
                maskedCovariateRaster = covariateRaster.where(covariateRaster != covariate_nodata, other=nodata_value)
                # Mask out template nodata → common nodata
                maskedCovariateRaster = maskedCovariateRaster.where(templateRaster != template_nodata, other=nodata_value)

                # Mask and set no data value
                # maskStack = xr.concat([covariateRaster, templateRaster], dim="band", join="override").chunk(
                #     {'band': -1, 'x': chunkDims, 'y': chunkDims})
                # maskedCovariateRaster = maskStack.map_blocks(mask, kwargs=dict(
                #     input_nodata=covariateRaster.rio.nodata, template_nodata=templateRaster.rio.nodata), template=covariateRaster)
                
                maskedCovariateRaster = maskedCovariateRaster.astype(signed_dtype)
                maskedCovariateRaster.rio.write_nodata(nodata_value, encoded=True, inplace=True)

                # write unique worker final (no overwrite)
                maskedCovariateRaster.rio.to_raster(
                    tmpFinal,
                    dtype=signed_dtype,
                    nodata=nodata_value,
                    tiled=True,
                    windowed=True,
                    overwrite=False,
                    compress='lzw'
                )

        return True  # success

    # Covariate processing with unique temps + locked writes
    tasks = []

    for i in range(len(covariateDataSheet.CovariatesID)):
        ps.environment.progress_bar("message", message="Processing Covariate: " + covariateDataSheet.CovariatesID[i] + " (" + str(
            i+1) + " of " + str(len(covariateDataSheet.CovariatesID)) + ") at " + datetime.now().strftime("%H:%M:%S"))
        ps.environment.update_run_log("Processing Covariate: " + covariateDataSheet.CovariatesID[i] + " (" + str(
            i+1) + " of " + str(len(covariateDataSheet.CovariatesID)) + ") at " + datetime.now().strftime("%H:%M:%S"))

        # Determine input and output file paths
        covId = covariateDataSheet.CovariatesID[i]
        covariatePath = covariateDataSheet.RasterFilePath[i]
        outputCovariatePath = os.path.join(
            ssimTempDir, "processed_" + os.path.basename(covariatePath))

        # check raster data type and set to signed type if necessary
        with rioxarray.open_rasterio(covariatePath, chunks=chunkDims) as covariateRaster:            
            src_dtype = covariateRaster.dtype
            if np.issubdtype(src_dtype, np.unsignedinteger):
                # Avoid overflow when introducing negative nodata
                if src_dtype == np.uint8:
                    signed_dtype = "int16"
                elif src_dtype == np.uint16:
                    signed_dtype = "int32"
                elif src_dtype == np.uint32:
                    signed_dtype = "int64"
                else:
                    # Fallback to float32 when integer range may exceed int64/driver limits
                    signed_dtype = "float32"
            else:
                signed_dtype = str(src_dtype)

            # Decide which resample method to use based on pixel size
            covPixelSize = np.abs(covariateRaster.rio.resolution()).prod()

        # if covariate resolution is finer then template use aggregate method (if coarser use resample method)
        if covPixelSize < templatePixelSize:
            rioResampleMethod = Resampling[covariateDataSheet.rioAggregate[i]]
            dropResAggCol = "ResampleMethod"
        else:
            rioResampleMethod = Resampling[covariateDataSheet.rioResample[i]]
            dropResAggCol = "AggregationMethod"
            
        # tmpUnmasked = uniqueTmp(ssimTempDir, f"tmp_unmasked_{covId}")
        tmpFinal = uniqueTmp(ssimTempDir, f"tmp_final_{covId}")

        pred = predictor_for_dtype(signed_dtype)

        task = delayed(process_one_covariate_worker_fast)(
            covariatePath, covId, tmpFinal,
            templatePath, templateBounds, templateCRS, templateTransform, templateHeight, templateWidth,
            rioResampleMethod, warpMemoryLimit,
            signed_dtype, nodata_value,
            blockDims=chunkDims, compress="lzw", predictor=pred, bigtiff="IF_NEEDED"
        )

         # bind values for later finalization
        def finalize(tmpFinal, covId, dropResAggCol, outputCovariatePath, dep=None):
            return {"CovariatesID": covId, "dropCol": dropResAggCol,
                    "finalTemp": tmpFinal, "finalPath": outputCovariatePath}

        tasks.append(delayed(finalize)(tmpFinal, covId, dropResAggCol, outputCovariatePath, task))

        # task = delayed(process_one_covariate_worker_fast)(
        #     covariatePath, covId, tmpUnmasked, tmpFinal,
        #     templatePath, templateBounds, templateCRS, templateTransform, templateHeight, templateWidth,
        #     rioResampleMethod, chunkDims, warpMemoryLimit,
        #     signed_dtype, nodata_value
        # )

        # define post-run function to clean temp and return mapping (also bound args!)
        # def post_cleanup(tmpUnmasked, tmpFinal, covId, dropResAggCol, outputCovariatePath, dep=None):
        #     # safeUnlink(tmpUnmasked)  # tolerate PermissionError now
        #     return {"CovariatesID": covId, "dropCol": dropResAggCol,
        #             "finalTemp": tmpFinal, "finalPath": outputCovariatePath}

        # tasks.append( delayed(post_cleanup)(tmpUnmasked, tmpFinal, covId, dropResAggCol, outputCovariatePath, task) )

    # Execute in parallel & collect results
    futures = client.compute(tasks)
    results = client.gather(futures)

    # Append to output datasheet
    for r in results:
        tmpFinal  = r["finalTemp"]
        finalPath = r["finalPath"]

        writeReplaceWithRetry(tmpFinal, finalPath)

        # update datasheet row to use final_path
        outputRow = covariateDataSheet.loc[covariateDataSheet.CovariatesID == r["CovariatesID"]].iloc[0, 0:4].copy()
        outputRow.RasterFilePath = finalPath
        outputRow[r["dropCol"]] = float('nan')
        outputCovariateSheet = pd.concat([outputCovariateSheet, pd.DataFrame([outputRow])], ignore_index=True)
        ps.environment.progress_bar()

    # Save updated covariate data to scenario
    ps.environment.progress_bar(
        "message", message="Saving covariate datasheet at " + datetime.now().strftime("%H:%M:%S"))
    ps.environment.update_run_log(
        "Saving covariate datasheet at " + datetime.now().strftime("%H:%M:%S"))
    myScenario.save_datasheet(
        name="wisdm_CovariateData", data=outputCovariateSheet)

    # Prepare restriction raster -------------------------------------------------- TO DO: update this section

    # check if restriction raster was provided
    if pd.notnull(restrictionPath):
        ps.environment.progress_bar(
            "message", message="Processing restriction raster at " + datetime.now().strftime("%H:%M:%S"))
        ps.environment.update_run_log(
            "Processing restriction raster at " + datetime.now().strftime("%H:%M:%S"))

        # Load raster, determine input and output paths
        outputRestrictionPath = os.path.join(
            ssimTempDir, "processed_" + os.path.basename(restrictionPath))

        with rioxarray.open_rasterio(restrictionPath, chunks=True) as restrictionRaster:
            # Check for valid crs
            if pd.isnull(restrictionRaster.rio.crs) or not restrictionRaster.rio.crs.is_valid:
                raise ValueError(
                    "The restriction rasters has an invalid or unknown CRS. Ensure that the raster has a valid CRS before continuing.")

            # check raster data type and set to signed type if necessary
            if np.issubdtype(restrictionRaster.dtype, np.unsignedinteger):
                src_dtype = restrictionRaster.dtype
                bitdepth = src_dtype.itemsize * 8
                signed_dtype = 'int16' if bitdepth <= 16 else 'int32'
            else:
                signed_dtype = str(restrictionRaster.dtype)

            # Decide which resample method to use based on pixel size
            resPixelSize = np.abs(restrictionRaster.rio.resolution()).prod()

            # if restriction resolution is finer then template use aggregate method (if coarser use resample method)
            if resPixelSize < templatePixelSize:
                rioResampleMethod = Resampling["nearest"]
            else:
                rioResampleMethod = Resampling["average"]

        try:

            # Setup client for processing
            with Client(threads_per_worker=num_threads, n_workers=1, processes=False) as client:

                # Connect to restriction raster file
                with rasterio.open(restrictionPath) as src:

                    # Calculate reprojection using a warped virtual layer
                    with WarpedVRT(src, resampling=rioResampleMethod, crs=templateCRS, transform=templateTransform, height=templateHeight, width=templateWidth, warp_mem_limit=warpMemoryLimit, warp_extras={'NUM_THREADS': num_threads}) as vrt:

                        # Convert to rioxarray
                        with rioxarray.open_rasterio(vrt, chunks={'x': chunkDims, 'y': chunkDims}, lock=False) as restrictionRaster:

                            # Check that restriction layer overlaps template
                            # - Note that bounds are ordered xmin, ymin, xmax, ymax
                            tb = templateBounds
                            rb = restrictionRaster.rio.bounds()
                            if rb[0] > tb[0] or rb[1] > tb[1] or rb[2] < tb[2] or rb[3] < tb[3]:
                                raise ValueError(
                                    "The extent of the restriction raster does not overlap the full extent of the template raster. "
                                    "Ensure the restriction raster overlaps the template extent before continuing.")

                            # Write reprojected and clipped layer to uncompresed temp file
                            # - Note! Some resampling algs fail when you try to write compressed, make sure you leave this temp file uncompresed
                            restrictionRaster.rio.to_raster(
                                unmaskedTempPath, tiled=True, lock=SerializableLock(), windowed=True, overwrite=True)

            # Restart client to avoid lock conflicts
            with Client(threads_per_worker=num_threads, n_workers=1, processes=False) as client:

                # Open template within client
                with rioxarray.open_rasterio(templatePath, chunks={'x': chunkDims, 'y': chunkDims}, lock=False) as templateRaster:

                    with rioxarray.open_rasterio(unmaskedTempPath, chunks={'x': chunkDims, 'y': chunkDims}, lock=False) as restrictionRaster:

                        # Mask and set no data value
                        maskedRestrictionRaster = xr.concat([restrictionRaster, templateRaster], "band").chunk({'band': -1, 'x': chunkDims, 'y': chunkDims}).map_blocks(
                            mask, kwargs=dict(input_nodata=restrictionRaster.rio.nodata, template_nodata=templateRaster.rio.nodata), template=restrictionRaster)
                        maskedRestrictionRaster.rio.write_nodata(
                            nodata_value, encoded=True, inplace=True)
                        maskedRestrictionRaster = maskedRestrictionRaster.astype(
                            signed_dtype)

                        # Write to disk
                        maskedRestrictionRaster.rio.to_raster(outputRestrictionPath, dtype=signed_dtype, nodata=nodata_value, tiled=True, lock=SerializableLock(
                        ), windowed=True, overwrite=True, compress='lzw')

        # Tornado's ioloop.py occassionally throws an attribute error looking for an f_code attribute, but the output is still produced correctly
        except AttributeError:
            pass

        # Update restriction raster in output dataframe
        restrictionRasterSheet.RasterFilePath = outputRestrictionPath

        # Save updated datasheet to scenario
        myScenario.save_datasheet(
            name="wisdm_RestrictionRaster", data=restrictionRasterSheet)

    # update progress bar
    ps.environment.progress_bar(report_type="end")
    ps.environment.progress_bar(
        "message", message="Done at " + datetime.now().strftime("%H:%M:%S"))
    ps.environment.update_run_log(
        "Done at " + datetime.now().strftime("%H:%M:%S"))


if __name__ == "__main__":
    prep_spatial_data()
