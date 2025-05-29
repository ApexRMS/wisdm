## ---------------------------------
## wisdm - spatial data preparation
## ApexRMS, March 2024
## ---------------------------------

# built under Python version 3.11.0 & SyncroSim version 3.0.0
# Script pulls in template and covariate rasters and projects, aggragates/resamples, and 
# clips (PARC) covariate rasters to match template; processes field data to ensure sites 
# are in the template CRS and extent; aggregates or weights sites by spatial distribution;
# extracts site-specific covaraite data and generates datasheet of covariate site data.

#def prep_spatial_data():
    #%% Source dependencies ----------------------------------------------------------

## Modify os path if multiple GDAL installations ----
def prep_spatial_data():

    #%%
    import os
    import glob
    from win32api import GetFileVersionInfo, LOWORD, HIWORD

    gdal_installations = []
    if "PATH" in os.environ:
        for p in os.environ["PATH"].split(os.pathsep):
            if p and glob.glob(os.path.join(p, "gdal*.dll")):
                gdal_installations.append(os.path.abspath(p))

    if len(gdal_installations) > 0:
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
    # from osgeo import gdal
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
    # import spatialUtils

    from dask.distributed import Client, Lock
    from dask.utils import SerializableLock
    from shapely.geometry import Point #, shape
    from rasterio.enums import Resampling #, MergeAlg
    from rasterio.vrt import WarpedVRT

    ps.environment.update_run_log('2 - Spatial Data Preparation => Begin')
    ps.environment.progress_bar("message", message = "Preparing inputs for spatial data prep...")

    # Suppress pandas Setting with Copy warning
    pd.options.mode.chained_assignment = None

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
    # ssimTempDir = myLibrary.info["Value"][myLibrary.info.Property == "Temporary files:"].item() # ps.runtime_temp_folder("DataTransfer")
    ssimTempDir = ps.runtime_temp_folder(os.path.join("DataTransfer", "Scenario-" + str(myScenario.sid)))

    # Load datasheets
    # inputs
    networkSheet = myScenario.library.datasheets("wisdm_Network")
    covariatesSheet = myScenario.project.datasheets("wisdm_Covariates")
    covariateDataSheet = myScenario.datasheets("wisdm_CovariateData", show_full_paths=True)
    templateRasterSheet = myScenario.datasheets("wisdm_TemplateRaster", show_full_paths=True)
    restrictionRasterSheet = myScenario.datasheets("wisdm_RestrictionRaster", show_full_paths=True)
    multiprocessingSheet = myScenario.datasheets("core_Multiprocessing")

    # outputs
    outputCovariateSheet = myScenario.datasheets("wisdm_CovariateData", empty = True)

    #%% Set progress bar ---------------------------------------------------------

    steps = 3 + len(covariateDataSheet.CovariatesID)
    ps.environment.progress_bar(report_type = "begin", total_steps = steps)

    #%% Set up dask client -------------------------------------------------------
    if multiprocessingSheet.EnableMultiprocessing.item() == "Yes":
        num_threads = multiprocessingSheet.MaximumJobs.item()
    else:
        num_threads = 1

    # Note: Follow link in output to view progress
    dask.config.set(**{'temporary-directory': os.path.join(ssimTempDir, 'dask-worker-space'),
                       'distributed.scheduler.worker-ttl': None}, 
                       scheduler='threads', serializers=['dask'], deserializers=['dask'])
    # client = Client(threads_per_worker = num_threads, n_workers = 1, processes=False)
    # client.dashboard_link

    #%% Check inputs and set defaults ---------------------------------------------

    # Set PROJ network connection
    if networkSheet.NetworkEnabled.item() == "No":
        pyproj.network.set_network_enabled(active=False)
        # os.environ["PROJ_NETWORK"] = "OFF"
        # pyproj.network.is_network_enabled()

    # Check that a template raster was provided
    if pd.isnull(templateRasterSheet.RasterFilePath.item()):
        raise ValueError("Template raster is missing.")

    # Identify categorical variables 
    catCovs = covariatesSheet.query('IsCategorical == "Yes"').CovariateName.tolist()

    # Ensure that covariate names do not contain spaces
    if any(covariateDataSheet.CovariatesID.str.contains(" ")):
        raise ValueError("Covariate names cannot contain spaces.")

    # Convert Resample and Aggregation method columns to strings
    covariateDataSheet = covariateDataSheet.astype({'ResampleMethod': 'string', 
                                                    "AggregationMethod" : "string"})

    # Set resample defaults
    for i in range(len(covariateDataSheet.ResampleMethod)):
        if pd.isnull(covariateDataSheet.ResampleMethod[i]):
            if covariateDataSheet.CovariatesID[i] in catCovs:
                covariateDataSheet.loc[i, "ResampleMethod"] = "Nearest Neighbor"
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
    resampleAggregateMethodName = ["Nearest Neighbor", "Bilinear", "Cubic", "Cubic Spline", "Lanczos", "Mean", "Min", "Max", "Majority"]
    resampleAggregateMethodCodes = ["nearest", "bilinear", "cubic", "cubic_spline", "lanczos", "average", "min", "max", "mode"]
    covariateDataSheet["rioResample"] = covariateDataSheet.ResampleMethod.replace(resampleAggregateMethodName, resampleAggregateMethodCodes)
    covariateDataSheet["rioAggregate"] = covariateDataSheet.AggregationMethod.replace(resampleAggregateMethodName, resampleAggregateMethodCodes)

    # check if field data was provided
    # if len(fieldDataSheet) == 0:
    #     # raise warning if field data was not provided
    #     ps.environment.update_run_log("Field data was not provided. Site data was not prepared, only the covaraite layers were prepared.")
        

    #%% Load template raster ----------------------------------------------------------------

    # set defualt chunk dimensions and no data value
    warpMemoryLimit = 8192 #MB
    chunkDims = 1024 #1024
    nodata_value = -9999

    templatePath = templateRasterSheet.RasterFilePath.item()
    templateRaster = rioxarray.open_rasterio(templatePath, chunks={'x': chunkDims, 'y': chunkDims}, lock = False)

    # Get information about template
    templateCRS = templateRaster.rio.crs
    if templateCRS.is_valid == False:
        raise ValueError("Template has an invalid CRS (authority code). See {documention} for a list of accepted authority codes.")
    templateTransform = templateRaster.rio.transform()
    templatePixelSize = np.abs(templateRaster.rio.resolution()).prod()
    templateBounds = templateRaster.rio.bounds()
    templateHeight = templateRaster.rio.height
    templateWidth = templateRaster.rio.width

    # update progress bar
    ps.environment.progress_bar()

    # Loop through and "PARC" covariate rasters -------------------------------------------

    # %%First check that all rasters have a valid crs and have valid dimensions
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
    if len(invalidCRS)>0:
        msg = "The following covariate rasters have an invalid or unknown CRS:\n" 
        for i in invalidCRS:
            msg += i + "\n"
        msg += "\nEnsure that the covariate rasters have a valid CRS before continuing."
        ps.environment.update_run_log(msg, type = "status")
        raise ValueError("Covariate rasters have invalid CRS.")     

    #%% Prep covariate rasters

    # Define masking function
    def mask(block, input_nodata, template_nodata):
        # Convert to numpy for faster masking
        npblock = block.to_numpy()

        # Mask and convert all nodata to common value
        temp = np.where(npblock[[0]] == input_nodata, nodata_value, npblock[[0]])
        masked = np.where(npblock[[1]] == template_nodata, nodata_value, temp)

        # Copy data back into appropriate xarray
        return block[[0]].copy(data = masked)

    unmaskedTempPath = os.path.join(ssimTempDir,  "temp.tif")

    # Load and process covariate rasters
    for i in range(len(covariateDataSheet.CovariatesID)):
        ps.environment.progress_bar("message", message = "Processing Covariate: " + covariateDataSheet.CovariatesID[i] + " (" + str(i+1) + " of " + str(len(covariateDataSheet.CovariatesID)) + ") at " +  datetime.now().strftime("%H:%M:%S"))
        ps.environment.update_run_log("Processing Covariate: " + covariateDataSheet.CovariatesID[i] + " (" + str(i+1) + " of " + str(len(covariateDataSheet.CovariatesID)) + ") at " +  datetime.now().strftime("%H:%M:%S"))

        # Determine input and output file paths
        covariatePath = covariateDataSheet.RasterFilePath[i]
        outputCovariatePath = os.path.join(ssimTempDir, "processed_" + os.path.basename(covariatePath))

        # Decide which resample method to use based on pixel size
        with rioxarray.open_rasterio(covariatePath, chunks = chunkDims) as covariateRaster:
            covPixelSize = np.abs(covariateRaster.rio.resolution()).prod()

            # if covariate resolution is finer then template use aggregate method (if coarser use resample method) 
            if covPixelSize < templatePixelSize:
                rioResampleMethod = Resampling[covariateDataSheet.rioAggregate[i]]
                dropResAggCol = "ResampleMethod"
            else:
                rioResampleMethod = Resampling[covariateDataSheet.rioResample[i]]
                dropResAggCol = "AggregationMethod"

        try:

            # Setup client for processing
            with Client(threads_per_worker = num_threads, n_workers = 1, processes=False) as client:

                # Connect to covariate raster file 
                with rasterio.open(covariatePath) as src:

                    # Calculate reprojection using a warped virtual layer
                    with WarpedVRT(src, resampling = rioResampleMethod, crs = templateCRS, transform = templateTransform, height = templateHeight, width = templateWidth, warp_mem_limit= warpMemoryLimit, warp_extras={'NUM_THREADS':num_threads}) as vrt:

                        # Convert to rioxarray
                        with rioxarray.open_rasterio(vrt, chunks = {'x': chunkDims, 'y': chunkDims}, lock = False) as covariateRaster:
                            
                            # Check that covariate layer overlaps template
                            # - Note that bounds are ordered xmin, ymin, xmax, ymax
                            tb = templateBounds
                            cb = covariateRaster.rio.bounds()
                            if cb[0] > tb[0] or cb[1] > tb[1] or cb[2] < tb[2] or cb[3] < tb[3]:
                                raise ValueError(print("The extent of the", covariateDataSheet.CovariatesID[i], 
                                                       "raster does not overlap the full extent of the template raster. Ensure all covariate rasters overlap the template extent before continuing."))

                            # Write reprojected and clipped layer to uncompresed temp file
                            # - Note! Some resampling algs fail when you try to write compressed, make sure you leave this temp file uncompresed
                            covariateRaster.rio.to_raster(unmaskedTempPath, tiled = True, lock = SerializableLock(), windowed=True, overwrite = True)

            # Restart client to avoid lock conflicts
            with Client(threads_per_worker = num_threads, n_workers = 1, processes=False) as client:

                # Open template within clien
                with rioxarray.open_rasterio(templatePath, chunks={'x': chunkDims, 'y': chunkDims}, lock = False) as templateRaster:

                    with rioxarray.open_rasterio(unmaskedTempPath, chunks={'x': chunkDims, 'y': chunkDims}, lock = False) as covariateRaster:
                        # Mask and set no data value
                        maskStack = xr.concat([covariateRaster, templateRaster], dim = "band", join="override").chunk({'band':-1, 'x':chunkDims, 'y':chunkDims})
                        maskedCovariateRaster = maskStack.map_blocks(mask, kwargs=dict(input_nodata=covariateRaster.rio.nodata, template_nodata=templateRaster.rio.nodata), template = covariateRaster)
                        maskedCovariateRaster.rio.write_nodata(nodata_value, encoded=True, inplace=True)

                        # Write to disk
                        maskedCovariateRaster.rio.to_raster(outputCovariatePath, tiled = True, lock = SerializableLock(), windowed=True, overwrite = True, compress = 'lzw')
    
        # Tornado's ioloop.py occassionally throws an attribute error looking for an f_code attribute, but the output is still produced correctly      
        except Exception as e:
            ps.environment.update_run_log("Caught exception processing covariate: " + covariateDataSheet.CovariatesID[i] + " -- " + str(type(e)) + " -- " + str(e))


        # Add covariate data to output dataframe
        outputRow = covariateDataSheet.iloc[i, 0:4]
        outputRow.RasterFilePath = outputCovariatePath
        outputRow[dropResAggCol] = float('nan')

        outputCovariateSheet = pd.concat([outputCovariateSheet, pd.DataFrame([outputRow])])
        
        # update progress bar
        ps.environment.progress_bar()

    #%% Save updated covariate data to scenario 

    ps.environment.progress_bar("message", message = "Saving covariate datasheet at " + datetime.now().strftime("%H:%M:%S"))
    ps.environment.update_run_log("Saving covariate datasheet at " + datetime.now().strftime("%H:%M:%S"))
    myScenario.save_datasheet(name="wisdm_CovariateData", data=outputCovariateSheet) 

    #%% Prepare restriction raster --------------------------------------------------

    # check if restriction raster was provided
    if len(restrictionRasterSheet.RasterFilePath) > 0:
        ps.environment.progress_bar("message", message = "Processing restriction raster at " + datetime.now().strftime("%H:%M:%S"))
        ps.environment.update_run_log("Processing restriction raster at " + datetime.now().strftime("%H:%M:%S"))

        # Load raster, determine input and output paths
        restrictionPath = restrictionRasterSheet.RasterFilePath.item()
        outputRestrictionPath = os.path.join(ssimTempDir, os.path.basename(restrictionRasterSheet.RasterFilePath.item()))
        
        with rioxarray.open_rasterio(restrictionPath, chunks = True) as restrictionRaster:
            # Check for valid crs
            if pd.isnull(restrictionRaster.rio.crs):
                raise ValueError(print("The restriction rasters has an invalid or unknown CRS. Ensure that the raster has a valid CRS before continuing."))  
            elif restrictionRaster.rio.crs.is_valid == False:
                raise ValueError(print("The restriction rasters has an invalid or unknown CRS. Ensure that the raster has a valid CRS before continuing."))  
            
            # Decide which resample method to use based on pixel size
            resPixelSize = np.abs(restrictionRaster.rio.resolution()).prod()

            # if restriction resolution is finer then template use aggregate method (if coarser use resample method) 
            if resPixelSize < templatePixelSize:
                rioResampleMethod = Resampling["nearest"]
            else:
                rioResampleMethod = Resampling["average"]

        try:

            # Setup client for processing
            with Client(threads_per_worker = num_threads, n_workers = 1, processes=False) as client:

                # Connect to restriction raster file
                with rasterio.open(restrictionPath) as src:

                    # Calculate reprojection using a warped virtual layer
                    with WarpedVRT(src, resampling = rioResampleMethod, crs = templateCRS, transform = templateTransform, height = templateHeight, width = templateWidth, warp_mem_limit= warpMemoryLimit, warp_extras={'NUM_THREADS':num_threads}) as vrt:

                        # Convert to rioxarray
                        with rioxarray.open_rasterio(vrt, chunks = {'x': chunkDims, 'y': chunkDims}, lock = False) as restrictionRaster:
                            
                            # Check that restriction layer overlaps template
                            # - Note that bounds are ordered xmin, ymin, xmax, ymax
                            tb = templateBounds
                            rb = restrictionRaster.rio.bounds()
                            if rb[0] > tb[0] or rb[1] > tb[1] or rb[2] < tb[2] or rb[3] < tb[3]:
                                raise ValueError(print("The extent of the restriction raster does not overlap the full extent of the template raster. Ensure the restiction raster overlaps the template extent before continuing."))
                            
                            # Write reprojected and clipped layer to uncompresed temp file
                            # - Note! Some resampling algs fail when you try to write compressed, make sure you leave this temp file uncompresed
                            restrictionRaster.rio.to_raster(unmaskedTempPath, tiled = True, lock = SerializableLock(), windowed=True, overwrite = True)
            
            # Restart client to avoid lock conflicts
            with Client(threads_per_worker = num_threads, n_workers = 1, processes=False) as client:

                # Open template within clien
                with rioxarray.open_rasterio(templatePath, chunks={'x': chunkDims, 'y': chunkDims}, lock = False) as templateRaster:

                    with rioxarray.open_rasterio(unmaskedTempPath, chunks={'x': chunkDims, 'y': chunkDims}, lock = False) as covariateRaster:

                            # Mask and set no data value
                            maskedRestrictionRaster = xr.concat([restrictionRaster, templateRaster], "band").chunk({'band':-1, 'x':chunkDims, 'y':chunkDims}).map_blocks(mask, kwargs=dict(input_nodata=restrictionRaster.rio.nodata, template_nodata=templateRaster.rio.nodata), template = restrictionRaster)
                            maskedRestrictionRaster.rio.write_nodata(nodata_value, encoded=True, inplace=True)

                                # Write to disk
                            maskedRestrictionRaster.rio.to_raster(outputRestrictionPath, tiled = True, lock = SerializableLock(), windowed=True, overwrite = True, compress = 'lzw')

        # Tornado's ioloop.py occassionally throws an attribute error looking for an f_code attribute, but the output is still produced correctly
        except AttributeError:
            pass
            
        # Update restriction raster in output dataframe
        restrictionRasterSheet.RasterFilePath = outputRestrictionPath

        # Save updated datasheet to scenario 
        myScenario.save_datasheet(name="wisdm_RestrictionRaster", data=restrictionRasterSheet)    

    # update progress bar
    ps.environment.progress_bar(report_type="end")
    ps.environment.progress_bar("message", message = "Done at " + datetime.now().strftime("%H:%M:%S"))
    ps.environment.update_run_log("Done at " + datetime.now().strftime("%H:%M:%S"))
    #%%
    # if __name__ == "__main__":
    #     prep_spatial_data()
    #%%

if __name__ == "__main__":
    prep_spatial_data()
    # prep_s