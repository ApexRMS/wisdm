## ---------------------------------
## wisdm - spatial data preparation
## ApexRMS, October 2022
## ---------------------------------

# built under Python version 3.10.6 & SyncroSim version 2.4.9
# Script pulls in template and covariate rasters and projects, aggragates/resamples, and 
# clips (PARC) covariate rasters to match template; processes field data to ensure sites 
# are in the template CRS and extent; aggregates or weights sites by spatial distribution;
# extracts site-specific covaraite data and generates datasheet of covariate site data.

def prep_spatial_data():
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
    import pysyncrosim as ps     
    import numpy as np          
    import pandas as pd          
    import rioxarray
    import xarray as xr
    import rasterio
    import geopandas as gpd
    import shapely
    import dask
    import pyproj
    from datetime import datetime
    # import spatialUtils

    from dask.distributed import Client, Lock
    from shapely.geometry import Point #, shape
    from rasterio.enums import Resampling #, MergeAlg
    from rasterio.vrt import WarpedVRT

    ps.environment.progress_bar("message", message = "Preparing inputs for spatial data prep...")

    # Suppress pandas Setting with Copy warning
    pd.options.mode.chained_assignment = None

    ## Modify the os PROJ path (when running with Conda) ----
    myLibrary = ps.Library()
    mySession = ps.Session()

    result = mySession._Session__call_console(["--conda", "--config"])
    conda_fpath = result.stdout.decode('utf-8').strip().split(": ")[1]
    if myLibrary.datasheets("core_Options").UseConda.item() == "Yes":
        os.environ["PROJ_DATA"] = os.path.join(conda_fpath , "envs\\wisdm\\wisdm-py-conda\\Library\\share\\proj")
        os.environ['PROJ_CURL_CA_BUNDLE'] = os.path.join(conda_fpath , "envs\\wisdm\\wisdm-py-conda\\Library\\ssl\\cacert.pem")

    # if myLibrary.datasheets("core_Options").UseConda.item() == "Yes":
    #    os.environ["PROJ_DATA"] = os.path.join(mySession.conda_filepath, "envs\\wisdm\\wisdm-py-conda\\Library\\share\\proj")
    #    os.environ['PROJ_CURL_CA_BUNDLE'] = os.path.join(mySession.conda_filepath, "envs\\wisdm\\wisdm-py-conda\\Library\\ssl\\cacert.pem")
        # pyproj.datadir.set_data_dir(os.path.join(mySession.conda_filepath, "envs\\wisdm\\wisdm-py-conda\\Library\\share\\proj"))
        # pyproj.network.set_ca_bundle_path(os.path.join(mySession.conda_filepath, "envs\\wisdm\\wisdm-py-conda\\Library\\ssl\\cacert.pem"))
        

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
    covariatesSheet = myScenario.project.datasheets("Covariates")
    covariateDataSheet = myScenario.datasheets("CovariateData", show_full_paths=False)
    # fieldDataSheet = myScenario.datasheets("FieldData")
    # fieldDataOptions = myScenario.datasheets("FieldDataOptions")
    templateRasterSheet = myScenario.datasheets("TemplateRaster", show_full_paths=False)
    restrictionRasterSheet = myScenario.datasheets("RestrictionRaster", show_full_paths=True)
    multiprocessingSheet = myScenario.datasheets("core_Multiprocessing")

    # outputs
    outputCovariateSheet = myScenario.datasheets("CovariateData", empty = True)

    #%% Set progress bar ---------------------------------------------------------

    steps = 3 + len(covariateDataSheet.CovariatesID)
    ps.environment.progress_bar(report_type = "begin", total_steps = steps)

    #%% Set up dask client -------------------------------------------------------
    if multiprocessingSheet.EnableMultiprocessing.item() == "Yes":
        num_threads = multiprocessingSheet.MaximumJobs.item()
    else:
        num_threads = 1

    # Note: Follow link in output to view progress
    dask.config.set(**{'temporary-directory': os.path.join(ssimTempDir, 'dask-worker-space')})

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

    # Set resample defaults
    for i in range(len(covariateDataSheet.ResampleMethod)):
        if pd.isnull(covariateDataSheet.ResampleMethod[i]):
            if covariateDataSheet.CovariatesID[i] in catCovs:
                covariateDataSheet.ResampleMethod[i] = "Nearest Neighbor"
            else:
                covariateDataSheet.ResampleMethod[i] = "Bilinear"  

    # Set aggregate defaults
    for i in range(len(covariateDataSheet.AggregationMethod)):
        if pd.isnull(covariateDataSheet.AggregationMethod[i]):
            if covariateDataSheet.CovariatesID[i] in catCovs:
                covariateDataSheet.AggregationMethod[i] = "Majority"
            else:
                covariateDataSheet.AggregationMethod[i] = "Mean" 

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
    warpMemoryLimit = 8192 # 1024 #  #MB
    chunkDims = 8192 # 1024 # 
    nodata_value = -9999

    templatePath = ssimInputDir + "\\wisdm_TemplateRaster\\" + templateRasterSheet.RasterFilePath.item()
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

    # update default chunk size based on template size
    # if templateHeight < 1000 or templateWidth < 1000:
    #     warpMemoryLimit = 512 #MB
    #     chunkDims = 512
    #     templateRaster = rioxarray.open_rasterio(templatePath, chunks={'x': chunkDims, 'y': chunkDims}, lock = False)
    
    # if templateHeight > 10000 or templateWidth > 10000:
    #     warpMemoryLimit = 8192 #MB
    #     chunkDims = 8192
    #     templateRaster = rioxarray.open_rasterio(templatePath, chunks={'x': chunkDims, 'y': chunkDims}, lock = False)

    # update progress bar
    ps.environment.progress_bar()

    # Loop through and "PARC" covariate rasters -------------------------------------------

    # %%First check that all rasters have a valid crs
    invalidCRS = []
    for i in range(len(covariateDataSheet.CovariatesID)):
        # Load covariate rasters
        covariatePath = ssimInputDir + "\\wisdm_CovariateData\\" + covariateDataSheet.RasterFilePath[i]
        covariateRaster = rioxarray.open_rasterio(covariatePath, chunks=True)
        if pd.isnull(covariateRaster.rio.crs):
            invalidCRS.append(covariateDataSheet.CovariatesID[i])
        elif covariateRaster.rio.crs.is_valid == False:
            invalidCRS.append(covariateDataSheet.CovariatesID[i]) 
    if len(invalidCRS)>0:
        raise ValueError(print("The following covariate rasters have an invalid or unknown CRS:", *invalidCRS, "Ensure that the covariate rasters have a valid CRS before continuing.", sep="\n") )      

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
        covariatePath = ssimInputDir + "\\wisdm_CovariateData\\" + covariateDataSheet.RasterFilePath[i]
        outputCovariatePath = os.path.join(ssimTempDir, os.path.basename(covariatePath))

        # Decide which resample method to use based on pixel size
        with rioxarray.open_rasterio(covariatePath, chunks = True) as covariateRaster:
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
                                raise ValueError(print("The extent of the", covariateDataSheet.CovariatesID[i], "raster does not overlap the full extent of the template raster. Ensure all covariate rasters overlap the template extent before continuing."))

                            # Write reprojected and clipped layer to uncompresed temp file
                            # - Note! Some resampling algs fail when you try to write compressed, make sure you leave this temp file uncompresed
                            covariateRaster.rio.write_nodata(nodata_value, inplace=True)
                            covariateRaster.rio.to_raster(unmaskedTempPath, tiled = True, lock = Lock(client = client), windowed=True, overwrite = True)

            # Restart client to avoid lock conflicts
            with Client(threads_per_worker = num_threads, n_workers = 1, processes=False) as client:

                # Open template within clien
                with rioxarray.open_rasterio(templatePath, chunks={'x': chunkDims, 'y': chunkDims}, lock = False) as templateRaster:

                    with rioxarray.open_rasterio(unmaskedTempPath, chunks={'x': chunkDims, 'y': chunkDims}, lock = False) as covariateRaster:
                        # Mask and set no data value
                        maskStack = xr.concat([covariateRaster, templateRaster], dim = "band", join="override").chunk({'band':-1, 'x':chunkDims, 'y':chunkDims})
                        maskedCovariateRaster = maskStack.map_blocks(mask, kwargs=dict(input_nodata=covariateRaster.rio.nodata, template_nodata=templateRaster.rio.nodata), template = covariateRaster)
                        maskedCovariateRaster.rio.write_nodata(nodata_value, inplace=True)

                        # Write to disk
                        maskedCovariateRaster.rio.to_raster(outputCovariatePath, tiled = True, lock = True, windowed=True, overwrite = True, compress = 'lzw')
        
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
    myScenario.save_datasheet(name="CovariateData", data=outputCovariateSheet) 

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
                            restrictionRaster.rio.to_raster(unmaskedTempPath, tiled = True, lock = Lock(client = client), windowed=True, overwrite = True)
            
            # Restart client to avoid lock conflicts
            with Client(threads_per_worker = num_threads, n_workers = 1, processes=False) as client:

                # Open template within clien
                with rioxarray.open_rasterio(templatePath, chunks={'x': chunkDims, 'y': chunkDims}, lock = False) as templateRaster:

                    with rioxarray.open_rasterio(unmaskedTempPath, chunks={'x': chunkDims, 'y': chunkDims}, lock = False) as covariateRaster:

                            # Mask and set no data value
                            maskedRestrictionRaster = xr.concat([restrictionRaster, templateRaster], "band").chunk({'band':-1, 'x':chunkDims, 'y':chunkDims}).map_blocks(mask, kwargs=dict(input_nodata=restrictionRaster.rio.nodata, template_nodata=templateRaster.rio.nodata), template = restrictionRaster)
                            maskedRestrictionRaster.rio.write_nodata(nodata_value, inplace=True)

                                # Write to disk
                            maskedRestrictionRaster.rio.to_raster(outputRestrictionPath, tiled = True, lock = True, windowed=True, overwrite = True, compress = 'lzw')

        # Tornado's ioloop.py occassionally throws an attribute error looking for an f_code attribute, but the output is still produced correctly
        except AttributeError:
            pass
            
        # Add covariate data to output dataframe
        restrictionRasterSheet.RasterFilePath = outputRestrictionPath

        # Save updated covariate data to scenario 
        myScenario.save_datasheet(name="RestrictionRaster", data=restrictionRasterSheet)     

    # update progress bar
    ps.environment.progress_bar(report_type="end")
    ps.environment.progress_bar("message", message = "Done at " + datetime.now().strftime("%H:%M:%S"))
    ps.environment.update_run_log("Done at " + datetime.now().strftime("%H:%M:%S"))
#%%
if __name__ == "__main__":
    prep_spatial_data()