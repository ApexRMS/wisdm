## ---------------------------------
## wisdm - spatial data preparation
## ApexRMS, October 2022
## ---------------------------------

# built under Python version 3.10.6 & SyncroSim version 2.4.9
# Script pulls in template and covariate rasters and projects, aggragates/resamples, and 
# clips (PARC) covariate rasters to match template; processes field data to ensure sites 
# are in the template CRS and extent; aggregates or weights sites by spatial distribution;
# extracts site-specific covaraite data and generates datasheet of covariate site data.

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
import pandas as pd          
import rioxarray
import xarray
import rasterio
import geopandas as gpd
import shapely
import dask
import pyproj
# import spatialUtils

from dask.distributed import Client, Lock
from shapely.geometry import Point #, shape
from rasterio.enums import Resampling #, MergeAlg

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
client = Client(threads_per_worker = num_threads, n_workers = 1, processes=False)
# client


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

# set defualt chunk dimensions
chunkDims = 4096 #1024

templatePath = ssimInputDir + "\\wisdm_TemplateRaster\\" + templateRasterSheet.RasterFilePath.item()
templateRaster = rioxarray.open_rasterio(templatePath, chunks={'x': chunkDims, 'y': chunkDims}, masked=True)
templateMask = templateRaster.to_masked_array().mask[0]

# Get information about template
templateCRS = templateRaster.rio.crs
if templateCRS.is_valid == False:
    raise ValueError("Template has an invalid CRS (authority code). See {documention} for a list of accepted authority codes.")
templateResolution = templateRaster.rio.resolution()
templateExtent = list(templateRaster.rio.bounds())
templateTransform = templateRaster.rio.transform()
templatePixelSize = templateResolution[0]*-templateResolution[1]

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
    else:
        if covariateRaster.rio.crs.is_valid == False:
           invalidCRS.append(covariateDataSheet.CovariatesID[i]) 
if len(invalidCRS)>0:
    raise ValueError(print("The following covariate rasters have an invalid or unknown CRS:", *invalidCRS, "Ensure that the covariate rasters have a valid CRS before continuing.", sep="\n") )      

#%% Prep covariate rasters
for i in range(len(covariateDataSheet.CovariatesID)):

    # Load and process covariate rasters
    covariatePath = ssimInputDir + "\\wisdm_CovariateData\\" + covariateDataSheet.RasterFilePath[i]
    covariateRaster = rioxarray.open_rasterio(covariatePath, chunks={'x': chunkDims, 'y': chunkDims})        
    outputCovariatePath = os.path.join(ssimTempDir, os.path.basename(covariatePath))

    # check if covariate extent fully overlaps template extent: [left, bottom, right, top]
    covExtent = list(rasterio.warp.transform_bounds(covariateRaster.rio.crs, templateCRS,*covariateRaster.rio.bounds()))
    
    overlap=[]
    for j in range(2):
        overlap.append(covExtent[j] <= templateExtent[j])
    for j in range(2,4):
        overlap.append(covExtent[j] >= templateExtent[j])
    if any(overlap) == False:
        raise ValueError(print("The extent of the", covariateDataSheet.CovariatesID[i], "raster does not overlap the full extent of the template raster. Ensure all covariate rasters overlap the template extent before continuing."))
    
    # Determine if raster needs to be resampled or aggregated
    # convert covariate resolution to template units (following code is faster then reprojecting for large rasters)
    covAffine = rasterio.warp.calculate_default_transform(covariateRaster.rio.crs, templateCRS, 
                                                            covariateRaster.rio.width, covariateRaster.rio.height,
                                                            *covariateRaster.rio.bounds())[0]
    
    covPixelSize = covAffine[0]*-covAffine[4]
        
    # if covariate resolution is finer then template use aggregate method (if coarser use resample method) 
    if covPixelSize < templatePixelSize:
        rioResampleMethod = covariateDataSheet.rioAggregate[i]
        dropResAggCol = "ResampleMethod"
    else:
        rioResampleMethod = covariateDataSheet.rioResample[i] 
        dropResAggCol = "AggregationMethod"

    # reproject covariate raster ## TO DO - build out memory save version of reproject-match
    # spatialUtils.parc(inputFile=covariatePath, 
    #                 templateRaster=templateRaster, 
    #                 outputFile=outputCovariatePath,
    #                 chunks = {'x': chunkDims, 'y': chunkDims},
    #                 resampling=Resampling[rioResampleMethod],
    #                 transform=covAffine,
    #                 mask = templatePolygons.geometry,
    #                 client=client)

    # reproject/resample/clip covariate raster to match template 
    covariateRaster = covariateRaster.rio.reproject_match(templateRaster,
                                                            resampling=Resampling[rioResampleMethod])
    
    #Add NoData mask to covariate raster
    covariateRaster = np.ma.masked_array(covariateRaster, templateMask)
    covariateRaster.fill_value = -9999
   
    # Set raster output parameters
    raster_params = {
        'driver': 'GTiff',
        'width': templateRaster.rio.width,
        'height': templateRaster.rio.height,
        'count': 1,
        'dtype': covariateRaster.dtype, 
        "nodata": -9999,
        'compress': 'lzw',
        'crs': templateCRS,
        'transform': templateTransform
    }
        
    # Write processed covariate raster to file with internal NoData mask 
    with rasterio.Env(GDAL_TIFF_INTERNAL_MASK=True):
        with rasterio.open(outputCovariatePath, mode="w", **raster_params,
                            masked= True, tiled=True, windowed=True, overwrite=True) as src:
            src.write(covariateRaster)
    
    # Add covariate data to output dataframe
    outputRow = covariateDataSheet.iloc[i, 0:4]
    outputRow.RasterFilePath = outputCovariatePath
    outputRow[dropResAggCol] = float('nan')

    outputCovariateSheet = pd.concat([outputCovariateSheet, pd.DataFrame([outputRow])])
    
    # update progress bar
    ps.environment.progress_bar()

#%% Save updated covariate data to scenario 
myScenario.save_datasheet(name="CovariateData", data=outputCovariateSheet) 

#%% Prepare restriction raster --------------------------------------------------

# check if restriction raster was provided
if len(restrictionRasterSheet.RasterFilePath) > 0:

    restrictionRaster = rioxarray.open_rasterio(restrictionRasterSheet.RasterFilePath.item(), chunks=True)
    outputPath = os.path.join(ssimTempDir, os.path.basename(restrictionRasterSheet.RasterFilePath.item()))
    
    if pd.isnull(restrictionRaster.rio.crs):
        raise ValueError(print("The restriction rasters has an invalid or unknown CRS. Ensure that the raster has a valid CRS before continuing."))  
    else:
        if restrictionRaster.rio.crs.is_valid == False:
           raise ValueError(print("The restriction rasters has an invalid or unknown CRS. Ensure that the raster has a valid CRS before continuing."))  

     # check if extent fully overlaps template extent: [left, bottom, right, top]
    resExtent = list(rasterio.warp.transform_bounds(restrictionRaster.rio.crs, templateCRS,*restrictionRaster.rio.bounds()))
    
    overlap=[]
    for j in range(2):
        overlap.append(resExtent[j] <= templateExtent[j])
    for j in range(2,4):
        overlap.append(resExtent[j] >= templateExtent[j])
    if any(overlap) == False:
        raise ValueError(print("The extent of the restriction raster does not overlap the full extent of the template raster. Ensure that the restriciton raster overlaps the template extent before continuing."))
    
    # Determine if raster needs to be resampled or aggregated
    # convert resolution to template units (following code is faster then reprojecting for large rasters)
    resAffine = rasterio.warp.calculate_default_transform(restrictionRaster.rio.crs, templateCRS, 
                                                            restrictionRaster.rio.width, restrictionRaster.rio.height,
                                                            *restrictionRaster.rio.bounds())[0]
    
    resPixelSize = resAffine[0]*-resAffine[4]
        
    # if resolution is finer then template use aggregate method (if coarser use resample method) 
    if resPixelSize < templatePixelSize:
        rioResampleMethod = "nearest"
    else:
        rioResampleMethod = "average"

    # reproject covariate raster ## TO DO - build out memory save version of reproject-match
    # spatialUtils.parc(inputFile=covariatePath, 
    #                 templateRaster=templateRaster, 
    #                 outputFile=outputCovariatePath,
    #                 chunks = {'x': chunkDims, 'y': chunkDims},
    #                 resampling=Resampling[rioResampleMethod],
    #                 transform=covAffine,
    #                 mask = templatePolygons.geometry,
    #                 client=client)

    # reproject/resample/clip raster to match template 
    restrictionRaster = restrictionRaster.rio.reproject_match(templateRaster,
                                                            resampling=Resampling[rioResampleMethod])
    
    #Add NoData mask to raster
    restrictionRaster = np.ma.masked_array(restrictionRaster, templateMask)
    restrictionRaster.fill_value = -9999
   
    # Set raster output parameters
    raster_params = {
        'driver': 'GTiff',
        'width': templateRaster.rio.width,
        'height': templateRaster.rio.height,
        'count': 1,
        'dtype': restrictionRaster.dtype, 
        "nodata": -9999,
        'compress': 'lzw',
        'crs': templateCRS,
        'transform': templateTransform
    }
        
    # Write processed raster to file with internal NoData mask 
    with rasterio.Env(GDAL_TIFF_INTERNAL_MASK=True):
        with rasterio.open(outputPath, mode="w", **raster_params,
                            masked= True, tiled=True, windowed=True, overwrite=True) as src:
            src.write(restrictionRaster)
    
    # Add covariate data to output dataframe
    restrictionRasterSheet.RasterFilePath = outputPath

    # Save updated covariate data to scenario 
    myScenario.save_datasheet(name="RestrictionRaster", data=restrictionRasterSheet)     

# update progress bar
ps.environment.progress_bar(report_type="end")
# %%
