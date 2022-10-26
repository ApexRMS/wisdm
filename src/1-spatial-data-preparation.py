#%% Load packages
import pysyncrosim as ps     
import numpy as np          
import pandas as pd          
import os
import rioxarray
import rasterio
import psutil

from rasterio.enums import Resampling

#%% Load the SyncroSim Scenario that is currently running
myScenario = ps.Scenario()  

# Create a temporary folder for storing rasters
ssimTempDir = ps.runtime_temp_folder("Data")

# Get path to scnario inputs 
ssimInputDir = myScenario.library.location + ".input\Scenario-" + str(myScenario.sid)

#%% Load datasheets
# inputs
covariatesSheet = myScenario.project.datasheets("Covariates")
covariateDataSheet = myScenario.datasheets("CovariateData")
# fieldDataSheet = myScenario.datasheets("FieldData")
# fieldDataOptions = myScenario.datasheets("FieldDataOptions")
templateRasterSheet = myScenario.datasheets("TemplateRaster")
# multiprocessingSheet = myScenario.datasheets("core_Multiprocessing")
# spatialMulitprocessingSheet = myScenario.datasheets("corestime_Multiprocessing")

# outputs
outputCovariateSheet = myScenario.datasheets("CovariateData", empty = True)
# siteDataSheet = myScenario.datasheets("SiteData")

#%% Check inputs and set defaults

# check that a template raster was provided
if pd.isnull(templateRasterSheet.RasterFilePath.item()):
    raise ValueError("Template raster is missing.")

# identify categorical variables 
catCovs = covariatesSheet.query('IsCategorical == "Yes"').CovariateName.tolist()

#%% resample defaults
for i in range(len(covariateDataSheet.ResampleMethod)):
    if pd.isnull(covariateDataSheet.ResampleMethod[i]):
        if covariateDataSheet.CovariatesID[i] in catCovs:
            covariateDataSheet.ResampleMethod[i] = "Nearest Neighbor"
        else:
            covariateDataSheet.ResampleMethod[i] = "Bilinear"  

#%% aggregate defaults
for i in range(len(covariateDataSheet.AggregationMethod)):
    if pd.isnull(covariateDataSheet.AggregationMethod[i]):
        if covariateDataSheet.CovariatesID[i] in catCovs:
            covariateDataSheet.AggregationMethod[i] = "Majority"
        else:
            covariateDataSheet.AggregationMethod[i] = "Mean" 

#%% define rio up/down-sample method
resampleAggregateMethodName = ["Nearest Neighbor", "Bilinear", "Cubic", "Cubic Spline", "Lanczos", "Mean", "Min", "Max", "Majority"]
resampleAggregateMethodCodes = ["nearest", "bilinear", "cubic", "cubic_spline", "lanczos", "average", "min", "max", "mode"]
covariateDataSheet["rioResample"] = covariateDataSheet.ResampleMethod.replace(resampleAggregateMethodName, resampleAggregateMethodCodes)
covariateDataSheet["rioAggregate"] = covariateDataSheet.AggregationMethod.replace(resampleAggregateMethodName, resampleAggregateMethodCodes)

#%% Load template raster 
templatePath = ssimInputDir + "\\wisdm_TemplateRaster\\" + templateRasterSheet.RasterFilePath.item()
templateRaster = rioxarray.open_rasterio(templatePath, chunks=True, masked=True)
templateMask = rasterio.open(templatePath).read(1, masked=True).mask
# templateMask = np.invert(templateMask)

#%% get information about template
templateCRS = templateRaster.rio.crs
templateResolution = templateRaster.rio.resolution()
templateExtent = list(templateRaster.rio.bounds())
templateTransform = templateRaster.rio.transform()
templatePixelSize = templateResolution[0]*-templateResolution[1]

#%% Loop through and "PARC" covariate rasters  
for i in range(len(covariateDataSheet.CovariatesID)):

    #%% Load and process covariate rasters
    covariatePath = ssimInputDir + "\\wisdm_CovariateData\\" + covariateDataSheet.RasterFilePath[i]
    covariateRaster = rioxarray.open_rasterio(covariatePath, chunks=True)
    # check if raster mem is larger then availble mem - if yes, assign mem appropriate chunk sizes
    if covariateRaster.nbytes > psutil.virtual_memory()[1]:
        if 1024*1024 < psutil.virtual_memory()[1]:
            chunkDims = 1024
        if 2290*2290 < psutil.virtual_memory()[1]:
            chunkDims = 2290
        covariateRaster = rioxarray.open_rasterio(covariatePath, chunks={'x': chunkDims, 'y': chunkDims})        
    outputCovariatePath = os.path.join(ssimTempDir, os.path.basename(covariatePath))

    #%% check if covariate extent fully overlaps template extent: [left, bottom, right, top]
    covExtent = list(rasterio.warp.transform_bounds(covariateRaster.rio.crs, templateCRS,*covariateRaster.rio.bounds()))
    
    overlap=[]
    for j in range(2):
        overlap.append(covExtent[j] <= templateExtent[j])
    for j in range(2,4):
        overlap.append(covExtent[j] >= templateExtent[j])
    if any(overlap) == False:
        raise ValueError(print("The extent of the", covariateDataSheet.CovariatesID[i], "raster does not overlap the full extent of the template raster. Ensure all covariate rasters overlap the template extent before continuing."))
    
    #%% Determine if raster needs to be resampled or aggregated
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

    #%% reproject/resample/clip covariate raster to match template 
    covariateRaster = covariateRaster.rio.reproject_match(templateRaster,
                                                            resampling=Resampling[rioResampleMethod])
    
    #%% clip raster to extent of valid data in template
    covariateRaster = np.ma.masked_array(covariateRaster, templateMask)
    # if covariateRaster.dtype in ["float_", "float16", "float32", "float64"]:
    #     covariateRaster.fill_value = float("nan")
    # else: 
    covariateRaster.fill_value = -9999
    #%% set raster output parameters
    raster_params = {
        'driver': 'GTiff',
        'width': templateRaster.rio.width,
        'height': templateRaster.rio.height,
        'count': 1,
        'dtype': covariateRaster.dtype,
        'compress': 'lzw',
        'crs': templateCRS,
        'transform': templateTransform
    }
    # %% Save outputs
    # save prepped covariate raster to temp file
    # covariateRaster.rio.write_nodata(-9999, inplace=True)
    # covariateRaster.rio.to_raster(outputCovariatePath, masked= True, tiled=True, windowed=True, overwrite=True, compress="lzw")
    
    #%% add nodata mask to covariate raster
    with rasterio.Env(GDAL_TIFF_INTERNAL_MASK=True):
        with rasterio.open(outputCovariatePath, mode="w", **raster_params,
                            masked= True, tiled=True, windowed=True, overwrite=True) as src:
            src.write(covariateRaster)
    
    #%% add covariate data to output dataframe
    outputRow = covariateDataSheet.iloc[i, 0:4]
    outputRow.RasterFilePath = outputCovariatePath
    outputRow[dropResAggCol] = float('nan')

    outputCovariateSheet = pd.concat([outputCovariateSheet, pd.DataFrame([outputRow])])

#%% save updated covariate data to scenario 
myScenario.save_datasheet(name="CovariateData", data=outputCovariateSheet) 

# %%
