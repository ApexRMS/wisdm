#%% Load packages
import pysyncrosim as ps     
import numpy as np          
import pandas as pd          
import os
import rioxarray
import rasterio

from rasterio.enums import Resampling

#%% Get the SyncroSim Scenario that is currently running
myScenario = ps.Scenario()  

#%% Create a temporary folder for storing rasters
ssimTempDir = ps.runtime_temp_folder("TempOutputs")

# Get path to scnarios inputs 
ssimInputDir = myScenario.library.location + ".input\Scenario-" + str(myScenario.sid)

#%% load datasheets
# inputs
covariatesSheet = myScenario.project.datasheets("Covariates")
covariateDataSheet = myScenario.datasheets("CovariateData")
templateRasterSheet = myScenario.datasheets("TemplateRaster")
multiprocessingSheet = myScenario.datasheets("core_Multiprocessing")
spatialMulitprocessingSheet = myScenario.datasheets("corestime_Multiprocessing")

# outputs
outputCovariateSheet = myScenario.datasheets("CovariateData", empty = True)

#%% Check inputs and set defaults

# check that a template raster was provided
if templateRasterSheet.RasterFilePath.isnull():
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
if covariateDataSheet.AggregationMethod.isnull().values.any():
    index1 = covariateDataSheet.AggregationMethod[covariateDataSheet.AggregationMethod.isnull()].index.tolist()
    covariateDataSheet.AggregationMethod[index1] = "Mean"
    index2 = covariateDataSheet.CovariatesID[covariateDataSheet.CovariatesID.isin(catCovs)].index.tolist()
    index3 = np.intersect1d(index1,index2)
    covariateDataSheet.AggregationMethod[index3] = "Majority"


#%% define rio up/down-sample method
resampleAggregateMethodName = ["Nearest Neighbor", "Bilinear", "Cubic", "Cubic Spline", "Lanczos", "Mean", "Min", "Max", "Majority"]
resampleAggregateMethodCodes = ["nearest", "bilinear", "cubic", "cubic_spline", "lanczos", "average", "min", "max", "mode"]
covariateDataSheet["rio_resample"] = covariateDataSheet.ResampleMethod.replace(resampleAggregateMethodName, resampleAggregateMethodCodes)
covariateDataSheet["rio_aggregate"] = covariateDataSheet.AggregationMethod.replace(resampleAggregateMethodName, resampleAggregateMethodCodes)

#%% Load template raster 
templatePath = ssimInputDir + "\\wisdm_TemplateRaster\\" + templateRasterSheet.RasterFilePath.item()
templateRaster = rioxarray.open_rasterio(templatePath, chunks=True)
# templateDS = gdal.Open(templatePath)

#%% get information about template
templateCRS = templateRaster.rio.crs
templateResolution = templateRaster.rio.resolution()
templateExtent = templateRaster.rio.bounds()
templateTransform = templateRaster.rio.transform()
templatePixelSize = templateResolution[0]*-templateResolution[1]

#%% Loop through and "PARC" covaraite rasters  
for i in range(len(covariateDataSheet.CovariatesID)):

    #%% Load and process covariate rasters
    covariatePath = ssimInputDir + "\\wisdm_CovariateData\\" + covariateDataSheet.RasterFilePath[i]
    covariateRaster = rioxarray.open_rasterio(covariatePath, chunks=True)
    outputCovariatePath = os.path.join(ssimTempDir, os.path.basename(covariatePath))

    #%% Determine if raster needs to be resampled or aggregated
    
    # convert covariate resolution to template units (following code is faster then reprojecting for large rasters)
    covAffine = rasterio.warp.calculate_default_transform(covariateRaster.rio.crs, templateCRS, covariateRaster.shape[0], covariateRaster.shape[1], *covariateRaster.rio.bounds())[0]
    covPixelSize = covAffine[0]*-covAffine[4]

    # if covariate resolution is finer then template use aggregate method (if coarser use resample method) 
    if covPixelSize < templatePixelSize:
        rioResampleMethod = covariateDataSheet.rio_aggregate[i]
        dropResAggCol = "ResampleMethod"
    else:
        rioResampleMethod = covariateDataSheet.rio_resample[i] 
        dropResAggCol = "AggregationMethod"

    #%% reproject/resample/clip covariate raster to match template 
    covariateRaster = covariateRaster.rio.reproject_match(templateRaster,
                                                            resampling=Resampling[rioResampleMethod])

    # %% Save outputs
    # save prepped covaraite raster to temp file
    covariateRaster.rio.write_nodata(-9999, inplace=True)
    covariateRaster.rio.to_raster(outputCovariatePath, tiled=True, windowed=True, overwrite=True, compress="lzw")
    
    #%% add covariate data to output dataframe
    outputRow = covariateDataSheet.iloc[i, 0:4]
    outputRow.RasterFilePath = outputCovariatePath
    outputRow[dropResAggCol] = float('nan')

    outputCovariateSheet = outputCovariateSheet.append(outputRow)

#%% save updated covariate data to scenario 
myScenario.save_datasheet(name="CovariateData", data=outputCovariateSheet) 