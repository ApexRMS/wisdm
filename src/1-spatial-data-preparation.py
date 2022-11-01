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

from typing_extensions import dataclass_transform
import pysyncrosim as ps     
import numpy as np          
import pandas as pd          
import os
import rioxarray
import rasterio
import psutil
import geopandas as gpd
import shapely

from shapely.geometry import Point #, shape
from rasterio.enums import Resampling #, MergeAlg

#%% Connect to SyncroSim library ------------------------------------------------

# Load current scenario
myScenario = ps.Scenario()  

# Create a temporary folder for storing rasters
ssimTempDir = ps.runtime_temp_folder("Data")

# Get path to scnario inputs 
ssimInputDir = myScenario.library.location + ".input\Scenario-" + str(myScenario.sid)

# Load datasheets
# inputs
covariatesSheet = myScenario.project.datasheets("Covariates")
covariateDataSheet = myScenario.datasheets("CovariateData")
fieldDataSheet = myScenario.datasheets("FieldData")
fieldDataOptions = myScenario.datasheets("FieldDataOptions")
templateRasterSheet = myScenario.datasheets("TemplateRaster")

# outputs
outputCovariateSheet = myScenario.datasheets("CovariateData", empty = True)

#%% Check inputs and set defaults ---------------------------------------------

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

#%% Load template raster ----------------------------------------------------------------

templatePath = ssimInputDir + "\\wisdm_TemplateRaster\\" + templateRasterSheet.RasterFilePath.item()
templateRaster = rioxarray.open_rasterio(templatePath, chunks=True, masked=True)
templateMask = rasterio.open(templatePath).read(1, masked=True).mask

# Get information about template
templateCRS = templateRaster.rio.crs
if templateCRS.is_valid == False:
    raise ValueError("Template has an invalid CRS (authority code). See {documention} for a list of accepted authority codes.")
templateResolution = templateRaster.rio.resolution()
templateExtent = list(templateRaster.rio.bounds())
templateTransform = templateRaster.rio.transform()
templatePixelSize = templateResolution[0]*-templateResolution[1]

if rasterio.dtypes.is_ndarray(templateMask):
    # Create template to polygon
    templatePolygons = []
    for shape, value in rasterio.features.shapes(templateMask.astype(np.int16), mask=np.invert(templateMask), transform=templateTransform):
        # shapely.geometry.shape(shape)
        templatePolygons.append(shapely.geometry.shape(shape))

    templatePolygons = shapely.geometry.MultiPolygon(templatePolygons)
    if not templatePolygons.is_valid:
        templatePolygons = templatePolygons.buffer(0)
        # Sometimes buffer() converts a simple Multipolygon to just a Polygon,
        # need to keep it a Multi throughout
        if templatePolygons.type == 'Polygon':
            templatePolygons = shapely.geometry.MultiPolygon([templatePolygons])
    # TO DO: save to geodataframe with mathcin transform etc as template....   
         
    templatePolygons = gpd.GeoDataFrame({'geometry':templatePolygons}, crs=templateCRS)
else:
    print('The template raster does not include a "No Data" mask.', 
    'All output covariate and site data will be clipped', 
    'to the full extent of the template raster.')
    templateMask = np.ones(templateRaster.shape, dtype=bool)
    templatePolygons = []


#%% Loop through and "PARC" covariate rasters -------------------------------------------

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
    
    #%% Add NoData mask to covariate raster
    covariateRaster = np.ma.masked_array(covariateRaster, templateMask)
    covariateRaster.fill_value = -9999
   
    #%% Set raster output parameters
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
        
    #%% Write processed covariate raster to file with internal NoData mask 
    with rasterio.Env(GDAL_TIFF_INTERNAL_MASK=True):
        with rasterio.open(outputCovariatePath, mode="w", **raster_params,
                            masked= True, tiled=True, windowed=True, overwrite=True) as src:
            src.write(covariateRaster)
    
    #%% Add covariate data to output dataframe
    outputRow = covariateDataSheet.iloc[i, 0:4]
    outputRow.RasterFilePath = outputCovariatePath
    outputRow[dropResAggCol] = float('nan')

    outputCovariateSheet = pd.concat([outputCovariateSheet, pd.DataFrame([outputRow])])

#%% Save updated covariate data to scenario 
myScenario.save_datasheet(name="CovariateData", data=outputCovariateSheet) 

#%% Prepare site data -------------------------------------------------------

nInitial = len(fieldDataSheet.SiteID)

# Create shapely points from the coordinate-tuple list
siteCoords = [Point(x, y) for x, y in zip(fieldDataSheet.X, fieldDataSheet.Y)]

# Define field data crs
if pd.isnull(fieldDataOptions.EPSG[0]):
    fieldDataCRS = templateCRS
else:
    fieldDataCRS = fieldDataOptions.EPSG[0]

# Convert shapely object to a geodataframe with a crs
sites = gpd.GeoDataFrame(fieldDataSheet, geometry=siteCoords, crs=fieldDataCRS)

# Reproject points if site crs differs from template crs
if sites.crs != templateCRS:
    sites = sites.to_crs(templateCRS)

#%% Clip sites to template extent
if rasterio.dtypes.is_ndarray(templatePolygons):
    sites = gpd.clip(sites,templatePolygons)
else:
    sites = gpd.clip(sites, templateExtent)
nFinal = len(sites.SiteID)

if nFinal < nInitial:
    print("\n", nInitial-nFinal, "sites out of", nInitial, 
        "total sites in the input field data were outside the template extent \nand were removed from the output.",
        nFinal, "sites were retained.\n")

#%% Create an index raster from template
baseArray = np.array(np.arange(0,templateRaster.rio.height*templateRaster.rio.width,1))
baseArray.shape = templateRaster.rio.height,templateRaster.rio.width 
band = baseArray
baseArray = band[np.newaxis,:,:]
indexRaster = templateRaster
indexRaster.values = baseArray

# Set raster output parameters
raster_params = {
    'driver': 'GTiff',
    'width': templateRaster.rio.width,
    'height': templateRaster.rio.height,
    'count': 1,
    'dtype': indexRaster.dtype, 
    'compress': 'lzw',
    'crs': templateCRS,
    'transform': templateTransform
    }
        
# Write raster to file
tempOutputPath = os.path.join(ssimTempDir, "indexRaster.tif")
with rasterio.open(tempOutputPath, mode="w", **raster_params,
                    masked= True, tiled=True, windowed=True, overwrite=True) as src:
    src.write(indexRaster)

#%% Extract raster ids for each point
with rasterio.open(tempOutputPath) as indexRast:
    rasterCellIDs = []
    for point in sites.geometry:
        x = point.xy[0][0]
        y = point.xy[1][0]
        row, col = indexRast.index(x,y)
        rasterCellIDs.append(indexRast.read(1)[row,col])
        
sites["RasterCellID"] = rasterCellIDs

#%% If there are multiple points per cell - Aggregate or Weight sites
if pd.notnull(fieldDataOptions.AggregateAndWeight[0]):
    if len(np.unique(sites.RasterCellID)) != len(sites.RasterCellID):
        # find duplicates
        seen = set()
        dupes = []
        for x in sites.RasterCellID.tolist():
            if x in seen:
                dupes.append(x)
            else:
                seen.add(x)
        dupes = np.unique(dupes).tolist()
        # if Aggregate sites is selected       
        if fieldDataOptions.AggregateAndWeight[0] == "Aggregate":
            # if presence absence data
            if all(sites.Response.isin([0,1])):
                for d in dupes:
                    sitesInd = sites.index[sites.RasterCellID == d].to_list() 
                    resp_d = sites.Response[sitesInd].to_list() 
                    if sum(resp_d) == 0 or np.mean(resp_d) == 1: # if all absence or all presence
                        sites.Response[sitesInd[1:]] = -9999 
                    else: # if response is mix of presence/absence
                        keep_d = (sites.Response[sitesInd] == 1).index[0] # keep a presence and convert rest of repeat sites to background points 
                        sitesInd.remove(keep_d)
                        sites.Response[sitesInd] = -9999 
            else: # if count data 
                for d in dupes:
                    sitesInd = sites.index[sites.RasterCellID == d].to_list() 
                    resp_d = sites.Response[sitesInd].to_list()
                    if sum(resp_d) == 0: # if all counts are zero
                        sites.Response[sitesInd[1:]] = -9999 
                    else: # if any counts are greater then zero
                        sites.Response[sitesInd[0]] = sum(resp_d)
                        sites.Response[sitesInd[1:]] = -9999 
        else: # if weight sites is selected
            if all(sites.Weight.isna()): # check for user defined weights;
                sites.Weight = 1
                for d in dupes:
                    sitesInd = sites.index[sites.RasterCellID == d].to_list()
                    weight_d = 1/len(sitesInd)
                    sites.Weights[sitesInd] = weight_d
            else: 
                print("\nWeights were already present in the field data, new weights were not assigned.\n")

    else: 
        print("\nOnly one field data observation present per pixel; no aggregation or weighting required.\n")

#%% Save updated field data to scenario 
outputFieldDataSheet = sites.iloc[:,0:7]
myScenario.save_datasheet(name="FieldData", data=outputFieldDataSheet) 

#%% Drop sites with repeat cell repeats  
dropInd = sites.index[sites.Response == -9999].tolist()
sites = sites.drop(dropInd)
sitesOut = sites[["SiteID", "RasterCellID"]]

#%% Extract covariate values for each site
# for i in range(len(covariateDataSheet.CovariatesID)):
#     # Load processed covariate rasters and extract site values
#     outputCovariatePath = os.path.join(ssimTempDir, covariateDataSheet.RasterFilePath[i])
#     with rasterio.open(outputCovariatePath) as r:
#         rastValues = []
#         for point in sites.geometry:
#             x = point.xy[0][0]
#             y = point.xy[1][0]
#             row, col = r.index(x,y)
#             rastValues.append(r.read(1)[row,col])
            
#     sites[covariateDataSheet.CovariatesID[i]] = rastValues
#
# siteData = sites.iloc[:,[0] + list(range(9,len(sites.columns)))]

#%% Create binary raster indicating raster cells with sites
binaryRaster = indexRaster.isin(rasterCellIDs).astype(int) 

# Write raster to file (for testing)
# tempOutputPath = os.path.join(ssimTempDir, "binarySitesRaster.tif")
# with rasterio.open(tempOutputPath, mode="w", **raster_params,
#                    masked= True, tiled=True, windowed=True, overwrite=True) as src:
#    src.write(binaryRaster)

#%% Extract covariate values for each cell with site data
siteData = binaryRaster*indexRaster
siteData = siteData.values.flatten()
siteData = pd.DataFrame(siteData, columns=["RasterCellID"])

for i in range(len(covariateDataSheet.CovariatesID)):
    # Load processed covariate rasters and extract site values
    outputCovariatePath = os.path.join(ssimTempDir, covariateDataSheet.RasterFilePath[i])
    covariateRaster = rioxarray.open_rasterio(outputCovariatePath, chunks=True)
    data_i = binaryRaster*covariateRaster
    data_i = data_i.values.flatten()
    siteData[covariateDataSheet.CovariatesID[i]] = data_i

#%% Merge and trim to Site data
siteData = pd.merge(sitesOut, siteData, on="RasterCellID")
siteData = siteData.drop(columns=["RasterCellID"])

#%% Convert site data to long format
siteData = pd.melt(siteData, id_vars= "SiteID", value_vars=siteData.columns[1:], var_name="CovariatesID", value_name="Value")

#%% Save site data to scenario 
myScenario.save_datasheet(name="SiteData", data=siteData) 


