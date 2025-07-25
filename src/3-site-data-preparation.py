## ---------------------------------
## wisdm - site data preparation
## ApexRMS, March 2024
## ---------------------------------

# built under Python version 3.11.0 & SyncroSim version 3.0.0
# Script pulls in covariate rasters and processes field data to ensure sites 
# are in the template CRS and extent; aggregates or weights sites by spatial distribution;
# extracts site-specific covariate data and generates datasheet of covariate site data.

# Source dependencies ----------------------------------------------------------

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
        filenames = [f for f in os.listdir(folder) if f.startswith("gdal") and f.endswith(".dll")]

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
import pandas as pd          
import rioxarray
import xarray
import geopandas as gpd
import shapely
import dask
import pyproj
# import spatialUtils

from dask.distributed import Client, Lock
from shapely.geometry import Point #, shape
from rasterio.enums import Resampling #, MergeAlg

pd.options.mode.chained_assignment = None

ps.environment.update_run_log('3 - Site Data Preparation => Begin')

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
       
# Connect to SyncroSim library ------------------------------------------------

# Load current scenario
myScenario = ps.Scenario()  

# Create a temporary folder for storing rasters
# ssimTempDir = myLibrary.info["Value"][myLibrary.info.Property == "Temporary files:"].item()
ssimTempDir = ps.runtime_temp_folder(os.path.join("DataTransfer", "Scenario-" + str(myScenario.sid)))

# Load datasheets
# inputs
networkSheet = myScenario.library.datasheets("wisdm_Network")
covariateDataSheet = myScenario.datasheets("wisdm_CovariateData", show_full_paths=True)
fieldDataSheet = myScenario.datasheets("wisdm_FieldData")
fieldDataOptions = myScenario.datasheets("wisdm_FieldDataOptions")
templateRasterSheet = myScenario.datasheets("wisdm_TemplateRaster", show_full_paths=True)
multiprocessingSheet = myScenario.datasheets("core_Multiprocessing")

# outputs
# outputCovariateSheet = myScenario.datasheets("CovariateData", empty = True)

# Set progress bar ---------------------------------------------------------

steps = 5 + len(covariateDataSheet.CovariatesID)
ps.environment.progress_bar(report_type = "begin", total_steps = steps)

# Set up dask client -------------------------------------------------------

if multiprocessingSheet.EnableMultiprocessing.item() == "Yes":
    num_threads = multiprocessingSheet.MaximumJobs.item()
else:
    num_threads = 1

# Note: Follow link in output to view progress
dask.config.set(**{'temporary-directory': os.path.join(ssimTempDir, 'dask-worker-space'),
                   'distributed.scheduler.worker-ttl': None},
                   scheduler='threads', serializers=['dask'], deserializers=['dask'])
# client = Client(threads_per_worker = num_threads, n_workers = 1, processes=False)
# client


# Check inputs and set defaults ---------------------------------------------

# Set PROJ network connection
if networkSheet.NetworkEnabled.item() == "No":
    pyproj.network.set_network_enabled(active=False)
    # os.environ["PROJ_NETWORK"] = "OFF"
    # pyproj.network.is_network_enabled()

# Check that a template raster was provided
if pd.isnull(templateRasterSheet.RasterFilePath.item()):
    raise ValueError("Template raster is missing.")

# check if field data was provided
if len(fieldDataSheet) == 0:
    # raise warning if field data was not provided
    raise ValueError("Field data was not provided. Please provide field data before continuing.")

# check that field data is between 0 and 1
if any(fieldDataSheet.Response)>1:
   # raise warning if field data includes counts
   raise ValueError("Field data contains counts in 'Response' column when occurrence data is expected. Please provide presence-(pseudo)absence data before continuing.")

# check if site ids were provided
if any(pd.isna(fieldDataSheet.SiteID)):
    fieldDataSheet.SiteID= range(1,len(fieldDataSheet)+1)
    
# check of field data options were provided
if len(fieldDataOptions) == 0:
    fieldDataOptions = pd.concat([fieldDataOptions, pd.DataFrame({"AggregateOrWeight":["None"]})], 
                                 ignore_index=True)

 
# Load template raster ----------------------------------------------------------------

# set defualt chunk dimensions
chunkDims = 4096 #1024

templatePath = templateRasterSheet.RasterFilePath.item()
templateRaster = rioxarray.open_rasterio(templatePath, chunks={'x': chunkDims, 'y': chunkDims}, masked=True)
templateMask = templateRaster.to_masked_array().mask[0]

# Get information about template
templateCRS = templateRaster.rio.crs
try: templateCRS.to_wkt()
except ValueError:
    raise ValueError("Template has an invalid CRS (authority code). See {documention} for a list of accepted authority codes.")
# templateResolution = templateRaster.rio.resolution()
templateExtent = list(templateRaster.rio.bounds())
templateTransform = templateRaster.rio.transform()
# templatePixelSize = templateResolution[0]*-templateResolution[1]

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
        if templatePolygons.geom_type == 'Polygon':
            templatePolygons = shapely.geometry.MultiPolygon([templatePolygons])
          
    templatePolygons = gpd.GeoDataFrame(geometry=[templatePolygons], crs=templateCRS)
else:
    ps.environment.update_run_log('The template raster does not include a "No Data" mask. ', 
    'All output covariate and site data will be clipped', 
    ' to the full extent of the template raster.')
    # templateMask = np.ones(templateRaster.shape, dtype=bool)
    templatePolygons = []

# update progress bar
ps.environment.progress_bar()
del templateMask

# Prepare site data -------------------------------------------------------

nInitial = len(fieldDataSheet.SiteID)

# Create shapely points from the coordinate-tuple list
siteCoords = [Point(x, y) for x, y in zip(fieldDataSheet.X, fieldDataSheet.Y)]
# gpd.GeoSeries(siteCoords).plot()
    
# Define field data crs
if pd.isnull(fieldDataOptions.EPSG[0]):
    fieldDataCRS = templateCRS
else:
    fieldDataCRS = fieldDataOptions.EPSG[0]

# Convert shapely object to a geodataframe with a crs
sites = gpd.GeoDataFrame(fieldDataSheet, geometry=siteCoords, crs=fieldDataCRS)
del fieldDataSheet, siteCoords

# Reproject points if site crs differs from template crs
if sites.crs != templateCRS:
    sites_reprojected = sites.to_crs(templateCRS)

    # Due to a bug in pyproj, the x and y coords may be infinite after the first reprojection.
    # If this is the case, then reproject again and it should be fixed
    # See bug report: https://github.com/geopandas/geopandas/issues/3433
    if np.isinf(sites_reprojected.geometry.x).any() and np.isinf(sites_reprojected.geometry.y).any():
        sites_reprojected = sites.to_crs(templateCRS)
        if np.isinf(sites_reprojected.geometry.x).all() and np.isinf(sites_reprojected.geometry.y).all():
            raise ValueError("Site coordinates are infinite after reprojection. Please check the input field data for errors.")
        if np.isinf(sites_reprojected.geometry.x).any() and np.isinf(sites_reprojected.geometry.y).any():
            if np.isinf(sites_reprojected.geometry.x).all() and np.isinf(sites_reprojected.geometry.y).all():
                msg = "All site coordinates are infinite after reprojection. "
                msg += "This could be due to restrictive firewall settings or an error in the input data. "
                msg += "\r\nPlease check the following input datasheets for errors: "
                msg += "\r\n - Field Data"
                msg += "\r\n - Field Data Options"
                msg += "\r\n - Template Raster"
                raise ValueError(msg)
            else:
                msg = "Some site coordinates are infinite after reprojection. "
                msg += "\r\nPlease check the following input datasheets for errors: "
                msg += "\r\n - Field Data"
                msg += "\r\n - Field Data Options"
                msg += "\r\n - Template Raster"
                msg += "\r\nThe following sites were not reprojected properly: "
                msg += str(sites_reprojected[sites_reprojected.geometry.x.isin([np.inf])].SiteID.tolist())
                raise ValueError(msg)
        
    sites = sites_reprojected

# Clip sites to template extent
if rasterio.dtypes.is_ndarray(templatePolygons):
    sites = gpd.clip(sites,templatePolygons)
else:
    sites = gpd.clip(sites, templateExtent)
nFinal = len(sites.SiteID)

if nFinal < nInitial:
    ps.environment.update_run_log(nInitial-nFinal, " sites out of ", nInitial, 
            " total sites in the input field data were outside the template extent and were removed from the output. ",
            nFinal, " sites were retained.")

# Update xy to match geometry
sites.X = sites.geometry.apply(lambda p: p.x)
sites.Y = sites.geometry.apply(lambda p: p.y)

# Extract raster ids for each point
rasterCellIDs = []
rasterRows = []
rasterCols = []
with rasterio.open(templatePath) as src:
    for point in sites.geometry:
        x = point.xy[0][0]
        y = point.xy[1][0]
        row, col = src.index(x,y)
        rasterCellIDs.append((row,col))
        rasterRows.append(row)
        rasterCols.append(col)

sites["RasterRow"] = rasterRows
sites["RasterCol"] = rasterCols
sites["RasterCellID"] = rasterCellIDs
    
# update progress bar
ps.environment.progress_bar()
del rasterRows, rasterCols, rasterCellIDs

# If there are multiple points per cell - Aggregate or Weight sites
if fieldDataOptions.AggregateAndWeight[0] != "None":
    if len(np.unique(sites.RasterCellID)) != len(sites.RasterCellID):
        # find duplicates
        seen = set()
        dupes = []
        for x in sites.RasterCellID.tolist():
            if x in seen:
                dupes.append(x)
            else:
                seen.add(x)
        dupes = list(set(dupes)) # get unsorted unique list of tuples

        # if Aggregate sites is selected       
        if fieldDataOptions.AggregateAndWeight[0] == "Aggregate":
            # if presence absence data
            if all(sites.Response.isin([0,1])):
                for d in dupes:
                    sitesInd = sites.index[sites.RasterCellID == d].to_list() 
                    resp_d = sites.Response[sitesInd].to_list() 
                    if sum(resp_d) == 0 or np.mean(resp_d) == 1: # if all absence or all presence
                        sites.loc[sitesInd[1:], "Response"] = -9999 
                    else: # if response is mix of presence/absence
                        keep_d = (sites.Response[sitesInd] == 1).index[0] # keep a presence and convert rest of repeat sites to background points 
                        sitesInd.remove(keep_d)
                        sites.loc[sitesInd, "Response"] = -9999
            else: # if count data 
                for d in dupes:
                    sitesInd = sites.index[sites.RasterCellID == d].to_list() 
                    resp_d = sites.Response[sitesInd].to_list()
                    if sum(resp_d) == 0: # if all counts are zero
                        sites.loc[sitesInd[1:], "Response"] = -9999
                    else: # if any counts are greater then zero
                        sites.loc[sitesInd[0], "Response"] = sum(resp_d)
                        sites.loc[sitesInd[1:], "Response"] = -9999
        else: # if weight sites is selected
            if all(sites.Weight.isna()): # check for user defined weights;
                sites.Weight = 1
                for d in dupes:
                    sitesInd = sites.index[sites.RasterCellID == d].to_list()
                    weight_d = np.float64(1/len(sitesInd))
                    sites = sites.astype({"Weight": float})
                    sites.loc[sitesInd, "Weight"] = weight_d
            else: 
                ps.environment.update_run_log("Weights were already present in the field data, new weights were not assigned.")

    else: 
        ps.environment.update_run_log("Only one field data observation present per pixel; no aggregation or weighting required.")

# Save updated field data to scenario 
outputFieldDataSheet = sites.iloc[:,0:7]
outputFieldDataSheet.SiteID = outputFieldDataSheet.SiteID.astype(int) # ensure SiteID is type integer
myScenario.save_datasheet(name="wisdm_FieldData", data=outputFieldDataSheet) 
  
# update progress bar
ps.environment.progress_bar()
del outputFieldDataSheet

# Drop sites with repeat cell repeats  
dropInd = sites.index[sites.Response == -9999].tolist()
sites = sites.drop(dropInd)

# Write sites to file (for testing)
# tempOutputPath = os.path.join(ssimTempDir, "sites.shp")
# sites.to_file(tempOutputPath)

# Create index arrays (note in xarray x=col and y=row from geodataframe)
yLoc = xarray.DataArray(sites.RasterRow, dims =["loc"])
xLoc = xarray.DataArray(sites.RasterCol, dims =["loc"])
sitesOut = sites[["SiteID"]] #, "RasterCellID"

# Extract covariate values for each site
for i in range(len(covariateDataSheet.CovariatesID)):
    # Load processed covariate rasters and extract site values
    outputCovariatePath = covariateDataSheet.RasterFilePath[i]
    covariateRaster = rioxarray.open_rasterio(outputCovariatePath, chunks=True)

    if covariateRaster.rio.width < templateRaster.rio.width or \
        covariateRaster.rio.height < templateRaster.rio.height:
        msg = f"Covariate raster dimensions are smaller than the template raster dimensions."
        msg += f"\nPlease ensure that all covariate rasters are at least the same size as the template raster."
        msg += f"\nCovariate name: {covariateDataSheet.CovariatesID[i]}."
        msg += f"\nCovariate raster dimensions: {covariateRaster.rio.width} x {covariateRaster.rio.height}. "
        raise ValueError(msg)

    sitesOut.loc[:,covariateDataSheet.CovariatesID[i]] = covariateRaster[0].isel(x=xLoc,y=yLoc).values.tolist() 

    # Replace no data values with -9999
    sitesOut.loc[sitesOut[covariateDataSheet.CovariatesID[i]] == covariateRaster.rio.nodata, covariateDataSheet.CovariatesID[i]] = None #-9999

    # update progress bar
    ps.environment.progress_bar()

# Convert site data to long format
siteData = pd.melt(sitesOut, id_vars= "SiteID", value_vars=sitesOut.columns[1:], var_name="CovariatesID", value_name="Value")
siteData.drop_duplicates(inplace=True)

# Save site data to scenario 
myScenario.save_datasheet(name="wisdm_SiteData", data=siteData)  

# update progress bar
ps.environment.progress_bar(report_type = "end")

