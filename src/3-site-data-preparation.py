# ---------------------------------
# wisdm - site data preparation
# ApexRMS, March 2024
# ---------------------------------

# built under Python version 3.11.0 & SyncroSim version 3.0.0
# Script pulls in covariate rasters and processes field data to ensure sites
# are in the template CRS and extent; aggregates or weights sites by spatial distribution;
# extracts site-specific covariate data and generates datasheet of covariate site data.

# Source dependencies ----------------------------------------------------------
# IMPORTANT: setup_functions must be imported before any non-stdlib packages.
# It removes user site-packages on import and setupCondaEnv() configures
# conda DLL paths. Linters must not reorder these imports.  # noqa: E402
import os  # noqa: E402
import sys  # noqa: E402
from setup_functions import (  # noqa: E402
    setupCondaEnv, checkGdalVersion, setupGdalProj,
    nodataValue, backgroundValue, defaultChunkDims, getNumThreads, signedDtype
)

setupCondaEnv()
checkGdalVersion()

# dependencies
import rasterio  # noqa: E402
import pysyncrosim as ps  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import rioxarray  # noqa: E402
import xarray  # noqa: E402
import geopandas as gpd  # noqa: E402
from shapely.geometry import Point, box  # noqa: E402
import dask  # noqa: E402
import pyproj  # noqa: E402
from rasterio.enums import Resampling  # noqa: E402

pd.options.mode.chained_assignment = None

ps.environment.update_run_log('3 - Site Data Preparation => Begin')

# Modify the os PROJ path (when running with Conda) ----
myLibrary = ps.Library()
setupGdalProj(myLibrary)

# Connect to SyncroSim library ------------------------------------------------
# Load current scenario
myScenario = ps.Scenario()

# Create a temporary folder for storing rasters
ssimTempDir = ps.runtime_temp_folder(os.path.join(
    "DataTransfer", "Scenario-" + str(myScenario.sid))
)

# Load datasheets
# inputs
networkSheet = myScenario.library.datasheets("wisdm_Network")
covariateDataSheet = myScenario.datasheets(
    "wisdm_CovariateData", show_full_paths=True
)
fieldDataSheet = myScenario.datasheets("wisdm_FieldData")
fieldDataOptions = myScenario.datasheets("wisdm_FieldDataOptions")
templateRasterSheet = myScenario.datasheets(
    "wisdm_TemplateRaster", show_full_paths=True
)
multiprocessingSheet = myScenario.datasheets("core_Multiprocessing")

# outputs
# outputCovariateSheet = myScenario.datasheets("CovariateData", empty = True)

# Set progress bar ---------------------------------------------------------
steps = 5 + len(covariateDataSheet.CovariatesID)
ps.environment.progress_bar(report_type="begin", total_steps=steps)

# Set up dask configuration -------------------------------------------------------
num_threads = getNumThreads(multiprocessingSheet)
dask.config.set(
    **{
        'temporary-directory': os.path.join(ssimTempDir, 'dask-worker-space'),
        'distributed.scheduler.worker-ttl': None
    },
    scheduler='threads',
    serializers=['dask'],
    deserializers=['dask']
)

# Check inputs and set defaults ---------------------------------------------
# Set PROJ network connection
if networkSheet.NetworkEnabled.item() == "No":
    pyproj.network.set_network_enabled(active=False)

# Check that a template raster was provided
if templateRasterSheet.empty or pd.isnull(templateRasterSheet.RasterFilePath.iloc[0]):
    raise ValueError("Template raster is missing.")

# check if field data was provided
if len(fieldDataSheet) == 0:
    raise ValueError(
        "Field data was not provided. Please provide field data before continuing."
    )

# check that response data is provided and between 0 and 1
if pd.isnull(fieldDataSheet.Response).any():
    raise ValueError(
        "Field data is missing values in the 'Response' column. Please provide presence-(pseudo)absence data before continuing."
    )
if (fieldDataSheet.Response > 1).any():
    raise ValueError(
        "Field data contains counts in 'Response' column when occurrence data is expected. Please provide presence-(pseudo)absence data before continuing."
    )

# check if site ids were provided
if any(pd.isna(fieldDataSheet.SiteID)):
    fieldDataSheet.SiteID = range(1, len(fieldDataSheet) + 1)

# check if field data options were provided
if len(fieldDataOptions) == 0:
    fieldDataOptions = pd.concat(
        [fieldDataOptions, pd.DataFrame({"AggregateAndWeight": ["None"]})],
        ignore_index=True
    )

# Load template raster ----------------------------------------------------------------
templatePath = templateRasterSheet.RasterFilePath.item()
templateRaster = rioxarray.open_rasterio(
    templatePath, chunks=defaultChunkDims, masked=True
)

# Get information about template
templateCRS = templateRaster.rio.crs
try:
    templateCRS.to_wkt()
except ValueError:
    raise ValueError(
        "Template has an invalid CRS (authority code). See documentation for a list of accepted authority codes."
    )

templateExtent = list(templateRaster.rio.bounds())
templateTransform = templateRaster.rio.transform()

# update progress bar
ps.environment.progress_bar()

# NEW: Build a shapely box for the raster extent (minx, miny, maxx, maxy) — replaces templatePolygons
extent_geom = box(*templateExtent)

# Prepare site data -------------------------------------------------------
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
del fieldDataSheet, siteCoords

# Convert background sites to 0 for processing; restore before saving (mirrors other transformers)
bgSiteIds = sites.SiteID[sites.Response == backgroundValue].tolist()
sites.loc[sites.Response == backgroundValue, "Response"] = 0

# Reproject points if site crs differs from template crs
if sites.crs != templateCRS:
    sites = sites.to_crs(templateCRS)

# Clip sites to template extent using bounding box (no polygons needed)
sites = gpd.clip(sites, extent_geom)
nFinal = len(sites.SiteID)

if nFinal < nInitial:
    ps.environment.update_run_log(
        nInitial - nFinal, " sites out of ", nInitial,
        " total sites in the input field data were outside the template extent and were removed from the output. ",
        nFinal, " sites were retained."
    )

# Update xy to match geometry
sites.X = sites.geometry.apply(lambda p: p.x)
sites.Y = sites.geometry.apply(lambda p: p.y)

# Extract raster ids for each point
rasterCellIDs = []
rasterRows = []
rasterCols = []
with rasterio.open(templatePath) as src:
    for point in sites.geometry:
        row, col = src.index(point.x, point.y)
        rasterCellIDs.append((row, col))
        rasterRows.append(row)
        rasterCols.append(col)

sites["RasterRow"] = rasterRows
sites["RasterCol"] = rasterCols
sites["RasterCellID"] = rasterCellIDs

# NEW: Filter out points that fall in NoData cells using the raster's dataset mask
with rasterio.open(templatePath) as src:
    data_mask = src.dataset_mask()  # 0 = NoData, non-zero = valid
valid_flags = [(data_mask[row, col] != 0) for (row, col) in sites["RasterCellID"]]
sites = sites.loc[valid_flags].copy()

# Optional run log note for removals due to NoData
nAfterMask = len(sites.SiteID)
removedByMask = nFinal - nAfterMask
if removedByMask > 0:
    ps.environment.update_run_log(
        removedByMask, " site(s) within the template extent fell in NoData cells and were removed."
    )

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
        dupes = list(set(dupes))  # get unsorted unique list of tuples

        # if Aggregate sites is selected
        if fieldDataOptions.AggregateAndWeight[0] == "Aggregate":
            # if presence absence data
            if all(sites.Response.isin([0, 1])):
                for d in dupes:
                    sitesInd = sites.index[sites.RasterCellID == d].to_list()
                    resp_d = sites.Response[sitesInd].to_list()
                    if sum(resp_d) == 0 or np.mean(resp_d) == 1:  # all absence or all presence
                        sites.loc[sitesInd[1:], "Response"] = nodataValue
                    else:  # mix of presence/absence
                        keep_d = sitesInd[resp_d.index(1)]
                        sitesInd.remove(keep_d)
                        sites.loc[sitesInd, "Response"] = nodataValue
            else:  # count data
                for d in dupes:
                    sitesInd = sites.index[sites.RasterCellID == d].to_list()
                    resp_d = sites.Response[sitesInd].to_list()
                    if sum(resp_d) == 0:  # all counts zero
                        sites.loc[sitesInd[1:], "Response"] = nodataValue
                    else:  # any counts > 0
                        sites.loc[sitesInd[0], "Response"] = sum(resp_d)
                        sites.loc[sitesInd[1:], "Response"] = nodataValue
        else:  # Weight sites
            if all(sites.Weight.isna()):  # user-defined weights absent
                sites.Weight = 1
                for d in dupes:
                    sitesInd = sites.index[sites.RasterCellID == d].to_list()
                    weight_d = np.float64(1 / len(sitesInd))
                    sites = sites.astype({"Weight": float})
                    sites.loc[sitesInd, "Weight"] = weight_d
            else:
                ps.environment.update_run_log(
                    "Weights were already present in the field data, new weights were not assigned."
                )
    else:
        ps.environment.update_run_log(
            "Only one field data observation present per pixel; no aggregation or weighting required."
        )

# Restore background sites to backgroundValue (only surviving sites remain as 0)
sites.loc[(sites.SiteID.isin(bgSiteIds)) & (sites.Response == 0), "Response"] = backgroundValue

# Save updated field data to scenario
outputFieldDataSheet = sites.iloc[:, 0:7]
outputFieldDataSheet.SiteID = outputFieldDataSheet.SiteID.astype(int)  # ensure SiteID is type integer
myScenario.save_datasheet(name="wisdm_FieldData", data=outputFieldDataSheet)

# update progress bar
ps.environment.progress_bar()
del outputFieldDataSheet

# Drop sites with repeat cell repeats
dropInd = sites.index[sites.Response == nodataValue].tolist()
sites = sites.drop(dropInd)

# Create index arrays (note in xarray x=col and y=row from geodataframe)
yLoc = xarray.DataArray(sites.RasterRow, dims=["loc"])
xLoc = xarray.DataArray(sites.RasterCol, dims=["loc"])
sitesOut = sites[["SiteID"]]  # , "RasterCellID"

# Extract covariate values for each site
for i in range(len(covariateDataSheet.CovariatesID)):
    # Load processed covariate rasters and extract site values
    outputCovariatePath = covariateDataSheet.RasterFilePath[i]
    covariateRaster = rioxarray.open_rasterio(
        outputCovariatePath, chunks=defaultChunkDims
    )

    # Ensure signed dtype (nodata value -9999)
    if covariateRaster.dtype != signedDtype(covariateRaster.dtype):
        raise ValueError(
            f"Covariate '{covariateDataSheet.CovariatesID[i]}' has an unsigned integer "
            f"dtype ({covariateRaster.dtype}), which is incompatible with the nodata "
            f"value -9999. Run Transformer 2 (Spatial Data Preparation) first to convert "
            f"covariate rasters to a signed type before continuing."
        )

    # Ensure coverage at least as large as template
    if covariateRaster.rio.width < templateRaster.rio.width or \
            covariateRaster.rio.height < templateRaster.rio.height:
        msg = f"Covariate raster dimensions are smaller than the template raster dimensions."
        msg += f"\nPlease ensure that all covariate rasters are at least the same size as the template raster."
        msg += f"\nCovariate name: {covariateDataSheet.CovariatesID[i]}."
        msg += f"\nCovariate raster dimensions: {covariateRaster.rio.width} x {covariateRaster.rio.height}. "
        raise ValueError(msg)

    # Extract values by row/col indices
    sitesOut.loc[:, covariateDataSheet.CovariatesID[i]] = (
        covariateRaster[0].isel(x=xLoc, y=yLoc).values.tolist()
    )

    # Replace no data values with None (to match SyncroSim NA behavior)
    sitesOut.loc[
        sitesOut[covariateDataSheet.CovariatesID[i]] == covariateRaster.rio.nodata,
        covariateDataSheet.CovariatesID[i]
    ] = None

    # update progress bar
    ps.environment.progress_bar()

# Convert site data to long format
siteData = pd.melt(
    sitesOut,
    id_vars="SiteID",
    value_vars=sitesOut.columns[1:],
    var_name="CovariatesID",
    value_name="Value"
)
siteData.drop_duplicates(inplace=True)

# Save site data to scenario
myScenario.save_datasheet(name="wisdm_SiteData", data=siteData)

# update progress bar
ps.environment.progress_bar(report_type="end")