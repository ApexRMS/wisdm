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
import dask  # noqa: E402
import shapely
from shapely.geometry.base import BaseGeometry
from shapely.geometry import Point, box  # noqa: E402
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
ps.environment.update_run_log(
    "is PROJ network enabled?", networkSheet.NetworkEnabled.item()
    )

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
# assuming template raster is single-band, open with rioxarray for CRS and extent info

# Get information about template

with rasterio.open(templatePath) as src:
    raster_crs = src.crs
    raster_epsg = raster_crs.to_epsg()

print("template CRS:", raster_crs)
print("template EPSG:", raster_epsg)

# # Get information about template
# templateCRS = CRS(templateRaster.rio.crs)
# try:
#     templateCRS.to_wkt()
# except ValueError:
#     raise ValueError(
#         "Template has an invalid CRS (authority code). See documentation for a list of accepted authority codes."
#     )

# update progress bar
ps.environment.progress_bar()

# NEW: Build a shapely box for the raster extent (minx, miny, maxx, maxy) — replaces templatePolygons
templateRaster = rioxarray.open_rasterio(templatePath, masked=True)
def raster_extent_box(rio_raster):
    return box(*rio_raster.rio.bounds())
templateExtent = raster_extent_box(templateRaster)

# Prepare site data ------------------------------------------------------------
nInitial = len(fieldDataSheet.SiteID)

# Create shapely points from the coordinate-tuple list
siteCoords = [Point(x, y) for x, y in zip(fieldDataSheet.X, fieldDataSheet.Y)]

# Define field data crs
fieldDataCRS = raster_crs if pd.isna(fieldDataOptions.EPSG[0]) else fieldDataOptions.EPSG[0]

# Convert shapely object to a geodataframe with a crs
sites = gpd.GeoDataFrame(fieldDataSheet, geometry=siteCoords, crs=fieldDataCRS)

# Reproject safely -------------------------------------------------------------
if sites.crs != raster_crs:
    sitesPRJ = sites.to_crs(raster_crs)
    del fieldDataSheet, siteCoords
    # Post-check for infinities (should be none)
    def has_inf_bounds(geom: BaseGeometry) -> bool:
        if geom is None or geom.is_empty:
            return True
        b = geom.bounds
        return any(np.isinf(v) or np.isnan(v) for v in b)
    bad_idx = sitesPRJ.index[sitesPRJ.geometry.apply(has_inf_bounds)]
    # if infinite geo, try again (GDAL connection)
    if len(bad_idx) > 0:
       sitesPRJ = sites.to_crs(raster_crs)
    else:
        ps.environment.update_run_log(
          "Reprojection successful: no inf/NaN bounds."
        )
    # if infinite geo again, report error
    bad_idx = sitesPRJ.index[sitesPRJ.geometry.apply(has_inf_bounds)]
    if len(bad_idx) > 0:
      raise ValueError(
            "Reprojection unsuccessful. Please check field (site) location data for: 1) coordinate values within known projection (CRS) range, 2) bad vertices (coordinates with NaN) "
      )
    # create copy
    sites = sitesPRJ.copy()

# Convert background sites to 0 for processing; restore before saving (mirrors other transformers)
bgSiteIds = sites.SiteID[sites.Response == backgroundValue].tolist()
sites.loc[sites.Response == backgroundValue, "Response"] = 0

# Clip sites to template extent using bounding box (no polygons needed)
sites = gpd.clip(sites, templateExtent)

# get final count of sites after clipping to template extent
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
xs = sites.geometry.x.values
ys = sites.geometry.y.values

with rasterio.open(templatePath) as src:
    rasterRows, rasterCols = src.index(xs, ys)
    data_mask = src.dataset_mask()  # 0 = NoData, non-zero = valid

sites["RasterRow"] = rasterRows
sites["RasterCol"] = rasterCols
sites["RasterCellID"] = list(zip(rasterRows, rasterCols))
del rasterRows, rasterCols

# NEW: Filter out points that fall in NoData cells using the raster's dataset mask
valid_flags = [(data_mask[row, col] != 0) for (row, col) in sites["RasterCellID"]]
sites = (
    sites
    .sort_values("SiteID")
    .loc[valid_flags]
)

# Optional run log note for removals due to NoData
nAfterMask = len(sites.SiteID)
removedByMask = nFinal - nAfterMask
if removedByMask > 0:
    ps.environment.update_run_log(
        removedByMask, " site(s) within the template extent fell in NoData cells and were removed."
    )

# update progress bar
ps.environment.progress_bar()

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
        outputCovariatePath)
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
        sitesOut[covariateDataSheet.CovariatesID[i]] == covariateRaster.rio.nodata, covariateDataSheet.CovariatesID[i]
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
# drop any duplicates
siteData.drop_duplicates(inplace=True)

# get freq of NoData values by covariate and report in log
def unique_site_nodata_by_cov(df):
    return (
    df[df["Value"].isna()]
    .groupby("CovariatesID")["SiteID"]
    .nunique()
    .reset_index(name="sites_with_nodata")
    )
freq_table = unique_site_nodata_by_cov(siteData)


# drop sites where covariate has NA values
siteData_filtered = siteData.groupby('SiteID').filter(lambda x: x['Value'].notna().all())
nInitial = siteData['SiteID'].nunique()
nFinal = siteData_filtered['SiteID'].nunique()
if nFinal < nInitial:
    ps.environment.update_run_log(
        nInitial - nFinal, " sites out of ", nInitial,
        " total sites in the input field data had NoData in 1 or more covariates and were removed. ",
        nFinal, " sites were retained. Please check covariate data for NoData values"
    )
    ps.environment.update_run_log(
        "Sites with No Data Values:", freq_table
    )


# Save site data to scenario
myScenario.save_datasheet(name="wisdm_SiteData", data=siteData_filtered)

# update progress bar
ps.environment.progress_bar(report_type="end")