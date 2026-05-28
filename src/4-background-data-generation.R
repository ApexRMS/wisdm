## -------------------------
## wisdm - background data generation
## ApexRMS, June 2025
## -------------------------

# built under R version 4.1.3, SyncroSim 3.1.10 & rsyncrosim 2.1.3
# Script generates background data (pseudo-absences) based on the defined
#  background data options and outputs updated field data and site data sheets

# source dependencies ----------------------------------------------------------

library(rsyncrosim) # install.packages("C:/GitHub/rsyncrosim", type="source", repos=NULL)
library(terra)
library(tidyr)
library(dplyr)
library(pander)

packageDir <- (Sys.getenv("ssim_package_directory"))
source(file.path(packageDir, "00-helper-functions.R"))
source(file.path(packageDir, "04-background-data-functions.R"))

updateRunLog('4 - Background Data Generation => Begin')

# Connect to library -----------------------------------------------------------

# Active scenario
myScenario <- scenario()

# temp directory
ssimTempDir <- ssimEnvironment()$TransferDirectory

# Read in datasheets
templateSheet <- datasheet(myScenario, "wisdm_TemplateRaster")
covariateDataSheet <- datasheet(
  myScenario,
  "wisdm_CovariateData",
  optional = T,
  lookupsAsFactors = F
)
fieldDataSheet <- datasheet(myScenario, "wisdm_FieldData", optional = T)
backgroundDataOptionsSheet <- datasheet(
  myScenario,
  "wisdm_BackgroundDataOptions",
  optional = T
)
siteDataSheet <- datasheet(
  myScenario,
  "wisdm_SiteData",
  optional = T,
  lookupsAsFactors = F
)

# Set progress bar -------------------------------------------------------------

steps <- 5 + length(covariateDataSheet$CovariatesID)
progressBar(type = "begin", totalSteps = steps)

# Generate and save random seed if not already set -----------------------------

if (nrow(validationDataSheet) < 1 || is.na(validationDataSheet$RandomSeed)) {
  if (nrow(validationDataSheet) < 1) {
    validationDataSheet <- safe_rbind(
      validationDataSheet,
      data.frame(RandomSeed = sample.int(.Machine$integer.max, 1))
    )
  } else {
    validationDataSheet$RandomSeed <- sample.int(.Machine$integer.max, 1)
  }
  saveDatasheet(myScenario, validationDataSheet, "wisdm_ValidationOptions")
}
set.seed(validationDataSheet$RandomSeed)

# Prep inputs ------------------------------------------------------------------

# drop no data (-9999) sites that resulted from spatial aggregation
fieldDataSheet <- fieldDataSheet[fieldDataSheet$Response != nodataValue, ]

#  Set defaults ----------------------------------------------------------------

## Background data options sheet
if (nrow(backgroundDataOptionsSheet) < 1) {
  backgroundDataOptionsSheet <- safe_rbind(
    backgroundDataOptionsSheet,
    data.frame(GenerateBackgroundSites = FALSE)
  )
}
if (is.na(backgroundDataOptionsSheet$GenerateBackgroundSites)) {
  backgroundDataOptionsSheet$GenerateBackgroundSites <- FALSE
}
if (backgroundDataOptionsSheet$GenerateBackgroundSites) {
  if (is.na(backgroundDataOptionsSheet$BackgroundSiteCount)) {
    backgroundDataOptionsSheet$BackgroundSiteCount <- sum(
      fieldDataSheet$Response
    )
  }
  if (is.na(backgroundDataOptionsSheet$BackgroundGenerationMethod)) {
    backgroundDataOptionsSheet$BackgroundGenerationMethod <- "Kernel Density Estimate (KDE)"
  }
  if (is.na(backgroundDataOptionsSheet$KDESurface)) {
    if (
      backgroundDataOptionsSheet$BackgroundGenerationMethod ==
      "Kernel Density Estimate (KDE)"
    ) {
      backgroundDataOptionsSheet$KDESurface <- "Continuous"
    }
  }
  if (is.na(backgroundDataOptionsSheet$Isopleth)) {
    if (
      backgroundDataOptionsSheet$KDESurface == "Binary" |
      backgroundDataOptionsSheet$BackgroundGenerationMethod ==
      "Minimum Convex Polygon (MCP)"
    ) {
      backgroundDataOptionsSheet$Isopleth <- 95
    }
  }
}

saveDatasheet(
  myScenario,
  backgroundDataOptionsSheet,
  "wisdm_BackgroundDataOptions"
)

# Error handling ---------------------------------------------------------------

if (nrow(fieldDataSheet) == 0L) {
  stop(
    "No Field Data found; please ensure that the Field Data datasheet is populated before continuing."
  )
}
if (backgroundDataOptionsSheet$GenerateBackgroundSites) {
  if (nrow(templateSheet) == 0L || is.na(templateSheet$RasterFilePath)) {
    stop(
      "No Template Raster found; please ensure that the Template Raster datasheet is populated before generating background sites."
    )
  }
  if (nrow(covariateDataSheet) == 0L) {
    stop(
      "No Covariate Data found; please ensure that the Covariate Data datasheet is populated before generating background sites."
    )
  }
}

# Generate pseudo-absences (if applicable) -------------------------------------

if (
  backgroundDataOptionsSheet$GenerateBackgroundSites &&
  any(fieldDataSheet$Response == backgroundValue)
) {
  updateRunLog(paste0(
    "\nWarning: ",
    sum(fieldDataSheet$Response == backgroundValue),
    " existing background site(s) found in ",
    "Field Data. These will be removed and regenerated based on the current Background Data Options. ",
    "To use the provided background data instead, disable Generate Background Sites.\n"
  ))
  fieldDataSheet <- fieldDataSheet[fieldDataSheet$Response != backgroundValue, ]
}

if (backgroundDataOptionsSheet$GenerateBackgroundSites) {
  if (
    backgroundDataOptionsSheet$BackgroundGenerationMethod ==
    "Kernel Density Estimate (KDE)"
  ) {
    methodInputs <- list(
      "method" = "kde",
      "surface" = backgroundDataOptionsSheet$KDESurface,
      "isopleth" = backgroundDataOptionsSheet$Isopleth
    )
  }
  if (
    backgroundDataOptionsSheet$BackgroundGenerationMethod ==
    "Minimum Convex Polygon (MCP)"
  ) {
    methodInputs <- list(
      "method" = "mcp",
      "surface" = NA,
      "isopleth" = backgroundDataOptionsSheet$Isopleth
    )
  }
  
  ### load template raster
  templateRaster <- rast(templateSheet$RasterFilePath)
  
  backgroundSurfacePointGeneration(
    sp = "species",
    n=backgroundDataOptionsSheet$BackgroundSiteCount+100,
    template = templateRaster,
    outputDir = ssimTempDir,
    dat = fieldDataSheet,
    method = methodInputs)
  
  progressBar()
  gc()
  

  # add background point to field data
  bgData <- data.table::fread(file.path(
    ssimTempDir,
    paste0("species_", methodInputs$method, "_bg_pts.csv")
  ))


  ### stop if no background pts were returned
  if (nrow(bgData) == 0) {
    updateRunLog(
      "\nWarning: No background sites remained after filtering against presence pixels. No background data was added to Field Data.\n"
    )
    stop("Warning: No background sites remained after filtering against presence pixels. No background data was added to Field Data.")
  }  
  
  progressBar()
  
  ### create SiteIDs for background data
  startId <- max(fieldDataSheet$SiteID, siteDataSheet$SiteID) + 1
  bgData$SiteID <- startId:(startId + nrow(bgData) - 1)

  ### create spatvect object for extract
  vbgPts <- vect(bgData, geom = c("X", "Y"), crs = crs(templateRaster))

  ### add columns for combining with field data. 
  bgData$UseInModelEvaluation <- NA
  bgData$ModelSelectionSplit <- NA
  bgData$Weight <- NA

  ### update field data weights. This is a legacy item, and needs to be revisted. Currently set the weights for ALL background data to 1 if any weights were supplied with field data. User options for setting other weights.
  if (any(!is.na(fieldDataSheet$Weight))) {
    bgData$Weight <- 1
  }

  ### subset cols
  sel_cols = names(fieldDataSheet)
  bgData = bgData[,..sel_cols]; rm(sel_cols)
  fieldDataSheet <- rbind(fieldDataSheet, bgData)
  fieldDataSheet$SiteID <- format(fieldDataSheet$SiteID, scientific = F)

  ## Extract covariate data for background sites  -----
  ### loop through each 
  for (i in 1:nrow(covariateDataSheet)) {
    ri = terra::rast(covariateDataSheet$RasterFilePath[i])
    vals = terra::extract(
      ri, vbgPts,
      method="simple",
      cells=FALSE,
      xy=FALSE,
      ID=TRUE)
    rm(ri); gc()
    bgData[,covariateDataSheet$CovariatesID[i]:=vals[,2]]
    # progressBar()
  }
  bgData[, SiteID := as.character(SiteID)]

  CovariatesID = covariateDataSheet$CovariatesID[!covariateDataSheet$CovariatesID %in% "mtpi_10m"] 

  ### trim cols to format and melt to long format
  keep_cols = c("SiteID", covariateDataSheet$CovariatesID)
  bgData = bgData[,..keep_cols]
  bgSiteData <- gather(
    data = bgData,
    key = CovariatesID,
    value = Value,
    -SiteID
  )

  siteDataSheet <- rbind(siteDataSheet, bgSiteData)
  siteDataSheet$SiteID <- format(siteDataSheet$SiteID, scientific = F)
  siteDataSheet <- siteDataSheet %>%
    distinct(SiteID, CovariatesID, .keep_all = T)
  rm(bgSiteData)
  gc()

  # save site data to scenario
  saveDatasheet(myScenario, siteDataSheet, "wisdm_SiteData")
  rm(siteDataSheet)
  gc()

  # save updated field data to scenario
  saveDatasheet(myScenario, fieldDataSheet, "wisdm_FieldData", append = F)
}

progressBar(type = "end")
