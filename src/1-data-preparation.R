## wisdm - data preparation
## ApexRMS, March 2022

# R version 4.1.1
# Script pulls in field data and covariate data and generates dataframe of site-specific covaraite data

# source dependencies ----------------------------------------------------------

packageDir <- (Sys.getenv("ssim_package_directory"))
source(file.path(packageDir, "0-dependencies.R"))
source(file.path(packageDir, "01-data-prep-functions.R"))

# Connect to library -----------------------------------------------------------

# Active project and scenario
myProject <- rsyncrosim::project()
myScenario <- scenario()

# Read in datasheets
covariatesSheet <- datasheet(myProject, "wisdm_Covariates", optional = T)
fieldDataSheet <- datasheet(myScenario, "wisdm_FieldData", optional = T)
fieldDataOptions <- datasheet(myScenario, "wisdm_FieldDataOptions")
covariateDataSheet <- datasheet(myScenario, "wisdm_CovariateData", lookupsAsFactors = F)
validationDataSheet <- datasheet(myScenario, "wisdm_ValidationOptions")
templateRasterSheet <- datasheet(myScenario, "wisdm_TemplateRaster")
# siteDataSheet <- datasheet(myScenario, "wisdm_SiteData", optional = T) # lookupsAsFactors = F,

#  Set defaults ----------------------------------------------------------------  

## field data options
if(nrow(fieldDataOptions)<1){
  fieldDataOptions <- addRow(fieldDataOptions, list(EPSG = NA,
                                                    AggregateAndWeight = FALSE))
  }

## Validation Sheet
if(nrow(validationDataSheet)<1){
  validationDataSheet <- addRow(validationDataSheet, list(SplitData = FALSE,
                                                          CrossValidate = FALSE))
}
if(is.na(validationDataSheet$SplitData)){validationDataSheet$SplitData <- FALSE}
if(is.na(validationDataSheet$CrossValidate)){validationDataSheet$CrossValidate <- FALSE}
if(validationDataSheet$CrossValidate){
  if(is.na(validationDataSheet$NumberOfFolds)){
    warning("Number of Folds not specified. Default value used: 10")
    validationDataSheet$NumberOfFolds <- 10
  }
  if(is.na(validationDataSheet$StratifyFolds)){validationDataSheet$StratifyFolds <- FALSE}
}
# Prep inputs ------------------------------------------------------------------

# identify categorical covariates
if(sum(covariatesSheet$IsCategorical, na.rm = T)>0){
  factorVars <- covariatesSheet$CovariateName[which(covariatesSheet$IsCategorical == T)]
} else { factorVars <- NULL }


# template raster
template <- rast(templateRasterSheet$RasterFilePath)
tempEPSG <- strsplit(strsplit(crs(template), 'EPSG\",')[[1]][2], "]")[[1]][1]
tempRes <- res(template)

# Prep field data --------------------------------------------------------------

# TO DO: add code for spatial multiprocessing 

# ptsMat <- pts %>% st_as_sf(coords = c("X", "Y")) # st_is_longlat(ptsMat)
# ptsObj <- st_set_crs(ptsMat, paste0("EPSG:",epsg)) # fieldDataEPSG$EPSG # st_is_longlat(ptsObj)
# ptsObj <- st_transform(ptsObj, paste0("EPSG:",tempEPSG)) # plot(ptsObj)

# check if EPGS matches template if not update coordinates  
if(fieldDataOptions$EPSG != tempEPSG){
  pts <- vect(fieldDataSheet, geom = c("X", "Y"), crs = paste0("EPSG:",fieldDataOptions$EPSG))
  pts <- terra::project(pts, template)
  geomPts <- geom(pts)
  fieldDataSheet$X <- geomPts[,"x"]
  fieldDataSheet$Y <- geomPts[,"y"]
} else {
  pts <- vect(fieldDataSheet, geom = c("X", "Y"), crs = paste0("EPSG:",tempEPSG))
}

# only keep site inside template extent
extPoly <- as.polygons(ext(template))
pts2 <- crop(pts, extPoly)
keepSites <- pts2$SiteID

warning(paste0(nrow(fieldDataSheet)-length(keepSites), " sites out of ", nrow(fieldDataSheet), 
               " total sites in the input Field Data sheet were outside the template extent and were REMOVED from the output Field Data Sheet."))
fieldDataSheet <- fieldDataSheet[fieldDataSheet$SiteID %in% keepSites,]

# rasterize field data
r <- rast(ext(template), resolution = tempRes, crs = paste0("EPSG:",tempEPSG))
rastPts <- rasterize(pts2, r)
matPts <- as.matrix(rastPts, wide=T)
keep <- which(!is.na(matPts))

rSites <- rasterize(pts2, r, field = "SiteID")
matSts <- as.matrix(rSites, wide=T)
SiteID <- matSts[keep]

# check if cells encompass multiple sites ## TO DO -- need to sort out weighting/selecting sites that occur within same pixel
if(fieldDataOptions$AggregateAndWeight){
  rIDs <- rast(r, vals = 1:(dim(r)[1]*dim(r)[2]))
  cellPerPt <- terra::extract(rIDs, pts2)
  cellPerPt$SiteID <- pts2$SiteID
  multiPtCells <- unique(cellPerPt$lyr.1[duplicated(cellPerPt$lyr.1)])
  multiSiteCells <- list()
  checkSites <- NULL
  for(i in multiPtCells){
    sites_i <- cellPerPt$SiteID[cellPerPt$lyr.1 == i]
    multiSiteCells[[as.character(i)]]$sites <- sites_i # pts2$Response[pts2$SiteID %in% sites_i]
    multiSiteCells[[as.character(i)]]$response <- pts2$Response[pts2$SiteID %in% sites_i]
    checkSites <- c(checkSites, sites_i)
  }
  
  repeatPts <- rasterize(pts2, r, fun = length)
  repeatPts <- rasterize(pts2, r, field = "SiteID", fun = function(x){unique(x)})
  repMat <- as.matrix(repeatPts, wide=T)
  # repMat[which(repMat>1)]
  repMat[keep]
}

  
# Prep covaraite data ----------------------------------------------------------

# TO DO: prepare raster layers to ensure layers match crs/res/extent of template layer
## for now code assumes layers are already processed
rStack <- NULL
rStack <- try(rast(covariateDataSheet$RasterFilePath))
if("try-error" %in% class(rStack)){
  stop("Covariate rasters do not match extent, CRS and/or resolution.")
}
remove(rStack)

# Extract covariate site data --------------------------------------------------

# TO DO: account for every site -- current approach extracts covariate data per raster pixel containing one or more sites,
# then stores data for only one of the sites in the pixel. Need to allow for all sites to be retained and either weight 
# sites with spatial repetition or change them to background points

siteData <- data.frame(SiteID = matSts[keep])
# checkSites %in% siteData$SiteID
# siteData$ptCount <- repMat[keep]

for(i in 1:nrow(covariateDataSheet)){
  ri <- rast(covariateDataSheet$RasterFilePath[i])
  mi <- as.matrix(ri, wide=TRUE)
  
  outMat <- mi*matPts
  vals <- outMat[keep]
  siteData[covariateDataSheet$CovariatesID[i]] <- vals
}

# convert site data to long format to save in datasheet
siteDataWide <- siteData
siteData <- gather(data = siteData, key = CovariatesID, value = Value, -SiteID)

# save site data to scenario
saveDatasheet(myScenario, siteData, "wisdm_SiteData")


# Split data for testing/training and validation -------------------------------

inputData <- merge(fieldDataSheet, siteDataWide)

# Define Train/Test Split (if specified)
if(validationDataSheet$SplitData){
  inputData <- testTrainSplit(inputData = inputData,
                               trainProp = validationDataSheet$ProportionTrainingData,
                               # ratioPresAbs = validationDataSheet$RatioPresenceAbsence,
                               factorVars = factorVars)
}




# Define Cross Validation folds (if specified) 
if(validationDataSheet$CrossValidate){
  outputData <- crossValidationSplit(inputData = inputData,
                                     factorVars = factorVars,
                                     nFolds = validationDataSheet$NumberOfFolds,
                                     stratify = validationDataSheet$StratifyFolds)
}

updateFieldData <- select(outputData, names(fieldDataSheet))


# save updated field data to scenario
saveDatasheet(myScenario, updateFieldData, "wisdm_FieldData", append = F)
