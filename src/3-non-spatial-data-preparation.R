## -------------------------
## wisdm - data preparation
## ApexRMS, March 2023
## -------------------------

# built under R version 4.1.3 & SyncroSim version 2.4.18
# Script pulls in field data and splits sites into test/train or CV groupings

# source dependencies ----------------------------------------------------------

library(rsyncrosim)
library(terra)
library(tidyr)
library(dplyr)
library(pander)

packageDir <- (Sys.getenv("ssim_package_directory"))
source(file.path(packageDir, "03-non-spatial-data-prep-functions.R"))
source(file.path(packageDir, "03-background-data-functions.R"))

# Connect to library -----------------------------------------------------------

# Active project and scenario
myProject <- rsyncrosim::project()
myScenario <- scenario()

# temp directory
ssimTempDir <- ssimEnvironment()$TransferDirectory 

# Read in datasheets
templateSheet <- datasheet(myScenario, "TemplateRaster")
covariatesSheet <- datasheet(myProject, "wisdm_Covariates", optional = T)
covariateDataSheet <- datasheet(myScenario, "wisdm_CovariateData", optional = T, lookupsAsFactors = F)
fieldDataSheet <- datasheet(myScenario, "wisdm_FieldData", optional = T)
backgroundDataOptionsSheet <- datasheet(myScenario, "wisdm_BackgroundDataOptions", optional = T)
validationDataSheet <- datasheet(myScenario, "wisdm_ValidationOptions")
siteDataSheet <- datasheet(myScenario, "wisdm_SiteData", optional = T, lookupsAsFactors = F)

# Prep inputs ------------------------------------------------------------------

# identify categorical covariates
if(sum(covariatesSheet$IsCategorical, na.rm = T)>0){
  factorVars <- covariatesSheet$CovariateName[which(covariatesSheet$IsCategorical == T & covariatesSheet$CovariateName %in% siteDataSheet$CovariatesID)]
  if(length(factorVars)<1){ factorVars <- NULL }
} else { factorVars <- NULL }

# drop no data (-9999) sites that resulted from spatial aggregation 
fieldDataSheet <- fieldDataSheet[!fieldDataSheet$Response == -9999,] 

#  Set defaults ----------------------------------------------------------------  

## Field data options sheet
if(nrow(backgroundDataOptionsSheet)<1){
  backgroundDataOptionsSheet <- addRow(backgroundDataOptionsSheet, list(GenerateBackgroundSites = FALSE))
}
if(is.na(backgroundDataOptionsSheet$GenerateBackgroundSites)){ backgroundDataOptionsSheet$GenerateBackgroundSites <- FALSE }
if(backgroundDataOptionsSheet$GenerateBackgroundSites){
  if(is.na(backgroundDataOptionsSheet$BackgroundSiteCount)){backgroundDataOptionsSheet$BackgroundSiteCount <- sum(fieldDataSheet$Response) }  
  if(is.na(backgroundDataOptionsSheet$BackgroundGenerationMethod)){backgroundDataOptionsSheet$BackgroundGenerationMethod <- "Kernel Density Estimate (KDE)"}
  if(is.na(backgroundDataOptionsSheet$KDESurface)){
    if(backgroundDataOptionsSheet$BackgroundGenerationMethod == "Kernel Density Estimate (KDE)"){
      backgroundDataOptionsSheet$KDESurface <- "Continuous"
    }}
  if(is.na(backgroundDataOptionsSheet$Isopleth)){
    if(backgroundDataOptionsSheet$KDESurface == "Binary" | backgroundDataOptionsSheet$BackgroundGenerationMethod == "Minimum Convex Polygon (MCP)"){
      backgroundDataOptionsSheet$Isopleth <- 95
    }
  }
}

saveDatasheet(myScenario, backgroundDataOptionsSheet, "wisdm_BackgroundDataOptions")

## Validation Sheet
if(nrow(validationDataSheet)<1){
  validationDataSheet <- addRow(validationDataSheet, list(SplitData = FALSE,
                                                          CrossValidate = FALSE))
}
if(is.na(validationDataSheet$SplitData)){validationDataSheet$SplitData <- FALSE}
if(validationDataSheet$SplitData){
  if(is.na(validationDataSheet$ProportionTrainingData)){
    validationDataSheet$ProportionTrainingData <- 0.5
    updateRunLog("\nTraining proportion not specified. Default value used: 0.5\n")
  }
}
if(is.na(validationDataSheet$CrossValidate)){validationDataSheet$CrossValidate <- FALSE}
if(validationDataSheet$CrossValidate){
  if(is.na(validationDataSheet$NumberOfFolds)){
    updateRunLog("\nNumber of Folds not specified. Default value used: 10\n")
    validationDataSheet$NumberOfFolds <- 10
  }
  if(is.na(validationDataSheet$StratifyFolds)){validationDataSheet$StratifyFolds <- FALSE}
}

saveDatasheet(myScenario, validationDataSheet, "wisdm_ValidationOptions")


# Generate pseudo-absences (if applicable) -------------------------------------

if(backgroundDataOptionsSheet$GenerateBackgroundSites){
  
  if(backgroundDataOptionsSheet$BackgroundGenerationMethod == "Kernel Density Estimate (KDE)"){ 
    methodInputs <- list("method"="kde", 
                   "surface"=backgroundDataOptionsSheet$KDESurface, 
                   "isopleth"=backgroundDataOptionsSheet$Isopleth)
    }
  if(backgroundDataOptionsSheet$BackgroundGenerationMethod == "Minimum Convex Polygon (MCP)"){ 
    methodInputs <- list("method"="mcp",
                         "surface" = NA,
                   "isopleth"=backgroundDataOptionsSheet$Isopleth)
    }

  templateRaster <- rast(templateSheet$RasterFilePath)
  
  # create mask polygon from extent of valid data in template 
  templateVector <- as.polygons(templateRaster)
  
  # generate background surface
  backgroundSurfaceGeneration(sp = "species",
                              template = templateRaster,
                              method = methodInputs, 
                              mask = templateVector,
                              outputDir = ssimTempDir,
                              dat = fieldDataSheet)
  
  # generate background (psuedo-absence) points
  backgroundPointGeneration(sp = "species",
                            outputDir = ssimTempDir,
                            n = backgroundDataOptionsSheet$BackgroundSiteCount+100,
                            method = methodInputs,
                            # target_file = backgroundDataOptionsSheet, # this is an external input... 
                            overwrite = T)
  
  # add background point to field data
  bgData <- read.csv(file.path(ssimTempDir, paste0("species_", methodInputs$method, "_bg_pts.csv")))
  
  startId <- max(fieldDataSheet$SiteID)+1
  bgData$SiteID <- startId:(startId+nrow(bgData)-1)
  bgData$UseInModelEvaluation <- NA
  bgData$ModelSelectionSplit <- NA
  bgData$Weight <- NA
  
  fieldData <- rbind(fieldDataSheet, bgData)
  
  ## Remove background points that occur in pixels with presence points -----
  
  # rasterize points data
  r <- rast(ext(templateRaster), resolution = res(templateRaster), crs = crs(templateRaster))
  pts <- vect(fieldData, geom = c("X", "Y"), crs = crs(templateRaster))
  rastPts <- rasterize(pts, r)
  matPts <- as.matrix(rastPts, wide=T)
  keep <- which(!is.na(matPts))
  
  rIDs <- rast(r, vals = 1:(dim(r)[1]*dim(r)[2]))
  cellPerPt <- terra::extract(rIDs, pts)
  cellPerPt$SiteID <- pts$SiteID
  names(cellPerPt)[2] <- "PixelID"
  pts <- merge(pts, cellPerPt[,c("PixelID", "SiteID")])
  
  # check for duplicate pixel ids
  if(any(duplicated(pts$PixelID))){
    dups <- pts[which(duplicated(pts$PixelID)),]
    dropSites <- dups$SiteID[which(dups$Response == -9998)]
    pts <- pts[!pts$SiteID %in% dropSites,]
  }
  
  bgPts <- pts[pts$Response == -9998]
  bgPts$PixelID <- NULL
  
  # remove extra bg sites
  if(nrow(bgPts)>backgroundDataOptionsSheet$BackgroundSiteCount){
    bgPts <- sample(x = bgPts, size = backgroundDataOptionsSheet$BackgroundSiteCount, replace = F)
  }
  if(nrow(bgPts)<backgroundDataOptionsSheet$BackgroundSiteCount){
    updateRunLog(paste0(nrow(bgPts), " psuedoabsence sites added to Field Data when ", backgroundDataOptionsSheet$BackgroundSiteCount, 
                        " were requested."))
    # backgroundDataOptionsSheet$BackgroundSiteCount <- nrow(bgPts)
  }
  
  ## Extract covariate data for background sites  -----
  
  # rasterize bg data
  rastPts <- rasterize(bgPts, r)
  matPts <- as.matrix(rastPts, wide=T)
  keep <- which(!is.na(matPts))
  
  rIDs <- rast(r, vals = 1:(dim(r)[1]*dim(r)[2]))
  cellPerPt <- terra::extract(rIDs, bgPts)
  cellPerPt$SiteID <- bgPts$SiteID
  names(cellPerPt)[2] <- "PixelID"
  bgPts <- merge(bgPts, cellPerPt[,c("PixelID", "SiteID")])
  
  rPixels <- rasterize(bgPts, r, field = "PixelID")
  matPixs <- as.matrix(rPixels, wide=T)
  PixelIDs <- matPixs[keep]
  
  PixelData <- data.frame(PixelID = PixelIDs)
  
  for(i in 1:nrow(covariateDataSheet)){
    ri <- rast(covariateDataSheet$RasterFilePath[i])
    mi <- as.matrix(ri, wide=TRUE)
    
    outMat <- mi*matPts
    vals <- outMat[keep]
    PixelData[covariateDataSheet$CovariatesID[i]] <- vals
  }
  
  # convert background site data to long format and add to existing site datasheet
  bgSiteData <- merge(cellPerPt[,c("PixelID", "SiteID")], PixelData)
  bgSiteData$PixelID <- NULL
  bgSiteData <- gather(data = bgSiteData, key = CovariatesID, value = Value, -SiteID)
  
  siteDataSheet <- rbind(siteDataSheet, bgSiteData)
  
  # save site data to scenario
  saveDatasheet(myScenario, siteDataSheet, "wisdm_SiteData")
  
  # update field data
  bgData <- bgData[which(bgData$SiteID %in% bgPts$SiteID),]
  if(any(!is.na(fieldDataSheet$Weight))){ bgData$Weight <- 1 }

  fieldDataSheet <- rbind(fieldDataSheet, bgData)

  # update response for background sites
  bgSiteIds <- fieldDataSheet$SiteID[fieldDataSheet$Response == -9998]
  fieldDataSheet$Response[which(fieldDataSheet$SiteID %in% bgSiteIds)] <- 0

}

# Split data for testing/training and validation -------------------------------

siteDataWide <- spread(data = siteDataSheet, key = CovariatesID, value = Value)

inputData <- left_join(fieldDataSheet, siteDataWide) # select(siteDataWide,-PixelID))

# Define Train/Test Split (if specified) 
if(validationDataSheet$SplitData){
  inputData <- testTrainSplit(inputData = inputData,
                               trainProp = validationDataSheet$ProportionTrainingData,
                               # ratioPresAbs = validationDataSheet$RatioPresenceAbsence,
                               factorVars = factorVars)
} else {
  # set all data to be used in model training
  inputData$UseInModelEvaluation <- FALSE
}


# Define Cross Validation folds (if specified) 
if(validationDataSheet$CrossValidate){
  
  inputData <- crossValidationSplit(inputData = inputData,
                                     factorVars = factorVars,
                                     nFolds = validationDataSheet$NumberOfFolds,
                                     stratify = validationDataSheet$StratifyFolds)
}

# Check categorical variables and update field data sheet (if no validation specified)
if(validationDataSheet$SplitData == F & validationDataSheet$CrossValidate == F){
  if(!is.null(factorVars)){
    for (i in 1:length(factorVars)){
      factor.table <- table(inputData[,factorVars[i]])
      if(any(factor.table<10)){
        updateRunLog(paste0("\nSome levels for the categorical predictor ",factorVars[i]," do not have at least 10 observations. ", 
                                                    "Consider removing or reclassifying this predictor before continuing.\n",
                                                    "Factors with few observations can cause failure in model fitting when the data is split and cannot be reliably used in training a model.\n"))
        factor.table <- as.data.frame(factor.table)
        colnames(factor.table)<-c("Factor Name","Factor Count")
        updateRunLog(paste("\n",factorVars[i],"\n"))
        updateRunLog(pander::pandoc.table.return(factor.table, style = "rmarkdown"))
      }
    }
  }
  # set all data to be used in model training
  inputData$UseInModelEvaluation <- FALSE
}  


# revert response for background sites
updateFieldData <- dplyr::select(inputData, names(fieldDataSheet))
if(backgroundDataOptionsSheet$GenerateBackgroundSites){
  updateFieldData$Response[which(updateFieldData$SiteID %in% bgSiteIds)] <- -9998
}

# save updated field data to scenario
saveDatasheet(myScenario, updateFieldData, "wisdm_FieldData", append = F)
