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
source(file.path(packageDir, "04-background-data-functions.R"))

updateRunLog('4 - Background Data Generation => Begin')

# Connect to library -----------------------------------------------------------

# Active project and scenario
myProject <- rsyncrosim::project()
myScenario <- scenario()

# temp directory
ssimTempDir <- ssimEnvironment()$TransferDirectory 

# Read in datasheets
templateSheet <- datasheet(myScenario, "wisdm_TemplateRaster")
covariateDataSheet <- datasheet(myScenario, "wisdm_CovariateData", optional = T, lookupsAsFactors = F)
fieldDataSheet <- datasheet(myScenario, "wisdm_FieldData", optional = T)
backgroundDataOptionsSheet <- datasheet(myScenario, "wisdm_BackgroundDataOptions", optional = T)
siteDataSheet <- datasheet(myScenario, "wisdm_SiteData", optional = T, lookupsAsFactors = F)

# Set progress bar -------------------------------------------------------------

steps <- 5 + length(covariateDataSheet$CovariatesID)
progressBar(type = "begin", totalSteps = steps)

# Prep inputs ------------------------------------------------------------------

# drop no data (-9999) sites that resulted from spatial aggregation 
fieldDataSheet <- fieldDataSheet[!fieldDataSheet$Response == -9999,] 

#  Set defaults ----------------------------------------------------------------  

## Background data options sheet
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
  
  progressBar()
  rm(templateVector); gc()
  
  # generate background (psuedo-absence) points
  backgroundPointGeneration(sp = "species",
                            outputDir = ssimTempDir,
                            n = backgroundDataOptionsSheet$BackgroundSiteCount+100,
                            method = methodInputs,
                            # target_file = backgroundDataOptionsSheet, # this is an external input... 
                            overwrite = T)
  
  progressBar()

  # add background point to field data
  bgData <- read.csv(file.path(ssimTempDir, paste0("species_", methodInputs$method, "_bg_pts.csv")))
  
  startId <- max(fieldDataSheet$SiteID, siteDataSheet$SiteID)+1
  bgData$SiteID <- startId:(startId+nrow(bgData)-1)
  bgData$UseInModelEvaluation <- NA
  bgData$ModelSelectionSplit <- NA
  bgData$Weight <- NA
  
  fieldData <- rbind(fieldDataSheet, bgData)
  
  ## Remove background points that occur in pixels with presence points -----
  
  # rasterize points data
  r <- rast(ext(templateRaster), resolution = res(templateRaster), crs = crs(templateRaster))
  pts <- vect(fieldData, geom = c("X", "Y"), crs = crs(templateRaster))
  rm(templateRaster, fieldData); gc()
  
  rIDs <- rast(r, vals = 1:(dim(r)[1]*dim(r)[2]))
  cellPerPt <- terra::extract(rIDs, pts)
  cellPerPt$SiteID <- pts$SiteID
  names(cellPerPt)[2] <- "PixelID"
  pts <- merge(pts, cellPerPt[,c("PixelID", "SiteID")])
  rm(rIDs, cellPerPt); gc()
  
  # check for duplicate pixel ids
  if(any(duplicated(pts$PixelID))){
    dups <- pts[which(duplicated(pts$PixelID)),]
    dropSites <- dups$SiteID[which(dups$Response == -9998)]
    pts <- pts[!pts$SiteID %in% dropSites,]
  }
  
  bgPts <- pts[pts$Response == -9998]
  bgPts$PixelID <- NULL
  rm(pts, dups, dropSites); gc()
  
  # remove extra bg sites
  if(nrow(bgPts)>backgroundDataOptionsSheet$BackgroundSiteCount){
    bgPts <- sample(x = bgPts, size = backgroundDataOptionsSheet$BackgroundSiteCount, replace = F)
  }
  if(nrow(bgPts)<backgroundDataOptionsSheet$BackgroundSiteCount){
    updateRunLog(paste0(nrow(bgPts), " psuedoabsence sites added to Field Data when ", backgroundDataOptionsSheet$BackgroundSiteCount, 
                        " were requested."))
    # backgroundDataOptionsSheet$BackgroundSiteCount <- nrow(bgPts)
  }
  
  progressBar()

  ## Extract covariate data for background sites  -----
  
  # rasterize bg data
  rastPts <- rasterize(bgPts, r)
  matPts <- as.matrix(rastPts, wide=T)
  keep <- which(!is.na(matPts))
  rm(rastPts); gc()
  
  rIDs <- rast(r, vals = 1:(dim(r)[1]*dim(r)[2]))
  cellPerPt <- terra::extract(rIDs, bgPts)
  cellPerPt$SiteID <- bgPts$SiteID
  names(cellPerPt)[2] <- "PixelID"
  bgPts <- merge(bgPts, cellPerPt[,c("PixelID", "SiteID")])
  rm(rIDs); gc()
  
  rPixels <- rasterize(bgPts, r, field = "PixelID")
  matPixs <- as.matrix(rPixels, wide=T)
  PixelIDs <- matPixs[keep]
  PixelData <- data.frame(PixelID = PixelIDs)
  rm(rPixels, matPixs, PixelIDs); gc()
  
  for(i in 1:nrow(covariateDataSheet)){
    ri <- rast(covariateDataSheet$RasterFilePath[i])
    mi <- as.matrix(ri, wide=TRUE)
    
    outMat <- mi*matPts
    vals <- outMat[keep]
    rm(ri, mi, outMat); gc()
    
    PixelData[covariateDataSheet$CovariatesID[i]] <- vals
    progressBar()
    
  }
  
  # convert background site data to long format and add to existing site datasheet
  bgSiteData <- merge(cellPerPt[,c("PixelID", "SiteID")], PixelData)
  bgSiteData$PixelID <- NULL
  bgSiteData <- gather(data = bgSiteData, key = CovariatesID, value = Value, -SiteID)
  
  siteDataSheet <- rbind(siteDataSheet, bgSiteData)
  siteDataSheet$SiteID <- format(siteDataSheet$SiteID, scientific = F)
  siteDataSheet <- siteDataSheet %>% distinct(SiteID, CovariatesID, .keep_all = T)
  rm(bgSiteData); gc()
  
  # save site data to scenario
  saveDatasheet(myScenario, siteDataSheet, "wisdm_SiteData")
  rm(siteDataSheet); gc()
  
  # update field data
  bgData <- bgData[which(bgData$SiteID %in% bgPts$SiteID),]
  if(any(!is.na(fieldDataSheet$Weight))){ bgData$Weight <- 1 }

  fieldDataSheet <- rbind(fieldDataSheet, bgData)
  fieldDataSheet$SiteID <- format(fieldDataSheet$SiteID, scientific = F)
  rm(bgData, bgPts); gc()
  
  # save updated field data to scenario
  saveDatasheet(myScenario, fieldDataSheet, "wisdm_FieldData", append = F)
}

progressBar(type = "end")

