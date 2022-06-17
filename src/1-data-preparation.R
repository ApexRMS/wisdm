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

# temp directory
ssimTempDir <- ssimEnvironment()$TransferDirectory 

# Read in datasheets
covariatesSheet <- datasheet(myProject, "wisdm_Covariates", optional = T)
fieldDataSheet <- datasheet(myScenario, "wisdm_FieldData", optional = T)
fieldDataOptions <- datasheet(myScenario, "wisdm_FieldDataOptions", optional = T)
covariateDataSheet <- datasheet(myScenario, "wisdm_CovariateData", lookupsAsFactors = F)
validationDataSheet <- datasheet(myScenario, "wisdm_ValidationOptions")
templateRasterSheet <- datasheet(myScenario, "wisdm_TemplateRaster")
# siteDataSheet <- datasheet(myScenario, "wisdm_SiteData", optional = T) # lookupsAsFactors = F,

#  Set defaults ----------------------------------------------------------------  

## field data options
if(nrow(fieldDataOptions)<1){
  fieldDataOptions <- addRow(fieldDataOptions, list(# Shapefile = NA,
                                                    EPSG = NA,
                                                    AggregateAndWeight = NA))
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

saveDatasheet(myScenario, validationDataSheet, "wisdm_ValidationOptions")

# Prep inputs ------------------------------------------------------------------

# identify categorical covariates
if(sum(covariatesSheet$IsCategorical, na.rm = T)>0){
  factorVars <- covariatesSheet$CovariateName[which(covariatesSheet$IsCategorical == T)]
} else { factorVars <- NULL }

# access crs database
# projDB <- system.file("proj/proj.db", package = "terra")
# crsTable <- sf::read_sf(projDB, "crs_view") 

# template raster
template <- rast(templateRasterSheet$RasterFilePath)
tempRes <- res(template)
tempExt <- ext(template)
tempCRS <- crs(template)

# Prep field data --------------------------------------------------------------

# TO DO: add code for spatial multiprocessing 

# ptsMat <- pts %>% st_as_sf(coords = c("X", "Y")) # st_is_longlat(ptsMat)
# st_crs(ptsMat)
# ptsObj <- st_set_crs(ptsMat, paste0("EPSG:",epsg)) # fieldDataEPSG$EPSG # st_is_longlat(ptsObj)
# ptsObj <- st_transform(ptsObj, paste0("EPSG:",tempEPSG)) # plot(ptsObj)

# proj_db <- system.file("proj/proj.db", package = "sf")
# crs_table <- sf::read_sf(proj_db, "crs_view") # extracts the "crs_view" table

# check if tabular or shapefile data was provided
# if(nrow(fieldDataSheet)>0){ # use tabular data
  
  # load points data and update coordinates if an authority code is defined  
  if(!is.na(fieldDataOptions$EPSG)){
    pts <- vect(fieldDataSheet, geom = c("X", "Y"), crs = fieldDataOptions$EPSG)
    pts <- terra::project(pts, template)
    geomPts <- geom(pts)
    fieldDataSheet$X <- geomPts[,"x"]
    fieldDataSheet$Y <- geomPts[,"y"]
  } else {
    pts <- vect(fieldDataSheet, geom = c("X", "Y"), crs = tempCRS)
  }


# } else { # use shapefile data
#   
#   unzip(fieldDataOptions$Shapefile, exdir = ssimTempDir)
#   tFiles <- list.files(ssimTempDir)
#   fileName <- tFiles[grepl("\\.shp$", tFiles)]
#   pts <- vect(file.path(ssimTempDir,fileName))
#   ptsEPSG <- strsplit(strsplit(crs(pts), 'EPSG\",')[[1]][2], "]")[[1]][1]
#
# df <- as.data.frame(pts[,c("X", "Y", "responseBi")])
# names(df) <- c("X", "Y", "Response")
# fieldDataSheet <- addRow(fieldDataSheet, df)
# fieldDataSheet$SiteID <- 1:nrow(fieldDataSheet)
# saveDatasheet(myScenario, fieldDataSheet, "wisdm_FieldData")
#
#   # check if EPGS matches template if not update coordinates  
#   if(!is.na(fieldDataOptions$EPSG)| ptsEPSG != tempEPSG){
#       pts <- terra::project(pts, template)
#       geomPts <- geom(pts)
#       fieldDataSheet$X <- geomPts[,"x"]
#       fieldDataSheet$Y <- geomPts[,"y"]
#       }
# }


# only keep site inside template extent
extPoly <- as.polygons(tempExt)
pts2 <- crop(pts, extPoly)
keepSites <- pts2$SiteID

if(length(keepSites)<nrow(fieldDataSheet)){
  warning(paste0(nrow(fieldDataSheet)-length(keepSites), " sites out of ", nrow(fieldDataSheet), 
               " total sites in the input Field Data were outside the template extent and were REMOVED from the output Field Data."))
  fieldDataSheet <- fieldDataSheet[fieldDataSheet$SiteID %in% keepSites,]
}

# rasterize field data
r <- rast(tempExt, resolution = tempRes, crs = tempCRS)
rastPts <- rasterize(pts2, r)
matPts <- as.matrix(rastPts, wide=T)
keep <- which(!is.na(matPts))

rIDs <- rast(r, vals = 1:(dim(r)[1]*dim(r)[2]))
cellPerPt <- terra::extract(rIDs, pts2)
cellPerPt$SiteID <- pts2$SiteID
names(cellPerPt)[2] <- "PixelID"

pts2 <- merge(pts2, cellPerPt[,c("PixelID", "SiteID")])
rPixels <- rasterize(pts2, r, field = "PixelID")
matPixs <- as.matrix(rPixels, wide=T)
PixelIDs <- matPixs[keep]

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

PixelData <- data.frame(PixelID = PixelIDs)

for(i in 1:nrow(covariateDataSheet)){
  ri <- rast(covariateDataSheet$RasterFilePath[i])
  mi <- as.matrix(ri, wide=TRUE)
  
  outMat <- mi*matPts
  vals <- outMat[keep]
  PixelData[covariateDataSheet$CovariatesID[i]] <- vals
}

# convert site data to long format to save in datasheet
siteData <- merge(cellPerPt[,c("PixelID", "SiteID")], PixelData)
siteDataWide <- siteData
siteData$PixelID <- NULL
siteData <- gather(data = siteData, key = CovariatesID, value = Value, -SiteID)

# save site data to scenario
saveDatasheet(myScenario, siteData, "wisdm_SiteData")

# Aggregate or Weight sites ----------------------------------------------------

repeatPIDs <- siteDataWide$PixelID[duplicated(siteDataWide$PixelID)]

if(!is.na(fieldDataOptions$AggregateAndWeight)){
  
  if(length(repeatPIDs)==0){
    message("Only one field data observation is represented per pixel; no aggregation or down-weighting required.")
    # fieldDataOptions$AggregateAndWeight <- NA
    # saveDatasheet(myScenario, fieldDataOptions, "wisdm_FieldDataOptions")
  } else {
  
  # if Aggregate sites is selected 
  if(fieldDataOptions$AggregateAndWeight == "Aggregate"){
    
    # if presence absence data
    if(all(na.omit(fieldDataSheet$Response) %in% 0:1)){ 
      for(i in repeatPIDs){
        sites_i <- siteDataWide$SiteID[siteDataWide$PixelID == i] 
        resp_i <- fieldDataSheet$Response[fieldDataSheet$SiteID %in% sites_i]
        if(sum(resp_i) == 0 | mean(resp_i) == 1){ # if all absence or all presence
          fieldDataSheet$Response[fieldDataSheet$SiteID %in% sites_i[2:length(sites_i)]] <- -9999 
        } else { # if response is mix of presence/absence 
          keep_i <- sites_i[resp_i == 1][1] # keep a presence and convert rest of repeat sites to background points 
          background_i <- sites_i[-which(sites_i == keep_i)]
          fieldDataSheet$Response[fieldDataSheet$SiteID %in% background_i] <- -9999 
        }  
      }
    } else { # if count data 
      for(i in repeatPIDs){
        sites_i <- siteDataWide$SiteID[siteDataWide$PixelID == i] 
        resp_i <- fieldDataSheet$Response[fieldDataSheet$SiteID %in% sites_i]
        if(sum(resp_i) == 0){ # if all counts are zero 
          fieldDataSheet$Response[fieldDataSheet$SiteID %in% sites_i[2:length(sites_i)]] <- -9999 
        } else if(any(resp_i > 0)){ # if any counts are greater then zero
          fieldDataSheet$Response[fieldDataSheet$SiteID %in% sites_i[1]] <- sum(resp_i)
          fieldDataSheet$Response[fieldDataSheet$SiteID %in% sites_i[2:length(sites_i)]] <- -9999 
        }
      }  
    }
  } # end site aggregation 
  
  if(fieldDataOptions$AggregateAndWeight == "Weight"){
    
    if(all(is.na(fieldDataSheet$Weight))){
      # set base weights
      fieldDataSheet$Weight <- 1 
    
      for(i in repeatPIDs){
        sites_i <- siteDataWide$SiteID[siteDataWide$PixelID == i] 
        weight_i <- 1/length(sites_i)
        fieldDataSheet$Weight[fieldDataSheet$SiteID %in% sites_i] <- weight_i
      } 
    } else { warning("Weights already provided in Field Data, new weights were NOT assigned.") }
  
  } # end site weighting
  } # end else
} # end aggregate and weight 

# Split data for testing/training and validation -------------------------------

inputData <- merge(fieldDataSheet, select(siteDataWide,-PixelID))

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
