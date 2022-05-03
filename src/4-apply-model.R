## --------------------
## wisdm - apply model
## ApexRMS, May 2022
## --------------------

# built under R version 4.1.1
# this transformer pulls in selected model objects and applies the models
# to specified spatial conditions 

# source dependencies ----------------------------------------------------------

packageDir <- Sys.getenv("ssim_package_directory")
# source(file.path(packageDir, "0-glm-constants.R"))
source(file.path(packageDir, "0-dependencies.R"))
source(file.path(packageDir, "0-helper-functions.R"))
source(file.path(packageDir, "0-apply-model-functions.R"))

# Connect to library -----------------------------------------------------------

# Active project and scenario
myLibrary <- ssimLibrary()
myProject <- rsyncrosim::project()
myScenario <- scenario()
# datasheet(myScenario)

# path to ssim directories
ssimTempDir <- paste0(Sys.getenv("ssim_temp_directory"),"\\Data\\")
ssimOutputDir <- Sys.getenv("ssim_output_directory")
resultScenario <- Sys.getenv("ssim_scenario_id")

# Read in datasheets
covariatesSheet <- datasheet(myProject, "Covariates", optional = T)
runControlSheet <- datasheet(myScenario, "RunControl", optional = T)
multiProcessingSheet <- datasheet(myScenario, "core_Multiprocessing")
# fieldDataSheet <- datasheet(myScenario, "wisdm_FieldData", optional = T)
# ValidationDataSheet <- datasheet(myScenario, "wisdm_ValidationOptions")
# reducedCovariatesSheet <- datasheet(myScenario, "wisdm_ReducedCovariates", lookupsAsFactors = F)
# siteDataSheet <- datasheet(myScenario, "wisdm_SiteData", lookupsAsFactors = F)
# GLMSheet <- datasheet(myScenario, "GLM")

covariateDataSheet <- datasheet(myScenario, "CovariateData", optional = T, lookupsAsFactors = F)
modelOutputsSheet <- datasheet(myScenario, "ModelOutputs", optional = T)
outputOptionsSheet <- datasheet(myScenario, "OutputOptions", optional = T)

spatialOutputsSheet <- datasheet(myScenario, "SpatialOutputs", optional = T)

# Set defaults -----------------------------------------------------------------

## run control 
runControlSheet <- addRow(runControlSheet, list(1,1,0,0))
saveDatasheet(myScenario, runControlSheet, "RunControl")

## output options
if(nrow(outputOptionsSheet)<1){
  outputOptionsSheet <- addRow(outputOptionsSheet, list(T,F,F))
}
if(any(is.na(outputOptionsSheet))){
  outputOptionsSheet[is.na(outputOptionsSheet)] <- F
}

# Load model object ------------------------------------------------------------

modType <- modelOutputsSheet$ModelType

mod <- readRDS(paste0(ssimOutputDir,"\\Scenario-", resultScenario,"\\wisdm_ModelOutputs\\",modelOutputsSheet$ModelRDS))

if(modType == "glm"){

  nVars <- length(attr(terms(formula(mod)),"term.labels"))
  modVars <- attr(terms(formula(mod)),"term.labels")
  # have to remove all the junk with powers and interactions for mess map production to work
  modVars <- unique(unlist(strsplit(gsub("I\\(","",gsub("\\^2)","",modVars)),":")))
  
  trainingData <-  mod$data 
}

# identify factor variables
factorVars <- covariatesSheet$CovariateName[which(covariatesSheet$IsCategorical == T)]
factorVars <- factorVars[factorVars %in% modVars]
if(length(factorVars)==0){ factorVars <- NULL }


# Predict model ----------------------------------------------------------------

  # get the paths off the raster files
  paths <- covariateDataSheet$RasterFilePath[which(covariateDataSheet$CovariatesID %in% modVars)]
  covData <- covariateDataSheet[which(covariateDataSheet$CovariatesID %in% modVars),]
  
  # check that all tifs are present
  if(all(file.exists(paths)) == F){
    missingTifs <- covData$CovariatesID[!file.exists(covData$RasterFilePath)]
    stop("ERROR: the following geotiff(s) are missing from Covariate Data:  ",
         paste(missingTifs, collapse=" ,"),sep="")
  }
  
  # check that we have read access to all tiffs
  if(sum(file.access(paths),mode=0)!=0){
    stop("ERROR: the following geotiff(s) are missing : ",
         paths[(file.access(paths)!=0),][1],sep="")
  }
  
  if(modType == "glm"){
    predictFct <- glm.predict
  }
  # if(modType == "rf")   {
  #   predictFct = rf.predict
  #   library(randomForest)
  # }
  # if(modType == "mars") {
  #   predictFct = mars.predict
  #   library(mda)
  # }
  # if(modType == "brt")  {
  #   model.covs<-levels(summary(out$mods$final.mod,plotit=FALSE)[,1])
  #   predictFct = brt.predict
  #   library(gbm)
  # }
  
  # create prediction maps
  proc.tiff(fit.model = mod,
            modType = modType,
            mod.vars = covData$CovariatesID,
            raster.files = covData$RasterFilePath,
            pred.fct = predictFct,
            output.options = outputOptionsSheet,
            factor.levels = factorVars,
            multiprocessing.cores = multiProcessingSheet$MaximumJobs,
            temp.directory = ssimTempDir,
            train.dat = trainingData, 
            tsize = 20,
            NAval = -3000)


# Save output maps -------------------------------------------------------------
  
  # possibleFolders <- c("ProbTiff", "MESSTiff", 'ModTiff', "ResidTiff")
  # tempFolders <- list.files(paste0(ssimTempDir))
  
  # add model Outputs to datasheet
  spatialOutputsSheet <- addRow(spatialOutputsSheet, 
                              list(Iteration = 1, Timestep = 0,
                                   ModelType = modType,
                                   ProbabilityRaster = paste0(ssimTempDir,"ProbTiff\\glm_prob_map.tif"),
                                   MessRaster = paste0(ssimTempDir,"MESSTiff\\glm_mess_map.tif")))
  
  saveDatasheet(myScenario, spatialOutputsSheet, "SpatialOutputs")
  
  
