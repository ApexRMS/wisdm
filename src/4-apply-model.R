## --------------------
## wisdm - apply model
## ApexRMS, May 2022
## --------------------

# built under R version 4.1.1
# this transformer pulls in selected model objects and applies the models
# to specified spatial conditions 

# source dependencies ----------------------------------------------------------

packageDir <- Sys.getenv("ssim_package_directory")
source(file.path(packageDir, "0-dependencies.R"))
source(file.path(packageDir, "0-helper-functions.R"))
source(file.path(packageDir, "04-apply-model-functions.R"))

# Connect to library -----------------------------------------------------------

# Active project and scenario
myLibrary <- ssimLibrary()
myProject <- rsyncrosim::project()
myScenario <- scenario()
# datasheet(myScenario)

# path to ssim directories
ssimTempDir <- ssimEnvironment()$TransferDirectory 
ssimOutputDir <- ssimEnvironment()$OutputDirectory                                 
resultScenario <- ssimEnvironment()$ScenarioId

# Read in datasheets
covariatesSheet <- datasheet(myProject, "Covariates", optional = T)
# runControlSheet <- datasheet(myScenario, "RunControl", optional = T)
multiprocessingSheet <- datasheet(myScenario, "core_Multiprocessing")
spatialMulitprocessingSheet <- datasheet(myScenario, "corestime_Multiprocessing")
covariateDataSheet <- datasheet(myScenario, "CovariateData", optional = T, lookupsAsFactors = F)
templateSheet <- datasheet(myScenario, "TemplateRaster")
modelOutputsSheet <- datasheet(myScenario, "ModelOutputs", optional = T)
outputOptionsSheet <- datasheet(myScenario, "OutputOptions", optional = T)
spatialOutputsSheet <- datasheet(myScenario, "SpatialOutputs", optional = T)

# Set defaults -----------------------------------------------------------------

## Run control sheet
# if(nrow(runControlSheet)<1){
#   runControlSheet <- addRow(runControlSheet, list(1,1,0,0))
#   saveDatasheet(myScenario, runControlSheet, "RunControl")
# }
## Output options sheet
if(nrow(outputOptionsSheet)<1){
  outputOptionsSheet <- addRow(outputOptionsSheet, list(T,F,F,F))
}
if(any(is.na(outputOptionsSheet))){
  outputOptionsSheet[is.na(outputOptionsSheet)] <- F
}

# set up spatial multiprocessing -----------------------------------------------

if(name(myLibrary) == "Partial"){
  
  # load multiprocessing raster
  maskFile <- rast(datasheet(myScenario, "corestime_Multiprocessing")$MaskFileName)
  # maskFile <- rast(file.choose())
  
  # identify multiprocessing tile
  tileID <- as.numeric(strsplit(strsplit(basename(ssimEnvironment()$LibraryFilePath), split = "-")[[1]][2],"\\.")[[1]][1])
  # tileID <- 1
  
  # define masking values
  maskValues <- unique(maskFile)[,1]
  if(length(maskValues)>1){
    maskValues <- maskValues[!maskValues == tileID]
  } else { maskValues <- NULL }
} else { maskValues <- NULL }

# Load model object ------------------------------------------------------------

modType <- modelOutputsSheet$ModelType

mod <- readRDS(paste0(ssimOutputDir,"\\Scenario-", resultScenario,"\\wisdm_ModelOutputs\\",modelOutputsSheet$ModelRDS))

if(modType == "glm"){
  
  nVars <- length(attr(terms(formula(mod)),"term.labels"))
  modVars <- attr(terms(formula(mod)),"term.labels")
  # have to remove all the junk with powers and interactions for mess map production to work
  modVars <- unique(unlist(strsplit(gsub("I\\(","",gsub("\\^2)","",modVars)),":")))
  
  trainingData <-  mod$data
  if(outputOptionsSheet$MakeResidualsMap){
    trainingData$predicted <- pred.fct(x=trainingData, mod=mod, modType=modType)
    if(max(trainingData$Response)>1){ modFamily <-"poisson" 
    } else { modFamily <- "binomial" }
  }
}

# identify factor variables
factorVars <- covariatesSheet$CovariateName[which(covariatesSheet$IsCategorical == T)]
factorVars <- factorVars[factorVars %in% modVars]
if(length(factorVars)==0){ factorVars <- NULL }


# Prepare inputs ---------------------------------------------------------------

# get the paths off the raster files
paths <- covariateDataSheet$RasterFilePath[which(covariateDataSheet$CovariatesID %in% modVars)]
covData <- covariateDataSheet[which(covariateDataSheet$CovariatesID %in% modVars),]

# check that all tifs are present
if(all(file.exists(paths)) == F){
  missingTifs <- covData$CovariatesID[!file.exists(covData$RasterFilePath)]
  stop("The following geotiff(s) are missing from Covariate Data:  ",
       paste(missingTifs, collapse=" ,"),sep="")
}

# check that we have read access to all tiffs
if(sum(file.access(paths),mode=0)!=0){
  stop("The following geotiff(s) are missing : ",
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


# prep prediction data ---------------------------------------------------------

if(nrow(templateSheet)<1){ stop("Template raster is missing") }

templateRaster <- rast(templateSheet$RasterFilePath)
  
if(!is.null(maskValues)){
  templateRaster <- mask(templateRaster, mask = maskFile, maskvalues = maskValues)
}
templateValues <- values(templateRaster)

if(all(is.na(templateValues))){  # if the template is completely NA values, don't read in any other data
  
  temp <- rep(NA, times = nrow(templateRaster) * ncol(templateRaster))
  
  } else {
  
    temp <- data.frame(matrix(ncol = nVars, nrow = nrow(templateRaster) * ncol(templateRaster)))
    
    for(k in 1:nVars){
      rast_k <- rast(covData$RasterFilePath[k])
      if(!is.null(maskValues)){
        rast_k <- mask(rast_k, mask = maskFile, maskvalues = maskValues)
      }
      temp[,k] <- as.vector(values(rast_k))
    }
  
    names(temp) <- covData$CovariatesID
    
  }

# Set predictor values to NA where ever the mask is NA
temp[is.na(templateValues),] <- NA 
# replace missing values with NA
for(m in covData$CovariatesID){ temp[is.nan(temp[,m]),m] <- NA }

if(!is.null(factorVars)){
  factor.cols <- match(names(factorVars), names(temp))
  if(sum(!is.na(factor.cols)) > 0){
    for(j in 1:length(factor.cols)){
      if(!is.na(factor.cols[j])){
        temp[,factor.cols[j]] <- factor(temp[, factor.cols[j]], levels = factorVars[[j]]$number, labels = factorVars[[j]]$class)
      }
    }
  }
}  

# create probability map -------------------------------------------------------
  
if(outputOptionsSheet$MakeProbabilityMap){
  
  preds <- t(matrix(predictFct(model = mod, x = temp), ncol = ncol(templateRaster), byrow = T))
  preds <- round((preds*100), 0)
  preds[is.na(preds)] <- -9999
  # typeof(preds)
  
  probRaster <- rast(templateRaster, vals = preds)
  # is.int(probRaster) 
  writeRaster(x = probRaster, 
              filename = file.path(ssimTempDir, paste0(modType,"_prob_map.tif")), 
              datatype = "INT4S",
              overwrite = TRUE)
}
  
# create MESS and MoD maps -----------------------------------------------------  
  
if(outputOptionsSheet$MakeMessMap){
  
  # order the training data so that we can consider the first and last row only in mess calculations
  train.dat <- trainingData[ ,match(modVars, names(trainingData))]
  for(k in 1:nVars){ train.dat[ ,k] <- sort(train.dat[ ,k]) } 
  # index <- names(train.dat)
  
  pred.rng <- rep(NA, nrow(temp))
  names(pred.rng) <- NA
    
  if(any(complete.cases(temp))){
      MessVals <- CalcMESS(rast = temp[complete.cases(temp), ], train.dat = train.dat)
      pred.rng[complete.cases(temp)] <- MessVals[ ,2]
      names(pred.rng)[complete.cases(temp)] <- MessVals[ ,1]
  }
  pred.rng <- round(pred.rng,0)
  pred.rng[is.na(pred.rng)] <- -9999
  # typeof(pred.rng)
  
  messRaster <- rast(templateRaster, vals = pred.rng)
  writeRaster(x = messRaster, 
                filename = file.path(ssimTempDir, paste0(modType,"_mess_map.tif")), 
                datatype = "INT4S",
                overwrite = TRUE) 
  # is.int(messRaster) 
  
  if(outputOptionsSheet$MakeModMap){ 
    
    if(is.null(names(pred.rng))){ names(pred.rng) <- NA }
    vals <- as.integer(names(pred.rng))
    vals[is.na(vals)] <- -9999
    # typeof(vals)
    
    modRaster <- rast(templateRaster, vals = vals)
    writeRaster(x = modRaster, 
                filename = file.path(ssimTempDir, paste0(modType,"_mod_map.tif")), 
                datatype = "INT4S", # "INT1U"
                overwrite = TRUE) 
    # is.int(modRaster) 
  }
}

  
# create residuals map ---------------------------------------------------------
  if(outputOptionsSheet$MakeResidualsMap){
    
    residSmooth <- readRDS(file.path(ssimEnvironment()$OutputDirectory, paste0("Scenario-",ssimEnvironment()$ScenarioId), "wisdm_ModelOutputs", modelOutputsSheet$ResidualSmoothRDS))

    predrast <- rast(probRaster)
    ablock <- 1:(ncol(probRaster) * nrow(probRaster))
    
    p <- xyFromCell(predrast, ablock) 
    p <- na.omit(p)
    
    blockvals <- data.frame(x=p[,1], y=p[,2])
  
    predv <- predict(residSmooth, blockvals)
    predv <- as.numeric(predv)
    predv[is.na(predv)] <- - 9999
    
    residRaster <- rast(templateRaster, vals = predv)
    writeRaster(x = residRaster, 
                filename = file.path(ssimTempDir, paste0(modType,"_resid_map.tif")), 
                datatype = "INT4S",
                overwrite = TRUE) 
    
}  


# Save output maps -------------------------------------------------------------
  
  # add model Outputs to datasheet
  spatialOutputsSheet <- addRow(spatialOutputsSheet, 
                              list( # Iteration = 1, Timestep = 0,
                                   ModelType = modType))
                              
  if(outputOptionsSheet$MakeProbabilityMap){
    spatialOutputsSheet$ProbabilityRaster <- file.path(ssimTempDir,"glm_prob_map.tif")
  }
  if(outputOptionsSheet$MakeMessMap){
    spatialOutputsSheet$MessRaster <- file.path(ssimTempDir,"glm_mess_map.tif")
    if(outputOptionsSheet$MakeModMap){
      spatialOutputsSheet$ModRaster <- file.path(ssimTempDir,"glm_mod_map.tif")
    }
  }
  if(outputOptionsSheet$MakeResidualsMap){
    spatialOutputsSheet$ResidualsRaster <- file.path(ssimTempDir,"glm_resid_map.tif")
  }
  
  # save outputs
  saveDatasheet(myScenario, spatialOutputsSheet, "SpatialOutputs")
    
