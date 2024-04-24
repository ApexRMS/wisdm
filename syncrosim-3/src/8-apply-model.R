## ---------------------
## wisdm - apply model
## ApexRMS, March 2024
## ---------------------

# built under R version 4.1.3 & SyncroSim version 3.0.0
# this transformer pulls in selected model objects and applies the models
# to specified spatial conditions to produce maps of occurrence probability, 
# multivariate environmental similarity surface, most dissimilar variable, 
# and residuals
# time tracking {dev code}

# Initialize first breakpoint for timing code
currentBreakPoint <- proc.time()

# source dependencies ----------------------------------------------------------

library(rsyncrosim)
library(terra)
library(tidyr)
library(dplyr)

packageDir <- Sys.getenv("ssim_package_directory")
source(file.path(packageDir, "00-helper-functions.R"))
source(file.path(packageDir, "08-apply-model-functions.R"))

# Connect to library -----------------------------------------------------------

# Active project and scenario
myLibrary <- ssimLibrary()
myScenario <- scenario()
# datasheet(myScenario)

# path to ssim directories
ssimTempDir <- ssimEnvironment()$TransferDirectory 
resultScenario <- ssimEnvironment()$ScenarioId

# Read in datasheets
covariatesSheet <- datasheet(myScenario, "wisdm_Covariates", optional = T, includeKey = T)
modelsSheet <- datasheet(myScenario, "wisdm_Models", includeKey = T)
multiprocessingSheet <- datasheet(myScenario, "core_Multiprocessing")
spatialMulitprocessingSheet <- datasheet(myScenario, "core_SpatialMultiprocessing")
covariateDataSheet <- datasheet(myScenario, "wisdm_CovariateData", optional = T, lookupsAsFactors = F)
templateSheet <- datasheet(myScenario, "wisdm_TemplateRaster")
restrictionSheet <- datasheet(myScenario, "wisdm_RestrictionRaster")
modelOutputsSheet <- datasheet(myScenario, "wisdm_OutputModel", optional = T, returnInvisible = T, lookupsAsFactors = F)
outputOptionsSheet <- datasheet(myScenario, "wisdm_OutputOptions", optional = T)

spatialOutputsSheet <- datasheet(myScenario, "wisdm_OutputSpatial", optional = T, lookupsAsFactors = T) %>% drop_na()
# spatialOutputsSheet$ModelsID <- as.numeric(spatialOutputsSheet$ModelsID)

# Set progress bar -------------------------------------------------------------

steps <- (nrow(modelOutputsSheet) * sum(outputOptionsSheet == T, na.rm = T))+2
progressBar(type = "begin", totalSteps = steps)

# Set defaults -----------------------------------------------------------------

## Covariates sheet
if(any(is.na(covariatesSheet$ID))){
  if (all(is.na(covariatesSheet$ID))){
    covariatesSheet$ID <- 1:nrow(covariatesSheet)
  } else {
    whichNA <- which(is.na(covariatesSheet$ID))
    maxID <- max(covariatesSheet$ID, na.rm = T)
    covariatesSheet$ID[whichNA] <- (maxID+1):(maxID+length(whichNA))
  }
 saveDatasheet(myScenario, covariatesSheet,  "wisdm_Covariates")  
}
 
## Output options sheet
if(nrow(outputOptionsSheet)<1){
  outputOptionsSheet <- addRow(outputOptionsSheet, list(T))
}
if(is.na(outputOptionsSheet$MakeBinaryMap)){ outputOptionsSheet$MakeBinaryMap <- F }
if(outputOptionsSheet$MakeBinaryMap){
  if(is.na(outputOptionsSheet$ThresholdOptimization)){
    outputOptionsSheet$ThresholdOptimization <- "Sensitivity equals specificity"
  }
}
if(is.na(outputOptionsSheet$MakeResidualsMap)){ outputOptionsSheet$MakeResidualsMap <- F }
if(is.na(outputOptionsSheet$MakeMessMap)){ outputOptionsSheet$MakeMessMap <- F }
if(is.na(outputOptionsSheet$MakeModMap)){ outputOptionsSheet$MakeModMap <- F }

saveDatasheet(myScenario, outputOptionsSheet,  "wisdm_OutputOptions") 

updateRunLog("Finished loading inputs in ", updateBreakpoint())

# set up spatial multiprocessing -----------------------------------------------

if(name(myLibrary) == "Partial"){
  
  # load multiprocessing raster
  maskFile <- rast(datasheet(myScenario, "core_SpatialMultiprocessing")$MaskFileName)
  maskExt <- ext(maskFile)
  
  # identify multiprocessing tile
  tileID <- as.numeric(strsplit(strsplit(basename(ssimEnvironment()$LibraryFilePath), split = "-")[[1]][2],"\\.")[[1]][1])
 
  # define masking values
  maskValues <- unique(maskFile)[,1]
  if(length(maskValues)>1){
    maskValues <- maskValues[!maskValues == tileID]
  } else { maskValues <- NULL }

  # trim tiling raster to extent of single tile
  maskFile <- mask(x = maskFile, mask = maskFile,  maskvalues = maskValues) # plot(maskFile)
  maskFile <- trim(x = maskFile)
 
} else { maskValues <- NULL }

# prepare restriction layer ---------------------------------------------------- 
  
# if(nrow((restrictionSheet))>0){
if(!is.na(restrictionSheet$RasterFilePath)){
  
   restrictRaster <- rast(restrictionSheet$RasterFilePath)
   
   if(!is.null(maskValues)){
     restrictRaster <- crop(restrictRaster, maskFile)
   }
   
   NAflag(restrictRaster) <- -9999
   
} else { restrictRaster <- NULL }

# Load model object ------------------------------------------------------------


for (i in 1:nrow(modelOutputsSheet)){ 
  
  
  modType <- modelsSheet$ModelType[modelsSheet$ModelName == modelOutputsSheet$ModelsID[i]]
  
  mod <- readRDS(modelOutputsSheet$ModelRDS[i])
  
  if(modType == "glm"){
    modVars <- attr(terms(formula(mod)),"term.labels")
    # have to remove all the junk with powers and interactions for mess map production to work
    modVars <- unique(unlist(strsplit(gsub("I\\(","",gsub("\\^2)","",modVars)),":")))
    trainingData <-  mod$data
    if(max(trainingData$Response)>1){ modFamily <-"poisson" 
    } else { modFamily <- "binomial" }
  }
  if(modType == "rf"){
    modVars <- rownames(mod$importance)
    trainingData <-  mod$trainingData
    if(max(trainingData$Response)>1){ modFamily <-"poisson" 
    } else { modFamily <- "binomial" }
  }
  if(modType == "maxent"){
    modVars <- mod$Raw.coef$V1
    trainingData <-  mod$trainingData
    mod$trainingData <- NULL
    modFamily <- "binomial"
  }
  if(modType == "brt"){
    modVars <- mod$contributions$var
    trainingData <-  mod$trainingData
    # mod$trainingData <- NULL
    # modFamily <- "binomial"
  }
  if(modType == "gam"){
    modVars <- attr(mod$terms, "term.labels")
    trainingData <-  mod$trainingData
    # mod$trainingData <- NULL
    # modFamily <- "binomial"
  }
  # identify factor variables
  factorVars <- covariatesSheet$CovariateName[which(covariatesSheet$IsCategorical == T)]
  factorVars <- factorVars[factorVars %in% modVars]
  if(length(factorVars)==0){ factorVars <- NULL }
  
  
  # Prepare inputs ---------------------------------------------------------------
  
  # get the paths off the raster files
  paths <- covariateDataSheet$RasterFilePath[which(covariateDataSheet$CovariatesID %in% modVars)]
  covData <- covariateDataSheet[which(covariateDataSheet$CovariatesID %in% modVars),]
  covData <- covData[!duplicated(covData),]
  
  # check that all tifs are present
  if(any(file.exists(paths)) == F){
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
  if(modType == "rf"){
    predictFct = rf.predict
    library(randomForest)
  }
  if(modType == "maxent"){
    predictFct = maxent.predict
  }
  if(modType == "brt"){
    predictFct = brt.predict
  }
  if(modType == "gam"){
    predictFct = gam.predict
  }
  # if(modType == "mars") {
  #   predictFct = mars.predict
  #   library(mda)
  # }
  
  # prep prediction data ---------------------------------------------------------
  
  # TO Do: Update scalability
  
  if(nrow(templateSheet)<1){ stop("Template raster is missing. Please provide a template raster before continuing.") }
  
  templateRaster <- rast(templateSheet$RasterFilePath)
    
  if(!is.null(maskValues)){
    templateRaster <- crop(templateRaster, maskFile)
   }
  templateValues <- values(templateRaster, dataframe=T)
  
  if(all(is.na(templateValues))){  # if the template is completely NA values, don't read in any other data
    
    temp <- rep(NA, times = nrow(templateRaster) * ncol(templateRaster))
    
    } else {
    
      temp <- data.frame(matrix(ncol = nrow(covData), nrow = nrow(templateRaster) * ncol(templateRaster)))
      
      for(k in 1:nrow(covData)){
        rast_k <- rast(covData$RasterFilePath[k])
        if(!is.null(maskValues)){
          rast_k <- crop(rast_k, maskFile)
          readStart(rast_k)
          temp[,k] <- readValues(rast_k) #, row = startRow, nrows = nrow(templateRaster))
          readStop(rast_k)
        } else {
          temp[,k] <- as.vector(values(rast_k))
        }
      }
      
      names(temp) <- covData$CovariatesID
      remove(rast_k)
    }
  
  # Set predictor values to NA where ever the mask is NA
  temp[is.na(templateValues),] <- NA 
  # replace missing values with NA
  for(m in covData$CovariatesID){ temp[is.nan(temp[,m]),m] <- NA }
  
  # set factor data
  if(!is.null(factorVars)){
    for(j in factorVars){
      lvls <- levels(mod$data[[j]])
      temp[,j] <- factor(temp[,j], levels = lvls)
    }
  }
  
  updateRunLog("Finished preparing mapping inputs in ", updateBreakpoint())
  
  # create probability map -----------------------------------------------------
    
  if(outputOptionsSheet$MakeProbabilityMap){ 
    preds <- matrix(predictFct(model = mod, x = temp), ncol = ncol(templateRaster), byrow = T)
    
    # apply restriction raster, if provided
    if(!is.null(restrictRaster)){
      readStart(restrictRaster)
      resVals <- readValues(restrictRaster, row = 1, nrows = nrow(templateRaster))
      readStop(restrictRaster)
      resVals <-  matrix(resVals, ncol = ncol(templateRaster), byrow = T)
      preds <- preds*resVals
    }
    preds <- round((preds*100), 0)
    # preds[is.na(preds)] <- -9999  # typeof(preds)
    
    probRaster <- rast(templateRaster, vals = preds) 
    
    if(!is.null(maskValues)){probRaster <- extend(probRaster, maskExt)}
    # is.int(probRaster) 
    writeRaster(x = probRaster, 
                filename = file.path(ssimTempDir, paste0(modType,"_prob_map.tif")), 
                datatype = "INT4S",
                NAflag = -9999,
                overwrite = TRUE)
    progressBar()
    updateRunLog("Finished Probability Map in ", updateBreakpoint())
  }
  
  # create binary map ---------------------------------------------------------- 
  
  if(outputOptionsSheet$MakeBinaryMap){
    
    if(!outputOptionsSheet$MakeProbabilityMap){
      updateRunLog("\nWarning: Binary map cannot be generated without generating the Probability map.\n")
    } else {
      
      thresholds <- mod$binThresholds
      names(thresholds) <- c("Max kappa", "Max sensitivity and specificity", "No omission", 
                           "Prevalence", "Sensitivity equals specificity")
      
      binThreshold <- as.numeric(thresholds[outputOptionsSheet$ThresholdOptimization])
      
      probVals <- preds/100
      probVals <- ifelse(probVals >= binThreshold, 1, 0)
      # probVals[is.na(probVals)] <- -9999
      binRaster <- rast(templateRaster, vals = probVals)
      
      if(!is.null(maskValues)){binRaster <- extend(binRaster, maskExt)}
      writeRaster(x = binRaster, 
                  filename = file.path(ssimTempDir, paste0(modType,"_bin_map.tif")), 
                  datatype = "INT4S", 
                  NAflag = -9999,
                  overwrite = TRUE) 
      progressBar()
      updateRunLog("Finished Binary Map in ", updateBreakpoint())
    
    }
  }
  
    
  # create MESS and MoD maps ---------------------------------------------------  
    
  if(outputOptionsSheet$MakeMessMap | outputOptionsSheet$MakeModMap){
    
    if(!is.null(factorVars)){ updateRunLog("\nWarning: MESS and MoD maps cannot be generated for models with categorical variables.\n")
      } else {
        
        # order the training data so that we can consider the first and last row only in mess calculations
        train.dat <- select(trainingData, all_of(modVars))
        for(k in 1:nrow(covData)){ train.dat[ ,k] <- sort(train.dat[ ,k]) } 
        ids <- NULL
        for(n in names(train.dat)){ ids[n] <- covariatesSheet$ID[covariatesSheet$CovariateName == n]}
        names(train.dat) <- ids
        
        pred.rng <- rep(NA, nrow(temp))
        names(pred.rng) <- NA
          
        if(any(complete.cases(temp))){
            MessVals <- CalcMESS(rast = data.frame(temp[complete.cases(temp),]), train.dat = train.dat)
            pred.rng[complete.cases(temp)] <- MessVals[ ,2]
            names(pred.rng)[complete.cases(temp)] <- MessVals[ ,1]
        }
        pred.rng <- round(pred.rng,0)
        pred.rng[is.na(pred.rng)] <- -9999  # typeof(pred.rng)
        
        if(outputOptionsSheet$MakeMessMap){ 
        messRaster <- rast(templateRaster, vals = pred.rng)
        if(!is.null(maskValues)){messRaster <- extend(messRaster, maskExt)}
        writeRaster(x = messRaster, 
                    filename = file.path(ssimTempDir, paste0(modType,"_mess_map.tif")), 
                    datatype = "INT4S",
                    NAflag = -9999,
                    overwrite = TRUE)  # is.int(messRaster) 
        
        remove(messRaster)
        progressBar()
        updateRunLog("Finished MESS Map in ", updateBreakpoint())
        }
        if(outputOptionsSheet$MakeModMap){ 
          
          if(is.null(names(pred.rng))){ names(pred.rng) <- NA }
          vals <- as.integer(names(pred.rng))
          vals[is.na(vals)] <- -9999  # typeof(vals)
          
          modRaster <- rast(templateRaster, vals = vals)
          if(!is.null(maskValues)){modRaster <- extend(modRaster, maskExt)}
          writeRaster(x = modRaster, 
                      filename = file.path(ssimTempDir, paste0(modType,"_mod_map.tif")), 
                      datatype = "INT4S", # "INT1U"
                      NAflag = -9999,
                      overwrite = TRUE)  # is.int(modRaster) 
          
          remove(modRaster)
          progressBar()
          updateRunLog("Finished MoD Map in ", updateBreakpoint())
        }
    }
  }
  
    
  # create residuals map ---------------------------------------------------------
    if(outputOptionsSheet$MakeResidualsMap){
      
      if(!outputOptionsSheet$MakeProbabilityMap){
        updateRunLog("\nWarning: Residuals map cannot be generated without generating the Probability map.\n")
      } else {
        
        residSmooth <- readRDS(modelOutputsSheet$ResidualSmoothRDS[i])
  
        predrast <- rast(probRaster)
        if(!is.null(maskValues)){ predrast <- crop(predrast, maskFile)}
        ablock <- 1:(ncol(probRaster) * nrow(probRaster))
        
        p <- xyFromCell(predrast, ablock) 
        p <- na.omit(p)
        
        blockvals <- data.frame(x=p[,1], y=p[,2])
      
        predv <- predict(residSmooth, blockvals)
        predv <- as.numeric(predv)
        predv[!complete.cases(temp)] <- -9999 
        predv[is.na(predv)] <- -9999
        
        residRaster <- rast(templateRaster, vals = predv)
        if(!is.null(maskValues)){residRaster <- extend(residRaster, maskExt)}
        writeRaster(x = residRaster, 
                    filename = file.path(ssimTempDir, paste0(modType,"_resid_map.tif")), 
                    datatype = "FLT8S", 
                    NAflag = -9999,
                    overwrite = TRUE) 
        progressBar()
        updateRunLog("Finished Residuals Map in ", updateBreakpoint())
      }
  }  
  
  # Save output maps -------------------------------------------------------------
    
    # add model Outputs to datasheet
    spatialOutputsSheet <- addRow(spatialOutputsSheet, 
                                  list( # Iteration = 1, Timestep = 0, ModelsID = modelsSheet$ModelsID[modelsSheet$ModelType == modType]))
                                    ModelsID = modelsSheet$ModelName[modelsSheet$ModelType == modType]))
  
    outputRow <- which(spatialOutputsSheet$ModelsID == modelsSheet$ModelName[modelsSheet$ModelType == modType])
                              
    if(outputOptionsSheet$MakeProbabilityMap){
      spatialOutputsSheet$ProbabilityRaster[outputRow] <- file.path(ssimTempDir,paste0(modType,"_prob_map.tif"))
    }
    if(outputOptionsSheet$MakeBinaryMap){
      spatialOutputsSheet$BinaryRaster[outputRow] <- file.path(ssimTempDir,paste0(modType,"_bin_map.tif"))
    }
    if(outputOptionsSheet$MakeMessMap & is.null(factorVars)){
      spatialOutputsSheet$MessRaster[outputRow] <- file.path(ssimTempDir,paste0(modType,"_mess_map.tif"))
    }
    if(outputOptionsSheet$MakeModMap & is.null(factorVars)){
      spatialOutputsSheet$ModRaster[outputRow] <- file.path(ssimTempDir,paste0(modType,"_mod_map.tif"))
    }
    if(outputOptionsSheet$MakeResidualsMap){
      spatialOutputsSheet$ResidualsRaster[outputRow] <- file.path(ssimTempDir,paste0(modType,"_resid_map.tif"))
    }
    
} # end modType loop

# save outputs
saveDatasheet(myScenario, spatialOutputsSheet, "wisdm_OutputSpatial")
progressBar(type = "end")    

