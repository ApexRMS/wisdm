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
source(file.path(packageDir, "0-apply-model-functions.R"))

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
runControlSheet <- datasheet(myScenario, "RunControl", optional = T)
multiprocessingSheet <- datasheet(myScenario, "core_Multiprocessing")
spatialMulitprocessingSheet <- datasheet(myScenario, "corestime_Multiprocessing")
covariateDataSheet <- datasheet(myScenario, "CovariateData", optional = T, lookupsAsFactors = F)
templateSheet <- datasheet(myScenario, "TemplateRaster")
modelOutputsSheet <- datasheet(myScenario, "ModelOutputs", optional = T)
outputOptionsSheet <- datasheet(myScenario, "OutputOptions", optional = T)
spatialOutputsSheet <- datasheet(myScenario, "SpatialOutputs", optional = T)

# Set defaults -----------------------------------------------------------------

## Run control sheet
if(nrow(runControlSheet)<1){
  runControlSheet <- addRow(runControlSheet, list(1,1,0,0))
  saveDatasheet(myScenario, runControlSheet, "RunControl")
}
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


# prep prediction data ---------------------------------------------------------

if(nrow(templateSheet)<1){ stop("ERROR: template raster is missing") }

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
  # typeof(preds)
  
  probRaster <- rast(templateRaster, vals = preds)
  # is.int(probRaster) 
  writeRaster(x = probRaster, 
              filename = file.path(ssimTempDir, paste0(modType,"_prob_map.tif")), 
              datatype = "INT4S",
              overwrite = TRUE)
  
  # writeStart(x = probRaster, 
  #            filename = file.path(ssimTempDir, paste0(modType,"_prob_map.tif")), 
  #            datatype = "INT4S",
  #            overwrite = TRUE) #, gdal=c("COMPRESS=LZW", "TILED=YES", "BIGTIFF=YES"))
  # 
  # writeValues(x = probRaster, v = preds, start = 1, nrows = nrow(templateRaster)) 
  # writeStop(probRaster) 
  
  # system(paste("gdalinfo", probRaster))
}
  
# create MESS and MoD maps -----------------------------------------------------  
  
if(outputOptionsSheet$MakeMessMap){
  
}  
  
# create residuals map ---------------------------------------------------------
  if(outputOptionsSheet$MakeResidualsMap){
    
    # ### Create residual surface of input data
    # 
    #   if(out$dat$split.label != "eval"){
    #     
    #    dev.contrib <- calc.deviance(obs = trainingData$Response, 
    #                                 preds = trainingData$predicted,
    #                                 weights = trainingData$Weight,
    #                                 family = modFamily,
    #                                 return.list = T)$dev.cont
    # 
    #     residual.smooth.fct <- resid.image(dev.contrib = dev.contrib,
    #                                        dat = trainingData, 
    #                                        
    #                                        pred = inlst$train$pred,
    #                                        raw.dat = inlst$train$dat$response, 
    #                                        x = inlst$train$XY$X, 
    #                                        y = inlst$train$XY$Y, 
    #                                        wgt = inlst$train$weight, out$input$model.family, out$input$output.dir, label = out$dat$split.label, out)
    #   } else {
    # 
    #     residual.smooth.fct <- resid.image(calc.dev(inlst$test$dat$response, inlst$test$pred, inlst$test$weight, family = out$input$model.family)$dev.cont, inlst$test$pred,
    #                                        inlst$test$dat$response, inlst$test$XY$X, inlst$test$XY$Y, inlst$test$weight, out$input$model.family, out$input$output.dir, label = out$dat$split.label, out)
    #   }

    Pred.Surface(object = rast(paste0(ssimTempDir,"ProbTiff\\glm_prob_map.tif")),
                 model = residSmooth, # residSmooth=out$mods$auc.output$residual.smooth.fct,
                 filename = paste0(ssimTempDir,"ResidTiff\\glm_resid_map.tif"),
                 NAval = -3000)
  }  


# Save output maps -------------------------------------------------------------
  
  # add model Outputs to datasheet
  spatialOutputsSheet <- addRow(spatialOutputsSheet, 
                              list(Iteration = 1, Timestep = 0,
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
    

  
  # old code --------------------------------
  # create probability/Mess/Mod maps
  # proc.tiff(fit.model = mod,
  #           modType = modType,
  #           mod.vars = covData$CovariatesID,
  #           raster.files = covData$RasterFilePath,
  #           pred.fct = predictFct,
  #           output.options = outputOptionsSheet,
  #           factor.levels = factorVars,
  #           multiprocessing.cores = multiprocessingSheet$MaximumJobs,
  #           temp.directory = ssimTempDir,
  #           train.dat = trainingData, 
  #           tsize = 20,
  #           NAval = -3000)



  
