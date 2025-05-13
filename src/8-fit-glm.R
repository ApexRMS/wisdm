## --------------------
## wisdm - fit glm
## ApexRMS, March 2024
## --------------------

# built under R version 4.1.3 & SyncroSim version 3.0.0
# script pulls in pre-processed field, site and covariate data; fits glm; builds
# model diagnostic and validation plots 

# source dependencies ----------------------------------------------------------
  
  library(rsyncrosim)
  library(tidyr)
  library(dplyr)

  packageDir <- Sys.getenv("ssim_package_directory")
  source(file.path(packageDir, "00-helper-functions.R"))
  source(file.path(packageDir, "08-fit-model-functions.R"))

# set progress bar -------------------------------------------------------------

steps <- 11
updateRunLog('8 - Generalized Linear Model => Begin')
progressBar(type = "begin", totalSteps = steps)

# Connect to library -----------------------------------------------------------

  # Active project and scenario
  myScenario <- scenario()
  # datasheet(myScenario)
  
  # Path to ssim temporary directory
  ssimTempDir <- ssimEnvironment()$TransferDirectory
  # ssimTempDir <- Sys.getenv("ssim_temp_directory")
  
  # Read in datasheets
  covariatesSheet <- datasheet(myScenario, "wisdm_Covariates", optional = T)
  modelsSheet <- datasheet(myScenario, "wisdm_Models")
  fieldDataSheet <- datasheet(myScenario, "wisdm_FieldData", optional = T)
  validationDataSheet <- datasheet(myScenario, "wisdm_ValidationOptions")
  retainedCovariatesSheet <- datasheet(myScenario, "wisdm_RetainedCovariates", lookupsAsFactors = F)
  siteDataSheet <- datasheet(myScenario, "wisdm_SiteData", lookupsAsFactors = F)
  GLMSheet <- datasheet(myScenario, "wisdm_GLM")
  modelOutputsSheet <- datasheet(myScenario, "wisdm_OutputModel", optional = T, returnInvisible = T, empty = T, lookupsAsFactors = F) %>% drop_na()

  progressBar()

# Error handling ---------------------------------------------------------------

  # check for both presence and absence data
  if(all(fieldDataSheet$Response == 1) | all(fieldDataSheet$Response == 0)){
    stop("GLM is a presence-absences method; please ensure that the Field Data includes both presence and absence (or pseudo-absence) data before continuing.")
  }
 
#  Set defaults ----------------------------------------------------------------  
  
  ## GLM Sheet
  if(nrow(GLMSheet)<1){
    GLMSheet <- addRow(GLMSheet, list(SelectBestPredictors = FALSE,
                                      SimplificationMethod = "AIC",
                                      ConsiderSquaredTerms = FALSE,
                                      ConsiderInteractions = FALSE))
  }
  if(is.na(GLMSheet$SelectBestPredictors)){GLMSheet$SelectBestPredictors <- FALSE}
  if(is.na(GLMSheet$SimplificationMethod)){GLMSheet$SimplificationMethod <- "AIC"}
  if(is.na(GLMSheet$ConsiderSquaredTerms)){GLMSheet$ConsiderSquaredTerms <- FALSE}
  if(is.na(GLMSheet$ConsiderInteractions)){GLMSheet$ConsiderInteractions <- FALSE}
  
  saveDatasheet(myScenario, GLMSheet, "wisdm_GLM")
  
  ## Validation Sheet
  if(nrow(validationDataSheet)<1){
    validationDataSheet <- addRow(validationDataSheet, list(SplitData = FALSE,
                                                            CrossValidate = FALSE))
  }
  if(is.na(validationDataSheet$CrossValidate)){validationDataSheet$CrossValidate <- FALSE}
  if(is.na(validationDataSheet$SplitData)){validationDataSheet$SplitData <- FALSE}
  progressBar()

# Prep data for model fitting --------------------------------------------------

  siteDataWide <- spread(siteDataSheet, key = CovariatesID, value = "Value")
  
  # remove variables dropped due to correlation
  siteDataWide <- siteDataWide[,c('SiteID', retainedCovariatesSheet$CovariatesID)]
  
  # merge field and site data
  siteDataWide <- merge(fieldDataSheet, siteDataWide, by = "SiteID")
  
  # remove sites with incomplete data 
  allCases <- nrow(siteDataWide)
  siteDataWide <- siteDataWide[complete.cases(subset(siteDataWide, select = c(-UseInModelEvaluation, -ModelSelectionSplit, -Weight))),]
  compCases <- nrow(siteDataWide)
  if(compCases/allCases < 0.9){updateRunLog(paste("\nWarning: ", round((1-compCases/allCases)*100,digits=2),"% of cases were removed because of missing values.\n",sep=""))}
  
  # set site weights to default of 1 if not already supplied
  if(all(is.na(siteDataWide$Weight))){siteDataWide$Weight <- 1}
  
  # ignore background data if present
  # siteDataWide <- siteDataWide[!siteDataWide$Response == -9999,]
  
  # set pseudo absences to zero 
  if(any(siteDataWide$Response == -9998)){pseudoAbs <- TRUE} else {pseudoAbs <- FALSE}
  siteDataWide$Response[siteDataWide$Response == -9998] <- 0
  
  # set categorical variable to factor
  factorInputVars <- covariatesSheet$CovariateName[which(covariatesSheet$IsCategorical == T & covariatesSheet$CovariateName %in% names(siteDataWide))]
  if(length(factorInputVars)>0){ 
    for (i in factorInputVars){
      siteDataWide[,i] <- factor(siteDataWide[,i])
    }
  }
 
  # identify training and testing sites 
  trainTestDatasets <- split(siteDataWide, f = siteDataWide[,"UseInModelEvaluation"], drop = T)
  trainingData <- trainTestDatasets$`FALSE`
  if(!validationDataSheet$CrossValidate){trainingData$ModelSelectionSplit <- FALSE}
  testingData <- trainTestDatasets$`TRUE`
  progressBar()

 # Model definitions ------------------------------------------------------------

  # create object to store intermediate model selection/evaluation inputs
  out <- list()
  
  ## Model type
  out$modType <- modType <- "glm"
  
  ## Model options
  out$modOptions <- GLMSheet
  out$modOptions$thresholdOptimization <- "Sens=Spec"	# To Do: link to defined Threshold Optimization Method in UI - currently set to default: sensitivity=specificity 
  
  ## Model family 
  # if response column contains only 1's and 0's response = presAbs
  if(max(fieldDataSheet$Response)>1){
     out$modelFamily <-"poisson" 
  } else { out$modelFamily <- "binomial" }
  
  ## Candidate variables 
  out$inputVars <- retainedCovariatesSheet$CovariatesID
  out$factorInputVars <- factorInputVars
  
  ## training data
  out$data$train <- trainingData
  
  ## testing data
  if(!is.null(testingData)){
    testingData$ModelSelectionSplit <- FALSE
  }
  out$data$test <- testingData
  
  ## pseudo absence  
  out$pseudoAbs <- pseudoAbs
  
  ## Validation options
  out$validationOptions <- validationDataSheet 
  
  ## path to temp ssim storage 
  out$tempDir <- ssimTempDir
  
# Create output text file ------------------------------------------------------

  capture.output(cat("Generalized Linear Model Results"), file=file.path(ssimTempDir,paste0(modType, "_output.txt"))) 
  on.exit(capture.output(cat("Model Failed\n\n"),file=file.path(ssimTempDir,paste0(modType, "_output.txt")),append=TRUE))  


# Review model data ------------------------------------------------------------

  if(nrow(trainingData)/(length(out$inputVars)-1)<10){
    updateRunLog(paste("\nYou have approximately ", round(nrow(trainingData)/(ncol(trainingData)-1),digits=1),
                             " observations for every predictor\n consider reducing the number of predictors before continuing\n",sep=""))
  }

  progressBar()

# Fit model --------------------------------------------------------------------

  finalMod <- fitModel(dat = trainingData, 
                      out = out)
  
  # save model to temp storage
  # saveRDS(finalMod, file = paste0(ssimTempDir,"\\Data\\", modType, "_model.rds"))
  
  # add relevant model details to out 
  out$finalMod <- finalMod
  out$finalVars <- attr(terms(formula(finalMod)),"term.labels")
  # have to remove all the junk with powers and interactions for mess map production to work
  out$finalVars <- unique(unlist(strsplit(gsub("I\\(","",gsub("\\^2)","",out$finalVars)),":")))
  out$nVarsFinal <- length(out$finalVars)
  
  # add relevant model details to text output 
  txt0 <- paste("\n\n","Settings:\n","\n\t model family:  ",out$modelFamily,
              "\n\t simplification method:  ",GLMSheet$SimplificationMethod,
              "\n\n\n","Results:\n\t ","number covariates in final model:  ",length(attr(terms(formula(finalMod)),"term.labels")),"\n",sep="")
  
  modSummary <- summary(finalMod)
  
  updateRunLog("\nSummary of Model:\n")
  coeftbl <- modSummary$coefficients
  coeftbl <- round(coeftbl, 6) 
  coeftbl <- cbind(rownames(coeftbl), coeftbl)
  rownames(coeftbl) <- NULL
  colnames(coeftbl) <- c("Variable", "Estimate", "Std. Error", "z value", "Pr(>|z|)")  
  updateRunLog(pander::pandoc.table.return(coeftbl, style = "simple", split.tables = 100))
  
  capture.output(cat(txt0),modSummary,file=file.path(ssimTempDir,paste0(modType, "_output.txt")), append=TRUE)
  
  if(length(coef(finalMod))==1) stop("Null model was selected. \nEvaluation metrics and plots will not be produced") 
  progressBar()

# Test model predictions -------------------------------------------------------
  
  out$data$train$predicted <- pred.fct(x=out$data$train, mod=finalMod, modType=modType)
  
  if(validationDataSheet$SplitData){
    out$data$test$predicted <- pred.fct(x=out$data$test, mod=finalMod, modType=modType)
  }
  progressBar()

# Evaluate thresholds (for use with binary output) ----------------------------
  
  predOcc <- out$data$train[out$data$train$Response >= 1 , "predicted"]
  predAbs <- out$data$train[out$data$train$Response == 0 , "predicted"]
  
  evalOut <- modelEvaluation(predOcc = predOcc, predAbs = predAbs)
  
  finalMod$binThresholds <- thresholds <- evalOut@thresholds
  
  names(thresholds) <- c("Max kappa", "Max sensitivity and specificity", "No omission", 
                         "Prevalence", "Sensitivity equals specificity")
  
  updateRunLog("\nThresholds:\n")
  tbl <- round(thresholds, 6) 
  updateRunLog(pander::pandoc.table.return(tbl, style = "simple", split.tables = 100))
  
  # save model info to temp storage
  saveRDS(finalMod, file = file.path(ssimTempDir, paste0(modType, "_model.rds")))
  
# Run Cross Validation (if specified) ------------------------------------------
  
  if(validationDataSheet$CrossValidate){
    
    out <- cv.fct(out = out,
                  nfolds = validationDataSheet$NumberOfFolds)
  }
  progressBar()

# Generate Model Outputs -------------------------------------------------------
 
  ## AUC/ROC - Residual Plots - Variable Importance -  Calibration - Confusion Matrix ##
  
  out <- suppressWarnings(makeModelEvalPlots(out=out))
  progressBar()
  ## Response Curves ##
  
  response.curves(out)
  progressBar()

# Save model outputs -----------------------------------------------------------

  tempFiles <- list.files(ssimTempDir)
  
  # add model Outputs to datasheet
  modelOutputsSheet <- addRow(modelOutputsSheet, 
                              list(ModelsID = modelsSheet$ModelName[modelsSheet$ModelType == modType],
                                   ModelRDS = file.path(ssimTempDir, paste0(modType, "_model.rds")),
                                   ResponseCurves = file.path(ssimTempDir, paste0(modType, "_ResponseCurves.png")),
                                   TextOutput = file.path(ssimTempDir, paste0(modType, "_output.txt")),
                                   ResidualSmoothPlot = file.path(ssimTempDir, paste0(modType, "_ResidualSmoothPlot.png")),
                                   ResidualSmoothRDS = file.path(ssimTempDir, paste0(modType, "_ResidualSmoothFunction.rds"))))
  
  
  if(out$modelFamily != "poisson"){
    if("glm_StandardResidualPlots.png" %in% tempFiles){ modelOutputsSheet$ResidualsPlot <- file.path(ssimTempDir, paste0(modType, "_StandardResidualPlots.png")) }
    modelOutputsSheet$ConfusionMatrix <- file.path(ssimTempDir, paste0(modType, "_ConfusionMatrix.png"))
    modelOutputsSheet$VariableImportancePlot <- file.path(ssimTempDir, paste0(modType, "_VariableImportance.png"))
    modelOutputsSheet$VariableImportanceData <- file.path(ssimTempDir, paste0(modType, "_VariableImportance.csv"))
    modelOutputsSheet$ROCAUCPlot <- file.path(ssimTempDir, paste0(modType, "_ROCAUCPlot.png"))
    modelOutputsSheet$CalibrationPlot <- file.path(ssimTempDir, paste0(modType, "_CalibrationPlot.png"))
  } else {
    modelOutputsSheet$ResidualsPlot <- file.path(ssimTempDir, paste0(modType, "_PoissonResidualPlots.png"))
  }
  
  if("glm_AUCPRPlot.png" %in% tempFiles){ modelOutputsSheet$AUCPRPlot <- file.path(ssimTempDir, paste0(modType, "_AUCPRPlot.png")) } 
  
  saveDatasheet(myScenario, modelOutputsSheet, "wisdm_OutputModel", append = T)
  progressBar(type = "end")
  