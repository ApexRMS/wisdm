## --------------------
## wisdm - fit brt
## ApexRMS, June 2025
## --------------------

# built under R version 4.1.3, SyncroSim 3.1.10 & rsyncrosim 2.1.3
# script pulls in pre-processed field, site and covariate data; fits boost 
# regression trees (brt) model; builds model diagnostic and validation plots 

# source dependencies ----------------------------------------------------------
  
  library(rsyncrosim) # install.packages("C:/GitHub/rsyncrosim", type="source", repos=NULL) 
  library(tidyr)
  library(dplyr)
  library(dismo)
  
  packageDir <- Sys.getenv("ssim_package_directory")
  source(file.path(packageDir, "00-helper-functions.R"))
  source(file.path(packageDir, "08-fit-model-functions.R"))

# Set progress bar -------------------------------------------------------------

steps <- 11
updateRunLog('8 - Boosted Regression Trees => Begin')
progressBar(type = "begin", totalSteps = steps)
  
# Connect to library -----------------------------------------------------------

  # Active project and scenario
  myScenario <- scenario()
  # datasheet(myScenario)
  
  # Path to ssim temporary directory
  ssimTempDir <- ssimEnvironment()$TransferDirectory
  
  # Read in datasheets
  covariatesSheet <- datasheet(myScenario, "wisdm_Covariates", optional = T)
  modelsSheet <- datasheet(myScenario, "wisdm_Models")
  fieldDataSheet <- datasheet(myScenario, "wisdm_FieldData", optional = T)
  validationDataSheet <- datasheet(myScenario, "wisdm_ValidationOptions")
  retainedCovariatesSheet <- datasheet(myScenario, "wisdm_RetainedCovariates", lookupsAsFactors = F)
  siteDataSheet <- datasheet(myScenario, "wisdm_SiteData", lookupsAsFactors = F)
  BRTSheet <- datasheet(myScenario, "wisdm_BRT")
  modelOutputsSheet <- datasheet(myScenario, "wisdm_OutputModel", optional = T, returnInvisible = T, empty = T, lookupsAsFactors = F) %>% drop_na()

  progressBar()

# Error handling ---------------------------------------------------------------

  # check for both presence and absence data
  if(all(fieldDataSheet$Response == 1) | all(fieldDataSheet$Response == 0)){
    stop("Random Forest is a presence-absences (or presence-background) method; please ensure that the Field Data includes both presence and absence (or pseudo-absence) data before continuing.")
  }

#  Set defaults ----------------------------------------------------------------  
  
  ## BRT Sheet
  if(nrow(BRTSheet)<1 | all(is.na(BRTSheet))){
    fitFromDefaults <- TRUE
    BRTSheet <- BRTSheet %>% drop_na()
    BRTSheet <- bind_rows(BRTSheet, list(FittingMethod = "Use defaults and tuning",
                                      LearningRate = 0.001,                     # learning.rate
                                      BagFraction = 0.75,                       # bag.fraction
                                      MaximumTrees = 10000,                     # max.trees
                                      NumberOfTrees = 50))                      # n.trees
  } else {
    if(BRTSheet$FittingMethod == "Use defaults and tuning"){ 
      fitFromDefaults <- TRUE
      BRTSheet$LearningRate <- 0.001
      BRTSheet$BagFraction <- 0.75
      BRTSheet$MaximumTrees <- 10000
      BRTSheet$NumberOfTrees <- 50
    } else {
      fitFromDefaults <- FALSE
      BRTSheet$FittingMethod <- "Use values provided below"
    }
  }
  
  if(is.na(BRTSheet$LearningRate)){BRTSheet$LearningRate <- 0.001}
  if(is.na(BRTSheet$BagFraction)){BRTSheet$BagFraction <- 0.75}
  if(is.na(BRTSheet$MaximumTrees)){BRTSheet$MaximumTrees <- 10000}
  if(is.na(BRTSheet$NumberOfTrees)){BRTSheet$NumberOfTrees <- 50}
  
  ## Validation Sheet
  if(nrow(validationDataSheet)<1){
    validationDataSheet <- bind_rows(validationDataSheet, list(SplitData = FALSE,
                                                            CrossValidate = FALSE))
  }
  if(is.na(validationDataSheet$CrossValidate)){validationDataSheet$CrossValidate <- FALSE}
  if(is.na(validationDataSheet$SplitData)){validationDataSheet$SplitData <- FALSE}
  progressBar()

# Prep data for model fitting --------------------------------------------------

  siteDataWide <- spread(siteDataSheet, key = CovariatesID, value = "Value")
  rm(siteDataSheet); gc()
  
  # remove variables dropped due to correlation
  siteDataWide <- siteDataWide[,c('SiteID', retainedCovariatesSheet$CovariatesID)]
  
  # merge field and site data
  siteDataWide <- merge(fieldDataSheet, siteDataWide, by = "SiteID")
  rm(fieldDataSheet); gc()
  
  # remove sites with incomplete data 
  allCases <- nrow(siteDataWide)
  siteDataWide <- siteDataWide[complete.cases(subset(siteDataWide, select = c(-UseInModelEvaluation, -ModelSelectionSplit, -Weight))),]
  compCases <- nrow(siteDataWide)
  if(compCases/allCases < 0.9){updateRunLog(paste("\nWarning: ", round((1-compCases/allCases)*100,digits=2),"% of cases were removed because of missing values.\n",sep=""))}
  
  # set site weights to default of 1 if not already supplied
  if(all(is.na(siteDataWide$Weight))){siteDataWide$Weight <- 1}
  
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
  rm(siteDataWide); gc()
  trainingData <- trainTestDatasets$`FALSE`
  if(!validationDataSheet$CrossValidate){trainingData$ModelSelectionSplit <- FALSE}
  testingData <- trainTestDatasets$`TRUE`
  rm(trainTestDatasets); gc()
  if(!is.null(testingData)){testingData$ModelSelectionSplit <- FALSE}
  progressBar()
  

# Model definitions ------------------------------------------------------------

  # create object to store intermediate model selection/evaluation inputs
  out <- list()
  
  ## Model type
  out$modType <- modType <- "brt"
  
  ## Model options
  out$modOptions <- BRTSheet
  # out$modOptions$stepSize <- out$modOptions$NumberOfTrees
  out$modOptions$thresholdOptimization <- "Sens=Spec"	# To Do: link to defined Threshold Optimization Method in UI - currently set to default: sensitivity=specificity 
  
  ## Model family 
  out$modelFamily <- "bernoulli" # "binomial"
  
  ## Candidate variables 
  out$inputVars <- retainedCovariatesSheet$CovariatesID
  out$factorInputVars <- factorInputVars
  
  ## training data
  out$data$train <- trainingData
 
  ## testing data
  out$data$test <- testingData
  
  ## pseudo absence  
  out$pseudoAbs <- pseudoAbs
  
  ## Validation options
  out$validationOptions <- validationDataSheet 
  
  ## path to temp ssim storage 
  out$tempDir <- ssimTempDir
  
# Create output text file ------------------------------------------------------

  capture.output(cat("Boosted Regression Tree Results"), file=file.path(ssimTempDir, paste0(modType, "_output.txt"))) 
  on.exit(capture.output(cat("Model Failed\n\n"),file=file.path(ssimTempDir, paste0(modType, "_output.txt"),append=TRUE)))  

# Review model data ------------------------------------------------------------

  if(nrow(trainingData)/(length(out$inputVars)-1)<10){
    updateRunLog(paste("\nYou have approximately ", round(nrow(trainingData)/(ncol(trainingData)-1),digits=1),
                             " observations for every predictor\n consider reducing the number of predictors before continuing\n",sep=""))
  }
  progressBar()
  
# Fit model --------------------------------------------------------------------
  
  if(fitFromDefaults){
      # tune learning rate to give approx 1000 trees in final model
      allLrOut <- est.lr(dat = trainingData, out = out)
      lrOut <- allLrOut[allLrOut$n.trees != allLrOut$trees.fit,]
      BRTSheet$LearningRate <- out$modOptions$LearningRate <- lrOut$lrs[1]
    }

  saveDatasheet(myScenario, BRTSheet, "wisdm_BRT")
  
  finalMod <- fitModel(dat = trainingData, 
                       out = out)

  if(is.null(finalMod)){
    # updateRunLog("Unable to fit model with defined parameters. Try setting a smaller learning rate or smaller step size (i.e., Number of trees add per stage).")
    stop("Unable to fit model with defined parameters. Try setting a smaller learning rate or smaller step size (i.e., Number of trees added per stage).")
  }
  
  # add relevant model details to out
  out$finalMod <- finalMod
  out$finalVars <- finalMod$contributions$var # brt doesn't drop variables
  out$nVarsFinal <- length(out$finalVars)

  # number of trees
  nTrees <- if (!is.null(finalMod$gbm.call$best.trees)) { # gbm.step stores optimal number of trees here
    finalMod$gbm.call$best.trees
  } else if (!is.null(finalMod$n.trees)) { # gbm stores the maximum number fitted
    finalMod$n.trees
  } else { NA }
  
  # number of folds
  cvFolds <- if (!is.null(finalMod$gbm.call$cv.folds)) { # gbm.step
    finalMod$gbm.call$cv.folds
  } else if (!is.null(finalMod$cv.folds)) { # gbm
    finalMod$cv.folds
  } else { NA }
  
  txt0 <- paste("\n\n","Settings:\n",              
                # "\n\trandom seed used             : ",out$input$seed,
                "\n\ttree complexity              : ",finalMod$interaction.depth,
                "\n\tlearning rate                : ",finalMod$shrinkage,
                "\n\tn(trees)                     : ",nTrees,
                "\n\tn folds                      : ",cvFolds,
                "\n\tn covariates in final model  : ",paste(out$finalVars, collapse = ", "),
                sep="")
  
  txt1 <- "\nRelative influence of predictors in final model:\n\n"
  # txt2 <- "\nImportant interactions in final model:\n\n"
  
  modelSummary <- finalMod$contributions
  modelSummary <- modelSummary[order(modelSummary[,2],decreasing=T),]
  row.names(modelSummary) <- NULL
  
  capture.output(cat(txt0),cat(txt1),print(modelSummary),file=file.path(ssimTempDir, paste0(modType, "_output.txt")),append=TRUE) 
  
  updateRunLog("\nSummary of Model:\n")
  coeftbl <- modelSummary
  coeftbl[,2] <- round(coeftbl[,2], 4) 
  colnames(coeftbl) <- c("Variable", "Relative Influence")  
  updateRunLog(pander::pandoc.table.return(coeftbl, style = "simple", split.tables = 100))
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
  finalMod$trainingData <- trainingData
  saveRDS(finalMod, file = file.path(ssimTempDir, paste0(modType, "_model.rds")))
  
## Run Cross Validation (if specified) -----------------------------------------
  
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
  modelOutputsSheet <- bind_rows(modelOutputsSheet, 
                              list(ModelsID = modelsSheet$ModelName[modelsSheet$ModelType == modType],
                                   ModelRDS = file.path(ssimTempDir, paste0(modType, "_model.rds")),
                                   ResponseCurves = file.path(ssimTempDir, paste0(modType, "_ResponseCurves.png")),
                                   TextOutput = file.path(ssimTempDir, paste0(modType, "_output.txt")),
                                   ResidualSmoothPlot = file.path(ssimTempDir, paste0(modType, "_ResidualSmoothPlot.png")),
                                   ResidualSmoothRDS = file.path(ssimTempDir, paste0(modType, "_ResidualSmoothFunction.rds")),
                                   ConfusionMatrix = file.path(ssimTempDir, paste0(modType, "_ConfusionMatrix.png")),
                                   VariableImportancePlot =  file.path(ssimTempDir, paste0(modType, "_VariableImportance.png")),
                                   VariableImportanceData =  file.path(ssimTempDir, paste0(modType, "_VariableImportance.csv")),
                                   ROCAUCPlot = file.path(ssimTempDir, paste0(modType, "_ROCAUCPlot.png")),
                                   CalibrationPlot = file.path(ssimTempDir, paste0(modType, "_CalibrationPlot.png"))))

  if("brt_StandardResidualPlots.png" %in% tempFiles){ modelOutputsSheet$ResidualsPlot <- file.path(ssimTempDir, paste0(modType, "_StandardResidualPlots.png")) }
  if("brt_AUCPRPlot.png" %in% tempFiles){ modelOutputsSheet$AUCPRPlot <- file.path(ssimTempDir, paste0(modType, "_AUCPRPlot.png")) } 
  
  saveDatasheet(myScenario, modelOutputsSheet, "wisdm_OutputModel", append = T)
  progressBar(type = "end")
  