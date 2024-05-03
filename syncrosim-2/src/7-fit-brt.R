## --------------------
## wisdm - fit brt
## ApexRMS, October 2023
## --------------------

# built under R version 4.1.3 & SyncroSim version 2.4.40
# script pulls in pre-processed field, site and covariate data; fits boost 
# regression trees (brt) model; builds model diagnostic and validation plots 

# source dependencies ----------------------------------------------------------
  
  library(rsyncrosim)
  library(tidyr)
  library(dplyr)
  library(dismo)
  
  packageDir <- Sys.getenv("ssim_package_directory")
  source(file.path(packageDir, "00-helper-functions.R"))
  source(file.path(packageDir, "07-fit-model-functions.R"))

# Set progress bar -------------------------------------------------------------

steps <- 11
progressBar(type = "begin", totalSteps = steps)
  
# Connect to library -----------------------------------------------------------

  # Active project and scenario
  myScenario <- scenario()
  # datasheet(myScenario)
  
  # Path to ssim temporary directory
  ssimTempDir <- Sys.getenv("ssim_temp_directory")
  
  # Read in datasheets
  covariatesSheet <- datasheet(myScenario, "wisdm_Covariates", optional = T)
  modelsSheet <- datasheet(myScenario, "wisdm_Models")
  fieldDataSheet <- datasheet(myScenario, "wisdm_FieldData", optional = T)
  validationDataSheet <- datasheet(myScenario, "wisdm_ValidationOptions")
  reducedCovariatesSheet <- datasheet(myScenario, "wisdm_ReducedCovariates", lookupsAsFactors = F)
  siteDataSheet <- datasheet(myScenario, "wisdm_SiteData", lookupsAsFactors = F)
  BRTSheet <- datasheet(myScenario, "wisdm_BRT")
  modelOutputsSheet <- datasheet(myScenario, "wisdm_ModelOutputs", optional = T, empty = T, lookupsAsFactors = F) #, returnInvisible = T

  progressBar()

# Error handling ---------------------------------------------------------------

  # check for both presence and absence data
  if(all(fieldDataSheet$Response == 1) | all(fieldDataSheet$Response == 0)){
    stop("Random Forest is a presence-absences (or presence-background) method; please ensure that the Field Data includes both presence and absence (or pseudo-absence) data before continuing.")
  }

#  Set defaults ----------------------------------------------------------------  
  
  ## BRT Sheet
  if(nrow(BRTSheet)<1){
    fitFromDefaults <- TRUE
    BRTSheet <- addRow(BRTSheet, list(FittingMethod = "Use defaults and tuning",
                                      LearningRate = 0.001,                     # learning.rate
                                      BagFraction = 0.75,                       # bag.fraction
                                      MaximumTrees = 10000,                     # max.trees
                                      NumberOfTrees = 50))                      # n.trees
  } else {
    if(BRTSheet$FittingMethod == "Use defaults and tuning"){ 
      fitFromDefaults <- TRUE
    } else {
      fitFromDefaults <- FALSE
      BRTSheet$FittingMethod <- "Use values provided below"
    }
  }
  
  if(is.na(BRTSheet$LearningRate)){BRTSheet$LearningRate <- 0.001}
  if(is.na(BRTSheet$BagFraction)){BRTSheet$BagFraction <- 0.75}
  if(is.na(BRTSheet$MaximumTrees)){BRTSheet$MaximumTrees <- 10000}
  if(is.na(BRTSheet$NumberOfTrees)){BRTSheet$NumberOfTrees <- 50}
  
  saveDatasheet(myScenario, BRTSheet, "wisdm_BRT")
  
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
  siteDataWide <- siteDataWide[,c('SiteID', reducedCovariatesSheet$CovariatesID)]
  
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
  out$modType <- modType <- "brt"
  
  ## Model options
  out$modOptions <- BRTSheet
  # out$modOptions$stepSize <- out$modOptions$NumberOfTrees
  out$modOptions$thresholdOptimization <- "Sens=Spec"	# To Do: link to defined Threshold Optimization Method in UI - currently set to default: sensitivity=specificity 
  
  ## Model family 
  # if response column contains only 1's and 0's response = presAbs
  if(max(fieldDataSheet$Response)>1){
     # stop("Random Forest not implemented for count data")
    out$modelFamily <-"poisson" 
  } else { out$modelFamily <- "bernoulli"} # "binomial"
  
  ## Candidate variables 
  out$inputVars <- reducedCovariatesSheet$CovariatesID
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
  out$tempDir <- file.path(ssimTempDir, "Data")
  
# Create output text file ------------------------------------------------------

  capture.output(cat("Boosted Regression Tree Results"), file=paste0(ssimTempDir,"\\Data\\", modType, "_output.txt")) 
  on.exit(capture.output(cat("Model Failed\n\n"),file=paste0(ssimTempDir,"\\Data\\", modType, "_output.txt"),append=TRUE))  


# Review model data ------------------------------------------------------------

  if(nrow(trainingData)/(length(out$inputVars)-1)<10){
    updateRunLog(paste("\nYou have approximately ", round(nrow(trainingData)/(ncol(trainingData)-1),digits=1),
                             " observations for every predictor\n consider reducing the number of predictors before continuing\n",sep=""))
  }
  progressBar()
  
# Fit model --------------------------------------------------------------------
  
  if(fitFromDefaults){
    
    # estimate learning rate
    if(validationDataSheet$CrossValidate){
      lr.list <- list()
      for(i in 1:validationDataSheet$NumberOfFolds){
        dat <- trainingData[trainingData$ModelSelectionSplit != i,]
        lr.i <- est.lr(dat = dat, out = out)
        if (is.null(lr.i)){
          stop("Unable to fit model using defaults, parameter tuning, and current data split. Try preparing data with fewer cross validation folds.")
        } else { lr.list[[i]] <-  lr.i}
      }
      BRTSheet$LearningRate <- out$modOptions$LearningRate <- mean(unlist(lapply(lr.list,function(lst){lst$lrs})))
      BRTSheet$NumberOfTrees <- out$modOptions$NumberOfTrees <- mean(unlist(lapply(lr.list,function(lst){lst$n.trees})))
    } else {
      lr.out <- est.lr(dat = trainingData, out = out)
      BRTSheet$LearningRate <- out$modOptions$LearningRate <- lr.out$lrs
      BRTSheet$NumberOfTrees <- out$modOptions$NumberOfTrees <- lr.out$n.trees
    }
  saveDatasheet(myScenario, BRTSheet, "wisdm_BRT")
  }
  
  finalMod <- fitModel(dat = trainingData, 
                       out = out)
  
  if(is.null(finalMod)){
    # updateRunLog("Unable to fit model with defined parameters. Try setting a smaller learning rate or smaller step size (i.e., Number of trees add per stage).")
    stop("Unable to fit model with defined parameters. Try setting a smaller learning rate or smaller step size (i.e., Number of trees added per stage).")
  }
  
  # if(is.null(finalMod)){
  #   repeat{
  #     # out$modOptions$stepSize <- out$modOptions$stepSize-10
  #     # out$modOptions$NumberOfTrees <- out$modOptions$NumberOfTrees-10
  #     out$modOptions$LearningRate <- out$modOptions$LearningRate/2
  #     
  #     finalMod <- fitModel(dat = trainingData,
  #                          out = out)
  #    if(!is.null(finalMod)) break
  #   }
  #   # updateRunLog(paste0("Smaller step size required for model fitting. 'Number of Trees' reduced from ", BRTSheet$NumberOfTrees, " to ", out$modOptions$NumberOfTrees, "."))
  #   # BRTSheet$NumberOfTrees <- out$modOptions$NumberOfTrees
  #   updateRunLog(paste0("Smaller learning rate required for model fitting. Learning rate reduced from ", BRTSheet$LearningRate, " to ", out$modOptions$LearningRate, "."))
  #   BRTSheet$LearningRate <- noquote(format(out$modOptions$LearningRate, scientific = F))
  #   
  #   saveDatasheet(myScenario, BRTSheet, "wisdm_BRT")
  # }
  
  finalMod$trainingData <- trainingData
  
  # add relevant model details to out 
  out$finalMod <- finalMod
  out$finalVars <- finalMod$contributions$var # brt doesn't drop variables
  out$nVarsFinal <- length(out$finalVars)
  
  txt0 <- paste("\n\n","Settings:\n",
                # if(out$pseudoAbs) "(Averaged across available splits)\n", 
                # "\n\trandom seed used             : ",out$input$seed,
                "\n\ttree complexity              : ",finalMod$interaction.depth,
                "\n\tlearning rate                : ",finalMod$shrinkage,
                "\n\tn(trees)                     : ",finalMod$n.trees,
                # "\n\tmodel simplification         : ",simp.method,
                "\n\tn folds                      : ",finalMod$gbm.call$cv.folds,
                "\n\tn covariates in final model  : ",paste(out$finalVars, collapse = ", "),
                sep="")
  txt1 <- "\nRelative influence of predictors in final model:\n\n"
  # txt2 <- "\nImportant interactions in final model:\n\n"
  
  modelSummary <- finalMod$contributions
  modelSummary <- modelSummary[order(modelSummary[,2],decreasing=T),]
  row.names(modelSummary) <- NULL
  
  capture.output(cat(txt0),cat(txt1),print(modelSummary),file=paste0(ssimTempDir,"\\Data\\", modType, "_output.txt"),append=TRUE) 
  
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
  saveRDS(finalMod, file = paste0(ssimTempDir,"\\Data\\", modType, "_model.rds"))
  
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

  tempFiles <- list.files(file.path(ssimTempDir, "Data"))
  
  # add model Outputs to datasheet
  modelOutputsSheet <- addRow(modelOutputsSheet, 
                              list(ModelsID = modelsSheet$ModelName[modelsSheet$ModelType == modType],
                                   ModelRDS = paste0(ssimTempDir,"\\Data\\", modType, "_model.rds"),
                                   ResponseCurves = paste0(ssimTempDir,"\\Data\\", modType, "_ResponseCurves.png"),
                                   TextOutput = paste0(ssimTempDir,"\\Data\\", modType, "_output.txt"),
                                   ResidualSmoothPlot = paste0(ssimTempDir,"\\Data\\", modType, "_ResidualSmoothPlot.png"),
                                   ResidualSmoothRDS = paste0(ssimTempDir,"\\Data\\", modType, "_ResidualSmoothFunction.rds")))
  
  
  if(out$modelFamily != "poisson"){
    if("brt_StandardResidualPlots.png" %in% tempFiles){ modelOutputsSheet$ResidualsPlot <- paste0(ssimTempDir,"\\Data\\", modType, "_StandardResidualPlots.png") }
    modelOutputsSheet$ConfusionMatrix <-  paste0(ssimTempDir,"\\Data\\", modType, "_ConfusionMatrix.png")
    modelOutputsSheet$VariableImportancePlot <-  paste0(ssimTempDir,"\\Data\\", modType, "_VariableImportance.png")
    modelOutputsSheet$VariableImportanceData <-  paste0(ssimTempDir,"\\Data\\", modType, "_VariableImportance.csv")
    modelOutputsSheet$ROCAUCPlot <- paste0(ssimTempDir,"\\Data\\", modType, "_ROCAUCPlot.png")
    modelOutputsSheet$CalibrationPlot <- paste0(ssimTempDir,"\\Data\\", modType, "_CalibrationPlot.png")
  } else {
    modelOutputsSheet$ResidualsPlot <- paste0(ssimTempDir,"\\Data\\", modType, "_PoissonResidualPlots.png")
  }
  
  if("brt_AUCPRPlot.png" %in% tempFiles){ modelOutputsSheet$AUCPRPlot <- paste0(ssimTempDir,"\\Data\\", modType, "_AUCPRPlot.png") } 
  
  saveDatasheet(myScenario, modelOutputsSheet, "wisdm_ModelOutputs", append = T)
  progressBar(type = "end")
  