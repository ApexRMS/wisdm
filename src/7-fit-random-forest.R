## --------------------
## wisdm - fit rf
## ApexRMS, October 2023
## --------------------

# built under R version 4.1.3 & SyncroSim version 2.4.40
# script pulls in pre-processed field, site and covariate data; fits random 
# forest (rf) model; builds model diagnostic and validation plots 

# source dependencies ----------------------------------------------------------
  
  library(rsyncrosim)
  library(tidyr)
  library(dplyr)
  library(randomForest)
  # library(splines)
  
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
  RFSheet <- datasheet(myScenario, "wisdm_RandomForest")
  modelOutputsSheet <- datasheet(myScenario, "wisdm_ModelOutputs", optional = T, empty = T, lookupsAsFactors = F)

  progressBar()

# Error handling ---------------------------------------------------------------

  # check for both presence and absence data
  if(all(fieldDataSheet$Response == 1) | all(fieldDataSheet$Response == 0)){
    stop("Random Forest is a presence-absences (or presence-background) method; please ensure that the Field Data includes both presence and absence (or pseudo-absence) data before continuing.")
  }

#  Set defaults ----------------------------------------------------------------  
  
  # if response column contains only 1's and 0's response = presAbs
  if(max(fieldDataSheet$Response)>1){
    modelFamily <-"poisson" 
  } else { modelFamily <- "binomial" }
  
  ## RF Sheet
  if(nrow(RFSheet)<1){
    RFSheet <- addRow(RFSheet, list(EvaluateCovariateImportance = TRUE,         # importance
                                    CalculateCasewiseImportance = FALSE,        # localImp
                                    NumberOfVariablesSampled = NA,              # mtry
                                    MaximumNodes = NA,                          # maxnodes
                                    NumberOfTrees = 1000,                       # n.trees
                                    NodeSize = NA,                              # nodesize
                                    NormalizeVotes = TRUE,                      # norm.votes
                                    CalculateProximity = FALSE,                 # proximity
                                    SampleWithReplacement = FALSE))             # samp.replace
  }
  
  if(is.na(RFSheet$EvaluateCovariateImportance)){RFSheet$EvaluateCovariateImportance <- TRUE}
  if(is.na(RFSheet$CalculateCasewiseImportance)){RFSheet$CalculateCasewiseImportance <- FALSE}
  if(is.na(RFSheet$NodeSize)){
    if(modelFamily == "poisson"){RFSheet$NodeSize <- 5}
    if(modelFamily == "binomial"){RFSheet$NodeSize <- 1}
  }
  if(is.na(RFSheet$NumberOfTrees)){RFSheet$NumberOfTrees <- 1000}
  if(is.na(RFSheet$NormalizeVotes)){RFSheet$NormalizeVotes <- TRUE}
  if(is.na(RFSheet$CalculateProximity)){RFSheet$CalculateProximity <- FALSE}
  if(is.na(RFSheet$SampleWithReplacement)){RFSheet$SampleWithReplacement <- FALSE}
  
  saveDatasheet(myScenario, RFSheet, "wisdm_RandomForest")
  
  if(!is.na(RFSheet$NumberOfVariablesSampled)){
    if(RFSheet$NumberOfVariablesSampled > nrow(reducedCovariatesSheet)){
      stop(paste0("The number of variables sampled must be between 1 and the total number of variables used in model fitting (in this case ", 
      nrow(reducedCovariatesSheet),"). If left blank, this input will be optimized by the tuneRF funciton."))
    }
  }
  
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
  out$modType <- modType <- "rf"
  
  ## Model options
  out$modOptions <- RFSheet
  out$modOptions$thresholdOptimization <- "Sens=Spec"	# To Do: link to defined Threshold Optimization Method in UI - currently set to default: sensitivity=specificity 
  
  ## Model family 
  # if response column contains only 1's and 0's response = presAbs
  if(max(fieldDataSheet$Response)>1){
     # stop("Random Forest not implemented for count data")
    out$modelFamily <-"poisson" 
  } else { out$modelFamily <- "binomial" }
  
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

  capture.output(cat("Random Forest Results"), file=paste0(ssimTempDir,"\\Data\\", modType, "_output.txt")) 
  on.exit(capture.output(cat("Model Failed\n\n"),file=paste0(ssimTempDir,"\\Data\\", modType, "_output.txt"),append=TRUE))  


# Review model data ------------------------------------------------------------

  if(nrow(trainingData)/(length(out$inputVars)-1)<10){
    updateRunLog(paste("\nYou have approximately ", round(nrow(trainingData)/(ncol(trainingData)-1),digits=1),
                             " observations for every predictor\n consider reducing the number of predictors before continuing\n",sep=""))
  }
  progressBar()
  
# Fit model --------------------------------------------------------------------

  finalMod <- fitModel(dat = trainingData, 
                       out = out)
  
  finalMod$trainingData <- trainingData
  
  # save model to temp storage
  saveRDS(finalMod, file = paste0(ssimTempDir,"\\Data\\", modType, "_model.rds"))
  
  # save output model options
  RFSheet$NumberOfVariablesSampled <- finalMod$mtry
  saveDatasheet(myScenario, RFSheet, "wisdm_RandomForest")
  
  # add relevant model details to out 
  out$finalMod <- finalMod
  out$finalVars <- out$inputVars # random forest doesn't drop variables
  out$nVarsFinal <- length(out$finalVars)
  
  txt0 <- paste("\n\n","Settings:",
                # "\n\trandom seed used                       : ",out$input$seed,
                "\n\tn covariates considered at each split  : ", RFSheet$NumberOfVariablesSampled,
                if(out$pseudoAbs==TRUE) "\n\t   (averaged over each used available split)\n",
                "\n\tn trees                                : ",RFSheet$NumberOfTrees,
                if(out$pseudoAbs==TRUE) "\n\t   (for each used available split)\n",
                sep="")
  txt1 <- "\n\nRelative performance of predictors in final model:\n\n"
  
  modelSummary <- finalMod$importance
  modelSummary <- modelSummary[order(modelSummary[,3],decreasing=T),]
  
  capture.output(cat(txt0),cat(txt1),print(round(modelSummary,4)),file=paste0(ssimTempDir,"\\Data\\", modType, "_output.txt"),append=TRUE) 
  
  updateRunLog("\nSummary of Model:\n")
  coeftbl <- modelSummary
  coeftbl <- round(coeftbl, 6) 
  coeftbl <- cbind(rownames(coeftbl), coeftbl)
  rownames(coeftbl) <- NULL
  colnames(coeftbl)[1] <- "Variable"  
  updateRunLog(pander::pandoc.table.return(coeftbl, style = "simple", split.tables = 100))
  progressBar()

# Test model predictions -------------------------------------------------------
  
  # For the training set for Random Forest take out of bag predictions rather than the regular predictions
  if(out$modelFamily == "poisson"){ out$data$train$predicted <- finalMod$predicted
  } else { out$data$train$predicted <- tweak.p(finalMod$votes[,2]) } # tweak predictions to remove 1/0 so that calc deviance doesn't produce NA/Inf values 
 
  if(validationDataSheet$SplitData){
    out$data$test$predicted <- pred.fct(x=out$data$test, mod=finalMod, modType=modType)
  }
  progressBar()

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
    if("rf_StandardResidualPlots.png" %in% tempFiles){ modelOutputsSheet$ResidualsPlot <- paste0(ssimTempDir,"\\Data\\", modType, "_StandardResidualPlots.png") }
    modelOutputsSheet$ConfusionMatrix <-  paste0(ssimTempDir,"\\Data\\", modType, "_ConfusionMatrix.png")
    modelOutputsSheet$VariableImportancePlot <-  paste0(ssimTempDir,"\\Data\\", modType, "_VariableImportance.png")
    modelOutputsSheet$VariableImportanceData <-  paste0(ssimTempDir,"\\Data\\", modType, "_VariableImportance.csv")
    modelOutputsSheet$ROCAUCPlot <- paste0(ssimTempDir,"\\Data\\", modType, "_ROCAUCPlot.png")
    modelOutputsSheet$CalibrationPlot <- paste0(ssimTempDir,"\\Data\\", modType, "_CalibrationPlot.png")
  } else {
    modelOutputsSheet$ResidualsPlot <- paste0(ssimTempDir,"\\Data\\", modType, "_PoissonResidualPlots.png")
  }
  
  if("rf_AUCPRPlot.png" %in% tempFiles){ modelOutputsSheet$AUCPRPlot <- paste0(ssimTempDir,"\\Data\\", modType, "_AUCPRPlot.png") } 
  
  saveDatasheet(myScenario, modelOutputsSheet, "wisdm_ModelOutputs", append = T)
  progressBar(type = "end")
  