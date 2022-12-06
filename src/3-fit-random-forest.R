## --------------------
## wisdm - fit rf
## ApexRMS, December 2022
## --------------------

# built under R version 4.1.3 & SyncroSim version 2.4.0
# script pulls in pre-processed field, site and covariate data; fits random 
# forest (rf) model; builds model diagnostic and validation plots 

# source dependencies ----------------------------------------------------------
  
  library(rsyncrosim)
  library(tidyr)
  library(dplyr)
  library(randomForest)
  library(PresenceAbsence)
  library(PRROC)
  library(ROCR)
  library(ggplot2)
  library(splines)
  
  packageDir <- Sys.getenv("ssim_package_directory")
  source(file.path(packageDir, "00-helper-functions.R"))
  source(file.path(packageDir, "03-fit-model-functions.R"))
  
# Connect to library -----------------------------------------------------------

  # Active project and scenario
  myLibrary <- ssimLibrary()
  myProject <- rsyncrosim::project()
  myScenario <- scenario()
  # datasheet(myScenario)
  
  # Path to ssim temporary directory
  ssimTempDir <- Sys.getenv("ssim_temp_directory")
  
  # Read in datasheets
  covariatesSheet <- datasheet(myProject, "wisdm_Covariates", optional = T)
  fieldDataSheet <- datasheet(myScenario, "wisdm_FieldData", optional = T)
  ValidationDataSheet <- datasheet(myScenario, "wisdm_ValidationOptions")
  reducedCovariatesSheet <- datasheet(myScenario, "wisdm_ReducedCovariates", lookupsAsFactors = F)
  siteDataSheet <- datasheet(myScenario, "wisdm_SiteData", lookupsAsFactors = F)
  RFSheet <- datasheet(myScenario, "wisdm_RandomForest")
  modelOutputsSheet <- datasheet(myScenario, "wisdm_ModelOutputs", optional = T)

  
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
                                    CalculateProximity = FALSE,                    # proximity
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
  if(nrow(ValidationDataSheet)<1){
    ValidationDataSheet <- addRow(ValidationDataSheet, list(SplitData = FALSE,
                                                            CrossValidate = FALSE))
  }
  if(is.na(ValidationDataSheet$CrossValidate)){ValidationDataSheet$CrossValidate <- FALSE}
  if(is.na(ValidationDataSheet$SplitData)){ValidationDataSheet$SplitData <- FALSE}

  
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
  siteDataWide <- siteDataWide[!siteDataWide$Response == -9999,]
  
  # set categorical variable to factor
  factorInputVars <- covariatesSheet$CovariateName[which(covariatesSheet$IsCategorical == T & covariatesSheet$CovariateName %in% names(siteDataWide))]
  if(length(factorInputVars)>0){ 
    for (i in factorInputVars){
      siteDataWide[,i] <- factor(siteDataWide[,i])
    }
  }
 
  # identify training and testing sites 
  # siteDataWide$UseInModelEvaluation[is.na(siteDataWide$UseInModelEvaluation)] <- FALSE
  trainTestDatasets <- split(siteDataWide, f = siteDataWide[,"UseInModelEvaluation"], drop = T)
  trainingData <- trainTestDatasets$`FALSE`
  testingData <- trainTestDatasets$`TRUE`
  
  # identify cross validation folds for model selection (if specified)
  if(ValidationDataSheet$CrossValidate){ 
    CVSplitsTrain <- list()
    CVSplitsTest <- list()
    for (i in 1:(ValidationDataSheet$NumberOfFolds)){
      folds <- 1:ValidationDataSheet$NumberOfFolds
      tfolds <- folds[-i]
      CVSplitsTrain[[i]] <- trainingData[trainingData$ModelSelectionSplit %in% tfolds,]
      CVSplitsTest[[i]] <- trainingData[trainingData$ModelSelectionSplit == i,]
    }
  } else {
    trainingData$ModelSelectionSplit <- FALSE
  }
  

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
  
  # CV splits
  if(ValidationDataSheet$CrossValidate){
    out$data$cvSplits$train <- CVSplitsTrain
    out$data$cvSplits$test <- CVSplitsTest
  }
  
  ## pseudo absence  
  out$pseudoAbs <- FALSE # To Do: build out call to pseudo absence 
  
  ## Validation options
  out$validationOptions <- ValidationDataSheet 
  
  ## path to temp ssim storage 
  out$tempDir <- paste0(ssimTempDir,"\\Data\\")
  
# Create output text file ------------------------------------------------------

  capture.output(cat("Random Forest Results"), file=paste0(ssimTempDir,"\\Data\\rf_output.txt")) 
  on.exit(capture.output(cat("Model Failed\n\n"),file=paste0(ssimTempDir,"\\Data\\rf_output.txt"),append=TRUE))  


# Review model data ------------------------------------------------------------

  if(nrow(trainingData)/(length(out$inputVars)-1)<10){
    updateRunLog(paste("\nYou have approximately ", round(nrow(trainingData)/(ncol(trainingData)-1),digits=1),
                             " observations for every predictor\n consider reducing the number of predictors before continuing\n",sep=""))
  }

  
# Fit model --------------------------------------------------------------------

  # source(file.path(packageDir, "fitModel.R")) ## GLM - code sourced and updated from model.fit.r 

  out <- fitModel(dat = trainingData, 
                  out = out,
                  modType = modType,
                  full.fit = T)
  
  finalMod <- out$finalMod
  finalMod[["trainingData"]] <- out$data$train
  
  # store output model option
  outModOptions <- out$modOptions
  outModOptions$thresholdOptimization <- NULL
  saveDatasheet(myScenario, outModOptions, "wisdm_RandomForest")
  
  out$modOptions$NumberOfVariablesSampled <- NA
  
  # save model to temp storage
  saveRDS(finalMod, file = paste0(ssimTempDir,"\\Data\\rf_model.rds"))
  finalMod[["trainingData"]] <- NULL
  
  # add relevant model details to out 
  out$finalVars <- out$inputVars # random forest doesn't drop variables
  out$nVarsFinal <- length(out$finalVars)
  
  txt0 <- paste("\n\n","Settings:",
                # "\n\trandom seed used                       : ",out$input$seed,
                "\n\tn covariates considered at each split  : ", mean(unlist(lapply(out$finalMod,"[",14))),
                if(out$pseudoAbs==TRUE) "\n\t   (averaged over each used available split)\n",
                "\n\tn trees                                : ",out$modOptions$NumberOfTrees,
                if(out$pseudoAbs==TRUE) "\n\t   (for each used available split)\n",
                sep="")
  txt1 <- "\n\nRelative performance of predictors in final model:\n\n"
  
  capture.output(cat(txt0),cat(txt1),print(round(out$modSummary,4)),file=paste0(ssimTempDir,"\\Data\\rf_output.txt"),append=TRUE) 
  
  updateRunLog("\nSummary of Model:\n")
  coeftbl <- out$modSummary
  coeftbl <- round(coeftbl, 6) 
  coeftbl <- cbind(rownames(coeftbl), coeftbl)
  rownames(coeftbl) <- NULL
  colnames(coeftbl)[1] <- "Variable"  
  updateRunLog(pander::pandoc.table.return(coeftbl, style = "simple", split.tables = 100))
  
# Test model predictions -------------------------------------------------------
  
  # Just for the training set for Random Forest we have to take out of bag predictions rather than the regular predictions
  out$data$train$predicted <- tweak.p(out$data$train$predicted) # out$trainModPredicted
  
  if(ValidationDataSheet$SplitData){
    out$data$test$predicted <- pred.fct(x=testingData, mod=finalMod, modType=modType)
  }
  
  if(ValidationDataSheet$CrossValidate){
     for(i in 1:length(out$data$cvSplits$test)){
       out$data$cvSplits$test[[i]]$predicted <- pred.fct(x=out$data$cvSplits$test[[i]], mod=finalMod, modType=modType)
     }
  }
 
# Run Cross Validation (if specified) ------------------------------------------
  
  if(ValidationDataSheet$CrossValidate){
    
    out <- cv.fct(out = out,
                  nfolds = ValidationDataSheet$NumberOfFolds)
  }
  
# Generate Model Outputs -------------------------------------------------------
 
  ## AUC/ROC - Residual Plots - Variable Importance -  Calibration - Confusion Matrix ##
  
  out <- suppressWarnings(makeModelEvalPlots(out=out))
  
  ## Response Curves ##
  
  response.curves(out)

# Save model outputs -----------------------------------------------------------

  tempFiles <- list.files(paste0(ssimTempDir,"\\Data"))
  
  # add model Outputs to datasheet
  modelOutputsSheet <- addRow(modelOutputsSheet, 
                              list(ModelType = modType, 
                                   ModelRDS = paste0(ssimTempDir,"\\Data\\rf_model.rds"),
                                   ResponseCurves = paste0(ssimTempDir,"\\Data\\rf_ResponseCurves.png"),
                                   TextOutput = paste0(ssimTempDir,"\\Data\\rf_output.txt"),
                                   ResidualSmoothPlot = paste0(ssimTempDir,"\\Data\\rf_ResidualSmoothPlot.png"),
                                   ResidualSmoothRDS = paste0(ssimTempDir,"\\Data\\rf_ResidualSmoothFunction.rds")))
  
  if(out$modelFamily != "poisson"){
    if("rf_StandardResidualPlots.png" %in% tempFiles){ modelOutputsSheet$ResidualsPlot <- paste0(ssimTempDir,"\\Data\\rf_StandardResidualPlots.png") }
    modelOutputsSheet$ConfusionMatrix <-  paste0(ssimTempDir,"\\Data\\rf_ConfusionMatrix.png")
    modelOutputsSheet$VariableImportancePlot <-  paste0(ssimTempDir,"\\Data\\rf_VariableImportance.png")
    modelOutputsSheet$VariableImportanceData <-  paste0(ssimTempDir,"\\Data\\rf_VariableImportance.csv")
    modelOutputsSheet$ROCAUCPlot <- paste0(ssimTempDir,"\\Data\\rf_ROCAUCPlot.png")
    modelOutputsSheet$CalibrationPlot <- paste0(ssimTempDir,"\\Data\\rf_CalibrationPlot.png")
  } else {
    modelOutputsSheet$ResidualsPlot <- paste0(ssimTempDir,"\\Data\\rf_PoissonResidualPlots.png")
  }
  
  if("rf_AUCPRPlot.png" %in% tempFiles){ modelOutputsSheet$AUCPRPlot <- paste0(ssimTempDir,"\\Data\\rf_AUCPRPlot.png") } 
  
  saveDatasheet(myScenario, modelOutputsSheet, "wisdm_ModelOutputs")
  
  