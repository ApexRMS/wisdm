## --------------------
## wisdm - parameter tuning
## ApexRMS, May 2025
## --------------------

# built under R version 4.1.3, SyncroSim 3.1.10 & rsyncrosim 2.1.2
# script pulls in pre-processed field, site and covariate data; fits multiple models 
# using defaults (or defined inputs) and provided tuning parameters; model diagnostics
# and validation plots are built for each parameter tuning combination and are displayed
# through and interactive shiny application; user selects a single model and the desired 
# tuning parameters saved to the results scenario; tuning input/output summaries are also
# saved to result scenario to preserve a record of parameter tuning results

# source dependencies ----------------------------------------------------------
  
  library(rsyncrosim)
  library(tidyr)
  library(dplyr)
  library(shiny)

  packageDir <- Sys.getenv("ssim_package_directory")
  source(file.path(packageDir, "00-helper-functions.R"))
  source(file.path(packageDir, "07-parameter-tuning-functions.R"))
  source(file.path(packageDir, "08-fit-model-functions.R"))
  

# set progress bar -------------------------------------------------------------

steps <- 30
updateRunLog('7 - Parameter Tuning => Begin')
progressBar(type = "begin", totalSteps = steps)

# Connect to library -----------------------------------------------------------

  # Active project and scenario
  myScenario <- scenario() # datasheet(myScenario)
  
  # Path to ssim temporary directory
  ssimTempDir <- ssimEnvironment()$TransferDirectory
  
  # Read in datasheets
  covariatesSheet <- datasheet(myScenario, "wisdm_Covariates", optional = T)
  modelsSheet <- datasheet(myScenario, "wisdm_Models")
  fieldDataSheet <- datasheet(myScenario, "wisdm_FieldData", optional = T)
  validationDataSheet <- datasheet(myScenario, "wisdm_ValidationOptions")
  retainedCovariatesSheet <- datasheet(myScenario, "wisdm_RetainedCovariates", lookupsAsFactors = F)
  siteDataSheet <- datasheet(myScenario, "wisdm_SiteData", lookupsAsFactors = F)
  tuningModelSheet <- datasheet(myScenario, "wisdm_TuningModel")
  
  modType <- modelsSheet$ModelType[modelsSheet$ModelName == tuningModelSheet$Model]
  if(modType == "rf"){ library(randomForest)
    tuningSheet <- datasheet(myScenario, "wisdm_rfTuning")
    modelSheet <- datasheet(myScenario, "wisdm_RF",returnScenarioInfo = F)
    modelArgs <- list(NumberOfVariablesSampled = "Number of variables sampled at split", MaximumNodes = "Maximum number of nodes", NumberOfTrees = "Number of trees", NodeSize = "Node size")
  } 
  if(modType == "brt"){ library(dismo)
    tuningSheet <- datasheet(myScenario, "wisdm_brtTuning")
    modelSheet <- datasheet(myScenario, "wisdm_BRT")
    modelArgs <- list(LearningRate = "Learning rate", NumberOfTrees = "Number of trees added per stage", BagFraction = "Bag fraction", MaximumTrees = "Maximum number of trees")
  }
  
  # output datasheets
  modelOutputsSheet <- datasheet(myScenario, "wisdm_OutputModel", optional = T, returnInvisible = T, empty = T, lookupsAsFactors = F)
  modelOutputParameterSheet <- datasheet(myScenario, "wisdm_OutputParameterTuning", optional = T, returnInvisible = T, empty = T, lookupsAsFactors = F)
  
  progressBar()

# Error handling ---------------------------------------------------------------

  # check for both presence and absence data
  if(all(fieldDataSheet$Response == 1) | all(fieldDataSheet$Response == 0)){
    stop("Please ensure that the Field Data includes both presence and absence (or pseudo-absence) data before continuing.")
  }
  
  # check that at least one tuning parameter is defined
  if(is.na(tuningSheet$Parameter1) | is.na(tuningSheet$Parameter1Values)){
    stop("No tuning parameters defined. Please select at least one parameter and provide tuning values before continuing.")
  }
 
#  Set defaults ----------------------------------------------------------------  
  
  # if response column contains only 1's and 0's response = presAbs
  if(max(fieldDataSheet$Response)>1){
    modelFamily <-"poisson" 
  } else { modelFamily <- "binomial" }
  
  ## Tuning Sheet
  tuningSheet$Parameter1 <- names(which(modelArgs == tuningSheet$Parameter1))
  v1 <- as.numeric(strsplit(tuningSheet$Parameter1Values, ",")[[1]])
  
  if(!is.na(tuningSheet$Parameter2)){
    tuningSheet$Parameter2 <- names(which(modelArgs == tuningSheet$Parameter2)) 
    v2 <- as.numeric(strsplit(tuningSheet$Parameter2Values, ",")[[1]])
    combos <- expand.grid(v1, v2)
    parameterNames <- names(combos) <- c(tuningSheet$Parameter1, tuningSheet$Parameter2)
  } else{
    combos <- expand.grid(v1)
    parameterNames <- names(combos) <- tuningSheet$Parameter1
  }
  
  # update progress bar
  steps <- 5+5*nrow(combos)
  progressBar(type = "begin", totalSteps = steps)
  
  ## Model Sheet [RF] 
  if(modType == "rf"){
    
    if(nrow(modelSheet)<1){ 
      nadf <- as.data.frame(matrix(ncol = ncol(modelSheet), nrow = 1, data = rep(NA,ncol(modelSheet))))
      names(nadf) <- names(modelSheet)
      modelSheet<- addRow(modelSheet, value = nadf) 
      }
    
    if(is.na(modelSheet$EvaluateCovariateImportance)){modelSheet$EvaluateCovariateImportance <- TRUE}
    if(is.na(modelSheet$CalculateCasewiseImportance)){modelSheet$CalculateCasewiseImportance <- FALSE}
    if(is.na(modelSheet$NodeSize)){
      if(modelFamily == "poisson"){modelSheet$NodeSize <- 5}
      if(modelFamily == "binomial"){modelSheet$NodeSize <- 1}
    }
    if(is.na(modelSheet$NumberOfTrees)){modelSheet$NumberOfTrees <- 1000}
    if(is.na(modelSheet$NormalizeVotes)){modelSheet$NormalizeVotes <- TRUE}
    if(is.na(modelSheet$CalculateProximity)){modelSheet$CalculateProximity <- FALSE}
    if(is.na(modelSheet$SampleWithReplacement)){modelSheet$SampleWithReplacement <- FALSE}
  }
  
  ## Model Sheet [BRT]
  if(modType == "brt"){
    
    modelSheet$FittingMethod <- "Use values provided below"
  
    if(is.na(modelSheet$LearningRate)){modelSheet$LearningRate <- 0.001}
    if(is.na(modelSheet$BagFraction)){modelSheet$BagFraction <- 0.75}
    if(is.na(modelSheet$MaximumTrees)){modelSheet$MaximumTrees <- 10000}
    if(is.na(modelSheet$NumberOfTrees)){modelSheet$NumberOfTrees <- 50}
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

 # Model definitions (general) -------------------------------------------------

  # create object to store intermediate model selection/evaluation inputs
  out <- list()
  
  ## Model type
  out$modType <- modType # <- modelsSheet$ModelType[modelsSheet$ModelName == tuningModelSheet$Model]
  
  ## Model family 
  out$modelFamily <- modelFamily
  
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

  capture.output(cat(as.character(tuningModelSheet$Model), "Results"), file=file.path(ssimTempDir,paste0(modType, "_output.txt"))) 
  on.exit(capture.output(cat("Model Failed\n\n"),file=file.path(ssimTempDir,paste0(modType, "_output.txt")),append=TRUE))  
  
# Review model data ------------------------------------------------------------

  if(nrow(trainingData)/(length(out$inputVars)-1)<10){
    updateRunLog(paste("\nYou have approximately ", round(nrow(trainingData)/(ncol(trainingData)-1),digits=1),
                             " observations for every predictor\n consider reducing the number of predictors before continuing\n",sep=""))
  }

  progressBar()

# Fit models -------------------------------------------------------------------
  
  outBase <- out
  outAll <- list()
  
  comboImgs <- combos
  comboImgs[,names(modelOutputsSheet)] <- NA
  
  p1Name <- paste0(modelArgs[parameterNames[1]], ": ", comboImgs[[parameterNames[1]]])
  if(length(parameterNames)>1){
    p2Name <- paste0(modelArgs[parameterNames[2]], ": ", comboImgs[[parameterNames[2]]])
    comboImgs$displayName <- displayNames <- paste0(p1Name, ", ", p2Name)
    } else { 
    comboImgs$displayName <- displayNames <- p1Name 
    }
  
  for (r in 1:nrow(combos)){ # loop fits a model for each parameter combo
    
    # update modelSheet args
    rVals <- combos[r,]
    modelSheet[names(rVals)] <- rVals
    
    # save out update model options
    out <- outBase
    out$modOptions <- modelSheet
    out$modOptions$thresholdOptimization <- "Sens=Spec"
    
    # create temp folder
    dir.create(file.path(ssimTempDir, r))
    out$tempDir <- file.path(ssimTempDir, r)
    
    if(modType == "rf"){ 
    if(!is.na(modelSheet$NumberOfVariablesSampled)){
      if(modelSheet$NumberOfVariablesSampled > nrow(retainedCovariatesSheet)){
        stop(paste0("The number of variables sampled must be between 1 and the total number of variables used in model fitting (in this case ",
                    nrow(retainedCovariatesSheet),")."))}}
    }
    # if(modType == "brt"){}
    
    # fit model 
    finalMod <- fitModel(dat = trainingData, 
                      out = out)
    
    if(is.null(finalMod)){
      
      buildFillerImage(outputPath = file.path(out$tempDir, "UnableToFitModel.png"))
      
      comboImgs[r, c("ResponseCurves", "ResidualsPlot", "ResidualSmoothPlot", "CalibrationPlot", 
                  "ROCAUCPlot", "AUCPRPlot", "ConfusionMatrix", "VariableImportancePlot")] <- file.path(out$tempDir, "UnableToFitModel.png")
      outAll[[comboImgs$displayName[r]]] <- out
      
    } else {
    
      finalMod$trainingData <- trainingData
      out$finalMod <- finalMod
      
      if(modType == "rf"){
        
        out$modOptions$NumberOfVariablesSampled <- finalMod$mtry 
        out$finalVars <- out$inputVars # random forest doesn't drop variables
        out$nVarsFinal <- length(out$finalVars)
        
        modelSummary <- finalMod$importance
        modelSummary <- modelSummary[order(modelSummary[,3],decreasing=T),]
        
        updateRunLog("\nSummary of Model:\n")
        coeftbl <- modelSummary
        coeftbl <- round(coeftbl, 6) 
        coeftbl <- cbind(rownames(coeftbl), coeftbl)
        rownames(coeftbl) <- NULL
        colnames(coeftbl)[1] <- "Variable"  
        updateRunLog(pander::pandoc.table.return(coeftbl, style = "simple", split.tables = 100))
        progressBar()
      }
      if(modType == "brt"){
        
        out$finalVars <- finalMod$contributions$var # brt doesn't drop variables
        out$nVarsFinal <- length(out$finalVars)
        
        modelSummary <- finalMod$contributions
        modelSummary <- modelSummary[order(modelSummary[,2],decreasing=T),]
        row.names(modelSummary) <- NULL
        
        updateRunLog("\nSummary of Model:\n")
        coeftbl <- modelSummary
        coeftbl[,2] <- round(coeftbl[,2], 4) 
        colnames(coeftbl) <- c("Variable", "Relative Influence")  
        updateRunLog(pander::pandoc.table.return(coeftbl, style = "simple", split.tables = 100))
        progressBar()
      }
      
      # Test model predictions -------------------------------------------------------
      
      # For the training set for Random Forest take out of bag predictions rather than the regular predictions
      if(modType == "rf"){
        if(out$modelFamily == "poisson"){ out$data$train$predicted <- finalMod$predicted
        } else { out$data$train$predicted <- tweak.p(finalMod$votes[,2]) } # tweak predictions to remove 1/0 so that calc deviance doesn't produce NA/Inf values 
      } else {
        out$data$train$predicted <- pred.fct(x=out$data$train, mod=finalMod, modType=modType)
      }
      
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
      saveRDS(finalMod, file = file.path(out$tempDir, paste0(modType, "_model.rds")))
      
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
      
      # save out tuning results
      
      tempFiles <- list.files(out$tempDir)
      
      comboImgs$ModelsID[r] <- modelsSheet$ModelName[modelsSheet$ModelType == modType]
      comboImgs$ModelRDS[r] <- file.path(out$tempDir, paste0(modType, "_model.rds"))
      comboImgs$ResponseCurves[r] <- file.path(out$tempDir, paste0(modType, "_ResponseCurves.png"))
      comboImgs$TextOutput[r] <- file.path(out$tempDir, paste0(modType, "_output.txt"))
      comboImgs$ResidualSmoothPlot[r] <- file.path(out$tempDir, paste0(modType, "_ResidualSmoothPlot.png"))
      comboImgs$ResidualSmoothRDS[r] <- file.path(out$tempDir, paste0(modType, "_ResidualSmoothFunction.rds"))
      
      if(out$modelFamily != "poisson"){
        if(paste0(modType, "_StandardResidualPlots.png") %in% tempFiles){ 
          comboImgs$ResidualsPlot[r] <- file.path(out$tempDir, paste0(modType, "_StandardResidualPlots.png")) 
          }
        comboImgs$ConfusionMatrix[r] <- file.path(out$tempDir, paste0(modType, "_ConfusionMatrix.png"))
        comboImgs$VariableImportancePlot[r] <- file.path(out$tempDir, paste0(modType, "_VariableImportance.png"))
        comboImgs$VariableImportanceData[r] <- file.path(out$tempDir, paste0(modType, "_VariableImportance.csv"))
        comboImgs$ROCAUCPlot[r] <- file.path(out$tempDir, paste0(modType, "_ROCAUCPlot.png"))
        comboImgs$CalibrationPlot[r] <- file.path(out$tempDir, paste0(modType, "_CalibrationPlot.png"))
      } else {
        comboImgs$ResidualsPlot[r] <- file.path(out$tempDir, paste0(modType, "_PoissonResidualPlots.png"))
      }
      
      if(paste0(modType, "_AUCPRPlot.png") %in% tempFiles){ 
        comboImgs$AUCPRPlot[r] <- file.path(out$tempDir, paste0(modType, "_AUCPRPlot.png")) 
      } 
      
    outAll[[comboImgs$displayName[r]]] <- out
    } # end else: model not NULL
  } # end tuning loop
  
# Generate comparison figures -------------------------------------------------- 
  
  
  buildTuningMatrices(modType = modType,
                      modelOutputsTable = comboImgs,
                      parameters = parameterNames,
                      outputPath = ssimTempDir)
  progressBar()
  
# Shiny App --------------------------------------------------------------------
  
  # TO DO: find better way to access default web app 
  browser.path <- NULL
  if(file.exists("C:/Program Files/Google/Chrome/Application/chrome.exe")){
    browser.path <- "C:/Program Files/Google/Chrome/Application/chrome.exe"
  } else if(file.exists("C:/Program Files(x86)/Google/Chrome/Application/chrome.exe")){
    browser.path <- "C:/Program Files(x86)/Google/Chrome/Application/chrome.exe"
  } else if(file.exists("C:/Program Files/Mozilla Firefox/firefox.exe")){
    browser.path <- "C:/Program Files/Mozilla Firefox/firefox.exe"
    # } else if(file.exists("C:/Program Files/Internet Explorer/iexplore.exe")){
    # browser.path <- "C:/Program Files/Internet Explorer/iexplore.exe"
  }
  
  if(is.null(browser.path)){
    runApp(appDir = file.path(packageDir, "07-parameter-tuning-app.R"),
           launch.browser = TRUE)  
  } else {
    runApp(appDir = file.path(packageDir, "07-parameter-tuning-app.R"),
           launch.browser = function(shinyurl) {
             system(paste0("\"", browser.path, "\" --app=", shinyurl, " -incognito"), wait = F)
           })
  }  
    
  selectedComboOutputs <- comboImgs[comboImgs$displayName == comboOut,]
  progressBar()

# Save model outputs -----------------------------------------------------------

  tempFiles<- list.files(ssimTempDir)
  
  # update model datasheet
  modelSheet <- outAll[[selectedComboOutputs$displayName]]$modOptions %>% select(-thresholdOptimization)
  if(modType == "rf"){ saveDatasheet(myScenario, modelSheet, "wisdm_RF")}
  if(modType == "brt"){ saveDatasheet(myScenario, modelSheet, "wisdm_BRT")}

  # save model outputs for selected model
  saveDatasheet(myScenario, selectedComboOutputs[,names(modelOutputsSheet)], "wisdm_OutputModel", append = T)
  
  # save tuning matrix outputs
  modelOutputParameterSheet <- addRow(modelOutputParameterSheet, 
                                      list(ModelsID = modelsSheet$ModelName[modelsSheet$ModelType == modType],
                                           ResponseCurves = file.path(ssimTempDir, "ResponseCurvesMatrix.png"),
                                           ResidualSmoothPlot = file.path(ssimTempDir, "ResidualSmoothPlotMatrix.png")))
  
  
  if(out$modelFamily != "poisson"){
    if("StandardResidualPlotsMatrix.png" %in% tempFiles){ modelOutputParameterSheet$ResidualsPlot <- file.path(ssimTempDir, "StandardResidualPlotsMatrix.png") }
    modelOutputParameterSheet$ConfusionMatrix <- file.path(ssimTempDir, "ConfusionMatrixMatrix.png")
    modelOutputParameterSheet$VariableImportancePlot <- file.path(ssimTempDir, "VariableImportanceMatrix.png")
    modelOutputParameterSheet$ROCAUCPlot <- file.path(ssimTempDir, "ROCAUCPlotMatrix.png")
    modelOutputParameterSheet$CalibrationPlot <- file.path(ssimTempDir, "CalibrationPlotMatrix.png")
  } else {
    modelOutputParameterSheet$ResidualsPlot <- file.path(ssimTempDir, "PoissonResidualPlotsMatrix.png")
  }
  
  if("AUCPRPlotMatrix.png" %in% tempFiles){ modelOutputParameterSheet$AUCPRPlot <- file.path(ssimTempDir, "AUCPRPlotMatrix.png") } 
  
  
  progressBar(type = "end")
  