## -------------------
## wisdm - fit glm
## ApexRMS, April 2022
## -------------------

# source dependencies ----------------------------------------------------------

  packageDir <- Sys.getenv("ssim_package_directory")
  source(file.path(packageDir, "0-dependencies.R"))
  source(file.path(packageDir, "0-helper-functions.R"))
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
  GLMSheet <- datasheet(myScenario, "wisdm_GLM")
  modelOutputsSheet <- datasheet(myScenario, "wisdm_ModelOutputs", optional = T)

  
#  Set defaults ----------------------------------------------------------------  
  
  ## GLM Sheet
  if(nrow(GLMSheet)<1){
    GLMSheet <- addRow(GLMSheet, list(SelectBestPredictors = FALSE,
                                      SimplificationMethod = "AIC",
                                      ConsiderSquaredTerms = FALSE,
                                      ConsiderInteractions = FALSE))
  }
  if(is.na(GLMSheet$SelectBestPredictors)){GLMSheet$SelectBestPredictors <- FALSE}
  if(is.na(GLMSheet$SimplificationMethod)){ValidationDataSheet$SplitData <- "AIC"}
  if(is.na(GLMSheet$ConsiderSquaredTerms)){GLMSheet$ConsiderSquaredTerms <- FALSE}
  if(is.na(GLMSheet$ConsiderInteractions)){GLMSheet$ConsiderInteractions <- FALSE}
  
  saveDatasheet(myScenario, GLMSheet, "wisdm_GLM")
  
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
  siteDataWide <- siteDataWide[complete.cases(subset(siteDataWide, select = c(-ModelSelectionSplit, -Weight))),]
  compCases <- nrow(siteDataWide)
  if(compCases/allCases < 0.9){updateRunLog(paste("Warning: ", round((1-compCases/allCases)*100,digits=2),"% of cases were removed because of missing values",sep=""))}
  
  # set site weights to default of 1 if not already supplied
  if(all(is.na(siteDataWide$Weight))){siteDataWide$Weight <- 1}
  
  # ignore background data if present
  siteDataWide <- siteDataWide[!siteDataWide$Response == -9999,]
  
  # set categorical variable to factor
  factorInputVars <- covariatesSheet$CovariateName[which(covariatesSheet$IsCategorical == T & covariatesSheet$CovariateName %in% names(siteDataWide))]
  if(length(factorInputVars)>0){ siteDataWide[,factorInputVars] <- factor(siteDataWide[,factorInputVars]) }
 
  # identify training and testing sites 
  siteDataWide$UseInModelEvaluation[is.na(siteDataWide$UseInModelEvaluation)] <- FALSE
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
  }
  

# Model definitions ------------------------------------------------------------

  # create object to store intermediate model selection/evaluation inputs
  out <-list()
  
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
  out$inputVars <- reducedCovariatesSheet$CovariatesID
  out$factorInputVars <- factorInputVars
  
  ## training data
  out$data$train <- trainingData
  
  ## testing data
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

  capture.output(cat("Generalized Linear Model Results"), file=paste0(ssimTempDir,"\\Data\\glm_output.txt")) 
  on.exit(capture.output(cat("Model Failed\n\n"),file=paste0(ssimTempDir,"\\Data\\glm_output.txt"),append=TRUE))  


# Review model data ------------------------------------------------------------

  if(nrow(trainingData)/(length(out$inputVars)-1)<10){
    updateRunLog(paste("\n Warning: You have approximately ", round(nrow(trainingData)/(ncol(trainingData)-1),digits=1),
                             " observations for every predictor\n consider reducing the number of predictors before continuing\n",sep=""))
  }

  
# Fit model --------------------------------------------------------------------

  # source(file.path(packageDir, "fitModel.R")) ## GLM - code sourced and updated from model.fit.r 

  finalMod <- fitModel(dat = trainingData, 
                      out = out,         
                      modType = modType)

  # save model to temp storage
  saveRDS(finalMod, file = paste0(ssimTempDir,"\\Data\\glm_model.rds"))
  progressBar(message = "Finished with stepwise GLM")
  
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
  
  finalMod$summary <- summary(finalMod)
  
  updateRunLog("Summary of Model:\n\n")
  updateRunLog(pander::pandoc.formula.return(finalMod$summary$call))
  coeftbl <- finalMod$summary$coefficients
  coeftbl <- round(coeftbl, 6) 
  coeftbl <- cbind(rownames(coeftbl), coeftbl)
  rownames(coeftbl) <- NULL
  colnames(coeftbl) <- c("Variable", "Estimate", "Std. Error", "z value", "Pr(>|z|)")  
  updateRunLog(pander::pandoc.table.return(coeftbl, style = "simple", split.tables = 100))
  
  capture.output(cat(txt0),finalMod$summary,file=paste0(ssimTempDir,"\\Data\\glm_output.txt"), append=TRUE)
  
  if(length(coef(finalMod))==1) stop("Null model was selected. \nEvaluation metrics and plots will not be produced") 

# Test model predictions -------------------------------------------------------
  
  out$data$train$predicted <- pred.fct(x=out$data$train, mod=finalMod, modType=modType)
  
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
                                   ModelRDS = paste0(ssimTempDir,"\\Data\\glm_model.rds"),
                                   ResponseCurves = paste0(ssimTempDir,"\\Data\\glm_ResponseCurves.png"),
                                   TextOutput = paste0(ssimTempDir,"\\Data\\glm_output.txt"),
                                   ResidualSmoothPlot = paste0(ssimTempDir,"\\Data\\glm_ResidualSmoothPlot.png"),
                                   ResidualSmoothRDS = paste0(ssimTempDir,"\\Data\\glm_ResidualSmoothFunction.rds")))
  
  if(out$modelFamily != "poisson"){
    if("glm_StandardResidualPlots.png" %in% tempFiles){ modelOutputsSheet$ResidualsPlot <- paste0(ssimTempDir,"\\Data\\glm_StandardResidualPlots.png") }
    modelOutputsSheet$ConfusionMatrix <-  paste0(ssimTempDir,"\\Data\\glm_ConfusionMatrix.png")
    modelOutputsSheet$VariableImportancePlot <-  paste0(ssimTempDir,"\\Data\\glm_VariableImportance.png")
    modelOutputsSheet$VariableImportanceData <-  paste0(ssimTempDir,"\\Data\\glm_VariableImportance.csv")
    modelOutputsSheet$ROCAUCPlot <- paste0(ssimTempDir,"\\Data\\glm_ROCAUCPlot.png")
    modelOutputsSheet$CalibrationPlot <- paste0(ssimTempDir,"\\Data\\glm_CalibrationPlot.png")
  } else {
    modelOutputsSheet$ResidualsPlot <- paste0(ssimTempDir,"\\Data\\glm_PoissonResidualPlots.png")
  }
  
  if("glm_AUCPRPlot.png" %in% tempFiles){ modelOutputsSheet$AUCPRPlot <- paste0(ssimTempDir,"\\Data\\glm_AUCPRPlot.png") } 
  
  saveDatasheet(myScenario, modelOutputsSheet, "wisdm_ModelOutputs")
  