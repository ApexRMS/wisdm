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
  siteDataWide <- siteDataWide[complete.cases(subset(siteDataWide, select = -ModelSelectionSplit)),]
  compCases <- nrow(siteDataWide)
  if(compCases/allCases < 0.9){warning(paste(round((1-compCases/allCases)*100,digits=2),"% of cases were removed because of missing values",sep=""))}
  
  # add site weights
  siteDataWide$Weight <- 1 # TO DO: build out call to site weighting
  
  # identify training and testing sites 
  siteDataWide$UseInModelEvaluation[is.na(siteDataWide$UseInModelEvaluation)] <- FALSE
  trainTestDatasets <- split(siteDataWide, f = siteDataWide[,"UseInModelEvaluation"], drop = T)
  trainingData <- trainTestDatasets$`FALSE`
  testingData <- trainTestDatasets$`TRUE`
  
  # identify cross validation folds for model selection (if specified)
  if(ValidationDataSheet$CrossValidate){ # TO Do: build out assigns nfold number to each training row
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
  out$factorInputVars <- covariatesSheet$CovariateName[which(covariatesSheet$IsCategorical == T)]
  
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

  ## check categorical variables 
  factor.names <- out$factorInputVars # out$inputVars[1]
  factor.levels <- list()
  if(length(factor.names)>0){
    for (i in 1:length(factor.names)){
      f.col <- factor.names[i]
      x <- table(trainingData[,f.col])
      if(nrow(x)<2){
        out$bad.factor.cols <- c(out$bad.factor.cols,factor.names[i])
      }
      if(any(x<10)) {
        warning(paste("Some levels for the categorical predictor ",factor.names[i]," do not have at least 10 observations.\n",
                      "You might want to consider removing or reclassifying this predictor before continuing.\n",
                      "Factors with few observations can cause failure in model fitting when the data is split and cannot be reilably used in training a model.",sep=""))
        factor.table <- as.data.frame(x)
        colnames(factor.table) <- c("Factor Name","Factor Count")
        cat(paste("\n",factor.names[i],"\n"))
        print(factor.table)
        cat("\n\n") }
      lc.levs <-  as.numeric(row.names(x))[x>0] # make sure there is at least one "available" observation at each level
      lc.levs <- data.frame(number=lc.levs,class=lc.levs)
      factor.levels[[i]] <- lc.levs
      
      trainingData[,f.col] <- factor(trainingData[,f.col],levels=lc.levs$number,labels=lc.levs$class)
    }
  }
  if(!is.null(out$bad.factor.cols)){
    capture.output(cat("\nWarning: the following categorical response variables were removed from consideration\n",
                       "because they had only one level:",paste(out$bad.factor.cols,collapse=","),"\n"),
                   file=paste0(ssimTempDir,"\\Data\\glm_output.txt"),append=T)
  }
  
  if(nrow(trainingData)/(length(out$inputVars)-1)<10){
    capture.output(cat(paste("\n Warning: You have approximately ", round(nrow(trainingData)/(ncol(trainingData)-1),digits=1),
                             " observations for every predictor\n consider reducing the number of predictors before continuing\n",sep="")),
                   file=paste0(ssimTempDir,"\\Data\\glm_output.txt"),append=T)
  }


# Fit model --------------------------------------------------------------------

  # source(file.path(packageDir, "fitModel.R")) ## GLM - code sourced and updated from model.fit.r 

  finalMod <- fitModel(dat = trainingData, 
                      out = out,         
                      modType = modType)

  # save model to temp storage
  saveRDS(finalMod, file = paste0(ssimTempDir,"\\Data\\glm_model.rds"))
  
  # add relevant model details to out 
  out$finalMod <- finalMod
  out$nVarsFinal <- length(attr(terms(formula(finalMod)),"term.labels"))
  out$finalVars <- attr(terms(formula(finalMod)),"term.labels")
  # have to remove all the junk with powers and interactions for mess map production to work
  out$finalVars <- unique(unlist(strsplit(gsub("I\\(","",gsub("\\^2)","",out$finalVars)),":")))
  
  # add relevant model details to text output 
  txt0 <- paste("\n\n","Settings:\n","\n\t model family:  ",out$modelFamily,
              "\n\t simplification method:  ",GLMSheet$SimplificationMethod,
              "\n\n\n","Results:\n\t ","number covariates in final model:  ",length(attr(terms(formula(finalMod)),"term.labels")),"\n",sep="")
  print(finalMod$summary <- summary(finalMod))
  
  capture.output(cat(txt0),finalMod$summary,file=paste0(ssimTempDir,"\\Data\\glm_output.txt"), append=TRUE)
  cat("\n","Finished with stepwise GLM","\n")
  cat("Summary of Model:","\n")
  
  if(length(coef(finalMod))==1) stop("Null model was selected.  \nEvaluation metrics and plots will not be produced") 

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
                                   CalibrationPlot = paste0(ssimTempDir,"\\Data\\glm_CalibrationPlot.png"), 
                                   ResidualSmoothPlot = paste0(ssimTempDir,"\\Data\\glm_ResidualSmoothPlot.png"),
                                   ResidualSmoothRDS = paste0(ssimTempDir,"\\Data\\glm_ResidualSmoothFunction.rds"),
                                   ROCAUCPlot = paste0(ssimTempDir,"\\Data\\glm_ROCAUCPlot.png"),
                                   # AUCPRPlot = paste0(ssimTempDir,"\\Data\\glm_AUCPRPlot.png"),
                                   ConfusionMatrix = paste0(ssimTempDir,"\\Data\\glm_ConfusionMatrix.png"),
                                   VariableImportancePlot = paste0(ssimTempDir,"\\Data\\glm_VariableImportance.png"),
                                   VariableImportanceData = paste0(ssimTempDir,"\\Data\\glm_VariableImportance.csv")))
  
  if(ValidationDataSheet$CrossValidate){
      modelOutputsSheet$AUCPRPlot <- paste0(ssimTempDir,"\\Data\\glm_AUCPRPlot.png")
  } 
  if("glm_StandardResidualPlots.png" %in% tempFiles){
    modelOutputsSheet$ResidualsPlot <- paste0(ssimTempDir,"\\Data\\glm_StandardResidualPlots.png")
  }


  saveDatasheet(myScenario, modelOutputsSheet, "wisdm_ModelOutputs")
  