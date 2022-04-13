## sdsim - fit glm
## ApexRMS, April 2022


# source dependencies ----------------------------------------------------------

  packageDir <- Sys.getenv("ssim_package_directory")
  source(file.path(packageDir, "0-glm-constants.R"))
  source(file.path(packageDir, "0-dependencies.R"))
  source(file.path(packageDir, "0-helper-functions.R"))

# Connect to library -----------------------------------------------------------

  # Active project and scenario
  myLibrary <- ssimLibrary()
  myProject <- rsyncrosim::project()
  myScenario <- scenario()
  # datasheet(myScenario)
  
  # path to ssim temporary directory
  ssimTempDir <- Sys.getenv("ssim_temp_directory")
  
  # Read in datasheets
  covariatesSheet <- datasheet(myProject, "sdsim_Covariates", optional = T)
  fieldDataSheet <- datasheet(myScenario, "sdsim_FieldData", optional = T)
  ValidationDataSheet <- datasheet(myScenario, "sdsim_ValidationOptions")
  reducedCovariatesSheet <- datasheet(myScenario, "sdsim_ReducedCovariates", lookupsAsFactors = F)
  siteDataSheet <- datasheet(myScenario, "sdsim_SiteData", lookupsAsFactors = F)
  GLMSheet <- datasheet(myScenario, "sdsim_GLM")
  modelOutputsSheet <- datasheet(myScenario, "sdsim_ModelOutputs", optional = T)

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
  if(ValidationDataSheet$CrossValidate == TRUE){ # TO Do: build out assigns nfold number to each training row
      trainingDataCVSplits <- split(trainingData, f = trainingData$ModelSelectionSplit, drop = T)
  }
  
  # identify cross validation folds for evaluation


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
  if(ValidationDataSheet$CrossValidate == TRUE){
    out$data$cvTrainSplits <- trainingDataCVSplits
  }
  
  ## testing data
  out$data$test <- testingData
  
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

  source(file.path(packageDir, "fitModel.R")) ## GLM - code sourced and updated from model.fit.r 

  finalMod <- fitModel(dat = trainingData, 
                      out = out,         
                      modType = modType)

  # save model to temp storage
  saveRDS(finalMod, file = paste0(ssimTempDir,"\\Data\\glm_model.rds"))
  
  # add model to Model Outputs
  modelOutputsSheet <- addRow(modelOutputsSheet, list(ModelType = modType, ModelRDS = paste0(ssimTempDir,"\\Data\\glm_model.rds")))
  # saveDatasheet(myScenario, modelOutputsSheet, "sdsim_ModelOutputs")
  
  # add relevant model details to out 
  out$finalMod <- finalMod
  out$nVarsFinal <- length(attr(terms(formula(finalMod)),"term.labels"))
  out$finalVars <- attr(terms(formula(finalMod)),"term.labels")
  # have to remove all the junk with powers and interactions for mess map production to work
  out$finalVars <- unique(unlist(strsplit(gsub("I\\(","",gsub("\\^2)","",out$finalVars)),":")))
  
  # add relevent model details to text output 
  txt0 <- paste("\n\n","Settings:\n","\n\t model family:  ",out$modelFamily,
              "\n\t simplification method:  ",GLMSheet$SimplificationMethod,
              "\n\n\n","Results:\n\t ","number covariates in final model:  ",length(attr(terms(formula(finalMod)),"term.labels")),"\n",sep="")
  # print(finalMod$summary <- summary(finalMod))
  
  capture.output(cat(txt0),finalMod$summary,file=paste0(ssimTempDir,"\\Data\\glm_output.txt"),append=TRUE)
  cat("\n","Finished with stepwise GLM","\n")
  cat("Summary of Model:","\n")
  
  if(length(coef(finalMod))==1) stop("Null model was selected.  \nEvaluation metrics and plots will not be produced") 

# Test model predictions -------------------------------------------------------
  
  out$data$train$predicted <- pred.fct(x=out$data$train, mod=finalMod, modType=modType)
  out$data$test$predicted <- pred.fct(x=testingData, mod=finalMod, modType=modType)
  
  if(ValidationDataSheet$CrossValidate == TRUE){
     for(i in 1:length(out$data$cvTrainSplits)){
       out$data$cvTrainSplits[[i]]$predicted <- pred.fct(x=out$data$cvTrainSplits[[i]], mod=finalMod, modType=modType)
     }
  }
 
  
  
# Run Cross Validation (if specified) ------------------------------------------
  
  if(ValidationDataSheet$CrossValidate == TRUE){
    
    out <- cv.fct(out = out,
                        nfolds = ValidationDataSheet$NumberOfFolds)
  }
  
# Generate Model Outputs -------------------------------------------------------
 
  # producing auc and residual plots model summary information and accross model evaluation metric
  out$mods$auc.output<-suppressWarnings(make.auc.plot.jpg(out=out))
  
  
## output plots
  

