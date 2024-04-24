## --------------------
## wisdm - fit maxent
## ApexRMS, April 2024
## --------------------

# built under R version 4.1.3 & SyncroSim version 3.0.0
# script pulls in pre-processed field, site and covariate data; 
# fits maxent model; builds model diagnostic and validation plots 

# source dependencies ----------------------------------------------------------

library(rsyncrosim) # install.packages("C:/GitHub/rsyncrosim", type="source", repos=NULL) 
library(tidyverse)
library(zip) # install.packages("zip")
# library(splines)

packageDir <- Sys.getenv("ssim_package_directory")
source(file.path(packageDir, "00-helper-functions.R"))
source(file.path(packageDir, "07-fit-model-functions.R"))

# disable scientific notation 
options(scipen = 999)

# Set progress bar -------------------------------------------------------------

steps <- 11
progressBar(type = "begin", totalSteps = steps)

# Connect to library -----------------------------------------------------------

# Active project and scenario
myScenario <- scenario() # datasheet(myScenario)

# Path to ssim temporary directory
ssimTempDir <- ssimEnvironment()$TransferDirectory

# Read in datasheets
covariatesSheet <- datasheet(myScenario, "wisdm_Covariates", optional = T)
modelsSheet <- datasheet(myScenario, "wisdm_Models")
fieldDataSheet <- datasheet(myScenario, "wisdm_FieldData", optional = T) %>% select(-ScenarioId)
validationDataSheet <- datasheet(myScenario, "wisdm_ValidationOptions")
retainedCovariatesSheet <- datasheet(myScenario, "wisdm_RetainedCovariates", lookupsAsFactors = F)
siteDataSheet <- datasheet(myScenario, "wisdm_SiteData", lookupsAsFactors = F) %>% select(-ScenarioId)
maxentSheet <- datasheet(myScenario, "wisdm_Maxent") %>% select(-ScenarioId)
mulitprocessingSheet <- datasheet(myScenario, "core_Multiprocessing")
modelOutputsSheet <- datasheet(myScenario, "wisdm_OutputModel", optional = T, returnInvisible = T, empty = T, lookupsAsFactors = F) %>% drop_na()

progressBar()

# Error handling ---------------------------------------------------------------

# check for for background data
if(!any(fieldDataSheet$Response == -9998)){
  if(any(siteDataSheet$Response == 0)){
    stop("Maxent is a presence-background method; please select a method suitable for presence-absense data before continuing.")
  }
  stop("Maxent is a presence-background method; please generate pseudo-absences before continuing. Note: Pseuso-absences can be generated in the '3 - Data Preperation (Non-Spatial)' Stage of the Pipeline. See 'Background Data Options' for more details.")
}

#  Set defaults ----------------------------------------------------------------  

## Validation Sheet
if(nrow(validationDataSheet)<1){
  validationDataSheet <- addRow(validationDataSheet, list(SplitData = FALSE,
                                                          CrossValidate = FALSE))
}
if(is.na(validationDataSheet$CrossValidate)){validationDataSheet$CrossValidate <- FALSE}
if(is.na(validationDataSheet$SplitData)){validationDataSheet$SplitData <- FALSE}

## Maxent Sheet 
if(nrow(maxentSheet)<1 | all(is.na(maxentSheet))){
  maxentSheet <- maxentSheet %>% drop_na()
  maxentSheet <- addRow(maxentSheet, list(MemoryLimit = 512,
                                          VisibleInterface = FALSE,
                                          MaximumBackgroundPoints = 10000,
                                          SaveMaxentFiles = FALSE))
}

if(is.na(maxentSheet$MemoryLimit)){maxentSheet$MemoryLimit <- 512}
if(is.na(maxentSheet$MaximumBackgroundPoints)){maxentSheet$MaximumBackgroundPoints <- 10000}
if(is.na(maxentSheet$VisibleInterface)){maxentSheet$VisibleInterface <- FALSE}
if(is.na(maxentSheet$SaveMaxentFiles)){maxentSheet$SaveMaxentFiles <- FALSE}

if(mulitprocessingSheet$EnableMultiprocessing){
  if(is.na(maxentSheet$MultiprocessingThreads)){
    maxentSheet$MultiprocessingThreads <- mulitprocessingSheet$MaximumJobs
    }
} else {
  maxentSheet$MultiprocessingThreads <- 1
  }
 
saveDatasheet(myScenario, maxentSheet, "wisdm_Maxent")   
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

# ignore repeat records data if present
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
if(!validationDataSheet$CrossValidate){ trainingData$ModelSelectionSplit <- FALSE }

testingData <- trainTestDatasets$`TRUE`
if(!is.null(testingData)){ testingData$ModelSelectionSplit <- FALSE }
progressBar()

# Model definitions ------------------------------------------------------------

# create object to store intermediate model selection/evaluation inputs
out <- list()

## Model type
out$modType <- modType <- "maxent"

## Model options
out$modOptions <- maxentSheet
out$modOptions$thresholdOptimization <- "Sens=Spec"	# To Do: link to defined Threshold Optimization Method in UI - currently set to default: sensitivity=specificity 

## Model family 
if(max(fieldDataSheet$Response)>1){
  stop("Maxent is a presence-background method; please select a method suitable for count data before continuing.")
  } else if(any(fieldDataSheet$Response==0)){
    stop("Maxent is a presence-background method; please select a method suitable for presence-absense data before continuing.")
    } else { out$modelFamily <- "binomial" }

## Candidate variables 
out$inputVars <- retainedCovariatesSheet$CovariatesID
out$factorInputVars <- factorInputVars

## pseudo absence
out$pseudoAbs <- pseudoAbs

## Validation options
out$validationOptions <- validationDataSheet 

## path to temp ssim storage 
out$tempDir <- ssimTempDir

# Format and save data for use in Maxent -------------------------------------------

# create temp maxent folders
dir.create(file.path(ssimTempDir, "Inputs"))
dir.create(file.path(ssimTempDir, "Outputs"))

## training data

out$swdPath <- swdPath <- file.path(ssimTempDir, "Inputs", "training-swd.csv", fsep = "\\")

trainingData %>%
  mutate(Species = case_when(Response == 1 ~ "species")) %>%
  drop_na(Species) %>%
  select(-SiteID, -Response, -UseInModelEvaluation, -ModelSelectionSplit, -Weight) %>%
  relocate(Species, .before = X) %>%
  write.csv(swdPath, row.names = F)

out$data$train <- trainingData

if(pseudoAbs){
  
  out$backgroundPath <- backgroundPath <- file.path(ssimTempDir, "Inputs", "background-swd.csv", fsep = "\\")
  
  trainingData %>%
    mutate(Species = case_when(Response != 1 ~ "background")) %>%
    drop_na(Species) %>%
    select(-SiteID, -Response, -UseInModelEvaluation, -ModelSelectionSplit, -Weight) %>%
    relocate(Species, .before = X) %>%
    write.csv(backgroundPath, row.names = F)
  }

## testing data

if(!is.null(testingData)){
  out$testDataPath <- testDataPath <- file.path(ssimTempDir, "Inputs", "testing-swd.csv", fsep = "\\")
  testingData %>%
    mutate(Species = case_when(Response == 1 ~ "species")) %>%
    drop_na(Species) %>%
    select(-SiteID, -Response, -UseInModelEvaluation, -ModelSelectionSplit, -Weight) %>%
    relocate(Species, .before = X) %>%
    write.csv(testDataPath, row.names = F)
  
  out$data$test <- testingData
  
  }

out$batchPath <- file.path(ssimTempDir, "Inputs", "runMaxent.bat", fsep = "\\")
# out$maxJobs <- mulitprocessingSheet$MaximumJobs

# Create output text file ------------------------------------------------------

capture.output(cat("Maxent Results"), file = file.path(ssimTempDir, paste0(modType, "_output.txt"))) 
on.exit(capture.output(cat("Model Failed\n\n"),file = file.path(ssimTempDir, paste0(modType, "_output.txt")),append=TRUE))  


# Review model data ------------------------------------------------------------

if(nrow(out$data$train)/(length(out$inputVars)-1)<10){
  updateRunLog(paste("\nYou have approximately ", round(nrow(out$data$train)/(ncol(out$data$train)-1),digits=1),
                     " observations for every predictor\n consider reducing the number of predictors before continuing\n",sep=""))
}
progressBar()

# Fit model --------------------------------------------------------------------

finalMod <- fitModel(dat = NULL, # maxent code pulls in data from csv files built/saved above  
                     out = out)

# finalMod$trainingData <- trainingData
# save model to temp storage
# saveRDS(finalMod, file = paste0(ssimTempDir,"\\Data\\", modType, "_model.rds"))
# finalMod$trainingData <- NULL

# add relevant model details to out 
out$finalMod <- finalMod
out$finalVars <- out$inputVars # maxent doesn't drop variables
out$nVarsFinal <- length(out$finalVars)

# load maxent run results
runReults <- read.csv(file.path(ssimTempDir, "Outputs", "maxentResults.csv"))

modSummary <- data.frame("Variabels" = gsub(".contribution", "", names(runReults)[grep("contribution",names(runReults))]))
modSummary$Contribution <- t(runReults[grep("contribution",names(runReults))])

updateRunLog("\nSummary of Model:\n")
coeftbl <- modSummary
rownames(coeftbl) <- NULL
updateRunLog(pander::pandoc.table.return(coeftbl, style = "simple", split.tables = 100))
capture.output(cat("\n\n"), modSummary, file = file.path(ssimTempDir, paste0(modType, "_output.txt")), append=TRUE)
progressBar()

# Test model predictions -------------------------------------------------------

out$data$train$predicted <- pred.fct(x=out$data$train, mod=finalMod, modType=modType) # alternatively should these be pulled from the samplePredictions.csv output by maxent??

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
finalMod$trainingData <- NULL

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
                                 ResidualSmoothRDS = file.path(ssimTempDir, paste0(modType, "_ResidualSmoothFunction.rds")),
                                 ConfusionMatrix = file.path(ssimTempDir, paste0(modType, "_ConfusionMatrix.png")),
                                 VariableImportancePlot = file.path(ssimTempDir, paste0(modType, "_VariableImportance.png")),
                                 VariableImportanceData = file.path(ssimTempDir, paste0(modType, "_VariableImportance.csv")),
                                 ROCAUCPlot = file.path(ssimTempDir, paste0(modType, "_ROCAUCPlot.png")),
                                 CalibrationPlot = file.path(ssimTempDir, paste0(modType, "_CalibrationPlot.png"))))


if("maxent_StandardResidualPlots.png" %in% tempFiles){ modelOutputsSheet$ResidualsPlot <- file.path(ssimTempDir, paste0(modType, "_StandardResidualPlots.png")) }
if("maxent_AUCPRPlot.png" %in% tempFiles){ modelOutputsSheet$AUCPRPlot <- file.path(ssimTempDir, paste0(modType, "_AUCPRPlot.png")) } 

# save maxent files to compressed zip
if(maxentSheet$SaveMaxentFiles){
  setwd(out$tempDir)
  
  # create zip file of all maxent specific inputs/outputs
  zip(zipfile = paste0("Maxent-",ssimEnvironment()$ScenarioId, ".zip"), files = c("Inputs", "Outputs"))
  
  # save zip to model outputs
  modelOutputsSheet$MaxentFiles <- filepath(ssimTempDir, paste0("Maxent-",ssimEnvironment()$ScenarioId, ".zip"))
}

# save outputs datasheet
saveDatasheet(myScenario, modelOutputsSheet, "wisdm_OutputModel", append = T)
progressBar(type = "end")
