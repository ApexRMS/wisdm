## --------------------
## wisdm - fit maxent
## ApexRMS, December 2022
## --------------------

# built under R version 4.1.3 & SyncroSim version 2.4.19
# script pulls in pre-processed field, site and covariate data; 
# fits maxent model; builds model diagnostic and validation plots 

# source dependencies ----------------------------------------------------------

library(rsyncrosim)
library(tidyr)
library(dplyr)
library(zip) # install.packages("zip")
# library(splines)

packageDir <- Sys.getenv("ssim_package_directory")
source(file.path(packageDir, "00-helper-functions.R"))
source(file.path(packageDir, "04-fit-model-functions.R"))

# Connect to library -----------------------------------------------------------

# Active project and scenario
myScenario <- scenario() # datasheet(myScenario)

# Path to ssim temporary directory
ssimTempDir <- Sys.getenv("ssim_temp_directory")
ssimInputDir <- ssimEnvironment()$InputDirectory

# Read in datasheets
covariatesSheet <- datasheet(myScenario, "wisdm_Covariates", optional = T)
modelsSheet <- datasheet(myScenario, "wisdm_Models")
fieldDataSheet <- datasheet(myScenario, "wisdm_FieldData", optional = T)
validationDataSheet <- datasheet(myScenario, "wisdm_ValidationOptions")
reducedCovariatesSheet <- datasheet(myScenario, "wisdm_ReducedCovariates", lookupsAsFactors = F)
siteDataSheet <- datasheet(myScenario, "wisdm_SiteData", lookupsAsFactors = F)
maxentSheet <- datasheet(myScenario, "wisdm_Maxent")
mulitprocessingSheet <- datasheet(myScenario, "core_Multiprocessing")
modelOutputsSheet <- datasheet(myScenario, "wisdm_ModelOutputs", optional = T, empty = T, lookupsAsFactors = F)

#  Set defaults ----------------------------------------------------------------  

## Validation Sheet
if(nrow(validationDataSheet)<1){
  validationDataSheet <- addRow(validationDataSheet, list(SplitData = FALSE,
                                                          CrossValidate = FALSE))
}
if(is.na(validationDataSheet$CrossValidate)){validationDataSheet$CrossValidate <- FALSE}
if(is.na(validationDataSheet$SplitData)){validationDataSheet$SplitData <- FALSE}

## Maxent Sheet 
if(nrow(maxentSheet)<1){
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

# ignore repeat records data if present
# siteDataWide <- siteDataWide[!siteDataWide$Response == -9999,]

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
out$inputVars <- reducedCovariatesSheet$CovariatesID
out$factorInputVars <- factorInputVars

## pseudo absence  
out$pseudoAbs <- pseudoAbs <- TRUE # To Do: build out call to pseudo absence 

## Validation options
out$validationOptions <- validationDataSheet 

## path to temp ssim storage 
out$tempDir <- file.path(ssimTempDir, "Data")

# Format and save data for use in Maxent -------------------------------------------

# create temp maxent folders
dir.create(file.path(ssimTempDir, "Data", "Inputs"))
dir.create(file.path(ssimTempDir, "Data", "Outputs"))

## training data

out$swdPath <- swdPath <- file.path(ssimTempDir, "Data", "Inputs", "training-swd.csv")

trainingData %>%
  mutate(Species = case_when(Response == 1 ~ "species")) %>%
  drop_na(Species) %>%
  select(-SiteID, -Response, -UseInModelEvaluation, -ModelSelectionSplit, -Weight) %>%
  relocate(Species, .before = X) %>%
  write.csv(swdPath, row.names = F)

out$data$train <- trainingData

if(pseudoAbs){
  
  out$backgroundPath <- backgroundPath <- file.path(ssimTempDir, "Data", "Inputs", "background-swd.csv")
  
  trainingData %>%
    mutate(Species = case_when(Response != 1 ~ "background")) %>%
    drop_na(Species) %>%
    select(-SiteID, -Response, -UseInModelEvaluation, -ModelSelectionSplit, -Weight) %>%
    relocate(Species, .before = X) %>%
    write.csv(backgroundPath, row.names = F)
  }

## testing data

if(!is.null(testingData)){
  out$testDataPath <- testDataPath <- file.path(ssimTempDir, "Data", "Inputs", "testing-swd.csv")
  testingData %>%
    mutate(Species = case_when(Response == 1 ~ "species")) %>%
    drop_na(Species) %>%
    select(-SiteID, -Response, -UseInModelEvaluation, -ModelSelectionSplit, -Weight) %>%
    relocate(Species, .before = X) %>%
    write.csv(testDataPath, row.names = F)
  }

out$batchPath <- file.path(ssimTempDir,"Data", "Inputs", "runMaxent.bat")
# out$maxJobs <- mulitprocessingSheet$MaximumJobs

# Create output text file ------------------------------------------------------

capture.output(cat("Maxent Results"), file=paste0(ssimTempDir,"\\Data\\", modType, "_output.txt")) 
on.exit(capture.output(cat("Model Failed\n\n"),file=paste0(ssimTempDir,"\\Data\\", modType, "_output.txt"),append=TRUE))  


# Review model data ------------------------------------------------------------

if(nrow(out$data$train)/(length(out$inputVars)-1)<10){
  updateRunLog(paste("\nYou have approximately ", round(nrow(out$data$train)/(ncol(out$data$train)-1),digits=1),
                     " observations for every predictor\n consider reducing the number of predictors before continuing\n",sep=""))
}

# Fit model --------------------------------------------------------------------

finalMod <- fitModel(dat = NULL, # maxent code pulls in data from csv files built/saved above  
                     out = out)

finalMod$trainingData <- trainingData
  
# save model to temp storage
saveRDS(finalMod, file = paste0(ssimTempDir,"\\Data\\", modType, "_model.rds"))
finalMod$trainingData <- NULL

# add relevant model details to out 
out$finalMod <- finalMod
out$finalVars <- out$inputVars # random forest doesn't drop variables
out$nVarsFinal <- length(out$finalVars)

# load maxent run results
runReults <- read.csv(file.path(ssimTempDir, "Data", "Outputs", "maxentResults.csv"))

modSummary <- data.frame("Variabels" = gsub(".contribution", "", names(runReults)[grep("contribution",names(runReults))]))
modSummary$Contribution <- t(runReults[grep("contribution",names(runReults))])

updateRunLog("\nSummary of Model:\n")
coeftbl <- modSummary
rownames(coeftbl) <- NULL
updateRunLog(pander::pandoc.table.return(coeftbl, style = "simple", split.tables = 100))
capture.output(cat("\n\n"), modSummary, file=paste0(ssimTempDir,"\\Data\\", modType, "_output.txt"), append=TRUE)

# Test model predictions -------------------------------------------------------

out$data$train$predicted <- pred.fct(x=out$data$train, mod=finalMod, modType=modType) # alternatively should these be pulled from the samplePredictions.csv output by maxent??

if(validationDataSheet$SplitData){
  out$data$test$predicted <- pred.fct(x=out$data$test, mod=finalMod, modType=modType)
}

# Run Cross Validation (if specified) ------------------------------------------
if(validationDataSheet$CrossValidate){
  
  out <- cv.fct(out = out,
                nfolds = validationDataSheet$NumberOfFolds)
}

# Generate Model Outputs -------------------------------------------------------

## AUC/ROC - Residual Plots - Variable Importance -  Calibration - Confusion Matrix ##

out <- suppressWarnings(makeModelEvalPlots(out=out))

## Response Curves ##

response.curves(out)

# Save model outputs -----------------------------------------------------------

tempFiles <- list.files(file.path(ssimTempDir,"Data"))

# add model Outputs to datasheet
modelOutputsSheet <- addRow(modelOutputsSheet, 
                            list(ModelsID = modelsSheet$ModelName[modelsSheet$ModelType == modType],
                                 ModelRDS = paste0(ssimTempDir,"\\Data\\", modType, "_model.rds"),
                                 ResponseCurves = paste0(ssimTempDir,"\\Data\\", modType, "_ResponseCurves.png"),
                                 TextOutput = paste0(ssimTempDir,"\\Data\\", modType, "_output.txt"),
                                 ResidualSmoothPlot = paste0(ssimTempDir,"\\Data\\", modType, "_ResidualSmoothPlot.png"),
                                 ResidualSmoothRDS = paste0(ssimTempDir,"\\Data\\", modType, "_ResidualSmoothFunction.rds"),
                                 ConfusionMatrix =  paste0(ssimTempDir,"\\Data\\", modType, "_ConfusionMatrix.png"),
                                 VariableImportancePlot = paste0(ssimTempDir,"\\Data\\", modType, "_VariableImportance.png"),
                                 VariableImportanceData =  paste0(ssimTempDir,"\\Data\\", modType, "_VariableImportance.csv"),
                                 ROCAUCPlot = paste0(ssimTempDir,"\\Data\\", modType, "_ROCAUCPlot.png"),
                                 CalibrationPlot = paste0(ssimTempDir,"\\Data\\", modType, "_CalibrationPlot.png")))


if("maxent_StandardResidualPlots.png" %in% tempFiles){ modelOutputsSheet$ResidualsPlot <- paste0(ssimTempDir,"\\Data\\", modType, "_StandardResidualPlots.png") }
if("maxent_AUCPRPlot.png" %in% tempFiles){ modelOutputsSheet$AUCPRPlot <- paste0(ssimTempDir,"\\Data\\", modType, "_AUCPRPlot.png") } 

# save maxent files to compressed zip
if(maxentSheet$SaveMaxentFiles){
  setwd(out$tempDir)
  
  # create zip file of all maxent specific inputs/outputs
  zip(zipfile = paste0("Maxent-",ssimEnvironment()$ScenarioId, ".zip"), files = c("Inputs", "Outputs"))
  
  # save zip to model outputs
  modelOutputsSheet$MaxentFiles <- paste0(ssimTempDir,"\\Data\\Maxent-",ssimEnvironment()$ScenarioId, ".zip")
}

# save outputs datasheet
saveDatasheet(myScenario, modelOutputsSheet, "wisdm_ModelOutputs", append = T)
