## -------------------------
## wisdm - data preparation
## ApexRMS, March 2022
## -------------------------

# built under R version 4.1.3 & SyncroSim version 2.4.0
# Script pulls in field data and splits sites into test/train or CV groupings

# source dependencies ----------------------------------------------------------

library(rsyncrosim)
library(tidyr)
library(dplyr)
library(pander)

packageDir <- (Sys.getenv("ssim_package_directory"))
source(file.path(packageDir, "03-non-spatial-data-prep-functions.R"))

# Connect to library -----------------------------------------------------------

# Active project and scenario
myProject <- rsyncrosim::project()
myScenario <- scenario()

# temp directory
ssimTempDir <- ssimEnvironment()$TransferDirectory 

# Read in datasheets
covariatesSheet <- datasheet(myProject, "wisdm_Covariates", optional = T)
fieldDataSheet <- datasheet(myScenario, "wisdm_FieldData", optional = T)
validationDataSheet <- datasheet(myScenario, "wisdm_ValidationOptions")
siteDataSheet <- datasheet(myScenario, "wisdm_SiteData", optional = T, lookupsAsFactors = F)

#  Set defaults ----------------------------------------------------------------  

## Validation Sheet
if(nrow(validationDataSheet)<1){
  validationDataSheet <- addRow(validationDataSheet, list(SplitData = FALSE,
                                                          CrossValidate = FALSE))
}
if(is.na(validationDataSheet$SplitData)){validationDataSheet$SplitData <- FALSE}
if(validationDataSheet$SplitData){
  if(is.na(validationDataSheet$ProportionTrainingData)){
    validationDataSheet$ProportionTrainingData <- 0.5
    updateRunLog("\nTraining proportion not specified. Default value used: 0.5\n")
  }
}
if(is.na(validationDataSheet$CrossValidate)){validationDataSheet$CrossValidate <- FALSE}
if(validationDataSheet$CrossValidate){
  if(is.na(validationDataSheet$NumberOfFolds)){
    updateRunLog("\nNumber of Folds not specified. Default value used: 10\n")
    validationDataSheet$NumberOfFolds <- 10
  }
  if(is.na(validationDataSheet$StratifyFolds)){validationDataSheet$StratifyFolds <- FALSE}
}

saveDatasheet(myScenario, validationDataSheet, "wisdm_ValidationOptions")

# Prep inputs ------------------------------------------------------------------

# identify categorical covariates
if(sum(covariatesSheet$IsCategorical, na.rm = T)>0){
  factorVars <- covariatesSheet$CovariateName[which(covariatesSheet$IsCategorical == T & covariatesSheet$CovariateName %in% siteDataSheet$CovariatesID)]
  if(length(factorVars)<1){ factorVars <- NULL }
} else { factorVars <- NULL }

# drop no data (-9999) sites that resulted from spatial aggregation 
fieldDataSheet <- fieldDataSheet[!fieldDataSheet$Response == -9999,] 

# Split data for testing/training and validation -------------------------------

siteDataWide <- spread(data = siteDataSheet, key = CovariatesID, value = Value)

inputData <- left_join(fieldDataSheet, siteDataWide) # select(siteDataWide,-PixelID))

# Define Train/Test Split (if specified)
if(validationDataSheet$SplitData){
  inputData <- testTrainSplit(inputData = inputData,
                               trainProp = validationDataSheet$ProportionTrainingData,
                               # ratioPresAbs = validationDataSheet$RatioPresenceAbsence,
                               factorVars = factorVars)
} else {
  # set all data to be used in model training
  inputData$UseInModelEvaluation <- FALSE
}


# Define Cross Validation folds (if specified) 
if(validationDataSheet$CrossValidate){
  
  inputData <- crossValidationSplit(inputData = inputData,
                                     factorVars = factorVars,
                                     nFolds = validationDataSheet$NumberOfFolds,
                                     stratify = validationDataSheet$StratifyFolds)
}

# Check categorical variables and update field data sheet (if no validation specified)
if(validationDataSheet$SplitData == F & validationDataSheet$CrossValidate == F){
  if(!is.null(factorVars)){
    for (i in 1:length(factorVars)){
      factor.table <- table(inputData[,factorVars[i]])
      if(any(factor.table<10)){
        updateRunLog(paste0("\nSome levels for the categorical predictor ",factorVars[i]," do not have at least 10 observations. ", 
                                                    "Consider removing or reclassifying this predictor before continuing.\n",
                                                    "Factors with few observations can cause failure in model fitting when the data is split and cannot be reliably used in training a model.\n"))
        factor.table <- as.data.frame(factor.table)
        colnames(factor.table)<-c("Factor Name","Factor Count")
        updateRunLog(paste("\n",factorVars[i],"\n"))
        updateRunLog(pander::pandoc.table.return(factor.table, style = "rmarkdown"))
      }
    }
  }
  # set all data to be used in model training
  inputData$UseInModelEvaluation <- FALSE
}  

# save updated field data to scenario
updateFieldData <- select(inputData, names(fieldDataSheet))
saveDatasheet(myScenario, updateFieldData, "wisdm_FieldData", append = F)
