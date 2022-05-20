## sdsim - variable reduction
## ApexRMS, May 2022

# R version 4.1.1
# Script pulls in site-specific covaraite data; calls shiny widget to display 
# correlation matrix and saves  reduced dataset of user selected covariates

# source dependencies ----------------------------------------------------------

packageDir <- Sys.getenv("ssim_package_directory")
source(file.path(packageDir, "0-dependencies.R"))
# source(file.path(packageDir, "0-helper-functions.R"))
# source(file.path(packageDir, "0-variable-reduction-functions.R"))

# Connect to library -----------------------------------------------------------

# Active project and scenario
myLibrary <- ssimLibrary()
myProject <- rsyncrosim::project()
myScenario <- scenario()
# datasheet(myScenario)

# path to ssim directories
ssimTempDir <- ssimEnvironment()$TransferDirectory 
ssimOutputDir <- ssimEnvironment()$OutputDirectory                                 
resultScenario <- ssimEnvironment()$ScenarioId

# Read in datasheets
covariatesSheet <- datasheet(myProject, "wisdm_Covariates", optional = T)
fieldDataSheet <- datasheet(myScenario, "wisdm_FieldData", optional = T)
siteDataSheet <- datasheet(myScenario, "wisdm_SiteData", lookupsAsFactors = F)
covariateSelectionSheet <- datasheet(myScenario, "wisdm_CovariateSelectionOptions", optional = T)
covariateCorrelationSheet <- datasheet(myScenario, "wisdm_CovariateCorrelationMatrix", optional = T)
reducedCovariatesSheet <- datasheet(myScenario, "wisdm_ReducedCovariates", lookupsAsFactors = F)


# Set defaults -----------------------------------------------------------------

## Covariate selection 
if(nrow(covariateSelectionSheet)<1){
  covariateSelectionSheet <- addRow(ValidationDataSheet, list(DisplayHighestCorrelations = TRUE,
                                                              CorrelationThreshold = 0.7, 
                                                              NumberOfPlots = 5))
}
if(is.na(covariateSelectionSheet$DisplayHighestCorrelations)){covariateSelectionSheet$DisplayHighestCorrelations <- TRUE}
if(is.na(covariateSelectionSheet$CorrelationThreshold)){covariateSelectionSheet$CorrelationThreshold <- 0.7}
if(is.na(covariateSelectionSheet$NumberOfPlots)){covariateSelectionSheet$NumberOfPlots <- 5}

# Prep inputs ------------------------------------------------------------------

# merge field and site data
siteDataWide <- spread(siteDataSheet, key = CovariatesID, value = "Value")
siteDataWide <- merge(fieldDataSheet, siteDataWide, by = "SiteID")
siteData <- select(siteDataWide, -c(SiteID, X, Y, UseInModelEvaluation, ModelSelectionSplit)) # ,Weights

# identify categorical covariates
if(sum(covariatesSheet$IsCategorical, na.rm = T)>0){
  factorVars <- covariatesSheet$CovariateName[which(covariatesSheet$IsCategorical == T)]
} else { factorVars <- NULL }

# model family 
# if response column contains only 1's and 0's response = presAbs
if(max(fieldDataSheet$Response)>1){
  modelFamily <-"poisson" 
} else { modelFamily <- "binomial" }


# run pairs explore ------------------------------------------------------------

pairsExplore(inputData = siteData,
             options = covariateSelectionSheet,
             factorVars = factorVars,
             family = modelFamily,
             outputFile = file.path(ssimTempDir, "CovariateCorrelationMatrix.png"))


