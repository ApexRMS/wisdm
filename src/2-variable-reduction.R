## sdsim - variable reduction
## ApexRMS, May 2022

# R version 4.1.1
# Script pulls in site-specific covaraite data; calls shiny widget to display 
# correlation matrix and saves  reduced dataset of user selected covariates

# source dependencies ----------------------------------------------------------

packageDir <- Sys.getenv("ssim_package_directory")
source(file.path(packageDir, "0-dependencies.R"))
source(file.path(packageDir, "02-variable-reduction-functions.R"))

# Connect to library -----------------------------------------------------------

# Active project and scenario
myProject <- rsyncrosim::project()
myScenario <- scenario()

# path to ssim directories
ssimTempDir <- ssimEnvironment()$TransferDirectory 

# Read in datasheets
covariatesSheet <- datasheet(myProject, "wisdm_Covariates", optional = T, includeKey = T)
fieldDataSheet <- datasheet(myScenario, "wisdm_FieldData", optional = T)
siteDataSheet <- datasheet(myScenario, "wisdm_SiteData", lookupsAsFactors = F)
covariateSelectionSheet <- datasheet(myScenario, "wisdm_CovariateSelectionOptions", optional = T)
covariateCorrelationSheet <- datasheet(myScenario, "wisdm_CovariateCorrelationMatrix", optional = T)


# Set defaults -----------------------------------------------------------------

## Covariate selection 
if(nrow(covariateSelectionSheet)<1){
  covariateSelectionSheet <- addRow(covariateSelectionSheet, list(DisplayHighestCorrelations = TRUE,
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

# pairsExplore(inputData = siteData,
#              options = covariateSelectionSheet,
#              selectedCovs = names(select(siteData, -Response)),
#              factorVars = factorVars,
#              family = modelFamily,
#              outputFile = file.path(ssimTempDir, "CovariateCorrelationMatrix.png"))
# 

# run shiny app ----------------------------------------------------------------

covData <- select(siteData, -Response)
options <- covariateSelectionSheet

# TO DO: find better way to access default web app 
browser.path <- NULL
if(file.exists("C:/Program Files/Google/Chrome/Application/chrome.exe")){
 browser.path <- "C:/Program Files/Google/Chrome/Application/chrome.exe"
} else if(file.exists("C:/Program Files(x86)/Google/Chrome/Application/chrome.exe")){
  browser.path <- "C:/Program Files(x86)/Google/Chrome/Application/chrome.exe"
} else if(file.exists("C:/Program Files/Mozilla Firefox/firefox.exe")){
  browser.path <- "C:/Program Files/Mozilla Firefox/firefox.exe"
} else if(file.exists("C:/Program Files/Internet Explorer/iexplore.exe")){
  browser.path <- "C:/Program Files/Internet Explorer/iexplore.exe"
}

# portable chrome - to large to store on git 
# browser.path = file.path(packageDir,"Apps/chrome/chrome.exe")

if(is.null(browser.path)){
  runApp(appDir = file.path(packageDir, "02-covariate-correlation-app.R"),
         launch.browser = TRUE)  
} else {
  runApp(appDir = file.path(packageDir, "02-covariate-correlation-app.R"),
       launch.browser = function(shinyurl) {
         system(paste0("\"", browser.path, "\" --app=", shinyurl, " -incognito"), wait = F)
        })
}


# save reduced covariate list and image file -----------------------------------

reducedCovariatesSheet <- data.frame(CovariatesID = SelectedCovariates)
saveDatasheet(myScenario, reducedCovariatesSheet, "wisdm_ReducedCovariates")

covariateCorrelationSheet <- addRow(covariateCorrelationSheet, data.frame(Filename = file.path(ssimTempDir, "CovariateCorrelationMatrix.png")))
saveDatasheet(myScenario, covariateCorrelationSheet, "wisdm_CovariateCorrelationMatrix")
