## ---------------------------
## wisdm - variable reduction
## ApexRMS, June 2025
## ---------------------------

# built under R version 4.1.3, SyncroSim 3.1.10 & rsyncrosim 2.1.3
# Script pulls in site-specific covariate data; runs automated variable reduction
# using USDM-VIF (vif-cor) or calls shiny widget to display interactive correlation
# tool; saves reduced dataset of auto/user selected covariates

# source dependencies ----------------------------------------------------------

library(rsyncrosim) # install.packages("C:/GitHub/rsyncrosim", type="source", repos=NULL) 
library(tidyr)
library(dplyr)
library(shiny)

packageDir <- Sys.getenv("ssim_package_directory")
source(file.path(packageDir, "06-variable-reduction-functions.R"))

# Set progress bar -------------------------------------------------------------

steps <- 6
updateRunLog('6 - Variable Reduction => Begin')
progressBar(type = "begin", totalSteps = steps)

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
retainedCovariatesSheet <- datasheet(myScenario, "wisdm_RetainedCovariates")
covariateCorrelationSheet <- datasheet(myScenario, "wisdm_OutputCovariateCorrelationMatrix", optional = T) %>% drop_na()

progressBar()

# Set defaults -----------------------------------------------------------------

## Covariate selection 
if(nrow(covariateSelectionSheet)<1){
  covariateSelectionSheet <- bind_rows(covariateSelectionSheet, list(SelectionMethod = "Interactive (Correlation Viewer)"))
}

if(is.na(covariateSelectionSheet$SelectionMethod)){covariateSelectionSheet$SelectionMethod <- "Interactive (Correlation Viewer)"}

## Retained covariates

if(any(is.na(retainedCovariatesSheet$CovariatesID))){
  retainedCovariatesSheet <- na.omit(retainedCovariatesSheet)
}

# Prep inputs ------------------------------------------------------------------

# merge field and site data
siteDataWide <- spread(siteDataSheet, key = CovariatesID, value = "Value")
siteDataWide <- merge(fieldDataSheet, siteDataWide, by = "SiteID")
siteData <- select(siteDataWide, -c(SiteID, X, Y, UseInModelEvaluation, ModelSelectionSplit, Weight)) # 
rm(siteDataSheet, fieldDataSheet, siteDataWide); gc()

# identify categorical covariates and drop any with a single level
if(sum(covariatesSheet$IsCategorical, na.rm = T)>0){
  factorVars <- covariatesSheet$CovariateName[which(covariatesSheet$IsCategorical == T & covariatesSheet$CovariateName %in% names(siteData))]
  if(length(factorVars)>0){
    badFactors <- NULL
    for (i in 1:length(factorVars)){
      factor.table <- table(siteData[,factorVars[i]])
      if(length(factor.table)<2){ badFactors <- c(badFactors, factorVars[i]) }
    }
    if(length(badFactors) > 0){
      factorVars <- factorVars[-which(factorVars %in% badFactors)]
      if(length(factorVars) == 0){ factorVars <- NULL }
      updateRunLog(paste0("\nThe following categorical response variables were removed from consideration\n",
                            "because they had only one level: ",paste(badFactors, collapse=","),"\n"))
    }
  } else { 
    badFactors <- NULL
    }
} else { 
  factorVars <- NULL
  badFactors <- NULL
  }

# model family 
modelFamily <- "binomial"

# Ignore background data if present
siteData <- siteData[!siteData$Response == -9999,]

# update response for pseudo-absence sites
siteData$Response[siteData$Response == -9998] <- 0

# covariate site data
covData <- select(siteData, -Response)

progressBar()

# Automatic selection (using VIF) ----------------------------------------------

if(covariateSelectionSheet$SelectionMethod == "Automatic (Variance Inflation Factor)"){
  
  ## Covariate selection 
  if(is.na(covariateSelectionSheet$CorrelationThreshold)){covariateSelectionSheet$CorrelationThreshold <- 0.7}
  if(is.na(covariateSelectionSheet$VIFThreshold)){covariateSelectionSheet$VIFThreshold <- 10}
  
  saveDatasheet(myScenario, covariateSelectionSheet, "wisdm_CovariateSelectionOptions")
  
  source(file.path(packageDir, "usdm-vif.R"))
  
  if(length(retainedCovariatesSheet$CovariatesID)>0){
    colin <- vifstep(x = covData, th = covariateSelectionSheet$VIFThreshold, keep = as.character(retainedCovariatesSheet$CovariatesID))
  } else{
    colin <- vifstep(x = covData, th = covariateSelectionSheet$VIFThreshold)
  }
  
  selectedCovariates <-  colin@results$Variables
  
  dropCovs <- colin@excluded
  allCovs <- colin@variables
  
  updateRunLog("\n",length(dropCovs), " of the ", length(allCovs), " input variables were removed from the covartite list due to collinarity:\n\n", paste(dropCovs, collapse = "\n"))
  updateRunLog("VIFs of remaining varaibles")
  updateRunLog(pander::pandoc.table.return(colin@results, style = "simple", split.tables = 100))
  
  corMat <- colin@corMatrix
  corMat[lower.tri(corMat, diag = TRUE)] <- NA
  corVals <- unique(corMat[corMat > covariateSelectionSheet$CorrelationThreshold])
  corVals <- corVals[!is.na(corVals)]
  
  if(length(corVals)>0){
    
    updateRunLog("\nThe following correlated covariates were retained:\n")
    
    for(i in corVals){
      row <- row.names(corMat)[as.vector(which(corMat==i,arr.ind=TRUE))[1]]
      col <- colnames(corMat)[as.vector(which(corMat==i,arr.ind=TRUE))[2]]
      updateRunLog(paste0(row, " ~ ", col, ": ", round(i,4)))
    }
  }
progressBar()
} # end if "Automatic"

# Interactive Selection (using viewer) ------------------------------------------

if(covariateSelectionSheet$SelectionMethod == "Interactive (Correlation Viewer)"){
  
  ## Covariate selection 
  if(is.na(covariateSelectionSheet$DisplayHighestCorrelations)){covariateSelectionSheet$DisplayHighestCorrelations <- TRUE}
  if(is.na(covariateSelectionSheet$CorrelationThreshold)){covariateSelectionSheet$CorrelationThreshold <- 0.7}
  if(is.na(covariateSelectionSheet$NumberOfPlots)){covariateSelectionSheet$NumberOfPlots <- 5}
  
  saveDatasheet(myScenario, covariateSelectionSheet, "wisdm_CovariateSelectionOptions")
  
  # prep deviance explained data
  devExp <- vector()
  for(i in (1:ncol(covData))){
    devExp[i] <- try(my.panel.smooth(x = covData[,i], 
                                     y = siteData$Response,
                                     plot.it=FALSE,
                                     family=modelFamily),silent=TRUE)
  }
  devExp <- round(devExp,2)
  devInfo <- as.data.frame(devExp)
  devInfo$covs <- names(covData)
  devInfo$covDE <- paste0(devInfo$covs, " (", devInfo$devExp, ")")
  covsDE <- devInfo$covs
  names(covsDE) <- devInfo$covDE
  
  # run pairs explore with all variables -----------------------------------------
  
  options <- covariateSelectionSheet
  options$NumberOfPlots <- ncol(select(siteData, -Response, -all_of(badFactors)))
  
  pairsExplore(inputData = siteData,
               options = options,
               selectedCovs = names(select(siteData, -Response, -all_of(badFactors))),
               factorVars = factorVars,
               family = modelFamily,
               outputFile = file.path(ssimTempDir, "InitialCovariateCorrelationMatrix.png"))
  
  # run shiny app ----------------------------------------------------------------
  
  # inputs
  covData <- select(siteData, -Response, -all_of(badFactors))
  selectedCovariates <- names(covData)
  # covsDE
  options <- covariateSelectionSheet
  
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
  
  # portable chrome - to large to store on git 
  # browser.path = file.path(packageDir,"Apps/chrome/chrome.exe")
  
  if(is.null(browser.path)){
    runApp(appDir = file.path(packageDir, "06-covariate-correlation-app.R"),
           launch.browser = TRUE)  
  } else {
    runApp(appDir = file.path(packageDir, "06-covariate-correlation-app.R"),
         launch.browser = function(shinyurl) {
           system(paste0("\"", browser.path, "\" --app=", shinyurl, " -incognito"), wait = F)
          })
  }
  
  # save image files
  covariateCorrelationSheet <- bind_rows(
    covariateCorrelationSheet, 
    data.frame(InitialMatrix = file.path(ssimTempDir, 
                                         "InitialCovariateCorrelationMatrix.png"), 
               SelectedMatrix = file.path(ssimTempDir, 
                                          "SelectedCovariateCorrelationMatrix.png")))
  saveDatasheet(myScenario, covariateCorrelationSheet, 
                "wisdm_OutputCovariateCorrelationMatrix")
  progressBar()
} # end if "Interactive" 

# save reduced covariate list --------------------------------------------------

selectedCovariates <- unique(c(as.character(retainedCovariatesSheet$CovariatesID), selectedCovariates))

retainedCovariatesSheet <- data.frame(CovariatesID = selectedCovariates)
saveDatasheet(myScenario, retainedCovariatesSheet, "wisdm_RetainedCovariates")
progressBar(type = "end")

