## sdsim - data preparation
## ApexRMS, March 2022

# R version 4.1.1
# Script pulls in field data and covariate data and generates dataframe of site-specific covaraite data

# source dependencies ----------------------------------------------------------

pkg_dir <- (Sys.getenv("ssim_package_directory"))
source(file.path(pkg_dir, "0-dependencies.R"))

# Connect to library -----------------------------------------------------------

# Active project and scenario
myProject <- project()
myScenario <- scenario()

# Read in datasheets
covariatesSheet <- datasheet(myProject, "sdsim_Covariates", optional = T)
fieldDataSheet <- datasheet(myScenario, "sdsim_FieldData", optional = T)
CovariateDataSheet <- datasheet(myScenario, "sdsim_CovariateData")
ValidationDataSheet <- datasheet(myScenario, "sdsim_ValidationOptions")
# siteDataSheet <- datasheet(myScenario, "sdsim_SiteData", lookupsAsFactors = F)


# Prep Field Data --------------------------------------------------------------

# update x coordinates to match template CRS

# Define Train/Test Split (if specified)
# siteSplits <- TestTrainSplit()
# fieldDataSheet$UseInModelEvaluation <- 


# convert 
siteDataWide <- spread(siteDataSheet, key = CovariatesID, value = "Value")

# merge field and site data
siteDataWide <- merge(fieldDataSheet, siteDataWide, by = "SiteID")


# Define Cross Validation folds (if specified) 
siteSplits <- CrossValidationSplit(data = siteDataWide, # site data
                                   covData = covariatesSheet,
                                   n.folds = ValidationDataSheet$NumberOfFolds,
                                   stratify = ValidationDataSheet$StratifyFolds)


# fieldDataSheet$ModelSelectionSplit <-  


