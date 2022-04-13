## sdsim - fit glm
## ApexRMS, March 2022



# source dependencies ----------------------------------------------------------

pkg_dir <- Sys.getenv("ssim_package_directory") # pkg_dir <- "C:/GitHub/sdsim/src"
source(file.path(pkg_dir, "0-dependencies.R"))

# Connect to library -----------------------------------------------------------

## for dev ##
ssimSession <- session("C:/SyncroSim Versions/2-3-11/")
myLibrary <- SsimLibrary("C:/GitHub/sdsim/demo/sdsim-demo.ssim", package = "sdsim", session = ssimSession)
myProject <- Project(myLibrary)
myScenario <- Scenario(myLibrary, 1)
# datasheet(myScenario)

covariateDataSheet <- read.csv("C:/GitHub/sdsim/demo/data/Covariate Data.csv")
siteDataSheet <- read.csv("C:/GitHub/sdsim/demo/data/Site Data.csv")
# saveDatasheet(myScenario, siteDataSheet,"sdsim_SiteData")



# Active project and scenario
myProject <- project()
myScenario <- scenario()

# Read in datasheets
# covariatesSheet <- datasheet(myProject, "sdsim_Covariates")
fieldDataSheet <- datasheet(myScenario, "sdsim_FieldData", optional = T)
# covariateDataSheet <- datasheet(myScenario, "sdsim_CovariateData", lookupsAsFactors = T)
# siteDataSheet <- datasheet(myScenario, "sdsim_SiteData")
GLMSheet <- datasheet(myScenario, "sdsim_GLM")
# outputOptionsSheet <- datasheet(myScenario, "sdsim_OutputOptions")


# Prep data for model fitting
library(tidyr)
siteDataWide <- spread(siteDataSheet, key = CovariatesID, value = "Value")

# remove variables dropped due to correlation
keepCovariates <- covariateDataSheet$CovariatesID[-which(covariateDataSheet$UseInModelSelection == "No")]
siteDataWide <- siteDataWide[,c('SiteID', keepCovariates)]

# remove testing sites 
trainingSites <- fieldDataSheet$SiteID[fieldDataSheet$UseInModelEvaluation == "No"]
trainingData <- merge(fieldDataSheet, siteDataWide, by = "SiteID")


