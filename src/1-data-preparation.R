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
fieldDataSheet <- datasheet()
CovariateDataSheet <- datasheet()



