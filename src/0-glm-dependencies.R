## sdsim - dependencies
## ApexRMS, March 2022

# R version 4.1.1 
# script checks for installed packages and installs missing ones

# check rsyncrosim version and install from Git if needed
if("rsyncrosim" %in% installed.packages()[,"Package"] & packageVersion("rsyncrosim") >= "1.3.10"){ 
  library(rsyncrosim) } else {
    install.packages("https://github.com/syncrosim/rsyncrosim/releases/download/1.3.10/rsyncrosim_1.3.10.tar.gz", repos = NULL)
    library(rsyncrosim)
  }

# check remaining packages
packagesToLoad <- c("")

# install missing packages
packagesToInstall <- packagesToLoad[!(packagesToLoad %in% installed.packages()[,"Package"])]
if(length(packagesToInstall)) install.packages(packagesToInstall)

# Load packages
lapply(packagesToLoad, library, character.only = TRUE)

