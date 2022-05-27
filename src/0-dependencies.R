## --------------------
## wisdm - dependencies
## ApexRMS, March 2022
## --------------------

# R version 4.1.1 
# script checks for installed packages and installs missing ones

# check rsyncrosim version and install from Git if needed
if("rsyncrosim" %in% installed.packages()[,"Package"]){
  if (packageVersion("rsyncrosim") >= "1.3.10"){
    library(rsyncrosim)
  } else {stop("rsyncrosim version 1.3.10 or higher is required. Installation options can be found at https://syncrosim.github.io/rsyncrosim/index.html")}
  } else {stop("rsyncrosim version 1.3.10 or higher is required. Installation options can be found at https://syncrosim.github.io/rsyncrosim/index.html")}
  #   # detach("package:rsyncrosim", unload=TRUE)
  #   install.packages("https://github.com/syncrosim/rsyncrosim/releases/download/1.3.10/rsyncrosim_1.3.10.tar.gz", repos = NULL)
  #   library(rsyncrosim)
  # }} else { 
  #   install.packages("https://github.com/syncrosim/rsyncrosim/releases/download/1.3.10/rsyncrosim_1.3.10.tar.gz", repos = NULL)
  #   library(rsyncrosim)}

# check remaining packages
packagesToLoad_02 <- c("shiny")
packagesToLoad_03 <- c("tidyr", "PresenceAbsence", "PRROC", "ROCR", "ggplot2", "dplyr", "splines")
packagesToLoad_04 <- c("terra", "gbm", "rgdal")

packagesToLoad <- c(packagesToLoad_02, packagesToLoad_03, packagesToLoad_04)
  
# install missing packages
packagesToInstall <- packagesToLoad[!(packagesToLoad %in% installed.packages()[,"Package"])]
if(length(packagesToInstall)) install.packages(packagesToInstall)

# Load packages
lapply(packagesToLoad, library, character.only = TRUE)

