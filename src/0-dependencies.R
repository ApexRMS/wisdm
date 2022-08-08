## --------------------
## wisdm - dependencies
## ApexRMS, Aug 2022
## --------------------

# R version 4.1.3 
# script checks for installed packages and installs missing ones

options(warn=-1)

# check rsyncrosim version and install from Git if needed
if("rsyncrosim" %in% installed.packages()[,"Package"]){
  if (packageVersion("rsyncrosim") >= "1.3.13"){
    library(rsyncrosim) 
    } else {
     remove.packages("rsyncrosim")
     install.packages("https://github.com/syncrosim/rsyncrosim/releases/download/1.3.13/rsyncrosim_1.3.13.tar.gz", repos = NULL)
     library(rsyncrosim)
  }} else { 
    install.packages("rsyncrosim")
    remove.packages("rsyncrosim")
    install.packages("https://github.com/syncrosim/rsyncrosim/releases/download/1.3.13/rsyncrosim_1.3.13.tar.gz", repos = NULL)
    library(rsyncrosim)}

# # check remaining packages
# packagesToLoad_01 <- c("terra", "sf", "pander") # "raster"
# packagesToLoad_02 <- c("shiny")
# packagesToLoad_03 <- c("tidyr", "PresenceAbsence", "PRROC", "ROCR", "ggplot2", "dplyr", "splines")
# packagesToLoad_04 <- c("terra", "gbm") #, "rgdal"
# 
# packagesToLoad <- c(packagesToLoad_02, packagesToLoad_03, packagesToLoad_04)
#   
# # install missing packages
# packagesToInstall <- packagesToLoad[!(packagesToLoad %in% installed.packages()[,"Package"])]
# if(length(packagesToInstall)) install.packages(packagesToInstall)
# 
# # Load packages
# lapply(packagesToLoad, FUN = function(x){library(x, character.only = T, warn.conflicts = F)}) 
