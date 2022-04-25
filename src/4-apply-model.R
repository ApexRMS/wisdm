## --------------------
## wisdm - apply model
## ApexRMS, April 2022
## --------------------

# built under R version 4.1.1
# this transformer pulls in selected model objects and applies the models
# to specified spatial conditions 

# source dependencies ----------------------------------------------------------

packageDir <- Sys.getenv("ssim_package_directory")
# source(file.path(packageDir, "0-glm-constants.R"))
source(file.path(packageDir, "0-dependencies.R"))
source(file.path(packageDir, "0-apply-model-functions.R"))

# Connect to library -----------------------------------------------------------

# Active project and scenario
myLibrary <- ssimLibrary()
myProject <- rsyncrosim::project()
myScenario <- scenario()
# datasheet(myScenario)

# path to ssim directories
ssimTempDir <- Sys.getenv("ssim_temp_directory")
ssimOutputDir <- Sys.getenv("ssim_output_directory")
resultScenario <- Sys.getenv("ssim_scenario_id")

# Read in datasheets
# covariatesSheet <- datasheet(myProject, "wisdm_Covariates", optional = T)
# fieldDataSheet <- datasheet(myScenario, "wisdm_FieldData", optional = T)
# ValidationDataSheet <- datasheet(myScenario, "wisdm_ValidationOptions")
# reducedCovariatesSheet <- datasheet(myScenario, "wisdm_ReducedCovariates", lookupsAsFactors = F)
# siteDataSheet <- datasheet(myScenario, "wisdm_SiteData", lookupsAsFactors = F)
# GLMSheet <- datasheet(myScenario, "GLM")

covariateDataSheet <- datasheet(myScenario, "CovariateData", optional = T, lookupsAsFactors = F)
modelOutputsSheet <- datasheet(myScenario, "ModelOutputs", optional = T)
outputOptionsSheet <- datasheet(myScenario, "OutputOptions", optional = T)

# Set defaults -----------------------------------------------------------------

if(nrow(outputOptionsSheet)<1){
  outputOptionsSheet <- addRow(outputOptionsSheet, list(T,F,F))
}
if(any(is.na(outputOptionsSheet))){
  outputOptionsSheet[is.na(outputOptionsSheet)] <- F
}


# Load model object ------------------------------------------------------------

modType <- modelOutputsSheet$ModelType

mod <- readRDS(paste0(ssimOutputDir,"\\Scenario-", resultScenario,"\\wisdm_ModelOutputs\\",modelOutputsSheet$ModelRDS))

if(modType == "glm"){

  nVars <- length(attr(terms(formula(mod)),"term.labels"))
  modVars <- attr(terms(formula(mod)),"term.labels")
  # have to remove all the junk with powers and interactions for mess map production to work
  modVars <- unique(unlist(strsplit(gsub("I\\(","",gsub("\\^2)","",modVars)),":")))
  
}


# Predict model ----------------------------------------------------------------

  # get the paths off the raster files
  paths <- covariateDataSheet$RasterFilePath[which(covariateDataSheet$CovariatesID %in% modVars)]
  covData <- covariateDataSheet[which(covariateDataSheet$CovariatesID %in% modVars),]
  
  # check that all tifs are present
  if(all(file.exists(paths))){
    missingTifs <- covData$CovariatesID[!file.exists(covData$RasterFilePath)]
    stop("ERROR: the following geotiff(s) are missing from Covariate Data:  ",
         paste(missingTifs, collapse=" ,"),sep="")
  }
  
  # check that we have read access to all tiffs
  if(sum(file.access(paths),mode=0)!=0){
    stop("ERROR: the following geotiff(s) are missing : ",
         paths[(file.access(paths)!=0),][1],sep="")
  }
  
  if(modType == "glm"){
    pred.fct <- glm.predict
  }
  # if(out$input$model.source.file=="rf.r")   {
  #   pred.fct = rf.predict
  #   library(randomForest)
  # }
  # if(out$input$model.source.file=="mars.r") {
  #   pred.fct=pred.mars
  #   library(mda)
  # }
  # if(out$input$model.source.file=="brt.r")  {
  #   pred.fct=brt.predict
  #   library(gbm)
  # }
  # if(out$input$model.source.file=="brt.r") {
  #   model.covs<-levels(summary(out$mods$final.mod,plotit=FALSE)[,1])
  #   pred.fct=brt.predict
  # }
  
  maps <- proc.tiff(fit.model = mod,
                    mod.vars = modVars,
                    raster.files = paths,
                    pred.fct = pred.fct,
                    output.options = outputOptionsSheet,
                    factor.levels=out$dat$ma$factor.levels,
                    make.binary.tif=make.binary.tif,
                    thresh=out$mods$auc.output$thresh,
                    make.p.tif=make.p.tif,outfile.p=file.path(out.dir,"prob_map.tif"),
                    outfile.bin=file.path(out.dir,"bin_map.tif"),tsize=50.0,NAval=-3000,logname=logname,out=out)
}


proc.tiff(model,
          mod.vars,
                        filenames=NULL,
                        factor.levels=NA,
                        make.binary.tif=F,
                        make.p.tif=T,
                        thresh=0.5,
                        # outfile.p="brt.prob.map.tif",
                        # outfile.bin="brt.bin.map.tif",
                        tsize=2.0,
                        NAval=-3000,
                        fnames=NULL,
                        out,
                        Model){


# Load spatial layers ----------------------------------------------------------

[covariateDataSheet$CovariatesID %in% modVars]
for()
raster_i <- rast(x=covariateDataSheet$RasterFilePath[i])

