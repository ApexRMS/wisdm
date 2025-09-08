## ---------------------
## wisdm - apply model
## ApexRMS, August 2025
## ---------------------

# built under R version 4.1.3, SyncroSim 3.1.10 & rsyncrosim 2.1.3
# this transformer pulls in selected model objects and applies the models
# to specified spatial conditions to produce maps of occurrence probability, 
# multivariate environmental similarity surface, most dissimilar variable, 
# and residuals
# time tracking {dev code}

# Initialize first breakpoint for timing code
currentBreakPoint <- proc.time()

# Source dependencies & helpers ------------------------------------------------

library(rsyncrosim)
library(terra)
library(tidyr)
library(dplyr)
library(xml2)

packageDir <- Sys.getenv("ssim_package_directory")
source(file.path(packageDir, "00-helper-functions.R"))
source(file.path(packageDir, "09-apply-model-functions.R"))

updateRunLog('9 - Apply Model => Begin')

# Connect to library & read datasheets -----------------------------------------

# Active project and scenario
myLibrary <- ssimLibrary()
myScenario <- scenario()
# datasheet(myScenario)
updateRunLog("check 0" )

# path to ssim directories
ssimTempDir <- ssimEnvironment()$TransferDirectory 
resultScenario <- ssimEnvironment()$ScenarioId
updateRunLog("check 0.5" )

# Read in datasheets
covariatesSheet <- datasheet(myScenario, "wisdm_Covariates", optional = T)
modelsSheet <- datasheet(myScenario, "wisdm_Models") 
updateRunLog("check 1" )
spatialMultiprocessingSheet <- datasheet(myScenario, "core_SpatialMultiprocessing")
updateRunLog("check 2" )
covariateDataSheet <- datasheet(myScenario, "wisdm_CovariateData", optional = T, lookupsAsFactors = F)
templateSheet <- datasheet(myScenario, "wisdm_TemplateRaster")
restrictionSheet <- datasheet(myScenario, "wisdm_RestrictionRaster")
modelOutputsSheet <- datasheet(myScenario, "wisdm_OutputModel", optional = T, returnInvisible = T, lookupsAsFactors = F)
outputOptionsSheet <- datasheet(myScenario, "wisdm_OutputOptions", optional = T)

spatialOutputsSheet <- datasheet(myScenario, "wisdm_OutputSpatial", optional = T, lookupsAsFactors = T)
updateRunLog("check 3" )
# Set progress bar -------------------------------------------------------------

steps <- (nrow(modelOutputsSheet) * sum(outputOptionsSheet == T, na.rm = T))+4
progressBar(type = "begin", totalSteps = steps)
updateRunLog("check 4" )

# Set defaults -----------------------------------------------------------------

## Covariates sheet
if(any(is.na(covariatesSheet$ID))){
  if (all(is.na(covariatesSheet$ID))){
    covariatesSheet$ID <- 1:nrow(covariatesSheet)
  } else {
    whichNA <- which(is.na(covariatesSheet$ID))
    maxID <- max(covariatesSheet$ID, na.rm = T)
    covariatesSheet$ID[whichNA] <- (maxID+1):(maxID+length(whichNA))
  }
  saveDatasheet(myScenario, covariatesSheet,  "wisdm_Covariates")  
}
updateRunLog("check 5" )

## Output options sheet
if(nrow(outputOptionsSheet)<1){
  updateRunLog("check 5.5" )
  outputOptionsSheet[1,] <- rbind(outputOptionsSheet, list(T, T, NA, T, T, T))
}
updateRunLog("check 6" )
if(is.na(outputOptionsSheet$MakeProbabilityMap)){ outputOptionsSheet$MakeProbabilityMap <- F }
if(is.na(outputOptionsSheet$MakeBinaryMap)){ outputOptionsSheet$MakeBinaryMap <- F }
if(outputOptionsSheet$MakeBinaryMap){
  if(is.na(outputOptionsSheet$ThresholdOptimization)){
    outputOptionsSheet$ThresholdOptimization <- "Sensitivity equals specificity"
  }
}
updateRunLog("check 7" )
if(is.na(outputOptionsSheet$MakeResidualsMap)){ outputOptionsSheet$MakeResidualsMap <- F }
if(is.na(outputOptionsSheet$MakeMessMap)){ outputOptionsSheet$MakeMessMap <- F }
if(is.na(outputOptionsSheet$MakeModMap)){ outputOptionsSheet$MakeModMap <- F }
updateRunLog("check 8" )
saveDatasheet(myScenario, outputOptionsSheet,  "wisdm_OutputOptions") 

updateRunLog("Finished loading inputs in ", updateBreakpoint())

# terra output options
wopt_int <- list(datatype = "INT4S", NAflag = -9999, gdal = c("COMPRESS=LZW", "PREDICTOR=2"))

# Set up spatial multiprocessing -----------------------------------------------

maskValues <- NULL
maskExt    <- NULL
maskFile   <- NULL

if (name(myLibrary) == "Partial") {
  
  # read the job.xml to get max jobs
  jobXML <- file.path(dirname(ssimEnvironment()$LibraryFilePath), "Jobs.xml")
  jobInfo <- xml2::read_xml(jobXML)
  maxJobs <- as.numeric(xml2::xml_text(xml2::xml_find_first(jobInfo, "//LibraryMPJobsConfig")))

  # Set terra tempdir + memory  
  totalMem <- detect_system_memory()$total_gb
  sessionDetails <- setup_session(ssim_temp_dir = ssimTempDir, concurrent_sessions = maxJobs, total_ram_gb = totalMem)
  updateRunLog("check 1")

  # load multiprocessing raster
  maskFile <- rast(datasheet(myScenario, "core_SpatialMultiprocessing")$MaskFileName)
  maskExt  <- ext(maskFile)

  # identify multiprocessing tile
  tileID <- as.numeric(strsplit(strsplit(basename(ssimEnvironment()$LibraryFilePath), "-")[[1]][2], "\\.")[[1]][1])

  # define masking values
  maskValues <- unique(maskFile)[,1]
  if (length(maskValues) > 1){
    maskValues <- maskValues[maskValues != tileID]
    } else { maskValues <- NULL }

  # trim tiling raster to extent of single tile
  temp_mask <- file.path(ssimTempDir, "mask_trimmed.tif")
  maskFile <- mask(x = maskFile, mask = maskFile, maskvalues = maskValues, 
                   filename = file.path(ssimTempDir, "mask.tif"), 
                   overwrite = TRUE,
                   wopt = wopt_int)
  maskFile <- trim(x = maskFile,
                   filename = temp_mask, 
                   overwrite = TRUE,
                   wopt = wopt_int)
  unlink(file.path(ssimTempDir, "mask.tif"), force = TRUE) # remove temp mask file

  } else { 
    
    # Set terra tempdir + memory (keeps all temps local so renames are atomic)  
    maxJobs <- 1
    totalMem <- detect_system_memory()$total_gb
    sessionDetails <- setup_session(ssim_temp_dir = ssimTempDir, concurrent_sessions = maxJobs, total_ram_gb = totalMem)
  
  }

updateRunLog(sprintf("Terra memmax set to %.2f GB for %d concurrent job(s)", sessionDetails$MEMMAX_GB, maxJobs))
progressBar()

# Restriction layer (align to template later) ----------------------------------

restrictRaster <- NULL

if(nrow((restrictionSheet))>0){
  
  restrictRaster <- rast(restrictionSheet$RasterFilePath)
  
  if(!is.null(maskValues)){
    temp_restrict <- file.path(ssimTempDir, "restrict_masked.tif")
    restrictRaster <- crop(restrictRaster, maskFile,
                            filename = temp_restrict,
                            overwrite = TRUE,
                            datatype = "INT2S", NAflag = -9999, gdal = c("COMPRESS=LZW", "PREDICTOR=2"))
  } else {
    NAflag(restrictRaster) <- -9999
  }
  
}

progressBar()

# Template raster --------------------------------------------------------------

# if (nrow(templateSheet) < 1){stop("Template raster is missing. Please provide a template raster before continuing.")}

# templateRaster <- rast(templateSheet$RasterFilePath)
# if (!is.null(maskValues)){
#   temp_temp <- file.path(ssimTempDir, "template_masked.tif")
#   templateRaster <- crop(templateRaster, maskFile,
#                           filename = temp_temp,
#                           overwrite = TRUE,
#                           NAflag = -9999, gdal = c("COMPRESS=LZW", "PREDICTOR=2"))
# }

# progressBar()

# Model loop -------------------------------------------------------------------

for (i in seq_len(nrow(modelOutputsSheet))) {
  
  modType <- modelsSheet$ModelType[modelsSheet$ModelName == modelOutputsSheet$ModelsID[i]]
  mod     <- readRDS(modelOutputsSheet$ModelRDS[i])
  
  # Model-specific vars & training data
  if (modType == "glm") {
    modVars <- attr(terms(formula(mod)), "term.labels")
    modVars <- unique(unlist(strsplit(gsub("I\\(", "", gsub("\\^2\\)", "", modVars)), ":")))
    trainingData <- mod$data
    predictFct   <- glm.predict
  } else if (modType == "rf") {
    modVars      <- rownames(mod$importance)
    trainingData <- mod$trainingData
    predictFct   <- rf.predict
    library(randomForest)
  } else if (modType == "maxent") {
    modVars      <- mod$Raw.coef$V1
    trainingData <- mod$trainingData
    mod$trainingData <- NULL
    predictFct   <- maxent.predict
  } else if (modType == "brt") {
    modVars      <- mod$contributions$var
    trainingData <- mod$trainingData
    predictFct   <- brt.predict
  } else if (modType == "gam") {
    modVars      <- attr(mod$terms, "term.labels")
    trainingData <- mod$trainingData
    predictFct   <- gam.predict
  } else {
    stop("Unsupported model type: ", modType)
  }
  
  # Factor variables used by this model
  factorVars <- intersect(covariatesSheet$CovariateName[covariatesSheet$IsCategorical == TRUE], modVars)
  factor_levels <- NULL
  if (length(factorVars) == 0L) {
    factorVars <- character(0)
  } else {
    # Precompute factor levels 
    for(v in factorVars){
      lvls <- levels(trainingData[[v]])
      factor_levels[[v]] <- lvls
    }
  }
  
  # Aligned covariates stack (only variables needed by this model)
  covs <- build_cov_stack(cov_data = covariateDataSheet, 
                          var_names = modVars, 
                          factor_vars = factorVars)
  if (!is.null(maskValues)){
    covs <- crop(covs, maskFile, 
                 filename = file.path(ssimTempDir, "covariates_masked.tif"), 
                 overwrite = TRUE,
                 NAflag = -9999, gdal = c("COMPRESS=LZW", "PREDICTOR=2"))
  }
  
  progressBar()
  updateRunLog("Finished preparing mapping inputs in ", updateBreakpoint())
  
  
  # ------------------------------------------------------------------
  # Probability map (single-pass write, temp files, then rename)
  # ------------------------------------------------------------------
  
  if (outputOptionsSheet$MakeProbabilityMap) {
    
    # Add restriction band if present so we can apply it inside predict()
    # covs_for_prob <- covs
    # restrict_name <- NULL
    # if (!is.null(restrictRaster)) {
    #   covs_for_prob <- c(covs_for_prob, restrictRaster)
    #   names(covs_for_prob)[nlyr(covs_for_prob)] <- "restrict"
    #   restrict_name <- "restrict"
    # }
   
    tmp_main  <- file.path(ssimTempDir, "prob_main.tif")
    r_main <- terra::predict(
      object = covs,
      model = mod,
      fun = predict_block_int,
      factor_levels = factor_levels,
      restrict = restrictRaster,
      predictFct = predictFct,
      filename = tmp_main,
      overwrite = TRUE,
      wopt = wopt_int) # compression may not be honored in terra 1.5-21
    
    final_tmp <- tmp_main
    tmp_ext1 <- NULL
    tmp_ext2 <- NULL
    
    if (!is.null(maskValues)) {
      
      tmp_ext1 <- file.path(ssimTempDir, "prob_temp_ext.tif")
      tmp_ext2 <- file.path(ssimTempDir, "prob_final_ext.tif")
      
      r_ext <- extend(r_main, maskExt, 
                      filename = tmp_ext1, 
                      overwrite = TRUE)
      
      r_ext <- writeRaster(r_ext, filename = tmp_ext2, 
                           datatype = "INT4S", NAflag = -9999, 
                           overwrite = TRUE, gdal = c("COMPRESS=LZW", "PREDICTOR=2"))
    
      final_tmp <- tmp_ext2
    }
    
    # Release file handles before rename/unlink (Windows)
    for (nm in c("r_main","r_ext")){ if (exists(nm, inherits = FALSE)){ rm(list = nm) }}
    gc()
    
    prob_path <- file.path(ssimTempDir, paste0(modType, "_prob_map.tif"))
    if (file.exists(prob_path)) { unlink(prob_path, force = TRUE) }
    file.rename(final_tmp, prob_path)
    
    # cleanup temps that exist
    for (p in c(tmp_main, tmp_ext1, tmp_ext2)){ if (!is.null(p) && file.exists(p) && !identical(p, prob_path)) { unlink(p, force = TRUE) }}
    
    progressBar()
    updateRunLog("Finished Probability Map in ", updateBreakpoint())
  }
    
  # ------------------------------------------------------------------
  # Residuals map (block write, temp -> rename)
  # ------------------------------------------------------------------
  # if (outputOptionsSheet$MakeResidualsMap) {
  #   if (!outputOptionsSheet$MakeProbabilityMap) {
  #     updateRunLog("\nWarning: Residuals map cannot be generated without generating the Probability map.\n")
  #   } else {
  #     residSmooth <- readRDS(modelOutputsSheet$ResidualSmoothRDS[i])
  #     
  #     predrast <- rast(file.path(ssimTempDir, paste0(modType, "_prob_map.tif")))
  #     if (!is.null(maskValues)) predrast <- crop(predrast, maskFile, mask = TRUE)
  #     
  #     out_path <- file.path(ssimTempDir, paste0(modType, "_resid_map.tif"))
  #     tmp_res  <- tempfile(paste0(modType, "_resid_main_"), tmpdir = ssimTempDir, fileext = ".tif")
  #     tmp_ext  <- if (!is.null(maskValues)) tempfile(paste0(modType, "_resid_ext_"), tmpdir = ssimTempDir, fileext = ".tif") else NULL
  #     
  #     r_out <- rast(predrast); NAflag(r_out) <- -9999
  #     r_out <- writeStart(r_out, filename = tmp_res, datatype = "FLT8S", NAflag = -9999, overwrite = TRUE)
  #     
  #     nrows <- nrow(predrast); ncols <- ncol(predrast)
  #     start_row <- 1L
  #     while (start_row <= nrows) {
  #       nrows_block <- min(1024L, nrows - start_row + 1L)
  #       rows <- start_row:(start_row + nrows_block - 1L)
  #       
  #       # cell indices for these rows (no cells() in 1.5-21)
  #       cell_idx <- unlist(lapply(rows, function(rw) {
  #         start_cell <- (rw - 1L) * ncols + 1L
  #         end_cell   <- rw * ncols
  #         start_cell:end_cell
  #       }), use.names = FALSE)
  #       
  #       coords <- xyFromCell(predrast, cell_idx)
  #       prob_vals <- values(predrast, row = start_row, nrows = nrows_block, dataframe = FALSE)
  #       if (!is.null(dim(prob_vals))) prob_vals <- as.vector(prob_vals)
  #       
  #       blockvals <- data.frame(x = coords[,1], y = coords[,2])
  #       predv <- as.numeric(predict(residSmooth, blockvals))
  #       predv[is.na(prob_vals)] <- NA
  #       
  #       writeValues(r_out, predv, start = ((start_row - 1L) * ncols + 1L))
  #       start_row <- start_row + nrows_block
  #     }
  #     
  #     r_out <- writeStop(r_out)
  #     
  #     final_tmp <- tmp_res
  #     if (!is.null(maskValues)) {
  #       r_ext <- extend(rast(tmp_res), maskExt, filename = tmp_ext, overwrite = TRUE)
  #       final_tmp <- tmp_ext
  #     }
  #     
  #     # release handles then rename & cleanup
  #     for (nm in c("r_out","r_ext","predrast")) if (exists(nm, inherits = FALSE)) rm(list = nm)
  #     gc()
  #     
  #     if (file.exists(out_path)) unlink(out_path, force = TRUE)
  #     file.rename(final_tmp, out_path)
  #     for (p in c(tmp_res, tmp_ext)) if (!is.null(p) && file.exists(p) && !identical(p, out_path)) unlink(p, force = TRUE)
  #     
  #     progressBar()
  #     updateRunLog("Finished Residuals Map (streamed) in ", updateBreakpoint())
  #   }
  # }
  # 
  # # ------------------------------------------------------------------
  # # Binary map (threshold from on-disk prob; temp -> rename)
  # # ------------------------------------------------------------------
  # if (outputOptionsSheet$MakeBinaryMap) {
  #   if (!outputOptionsSheet$MakeProbabilityMap) {
  #     updateRunLog("\nWarning: Binary map cannot be generated without generating the Probability map.\n")
  #   } else {
  #     thresholds <- mod$binThresholds
  #     names(thresholds) <- c("Max kappa", "Max sensitivity and specificity", "No omission",
  #                            "Prevalence", "Sensitivity equals specificity")
  #     binThreshold <- as.numeric(thresholds[outputOptionsSheet$ThresholdOptimization])
  #     thr_int <- as.integer(round(binThreshold * 100))
  #     
  #     prob_r <- rast(file.path(ssimTempDir, paste0(modType, "_prob_map.tif")))
  #     bin_r  <- terra::ifel(prob_r >= thr_int, 1L, 0L)  # preserves NA where prob is NA
  #     
  #     bin_path <- file.path(ssimTempDir, paste0(modType, "_bin_map.tif"))
  #     tmp_bin  <- tempfile(paste0(modType, "_bin_main_"),  tmpdir = ssimTempDir, fileext = ".tif")
  #     tmp_ext  <- if (!is.null(maskValues)) tempfile(paste0(modType, "_bin_ext_"),  tmpdir = ssimTempDir, fileext = ".tif") else NULL
  #     
  #     bin_r <- writeRaster(bin_r, filename = tmp_bin, datatype = "INT2S", NAflag = -9999, overwrite = TRUE)
  #     
  #     final_tmp <- tmp_bin
  #     if (!is.null(maskValues)) {
  #       r_ext <- extend(rast(tmp_bin), maskExt, filename = tmp_ext, overwrite = TRUE)
  #       final_tmp <- tmp_ext
  #     }
  #     
  #     for (nm in c("prob_r","bin_r","r_ext")) if (exists(nm, inherits = FALSE)) rm(list = nm)
  #     gc()
  #     
  #     if (file.exists(bin_path)) unlink(bin_path, force = TRUE)
  #     file.rename(final_tmp, bin_path)
  #     for (p in c(tmp_bin, tmp_ext)) if (!is.null(p) && file.exists(p) && !identical(p, bin_path)) unlink(p, force = TRUE)
  #     
  #     progressBar()
  #     updateRunLog("Finished Binary Map (streamed) in ", updateBreakpoint())
  #   }
  # }
  # 
  # # ------------------------------------------------------------------
  # # MESS & MoD maps (temp -> rename)
  # # ------------------------------------------------------------------
  # if (outputOptionsSheet$MakeMessMap | outputOptionsSheet$MakeModMap) {
  #   if (length(factorVars)) {
  #     updateRunLog("\nWarning: MESS and MoD maps cannot be generated for models with categorical variables.\n")
  #   } else {
  #     train.dat <- select(trainingData, all_of(modVars))
  #     for (k in modVars) train.dat[, k] <- sort(train.dat[, k])
  #     prep <- prepare_mess(train.dat, var_names = modVars)
  #     
  #     if (outputOptionsSheet$MakeMessMap) {
  #       mess_path <- file.path(ssimTempDir, paste0(modType, "_mess_map.tif"))
  #       tmp_mess  <- tempfile(paste0(modType, "_mess_main_"), tmpdir = ssimTempDir, fileext = ".tif")
  #       tmp_ext   <- if (!is.null(maskValues)) tempfile(paste0(modType, "_mess_ext_"), tmpdir = ssimTempDir, fileext = ".tif") else NULL
  #       
  #       mess_fun <- function(prep, data, ...) {
  #         round(CalcMESS_fast(data, prep)$output)
  #       }
  #       
  #       r_main <- terra::predict(
  #         object    = covs[[modVars]],
  #         model     = prep,
  #         fun       = mess_fun,
  #         filename  = tmp_mess,
  #         datatype  = "INT2S",
  #         NAflag    = -9999,
  #         overwrite = TRUE
  #       )
  #       
  #       final_tmp <- tmp_mess
  #       if (!is.null(maskValues)) {
  #         r_ext <- extend(r_main, maskExt, filename = tmp_ext, overwrite = TRUE)
  #         final_tmp <- tmp_ext
  #       }
  #       
  #       for (nm in c("r_main","r_ext")) if (exists(nm, inherits = FALSE)) rm(list = nm)
  #       gc()
  #       
  #       if (file.exists(mess_path)) unlink(mess_path, force = TRUE)
  #       file.rename(final_tmp, mess_path)
  #       for (p in c(tmp_mess, tmp_ext)) if (!is.null(p) && file.exists(p) && !identical(p, mess_path)) unlink(p, force = TRUE)
  #       
  #       progressBar()
  #       updateRunLog("Finished MESS Map (streamed) in ", updateBreakpoint())
  #     }
  #     
  #     if (outputOptionsSheet$MakeModMap) {
  #       mod_path  <- file.path(ssimTempDir, paste0(modType, "_mod_map.tif"))
  #       tmp_mod   <- tempfile(paste0(modType, "_mod_main_"), tmpdir = ssimTempDir, fileext = ".tif")
  #       tmp_ext   <- if (!is.null(maskValues)) tempfile(paste0(modType, "_mod_ext_"), tmpdir = ssimTempDir, fileext = ".tif") else NULL
  #       name_to_id <- setNames(covariatesSheet$ID, covariatesSheet$CovariateName)
  #       
  #       mod_fun <- function(prep, data, name_to_id, ...) {
  #         res <- CalcMESS_fast(data, prep)
  #         ids <- unname(as.integer(name_to_id[res$Indx]))
  #         ids[is.na(ids)] <- -9999L
  #         ids
  #       }
  #       
  #       r_main <- terra::predict(
  #         object     = covs[[modVars]],
  #         model      = prep,
  #         fun        = mod_fun,
  #         name_to_id = name_to_id,
  #         filename   = tmp_mod,
  #         datatype   = "INT2S",
  #         NAflag     = -9999,
  #         overwrite  = TRUE
  #       )
  #       
  #       final_tmp <- tmp_mod
  #       if (!is.null(maskValues)) {
  #         r_ext <- extend(r_main, maskExt, filename = tmp_ext, overwrite = TRUE)
  #         final_tmp <- tmp_ext
  #       }
  #       
  #       for (nm in c("r_main","r_ext")) if (exists(nm, inherits = FALSE)) rm(list = nm)
  #       gc()
  #       
  #       if (file.exists(mod_path)) unlink(mod_path, force = TRUE)
  #       file.rename(final_tmp, mod_path)
  #       for (p in c(tmp_mod, tmp_ext)) if (!is.null(p) && file.exists(p) && !identical(p, mod_path)) unlink(p, force = TRUE)
  #       
  #       progressBar()
  #       updateRunLog("Finished MoD Map (streamed) in ", updateBreakpoint())
  #     }
  #   }
  # }

  # Release covs handles and temporary files (Windows)
  for (nm in c("covs")){ if (exists(nm, inherits = FALSE)){ rm(list = nm) }}
  gc()

  # ------------------------------------------------------------------
  # Save output paths to datasheet
  # ------------------------------------------------------------------
  spatialOutputsSheet <- addRow(
    spatialOutputsSheet,
    list(ModelsID = modelsSheet$ModelName[modelsSheet$ModelType == modType])
  )
  outputRow <- which(spatialOutputsSheet$ModelsID == modelsSheet$ModelName[modelsSheet$ModelType == modType])
  
  if (outputOptionsSheet$MakeProbabilityMap){ spatialOutputsSheet$ProbabilityRaster[outputRow] <- file.path(ssimTempDir, paste0(modType, "_prob_map.tif")) }
  if (outputOptionsSheet$MakeBinaryMap){ spatialOutputsSheet$BinaryRaster[outputRow] <- file.path(ssimTempDir, paste0(modType, "_bin_map.tif")) }
  if (outputOptionsSheet$MakeMessMap && length(factorVars) == 0L){ spatialOutputsSheet$MessRaster[outputRow] <- file.path(ssimTempDir, paste0(modType, "_mess_map.tif")) }
  if (outputOptionsSheet$MakeModMap && length(factorVars) == 0L){ spatialOutputsSheet$ModRaster[outputRow] <- file.path(ssimTempDir, paste0(modType, "_mod_map.tif")) }
  if (outputOptionsSheet$MakeResidualsMap){ spatialOutputsSheet$ResidualsRaster[outputRow] <- file.path(ssimTempDir, paste0(modType, "_resid_map.tif")) }
  
  # occasional GC after heavy writes
  if (i %% 3L == 0L) gc()
}

# Finalize
saveDatasheet(myScenario, spatialOutputsSheet, "wisdm_OutputSpatial")
progressBar(type = "end")
