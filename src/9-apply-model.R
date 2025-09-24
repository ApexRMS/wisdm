
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

# path to ssim directories
ssimTempDir <- ssimEnvironment()$TransferDirectory 
resultScenario <- ssimEnvironment()$ScenarioId

# Read in datasheets
covariatesSheet <- datasheet(myScenario, "wisdm_Covariates", optional = T)
modelsSheet <- datasheet(myScenario, "wisdm_Models") 
spatialMultiprocessingSheet <- datasheet(myScenario, "core_SpatialMultiprocessing")
covariateDataSheet <- datasheet(myScenario, "wisdm_CovariateData", optional = T, lookupsAsFactors = F)
restrictionSheet <- datasheet(myScenario, "wisdm_RestrictionRaster")
modelOutputsSheet <- datasheet(myScenario, "wisdm_OutputModel", optional = T, returnInvisible = T, lookupsAsFactors = F)
outputOptionsSheet <- datasheet(myScenario, "wisdm_OutputOptions", optional = T)

spatialOutputsSheet <- datasheet(myScenario, "wisdm_OutputSpatial", optional = T, lookupsAsFactors = T)

# Set progress bar -------------------------------------------------------------

steps <- (nrow(modelOutputsSheet) * sum(outputOptionsSheet == T, na.rm = T))+4
if(name(myLibrary) == "Partial"){ steps <- steps + 2 } # extra steps for partial library setup
progressBar(type = "begin", totalSteps = steps)

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

## Output options sheet
if(nrow(outputOptionsSheet)<1){
  outputOptionsSheet <- safe_rbind(outputOptionsSheet, data.frame(
    MakeProbabilityMap = T, MakeBinaryMap = F, ThresholdOptimization = NA, 
    MakeResidualsMap = F, MakeMessMap = F, MakeModMap = F))
}
if(is.na(outputOptionsSheet$MakeProbabilityMap)){ outputOptionsSheet$MakeProbabilityMap <- F }
if(is.na(outputOptionsSheet$MakeBinaryMap)){ outputOptionsSheet$MakeBinaryMap <- F }
if(outputOptionsSheet$MakeBinaryMap){
  if(is.na(outputOptionsSheet$ThresholdOptimization)){
    outputOptionsSheet$ThresholdOptimization <- "Sensitivity equals specificity"
  }
}
if(is.na(outputOptionsSheet$MakeResidualsMap)){ outputOptionsSheet$MakeResidualsMap <- F }
if(is.na(outputOptionsSheet$MakeMessMap)){ outputOptionsSheet$MakeMessMap <- F }
if(is.na(outputOptionsSheet$MakeModMap)){ outputOptionsSheet$MakeModMap <- F }
saveDatasheet(myScenario, outputOptionsSheet,  "wisdm_OutputOptions") 

updateRunLog("Finished loading inputs in ", updateBreakpoint())

# terra output options
wopt_int <- list(datatype = "INT4S", NAflag = -9999, gdal = c("COMPRESS=LZW", "PREDICTOR=2"))
wopt_flt <- list(datatype = "FLT4S", NAflag = -9999, gdal = c("COMPRESS=LZW"))

# Set up spatial multiprocessing -----------------------------------------------

maskValues <- NULL
fullExt    <- NULL
fullMaskFile   <- NULL

if (name(myLibrary) == "Partial") {
  
  # read the job.xml to get max jobs and total jobs (i.e., tile ids)
  jobXML <- file.path(dirname(ssimEnvironment()$LibraryFilePath), "Jobs.xml")
  jobInfo <- xml2::read_xml(jobXML)
  maxJobs <- as.numeric(xml2::xml_text(xml2::xml_find_first(jobInfo, "//LibraryMPJobsConfig")))
  allJobs <- as.numeric(xml2::xml_text(xml2::xml_find_all(jobInfo, "//JobId")))

  # Set terra tempdir + memory  
  totalMem <- detect_system_memory()$total_gb
  sessionDetails <- setup_session(ssim_temp_dir = ssimTempDir, concurrent_sessions = maxJobs, total_ram_gb = totalMem)
  progressBar()

  # load multiprocessing raster
  fullMaskFile <- rast(spatialMultiprocessingSheet$MaskFileName)
  fullExt  <- ext(fullMaskFile)

  # identify multiprocessing tile
  tileID <- as.numeric(strsplit(strsplit(basename(ssimEnvironment()$LibraryFilePath), "-")[[1]][2], "\\.")[[1]][1])

  # define masking values
  maskValues <- allJobs
  if (length(maskValues) > 1){
    maskValues <- maskValues[maskValues != tileID]
    } else { maskValues <- NULL }

  # trim tiling raster to extent of single tile
  maskedGrid <- mask(x = fullMaskFile, mask = fullMaskFile, maskvalues = maskValues, 
                   filename = file.path(ssimTempDir, "full_masked.tif"), 
                   overwrite = TRUE,
                   wopt = wopt_int)
  progressBar() 
  updateRunLog("Finished tile mask in ", updateBreakpoint()) 
  tileExt <- non_NA_extent(x = maskedGrid)
  updateRunLog("Finished tile extent calculation in ", updateBreakpoint())
  # temp_mask <- file.path(ssimTempDir, "mask_trimmed.tif")
  # trimmedGrid <- trim(x = maskedGrid) # , filename = temp_mask, overwrite = TRUE, wopt = wopt_int)
  # updateRunLog("Finished trimming tile extent in ", updateBreakpoint())
  # trimmedGrid  <- writeRaster(trimmedGrid, filename = temp_mask, overwrite = TRUE)                 
  unlink(file.path(ssimTempDir, "full_masked.tif"), force = TRUE) # remove temp mask file

  } else { 
    
    # Set terra tempdir + memory (keeps all temps local so renames are atomic)  
    maxJobs <- 1
    totalMem <- detect_system_memory()$total_gb
    sessionDetails <- setup_session(ssim_temp_dir = ssimTempDir, concurrent_sessions = maxJobs, total_ram_gb = totalMem)
  
  }

updateRunLog(sprintf("Terra memmax set to %.2f GB for %d concurrent job(s)", sessionDetails$MEMMAX_GB, maxJobs))
progressBar()

# Restriction layer -------------------------------------------------------------

restrictRaster <- NULL
tmp_restrict <- NULL

if(nrow((restrictionSheet))>0){
  
  restrictRaster <- rast(restrictionSheet$RasterFilePath)
  
  if(!is.null(maskValues)){
    tmp_restrict <- file.path(ssimTempDir, "restrict_cropped.tif")
    restrictRaster <- crop(restrictRaster, tileExt,
                            filename = tmp_restrict,
                            overwrite = TRUE,
                            datatype = "INT2S", NAflag = -9999, gdal = c("COMPRESS=LZW", "PREDICTOR=2"))
  } else {
    NAflag(restrictRaster) <- -9999
  }
  
}

progressBar()

# Model loop -------------------------------------------------------------------

for (i in seq_len(nrow(modelOutputsSheet))) {
  
  modType <- modelsSheet$ModelType[modelsSheet$ModelName == modelOutputsSheet$ModelsID[i]]
  mod     <- readRDS(modelOutputsSheet$ModelRDS[i])
  
  # Model-specific vars & training data
  if (modType == "glm") {
    modVars <- attr(terms(formula(mod)), "term.labels")
    modVars <- unique(unlist(strsplit(gsub("I\\(", "", gsub("\\^2\\)", "", modVars)), ":")))
    trainingData <- mod$data
    predictFct <- glm.predict
  } else if (modType == "rf") {
    modVars <- rownames(mod$importance)
    trainingData <- mod$trainingData
    predictFct <- rf.predict
    library(randomForest)
  } else if (modType == "maxent") {
    modVars <- mod$Raw.coef$V1
    trainingData <- mod$trainingData
    mod$trainingData <- NULL
    predictFct <- maxent.predict
  } else if (modType == "brt") {
    modVars <- mod$contributions$var
    trainingData <- mod$trainingData
    predictFct <- brt.predict
  } else if (modType == "gam") {
    modVars <- attr(mod$terms, "term.labels")
    trainingData <- mod$trainingData
    predictFct <- gam.predict
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
    tmp_covs <- file.path(ssimTempDir, "covariates_temp.tif")
    covs <- crop(covs, tileExt, 
                 filename = tmp_covs, 
                 overwrite = TRUE,
                 NAflag = -9999, gdal = c("COMPRESS=LZW", "PREDICTOR=2"))
  }
  
  progressBar()
  updateRunLog("Finished preparing mapping inputs for ", modType, " in ", updateBreakpoint())

  # ------------------------------------------------------------------
  # Probability map 
  # ------------------------------------------------------------------
  
  if (outputOptionsSheet$MakeProbabilityMap) {
    
    # Add restriction band if present so we can apply it inside predict()
    covs_for_prob <- covs
    restrict_name <- NULL
    if (!is.null(restrictRaster)) {
      covs_for_prob <- c(covs_for_prob, restrictRaster)
      names(covs_for_prob)[nlyr(covs_for_prob)] <- "restrict"
      restrict_name <- "restrict"
    }
   
    tmp_main  <- file.path(ssimTempDir, "prob_main.tif")
    r_main <- terra::predict(
      object = covs_for_prob,
      model = mod,
      fun = predict_block_int,
      factor_levels = factor_levels,
      restrict_col = restrict_name,
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
      
      r_ext <- extend(r_main, fullExt, 
                      filename = tmp_ext1, 
                      overwrite = TRUE)
      
      r_ext <- writeRaster(r_ext, filename = tmp_ext2, 
                           datatype = "INT4S", NAflag = -9999, 
                           overwrite = TRUE, gdal = c("COMPRESS=LZW", "PREDICTOR=2"))
    
      final_tmp <- tmp_ext2
    }
    
    # Release file handles before rename/unlink (Windows)
    for (nm in c("r_main","r_ext", "covs_for_prob")){ if (exists(nm, inherits = FALSE)){ rm(list = nm) }}
    gc()
    
    prob_path <- file.path(ssimTempDir, paste0(modType, "_prob_map.tif"))
    if (file.exists(prob_path)) { unlink(prob_path, force = TRUE) }
    file.rename(final_tmp, prob_path)
    
    # cleanup temps that exist
    for (p in c(tmp_main, tmp_ext1, tmp_ext2)){ if (!is.null(p) && file.exists(p) && !identical(p, prob_path)) { unlink(p, force = TRUE) }}
    
    progressBar()
    updateRunLog("Finished Probability Map for ", modType, " in ", updateBreakpoint())
  }
    
  # ------------------------------------------------------------------
  # Residuals map 
  # ------------------------------------------------------------------

  if (outputOptionsSheet$MakeResidualsMap) {

    if (!outputOptionsSheet$MakeProbabilityMap) {
      updateRunLog("\nWarning: Residuals map cannot be generated without generating the Probability map.\n")
    
    } else {
      residSmooth <- readRDS(modelOutputsSheet$ResidualSmoothRDS[i])

      prob_path <- file.path(ssimTempDir, paste0(modType, "_prob_map.tif"))
      prob_r <- terra::rast(prob_path)
      if (!is.null(maskValues)){ prob_r <- terra::crop(prob_r, tileExt) }

      out_path <- file.path(ssimTempDir, paste0(modType, "_resid_map.tif"))
      tmp_res <- file.path(ssimTempDir, "resid_stage1.tif")
      tmp_ext <- if (!is.null(maskValues)){ file.path(ssimTempDir, "resid_stage2ext.tif") } else { NULL }
      tmp_fin <- file.path(ssimTempDir, "resid_finalized.tif")

      nrows <- terra::nrow(prob_r)
      ncols <- terra::ncol(prob_r)

      # Precompute column/row coordinates
      x_cols <- terra::xFromCol(prob_r, 1:ncols)     # length = ncols
      y_rows <- terra::yFromRow(prob_r, 1:nrows)     # length = nrows

      # Stage-1: open a streaming writer (no datatype/compression here)
      wr <- terra::rast(prob_r)
      terra::writeStart(wr, filename = tmp_res, overwrite = TRUE)

      block_rows <- if (!is.null(sessionDetails$BLOCK_ROWS)){ sessionDetails$BLOCK_ROWS } else { 384L }
      start_row  <- 1L
      while (start_row <= nrows) {
        nrb <- min(block_rows, nrows - start_row + 1L)
        rows <- start_row:(start_row + nrb - 1L)

        # Arithmetic coords for this block
        xvec <- rep.int(x_cols, times = nrb)         # length = ncols * nrb
        yvec <- rep(y_rows[rows], each = ncols)

        # Predict residual surface for this block
        predv <- as.numeric(stats::predict(residSmooth, data.frame(x = xvec, y = yvec)))

        # Apply NA mask from probability map (preserve NA semantics)
        pv <- terra::readValues(prob_r, row = start_row, nrows = nrb, dataframe = FALSE)
        if (!is.null(dim(pv))){ pv <- as.vector(pv) }
        predv[is.na(pv)] <- NA_real_

        # Write this block
        terra::writeValues(wr, predv, start = start_row, nrows = nrb)

        start_row <- start_row + nrb
      }

      terra::writeStop(wr)  # flush tmp_res

      # Extend to full extent (if needed)
      src_for_final <- tmp_res
      if (!is.null(maskValues)) {
        r_ext <- terra::extend(terra::rast(tmp_res), fullExt, filename = tmp_ext, overwrite = TRUE)
        rm(r_ext)
        src_for_final <- tmp_ext
      }

      # Enforce float type + NAflag + compression
      src_r  <- terra::rast(src_for_final)
      r_fin <- terra::writeRaster(src_r, filename = tmp_fin,
                                  overwrite = TRUE, wopt = wopt_flt)
      rm(r_fin, src_r); gc()

      # Clean temps
      if (file.exists(out_path)) { unlink(out_path, force = TRUE) }
      file.rename(tmp_fin, out_path)
      
      for (p in c(tmp_res, if (!is.null(maskValues)){ src_for_final })) {
        if (!is.null(p) && file.exists(p) && !identical(p, out_path)) { unlink(p, force = TRUE) }
      }

      for (nm in c("prob_r","wr","x_cols","y_rows")) if (exists(nm, inherits = FALSE)){ rm(list = nm) }
      gc()

      progressBar()
      updateRunLog("Finished Residuals Map for ", modType, " in ", updateBreakpoint())
    }
  }
  # ------------------------------------------------------------------
  # Binary map 
  # ------------------------------------------------------------------

  if (outputOptionsSheet$MakeBinaryMap) {

    if (!outputOptionsSheet$MakeProbabilityMap) {
      updateRunLog("\nWarning: Binary map cannot be generated without generating the Probability map.\n")
    
    } else {
      thresholds <- mod$binThresholds
      names(thresholds) <- c("Max kappa", "Max sensitivity and specificity", "No omission",
                            "Prevalence", "Sensitivity equals specificity")
      binThreshold <- as.numeric(thresholds[outputOptionsSheet$ThresholdOptimization])
      thr_int <- as.integer(round(binThreshold * 100))

      prob_r  <- rast(file.path(ssimTempDir, paste0(modType, "_prob_map.tif")))
      if (!is.null(maskValues)) { prob_r <- crop(prob_r, tileExt) }

      bin_r   <- terra::ifel(prob_r >= thr_int, 1L, 0L)  # preserves NA from prob
      bin_path <- file.path(ssimTempDir, paste0(modType, "_bin_map.tif"))

      tmp_bin <- file.path(ssimTempDir, "bin_stage1.tif")
      tmp_ext <- if (!is.null(maskValues)){ file.path(ssimTempDir, "bin_stage2ext.tif") } else { NULL }
      tmp_fin <- file.path(ssimTempDir, "bin_finalized.tif")

      # Stage-1 write (datatype later)
      bin_r <- writeRaster(bin_r, filename = tmp_bin, overwrite = TRUE)

      src_for_final <- tmp_bin
      if (!is.null(maskValues)) {
        r_ext <- extend(rast(tmp_bin), fullExt, filename = tmp_ext, overwrite = TRUE)
        src_for_final <- tmp_ext
        if (exists("r_ext", inherits = FALSE)) rm(r_ext)
      }

      # Stage-2 enforce INT4S (+ NAflag/compression)
      r_fin <- writeRaster(rast(src_for_final), filename = tmp_fin, overwrite = TRUE, wopt = wopt_int)
      rm(r_fin)

      if (file.exists(bin_path)) { unlink(bin_path, force = TRUE) }
      file.rename(tmp_fin, bin_path)
      for (p in c(tmp_bin, tmp_ext)) if (!is.null(p) && file.exists(p) && !identical(p, bin_path)) { unlink(p, force = TRUE) }

      for (nm in c("prob_r","bin_r")) if (exists(nm, inherits = FALSE)) { rm(list = nm) }
      gc()

      progressBar()
      updateRunLog("Finished Binary Map for ", modType, " in ", updateBreakpoint())
    }
  }

  # ------------------------------------------------------------------
  # MESS & MoD maps 
  # ------------------------------------------------------------------

  if (outputOptionsSheet$MakeMessMap | outputOptionsSheet$MakeModMap) {
    
    if (length(factorVars)) {
      updateRunLog("\nWarning: MESS and MoD maps cannot be generated for models with categorical variables.\n")
    
    } else {
      # Prep training data order & precompute structures
      train.dat <- dplyr::select(trainingData, dplyr::all_of(modVars))
      for (k in modVars){ train.dat[, k] <- sort(train.dat[, k]) }
      prep <- prepare_mess(train.dat, var_names = modVars)

      # ----- MESS -----
      if (outputOptionsSheet$MakeMessMap) {
        mess_path <- file.path(ssimTempDir, paste0(modType, "_mess_map.tif"))
        tmp_mess <- file.path(ssimTempDir, "mess_stage1.tif")
        tmp_ext <- if (!is.null(maskValues)){ file.path(ssimTempDir, "mess_stage2ext.tif") } else { NULL }
        tmp_fin <- file.path(ssimTempDir, "mess_finalized.tif")

        r_main <- terra::predict(
          object    = covs,
          model     = prep,
          fun       = mess_fun,
          filename  = tmp_mess,
          overwrite = TRUE
        )

        src_for_final <- tmp_mess
        if (!is.null(maskValues)) {
          r_ext <- extend(r_main, fullExt, filename = tmp_ext, overwrite = TRUE)
          src_for_final <- tmp_ext
          if (exists("r_ext", inherits = FALSE)) { rm(r_ext) }
        }

        r_fin <- writeRaster(rast(src_for_final), filename = tmp_fin, overwrite = TRUE, wopt = wopt_int)
        rm(r_fin)

        if (file.exists(mess_path)) { unlink(mess_path, force = TRUE) }
        file.rename(tmp_fin, mess_path)
        for (p in c(tmp_mess, tmp_ext)) if (!is.null(p) && file.exists(p) && !identical(p, mess_path)) { unlink(p, force = TRUE) }

        for (nm in c("r_main")) if (exists(nm, inherits = FALSE)) { rm(list = nm) }
        gc()

        progressBar()
        updateRunLog("Finished MESS Map for ", modType, " in ", updateBreakpoint())
      }

      # ----- MoD -----
      if (outputOptionsSheet$MakeModMap) {
        mod_path <- file.path(ssimTempDir, paste0(modType, "_mod_map.tif"))
        tmp_mod <- file.path(ssimTempDir, "mod_stage1.tif")
        tmp_ext <- if (!is.null(maskValues)){ file.path(ssimTempDir, "mod_stage2ext.tif") } else { NULL }
        tmp_fin <- file.path(ssimTempDir, "mod_finalized.tif")
        name_to_id <- setNames(covariatesSheet$ID, covariatesSheet$CovariateName)

        r_main <- terra::predict(
          object     = covs,
          model      = prep,
          fun        = mod_fun,
          name_to_id = name_to_id,
          filename   = tmp_mod,
          overwrite  = TRUE
        )

        src_for_final <- tmp_mod
        if (!is.null(maskValues)) {
          r_ext <- extend(r_main, fullExt, filename = tmp_ext, overwrite = TRUE)
          src_for_final <- tmp_ext
          if (exists("r_ext", inherits = FALSE)) { rm(r_ext) }
        }

        r_fin <- writeRaster(rast(src_for_final), filename = tmp_fin, overwrite = TRUE, wopt = wopt_int)
        rm(r_fin)

        if (file.exists(mod_path)) { unlink(mod_path, force = TRUE) }
        file.rename(tmp_fin, mod_path)
        for (p in c(tmp_mod, tmp_ext)) if (!is.null(p) && file.exists(p) && !identical(p, mod_path)) { unlink(p, force = TRUE) }

        for (nm in c("r_main")) if (exists(nm, inherits = FALSE)) { rm(list = nm) }
        gc()

        progressBar()
        updateRunLog("Finished MoD Map for ", modType, " in ", updateBreakpoint())
      }
    }
  }

  # Release covs handles and temporary files (Windows)
  for (nm in c("covs")){ if (exists(nm, inherits = FALSE)){ rm(list = nm) }}
  for (p in c(tmp_covs)){ unlink(p, force = TRUE) }
  gc()
  
  # ------------------------------------------------------------------
  # Save output paths to datasheet
  # ------------------------------------------------------------------
  spatialOutputsSheet <- safe_rbind(spatialOutputsSheet, 
                                  data.frame(ModelsID = modelsSheet$ModelName[modelsSheet$ModelType == modType]))
  outputRow <- nrow(spatialOutputsSheet)
  
  if (outputOptionsSheet$MakeProbabilityMap){ spatialOutputsSheet$ProbabilityRaster[outputRow] <- file.path(ssimTempDir, paste0(modType, "_prob_map.tif")) }
  if (outputOptionsSheet$MakeBinaryMap){ spatialOutputsSheet$BinaryRaster[outputRow] <- file.path(ssimTempDir, paste0(modType, "_bin_map.tif")) }
  if (outputOptionsSheet$MakeMessMap && length(factorVars) == 0L){ spatialOutputsSheet$MessRaster[outputRow] <- file.path(ssimTempDir, paste0(modType, "_mess_map.tif")) }
  if (outputOptionsSheet$MakeModMap && length(factorVars) == 0L){ spatialOutputsSheet$ModRaster[outputRow] <- file.path(ssimTempDir, paste0(modType, "_mod_map.tif")) }
  if (outputOptionsSheet$MakeResidualsMap){ spatialOutputsSheet$ResidualsRaster[outputRow] <- file.path(ssimTempDir, paste0(modType, "_resid_map.tif")) }
  
  # occasional GC after heavy writes
  if (i %% 3L == 0L) gc()
}

# release/remove any remaining temps
for (nm in c("restrictRaster")){ if (exists(nm, inherits = FALSE)){ rm(list = nm) }}
if (!is.null(tmp_restrict)) {
  if (file.exists(tmp_restrict)) {unlink(tmp_restrict, force = TRUE)}
}
gc()

# Finalize
saveDatasheet(myScenario, spatialOutputsSheet, "wisdm_OutputSpatial")
progressBar(type = "end")
