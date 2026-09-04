## Integration test B -- builds a wisdm library entirely from scratch (no
## hosted template, no network access) with a tiny synthetic template raster,
## two covariate rasters, and a handful of field data points, then runs a
## fast (GLM-only) pipeline end-to-end.
##
## This is the environment-independent counterpart to
## test-integration-template-library.R: it only needs SyncroSim + the wisdm
## package (and its conda environments) to be installed locally, with no
## dependency on the "WISDM Example" online template being reachable/current.
##
## Requires SyncroSim Studio + the wisdm package to be installed and
## registered with rsyncrosim; self-skips otherwise.

test_that("a minimal synthetic library runs to completion for a fast (GLM-only) pipeline", {
  skip_if_not(syncrosim_available(), "SyncroSim + the wisdm package are not available")
  skip_if_not_installed("terra")

  library(rsyncrosim)
  library(terra)

  workDir <- file.path(tempdir(), paste0("wisdm-synthetic-test-", Sys.getpid()))
  dir.create(workDir, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(workDir, recursive = TRUE, force = TRUE), add = TRUE)

  # --- Build a tiny synthetic landscape -----------------------------------
  # A 30x30 cell template raster and two continuous covariates. Covariate 1
  # carries the (synthetic) habitat signal; covariate 2 is noise.

  set.seed(1)
  templateRast <- rast(nrows = 30, ncols = 30, xmin = 0, xmax = 30, ymin = 0, ymax = 30, crs = "EPSG:32610")
  values(templateRast) <- 1L

  covariate1 <- templateRast
  values(covariate1) <- as.vector(t(matrix(1:30, nrow = 30, ncol = 30))) + rnorm(900, sd = 0.5)
  covariate2 <- templateRast
  values(covariate2) <- rnorm(900)

  templatePath <- file.path(workDir, "template.tif")
  covariate1Path <- file.path(workDir, "covariate1.tif")
  covariate2Path <- file.path(workDir, "covariate2.tif")
  writeRaster(templateRast, templatePath, overwrite = TRUE, datatype = "INT2S")
  writeRaster(covariate1, covariate1Path, overwrite = TRUE)
  writeRaster(covariate2, covariate2Path, overwrite = TRUE)

  # Sample 30 "presence" cells where covariate1 is high and 30 "absence"
  # cells where covariate1 is low, so GLM has a learnable signal.
  covVals <- values(covariate1)[, 1]
  presenceCells <- order(covVals, decreasing = TRUE)[1:30]
  absenceCells <- order(covVals, decreasing = FALSE)[1:30]
  presenceXY <- xyFromCell(covariate1, presenceCells)
  absenceXY <- xyFromCell(covariate1, absenceCells)

  fieldData <- data.frame(
    X = c(presenceXY[, "x"], absenceXY[, "x"]),
    Y = c(presenceXY[, "y"], absenceXY[, "y"]),
    Response = c(rep(1L, 30), rep(0L, 30))
  )

  # --- Build the library ---------------------------------------------------

  mySession <- session()
  myLibrary <- ssimLibrary(
    name = file.path(workDir, "wisdm-synthetic-test.ssim"),
    session = mySession,
    packages = "wisdm",
    overwrite = TRUE
  )
  myProject <- rsyncrosim::project(myLibrary, project = "Definitions")

  covariatesSheet <- datasheet(myProject, "wisdm_Covariates")
  covariatesSheet <- addRow(
    covariatesSheet,
    data.frame(CovariateName = c("covariate1", "covariate2"), IsCategorical = c(FALSE, FALSE))
  )
  saveDatasheet(myProject, covariatesSheet, "wisdm_Covariates")

  myScenario <- scenario(myProject, scenario = "Synthetic test scenario")

  saveDatasheet(myScenario, data.frame(RasterFilePath = templatePath), "wisdm_TemplateRaster")
  saveDatasheet(myScenario, fieldData, "wisdm_FieldData")

  covariateDataSheet <- datasheet(myScenario, "wisdm_CovariateData", lookupsAsFactors = FALSE)
  covariateDataSheet <- addRow(
    covariateDataSheet,
    data.frame(
      CovariatesID = c("covariate1", "covariate2"),
      RasterFilePath = c(covariate1Path, covariate2Path)
    )
  )
  saveDatasheet(myScenario, covariateDataSheet, "wisdm_CovariateData")

  saveDatasheet(
    myScenario,
    data.frame(GenerateBackgroundSites = FALSE),
    "wisdm_BackgroundDataOptions"
  )
  saveDatasheet(
    myScenario,
    data.frame(SplitData = FALSE, CrossValidate = FALSE),
    "wisdm_ValidationOptions"
  )
  saveDatasheet(
    myScenario,
    data.frame(
      SelectionMethod = 1, # Automatic (Variance Inflation Factor) -- avoids interactive plotting
      DisplayHighestCorrelations = FALSE,
      CorrelationThreshold = 0.7,
      NumberOfPlots = 2,
      VIFThreshold = 10
    ),
    "wisdm_CovariateSelectionOptions"
  )
  saveDatasheet(
    myScenario,
    data.frame(
      SelectBestPredictors = FALSE,
      SimplificationMethod = 0,
      ConsiderSquaredTerms = FALSE,
      ConsiderInteractions = FALSE
    ),
    "wisdm_GLM"
  )

  set_wisdm_pipeline(myScenario, wisdm_fast_pipeline_stages)

  # --- Run and verify --------------------------------------------------------

  myResultScenario <- run(myScenario)

  outputModel <- datasheet(myResultScenario, "wisdm_OutputModel")
  expect_gte(nrow(outputModel), 1)

  outputSpatial <- datasheet(myResultScenario, "wisdm_OutputSpatial")
  expect_gte(nrow(outputSpatial), 1)
  for (rasterPath in na.omit(outputSpatial$ProbabilityRaster)) {
    expect_true(file.exists(rasterPath))
    outRast <- rast(rasterPath)
    expect_false(all(is.na(values(outRast))))
  }
})
