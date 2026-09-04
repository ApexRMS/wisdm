## Integration test A -- builds a library from the wisdm package's hosted
## "WISDM Example" online template (see docs/getting_started.md, Step 2) and
## runs its pre-configured "Brewer's Sparrow" scenario end-to-end.
##
## Requires SyncroSim Studio + the wisdm package (and its conda environments)
## to be installed and registered with rsyncrosim; self-skips otherwise.
##
## NOTE: confirm `ssimLibrary(..., template = ...)` against the rsyncrosim
## version actually installed (src/wisdm-r-conda.yml pins r-rsyncrosim) --
## some releases exposed this as `library()` instead. Run `?ssimLibrary` /
## `?library` in that environment if this test errors on an unknown argument.

test_that("the 'WISDM Example' template runs to completion for a fast (GLM-only) pipeline", {
  skip_if_not(syncrosim_available(), "SyncroSim + the wisdm package are not available")

  library(rsyncrosim)

  workDir <- file.path(tempdir(), paste0("wisdm-template-test-", Sys.getpid()))
  dir.create(workDir, showWarnings = FALSE, recursive = TRUE)
  on.exit(unlink(workDir, recursive = TRUE, force = TRUE), add = TRUE)

  mySession <- session()
  myLibrary <- ssimLibrary(
    name = file.path(workDir, "wisdm-template-test.ssim"),
    session = mySession,
    packages = "wisdm",
    template = "WISDM Example",
    overwrite = TRUE
  )
  myProject <- rsyncrosim::project(myLibrary, project = "Definitions")
  myScenario <- scenario(myProject, scenario = "Brewer's Sparrow")

  set_wisdm_pipeline(myScenario, wisdm_fast_pipeline_stages)

  myResultScenario <- run(myScenario)

  outputModel <- datasheet(myResultScenario, "wisdm_OutputModel")
  expect_gte(nrow(outputModel), 1)

  outputSpatial <- datasheet(myResultScenario, "wisdm_OutputSpatial")
  expect_gte(nrow(outputSpatial), 1)
  for (rasterPath in na.omit(outputSpatial$ProbabilityRaster)) {
    expect_true(file.exists(rasterPath))
  }
})
