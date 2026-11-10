## -------------------------
## wisdm - test setup
## ApexRMS, October 2026
## -------------------------

# testthat runs any file starting with "setup" before the tests, with the
# working directory set to tests/r. The helper functions read the package
# folder from an environment variable that SyncroSim normally sets, so we set
# it here to point at src/.

srcDir <- normalizePath(file.path("..", "..", "src"))
Sys.setenv(ssim_package_directory = srcDir)

source(file.path(srcDir, "00-helper-functions.R"))
