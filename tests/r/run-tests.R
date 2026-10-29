## -------------------------
## wisdm - run R tests
## ApexRMS, October 2026
## -------------------------

# Run from the repository root:
#   Rscript tests/r/run-tests.R

testDir <- file.path("tests", "r")

if (!dir.exists(testDir)) {
  stop("Run this script from the wisdm repository root folder.")
}

testthat::test_dir(testDir, stop_on_failure = TRUE)
