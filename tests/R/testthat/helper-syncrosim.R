## -------------------------------------------------------------------------
## Shared test helpers: locating src/*.R and detecting a SyncroSim install
## -------------------------------------------------------------------------

# Walk upward from `start` until a directory containing src/package.xml is
# found. Lets tests locate the repo root regardless of the working directory
# testthat happens to run from.
find_wisdm_repo_root <- function(start = getwd()) {
  path <- normalizePath(start, mustWork = TRUE)
  repeat {
    if (file.exists(file.path(path, "src", "package.xml"))) {
      return(path)
    }
    parent <- dirname(path)
    if (parent == path) {
      stop("Could not locate the wisdm repo root (src/package.xml) above ", start)
    }
    path <- parent
  }
}

wisdm_src_dir <- file.path(find_wisdm_repo_root(), "src")

# src/00-helper-functions.R itself does
# `source(file.path(Sys.getenv("ssim_package_directory"), "00-constants.R"))`
# at its top level -- exactly as SyncroSim sets it up when invoking a
# transformer -- so this env var must be set for that (and any other
# src/*.R file that reads it) to source correctly under testthat.
Sys.setenv(ssim_package_directory = wisdm_src_dir)

# Source one of the wisdm src/*.R files into `env` for unit testing.
#
# - Ensures 00-constants.R (nodataValue, backgroundValue, ...) has been
#   loaded, since several src/*.R files rely on it being defined globally
#   rather than sourcing it themselves.
# - Stubs `updateRunLog()` (normally supplied by rsyncrosim at runtime) so
#   functions that log warnings can be exercised without a live scenario.
source_wisdm_r <- function(filename, env = parent.frame()) {
  if (!exists("nodataValue", envir = globalenv(), inherits = FALSE)) {
    sys.source(file.path(wisdm_src_dir, "00-constants.R"), envir = globalenv())
  }
  if (!exists("updateRunLog", envir = env, inherits = FALSE)) {
    assign("updateRunLog", function(...) invisible(NULL), envir = env)
  }
  sys.source(file.path(wisdm_src_dir, filename), envir = env)
}

# TRUE if a SyncroSim session is reachable and the wisdm package is
# registered with it, so integration tests can self-skip on a machine
# without SyncroSim Studio + the wisdm package installed.
syncrosim_available <- function() {
  if (!requireNamespace("rsyncrosim", quietly = TRUE)) {
    return(FALSE)
  }
  tryCatch(
    {
      mySession <- rsyncrosim::session()
      pkgs <- rsyncrosim::package(mySession)
      "wisdm" %in% pkgs$name
    },
    error = function(e) FALSE
  )
}

# Restrict a scenario's pipeline to a specific, ordered set of transformer
# display names (see the `displayName` attribute of each <transformer> in
# src/package.xml) -- used by the integration tests to skip slow stages
# (hyperparameter tuning, the non-GLM model fits, and the ensemble step).
set_wisdm_pipeline <- function(myScenario, stageDisplayNames) {
  pipeline <- data.frame(
    StageNameId = stageDisplayNames,
    RunOrder = seq_along(stageDisplayNames)
  )
  rsyncrosim::saveDatasheet(myScenario, pipeline, "core_Pipeline")
}

# The fast subset of the wisdm pipeline used by the integration tests:
# data prep through a single (GLM) model fit, then spatial prediction.
# Skips HyperparameterTuning (stage 7, only needed for RF/BRT) and
# EnsembleModel (stage 10, needs more than one fitted model).
wisdm_fast_pipeline_stages <- c(
  "1 - Prepare Multiprocessing",
  "2 - Spatial Data Preparation",
  "3 - Site Data Preparation",
  "4 - Background Data Generation",
  "5 - Prepare Training/Testing Data",
  "6 - Variable Reduction",
  "8 - Generalized Linear Model",
  "9 - Apply Model"
)
